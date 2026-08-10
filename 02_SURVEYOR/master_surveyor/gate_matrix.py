"""Audited five-gate selection matrix for trial_failure_candidate docking pairs.

This supersedes the earlier six-gate framework (recorded in the
project_ds_docking_selection memory). Two structural changes and four bug
fixes separate the two; both are load-bearing, so the old funnel
(55 -> 34 -> 28 -> 7 -> 2 -> 2 -> 2, survivors IFNAR1 + MLST8) is NOT
reproducible from this module and should not be treated as a regression.

Structural changes
------------------
1. Substrate is the UNCOLLAPSED set of distinct
   (canonical_protein_seq, alt_protein_seq) pairs -- 94, not 55.
   m0_select.representative_row() keeps one alt per hit, which is correct for
   the export path but hides genuine docking candidates when several alts of
   the same canonical differ. Cell type is still not a docking unit, so rows
   sharing a protein pair are collapsed carrying the MAX alt_usage_delta, and
   the canonical's own row is dropped (see build_gate_matrix).
2. The old gate E ("touched domain still detected in the alt") is dropped and
   structures-ready is renumbered F -> E. Domain retention is still reported
   as domains_touched_kept / domain_disruption -- context for reading the
   shortlist, not a filter.

Bug fixes vs the six-gate version
---------------------------------
* Gate D used changed_aa_start/end, an ENVELOPE that marks every domain as
  affected on any N-truncation. It now re-runs the edlib alignment to recover
  the true changed-residue set and intersects that with real Pfam coordinates.
  14 of 56 D-passers under the old test spanned >0.8 of the protein.
* Gate B now requires an Ensembl-annotated coding CDS (PC_CDS). This retires
  three separate defects at once: premature_stop fires on 1/163 rows and is
  effectively inert; PC_CDS_ND carries SQANTI-PREDICTED ORFs (IFNAR1-216,
  POMT2-240); and PC_UTR conflates "no CDS difference" with "no CDS data"
  (ZDHHC3-202).
* Gate C tested for drug NAME strings, which is exactly equivalent to
  phase >= 1 across all 149 pairs -- and is applied to a universe J4 already
  filtered on DGIdb (j4_gate._has_drug_evidence ORs it in; 106/149 pairs enter
  ONLY via that clause). It is relabelled to say what it measures, not
  rescored. DGIdb is legitimate for RECRUITMENT but unusable as a ligand list.
* Dedup across cell types took whichever row came first, which is
  order-dependent -- CHD1-205 straddled the gate A threshold on row order
  alone. It now takes the max.

Deliberately NOT changed: gate D stays "any true changed residue lands in a
Pfam domain" rather than a fraction-of-domain threshold. That metric is
bimodal -- a domain is either obliterated by a truncation (~0.95) or nicked by
a local indel (~0.03), nothing between -- so any cutoff silently becomes the
real selector. It is also the wrong proxy for pocket relevance: PRKG2's
29-residue indel is 7% of its kinase domain and could still reshape the ATP
site. Pocket overlap is the right test and needs m2c_pocket, not this table.

Requires edlib (present in the oneash_dtu env, same as j2_protein_diff).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import edlib
import numpy as np
import pandas as pd

from master_surveyor import m2_structures
from master_surveyor.config import HITS_CSV, REPO_ROOT, STRUCTURE_CACHE_DIR

PFAM_JSON = REPO_ROOT / "outputs/junior_surveyor/cache/j2_pfam_hits.json"

GATES = list("ABCDE")
MIN_USAGE_DELTA = 0.10        # gate A
MIN_KEPT_FRAC = 0.50          # gate B

# (letter, short name, one-line definition) -- shared with plot_gate_overlap
GATE_LABEL = {
    "A": ("A", "coherent switch", "alt rises ≥ 0.10 usage (max over cell types)"),
    "B": ("B", "intact protein", "≥ 50% kept, Ensembl-annotated CDS"),
    "C": ("C", "curated drug", "ChEMBL or OpenTargets phase ≥ 1"),
    "D": ("D", "domain hit", "changed residues land in a Pfam domain"),
    "E": ("E", "structures ready", "canonical + alt folded (was gate F)"),
}

# Mechanism-only working set: every gate except the drug-evidence one. Gate C
# is a literature-coverage filter, not a property of the protein, so this is
# the set worth carrying into pocket analysis.
MECHANISM_GATES = "ABDE"

_CIGAR_RE = re.compile(r"(\d+)([=XID])")

DETAIL_COLUMNS = [
    "gene_name", "cell_type", "n_cell_types", "transcript_name",
    "alt_transcript_name", "alt_rank",
    # usage: `delta_usage` / AD / Control describe the CANONICAL transcript,
    # `alt_usage_*` the alt. Every trial_failure hit is CT_enriched, so the
    # canonical's delta_usage is negative by construction.
    "delta_usage", "AD", "Control",
    "alt_usage_pct_AD", "alt_usage_pct_control", "alt_usage_delta",
    "alt_biotype_class", "alt_cds_source", "protein_change_type",
    "can_aa_len", "alt_aa_len", "pct_identity", "protein_length_diff",
    "changed_aa_start", "changed_aa_end", "changed_aa_fraction",
    "n_changed", "true_changed_frac", "top_domain", "domain_disruption",
    "domains_touched", "domains_touched_kept", "n_touched", "n_touched_kept",
    "affected_domain", "chembl_bioactive_compounds", "chembl_best_pchembl",
    "dgidb_interactions", "chembl_max_phase", "ot_max_phase", "uniprot_acc",
]


def changed_positions(canonical: str, alt: str) -> set[int]:
    """True set of 1-based canonical residue positions altered in `alt`.

    Mirrors j2_protein_diff._align_proteins' CIGAR walk exactly (query=canonical,
    '=' match, 'X' mismatch, 'I' present in canonical only, 'D' present in alt
    only anchored to the current canonical position clamped to >= 1).
    """
    if not canonical or not alt or canonical == alt:
        return set()
    ops = _CIGAR_RE.findall(edlib.align(canonical, alt, mode="NW",
                                        task="path")["cigar"])
    can_pos, changed = 0, set()
    for n, op in ops:
        n = int(n)
        if op == "=":
            can_pos += n
        elif op in ("X", "I"):
            for _ in range(n):
                can_pos += 1
                changed.add(can_pos)
        elif op == "D":
            changed.add(max(can_pos, 1))
    return changed


def _canonical_folded(seq: str) -> bool:
    out_dir = STRUCTURE_CACHE_DIR / m2_structures.seq_hash(str(seq))
    if (out_dir / "canonical_afdb.pdb").exists():
        return True
    cf_dir = out_dir / "canonical_colabfold"
    return bool(cf_dir.exists() and list(cf_dir.glob("*_unrelaxed_rank_001_*.pdb")))


def _alt_folded(seq: str) -> bool:
    cf_dir = STRUCTURE_CACHE_DIR / m2_structures.seq_hash(str(seq)) / "alt_colabfold"
    return bool(cf_dir.exists() and list(cf_dir.glob("*_unrelaxed_rank_*.pdb")))


def _domain_overlap(row: pd.Series, pfam: dict) -> dict:
    """True changed-residue overlap against real Pfam coordinates."""
    chg = changed_positions(row.canonical_protein_seq, row.alt_protein_seq)
    canon_doms = pfam.get(f"{row.canonical_enst}__canonical", [])
    alt_accs = {d["acc"] for d in pfam.get(f"{row.ENST_ID}__alt", [])}

    touched, touched_kept = [], []
    best_name, best_disr = "", 0.0
    for d in canon_doms:
        span = set(range(d["start"], d["end"] + 1))
        n_hit = len(span & chg)
        if not span or not n_hit:
            continue
        touched.append(d["name"])
        if d["acc"] in alt_accs:
            touched_kept.append(d["name"])
        if n_hit / len(span) > best_disr:
            best_disr, best_name = n_hit / len(span), d["name"]

    n_can = len(row.canonical_protein_seq)
    return {
        "n_changed": len(chg),
        "true_changed_frac": len(chg) / n_can if n_can else 0.0,
        "domains_touched": ", ".join(dict.fromkeys(touched)),
        "domains_touched_kept": ", ".join(dict.fromkeys(touched_kept)),
        "n_touched": len(set(touched)),
        "n_touched_kept": len(set(touched_kept)),
        "top_domain": best_name,
        "domain_disruption": round(best_disr, 4),
        "n_canon_domains": len(canon_doms),
    }


def build_gate_matrix() -> pd.DataFrame:
    """One row per distinct protein pair, with the five gates as bool columns.

    Gate E re-scans STRUCTURE_CACHE_DIR on every call, so it tracks whatever
    ColabFold has finished at the moment of the run.
    """
    pfam = json.loads(PFAM_JSON.read_text())
    hits = pd.read_csv(HITS_CSV, low_memory=False)
    tf = hits[hits["master_group"] == "trial_failure_candidate"].copy()
    tf["canonical_protein_seq"] = tf["canonical_protein_seq"].fillna("")
    tf["alt_protein_seq"] = tf["alt_protein_seq"].fillna("")
    tf = tf[(tf.canonical_protein_seq != "") & (tf.alt_protein_seq != "")]

    # The canonical appears in its own alt-ranking table (is_canonical, with
    # alt_ENST_ID == canonical_enst and protein_change_type "identical"). A
    # protein compared against itself is not a docking pair: 60 such rows
    # collapse to exactly one self-pair per gene and would inflate the
    # denominator by 55. They never survive -- they fail gate D (nothing
    # changed) and gate A (the canonical is CT_enriched by construction, so
    # its own delta is negative) -- so dropping them leaves the funnel
    # identical while correcting the per-gate totals for B, C and E.
    tf = tf[~(tf["is_canonical"] == True)]  # noqa: E712 -- NaN-safe

    # One row per protein pair, carrying the MAX usage delta over cell types.
    key = ["canonical_protein_seq", "alt_protein_seq"]
    tf = tf.sort_values("alt_usage_delta", ascending=False)
    pairs = tf.drop_duplicates(key).copy()
    n_ct = tf.groupby(key)["cell_type"].nunique()
    pairs["n_cell_types"] = pairs.set_index(key).index.map(n_ct)

    overlap = pd.DataFrame([_domain_overlap(r, pfam) for _, r in pairs.iterrows()])
    pairs = pd.concat([pairs.reset_index(drop=True), overlap], axis=1)
    pairs["can_aa_len"] = pairs["canonical_protein_seq"].str.len()
    pairs["alt_aa_len"] = pairs["alt_protein_seq"].str.len()

    def num(col: str) -> pd.Series:
        return pd.to_numeric(pairs[col], errors="coerce").fillna(0)

    kept_frac = (pairs["alt_protein_seq"].str.len()
                 / pairs["canonical_protein_seq"].str.len().replace(0, np.nan))

    gates = pd.DataFrame({
        "A": pairs["alt_usage_delta"] >= MIN_USAGE_DELTA,
        "B": ((kept_frac >= MIN_KEPT_FRAC)
              & (~(pairs["premature_stop"] == True))  # noqa: E712 -- NaN-safe
              & (pairs["alt_biotype_class"] == "PC_CDS")),
        "C": (num("chembl_max_phase") >= 1) | (num("ot_max_phase") >= 1),
        "D": pairs["n_touched"] > 0,
        "E": (pairs["canonical_protein_seq"].map(_canonical_folded)
              & pairs["alt_protein_seq"].map(_alt_folded)),
    }).fillna(False).astype(bool)

    return pd.concat([pairs[DETAIL_COLUMNS].reset_index(drop=True), gates], axis=1)


def funnel(gates: pd.DataFrame) -> tuple[list[str], list[int]]:
    """Cumulative A->E survivor counts, plus the ungated total."""
    labels, counts = ["all trial_failure pairs"], [len(gates)]
    mask = pd.Series(True, index=gates.index)
    for g in GATES:
        mask &= gates[g]
        labels.append("A" if g == "A" else f"A–{g}")
        counts.append(int(mask.sum()))
    return labels, counts


if __name__ == "__main__":
    m = build_gate_matrix()
    print(f"pairs={len(m)} genes={m.gene_name.nunique()}")
    print("per-gate:", {g: int(m[g].sum()) for g in GATES})
    labels, counts = funnel(m)
    print("funnel:", " → ".join(str(c) for c in counts))
    print("all five:", sorted(set(m.loc[m[GATES].all(axis=1), "gene_name"])) or "(none)")
    mech = m[m[list(MECHANISM_GATES)].all(axis=1)]
    print(f"mechanism-only ({'+'.join(MECHANISM_GATES)}):",
          sorted(set(mech.gene_name)) or "(none)")
    waiting = m[m[["A", "B", "D"]].all(axis=1) & ~m["E"]]
    print("waiting only on folding:", sorted(set(waiting.gene_name)) or "(none)")
