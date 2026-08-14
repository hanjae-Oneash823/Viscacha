"""Audited four-gate selection matrix for trial_failure_candidate docking pairs.

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
2. Domain disruption is reported as context (domains_touched_kept /
   domain_disruption), not used as a gate.  Structure availability is Gate D.

Bug fixes vs the six-gate version
---------------------------------
* Domain context uses the true changed-residue set rather than the broad
  changed_aa_start/end envelope, which can mark every domain as affected on
  an N-truncation.
* Gate B briefly required an Ensembl-annotated coding CDS (PC_CDS). That has
  since been dropped again -- all trial_failure alts are protein-coding
  (PC_CDS/PC_CDS_ND/PC_UTR are all coding, just differing in CDS-annotation
  provenance), so the biotype clause was excluding real candidates rather
  than filtering noise. Gate B is now length-retention + no-premature-stop
  only; premature_stop still fires on only 1/163 rows and is effectively
  inert on its own.
* Gate C tested for drug NAME strings, which is exactly equivalent to
  phase >= 1 across all 149 pairs -- and is applied to a universe J4 already
  filtered on DGIdb (j4_gate._has_drug_evidence ORs it in; 106/149 pairs enter
  ONLY via that clause). It is relabelled to say what it measures, not
  rescored. DGIdb is legitimate for RECRUITMENT but unusable as a ligand list.
* Dedup across cell types took whichever row came first, which is
  order-dependent -- CHD1-205 straddled the gate A threshold on row order
  alone. It now takes the max.

Domain disruption is not a selection criterion: the fraction affected is
bimodal, and any cutoff would silently become the real selector. Pocket
overlap is the right test and needs m2c_pocket, not this table.

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

# Pfam overlap is descriptive context only; Gate D is structure availability.
GATES = list("ABCD")
MIN_USAGE_DELTA = 0.10        # gate A -- rise
MIN_ALT_LEVEL_AD = 0.25       # gate A -- OR absolute AD level
MIN_KEPT_FRAC = 0.50          # gate B

# (letter, short name, one-line definition) -- shared with plot_gate_overlap
GATE_LABEL = {
    "A": ("A", "coherent switch",
          "alt rises ≥ 0.10 usage\nor is > 0.25 of AD usage"),
    "B": ("B", "intact protein", "≥ 50% kept, no premature stop"),
    "C": ("C", "curated drug", "ChEMBL or OpenTargets phase ≥ 1"),
    "D": ("D", "structures ready", "canonical + alt folded"),
}

# Mechanism-only working set: every gate except the drug-evidence one. Gate C
# is a literature-coverage filter, not a property of the protein, so this is
# the set worth carrying into pocket analysis.
MECHANISM_GATES = "ABD"

_CIGAR_RE = re.compile(r"(\d+)([=XID])")

DETAIL_COLUMNS = [
    "gene_name", "cell_type", "n_cell_types", "transcript_name",
    "alt_transcript_name", "alt_rank",
    # Retained for downstream, non-gating structural analyses.  Gate-matrix
    # plots do not display these long sequences, but cached model lookup is
    # sequence-hash based so an audit table must be able to recover them.
    "canonical_protein_seq", "alt_protein_seq",
    # usage: `delta_usage` / AD / Control describe the CANONICAL transcript,
    # `alt_usage_*` the alt. Every trial_failure hit is CT_enriched, so the
    # canonical's delta_usage is negative by construction.
    "delta_usage", "AD", "Control",
    "alt_usage_pct_AD", "alt_usage_pct_AD_max", "alt_usage_pct_control", "alt_usage_delta",
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
    """One row per distinct protein pair, with the four active gates as bool columns.

    Gate D re-scans STRUCTURE_CACHE_DIR on every call, so it tracks whatever
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
    # denominator by 55. They are not alternative-protein candidates and the
    # canonical is CT-enriched by construction, so dropping them corrects the
    # per-gate totals without discarding a meaningful comparison.
    tf = tf[~(tf["is_canonical"] == True)]  # noqa: E712 -- NaN-safe

    # One row per protein pair, carrying the MAX usage delta over cell types.
    key = ["canonical_protein_seq", "alt_protein_seq"]
    tf = tf.sort_values("alt_usage_delta", ascending=False)
    pairs = tf.drop_duplicates(key).copy()
    n_ct = tf.groupby(key)["cell_type"].nunique()
    pairs["n_cell_types"] = pairs.set_index(key).index.map(n_ct)
    # Gate A's OR clause needs the max AD level over cell types too, not just
    # whichever cell type happened to carry the max delta (the row `pairs`
    # already collapsed to) -- a pair can peak on level in one cell type and
    # on delta in another.
    max_level_ad = tf.groupby(key)["alt_usage_pct_AD"].max()
    pairs["alt_usage_pct_AD_max"] = pairs.set_index(key).index.map(max_level_ad)

    overlap = pd.DataFrame([_domain_overlap(r, pfam) for _, r in pairs.iterrows()])
    pairs = pd.concat([pairs.reset_index(drop=True), overlap], axis=1)
    pairs["can_aa_len"] = pairs["canonical_protein_seq"].str.len()
    pairs["alt_aa_len"] = pairs["alt_protein_seq"].str.len()

    def num(col: str) -> pd.Series:
        return pd.to_numeric(pairs[col], errors="coerce").fillna(0)

    kept_frac = (pairs["alt_protein_seq"].str.len()
                 / pairs["canonical_protein_seq"].str.len().replace(0, np.nan))

    gates = pd.DataFrame({
        "A": (pairs["alt_usage_delta"] >= MIN_USAGE_DELTA)
             | (pairs["alt_usage_pct_AD_max"] > MIN_ALT_LEVEL_AD),
        "B": ((kept_frac >= MIN_KEPT_FRAC)
              & (~(pairs["premature_stop"] == True))),  # noqa: E712 -- NaN-safe
        "C": (num("chembl_max_phase") >= 1) | (num("ot_max_phase") >= 1),
        "D": (pairs["canonical_protein_seq"].map(_canonical_folded)
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
    print("all active gates:", sorted(set(m.loc[m[GATES].all(axis=1), "gene_name"])) or "(none)")
    mech = m[m[list(MECHANISM_GATES)].all(axis=1)]
    print(f"mechanism-only ({'+'.join(MECHANISM_GATES)}):",
          sorted(set(mech.gene_name)) or "(none)")
    waiting = m[m[["A", "B"]].all(axis=1) & ~m["D"]]
    print("waiting only on folding:", sorted(set(waiting.gene_name)) or "(none)")
