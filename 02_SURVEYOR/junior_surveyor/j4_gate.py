"""J4_GATE — unified selection for master_surveyor, branch-specific criteria.

new_target_candidate (one alt_rank=0 row per hit):
    passes_to_master = protein_change(alt_rank==0)
    (no drug-evidence requirement -- the downstream plan screens the
    AD-enriched isoform against an external drug panel regardless of whether
    this gene already has a known ChEMBL/DGIdb/Open Targets binder; requiring
    prior drug evidence would discard the most interesting "genuinely novel
    target" cases. chembl_max_phase/dgidb_interactions/ot_max_phase/drug_names
    are still attached as informational columns for tiering, just not gated on.)

trial_failure_candidate (alt_rank=1..N ranked alternates surviving J1c's
    biotype + >=5% AD-usage filter, THEN J2's protein-identical-decoy filter
    -- J1b itself ranks and fetches every isoform with no cutoff; J1c drops
    low-quality/low-usage rows; J2's _drop_protein_identical_alts drops
    non-canonical rows that turn out protein-identical to canonical and any
    hit left with none remaining):
    passes_to_master = ANY surviving ranked alternate has a real protein change
                        AND drug_evidence(gene)
    (the protein-change half is trivially true by construction for every hit
    that reaches J4 -- J2 already removed the identical/no_sequence rows and
    the hits with nothing left -- so in practice this reduces to
    drug_evidence(gene). Kept as an explicit .any() rather than assumed True
    so the gate stays correct even if a future no-J2 code path skips that
    filter. drug_evidence stays required here -- the downstream plan needs an
    actual known drug (drug_names) to test against the AD-shifted alternate(s).)

Gate is computed once per hit (gene_name, cell_type) and applied to every
alt_rank row of that hit, so a downstream consumer can filter on one boolean
column regardless of branch.

`selected_for_next_stage` is kept as an alias of `passes_to_master` for
backward compatibility with master_surveyor, which has not yet been updated
for the long-format (multi-row-per-hit) trial_failure shape.

master_group (3-way split of passes_to_master hits, for master_surveyor):
    trial_failure_candidate   -- the existing trial_failure branch, unchanged
                                  (drug_evidence already required by its gate)
    drug_repurposing_candidate -- new_target hits whose gene has ANY existing
                                  chemical-matter evidence (broader than the
                                  gate's own drug_evidence: also counts a bare
                                  ChEMBL target entry and a Pharos ligand, not
                                  just a clinical-phase drug/DGIdb hit -- the
                                  point is "something to repurpose/screen from",
                                  not "already an approved drug")
    novel_target_candidate    -- new_target hits with no chemical-matter
                                  evidence at all: genuinely undrugged targets
    "" (empty)                 -- dropped hits (passes_to_master == False)
"""

from __future__ import annotations

import pandas as pd

from junior_surveyor.config import NULL_PROTEIN_CHANGE_TYPES

_HIT_KEY = ["gene_name", "cell_type"]


def _has_protein_change(row: pd.Series) -> bool:
    return row.get("protein_change_type", "no_sequence") not in NULL_PROTEIN_CHANGE_TYPES


def _has_drug_evidence(row: pd.Series) -> bool:
    phase    = int(row.get("chembl_max_phase", 0) or 0)
    dgidb    = int(row.get("dgidb_interactions", 0) or 0)
    ot_phase = int(row.get("ot_max_phase", 0) or 0)
    return phase >= 1 or dgidb > 0 or ot_phase >= 1


def _has_any_chemical_evidence(row: pd.Series) -> bool:
    """Broader than _has_drug_evidence: any sign the gene is chemically
    tractable at all -- a bare ChEMBL target entry or a Pharos-listed ligand
    counts, not just a clinical-phase drug or curated DGIdb interaction.
    Used only for master_group's drug_repurposing/novel_target split, not
    for the passes_to_master gate itself.
    """
    is_druggable = bool(row.get("is_druggable", False))
    dgidb        = int(row.get("dgidb_interactions", 0) or 0)
    ot_phase     = int(row.get("ot_max_phase", 0) or 0)
    n_ligands    = int(row.get("pharos_n_ligands", 0) or 0)
    return is_druggable or dgidb > 0 or ot_phase >= 1 or n_ligands > 0


def run(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["alt_is_gate_driver"] = df.apply(_has_protein_change, axis=1)
    drug_evidence = df.apply(_has_drug_evidence, axis=1)

    nt_mask = df["candidate_group"] == "new_target_candidate"
    tf_mask = df["candidate_group"] == "trial_failure_candidate"

    # new_target: the single alt_rank=0 row IS the gate check.
    # No drug-evidence requirement -- see module docstring.
    nt = df[nt_mask & (df["alt_rank"] == 0)].copy()
    nt["_pass"] = nt["alt_is_gate_driver"]
    nt_pass_keys = set(map(tuple, nt.loc[nt["_pass"], _HIT_KEY].values))

    # trial_failure: pass if ANY ranked alternate is a gate driver AND the
    # gene has drug evidence (drug evidence is gene-level, identical across
    # all alt_rank rows of a hit, so .any() == the row-level value)
    tf = df[tf_mask].copy()
    tf["_drug"] = drug_evidence.loc[tf.index]
    tf_grp = tf.groupby(_HIT_KEY)
    tf_pass = (tf_grp["alt_is_gate_driver"].any() & tf_grp["_drug"].any())
    tf_pass_keys = set(tf_pass[tf_pass].index)

    def _mark(row: pd.Series) -> bool:
        key = (row["gene_name"], row["cell_type"])
        if row["candidate_group"] == "new_target_candidate":
            return key in nt_pass_keys
        if row["candidate_group"] == "trial_failure_candidate":
            return key in tf_pass_keys
        return False

    def _reason(row: pd.Series, passed: bool) -> str:
        if not passed:
            return ""
        if row["candidate_group"] == "new_target_candidate":
            return "new_target: protein_change (structural difference from canonical)"
        if row["candidate_group"] == "trial_failure_candidate":
            return "trial_failure: dominant_alt_structurally_different + drug_evidence"
        return ""

    df["passes_to_master"] = df.apply(_mark, axis=1)
    df["gate_branch"] = df["candidate_group"]
    df["gate_reason"] = df.apply(lambda r: _reason(r, r["passes_to_master"]), axis=1)

    # Backward-compat alias for master_surveyor (not yet updated for the
    # long-format trial_failure shape — see junior_surveyor redesign notes).
    df["selected_for_next_stage"] = df["passes_to_master"]

    # ── master_group: 3-way split of passes_to_master hits ──────────────────
    # trial_failure stays as-is; new_target splits into drug_repurposing vs
    # novel_target by chemical-tractability evidence (gene-level, so it's
    # identical across every alt_rank row of a hit -- see module docstring).
    chem_evidence = df.apply(_has_any_chemical_evidence, axis=1)

    def _master_group(row: pd.Series, has_chem: bool) -> str:
        if not row["passes_to_master"]:
            return ""
        if row["candidate_group"] == "trial_failure_candidate":
            return "trial_failure_candidate"
        if row["candidate_group"] == "new_target_candidate":
            return "drug_repurposing_candidate" if has_chem else "novel_target_candidate"
        return ""

    df["master_group"] = [
        _master_group(row, has_chem)
        for (_, row), has_chem in zip(df.iterrows(), chem_evidence)
    ]

    # NOTE: totals per branch use (candidate_group, gene_name, cell_type) so a
    # gene/cell_type appearing in both branches (different transcripts) isn't
    # collapsed into one when summing across branches.
    branch_key = ["candidate_group"] + _HIT_KEY
    n_hits_total = df[branch_key].drop_duplicates().shape[0]
    n_hits_pass  = df.loc[df["passes_to_master"], branch_key].drop_duplicates().shape[0]
    print(f"  [j4_gate] passes_to_master: {n_hits_pass:,}/{n_hits_total:,} hits")
    for branch in ["new_target_candidate", "trial_failure_candidate"]:
        b = df[df["candidate_group"] == branch]
        b_hits = b[_HIT_KEY].drop_duplicates().shape[0]
        b_pass = b.loc[b["passes_to_master"], _HIT_KEY].drop_duplicates().shape[0]
        print(f"  [j4_gate]   {branch}: {b_pass:,}/{b_hits:,} hits")

    print(f"  [j4_gate] master_group (passes_to_master hits only):")
    for group in ["trial_failure_candidate", "drug_repurposing_candidate", "novel_target_candidate"]:
        g_hits = df.loc[df["master_group"] == group, _HIT_KEY].drop_duplicates().shape[0]
        print(f"  [j4_gate]   {group}: {g_hits:,} hits")

    return df
