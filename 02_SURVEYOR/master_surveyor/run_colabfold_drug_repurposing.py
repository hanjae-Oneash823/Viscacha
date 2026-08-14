"""ColabFold pass for the drug_repurposing_candidate group (168 genes / 175
hits, already one row per hit -- see m0_select.representative_row, unlike
trial_failure_candidate's long-format alt_rank rows).

Canonical sequences: AlphaFold DB REST lookup first (get_canonical_structure
checks uniprot_acc against AFDB, free/instant/no GPU) -- these are the plain
UniProt/MANE sequence by construction, so AFDB is expected to already have
almost all of them. ColabFold only runs locally on an AFDB miss.

Alt/isoform sequences: always folded locally on GPU1, same no-template,
elevated-recycle (20), multi-seed (5) ensemble protocol as
run_colabfold_all_alts.py used for trial_failure_candidate
(config.COLABFOLD_ALT_ARGS / CUDA_VISIBLE_DEVICES).

Dedup is by sequence, not by hit -- a handful of genes recur across
multiple cell types with an identical canonical/alt sequence, and
get_alt_structure()/get_canonical_structure() cache by seq_hash anyway, so
this just avoids redundant fasta writes/log lines for those.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from master_surveyor import m2_structures
from master_surveyor.config import HITS_CSV


def _log(msg: str) -> None:
    print(f"[run_colabfold_drug_repurposing] {msg}", file=sys.stderr, flush=True)


def main() -> None:
    df = pd.read_csv(HITS_CSV)
    dr = df[df["master_group"] == "drug_repurposing_candidate"]
    n_hits = dr.groupby(["gene_name", "cell_type", "hit_ENST_ID"]).ngroups
    _log(f"{len(dr):,} rows across {n_hits} drug_repurposing hits")

    canon = dr[["gene_name", "uniprot_acc", "canonical_protein_seq"]].drop_duplicates(subset=["canonical_protein_seq"])
    alts = dr[["gene_name", "alt_rank", "alt_protein_seq"]].drop_duplicates(subset=["alt_protein_seq"])
    _log(f"{len(canon):,} unique canonical sequences, {len(alts):,} unique alt sequences to process")

    _log("--- canonical structures (AlphaFold DB lookup; ColabFold only on a miss) ---")
    canon_done, canon_failed, canon_afdb, canon_folded = 0, 0, 0, 0
    for i, (_, row) in enumerate(canon.iterrows(), 1):
        gene, acc = row["gene_name"], row["uniprot_acc"]
        t0 = time.monotonic()
        try:
            result = m2_structures.get_canonical_structure(acc, row["canonical_protein_seq"])
            elapsed = time.monotonic() - t0
            if result.source == "afdb":
                canon_afdb += 1
            else:
                canon_folded += 1
            _log(f"  [{i}/{len(canon)}] {gene} ({acc}) ok in {elapsed:.0f}s -- source={result.source}")
            canon_done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(canon)}] {gene} ({acc}) FAILED after {elapsed:.0f}s: {exc}")
            canon_failed += 1
    _log(
        f"canonical done: {canon_done:,} ok ({canon_afdb:,} from AFDB, {canon_folded:,} folded locally), "
        f"{canon_failed:,} failed"
    )

    _log("--- alt/isoform structures (always local ColabFold, GPU1) ---")
    alt_done, alt_failed = 0, 0
    for i, (_, row) in enumerate(alts.iterrows(), 1):
        seq = row["alt_protein_seq"]
        gene = row["gene_name"]
        rank = row["alt_rank"]
        t0 = time.monotonic()
        try:
            result = m2_structures.get_alt_structure(seq)
            elapsed = time.monotonic() - t0
            _log(
                f"  [{i}/{len(alts)}] {gene} rank{rank} ({len(seq)} aa) ok in {elapsed:.0f}s -- "
                f"{len(result.seed_models)} seed models, esmfold_error={result.esmfold_error}"
            )
            alt_done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(alts)}] {gene} rank{rank} ({len(seq)} aa) FAILED after {elapsed:.0f}s: {exc}")
            alt_failed += 1

    _log(f"alt done: {alt_done:,} ok, {alt_failed:,} failed")


if __name__ == "__main__":
    main()
