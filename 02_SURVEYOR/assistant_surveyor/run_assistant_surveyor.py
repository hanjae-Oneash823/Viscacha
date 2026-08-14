"""ASSISTANT_SURVEYOR — orchestrator.

Run from 02_SURVEYOR/:
    /home/welcome3/anaconda3/envs/oneash_dtu/bin/python -m assistant_surveyor.run_assistant_surveyor

Options:
    --no-cache   ignore cached API responses (re-fetch everything)
    --offline    use cache only; fail if cache absent for any batch
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

# Ensure 02_SURVEYOR (this package's parent) and repo root (for
# classify_hit_scenarios_mane, imported by l1_biotype) are on path regardless
# of invocation cwd.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from assistant_surveyor.config import CACHE_DIR, HITS_CSV, OUT_DIR
from assistant_surveyor import l1_biotype, l2_ad_prior, l3_opentargets, l4_uniprot, l5_consequences


def _ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


def _log(msg: str) -> None:
    print(f"[{_ts()}] {msg}", flush=True)


def _gene_summary(hits: pd.DataFrame) -> pd.DataFrame:
    """Collapse enriched hits to one row per gene."""
    proxy_order = {"C": 0, "D": 1, "frame": 2, "NMD": 3, "N": 4}

    def best_proxy(series):
        vals = set(series.dropna())
        for pt in ("C", "D", "frame", "NMD", "N"):
            if pt in vals:
                return pt
        return "N"

    agg = hits.groupby("gene_name").agg(
        ENSG_ID                = ("ENSG_ID", "first"),
        n_hit_transcripts      = ("transcript_name", "nunique"),
        n_hit_cell_types       = ("cell_type", "nunique"),
        max_abs_delta_usage    = ("delta_usage", lambda x: x.abs().max()),
        min_permutation_pval   = ("permutation_pval", "min"),
        best_proxy_type        = ("proxy_type", best_proxy),
        ad_prior_category      = ("ad_prior_category", "first"),
        ad_ot_score            = ("ad_ot_score", "max"),
        ad_ot_label            = ("ad_ot_label", lambda x: x.iloc[0]),
        has_structural_feat    = ("has_structural_feat", "any"),
        domain_names           = ("domain_names", "first"),
        any_junior_pass        = ("junior_pass", "any"),
    ).reset_index()

    return agg.sort_values(["any_junior_pass", "max_abs_delta_usage"], ascending=[False, False])


def _write_summary(hits: pd.DataFrame, gene_sum: pd.DataFrame,
                   elapsed: float, out_path: Path) -> None:
    passed  = hits[hits["junior_pass"]]
    dropped = hits[~hits["junior_pass"]]

    top10 = (
        passed.assign(_abs_du=passed["delta_usage"].abs())
              .drop_duplicates("gene_name")
              .nlargest(10, "_abs_du")[
                  ["gene_name", "transcript_name", "cell_type",
                   "delta_usage", "biotype_class", "proxy_type",
                   "ad_ot_score", "ad_ot_label", "ad_prior_category"]
              ]
    )

    lines = [
        "# ASSISTANT_SURVEYOR — Run Summary",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Runtime: {elapsed:.1f} s",
        "",
        "## Input",
        f"- Hits: {len(hits):,} rows ({hits['gene_name'].nunique():,} genes, "
        f"{hits['transcript_name'].nunique():,} transcripts)",
        "",
        "## Junior-Layer Gate (by biotype_class)",
        "| Outcome | Hits | Genes |",
        "|---------|------|-------|",
        f"| Pass | {len(passed)} | {passed['gene_name'].nunique()} |",
        f"| Drop | {len(dropped)} | {dropped['gene_name'].nunique()} |",
        "",
        "## biotype_class Distribution",
        "| biotype_class | Count | Junior Pass? |",
        "|---------------|-------|--------------|",
    ]
    for bc, cnt in hits["biotype_class"].value_counts().items():
        passes = "yes" if hits.loc[hits["biotype_class"] == bc, "junior_pass"].iloc[0] else "no"
        lines.append(f"| {bc} | {cnt} | {passes} |")

    lines += [
        "",
        "## proxy_type Distribution",
        "| proxy_type | Count |",
        "|------------|-------|",
    ]
    for pt, cnt in hits["proxy_type"].value_counts().items():
        lines.append(f"| {pt} | {cnt} |")

    lines += [
        "",
        "## OT Label Distribution",
        "| OT Label | Count |",
        "|----------|-------|",
    ]
    for lbl, cnt in hits["ad_ot_label"].value_counts().items():
        lines.append(f"| {lbl} | {cnt} |")

    lines += [
        "",
        "## AD Prior Hits",
        "| Category | Count |",
        "|----------|-------|",
    ]
    for cat, cnt in hits["ad_prior_category"].value_counts().items():
        lines.append(f"| {cat} | {cnt} |")

    lines += [
        "",
        "## Top 10 Junior-Pass Candidates (by |Δ usage|)",
        "| gene | transcript | cell_type | dPSI | biotype | proxy | OT score | OT label | AD prior |",
        "|------|-----------|-----------|------|---------|-------|----------|----------|----------|",
    ]
    for _, r in top10.iterrows():
        lines.append(
            f"| {r['gene_name']} | {r['transcript_name']} | {r['cell_type']} "
            f"| {r['delta_usage']:+.3f} | {r['biotype_class']} | {r['proxy_type']} "
            f"| {r['ad_ot_score']:.3f} | {r['ad_ot_label']} "
            f"| {r['ad_prior_category']} |"
        )

    lines += [
        "",
        "## Cell-Type Breakdown (Junior Pass)",
        "| Cell type | Pass hits |",
        "|-----------|-----------|",
    ]
    for ct, cnt in passed["cell_type"].value_counts().items():
        lines.append(f"| {ct} | {cnt} |")

    out_path.write_text("\n".join(lines) + "\n")


def main(no_cache: bool = False) -> None:
    t0 = time.time()

    if no_cache:
        _log("--no-cache: clearing existing cache files")
        for f in CACHE_DIR.glob("*.json"):
            f.unlink()
        for f in CACHE_DIR.glob("*.tsv"):
            f.unlink()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    _log("Loading hits CSV ...")
    hits_raw = pd.read_csv(HITS_CSV)
    _log(f"Loaded {len(hits_raw):,} hits")

    _log("=== L1: Biotype & CDS Classification ===")
    hits = l1_biotype.run(hits_raw)

    _log("=== L2: AD Gene Prior ===")
    hits = l2_ad_prior.run(hits)

    _log("=== L3: OpenTargets Batch Query ===")
    hits = l3_opentargets.run(hits)

    _log("=== L4: UniProt Domain Presence ===")
    hits = l4_uniprot.run(hits)

    _log("=== L5: Consequence Classification ===")
    hits = l5_consequences.run(hits)

    _log("=== Writing outputs ===")
    enriched_path = OUT_DIR / "hits_enriched.csv"
    hits.to_csv(enriched_path, index=False)
    _log(f"  hits_enriched.csv ({len(hits):,} rows, {len(hits.columns)} columns)")

    junior_pass_path = OUT_DIR / "junior_pass_candidates.csv"
    passed = hits[hits["junior_pass"]].copy()
    passed.to_csv(junior_pass_path, index=False)
    _log(f"  junior_pass_candidates.csv ({len(passed):,} rows)")

    gene_sum = _gene_summary(hits)
    gene_sum_path = OUT_DIR / "gene_summary.csv"
    gene_sum.to_csv(gene_sum_path, index=False)
    _log(f"  gene_summary.csv ({len(gene_sum):,} rows)")

    elapsed = time.time() - t0
    summary_path = OUT_DIR / "run_summary.md"
    _write_summary(hits, gene_sum, elapsed, summary_path)
    _log(f"  run_summary.md")

    _log(f"=== ASSISTANT_SURVEYOR complete in {elapsed:.1f} s ===")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ASSISTANT_SURVEYOR pipeline")
    parser.add_argument("--no-cache", action="store_true",
                        help="Ignore cached API responses and re-fetch")
    args = parser.parse_args()
    main(no_cache=args.no_cache)
