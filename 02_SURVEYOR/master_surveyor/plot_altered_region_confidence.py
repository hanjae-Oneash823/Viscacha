"""Altered-region structural confidence, intentionally outside the gate system.

Gate D only records whether canonical and alternative structures are present.
This module reuses those cached models to quantify how confidently a local
structural difference can be interpreted.  It neither changes nor contributes
to any gate.

The main view places model confidence for the altered span (mean pLDDT) against
the local displacement of that span after canonical/alternative superposition.
Large displacement is an observation, not a pass/fail rule; the plot makes it
possible to distinguish a strong structural signal from one in an uncertain
region of a model.
"""

from __future__ import annotations

from pathlib import Path
import sys
import gc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

from master_surveyor import m2_structures, m2b_structure_qc
from master_surveyor.config import OUT_DIR, STRUCTURE_CACHE_DIR
from master_surveyor.gate_matrix import build_gate_matrix


QC_CSV = OUT_DIR.parent / "altered_region_structure_confidence.csv"
PLOT_PATH = OUT_DIR / "10e_altered_region_structure_confidence.png"

VERDICT_COLORS = {"high": "#197a62", "medium": "#d17b28", "low": "#b63b50", "unknown": "#8d96a3"}
MAX_QC_SEQUENCE_LENGTH = 3000


def _canonical_path(sequence: str) -> Path | None:
    folder = STRUCTURE_CACHE_DIR / m2_structures.seq_hash(sequence)
    afdb = folder / "canonical_afdb.pdb"
    if afdb.exists():
        return afdb
    models = sorted((folder / "canonical_colabfold").glob("*_unrelaxed_rank_001_*.pdb"))
    return models[0] if models else None


def _cached_alt_result(sequence: str) -> m2_structures.AltStructureResult | None:
    folder = STRUCTURE_CACHE_DIR / m2_structures.seq_hash(sequence)
    models = sorted((folder / "alt_colabfold").glob("*_unrelaxed_rank_*.pdb"))
    if not models:
        return None
    scores = [m2_structures._scores_path_for(model) for model in models]
    esmfold = folder / "esmfold.pdb"
    return m2_structures.AltStructureResult(
        seed_models=models,
        scores_json=[path for path in scores if path is not None],
        esmfold_model=esmfold if esmfold.exists() else None,
    )


def build_qc_table() -> pd.DataFrame:
    """Assess every cached, Gate-E-ready pair without requesting new models."""
    pairs = build_gate_matrix()
    pairs = pairs[pairs["D"]].copy()
    records: list[dict] = []

    for _, pair in pairs.iterrows():
        base = {
            "gene_name": pair["gene_name"],
            "alt_transcript_name": pair["alt_transcript_name"],
            "cell_type": pair["cell_type"],
            "n_changed": pair["n_changed"],
            "changed_aa_start": pair["changed_aa_start"],
            "changed_aa_end": pair["changed_aa_end"],
        }
        # Parsing five full-atom ColabFold models plus the reference structure
        # for a multi-thousand-residue protein exceeds the available worker
        # memory.  Preserve the pair in the audit table, but do not pretend
        # it has a confidence score until it is handled in a larger-memory job.
        if max(len(pair["canonical_protein_seq"]), len(pair["alt_protein_seq"])) > MAX_QC_SEQUENCE_LENGTH:
            records.append({
                **base,
                "assessment_status": "not scored: model exceeds 3,000 aa QC limit",
            })
            continue
        canonical = _canonical_path(pair["canonical_protein_seq"])
        alt = _cached_alt_result(pair["alt_protein_seq"])
        if canonical is None or alt is None:
            records.append({**base, "assessment_status": "model files unavailable"})
            continue
        try:
            qc = m2b_structure_qc.assess(
                alt, canonical,
                pair["canonical_protein_seq"], pair["alt_protein_seq"],
                int(pair["changed_aa_start"]), int(pair["changed_aa_end"]),
            )
            region = qc.get("region_confidence", {})
            superposition = qc.get("superposition", {})
            records.append({
                **base,
                "assessment_status": "scored",
                "verdict": qc.get("verdict", "unknown"),
                "region_mean_plddt": region.get("region_mean_plddt"),
                "whole_protein_mean_plddt": region.get("whole_protein_mean_plddt"),
                "region_vs_rest_mean_pae": region.get("region_vs_rest_mean_pae"),
                "region_length": region.get("region_length"),
                "ensemble_spread_rmsd": qc.get("ensemble_spread_rmsd"),
                "colabfold_vs_esmfold_rmsd": qc.get("colabfold_vs_esmfold_rmsd"),
                "n_ensemble_models": qc.get("n_ensemble_models"),
                "changed_region_mean_local_dist": superposition.get("changed_region_mean_local_dist"),
                "changed_region_max_local_dist": superposition.get("changed_region_max_local_dist"),
                "outside_span_mean_local_dist": superposition.get("outside_span_mean_local_dist"),
                "anchor_fit_rmsd": superposition.get("anchor_fit_rmsd"),
                "n_anchor_residues": superposition.get("n_anchor_residues"),
                "error": qc.get("error") or superposition.get("error"),
            })
        except Exception as exc:  # one malformed model must not hide the rest
            records.append({**base, "assessment_status": "QC error", "error": str(exc)})
        finally:
            gc.collect()

    table = pd.DataFrame(records)
    QC_CSV.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(QC_CSV, index=False)
    return table


def plot(table: pd.DataFrame) -> None:
    scored = table[
        (table["assessment_status"] == "scored")
        & table["region_mean_plddt"].notna()
        & table["changed_region_mean_local_dist"].notna()
    ].copy()
    if scored.empty:
        raise RuntimeError("No pairs could be scored; see the QC table for details.")

    fig, ax = plt.subplots(figsize=(4.9, 4.35))
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.16, top=0.87)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axvspan(0, 50, color="#f5f6f8", zorder=0)
    ax.axvline(50, color="#a7adb5", linewidth=0.8, zorder=1)
    ax.axvline(70, color="#20242a", linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)

    sizes = 28 + 8 * np.sqrt(pd.to_numeric(scored["n_changed"], errors="coerce").fillna(1).clip(1, 100))
    for verdict in ["high", "medium", "low", "unknown"]:
        subset = scored[scored["verdict"].fillna("unknown") == verdict]
        if subset.empty:
            continue
        ax.scatter(
            subset["region_mean_plddt"], subset["changed_region_mean_local_dist"],
            s=sizes.loc[subset.index], c=VERDICT_COLORS[verdict], label=verdict.capitalize(),
            edgecolors="black", linewidths=0.55, alpha=0.9, zorder=3,
        )

    # Labels are reserved for the most interpretable candidates: confident
    # altered spans with the largest local displacement.  Labeling every pair
    # makes this plot unreadable and implies a false ranking for the rest.
    annotations = (
        scored[scored["region_mean_plddt"] >= 70]
        .sort_values("changed_region_mean_local_dist", ascending=False)
        .drop_duplicates("gene_name")
        .head(5)
    )
    texts = []
    for _, row in annotations.iterrows():
        texts.append(ax.text(
            row["region_mean_plddt"] + 1.1,
            row["changed_region_mean_local_dist"] + 1.2,
            row["gene_name"],
            fontsize=7.2, color="#15181c", zorder=4,
        ))

    # Use leader lines only when a label needs to move.  This retains a clean
    # compact panel while preventing the high-confidence candidates from
    # obscuring one another or their points.
    adjust_text(
        texts, ax=ax,
        expand=(1.18, 1.30), force_text=(0.28, 0.42), force_points=(0.22, 0.32),
        arrowprops={"arrowstyle": "-", "color": "#6b727c", "lw": 0.55},
    )

    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Altered-region mean pLDDT", fontsize=9.5)
    ax.set_ylabel("Mean altered-region displacement (Å)", fontsize=9.5)
    ax.set_title("Altered-region structural confidence", loc="left", fontsize=12.5, fontweight="bold", pad=8)
    ax.text(0.01, 1.01, f"{len(scored)} comparable pairs · not a gate", transform=ax.transAxes,
            fontsize=7.8, color="#5b6470", va="bottom")
    ax.text(4, ax.get_ylim()[1] * 0.96, "interpret cautiously", fontsize=8, color="#707780")
    ax.grid(axis="y", color="#e3e6ea", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("black")
    ax.tick_params(colors="#20242a", labelsize=8)
    legend = ax.legend(title="QC verdict", frameon=False, fontsize=7.5, title_fontsize=7.5,
                       loc="upper right", handletextpad=0.5)
    for handle in legend.legend_handles:
        handle.set_sizes([40])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    table = build_qc_table()
    plot(table)
    scored = (table["assessment_status"] == "scored").sum()
    print(f"Wrote {QC_CSV} ({scored}/{len(table)} scored)")
    print(f"Wrote {PLOT_PATH}")


if __name__ == "__main__":
    main()
