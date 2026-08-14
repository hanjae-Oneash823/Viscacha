"""Gate C -- ChEMBL and Open Targets clinical-phase matrix."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from master_surveyor.config import HITS_CSV
from master_surveyor.fda_approval import approved_names
from master_surveyor.gate_matrix import build_gate_matrix
from master_surveyor.plot_results import BG, FG, GRID, _style, save

CHEMBL = "#e5670a"
OT = "#5726b8"
NO_EVIDENCE = "#d9dde3"
MUTED = "#6b6b6b"


def _gene_phases(pairs: pd.DataFrame) -> pd.DataFrame:
    """One Gate-C evidence row per gene, including drug-record status counts."""
    d = pairs.copy()
    for col in ("chembl_max_phase", "ot_max_phase"):
        d[col] = pd.to_numeric(d[col], errors="coerce").fillna(0).astype(int)
    d = (d.groupby("gene_name", as_index=False)[["chembl_max_phase", "ot_max_phase"]]
           .max())
    d["max_phase"] = d[["chembl_max_phase", "ot_max_phase"]].max(axis=1)
    d["n_alt_proteins"] = pairs.groupby("gene_name").size().reindex(d["gene_name"]).values

    source = pd.read_csv(HITS_CSV, usecols=["gene_name", "drug_records_json"])
    source = source.dropna(subset=["drug_records_json"]).drop_duplicates("gene_name")

    direct_records: dict[str, list[dict]] = {}
    for _, row in source.iterrows():
        try:
            records = json.loads(row["drug_records_json"])
        except (TypeError, json.JSONDecodeError):
            records = []
        direct_records[row["gene_name"]] = [
            r for r in records if {"chembl", "ot"} & set(r.get("sources", []))
        ]

    # Phase four is a candidate approval signal, but FDA Drugs@FDA supplies
    # the independent US-approval verification used for the plotted count.
    fda_candidates = {
        r.get("name", "") for records in direct_records.values() for r in records
        if r.get("phase", 0) >= 4 and r.get("name")
    }
    fda_approved = approved_names(fda_candidates)

    def status_counts(gene: str) -> pd.Series:
        records = direct_records.get(gene, [])
        try:
            approved = sum(r.get("name") in fda_approved for r in records)
        except TypeError:
            approved = 0
        return pd.Series({
            "n_fda_approved": approved,
            "n_in_trials": sum(r.get("status") == "in_trials" for r in records),
            "n_failed_trial": sum(r.get("status") == "failed_trial" for r in records),
        })

    statuses = pd.DataFrame({gene: status_counts(gene) for gene in direct_records}).T
    d = d.join(statuses, on="gene_name").fillna({
        "n_fda_approved": 0, "n_in_trials": 0, "n_failed_trial": 0,
    })
    return d.sort_values(["max_phase", "ot_max_phase", "chembl_max_phase", "gene_name"],
                         ascending=[False, False, False, True]).reset_index(drop=True)


def plot_gate_c_evidence(pairs: pd.DataFrame) -> None:
    d = _gene_phases(pairs)
    d = d[d["max_phase"] >= 1].copy()
    y = np.arange(len(d))[::-1]

    fig = plt.figure(figsize=(6.9, 5.9), facecolor=BG)
    grid = fig.add_gridspec(
        1, 3, width_ratios=[0.72, 0.72, 0.62], wspace=0.13,
        left=0.16, right=0.99, top=0.84, bottom=0.17,
    )
    ax_chembl = fig.add_subplot(grid[0, 0])
    ax_ot = fig.add_subplot(grid[0, 1])
    ax_status = fig.add_subplot(grid[0, 2])
    axes = [ax_chembl, ax_ot, ax_status]
    for ax in axes:
        _style(ax)
        ax.set_ylim(-0.60, len(d) - 0.40)
        ax.set_yticks(y)
        ax.tick_params(axis="y", length=0)
        for side in ax.spines.values():
            side.set_color(FG)
            side.set_linewidth(0.9)

    for ax, col, title, color in (
        (ax_chembl, "chembl_max_phase", "ChEMBL phase", CHEMBL),
        (ax_ot, "ot_max_phase", "OT phase", OT),
    ):
        present = d[col] >= 1
        ax.scatter(d.loc[~present, col], y[~present], marker="o", s=24,
                   color=BG, edgecolor=FG, linewidth=0.8, zorder=2)
        for xi, yi in zip(d.loc[present, col], y[present]):
            ax.plot([0, xi], [yi, yi], color=color, alpha=0.65, linewidth=1.2, zorder=2,
                    solid_capstyle="round")
        ax.scatter(d.loc[present, col], y[present], marker="s" if col.startswith("chembl") else "o",
                   s=70, color=color, edgecolor=BG, linewidth=0.8, zorder=3)
        ax.set_xlim(-0.35, 4.45)
        ax.set_xticks([0, 1, 2, 3, 4])
        ax.set_title(title, fontsize=10, color=FG, pad=7)
        if ax is ax_chembl:
            ax.set_yticklabels(d["gene_name"], fontsize=8.8, color=FG)
            for label in ax.get_yticklabels():
                label.set_fontweight("bold")
        else:
            ax.set_yticklabels([])
        ax.xaxis.grid(True, color=GRID, lw=0.7)

    status_specs = [
        ("n_fda_approved", "#167a5b"),
        ("n_in_trials", "#246fa5"),
        ("n_failed_trial", "#b43a4e"),
    ]
    for yi, (_, row) in zip(y, d.iterrows()):
        for xpos, (col, color) in zip((-0.25, 0, 0.25), status_specs):
            value = int(row[col])
            ax_status.text(xpos, yi, str(value) if value else "—", ha="center", va="center",
                           fontsize=9, color=color if value else "#b7bcc4",
                           fontweight="bold" if value else "normal", zorder=3)
    ax_status.set_xlim(-0.42, 0.42)
    ax_status.set_xticks([])
    ax_status.set_title("Drug status", fontsize=10, color=FG, pad=7)
    ax_status.set_xlabel("FDA  ·  trial  ·  failed", fontsize=8.2, color=FG, labelpad=13)
    ax_status.set_yticklabels([])
    ax_status.xaxis.grid(False)

    fig.suptitle("Gate C — clinical drug evidence", fontsize=14, color=FG,
                 fontweight="bold", x=0.012, ha="left", y=0.99)
    fig.text(0.50, 0.075, "maximum clinical phase", fontsize=9, color=FG,
             ha="center", va="center")
    save(fig, "10c_gate_c_evidence_matrix.png")


def main() -> None:
    plot_gate_c_evidence(build_gate_matrix())


if __name__ == "__main__":
    main()
