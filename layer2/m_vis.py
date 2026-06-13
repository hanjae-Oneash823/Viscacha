"""
M_VIS — Visualization Generation (plan 5, M_VIS).

Generates the inline SVGs embedded in the M08 HTML report, from M07's per-donor
PSI output. Vis 2 (gene expression heatmap) is intentionally skipped per the
build decision; Vis 1 and Vis 3 are produced per significant cell type.

  Vis 1 — per-donor PSI strip plot (focal AD-enriched isoform), points by
          condition with mean +/- SD overlay.
  Vis 3 — isoform proportion stacked bar (mean PSI per transcript per condition,
          stacked to 100% with an 'other' remainder).

Returns a dict {name: svg_string} consumed by M08.
"""

from __future__ import annotations

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from layer2.m07_expression import M07Result

COND_COLORS = {"Control": "#2ca02c", "AD": "#d62728", "Active control": "#ff7f0e"}
COND_ORDER = ["Control", "Active control", "AD"]
ISO_PALETTE = ["#1b7837", "#762a83", "#e08214", "#4575b4", "#d73027", "#01665e"]


def _fig_to_svg(fig) -> str:
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    svg = buf.getvalue()
    return svg[svg.index("<svg"):]      # drop XML/doctype header for inline embed


def _strip_plot(m07: M07Result) -> str | None:
    focal = next((t for t in m07.transcript_psi if t.role == "ad_enriched"),
                 m07.transcript_psi[0] if m07.transcript_psi else None)
    if not focal or not focal.donor_points:
        return None
    conds = [c for c in COND_ORDER if any(p[0] == c for p in focal.donor_points)]
    xpos = {c: i for i, c in enumerate(conds)}

    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    import numpy as np
    rng = np.random.default_rng(0)
    for c in conds:
        ys = [p[1] for p in focal.donor_points if p[0] == c]
        xs = xpos[c] + rng.uniform(-0.08, 0.08, len(ys))
        ax.scatter(xs, ys, color=COND_COLORS.get(c, "#777"), alpha=0.8, s=28,
                   edgecolor="white", linewidth=0.5, zorder=3)
        cp = next((x for x in focal.per_condition if x.condition == c), None)
        if cp and cp.n:
            ax.hlines(cp.mean, xpos[c] - 0.2, xpos[c] + 0.2, color="#222", lw=2, zorder=4)
            ax.vlines(xpos[c], cp.mean - cp.sd, cp.mean + cp.sd, color="#222", lw=1, zorder=4)

    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels(conds, rotation=15)
    ax.set_ylim(0, 1)
    ax.set_ylabel(f"PSI — {focal.transcript_name}")
    ax.set_title(f"Per-donor isoform usage\n{m07.gene_id} · {m07.cell_type.replace('_',' ')}",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    return _fig_to_svg(fig)


def _stacked_bar(m07: M07Result) -> str | None:
    if not m07.transcript_psi:
        return None
    conds, present = [], set()
    for t in m07.transcript_psi:
        for c in t.per_condition:
            present.add(c.condition)
    conds = [c for c in COND_ORDER if c in present]
    if not conds:
        return None

    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    bottoms = [0.0] * len(conds)
    for i, t in enumerate(m07.transcript_psi):
        means = []
        for c in conds:
            cp = next((x for x in t.per_condition if x.condition == c), None)
            means.append(cp.mean if cp and cp.mean == cp.mean else 0.0)
        ax.bar(range(len(conds)), means, bottom=bottoms,
               color=ISO_PALETTE[i % len(ISO_PALETTE)], label=t.transcript_name,
               width=0.6, edgecolor="white")
        bottoms = [b + m for b, m in zip(bottoms, means)]
    # 'other' remainder to 100%
    other = [max(0.0, 1.0 - b) for b in bottoms]
    if any(o > 0.01 for o in other):
        ax.bar(range(len(conds)), other, bottom=bottoms, color="#dddddd",
               label="other", width=0.6, edgecolor="white")

    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels(conds, rotation=15)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Mean PSI")
    ax.set_title(f"Isoform proportions\n{m07.gene_id} · {m07.cell_type.replace('_',' ')}",
                 fontsize=10)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    return _fig_to_svg(fig)


def run(m07: M07Result) -> dict[str, str]:
    svgs = {}
    s1 = _strip_plot(m07)
    if s1:
        svgs["psi_strip"] = s1
    s3 = _stacked_bar(m07)
    if s3:
        svgs["isoform_proportions"] = s3
    return svgs
