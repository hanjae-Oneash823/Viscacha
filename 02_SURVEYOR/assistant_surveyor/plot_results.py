"""
ASSISTANT_SURVEYOR visualization suite — 4 diagnostic plots.

Run from the repo root:
    python assistant_surveyor/plot_results.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

# ── paths ──────────────────────────────────────────────────────────────────
HITS_CSV = Path("outputs/assistant_surveyor/hits_enriched.csv")
OUT_DIR  = Path("outputs/assistant_surveyor/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── neon colormap (same as expression_vs_psi plots) ───────────────────────
NEON_CMAP = LinearSegmentedColormap.from_list(
    "neon", ["#0a0020", "#5522cc", "#cc1188", "#cc5500", "#bbbb00"]
)

def _neon(t: float) -> str:
    """Sample a hex color from the neon cmap at position t ∈ [0, 1]."""
    return mcolors.to_hex(NEON_CMAP(t))

# ── consistent colour palettes derived from the neon cmap ─────────────────
PROXY_COLORS = {
    "C":   "#cc1155",   # crimson-rose      ~330° — CDS + domain
    "D":   "#0099bb",   # teal-cyan         ~195° — CDS, no domain
    "NMD": "#cc8800",   # amber             ~ 40° — NMD isoform
    "N":   "#999999",   # neutral grey              — no protein change
}
GATE_COLORS  = {True: "#17a2a8", False: "#cccccc"}   # junior_pass: pass (turquoise) / drop
BIO_COLORS   = {
    "PC_CDS":    "#146b3a",   # dark green   — protein-coding, CDS-altering
    "PC_UTR":    "#3fa66d",   # medium green — protein-coding, UTR-only
    "PC_CDS_ND": "#8fd9b6",   # light green  — protein-coding, CDS not defined
    "novel":     _neon(0.75), # neon orange  — unchanged
    "RI":        "#4a1f7a",   # dark purple  — retained intron
    "NMD":       "#7c4aa8",   # medium purple — nonsense-mediated decay
    "TEC":       "#a97fc9",   # light purple — to be experimentally confirmed
    "other":     "#d4bce6",   # lightest purple — dropped, minor GENCODE biotypes
}
BIO_ORDER = ["PC_CDS", "PC_UTR", "PC_CDS_ND", "novel", "TEC", "RI", "NMD", "other"]
GROUP_COLORS = {
    "trial_failure_candidate": "#2a78d6",   # blue — dominant/canonical, down in AD
    "new_target_candidate":    "#e34948",   # red  — minor/alternate, up in AD
}

ZERO_COLOR   = "black"

CELL_ORDER = [
    "Excitatory_neuron", "Inhibitory_neuron", "Astrocyte",
    "Oligodendrocyte", "OPC", "Microglia", "Vascular_cell", "Lymphocyte",
]

BG   = "#ffffff"
FG   = "#1a1a1a"
GRID = "#e0e0e0"

plt.rcParams.update({
    "font.family":          "sans-serif",
    "font.size":            10,
    "text.color":           FG,
    "axes.facecolor":       BG,
    "figure.facecolor":     BG,
    "axes.edgecolor":       "#cccccc",
    "axes.labelcolor":      FG,
    "xtick.color":          FG,
    "ytick.color":          FG,
    "axes.spines.top":      False,
    "axes.spines.right":    False,
    "axes.spines.left":     True,
    "axes.spines.bottom":   True,
    "axes.grid":            False,
    "legend.facecolor":     BG,
    "legend.edgecolor":     "#cccccc",
    "legend.labelcolor":    FG,
    "figure.dpi":           150,
})


def save(fig: plt.Figure, name: str) -> None:
    path = OUT_DIR / name
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  saved → {path.name}")


def load() -> pd.DataFrame:
    df = pd.read_csv(HITS_CSV)
    df["abs_delta"]      = df["delta_usage"].abs()
    df["neg_log10_pval"] = -np.log10(df["permutation_pval"].clip(lower=1e-4))
    return df


def push_apart_y(ys: list[float], min_gap: float = 0.25,
                 ylim: tuple[float, float] = (0.0, 4.5)) -> list[float]:
    """Iteratively push label y-positions apart so annotations don't overlap."""
    ys = list(ys)
    for _ in range(300):
        moved = False
        order = sorted(range(len(ys)), key=lambda i: ys[i])
        for k in range(len(order) - 1):
            i, j = order[k], order[k + 1]
            gap = ys[j] - ys[i]
            if gap < min_gap:
                push   = (min_gap - gap) / 2
                ys[i] -= push
                ys[j] += push
                moved  = True
        if not moved:
            break
    margin = 0.05
    return [max(ylim[0] + margin, min(ylim[1] - margin, y)) for y in ys]


def plot_biotype_by_celltype(df: pd.DataFrame) -> None:
    """01 — Stacked bar: biotype_class per cell type."""
    bio_order = BIO_ORDER
    cells     = [c for c in CELL_ORDER if c in df["cell_type"].unique()]

    counts = (
        df.groupby(["cell_type", "biotype_class"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=cells, columns=bio_order, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    bottoms = np.zeros(len(cells))
    for bc in bio_order:
        vals = counts[bc].values
        ax.bar(cells, vals, bottom=bottoms, color=BIO_COLORS[bc],
               label=bc, edgecolor=BG, linewidth=0.5)
        bottoms += vals

    ax.set_xticklabels(cells, rotation=35, ha="right")
    ax.set_ylabel("Hits")
    ax.set_title("Biotype class distribution per cell type", fontweight="bold")
    ax.legend(title="biotype_class", bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, "01_biotype_by_celltype.png")


def plot_usage_direction(df: pd.DataFrame) -> None:
    """02 — Diverging bar: AD-enriched up, CT-enriched down, per cell type."""
    cells     = [c for c in CELL_ORDER if c in df["cell_type"].unique()]
    ad_counts = (df[df["usage_direction"] == "AD_enriched"]
                 .groupby("cell_type").size().reindex(cells, fill_value=0))
    ct_counts = (df[df["usage_direction"] == "CT_enriched"]
                 .groupby("cell_type").size().reindex(cells, fill_value=0))

    AD_COLOR = _neon(0.50)   # muted magenta
    CT_COLOR = _neon(0.22)   # muted violet

    x    = np.arange(len(cells))
    ymax = max(ad_counts.max(), ct_counts.max()) * 1.18

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.bar(x,  ad_counts.values,  color=AD_COLOR, edgecolor=BG, label="AD-enriched",  width=0.6)
    ax.bar(x, -ct_counts.values,  color=CT_COLOR, edgecolor=BG, label="CT-enriched",  width=0.6)

    ax.axhline(0, color=ZERO_COLOR, lw=0.8, alpha=0.5)

    # count labels
    for i, (ad, ct) in enumerate(zip(ad_counts.values, ct_counts.values)):
        if ad > 0:
            ax.text(i, ad + ymax * 0.02, str(ad), ha="center", va="bottom",
                    fontsize=8, color=FG)
        if ct > 0:
            ax.text(i, -(ct + ymax * 0.02), str(ct), ha="center", va="top",
                    fontsize=8, color=FG)

    # direction labels in chart area
    ax.text(len(cells) - 0.5, ymax * 0.88, "AD-enriched ↑",
            ha="right", color=AD_COLOR, fontsize=10, fontweight="bold")
    ax.text(len(cells) - 0.5, -ymax * 0.88, "↓ CT-enriched",
            ha="right", color=CT_COLOR, fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(cells, rotation=35, ha="right")
    ax.set_ylim(-ymax, ymax)

    # show absolute values on y axis
    yticks = [t for t in ax.get_yticks() if abs(t) <= ymax]
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(abs(int(t))) if t == int(t) else "" for t in yticks])
    ax.set_ylabel("Hits")

    ax.set_title("Isoform usage direction per cell type", fontweight="bold")
    save(fig, "02_usage_direction_by_celltype.png")


def plot_ot_specificity(df: pd.DataFrame) -> None:
    """03 — OT score vs AD specificity. Top-4 OT (navy) and top-4 specificity (teal) annotated."""
    SPEC_CAP   = 35.0
    ANNOT_OT   = "#003399"   # navy  — top OT score
    ANNOT_SPEC = "#007755"   # teal  — top AD specificity

    df_ot = df[df["ad_ot_score"] > 0].copy()
    df_ot["spec_clipped"] = df_ot["ad_specificity"].clip(upper=SPEC_CAP)

    passing  = df[df["junior_pass"]].drop_duplicates("gene_name")
    top_ot   = passing.nlargest(4, "ad_ot_score").copy()
    top_spec = passing[passing["ad_ot_score"] > 0].nlargest(4, "ad_specificity").copy()
    top_ot  ["spec_clipped"] = top_ot  ["ad_specificity"].clip(upper=SPEC_CAP)
    top_spec["spec_clipped"] = top_spec["ad_specificity"].clip(upper=SPEC_CAP)

    fig, ax = plt.subplots(figsize=(9, 7))

    # Layer 1 — N proxy_type: grey background
    sub_n = df_ot[df_ot["proxy_type"] == "N"]
    ax.scatter(
        sub_n["ad_ot_score"], sub_n["spec_clipped"],
        c="#cccccc", s=6, alpha=0.35, linewidths=0,
        rasterized=True, zorder=1,
    )

    # Layer 2 — NMD / D / C (not junior-pass): coloured, no edge
    for pt in ["NMD", "D", "C"]:
        sub = df_ot[(df_ot["proxy_type"] == pt) & (~df_ot["junior_pass"])]
        ax.scatter(
            sub["ad_ot_score"], sub["spec_clipped"],
            c=PROXY_COLORS[pt], s=14, alpha=0.75, linewidths=0,
            rasterized=True, zorder=2,
        )

    # Layer 3 — junior pass: coloured, larger, black edge
    sub_t1 = df_ot[df_ot["junior_pass"]]
    ax.scatter(
        sub_t1["ad_ot_score"], sub_t1["spec_clipped"],
        c=[PROXY_COLORS[pt] for pt in sub_t1["proxy_type"]],
        s=38, alpha=1.0, linewidths=0.7, edgecolors="black",
        rasterized=True, zorder=3,
    )

    # Reference line: specificity = 1
    ax.axhline(1.0, color=ZERO_COLOR, lw=0.6, ls="--", alpha=0.35)

    def _annotate_group(rows: pd.DataFrame, color: str,
                        min_gap: float, dx_sign: int) -> None:
        pts_y    = rows["spec_clipped"].tolist()
        label_ys = push_apart_y(pts_y, min_gap=min_gap, ylim=(0.0, SPEC_CAP))
        for (_, r), ly in zip(rows.iterrows(), label_ys):
            px = r["ad_ot_score"]
            dx = dx_sign * 0.04
            ax.annotate(
                r["gene_name"],
                xy=(px, r["spec_clipped"]), xytext=(px + dx, ly),
                fontsize=7.5, color=color,
                arrowprops=dict(arrowstyle="-", color=color, lw=0.7),
                ha="right" if dx < 0 else "left",
                va="center", zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor=color, linewidth=0.9, alpha=0.92),
            )

    # Top-4 OT score: labels go left (high x region)
    _annotate_group(top_ot,   ANNOT_OT,   min_gap=1.2, dx_sign=-1)
    # Top-4 specificity: labels go right (moderate x region)
    _annotate_group(top_spec, ANNOT_SPEC, min_gap=1.8, dx_sign=1)

    ax.set_xlim(-0.02, 0.60)
    ax.set_ylim(-0.8, SPEC_CAP * 1.05)
    ax.set_xlabel("AD OpenTargets score", fontsize=9)
    ax.set_ylabel("AD specificity  (AD / mean(PD, FTD, ALS))", fontsize=9)
    ax.set_title(
        "OT association vs disease specificity  |  colour = proxy type  |  black edge = junior pass",
        fontsize=10, pad=6,
    )
    ax.tick_params(labelsize=8)

    proxy_patches = [mpatches.Patch(color=PROXY_COLORS[pt], label=pt)
                     for pt in ["C", "D", "NMD", "N"]]
    pass_patch    = mpatches.Patch(facecolor="white", edgecolor="black",
                                   label="Junior pass", linewidth=0.8)
    ot_patch      = mpatches.Patch(facecolor="white", edgecolor=ANNOT_OT,
                                   label="top-4 OT score", linewidth=0.9)
    spec_patch    = mpatches.Patch(facecolor="white", edgecolor=ANNOT_SPEC,
                                   label="top-4 AD specificity", linewidth=0.9)
    ax.legend(handles=proxy_patches + [pass_patch, ot_patch, spec_patch],
              bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title_fontsize=8)

    save(fig, "03_ot_specificity.png")


# ═══════════════════════════════════════════════════════════════════════════
# SANKEY
# ═══════════════════════════════════════════════════════════════════════════

def plot_sankey(df: pd.DataFrame) -> None:
    """04 — Sankey: candidate group → biotype_class → proxy_type → junior gate."""
    import plotly.graph_objects as go

    def _rgba(hex_color: str, alpha: float = 0.45) -> str:
        r = int(hex_color[1:3], 16)
        g = int(hex_color[3:5], 16)
        b = int(hex_color[5:7], 16)
        return f"rgba({r},{g},{b},{alpha})"

    group_order = ["trial_failure_candidate", "new_target_candidate"]
    bio_order   = BIO_ORDER
    proxy_order = ["C", "D", "NMD", "N"]
    gate_order  = [True, False]
    gate_labels = {True: "Pass", False: "Drop"}

    n_grp = len(group_order)
    n_bio = len(bio_order)
    n_pxy = len(proxy_order)

    group_idx = {g: i                        for i, g in enumerate(group_order)}
    bio_idx   = {b: n_grp + i                for i, b in enumerate(bio_order)}
    proxy_idx = {p: n_grp + n_bio + i        for i, p in enumerate(proxy_order)}
    gate_idx  = {t: n_grp + n_bio + n_pxy + i for i, t in enumerate(gate_order)}

    group_labels = {
        "trial_failure_candidate": "Trial-failure candidate",
        "new_target_candidate":    "New-target candidate",
    }
    node_labels = (
        [group_labels[g] for g in group_order]
        + bio_order
        + proxy_order
        + [gate_labels[t] for t in gate_order]
    )
    node_colors = (
        [GROUP_COLORS[g]      for g in group_order]
        + [BIO_COLORS[b]        for b in bio_order]
        + [PROXY_COLORS[p]    for p in proxy_order]
        + [GATE_COLORS[t]     for t in gate_order]
    )

    sources, targets, values, link_colors = [], [], [], []

    # candidate group → biotype
    for grp in group_order:
        for bio in bio_order:
            n = int(((df["candidate_group"] == grp) & (df["biotype_class"] == bio)).sum())
            if n:
                sources.append(group_idx[grp])
                targets.append(bio_idx[bio])
                values.append(n)
                link_colors.append(_rgba(GROUP_COLORS[grp], 0.35))

    # biotype → proxy
    for bio in bio_order:
        for pxy in proxy_order:
            n = int(((df["biotype_class"] == bio) & (df["proxy_type"] == pxy)).sum())
            if n:
                sources.append(bio_idx[bio])
                targets.append(proxy_idx[pxy])
                values.append(n)
                link_colors.append(_rgba(BIO_COLORS[bio], 0.35))

    # proxy → junior gate
    for pxy in proxy_order:
        for t in gate_order:
            n = int(((df["proxy_type"] == pxy) & (df["junior_pass"] == t)).sum())
            if n:
                sources.append(proxy_idx[pxy])
                targets.append(gate_idx[t])
                values.append(n)
                link_colors.append(_rgba(PROXY_COLORS[pxy], 0.35))

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=18,
            thickness=22,
            line=dict(color="#cccccc", width=0.5),
            label=node_labels,
            color=node_colors,
            hovertemplate="%{label}: %{value}<extra></extra>",
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors,
        ),
    ))

    fig.update_layout(
        title=dict(
            text="Candidate group → Biotype class → Proxy type → Junior gate",
            font=dict(size=15, family="sans-serif", color=FG),
            x=0.5, xanchor="center",
        ),
        font=dict(family="sans-serif", size=12, color=FG),
        paper_bgcolor=BG,
        plot_bgcolor=BG,
        width=1100,
        height=600,
        margin=dict(l=20, r=20, t=60, b=20),
    )

    out = OUT_DIR / "04_sankey_group_biotype_proxy_gate.png"
    fig.write_image(str(out), scale=2)
    print(f"  saved → {out.name}")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main() -> None:
    print("Loading hits_enriched.csv ...")
    df = load()
    print(f"  {len(df):,} hits  |  {df['gene_name'].nunique():,} genes\n")

    print("Rendering plots ...")
    plot_biotype_by_celltype(df)
    plot_usage_direction(df)
    plot_ot_specificity(df)
    plot_sankey(df)

    n_saved = len(list(OUT_DIR.glob("*.png")))
    print(f"\nDone — {n_saved} plots in {OUT_DIR}")


if __name__ == "__main__":
    main()
