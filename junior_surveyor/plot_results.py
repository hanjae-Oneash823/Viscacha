"""JUNIOR_SURVEYOR — visualization suite."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR   = REPO_ROOT / "outputs/junior_surveyor/plots"

# ---------------------------------------------------------------------------
# Design system
# ---------------------------------------------------------------------------
NEON_CMAP = LinearSegmentedColormap.from_list(
    "neon", ["#0a0020", "#5522cc", "#cc1188", "#cc5500", "#bbbb00"]
)

def _neon(t: float) -> str:
    return mcolors.to_hex(NEON_CMAP(t))

CHANGE_COLORS = {
    "frameshift_stop":    _neon(0.95),
    "N_truncation":       "#E5670A",   # burnt orange  — clearly distinct from C
    "C_truncation":       "#C41B3A",   # crimson       — clearly distinct from N
    "internal_indel":     _neon(0.50),
    "substitution":       _neon(0.42),
    "N_extension":        _neon(0.35),
    "internal_insertion": _neon(0.28),
    "C_extension":        _neon(0.20),
    "identical":          "#cccccc",
    "no_sequence":        "#eeeeee",
}

_CHANGE_ORDER = [
    "frameshift_stop", "N_truncation", "C_truncation",
    "internal_indel", "substitution",
    "N_extension", "internal_insertion", "C_extension",
    "identical", "no_sequence",
]

PROXY_COLORS = {"C": "#cc1155", "D": "#0099bb", "NMD": "#cc8800", "N": "#999999"}

# Mirrors assistant_surveyor's plot_sankey palette so the two Sankeys
# (biotype/proxy/junior-gate and its J4 extension here) read as one visual system.
BIO_COLORS = {
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
JUNIOR_GATE_COLORS = {True: "#17a2a8", False: "#cccccc"}   # junior_pass: pass (turquoise) / drop
GROUP_COLORS = {
    "trial_failure_candidate": "#2a78d6",   # blue — dominant/canonical, down in AD
    "new_target_candidate":    "#e34948",   # red  — minor/alternate, up in AD
    "other":                   "#c3c2b7",   # muted — dead end, not passed to assistant_surveyor
    "no_MANE_coverage":        "#898781",   # muted (darker) — dead end, no MANE Select entry
}
HITS_NODE_COLOR = "#1a1a1a"
GATE_COLORS = {
    "Both pass":   "#2e9e44",   # green      — selected_for_next_stage
    "Gate 1 only": "#7fb8e0",   # light blue — only one J4 sub-gate passed
    "Gate 2 only": "#7fb8e0",   # light blue — only one J4 sub-gate passed
    "Neither":     "#999999",   # grey       — dead end
}

BG         = "#ffffff"
FG         = "#1a1a1a"
GRID       = "#e0e0e0"
ANNOT_NAVY = "#003399"
ANNOT_TEAL = "#007755"

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _push_apart(ys: list[float], min_gap: float,
                lo: float, hi: float) -> list[float]:
    """Iteratively separate label positions to avoid overlap."""
    ys = list(ys)
    for _ in range(400):
        moved = False
        order = sorted(range(len(ys)), key=lambda i: ys[i])
        for k in range(len(order) - 1):
            i, j = order[k], order[k + 1]
            gap = ys[j] - ys[i]
            if gap < min_gap:
                push   = (min_gap - gap) / 2
                ys[i] -= push
                ys[j] += push
                moved   = True
        if not moved:
            break
    return [max(lo, min(hi, y)) for y in ys]


def _style(ax: plt.Axes) -> None:
    ax.set_facecolor(BG)
    ax.tick_params(colors=FG, labelsize=8)
    ax.xaxis.label.set_color(FG)
    ax.yaxis.label.set_color(FG)
    ax.title.set_color(FG)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID)
    ax.set_axisbelow(True)


def _grid(ax: plt.Axes) -> None:
    ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)


def save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    fig.patch.set_facecolor(BG)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"  saved → {path.name}")


def load() -> pd.DataFrame:
    csv = REPO_ROOT / "outputs/junior_surveyor/hits_deep.csv"
    df = pd.read_csv(csv)
    df["protein_change_type"] = df["protein_change_type"].fillna("no_sequence")
    df["affected_domain"]     = df["affected_domain"].fillna("none")
    df["is_druggable"]        = df["is_druggable"].fillna(False).astype(bool)
    df["premature_stop"]      = df["premature_stop"].fillna(False).astype(bool)
    df["chembl_max_phase"]    = df["chembl_max_phase"].fillna(0).astype(int)
    df["dgidb_interactions"]  = df["dgidb_interactions"].fillna(0).astype(int)
    df["abs_delta"]           = df["delta_usage"].abs()
    df["canon_len"]           = df["canonical_protein_seq"].fillna("").str.len()
    df["alt_len"]             = df["alt_protein_seq"].fillna("").str.len()
    return df


# ---------------------------------------------------------------------------
# Plot 01 — canonical vs alt protein length
# ---------------------------------------------------------------------------

def plot_protein_lengths(df: pd.DataFrame) -> None:
    """01 — Scatter: canonical protein length × alt protein length, colour = change type."""
    sub = df[
        (df["canon_len"] > 0) &
        (df["alt_len"] > 0) &
        (df["protein_change_type"] != "no_sequence")
    ].copy()

    lim    = max(sub["canon_len"].max(), sub["alt_len"].max()) * 1.35
    lo     = 40
    xs_ref = np.logspace(np.log10(lo), np.log10(lim), 300)

    fig, ax = plt.subplots(figsize=(8, 8))

    # ── Diagonal band: ±20% length change ────────────────────────────────────
    ax.fill_between(xs_ref, xs_ref * 0.80, xs_ref * 1.20,
                    color=GRID, alpha=0.18, linewidth=0, zorder=0)

    # Exact diagonal
    ax.plot(xs_ref, xs_ref, color="#aaaaaa", lw=0.8, ls="--", zorder=1)

    # ── Layer 1 — not selected ───────────────────────────────────────────────
    order = [ct for ct in reversed(_CHANGE_ORDER)
             if ct not in ("no_sequence", "frameshift_stop")]
    for ct in order:
        s = sub[(sub["protein_change_type"] == ct) & (~sub["selected_for_next_stage"])]
        if s.empty:
            continue
        ax.scatter(s["canon_len"], s["alt_len"],
                   c=CHANGE_COLORS[ct], s=11, alpha=0.40,
                   linewidths=0, rasterized=True, zorder=2)

    # ── Layer 2 — selected for next stage ────────────────────────────────────
    sel_all = sub[sub["selected_for_next_stage"]].copy()

    for ct in order:
        s = sel_all[sel_all["protein_change_type"] == ct]
        if s.empty:
            continue
        ax.scatter(s["canon_len"], s["alt_len"],
                   c=CHANGE_COLORS[ct], s=70,
                   alpha=1.0, linewidths=0.7, edgecolors=FG,
                   rasterized=True, zorder=3)

    # ── Annotations ──────────────────────────────────────────────────────────
    # Select by most extreme length change — these are furthest from diagonal
    # so arrows stay short and point to spatially prominent dots.
    sel_genes = sel_all.drop_duplicates("gene_name").copy()
    sel_genes["len_loss"] = sel_genes["canon_len"] - sel_genes["alt_len"]
    sel_genes["len_gain"] = sel_genes["alt_len"]   - sel_genes["canon_len"]

    _TRUNC_TYPES = {"N_truncation", "C_truncation", "internal_indel"}
    _EXT_TYPES   = {"N_extension", "C_extension", "internal_insertion"}

    top_trunc = (sel_genes[sel_genes["protein_change_type"].isin(_TRUNC_TYPES)]
                 .nlargest(5, "len_loss"))
    top_ext   = (sel_genes[sel_genes["protein_change_type"].isin(_EXT_TYPES)]
                 .nlargest(4, "len_gain"))

    # Use offset points (display coords) for arrows — avoids log-scale
    # coordinate transform bugs that break bbox_inches="tight".
    _TRUNC_OFFSETS = [(-65, -25), (-50, -42), (-80, -10), (-45, -55), (-70, -38)]
    _EXT_OFFSETS   = [(65, 25),   (50, 42),   (75, 10),   (80, 40)]

    def _ann(row: pd.Series, offsets: list, color: str, idx: int) -> None:
        dx, dy = offsets[idx % len(offsets)]
        ax.annotate(
            row["gene_name"],
            xy=(row["canon_len"], row["alt_len"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.5, color=color, va="center", ha="center",
            arrowprops=dict(arrowstyle="-", color="#888888",
                            lw=0.7, shrinkA=6, shrinkB=3),
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor=color, linewidth=0.9, alpha=0.92),
            zorder=6,
        )

    for i, (_, r) in enumerate(top_trunc.iterrows()):
        _ann(r, _TRUNC_OFFSETS, ANNOT_NAVY, i)
    for i, (_, r) in enumerate(top_ext.iterrows()):
        _ann(r, _EXT_OFFSETS, ANNOT_TEAL, i)

    # ── Axes ─────────────────────────────────────────────────────────────────
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, lim)
    ax.set_ylim(lo, lim)

    # Human-readable tick labels (100, 1000, etc.)
    from matplotlib.ticker import FuncFormatter
    fmt = FuncFormatter(lambda x, _: f"{int(x):,}")
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)

    _style(ax)
    _grid(ax)
    ax.set_xlabel("canonical protein length (aa)", fontsize=9)
    ax.set_ylabel("alt protein length (aa)", fontsize=9)
    ax.set_title(
        "Canonical vs alt isoform protein length  |  band = ±20%  |  edge = selected for next stage",
        fontsize=10, pad=8,
    )

    # ── Legend ────────────────────────────────────────────────────────────────
    ct_handles = [
        mpatches.Patch(color=CHANGE_COLORS[ct], label=ct)
        for ct in _CHANGE_ORDER
        if ct not in ("no_sequence", "frameshift_stop")
        and ct in sub["protein_change_type"].values
    ]
    tier1_h = mpatches.Patch(facecolor="white", edgecolor=FG,
                              linewidth=0.8, label="Deep Tier 1")
    trunc_h = mpatches.Patch(facecolor="white", edgecolor=ANNOT_NAVY,
                              linewidth=0.9, label="largest truncations")
    ext_h   = mpatches.Patch(facecolor="white", edgecolor=ANNOT_TEAL,
                              linewidth=0.9, label="largest extensions")
    ax.legend(
        handles=ct_handles + [tier1_h, trunc_h, ext_h],
        bbox_to_anchor=(1.02, 1), loc="upper left",
        fontsize=8, title="change type", title_fontsize=8,
        frameon=True, framealpha=1, edgecolor=GRID,
    )

    save(fig, "01_protein_lengths.png")


# ---------------------------------------------------------------------------
# Plot 02 — Pfam domain landscape (canonical transcripts)
# ---------------------------------------------------------------------------

_PFAM_CACHE = REPO_ROOT / "outputs/junior_surveyor/cache/j2_pfam_hits.json"

# Generic positional descriptors — not functional domain families
_PFAM_FILTER: set[str] = {
    "C-terminal", "C-terminal domain", "N-terminal", "N-terminal domain",
    "N-terminus", "C-terminus", "catalytic domain", "middle domain",
    "Lobe A", "Lobe B", "Lobe C",
}

# Collapse redundant Pfam entries to a single canonical family name.
# Keys are the EXACT strings stored in j2_pfam_hits.json (hit.name from pyhmmer).
# Some Pfam names legitimately contain commas — they appear here as single entries.
_PFAM_NORM: dict[str, str] = {
    # Zinc finger C2H2 — multiple Pfam entries describe the same structural fold
    "Zinc finger, C2H2 type":                              "Zinc finger (C2H2)",
    "C2H2-type zinc finger":                               "Zinc finger (C2H2)",
    "Zinc-finger of C2H2 type":                            "Zinc finger (C2H2)",
    "Zinc finger, C2H2 type, C2H2-type zinc finger":       "Zinc finger (C2H2)",
    "Zinc finger C-x8-C-x5-C-x3-H type (and similar)":   "Zinc finger (C2H2)",
    "Zinc finger protein 462, seventh C2H2 zinc finger":   "Zinc finger (C2H2)",
    "Zinc finger protein 462, C2H2 zinc finger":           "Zinc finger (C2H2)",
    # Protein kinase
    "Protein kinase domain":                               "Protein kinase",
    "Protein tyrosine and serine/threonine kinase":        "Protein kinase",
    "ABC1 atypical kinase-like domain":                    "Protein kinase",
    "Kinase-like":                                         "Protein kinase",
    "Protein kinase C terminal domain":                    "Protein kinase",
    # WD40 — pyhmmer stores as single compound entry
    "WD domain, G-beta repeat":                            "WD40 repeat",
    # Immunoglobulin
    "Immunoglobulin domain":                               "Immunoglobulin-like",
    "Immunoglobulin I-set domain":                         "Immunoglobulin-like",
    "Immunoglobulin V-set domain":                         "Immunoglobulin-like",
    "immunoglobulin-like domain":                          "Immunoglobulin-like",
    # SH3
    "Variant SH3 domain":                                  "SH3 domain",
    # Ankyrin
    "Ankyrin repeats (3 copies)":                          "Ankyrin repeat",
    "Ankyrin repeats (many copies)":                       "Ankyrin repeat",
    # EGF-like
    "Calcium-binding EGF domain":                          "EGF-like domain",
    "Human growth factor-like EGF":                        "EGF-like domain",
    "Laminin EGF domain":                                  "EGF-like domain",
    "EGF domain":                                          "EGF-like domain",
    # Leucine-rich repeat
    "Leucine rich repeat":                                 "Leucine-rich repeat",
    "Leucine Rich Repeat":                                 "Leucine-rich repeat",
    "Leucine Rich repeats (2 copies)":                     "Leucine-rich repeat",
    "Leucine Rich repeat":                                 "Leucine-rich repeat",
    # RING finger
    "Zinc finger, C3HC4 type (RING finger)":               "RING finger",
    "Ring finger domain":                                  "RING finger",
    "RING-type zinc-finger":                               "RING finger",
    "RING-like zinc finger":                               "RING finger",
    # PH domain
    "Pleckstrin homology domain":                          "PH domain",
    # HEAT / TPR
    "HEAT repeats":                                        "HEAT repeat",
    "Tetratricopeptide repeat":                            "TPR repeat",
    "TPR repeats":                                         "TPR repeat",
    "TPR domain":                                          "TPR repeat",
    # Helicase
    "Helicase conserved C-terminal domain":                "Helicase",
    "DEAD/DEAH box helicase":                              "Helicase",
    "Helicase associated domain (HA2)":                    "Helicase",
    # AAA ATPase
    "AAA+ ATPase lid domain":                              "AAA ATPase",
    "AAA domain (dynein-related subfamily)":               "AAA ATPase",
    "AAA domain":                                          "AAA ATPase",
    # Cadherin
    "Cadherin-like":                                       "Cadherin domain",
    "Cadherin C-terminal cytoplasmic tail":                "Cadherin domain",
    "catenin-binding region":                              "Cadherin domain",
    # PAS domain
    "PAS fold":                                            "PAS domain",
    # P-type ATPase
    "P-type ATPase, cytoplasmic domain N":                 "P-type ATPase",
    "P-type ATPase actuator domain":                       "P-type ATPase",
    "Cation transporter/ATPase":                           "P-type ATPase",
    "Cation transporter/ATPase, N-terminus":               "P-type ATPase",
    "Cation transporting ATPase":                          "P-type ATPase",
    # Roc/COR (LRRK/DAPK family — stored as compound entry)
    "Ras of Complex, Roc, domain of DAPkinase":            "Roc/COR domain",
    # DENN
    "uDENN domain":                                        "DENN domain",
    "dDENN domain":                                        "DENN domain",
    "DENN (AEX-3) domain":                                 "DENN domain",
    # C1 domain
    "Phorbol esters/diacylglycerol binding domain (C1 domain)": "C1 domain",
    # PHD finger
    "PHD-zinc-finger like domain":                         "PHD finger",
    "PHD-like zinc-binding domain":                        "PHD finger",
    # FERM
    "FERM central domain":                                 "FERM domain",
    "FERM N-terminal domain":                              "FERM domain",
    # Ubiquitin
    "Ubiquitin family":                                    "Ubiquitin-like",
    "ubiquitin-like domain":                               "Ubiquitin-like",
    "Ubiquitin carboxyl-terminal hydrolase":               "Ubiquitin-like",
    # Fibronectin
    "Fibronectin type III domain":                         "Fibronectin type III",
    # Reprolysin metalloprotease
    "Reprolysin (M12B) family zinc metalloprotease":       "Reprolysin (M12B)",
    "Metallo-peptidase family M12B Reprolysin-like":       "Reprolysin (M12B)",
    # F-box
    "F-box-like":                                          "F-box domain",
    # RRM — Pfam full name contains commas; stored as single entry in cache
    "RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)": "RNA recognition motif",
    # TOG domain
    "XMAP215/Dis1/CLASP, TOG domain":                      "TOG domain",
    # LDL receptor repeats
    "Low-density lipoprotein receptor domain class A":     "LDL receptor repeat",
    "Low-density lipoprotein receptor repeat class B":     "LDL receptor repeat",
    # Type III restriction enzyme — spurious annotation in human proteins, exclude
    "Type III restriction enzyme, res subunit":            None,
}


def _norm_domain(name: str) -> str | None:
    if name in _PFAM_FILTER:
        return None
    mapped = _PFAM_NORM.get(name, name)
    return mapped  # may be None for explicitly excluded entries


def plot_pfam_domains(df: pd.DataFrame) -> None:
    """02 — Pfam domain landscape of canonical transcripts (pyhmmer source)."""
    with open(_PFAM_CACHE) as fh:
        cache: dict = json.load(fh)

    genes = df.drop_duplicates("canonical_enst")[
        ["canonical_enst", "pfam_n_domains"]
    ].copy()

    # ── n_domains distribution ────────────────────────────────────────────────
    ndist = genes["pfam_n_domains"].value_counts().sort_index()

    # ── Domain family counts (one count per gene, not per hit) ────────────────
    domain_counts: Counter = Counter()
    for enst in genes["canonical_enst"].dropna():
        key   = f"{enst}__canonical"
        normed = {_norm_domain(h["name"]) for h in cache.get(key, [])}
        normed.discard(None)
        domain_counts.update(normed)

    top = pd.Series(dict(domain_counts)).sort_values(ascending=True).tail(25)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1,
        figsize=(9, 11),
        gridspec_kw={"height_ratios": [1.6, 5], "hspace": 0.38},
    )

    # ── Panel A: n_domains distribution ──────────────────────────────────────
    xs = ndist.index.tolist()
    ys = ndist.values.tolist()
    ax_top.bar(xs, ys, color=ANNOT_TEAL, alpha=0.82, width=0.65, zorder=2)
    for x, y in zip(xs, ys):
        if y > 0:
            ax_top.text(x, y + 1.5, str(y), ha="center", va="bottom",
                        fontsize=6.5, color=FG)
    ax_top.set_xlabel("Pfam domains per canonical transcript", fontsize=8)
    ax_top.set_ylabel("genes", fontsize=8)
    ax_top.set_title(
        "A  —  Pfam domain count per canonical transcript  "
        f"({int((genes['pfam_n_domains'] > 0).sum()):,} / {len(genes):,} genes have ≥1 domain)",
        fontsize=9, pad=7, loc="left",
    )
    ax_top.set_xticks(range(0, max(xs) + 1))
    _style(ax_top)
    _grid(ax_top)
    ax_top.set_ylim(0, max(ys) * 1.15)

    # ── Panel B: top 25 domain families (lollipop) ───────────────────────────
    ypos = range(len(top))
    ax_bot.hlines(ypos, 0, top.values, color=GRID, lw=1.4, zorder=1)
    ax_bot.scatter(top.values, list(ypos),
                   color=ANNOT_TEAL, s=60, zorder=3, linewidths=0)
    for val, y in zip(top.values, ypos):
        ax_bot.text(val + 0.6, y, str(val),
                    va="center", fontsize=7, color=FG)
    ax_bot.set_yticks(list(ypos))
    ax_bot.set_yticklabels(top.index.tolist(), fontsize=7.5)
    ax_bot.set_xlabel("number of canonical transcripts", fontsize=8)
    ax_bot.set_title(
        "B  —  Top 25 Pfam domain families detected by pyhmmer\n"
        "      (counted per gene; redundant sub-entries collapsed; positional terms excluded)",
        fontsize=9, pad=7, loc="left",
    )
    ax_bot.set_xlim(0, top.max() * 1.14)
    ax_bot.yaxis.grid(False)
    ax_bot.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    _style(ax_bot)
    # remove left spine for cleaner look
    ax_bot.spines["left"].set_visible(False)
    ax_bot.tick_params(axis="y", length=0)

    save(fig, "02_pfam_domains.png")


# ---------------------------------------------------------------------------
# Plot 03 — Heatmap: affected domain family × protein change type
# ---------------------------------------------------------------------------

_CT_ORDER = [
    "N_truncation", "C_truncation", "internal_indel",
    "N_extension", "C_extension", "internal_insertion", "substitution",
]

_CT_LABELS = {
    "N_truncation":      "N truncation",
    "C_truncation":      "C truncation",
    "internal_indel":    "internal indel",
    "N_extension":       "N extension",
    "C_extension":       "C extension",
    "internal_insertion":"internal insertion",
    "substitution":      "substitution",
}


def _affected_domain_events(df: pd.DataFrame, cache: dict) -> pd.DataFrame:
    """Return long-form table of (change_type, normalized_domain_family) events."""
    rows = []
    skip = {"no_sequence", "identical", "frameshift_stop"}
    for _, row in df.iterrows():
        ct    = row.get("protein_change_type", "")
        start = row.get("changed_aa_start", 0)
        end   = row.get("changed_aa_end",   0)
        if ct in skip or not start:
            continue
        key = f"{row.get('canonical_enst', '')}__canonical"
        for h in cache.get(key, []):
            if h["start"] <= end and h["end"] >= start:
                normed = _norm_domain(h["name"])
                if normed:
                    rows.append({"change_type": ct, "domain": normed})
    return pd.DataFrame(rows)


def plot_affected_domain_heatmap(df: pd.DataFrame) -> None:
    """03 — Heatmap: which domain families are disrupted by each change type."""
    with open(_PFAM_CACHE) as fh:
        cache: dict = json.load(fh)

    events = _affected_domain_events(df, cache)

    # Count matrix: domain (rows) × change_type (cols)
    mat = (
        events.groupby(["domain", "change_type"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=[ct for ct in _CT_ORDER if ct in events["change_type"].unique()],
                 fill_value=0)
    )
    # Keep top 20 domains by total events
    mat["_total"] = mat.sum(axis=1)
    mat = mat.nlargest(20, "_total").drop(columns="_total")
    mat = mat.sort_values(mat.columns.tolist(), ascending=True)  # sort by first col for readability
    # Sort by total descending (top domain at top of plot)
    mat = mat.iloc[::-1]

    n_rows, n_cols = mat.shape

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 8))

    # Square-root colorscale so large truncation counts don't wash out others
    data   = mat.values.astype(float)
    data_s = np.sqrt(data)
    cmap   = LinearSegmentedColormap.from_list("teal_seq", ["#f0f9f5", ANNOT_TEAL])
    im     = ax.imshow(data_s, aspect="auto", cmap=cmap,
                       vmin=0, vmax=data_s.max())

    # Cell annotations: show raw count; blank for zero
    for r in range(n_rows):
        for c in range(n_cols):
            val = int(data[r, c])
            if val == 0:
                continue
            brightness = data_s[r, c] / data_s.max()
            txt_color  = "white" if brightness > 0.55 else FG
            ax.text(c, r, str(val), ha="center", va="center",
                    fontsize=7.5, color=txt_color, fontweight="normal")

    # Axes
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(
        [_CT_LABELS.get(c, c) for c in mat.columns],
        fontsize=8, rotation=35, ha="right",
    )
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(mat.index.tolist(), fontsize=7.5)
    ax.set_title(
        "Pfam domain families disrupted by protein change type\n"
        "(count of transcripts per cell; colour scale = √count)",
        fontsize=9, pad=10,
    )

    # Colorbar
    cb = fig.colorbar(im, ax=ax, shrink=0.55, pad=0.02)
    cb.set_label("√(transcript count)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    # Grid lines between cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    _style(ax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    save(fig, "03_affected_domain_heatmap.png")


# ---------------------------------------------------------------------------
# Plots 04 / 05 — Domain loss vs gain (diverging bar)
# ---------------------------------------------------------------------------

def _domain_loss_gain_events(df: pd.DataFrame, cache: dict) -> pd.DataFrame:
    """Long-form table of domain loss/gain events with normalized domain names."""
    skip = {"no_sequence", "frameshift_stop"}
    rows = []
    for _, row in df.iterrows():
        ct = row.get("protein_change_type", "")
        if ct in skip:
            continue
        canon_key = f"{row.get('canonical_enst', '')}__canonical"
        alt_key   = f"{row.get('ENST_ID', '')}__alt"
        canon_acc = {h["acc"]: h["name"] for h in cache.get(canon_key, [])}
        alt_acc   = {h["acc"]: h["name"] for h in cache.get(alt_key,   [])}
        for acc, name in canon_acc.items():
            if acc not in alt_acc:
                normed = _norm_domain(name)
                if normed:
                    rows.append({"side": "loss", "change_type": ct, "domain": normed})
        for acc, name in alt_acc.items():
            if acc not in canon_acc:
                normed = _norm_domain(name)
                if normed:
                    rows.append({"side": "gain", "change_type": ct, "domain": normed})
    return pd.DataFrame(rows)


def _build_stacked_matrix(
    events: pd.DataFrame,
    side: str,
    domains: list[str],
    group_map: dict[str, str],
    group_order: list[str],
) -> pd.DataFrame:
    """Return (domains × groups) count matrix for one side of the diverging bar."""
    sub = events[events["side"] == side].copy()
    sub["group"] = sub["change_type"].map(lambda x: group_map.get(x, group_order[-1]))
    mat = (
        sub.groupby(["domain", "group"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=domains, fill_value=0)
    )
    for g in group_order:
        if g not in mat.columns:
            mat[g] = 0
    return mat[group_order]


def _draw_diverging(
    ax: plt.Axes,
    domains: list[str],
    mat_loss: pd.DataFrame,
    mat_gain: pd.DataFrame,
    loss_colors: list[str],
    gain_colors: list[str],
    loss_labels: list[str],
    gain_labels: list[str],
    title: str,
) -> None:
    """Core diverging stacked bar renderer."""
    n    = len(domains)
    ypos = np.arange(n)

    # Loss bars (extend left, negative x)
    left_loss = np.zeros(n)
    loss_handles = []
    for col, color, label in zip(mat_loss.columns, loss_colors, loss_labels):
        vals = mat_loss[col].values.astype(float)
        bars = ax.barh(ypos, -vals, left=-left_loss,
                       color=color, height=0.62, zorder=2)
        loss_handles.append(mpatches.Patch(color=color, label=label))
        left_loss += vals

    # Gain bars (extend right, positive x)
    left_gain = np.zeros(n)
    gain_handles = []
    for col, color, label in zip(mat_gain.columns, gain_colors, gain_labels):
        vals = mat_gain[col].values.astype(float)
        ax.barh(ypos, vals, left=left_gain,
                color=color, height=0.62, zorder=2, alpha=0.88)
        gain_handles.append(mpatches.Patch(color=color, label=label))
        left_gain += vals

    # Total count annotations
    total_loss = mat_loss.sum(axis=1).values
    total_gain = mat_gain.sum(axis=1).values
    x_max = max(total_loss.max(), 1)

    for i, (tl, tg) in enumerate(zip(total_loss, total_gain)):
        if tl > 0:
            ax.text(-tl - x_max * 0.012, i, str(int(tl)),
                    ha="right", va="center", fontsize=6.5, color=FG)
        if tg > 0:
            ax.text(tg + x_max * 0.012, i, str(int(tg)),
                    ha="left",  va="center", fontsize=6.5, color=FG)

    # Centre line
    ax.axvline(0, color=FG, lw=0.8, zorder=3)

    # Axes styling
    ax.set_yticks(ypos)
    ax.set_yticklabels(domains, fontsize=7.5)
    ax.set_xlabel("← domain loss                domain gain →", fontsize=8)
    ax.set_title(title, fontsize=9, pad=8, loc="left")

    # Symmetric-ish x-axis; let loss dominate naturally
    x_edge = total_loss.max() * 1.18
    ax.set_xlim(-x_edge, max(total_gain.max() * 3.5, x_edge * 0.15))
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: str(int(abs(x))))
    )

    _style(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(False)

    # Legend: loss (left column) | gain (right column)
    leg_loss = ax.legend(handles=loss_handles, loc="lower left",
                         bbox_to_anchor=(0.01, -0.22), ncol=len(loss_handles),
                         fontsize=7, title="loss by:", title_fontsize=7,
                         frameon=True, framealpha=1, edgecolor=GRID)
    ax.add_artist(leg_loss)
    ax.legend(handles=gain_handles, loc="lower right",
              bbox_to_anchor=(0.99, -0.22), ncol=len(gain_handles),
              fontsize=7, title="gain by:", title_fontsize=7,
              frameon=True, framealpha=1, edgecolor=GRID)


# Shared colour palette for change types
_LOSS_COLORS_A = [CHANGE_COLORS["N_truncation"],
                  CHANGE_COLORS["C_truncation"],
                  "#cccccc"]
_LOSS_LABELS_A = ["N truncation", "C truncation", "other"]

_GAIN_COLORS_A = [ANNOT_TEAL]
_GAIN_LABELS_A = ["all gains"]

_LOSS_COLORS_C = [CHANGE_COLORS["N_truncation"],
                  CHANGE_COLORS["C_truncation"],
                  "#cccccc"]
_LOSS_LABELS_C = ["N truncation", "C truncation", "other"]

_GAIN_COLORS_C = [CHANGE_COLORS["C_extension"],
                  CHANGE_COLORS["N_extension"],
                  "#cccccc"]
_GAIN_LABELS_C = ["C extension", "N extension", "other"]


def plot_domain_loss_gain(df: pd.DataFrame) -> None:
    """04 — Diverging bar: domain loss (stacked N/C-trunc) vs gain (single colour)."""
    with open(_PFAM_CACHE) as fh:
        cache: dict = json.load(fh)

    events  = _domain_loss_gain_events(df, cache)
    domains = (events[events["side"] == "loss"]["domain"]
               .value_counts().head(15).index.tolist())
    domains = domains[::-1]  # highest at top

    loss_map = {"N_truncation": "N_truncation", "C_truncation": "C_truncation"}
    gain_map: dict[str, str] = {}

    mat_loss = _build_stacked_matrix(
        events, "loss", domains, loss_map,
        ["N_truncation", "C_truncation", "other"])
    mat_gain = _build_stacked_matrix(
        events, "gain", domains, gain_map, ["other"])

    fig, ax = plt.subplots(figsize=(9, 7))
    _draw_diverging(
        ax, domains, mat_loss, mat_gain,
        _LOSS_COLORS_A, _GAIN_COLORS_A,
        _LOSS_LABELS_A, _GAIN_LABELS_A,
        "A  —  Domain loss vs gain per Pfam family\n"
        "      (loss stacked by truncation direction; gain = all types combined)",
    )
    fig.subplots_adjust(bottom=0.18)
    save(fig, "04_domain_loss_gain.png")


def plot_domain_loss_gain_detail(df: pd.DataFrame) -> None:
    """05 — Diverging bar: full change-type breakdown on both loss and gain sides."""
    with open(_PFAM_CACHE) as fh:
        cache: dict = json.load(fh)

    events  = _domain_loss_gain_events(df, cache)
    domains = (events[events["side"] == "loss"]["domain"]
               .value_counts().head(15).index.tolist())
    domains = domains[::-1]

    loss_map = {"N_truncation": "N_truncation", "C_truncation": "C_truncation"}
    gain_map = {"C_extension": "C_extension",   "N_extension":  "N_extension"}

    mat_loss = _build_stacked_matrix(
        events, "loss", domains, loss_map,
        ["N_truncation", "C_truncation", "other"])
    mat_gain = _build_stacked_matrix(
        events, "gain", domains, gain_map,
        ["C_extension", "N_extension", "other"])

    fig, ax = plt.subplots(figsize=(9, 7))
    _draw_diverging(
        ax, domains, mat_loss, mat_gain,
        _LOSS_COLORS_C, _GAIN_COLORS_C,
        _LOSS_LABELS_C, _GAIN_LABELS_C,
        "B  —  Domain loss vs gain — full change-type breakdown\n"
        "      (gain side: C extension / N extension / other)",
    )
    fig.subplots_adjust(bottom=0.18)
    save(fig, "05_domain_loss_gain_detail.png")


# ---------------------------------------------------------------------------
# Plots 06 / 07 / 08 — J3 drug-target landscape
# ---------------------------------------------------------------------------

_PHASE_COLORS = {
    0: "#cccccc",   # no clinical data
    1: "#8ecae6",   # phase 1
    2: "#219ebc",   # phase 2
    3: "#fb8500",   # phase 3
    4: "#e63946",   # approved
}
_PHASE_LABELS = {
    0: "Phase 0 / preclinical",
    1: "Phase 1",
    2: "Phase 2",
    3: "Phase 3",
    4: "Phase 4 (approved)",
}

# Druggability tier assignment per gene row
def _drug_tier(row: pd.Series) -> str:
    phase = int(row.get("chembl_max_phase", 0) or 0)
    if phase == 4:
        return "phase4"
    if phase in (1, 2, 3):
        return "phase13"
    if row.get("is_druggable") or int(row.get("dgidb_interactions", 0) or 0) > 0:
        return "evidence"
    return "none"


def plot_druggability_funnel(df: pd.DataFrame) -> None:
    """06 — Horizontal funnel: how many genes pass each druggability threshold."""
    genes = df.drop_duplicates("gene_name")
    n_total = len(genes)

    steps = [
        ("All cohort genes",             n_total),
        ("ChEMBL target entry",          int(genes["is_druggable"].sum())),
        ("DGIdb interactions (any)",     int((genes["dgidb_interactions"] > 0).sum())),
        ("Clinical trial (Phase ≥ 1)",   int((genes["chembl_max_phase"] >= 1).sum())),
        ("Approved drug (Phase 4)",      int((genes["chembl_max_phase"] == 4).sum())),
    ]
    labels = [s[0] for s in steps]
    counts = [s[1] for s in steps]
    colors = ["#d0d0d0", "#a8dadc", ANNOT_TEAL, "#fb8500", "#e63946"]

    fig, ax = plt.subplots(figsize=(8, 5))

    ypos = np.arange(len(steps))
    max_c = counts[0]

    for i, (label, count, color) in enumerate(zip(labels, counts, colors)):
        half = count / 2
        ax.barh(i, count, left=-half, height=0.58, color=color, zorder=2)
        pct = count / max_c * 100
        # Count left of bar, percentage right
        ax.text(-half - max_c * 0.01, i, f"{count:,}",
                ha="right", va="center", fontsize=8, color=FG, fontweight="bold")
        ax.text(half + max_c * 0.01, i, f"{pct:.0f}%",
                ha="left",  va="center", fontsize=8, color=FG)

    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlim(-max_c * 0.68, max_c * 0.68)
    ax.set_xlabel("Gene count (centered)", fontsize=8)
    ax.set_title("Druggability funnel — Junior Surveyor cohort\n"
                 "Each bar is an independent criterion (not a strict subset)", fontsize=9, loc="left")
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: str(int(abs(x)))))
    ax.axvline(0, color=FG, lw=0.6, zorder=3)

    _style(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(False)

    fig.tight_layout()
    save(fig, "06_druggability_funnel.png")


def plot_drug_scatter(df: pd.DataFrame) -> None:
    """07 — Scatter: ChEMBL n_drugs vs DGIdb interactions, coloured by max clinical phase."""
    genes = df.drop_duplicates("gene_name").copy()
    # Include all genes with any drug evidence
    sub = genes[(genes["is_druggable"]) | (genes["dgidb_interactions"] > 0)].copy()

    rng = np.random.default_rng(42)
    jitter_x = rng.uniform(-0.35, 0.35, len(sub))
    jitter_y = rng.uniform(-1.2,   1.2,  len(sub))

    x = sub["n_drugs"].values.astype(float)          + jitter_x
    y = sub["dgidb_interactions"].values.astype(float) + jitter_y
    phases = sub["chembl_max_phase"].fillna(0).astype(int).values
    colors_pts = [_PHASE_COLORS.get(p, "#cccccc") for p in phases]
    sizes = [15 if p == 0 else 35 for p in phases]

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(x, y, c=colors_pts, s=sizes, alpha=0.75, linewidths=0.3,
               edgecolors="#555555", zorder=2)

    # Annotate standouts
    ann_genes = sub[(sub["n_drugs"] >= 4) | (sub["dgidb_interactions"] >= 50)]
    _offsets_scatter = [
        (22, 8), (-22, 8), (22, -10), (-22, -10), (15, 15),
        (-30, 12), (25, -15), (-15, 20), (30, 5),
    ]
    for i, (_, row) in enumerate(ann_genes.iterrows()):
        gx = row["n_drugs"] + jitter_x[sub.index.get_loc(row.name)]
        gy = row["dgidb_interactions"] + jitter_y[sub.index.get_loc(row.name)]
        dx, dy = _offsets_scatter[i % len(_offsets_scatter)]
        ax.annotate(
            row["gene_name"],
            xy=(gx, gy), xytext=(dx, dy), textcoords="offset points",
            fontsize=7, color=FG, va="center", ha="center",
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7,
                            shrinkA=6, shrinkB=3),
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="#bbbbbb", linewidth=0.8, alpha=0.92),
            zorder=5,
        )

    # Legend by phase
    handles = []
    for phase in [0, 1, 2, 3, 4]:
        if phase in phases:
            handles.append(mpatches.Patch(
                color=_PHASE_COLORS[phase], label=_PHASE_LABELS[phase]))
    ax.legend(handles=handles, fontsize=7.5, title="Max clinical phase",
              title_fontsize=7.5, frameon=True, framealpha=1,
              edgecolor=GRID, loc="upper right")

    ax.set_xlabel("ChEMBL drug molecules with binding data  (+ jitter)", fontsize=8)
    ax.set_ylabel("DGIdb drug–gene interactions  (+ jitter)", fontsize=8)
    ax.set_title("Drug target evidence: ChEMBL vs DGIdb\n"
                 "Druggable genes only  |  dot size ∝ clinical phase", fontsize=9, loc="left")
    _style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)

    fig.tight_layout()
    save(fig, "07_drug_scatter.png")


def plot_druggability_by_changetype(df: pd.DataFrame) -> None:
    """08 — Stacked bar: druggability tier breakdown per protein change type."""
    _CT_DISPLAY = [
        "N_truncation", "C_truncation", "internal_indel",
        "N_extension", "C_extension", "internal_insertion",
        "substitution", "identical",
    ]
    _CT_LABELS = {
        "N_truncation":       "N trunc",
        "C_truncation":       "C trunc",
        "internal_indel":     "internal\nindel",
        "N_extension":        "N ext",
        "C_extension":        "C ext",
        "internal_insertion": "internal\nins",
        "substitution":       "substitution",
        "identical":          "identical",
    }

    _TIER_ORDER  = ["phase4", "phase13", "evidence", "none"]
    _TIER_COLORS = {
        "phase4":   "#e63946",
        "phase13":  "#fb8500",
        "evidence": ANNOT_TEAL,
        "none":     "#e0e0e0",
    }
    _TIER_LABELS = {
        "phase4":   "Phase 4 (approved)",
        "phase13":  "Phase 1–3 (clinical)",
        "evidence": "Target / interactions (no Rx)",
        "none":     "No drug evidence",
    }

    df2 = df.copy()
    df2["tier"] = df2.apply(_drug_tier, axis=1)
    df2 = df2[df2["protein_change_type"].isin(_CT_DISPLAY)]

    mat = (
        df2.groupby(["protein_change_type", "tier"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=_CT_DISPLAY, fill_value=0)
    )
    for t in _TIER_ORDER:
        if t not in mat.columns:
            mat[t] = 0
    mat = mat[_TIER_ORDER]

    fig, ax = plt.subplots(figsize=(10, 6))
    xpos   = np.arange(len(_CT_DISPLAY))
    bottom = np.zeros(len(_CT_DISPLAY))
    handles = []

    for tier in _TIER_ORDER:
        vals = mat[tier].values.astype(float)
        ax.bar(xpos, vals, bottom=bottom, color=_TIER_COLORS[tier],
               width=0.65, zorder=2, label=_TIER_LABELS[tier])
        handles.append(mpatches.Patch(color=_TIER_COLORS[tier],
                                      label=_TIER_LABELS[tier]))
        bottom += vals

    # Add druggable-fraction label above each bar
    totals = mat.sum(axis=1).values
    druggable = (mat["phase4"] + mat["phase13"] + mat["evidence"]).values
    for i, (tot, drg) in enumerate(zip(totals, druggable)):
        if tot > 0:
            pct = drg / tot * 100
            ax.text(i, tot + totals.max() * 0.01, f"{pct:.0f}%",
                    ha="center", va="bottom", fontsize=7, color=FG)

    ax.set_xticks(xpos)
    ax.set_xticklabels([_CT_LABELS[ct] for ct in _CT_DISPLAY], fontsize=8)
    ax.set_ylabel("Transcript count", fontsize=8)
    ax.set_title("Druggability tier per protein change type\n"
                 "% = fraction of transcripts affecting a gene with any drug evidence",
                 fontsize=9, loc="left")
    ax.legend(handles=handles[::-1], fontsize=7.5, frameon=True,
              framealpha=1, edgecolor=GRID, loc="upper right")

    _style(ax)
    ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.xaxis.grid(False)

    fig.tight_layout()
    save(fig, "08_druggability_by_changetype.png")


# ---------------------------------------------------------------------------
# Plots 09 / 10 / 11 — Cross-module: J1 × J2 × J3
# ---------------------------------------------------------------------------

_PROXY_COLORS = {
    "C":   ANNOT_TEAL,
    "D":   "#fb8500",
    "NMD": "#e63946",
}
_PROXY_LABELS = {
    "C":   "C — conservation",
    "D":   "D — differential expression",
    "NMD": "NMD — nonsense-mediated decay",
}

_TIER_COLORS_CROSS = {
    "Phase 4":     "#e63946",
    "Phase 1–3":   "#fb8500",
    "Target/DGIdb":"#007755",
    "None":        "#e0e0e0",
}


def _drug_tier_label(row: pd.Series) -> str:
    phase = int(row.get("chembl_max_phase", 0) or 0)
    if phase == 4:
        return "Phase 4"
    if phase in (1, 2, 3):
        return "Phase 1–3"
    if row.get("is_druggable") or int(row.get("dgidb_interactions", 0) or 0) > 0:
        return "Target/DGIdb"
    return "None"


def plot_proxy_changetype(df: pd.DataFrame) -> None:
    """09 — 100% stacked bar: which protein consequences does each proxy type detect?"""
    _CT_PLOT = [
        "C_truncation", "N_truncation",
        "internal_indel", "internal_insertion",
        "N_extension", "C_extension",
        "identical", "no_sequence",
    ]
    _CT_PLOT_LABELS = {
        "C_truncation":       "C truncation",
        "N_truncation":       "N truncation",
        "internal_indel":     "internal indel",
        "internal_insertion": "internal insertion",
        "N_extension":        "N extension",
        "C_extension":        "C extension",
        "identical":          "identical",
        "no_sequence":        "no sequence",
    }

    proxy_order = ["C", "D", "NMD"]
    counts = pd.crosstab(df["proxy_type"], df["protein_change_type"])
    counts = counts.reindex(index=proxy_order, fill_value=0)
    for ct in _CT_PLOT:
        if ct not in counts.columns:
            counts[ct] = 0
    counts = counts[_CT_PLOT]
    pct = counts.div(counts.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(8, 5))
    xpos    = np.arange(len(proxy_order))
    bottom  = np.zeros(len(proxy_order))
    handles = []

    for ct in _CT_PLOT:
        vals  = pct[ct].values
        color = CHANGE_COLORS.get(ct, "#cccccc")
        ax.bar(xpos, vals, bottom=bottom, color=color,
               width=0.55, zorder=2)
        handles.append(mpatches.Patch(color=color, label=_CT_PLOT_LABELS[ct]))
        # Label segments ≥ 5 % wide
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 0.05:
                ax.text(i, b + v / 2, f"{v*100:.0f}%",
                        ha="center", va="center", fontsize=7.5,
                        color="white" if v > 0.12 else FG, fontweight="bold")
        bottom += vals

    # Total-count annotation above bar
    totals = counts.sum(axis=1).values
    for i, n in enumerate(totals):
        ax.text(i, 1.02, f"n={n:,}", ha="center", va="bottom", fontsize=8, color=FG)

    ax.set_xticks(xpos)
    ax.set_xticklabels(
        [_PROXY_LABELS[p] for p in proxy_order], fontsize=8.5
    )
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Fraction of transcripts", fontsize=8)
    ax.set_title("Protein consequence profile per scoring proxy  (J1 × J2)\n"
                 "Each proxy detects a distinct structural biology",
                 fontsize=9, loc="left")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y*100:.0f}%"))
    ax.legend(handles=handles[::-1], fontsize=7, ncol=2,
              loc="upper right", frameon=True, framealpha=1, edgecolor=GRID)

    _style(ax)
    ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.xaxis.grid(False)
    ax.spines["bottom"].set_visible(True)

    fig.tight_layout()
    save(fig, "09_proxy_changetype.png")


def plot_deltausage_vs_structural(df: pd.DataFrame) -> None:
    """10 — Scatter: |delta_usage| vs changed_aa_fraction — decoupling expression and structure."""
    sub = df[df["protein_change_type"].isin(
        ["N_truncation", "C_truncation", "internal_indel",
         "N_extension", "C_extension", "internal_insertion", "substitution"]
    )].copy()

    rng = np.random.default_rng(0)
    jy  = rng.uniform(-0.005, 0.005, len(sub))

    fig, ax = plt.subplots(figsize=(9, 6))

    for proxy, grp in sub.groupby("proxy_type"):
        idx = grp.index
        jitter_slice = jy[sub.index.get_indexer(idx)]
        ax.scatter(
            grp["delta_usage"].abs(),
            grp["changed_aa_fraction"] + jitter_slice,
            c=_PROXY_COLORS.get(proxy, "#cccccc"),
            s=12, alpha=0.45, linewidths=0,
            label=_PROXY_LABELS.get(proxy, proxy), zorder=2,
        )

    # Per-proxy horizontal mean line for changed_aa_fraction
    for proxy, grp in sub.groupby("proxy_type"):
        mean_caf = grp["changed_aa_fraction"].mean()
        color = _PROXY_COLORS.get(proxy, "#cccccc")
        ax.axhline(mean_caf, color=color, lw=1.2,
                   ls="--", alpha=0.7, zorder=3)
        ax.text(0.97, mean_caf + 0.01, f"{proxy} mean",
                ha="right", va="bottom", fontsize=7,
                color=color, transform=ax.get_yaxis_transform())

    # Pearson r on non-identical transcripts
    r_val = sub[["delta_usage", "changed_aa_fraction"]].corr().iloc[0, 1]
    ax.text(0.97, 0.96, f"Pearson r = {r_val:.2f}\n(structural change vs. Δusage)",
            ha="right", va="top", transform=ax.transAxes,
            fontsize=8, color=FG,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor=GRID, linewidth=0.8))

    ax.set_xlabel("|Δ usage|  (AD − Control)", fontsize=8)
    ax.set_ylabel("Changed amino-acid fraction  (J2)", fontsize=8)
    ax.set_title("Differential splicing magnitude vs protein structural impact  (J1 × J2)\n"
                 "Near-zero correlation: expression shift and structural severity are independent",
                 fontsize=9, loc="left")
    ax.legend(fontsize=7.5, frameon=True, framealpha=1,
              edgecolor=GRID, loc="upper left")

    _style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)

    fig.tight_layout()
    save(fig, "10_deltausage_vs_structural.png")


def plot_proxy_druggability(df: pd.DataFrame) -> None:
    """11 — Stacked bar: druggability tier breakdown per scoring proxy (J1 × J3)."""
    tier_order  = ["Phase 4", "Phase 1–3", "Target/DGIdb", "None"]
    proxy_order = ["C", "D", "NMD"]

    df2 = df.copy()
    df2["drug_tier"] = df2.apply(_drug_tier_label, axis=1)

    counts = pd.crosstab(df2["proxy_type"], df2["drug_tier"])
    counts = counts.reindex(index=proxy_order, fill_value=0)
    for t in tier_order:
        if t not in counts.columns:
            counts[t] = 0
    counts = counts[tier_order]

    fig, (ax_abs, ax_pct) = plt.subplots(1, 2, figsize=(11, 5))

    for ax, data, ylabel, title_suffix in [
        (ax_abs, counts,
         "Transcript count", "absolute counts"),
        (ax_pct, counts.div(counts.sum(axis=1), axis=0),
         "Fraction of transcripts", "fractions"),
    ]:
        xpos   = np.arange(len(proxy_order))
        bottom = np.zeros(len(proxy_order))
        for tier in tier_order:
            vals  = data[tier].values.astype(float)
            color = _TIER_COLORS_CROSS[tier]
            ax.bar(xpos, vals, bottom=bottom, color=color, width=0.5, zorder=2)
            # Label visible segments
            threshold = (data.sum(axis=1).max() if ax is ax_abs else 1) * 0.04
            for i, (v, b) in enumerate(zip(vals, bottom)):
                if v >= threshold:
                    label = str(int(v)) if ax is ax_abs else f"{v*100:.0f}%"
                    ax.text(i, b + v / 2, label,
                            ha="center", va="center", fontsize=7.5,
                            color="white" if (ax is ax_pct and tier != "None") else FG,
                            fontweight="bold")
            bottom += vals

        ax.set_xticks(xpos)
        ax.set_xticklabels([_PROXY_LABELS[p] for p in proxy_order],
                           fontsize=8, rotation=15, ha="right")
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_title(title_suffix.capitalize(), fontsize=8.5)
        if ax is ax_pct:
            ax.set_ylim(0, 1.08)
            ax.yaxis.set_major_formatter(
                plt.FuncFormatter(lambda y, _: f"{y*100:.0f}%"))
        _style(ax)
        ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)
        ax.xaxis.grid(False)

    # Shared legend
    handles = [mpatches.Patch(color=_TIER_COLORS_CROSS[t], label=t)
               for t in tier_order]
    fig.legend(handles=handles[::-1], fontsize=8, loc="upper right",
               bbox_to_anchor=(0.99, 0.98), frameon=True, framealpha=1,
               edgecolor=GRID, title="Druggability tier", title_fontsize=8)

    fig.suptitle("Druggability profile per scoring proxy  (J1 × J3)\n"
                 "Conservation proxy (C) finds proportionally more drug targets",
                 fontsize=9, x=0.42, y=1.02)
    fig.tight_layout()
    save(fig, "11_proxy_druggability.png")


# ---------------------------------------------------------------------------
# Plot 12 — Sankey: biotype_class -> proxy_type -> junior gate -> J4 gate outcome
# ---------------------------------------------------------------------------

def plot_sankey_full(df: pd.DataFrame) -> None:
    """12 — Sankey: significant hits -> candidate group -> biotype_class -> junior gate -> J4 gate outcome.

    Stage 0->1 (significant hits -> candidate group) uses ALL 2,599
    permutation-significant DTU hits (config.ALL_GROUPS_CSV, written by
    initial_filter.py) — includes the two dead-end groups ("other" and
    "no_MANE_coverage") that never reached assistant_surveyor at all, so
    they have inbound flow only, no outbound.

    Stages 1->3 (candidate group -> biotype -> junior gate) are computed
    from assistant_surveyor's full hits_enriched.csv (only the two real
    groups that were actually passed to assistant_surveyor -- "other" and
    "no_MANE_coverage" never appear here, so those links resolve to zero
    automatically), not junior_surveyor's already-J1-filtered df. Stage 4
    (junior gate -> J4 gate outcome) uses df, since gate_outcome only exists
    for hits that actually reached J1-J4 (i.e. passed the junior gate).
    """
    import plotly.graph_objects as go

    from junior_surveyor.config import ALL_GROUPS_CSV, HITS_CSV, NULL_PROTEIN_CHANGE_TYPES

    def _rgba(hex_color: str, alpha: float = 0.35) -> str:
        r = int(hex_color[1:3], 16)
        g = int(hex_color[3:5], 16)
        b = int(hex_color[5:7], 16)
        return f"rgba({r},{g},{b},{alpha})"

    all_groups_df = pd.read_csv(ALL_GROUPS_CSV, usecols=["candidate_group"])
    full_df = pd.read_csv(HITS_CSV, usecols=["candidate_group", "biotype_class", "junior_pass"])

    has_protein_change = ~df["protein_change_type"].isin(NULL_PROTEIN_CHANGE_TYPES)
    has_drug_evidence   = (df["chembl_max_phase"] >= 1) | (df["dgidb_interactions"] > 0)

    df = df.copy()
    df["gate_outcome"] = np.select(
        [has_protein_change & has_drug_evidence,
         has_protein_change & ~has_drug_evidence,
         ~has_protein_change & has_drug_evidence],
        ["Both pass", "Gate 1 only", "Gate 2 only"],
        default="Neither",
    )

    group_order = ["trial_failure_candidate", "new_target_candidate", "other", "no_MANE_coverage"]
    # Drop biotype nodes with zero hits (e.g. TEC when none are present). A
    # disconnected zero-flow node in a middle column corrupts Plotly's
    # arrangement="fixed" layout and shoves the downstream "Junior pass" node
    # back into the biotype column, so it must be excluded, not just hidden.
    bio_order   = [b for b in BIO_ORDER if int((full_df["biotype_class"] == b).sum()) > 0]
    jgate_order = [True, False]
    jgate_labels = {True: "Junior pass", False: "Junior drop"}
    gate_order  = ["Both pass", "Gate 1 only", "Gate 2 only", "Neither"]

    n_hits_node = 1
    n_grp, n_bio, n_jgate = len(group_order), len(bio_order), len(jgate_order)

    hits_idx  = 0
    group_idx = {g: n_hits_node + i                        for i, g in enumerate(group_order)}
    bio_idx   = {b: n_hits_node + n_grp + i                for i, b in enumerate(bio_order)}
    jgate_idx = {t: n_hits_node + n_grp + n_bio + i        for i, t in enumerate(jgate_order)}
    gate_idx  = {g: n_hits_node + n_grp + n_bio + n_jgate + i for i, g in enumerate(gate_order)}

    group_labels = {
        "trial_failure_candidate": "Trial-failure candidate",
        "new_target_candidate":    "New-target candidate",
        "other":                   "Other",
        "no_MANE_coverage":        "No MANE annotation",
    }

    # True per-node counts for the displayed label. Candidate group uses the
    # FULL 2,599-hit universe (all_groups_df); biotype/junior-gate use the
    # FULL assistant_surveyor universe (full_df, 1,201 hits); J4 gate outcome
    # uses junior's df.
    n_total_hits = len(all_groups_df)
    true_group = {g: int((all_groups_df["candidate_group"] == g).sum()) for g in group_order}
    true_bio   = {b: int((full_df["biotype_class"] == b).sum()) for b in bio_order}
    true_jgate = {t: int((full_df["junior_pass"] == t).sum())   for t in jgate_order}
    true_gate  = {g: int((df["gate_outcome"] == g).sum())       for g in gate_order}
    true_counts = (
        [n_total_hits]
        + [true_group[g] for g in group_order]
        + [true_bio[b] for b in bio_order]
        + [true_jgate[t] for t in jgate_order]
        + [true_gate[g] for g in gate_order]
    )

    node_labels = [
        f"{lbl} ({n:,})" for lbl, n in zip(
            ["Significant hits"]
            + [group_labels[g] for g in group_order]
            + bio_order + [jgate_labels[t] for t in jgate_order] + gate_order,
            true_counts,
        )
    ]
    node_colors = (
        [HITS_NODE_COLOR]
        + [GROUP_COLORS[g]       for g in group_order]
        + [BIO_COLORS[b]         for b in bio_order]
        + [JUNIOR_GATE_COLORS[t] for t in jgate_order]
        + [GATE_COLORS[g]        for g in gate_order]
    )

    sources, targets, values, link_colors = [], [], [], []

    # significant hits -> candidate group  (full 2,599-hit universe)
    for grp in group_order:
        n = true_group[grp]
        if n:
            sources.append(hits_idx)
            targets.append(group_idx[grp])
            values.append(n)
            link_colors.append(_rgba(GROUP_COLORS[grp]))

    # candidate group -> biotype  (assistant_surveyor universe — "other" and
    # "no_MANE_coverage" never appear in full_df, so these resolve to zero
    # and are true dead ends, not a rendering artifact)
    for grp in group_order:
        for bio in bio_order:
            n = int(((full_df["candidate_group"] == grp) & (full_df["biotype_class"] == bio)).sum())
            if n:
                sources.append(group_idx[grp])
                targets.append(bio_idx[bio])
                values.append(n)
                link_colors.append(_rgba(GROUP_COLORS[grp]))

    # biotype -> junior gate  (direct — proxy_type stage removed)
    for bio in bio_order:
        for t in jgate_order:
            n = int(((full_df["biotype_class"] == bio) & (full_df["junior_pass"] == t)).sum())
            if n:
                sources.append(bio_idx[bio])
                targets.append(jgate_idx[t])
                values.append(n)
                link_colors.append(_rgba(BIO_COLORS[bio]))

    # junior gate -> J4 gate outcome  (J1-J4 survivors only — Junior drop has
    # none: real dead end, not a rendering artifact.)
    for t in jgate_order:
        for g in gate_order:
            n = int(((df["junior_pass"] == t) & (df["gate_outcome"] == g)).sum())
            if n:
                sources.append(jgate_idx[t])
                targets.append(gate_idx[g])
                values.append(n)
                link_colors.append(_rgba(JUNIOR_GATE_COLORS[t]))

    # Pin every node's position explicitly. Plotly's default "snap"
    # arrangement (like d3-sankey's justify) auto-places nodes by graph
    # depth and shoves any node with no outgoing links — Junior drop, a true
    # dead end — flush into the rightmost column, mixing it visually in
    # with the gate outcomes. "perpendicular" (x fixed, y free) turned out
    # not to hold x either, so lay out both axes by hand under "fixed".
    #
    # Node *thickness* is rendered by Plotly from one shared value-per-pixel
    # scale across the whole figure (so a link's ribbon width matches at
    # both ends) — not independently per column. Normalizing each column's
    # y-positions against its own total (as before) assumes a per-column
    # scale that Plotly doesn't actually use, which desyncs our y-centers
    # from the real rendered thickness: columns misalign at the top and
    # gaps come out uneven. Using one GLOBAL_TOTAL (the gate column has
    # fewer hits than bio/proxy/junior-gate — J1 already dropped the rest)
    # and a constant pad for every column fixes both: every column's first
    # node starts at y=0 (tops aligned), and equal pad means equal gaps.
    GLOBAL_TOTAL = max(
        n_total_hits, sum(true_group.values()), sum(true_bio.values()),
        sum(true_jgate.values()), sum(true_gate.values()),
    )
    max_nodes = max(n_hits_node, n_grp, n_bio, n_jgate, len(gate_order))
    PAD = 0.02
    avail = 1.0 - PAD * (max_nodes - 1)
    scale = avail / GLOBAL_TOTAL

    def _column_y(order, value_of):
        cursor, ys = 0.0, []
        for k in order:
            h = value_of(k) * scale
            ys.append(cursor + h / 2)
            cursor += h + PAD
        return ys

    col_x = [0.001, 0.2505, 0.5, 0.7495, 0.999]
    node_x = (
        [col_x[0]] * n_hits_node
        + [col_x[1]] * n_grp
        + [col_x[2]] * n_bio
        + [col_x[3]] * n_jgate
        + [col_x[4]] * len(gate_order)
    )
    node_y = (
        _column_y(["Significant hits"], lambda k: n_total_hits)
        + _column_y(group_order, lambda k: true_group[k])
        + _column_y(bio_order, lambda k: true_bio[k])
        + _column_y(jgate_order, lambda k: true_jgate[k])
        + _column_y(gate_order, lambda k: true_gate[k])
    )

    fig = go.Figure(go.Sankey(
        arrangement="fixed",
        node=dict(
            pad=18,
            thickness=22,
            line=dict(color="#cccccc", width=0.5),
            label=node_labels,
            color=node_colors,
            x=node_x,
            y=node_y,
            hovertemplate="%{label}<extra></extra>",
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
            text="Significant hits → Candidate group → Biotype class → Junior gate → J4 gate outcome",
            font=dict(size=15, family="sans-serif", color=FG),
            x=0.5, xanchor="center",
        ),
        font=dict(family="sans-serif", size=12, color=FG),
        paper_bgcolor=BG,
        plot_bgcolor=BG,
        width=1450,
        height=650,
        margin=dict(l=20, r=20, t=60, b=20),
    )

    out = OUT_DIR / "12_sankey_hits_group_biotype_juniorgate_j4gate.png"
    fig.write_image(str(out), scale=2)
    print(f"  saved → {out.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading hits_deep.csv …")
    df = load()
    print(f"  {len(df):,} hits  |  {df['gene_name'].nunique():,} genes\n")

    print("Rendering plots …")
    plot_protein_lengths(df)
    plot_pfam_domains(df)
    plot_affected_domain_heatmap(df)
    plot_domain_loss_gain(df)
    plot_domain_loss_gain_detail(df)
    plot_druggability_funnel(df)
    plot_drug_scatter(df)
    plot_druggability_by_changetype(df)
    plot_proxy_changetype(df)
    plot_deltausage_vs_structural(df)
    plot_proxy_druggability(df)
    plot_sankey_full(df)

    n = len(list(OUT_DIR.glob("*.png")))
    print(f"\nDone — {n} plot(s) in {OUT_DIR}")


if __name__ == "__main__":
    main()
