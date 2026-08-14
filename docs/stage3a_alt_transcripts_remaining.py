"""Stage 3A summary plot — surviving alt transcripts per trial_failure hit.

One-off analysis script (not part of the junior_surveyor pipeline), matching
the color/style conventions of junior_surveyor/plot_results.py.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BG, FG, GRID = "#ffffff", "#1a1a1a", "#e0e0e0"
TF_BLUE = "#2a78d6"

df = pd.read_csv("outputs/junior_surveyor/hits_deep.csv")
tf = df[df["candidate_group"] == "trial_failure_candidate"]
tf = tf[tf["has_viable_alt"]]  # exclude the 127 hits with no surviving non-canonical alt
tf_alt = tf[~tf["is_canonical"].fillna(False).astype(bool)]

counts = tf_alt.groupby(["gene_name", "cell_type"]).size()
canon_hits = tf[tf["is_canonical"] == True][["gene_name", "cell_type"]].drop_duplicates()
counts_full = counts.reindex(
    canon_hits.set_index(["gene_name", "cell_type"]).index, fill_value=0
)

dist = counts_full.value_counts().sort_index()
n_hits = len(counts_full)

fig, ax = plt.subplots(figsize=(7.5, 5))
bars = ax.bar(dist.index.astype(str), dist.values, color=TF_BLUE, width=0.62, zorder=2)

for rect, v in zip(bars, dist.values):
    pct = v / n_hits * 100
    ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + n_hits * 0.008,
            f"{v}\n({pct:.0f}%)", ha="center", va="bottom", fontsize=8.5, color=FG)

ax.set_facecolor(BG)
fig.patch.set_facecolor(BG)
ax.set_xlabel("surviving alt transcripts per hit\n(after alt_biotype_class + 3% AD-usage filter)", fontsize=9)
ax.set_ylabel("number of trial_failure hits", fontsize=9)
ax.set_title(
    f"Stage 3A: alt transcripts remaining per trial_failure hit\n"
    f"n = {n_hits} hits with ≥1 viable alt (127 hits with none excluded)\n"
    f"mean = {counts_full.mean():.2f}   median = {counts_full.median():.0f}   max = {counts_full.max()}",
    fontsize=10.5, pad=10,
)
ax.tick_params(colors=FG, labelsize=9)
ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)
ax.xaxis.grid(False)
ax.set_ylim(0, dist.values.max() * 1.18)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["left", "bottom"]:
    ax.spines[spine].set_color(GRID)

fig.tight_layout()
fig.savefig("docs/stage3a_alt_transcripts_remaining.png", dpi=150, facecolor=BG)
print("saved -> docs/stage3a_alt_transcripts_remaining.png")
print(dist)
print(f"mean={counts_full.mean():.3f} median={counts_full.median()} max={counts_full.max()}")
