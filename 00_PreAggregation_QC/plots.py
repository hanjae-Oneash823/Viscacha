"""
00_PreAggregation_QC visualizations — styled to match project notebook conventions.
All figures saved to outputs/00_PreAggregation_QC/plots/.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import seaborn as sns
import anndata as ad
from pathlib import Path

from .config import OUT_PLOTS, COND_AD, COND_CTRL, COND_ACTIVE

# --- Project-standard style ---
COND_ORDER  = [COND_CTRL, COND_ACTIVE, COND_AD]
COND_COLORS = {COND_AD: '#d62728', COND_CTRL: '#2ca02c', COND_ACTIVE: '#ff7f0e'}
BOX_OVERLAY = dict(width=0.08, color='#2C3E50', fliersize=0,
                   boxprops={'zorder': 3}, linewidth=1.5)
SEX_MARKERS = {'M': '^', 'F': 'o'}


def _savefig(fig: plt.Figure, name: str) -> Path:
    OUT_PLOTS.mkdir(parents=True, exist_ok=True)
    path = OUT_PLOTS / name
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  [plot] saved {path.name}')
    return path


# ---------------------------------------------------------------------------
# Step 0.1 — Cohort overview
# ---------------------------------------------------------------------------

def plot_step01(meta: pd.DataFrame) -> None:
    """
    01a_cohort_counts.png   — donor counts by condition, stacked by source (SMC / PO)
    01b_age_distribution.png — age by condition, markers = sex
    """
    df = meta.reset_index().copy()
    df['source'] = df['donor_id'].str[:2].map({'SM': 'SMC', 'PO': 'PO'})

    # --- 01a ---
    fig, ax = plt.subplots(figsize=(6, 5))
    grp = (df.groupby(['condition', 'source'])
             .size()
             .unstack(fill_value=0)
             .reindex([c for c in COND_ORDER if c in df['condition'].values]))
    grp.plot(kind='bar', ax=ax,
             color=['#5D6D7E', '#A9CCE3'], edgecolor='white', width=0.6)
    ax.set_xlabel('')
    ax.set_ylabel('Number of donors', labelpad=10)
    ax.set_title('Cohort composition', weight='semibold', pad=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.legend(title='Source', frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, '01a_cohort_counts.png')

    # --- 01b ---
    df_v = df.dropna(subset=['age', 'condition', 'sex'])
    fig, ax = plt.subplots(figsize=(7, 6))
    for cond_i, cond in enumerate(COND_ORDER):
        for sex, marker in SEX_MARKERS.items():
            pts = df_v[(df_v['condition'] == cond) & (df_v['sex'] == sex)]['age']
            jitter = np.random.default_rng(42).uniform(-0.14, 0.14, len(pts))
            ax.scatter(cond_i + jitter, pts,
                       color=COND_COLORS.get(cond, 'gray'),
                       marker=marker, s=65, alpha=0.85, zorder=3,
                       label=sex if cond_i == 0 else '_nolegend_')
    bp_data = [df_v[df_v['condition'] == c]['age'].dropna().values for c in COND_ORDER]
    ax.boxplot(bp_data, positions=range(len(COND_ORDER)),
               widths=0.3, patch_artist=False,
               medianprops=dict(color='black', linewidth=2),
               whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1),
               flierprops=dict(marker=''))
    ax.set_xticks(range(len(COND_ORDER)))
    ax.set_xticklabels(COND_ORDER, rotation=20, ha='right')
    ax.set_ylabel('Age at death (years)', labelpad=10)
    ax.set_title('Age distribution by condition', weight='semibold', pad=12)
    ax.legend(title='Sex', frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, '01b_age_distribution.png')


# ---------------------------------------------------------------------------
# Step 0.3 — Barcode confidence filter
# ---------------------------------------------------------------------------

def plot_step03(adata_pre: ad.AnnData, adata_post: ad.AnnData) -> None:
    """
    03a_doublet_score_violin.png   — doublet score per cell type with threshold
    03b_barcode_dropout_heatmap.png — % barcodes dropped per donor × cell type
    """
    obs = adata_pre.obs.copy()
    obs['kept'] = obs.index.isin(adata_post.obs_names)

    cell_types_sorted = (
        obs.groupby('cell_type')['doublet_score']
           .median()
           .sort_values()
           .index.tolist()
    )

    # --- 03a: violin + narrow boxplot overlay ---
    plot_data = [obs[obs['cell_type'] == ct]['doublet_score'].dropna().values
                 for ct in cell_types_sorted]
    fig, ax = plt.subplots(figsize=(11, 5))
    parts = ax.violinplot(plot_data, positions=range(len(cell_types_sorted)),
                          showmedians=False, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor('#AED6F1')
        pc.set_alpha(0.75)
        pc.set_edgecolor('#2980B9')
        pc.set_linewidth(0.8)
    # narrow boxplot overlay (project style)
    ax.boxplot(plot_data, positions=range(len(cell_types_sorted)),
               widths=0.05, patch_artist=True,
               boxprops=dict(facecolor='#2C3E50', zorder=3),
               medianprops=dict(color='white', linewidth=1.5),
               whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1),
               flierprops=dict(marker=''))
    ax.axhline(0.3, color='#d62728', linestyle='--', linewidth=1.4,
               label='threshold (0.3)')
    ax.set_xticks(range(len(cell_types_sorted)))
    ax.set_xticklabels(cell_types_sorted, rotation=30, ha='right')
    ax.set_ylabel('Doublet score', labelpad=10)
    ax.set_title('Doublet score distribution per cell type\n(sorted by median)',
                 weight='semibold', pad=12)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, '03a_doublet_score_violin.png')

    # --- 03b: dropout heatmap ---
    donors     = sorted(obs['donor'].unique())
    cond_map   = obs.drop_duplicates('donor').set_index('donor')['condition'].to_dict()
    cond_abbr  = {COND_AD: 'AD', COND_CTRL: 'Ct', COND_ACTIVE: 'Ac'}

    donor_order = (
        [d for d in donors if cond_map.get(d) == COND_CTRL] +
        [d for d in donors if cond_map.get(d) == COND_ACTIVE] +
        [d for d in donors if cond_map.get(d) == COND_AD]
    )
    row_labels = [f'{d} [{cond_abbr.get(cond_map.get(d,""), "?")}]'
                  for d in donor_order]

    matrix = pd.DataFrame(index=donor_order, columns=cell_types_sorted, dtype=float)
    for ct in cell_types_sorted:
        for donor in donor_order:
            mask = (obs['cell_type'] == ct) & (obs['donor'] == donor)
            n_total = mask.sum()
            if n_total == 0:
                matrix.loc[donor, ct] = np.nan
            else:
                matrix.loc[donor, ct] = (~obs.loc[mask, 'kept']).sum() / n_total * 100

    fig, ax = plt.subplots(figsize=(11, 9))
    sns.heatmap(matrix.astype(float), ax=ax,
                cmap='YlOrRd', vmin=0, vmax=50,
                linewidths=0.3, linecolor='white',
                yticklabels=row_labels,
                xticklabels=[c.replace(' ', '\n') for c in cell_types_sorted],
                cbar_kws={'label': '% barcodes dropped'},
                annot=True, fmt='.0f', annot_kws={'size': 7})
    n_ctrl  = sum(1 for d in donor_order if cond_map.get(d) == COND_CTRL)
    n_ac    = sum(1 for d in donor_order if cond_map.get(d) == COND_ACTIVE)
    ax.axhline(n_ctrl,        color='black', linewidth=1.5)
    ax.axhline(n_ctrl + n_ac, color='black', linewidth=1.5)
    ax.set_title('Barcode dropout % per donor × cell type  (Step 0.3 filter)',
                 weight='semibold', pad=12)
    fig.tight_layout()
    _savefig(fig, '03b_barcode_dropout_heatmap.png')


# ---------------------------------------------------------------------------
# Step 0.4 — Prevalence filter
# ---------------------------------------------------------------------------

def plot_step04(prevalence_masks: dict, n_total: int) -> None:
    """
    04_transcript_prevalence_filter.png — transcripts kept / dropped per cell type
    """
    cell_types = list(prevalence_masks.keys())
    n_kept    = [len(prevalence_masks[ct]) for ct in cell_types]
    n_dropped = [n_total - k for k in n_kept]
    pct_kept  = [k / n_total * 100 for k in n_kept]

    order = sorted(range(len(cell_types)), key=lambda i: n_kept[i], reverse=True)
    cell_types = [cell_types[i] for i in order]
    n_kept     = [n_kept[i]    for i in order]
    n_dropped  = [n_dropped[i] for i in order]
    pct_kept   = [pct_kept[i]  for i in order]

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(cell_types))
    ax.bar(x, n_kept,    color='steelblue',   label='Kept',    width=0.6, edgecolor='white')
    ax.bar(x, n_dropped, bottom=n_kept,
           color='#DCDCDC', label='Dropped',  width=0.6, edgecolor='white')
    for i, (k, pct) in enumerate(zip(n_kept, pct_kept)):
        ax.text(i, k + 400, f'{pct:.0f}%', ha='center', va='bottom', fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=30, ha='right')
    ax.set_ylabel('Number of transcripts', labelpad=10)
    ax.set_title(f'Transcript prevalence filter  (≥40% donors)\nn_total = {n_total:,}',
                 weight='semibold', pad=12)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, '04_transcript_prevalence_filter.png')


# ---------------------------------------------------------------------------
# Step 0.6 — Pseudo-bulk outputs
# ---------------------------------------------------------------------------

def plot_step06(adata: ad.AnnData) -> None:
    """
    06a_ncells_heatmap.png    — n_cells per donor × cell type
    06b_median_pct_mt.png     — median_pct_mt violin per cell type by condition
    06c_library_size.png      — library size violin per cell type by condition
    """
    from .config import OUT_PB, CELL_TYPES

    obs = adata.obs.copy()
    cell_types = [ct for ct in CELL_TYPES if ct in obs['cell_type'].values]
    donors     = sorted(obs['donor'].unique())
    cond_map   = obs.drop_duplicates('donor').set_index('donor')['condition'].to_dict()
    cond_abbr  = {COND_AD: 'AD', COND_CTRL: 'Ct', COND_ACTIVE: 'Ac'}

    donor_order = (
        [d for d in donors if cond_map.get(d) == COND_CTRL] +
        [d for d in donors if cond_map.get(d) == COND_ACTIVE] +
        [d for d in donors if cond_map.get(d) == COND_AD]
    )
    row_labels = [f'{d} [{cond_abbr.get(cond_map.get(d,""), "?")}]'
                  for d in donor_order]

    # --- 06a: n_cells heatmap ---
    ncells_mat = pd.DataFrame(0, index=donor_order, columns=cell_types, dtype=float)
    for ct in cell_types:
        for donor in donor_order:
            n = ((obs['cell_type'] == ct) & (obs['donor'] == donor)).sum()
            ncells_mat.loc[donor, ct] = float(n) if n > 0 else np.nan

    fig, ax = plt.subplots(figsize=(11, 9))
    im = ax.imshow(ncells_mat.values, cmap='Blues', aspect='auto',
                   vmin=0, vmax=np.nanmax(ncells_mat.values))
    ax.set_xticks(range(len(cell_types)))
    ax.set_yticks(range(len(donor_order)))
    ax.set_xticklabels([c.replace(' ', '\n') for c in cell_types], fontsize=9)
    ax.set_yticklabels(row_labels, fontsize=8)
    # annotate cells
    for i in range(len(donor_order)):
        for j in range(len(cell_types)):
            val = ncells_mat.iloc[i, j]
            if not np.isnan(val):
                ax.text(j, i, f'{int(val):,}', ha='center', va='center',
                        fontsize=6.5, color='white' if val > ncells_mat.values[~np.isnan(ncells_mat.values)].max() * 0.55 else '#333333')
    n_ctrl = sum(1 for d in donor_order if cond_map.get(d) == COND_CTRL)
    n_ac   = sum(1 for d in donor_order if cond_map.get(d) == COND_ACTIVE)
    ax.axhline(n_ctrl - 0.5,        color='white', linewidth=2)
    ax.axhline(n_ctrl + n_ac - 0.5, color='white', linewidth=2)
    plt.colorbar(im, ax=ax, label='n cells', shrink=0.7)
    ax.set_title('Cells per donor × cell type  (after Step 0.3 filter)',
                 weight='semibold', pad=12)
    fig.tight_layout()
    _savefig(fig, '06a_ncells_heatmap.png')

    # --- 06b: pct_mt violin + narrow boxplot overlay ---
    records = []
    for ct in cell_types:
        for donor in donors:
            sub = obs[(obs['cell_type'] == ct) & (obs['donor'] == donor)]
            if len(sub) == 0:
                continue
            records.append({
                'cell_type': ct,
                'condition': sub['condition'].iloc[0],
                'median_pct_mt': sub['pct_counts_mt'].median(),
            })
    mt_df = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.violinplot(
        data=mt_df, x='cell_type', y='median_pct_mt', hue='condition',
        hue_order=COND_ORDER, palette=COND_COLORS,
        order=cell_types, ax=ax,
        inner=None, cut=0, linewidth=0.8,
    )
    sns.stripplot(
        data=mt_df, x='cell_type', y='median_pct_mt', hue='condition',
        hue_order=COND_ORDER, palette=COND_COLORS,
        order=cell_types, ax=ax,
        dodge=True, size=4, alpha=0.65, color='black', legend=False,
    )
    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=30, ha='right')
    ax.set_xlabel('')
    ax.set_ylabel('Median % mitochondrial reads', labelpad=10)
    ax.set_title('RNA quality proxy  (pct_counts_mt)  per donor × cell type',
                 weight='semibold', pad=12)
    handles = [mpatches.Patch(color=COND_COLORS[c], label=c) for c in COND_ORDER]
    ax.legend(handles=handles, title='Condition', frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, '06b_median_pct_mt.png')

    # --- 06c: library size violin + strip ---
    lib_records = []
    for ct in cell_types:
        stem = ct.replace(' ', '_')
        for suffix in ['', '_sensitivity']:
            path = OUT_PB / f'counts_{stem}{suffix}.csv'
            if not path.exists():
                continue
            df = pd.read_csv(path, index_col=0)
            for donor in df.index:
                lib_records.append({
                    'cell_type': ct,
                    'condition': cond_map.get(donor, 'Unknown'),
                    'library_size': int(df.loc[donor].sum()),
                })
    lib_df = pd.DataFrame(lib_records)

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.violinplot(
        data=lib_df, x='cell_type', y='library_size', hue='condition',
        hue_order=COND_ORDER, palette=COND_COLORS,
        order=cell_types, ax=ax,
        inner=None, cut=0, linewidth=0.8,
    )
    sns.stripplot(
        data=lib_df, x='cell_type', y='library_size', hue='condition',
        hue_order=COND_ORDER, palette=COND_COLORS,
        order=cell_types, ax=ax,
        dodge=True, size=4, alpha=0.65, color='black', legend=False,
    )
    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=30, ha='right')
    ax.set_xlabel('')
    ax.set_ylabel('Total transcript counts  (library size)', labelpad=10)
    ax.set_title('Pseudo-bulk library size per cell type by condition',
                 weight='semibold', pad=12)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    handles = [mpatches.Patch(color=COND_COLORS[c], label=c) for c in COND_ORDER]
    ax.legend(handles=handles, title='Condition', frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, '06c_library_size.png')
