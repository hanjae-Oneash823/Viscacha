#!/usr/bin/env python3
"""
Export permutation_significant=True DIU hits to a single combined CSV.

Runs before INITIAL_FILTER (repo root) -> ASSISTANT_SURVEYOR.

Columns added from h5ad pseudo-bulk:
  total_counts_AD / total_counts_CT  — raw UMI sums across all cells of that cell type
  mean_CPM_AD    / mean_CPM_CT       — mean CPM across donors (used for log2FC)
  log2FC                             — log2(mean_CPM_AD + 1) - log2(mean_CPM_CT + 1)

Active control donors excluded. One row per significant transcript × cell type.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python 02_SURVEYOR/export_significant_hits.py
"""

import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

DATA_DIR   = Path('/node212data/welcome3/Grad_proj_2026/DATA')
DIU_DIR    = DATA_DIR / 'DIU_result_with_permutation'
OUTPUT_DIR = Path('/home/welcome3/Viscacha_pipeline/outputs/DIU_significant_hits')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CT_MAP = {
    'Astrocyte':         'Astrocyte',
    'Excitatory neuron': 'Excitatory_neuron',
    'Inhibitory neuron': 'Inhibitory_neuron',
    'Lymphocyte':        'Lymphocyte',
    'Microglia':         'Microglia',
    'Oligodendrocyte':   'Oligodendrocyte',
    'OPC':               'OPC',
    'Vascular cell':     'Vascular_cell',
}

# ── 1. Load h5ad ──────────────────────────────────────────────────────────────
print("Loading h5ad...", flush=True)
adata = ad.read_h5ad(DATA_DIR / 'adata_transcript_loose_filtering_for_bulk_analysis.h5ad')
adata = adata[adata.obs['condition'].isin(['AD', 'Control'])].copy()
print(f"  {adata.n_obs} cells × {adata.n_vars} transcripts (AD + Control only)", flush=True)

# ── 2. Pseudo-bulk per cell type ──────────────────────────────────────────────
print("Computing pseudo-bulk per cell type...", flush=True)
expr_dict = {}   # ct_csv → DataFrame indexed by transcript_name

for ct_h5ad, ct_csv in CT_MAP.items():
    adata_ct  = adata[adata.obs['cell_type'] == ct_h5ad]
    donor_arr = adata_ct.obs['donor'].values
    cond_arr  = adata_ct.obs['condition'].values
    X         = adata_ct.X

    donors_ad = adata_ct.obs.loc[cond_arr == 'AD',      'donor'].unique()
    donors_ct = adata_ct.obs.loc[cond_arr == 'Control', 'donor'].unique()

    pb_ad = np.zeros((len(donors_ad), adata_ct.n_vars))
    pb_ct = np.zeros((len(donors_ct), adata_ct.n_vars))

    for i, d in enumerate(donors_ad):
        pb_ad[i] = np.asarray(X[donor_arr == d].sum(axis=0)).flatten()
    for i, d in enumerate(donors_ct):
        pb_ct[i] = np.asarray(X[donor_arr == d].sum(axis=0)).flatten()

    # Raw totals across all cells in each condition
    total_ad = pb_ad.sum(axis=0)
    total_ct = pb_ct.sum(axis=0)

    # CPM per donor then mean across donors
    row_sum_ad = pb_ad.sum(axis=1, keepdims=True); row_sum_ad[row_sum_ad == 0] = 1
    row_sum_ct = pb_ct.sum(axis=1, keepdims=True); row_sum_ct[row_sum_ct == 0] = 1
    mean_cpm_ad = (pb_ad / row_sum_ad * 1e6).mean(axis=0)
    mean_cpm_ct = (pb_ct / row_sum_ct * 1e6).mean(axis=0)

    log2fc = np.log2(mean_cpm_ad + 1) - np.log2(mean_cpm_ct + 1)

    expr_dict[ct_csv] = pd.DataFrame({
        'total_counts_AD': total_ad,
        'total_counts_CT': total_ct,
        'mean_CPM_AD':     mean_cpm_ad,
        'mean_CPM_CT':     mean_cpm_ct,
        'log2FC':          log2fc,
    }, index=adata_ct.var_names)

    n_ad_donors = len(donors_ad)
    n_ct_donors = len(donors_ct)
    print(f"  {ct_h5ad}: AD donors={n_ad_donors}, Control donors={n_ct_donors}", flush=True)

# ── 3. Load DIU CSVs, filter, attach expression columns ───────────────────────
print("Building combined hits table...", flush=True)
chunks = []

for ct_csv in CT_MAP.values():
    fpath = DIU_DIR / f'DIU_result_with_Permutation_10000_donor_umi_cutoff_{ct_csv}.csv'
    df = pd.read_csv(fpath)
    df = df[df['permutation_significant']].copy()

    expr = expr_dict[ct_csv]
    df['total_counts_AD'] = df['transcript_name'].map(expr['total_counts_AD'])
    df['total_counts_CT'] = df['transcript_name'].map(expr['total_counts_CT'])
    df['mean_CPM_AD']     = df['transcript_name'].map(expr['mean_CPM_AD'])
    df['mean_CPM_CT']     = df['transcript_name'].map(expr['mean_CPM_CT'])
    df['log2FC']          = df['transcript_name'].map(expr['log2FC'])
    df.insert(0, 'cell_type', ct_csv)

    chunks.append(df)
    print(f"  {ct_csv}: {len(df)} hits", flush=True)

combined = pd.concat(chunks, ignore_index=True)

# Round float columns for readability
for col in ['mean_CPM_AD', 'mean_CPM_CT', 'log2FC']:
    combined[col] = combined[col].round(4)
combined['total_counts_AD'] = combined['total_counts_AD'].round(0).astype(int)
combined['total_counts_CT'] = combined['total_counts_CT'].round(0).astype(int)

# ── 4. Save ───────────────────────────────────────────────────────────────────
out = OUTPUT_DIR / 'DIU_significant_hits_combined.csv'
combined.to_csv(out, index=False)
print(f"\nSaved {len(combined)} rows × {len(combined.columns)} columns → {out}", flush=True)
print(f"Columns: {list(combined.columns)}", flush=True)
