#!/usr/bin/env python
# ============================================================
# Layer 1 helper: gene-level expressing-cell counts per donor
# ============================================================
# The pseudobulk metadata's `n_cells` is the TOTAL number of cells of that
# broad cell type contributed by a donor (the pseudobulk denominator) -- the
# same value for every gene in that cell type. This script computes the more
# specific number per (cell_type, gene, donor): how many of those cells have
# nonzero expression of THIS gene (summed across all its isoform transcripts,
# using the same "GENE-" prefix grouping the R pipeline uses).
#
# Only computed for genes that are significant DTU hits (read from
# dtu_significant_all_celltypes.csv), since that's a short, known list and
# the source AnnData is otherwise too large to scan gene-by-gene.
#
# Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python 01_ViscachaDTU_Analysis/extract_gene_expressing_cells.py
#      (from anywhere -- paths below are absolute)
#      Must run AFTER run_layer1.R has written dtu_significant_all_celltypes.csv
# Output: outputs/01_ViscachaDTU_Analysis/gene_expressing_cells.csv
#   columns: cell_type, gene_id, donor, n_cells_expressing_gene, n_cells_total_donor
# ============================================================

import anndata as ad
import numpy as np
import pandas as pd

ADATA_PATH    = "/home/welcome3/Viscacha_pipeline/outputs/00_PreAggregation_QC/filtered_adata/adata_tx_step03.h5ad"
SIG_HITS_PATH = "/home/welcome3/Viscacha_pipeline/outputs/01_ViscachaDTU_Analysis/dtu_significant_all_celltypes.csv"
OUT_PATH      = "/home/welcome3/Viscacha_pipeline/outputs/01_ViscachaDTU_Analysis/gene_expressing_cells.csv"


def main():
    sig = pd.read_csv(SIG_HITS_PATH)
    pairs = sig[["cell_type", "gene_id"]].drop_duplicates()
    if pairs.empty:
        print("No significant hits found -- nothing to extract.")
        return

    adata = ad.read_h5ad(ADATA_PATH)
    # pseudobulk cell-type labels use underscores; AnnData obs uses spaces
    adata.obs["cell_type_us"] = adata.obs["cell_type"].str.replace(" ", "_")

    rows = []
    for ct in pairs["cell_type"].unique():
        sub = adata[adata.obs["cell_type_us"] == ct]
        genes = pairs.loc[pairs["cell_type"] == ct, "gene_id"].unique()
        for gn in genes:
            # Same prefix-grouping convention the R side uses for isoforms of a gene
            tx_cols = [v for v in sub.var_names if v.startswith(gn + "-")]
            if not tx_cols:
                print(f"  WARNING: no transcripts found for {ct}/{gn} -- skipping")
                continue

            gene_counts = np.asarray(sub[:, tx_cols].X.sum(axis=1)).ravel()
            donor = sub.obs["donor"].values

            df = pd.DataFrame({"donor": donor, "gene_count": gene_counts})
            grp = df.groupby("donor")["gene_count"].agg(
                n_cells_expressing_gene=lambda x: int((x > 0).sum()),
                n_cells_total_donor="size",
            ).reset_index()
            grp["cell_type"] = ct
            grp["gene_id"] = gn
            rows.append(grp)

    out = pd.concat(rows, ignore_index=True)
    out = out[["cell_type", "gene_id", "donor", "n_cells_expressing_gene", "n_cells_total_donor"]]
    out = out.sort_values(["cell_type", "gene_id", "donor"])
    out.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(out)} rows to {OUT_PATH}")


if __name__ == "__main__":
    main()
