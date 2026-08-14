"""ASSISTANT_SURVEYOR — central configuration."""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

HITS_CSV   = REPO_ROOT / "outputs/DIU_significant_hits/DIU_significant_hits_initial_filter.csv"
TX_ID_MAP  = REPO_ROOT / "outputs/annotation/tx_id_map.csv"
ANNOT_TSV  = Path(
    "/node212data/welcome3/Grad_proj_2026/DATA/"
    "extended_annotation_including_refTSS_umi10_donor3_supported_"
    "novel_tx_with_gene_name_for_novel_tx.tsv"
)

OUT_DIR    = REPO_ROOT / "outputs/assistant_surveyor"
CACHE_DIR  = OUT_DIR / "cache"

# ---------------------------------------------------------------------------
# API endpoints & request settings
# ---------------------------------------------------------------------------
OT_ENDPOINT         = "https://api.platform.opentargets.org/api/v4/graphql"
OT_BATCH_SIZE       = 200
OT_TIMEOUT          = 45
OT_RETRIES          = 3
OT_BACKOFF          = 5.0

UNIPROT_ENDPOINT    = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_BATCH_SIZE  = 100
UNIPROT_TIMEOUT     = 30
UNIPROT_RETRIES     = 3
UNIPROT_BACKOFF     = 3.0

# ---------------------------------------------------------------------------
# OpenTargets disease EFO IDs
# ---------------------------------------------------------------------------
OT_EFO = {
    "AD":  "MONDO_0004975",
    "PD":  "MONDO_0005180",
    "FTD": "MONDO_0017276",
    "ALS": "MONDO_0004976",
}

# OT score label thresholds
OT_LABEL_SUPPORTED = 0.20
OT_LABEL_EMERGING  = 0.05

# ---------------------------------------------------------------------------
# Junior-layer gate — biotype_class values that pass a hit through to
# junior_surveyor. Everything else (RI, NMD, other) is dropped here.
# ---------------------------------------------------------------------------
JUNIOR_PASS_BIOTYPES: set[str] = {
    "PC_CDS", "PC_UTR", "PC_CDS_ND", "novel", "TEC",
}

# ---------------------------------------------------------------------------
# AD gene lists (L2)
# ---------------------------------------------------------------------------
AD_CAUSAL: set[str] = {
    "APP", "PSEN1", "PSEN2", "TREM2",
}

AD_GWAS: set[str] = {
    # Bellenguez et al. 2022 (Nature Genetics) + Lambert 2013 + meta-analyses
    "BIN1", "CLU", "ABCA7", "CR1", "PICALM", "MS4A4A", "MS4A6A",
    "EPHA1", "CD33", "CD2AP", "SORL1", "FERMT2", "HLA-DRB1", "INPP5D",
    "MEF2C", "CELF1", "ZCWPW1", "SLC24A4", "NME8", "ECHDC3", "APH1B",
    "CASS4", "PILRA", "ADAM10", "SHARPIN", "PLCG2", "SPI1", "SCIMP",
    "HAVCR2", "KAT8", "ACE", "APOE", "IQCK", "CLNK", "RIN3", "UMAD1",
    "HESX1", "ALPK2", "WDR12", "FAM193B", "PTK2B", "CELF2",
}

AD_PATHWAY: set[str] = {
    "MAPT", "GSK3B", "CDK5", "BACE1", "BACE2", "GBA", "LRRK2",
    "ADAM17", "APBB1", "APBB2", "FE65", "SORLA",
}
