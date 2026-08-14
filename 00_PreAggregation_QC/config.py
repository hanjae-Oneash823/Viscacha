import glob
from pathlib import Path

def _glob_one(pattern: str) -> Path:
    """Resolve a glob pattern to a single Path; raises if not found."""
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No file matched: {pattern}")
    return Path(matches[0])

# --- Directories ---
DATA_DIR  = Path("/node212data/welcome3/Grad_proj_2026/DATA")
META_DIR  = Path("/node212data/welcome3/Grad_proj_2026/sample_patient_metadata")
OUT_DIR   = Path("/home/welcome3/Viscacha_pipeline/outputs/00_PreAggregation_QC")

# --- Input files ---
SR_PATH   = DATA_DIR / "adata_sr.h5ad"
TX_PATH   = DATA_DIR / "adata_transcript_loose_filtering_for_bulk_analysis.h5ad"
GX_PATH   = DATA_DIR / "adata_gene_loose_filtering_for_bulk_analysis.h5ad"
# Korean filenames resolved via glob to avoid NFC/NFD encoding mismatch
SMC_XLSX  = _glob_one(str(META_DIR / "*.xlsx"))
PO_PPTX   = _glob_one(str(META_DIR / "*.pptx"))

# --- Output directories ---
OUT_META    = OUT_DIR / "metadata"
OUT_PB      = OUT_DIR / "pseudobulk"
OUT_ADATA   = OUT_DIR / "filtered_adata"
OUT_QC      = OUT_DIR / "qc_logs"
OUT_PLOTS   = OUT_DIR / "plots"

# --- Fields transferred from adata_sr.obs into adata_tx.obs ---
# cell_type is already present in adata_tx; only these two need transfer
TRANSFER_COLS = ["doublet_score", "pct_counts_mt"]

# --- Unified metadata output columns ---
UNIFIED_META_COLS = [
    "donor_id", "condition", "age", "sex",
    "braak_stage", "thal_phase", "cerad_score", "apoe", "adnc_grade",
]

# --- Cell types (must match adata.obs['cell_type'] values exactly — spaces, not underscores) ---
CELL_TYPES = [
    "Excitatory neuron",
    "Inhibitory neuron",
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "Microglia",
    "Vascular cell",
    "Lymphocyte",
]

# --- Thresholds ---
DOUBLET_THRESHOLD   = 0.3
MIN_PREVALENCE      = 0.30
MIN_CELLS_PB        = 10
DROPOUT_FLAG_PCT    = 0.20   # flag donor×cell_type losing >20% barcodes

# --- Condition labels (canonical strings matching AnnData) ---
COND_AD      = "AD"
COND_CTRL    = "Control"
COND_ACTIVE  = "Active control"   # lowercase 'c' — must match AnnData exactly

# --- Braak stage: B-score tier → numeric midpoint ---
B_TO_BRAAK = {"B0": 0, "B1": 2, "B2": 4, "B3": 6}

# --- Roman numeral → integer (for pptx Thal phase and Braak stage) ---
ROMAN_TO_INT = {"0": 0, "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}

# --- CERAD text → ordinal ---
CERAD_TEXT_TO_INT = {"none": 0, "sparse": 1, "moderate": 2, "frequent": 3}
