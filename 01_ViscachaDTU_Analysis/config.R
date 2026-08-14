# ============================================================
# Viscacha Layer 1 — configuration
# ============================================================

# --- Directories ---
IN_DIR  <- "/home/welcome3/Viscacha_pipeline/outputs/00_PreAggregation_QC/pseudobulk"
OUT_DIR <- "/home/welcome3/Viscacha_pipeline/outputs/01_ViscachaDTU_Analysis"

# --- Cell types to analyze (Lymphocyte skipped: only 2 donors pass MIN_CELLS_PB) ---
CELL_TYPES <- c(
  "Excitatory_neuron",
  "Inhibitory_neuron",
  "Oligodendrocyte",
  "OPC",
  "Astrocyte",
  "Microglia",
  "Vascular_cell"
)

# --- Model formulas (Option C) ---
# Primary:     main results, no braak_stage (avoids r=0.914 collinearity with condition)
# Braak:       sensitivity check — same model + braak_stage
FORMULA_PRIMARY <- ~ condition + age + sex + median_pct_mt
FORMULA_BRAAK   <- ~ condition + age + sex + median_pct_mt + braak_stage

# --- Fit-set expression filter (step11's filter_counts_for_fit) ---
# Restricts which transcripts are fit/tested by satuRn+stageR, to keep the
# multiple-testing universe sized to testable signal. PSI denominators are
# unaffected — those always use the full structural set from filter_counts().
MIN_TX_COUNT   <- 5    # transcript must have count >= this in >= MIN_SAMPS_FRAC of samples
MIN_SAMPS_FRAC <- 0.30 # fraction of samples that must pass the count threshold

# --- DTU test parameters ---
ALPHA   <- 0.05
N_CORES <- 4

# --- Condition levels (Control is reference so conditionAD is the tested coefficient) ---
COND_LEVELS <- c("Control", "AD")
