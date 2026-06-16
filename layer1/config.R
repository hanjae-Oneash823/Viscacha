# ============================================================
# Viscacha Layer 1 — configuration
# ============================================================

# --- Directories ---
IN_DIR  <- "/home/welcome3/Viscacha_pipeline/outputs/layer0/pseudobulk"
OUT_DIR <- "/home/welcome3/Viscacha_pipeline/outputs/layer1"

# --- Cell types to analyze (Lymphocyte skipped: only 6 transcripts pass filter) ---
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

# --- Count-based filter thresholds (applied per cell type in step11) ---
MIN_TX_COUNT         <- 5    # transcript must have count >= this in >= MIN_SAMPS_FRAC of samples
MIN_GENE_COUNT       <- 10   # gene total must be >= this in >= MIN_SAMPS_FRAC of samples
MIN_SAMPS_FRAC       <- 0.30 # fraction of samples that must pass the count threshold

# --- DTU test parameters ---
ALPHA   <- 0.05
N_CORES <- 8

# --- Condition levels (Control is reference so conditionAD is the tested coefficient) ---
COND_LEVELS <- c("Control", "AD")
