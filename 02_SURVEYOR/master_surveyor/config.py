"""MASTER_SURVEYOR — central configuration."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

HITS_CSV = REPO_ROOT / "outputs/junior_surveyor/hits_deep.csv"
OUT_DIR  = REPO_ROOT / "outputs/master_surveyor/plots"

# ---------------------------------------------------------------------------
# Scope (docs/MASTER_SURVEYOR_plan.md) -- novel_target_candidate is excluded.
# ---------------------------------------------------------------------------
MASTER_SURVEYOR_GROUPS = ["trial_failure_candidate", "drug_repurposing_candidate"]

BASE_DIR        = REPO_ROOT / "outputs/master_surveyor"
CACHE_DIR       = BASE_DIR / "cache"
STRUCTURE_CACHE_DIR = CACHE_DIR / "structures"   # <seq_hash>/ subdirs
FASTA_CACHE_DIR = CACHE_DIR / "fasta"
JOBS_DIR        = BASE_DIR / "jobs"
DOCKING_DIR     = BASE_DIR / "docking"
SHORTLIST_CSV   = BASE_DIR / "shortlist.csv"
CART_JSON       = BASE_DIR / "cart.json"

# ---------------------------------------------------------------------------
# m0_select — deliberately permissive default thresholds. Passed as a dict
# (not read as fixed module constants) so dossier_server can override them
# per-request from live UI controls without any code change.
# ---------------------------------------------------------------------------
DEFAULT_THRESHOLDS = {
    "tf_min_abs_delta_usage":    0.0,     # trial_failure: |delta_usage| floor
    "tf_require_domain_overlap": False,   # trial_failure: affected_domain != "none"
    "dr_min_chembl_or_ot_phase": 0,       # drug_repurposing: clinical-phase floor
    "dr_require_structural_change": False,  # drug_repurposing: exclude "substitution"-only
    "min_structure_confidence": None,     # None = no filter; else "low"/"medium"/"high" floor
}

# ---------------------------------------------------------------------------
# ChEMBL REST (ligand SMILES lookup) -- same endpoint junior_surveyor's J3
# already uses for target/drug lookups.
# ---------------------------------------------------------------------------
CHEMBL_BASE    = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_TIMEOUT = 20
CHEMBL_RETRIES = 3
CHEMBL_BACKOFF = 3.0

# ---------------------------------------------------------------------------
# AlphaFold DB REST (canonical-sequence structure lookup)
# ---------------------------------------------------------------------------
AFDB_API     = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
AFDB_TIMEOUT = 20

# ---------------------------------------------------------------------------
# ColabFold (local, GPU-pinned) -- installed under /home/welcome3/tools per
# docs/MASTER_SURVEYOR_plan.md's infra notes (root disk too small for the
# AF2 params + conda env; /home has the space).
# ---------------------------------------------------------------------------
COLABFOLD_BIN       = "/home/welcome3/tools/localcolabfold/colabfold-conda/bin/colabfold_batch"
CUDA_VISIBLE_DEVICES = "1"

# Canonical-sequence fallback (AFDB miss): standard settings, templates ON --
# this is the well-studied fold, a template hit is expected/desirable.
COLABFOLD_CANONICAL_ARGS: list[str] = []

# Alt/isoform sequence: ALWAYS run this way, never AFDB. Templates OFF (a
# template hit here would just be the canonical protein's own structure,
# biasing the model toward the fold this analysis exists to question) and
# recycles well above the AF2/AFDB default of 3, following the published
# no-template + elevated-recycle protocol used for AD-variant (Presenilin-1)
# structure prediction. --num-seeds gives a small ensemble so model-to-model
# spread over the altered region is itself a confidence signal (see
# m2b_structure_qc.py), not just the top-ranked model.
COLABFOLD_ALT_NUM_RECYCLE = 20
COLABFOLD_ALT_NUM_SEEDS   = 5
COLABFOLD_ALT_ARGS = [
    "--num-recycle", str(COLABFOLD_ALT_NUM_RECYCLE),
    "--num-seeds", str(COLABFOLD_ALT_NUM_SEEDS),
    "--num-models", "1",
]
# --num-models 1 pins to a single AF2-ptm checkpoint (model_1) so
# num_seeds x num_models = exactly 5 predictions, one per seed -- verified
# empirically that colabfold_batch ranks ALL model x seed combinations
# together into one global rank_001..rank_NNN list (NOT a separate rank_001
# per seed), so without this, --num-models's default of 5 would produce 25
# predictions with no clean "one per seed" file to glob for the ensemble.
# --templates is opt-in in colabfold_batch (omitting the flag means no
# templates are used), so COLABFOLD_ALT_ARGS deliberately does NOT include it.

COLABFOLD_TIMEOUT_S = 14400   # one alt-sequence ensemble prediction, wall clock

# ---------------------------------------------------------------------------
# ESMFold (local, GPU) -- HuggingFace transformers port (facebook/esmfold_v1),
# not fair-esm+openfold: the HF port vendors the folding trunk in pure
# PyTorch (transformers.models.esm.openfold_utils) so it needs no separate
# openfold pip package or CUDA-kernel compilation step, which matters a lot
# for install reliability on a shared box. Used as an MSA-free cross-check
# for the isoform-altered span specifically, since MMseqs2 pulls homologs of
# the CANONICAL sequence and a truncated/frameshifted/novel alt region can
# have thin-to-zero aligned homology -- exactly where AF2's co-evolutionary
# signal runs out.
# ---------------------------------------------------------------------------
ESMFOLD_ENV_PYTHON = "/home/welcome3/anaconda3/envs/esmfold/bin/python"
ESMFOLD_MODEL_NAME  = "facebook/esmfold_v1"
ESMFOLD_CHUNK_SIZE  = 128   # trades speed for memory if a sequence OOMs; None = full

# Warm server (scripts/esmfold_server.py) -- keeps the model resident on GPU
# instead of reloading ~2GB per call (measured ~45s of pure overhead per
# call via the subprocess path). m2_structures.py tries this first and
# falls back to the subprocess script if it's not running, so starting this
# server is an optional speedup, not a hard dependency.
ESMFOLD_SERVER_URL = "http://127.0.0.1:5058"
ESMFOLD_SERVER_HEALTHCHECK_TIMEOUT_S = 2
ESMFOLD_SERVER_FOLD_TIMEOUT_S = 1800
