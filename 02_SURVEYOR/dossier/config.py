"""DOSSIER — central configuration.

Standalone per-candidate report generator. Reads outputs already produced by
junior_surveyor / 01_ViscachaDTU_Analysis; writes nothing back into those
pipelines. Not gated by master_surveyor's TF+DR-only scope -- any (gene,
cell_type) present in hits_deep.csv can be rendered, including
novel_target_candidate hits.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

HITS_CSV      = REPO_ROOT / "outputs/junior_surveyor/hits_deep.csv"
TX_ID_MAP_CSV = REPO_ROOT / "outputs/annotation/tx_id_map.csv"
PFAM_CACHE    = REPO_ROOT / "outputs/junior_surveyor/cache/j2_pfam_hits.json"
PSEUDOBULK_DIR = REPO_ROOT / "outputs/00_PreAggregation_QC/pseudobulk"

OUT_DIR   = REPO_ROOT / "outputs/dossier"
ASSETS_DIR = Path(__file__).resolve().parent / "assets"

# Human-readable label for each master_group value, shown as the summary badge.
CANDIDATE_TYPE_LABELS = {
    "trial_failure_candidate":   "may explain drug trial failure",
    "drug_repurposing_candidate": "may be a target for drug repurposing",
    "novel_target_candidate":    "novel drug target",
}

CELL_TYPES = [
    "Astrocyte", "Excitatory_neuron", "Inhibitory_neuron",
    "Microglia", "Oligodendrocyte", "OPC", "Vascular_cell", "Lymphocyte",
]

# master_surveyor's own scope (docs/MASTER_SURVEYOR_plan.md) -- novel_target_candidate
# is explicitly excluded there, so every batch/index tool matches it here.
MASTER_SURVEYOR_GROUPS = ["trial_failure_candidate", "drug_repurposing_candidate"]

# ---------------------------------------------------------------------------
# Canonical structure images (fetch_structures.py) -- for visual purposes
# only in the dossier header, not used for any docking/structural analysis.
# ---------------------------------------------------------------------------
STRUCTURE_CACHE_DIR = OUT_DIR / "cache" / "structures"
AFDB_API      = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
AFDB_TIMEOUT  = 20
CHROME_BIN    = "/home/welcome3/.cache/puppeteer/chrome/linux-150.0.7871.24/chrome-linux64/chrome"
