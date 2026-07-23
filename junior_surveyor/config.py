"""JUNIOR_SURVEYOR — central configuration."""

from pathlib import Path

from assistant_surveyor.config import JUNIOR_PASS_BIOTYPES

REPO_ROOT = Path(__file__).resolve().parent.parent

HITS_CSV  = REPO_ROOT / "outputs/assistant_surveyor/hits_enriched.csv"
ALL_GROUPS_CSV = REPO_ROOT / "outputs/layer1_1/DIU_significant_hits_all_candidate_groups.csv"
OUT_DIR   = REPO_ROOT / "outputs/junior_surveyor"
CACHE_DIR = OUT_DIR / "cache"

# ---------------------------------------------------------------------------
# Ensembl REST
# ---------------------------------------------------------------------------
ENSEMBL_BASE    = "https://rest.ensembl.org"
ENSEMBL_BATCH   = 50      # max IDs per POST
ENSEMBL_TIMEOUT = 30
ENSEMBL_RETRIES = 4
ENSEMBL_BACKOFF = 4.0

# ---------------------------------------------------------------------------
# ChEMBL REST
# ---------------------------------------------------------------------------
CHEMBL_BASE    = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_TIMEOUT = 20
CHEMBL_RETRIES = 3
CHEMBL_BACKOFF = 3.0

# ---------------------------------------------------------------------------
# DGIdb GraphQL
# ---------------------------------------------------------------------------
DGIDB_ENDPOINT = "https://dgidb.org/api/graphql"
DGIDB_BATCH    = 50
DGIDB_TIMEOUT  = 30
DGIDB_RETRIES  = 3
DGIDB_BACKOFF  = 3.0

# ---------------------------------------------------------------------------
# Pfam / HMMER (local scan — replaces name-only UniProt domain guessing)
# ---------------------------------------------------------------------------
PFAM_HMM     = REPO_ROOT / "outputs/layer1_1/reference/pfam/Pfam-A.hmm"
HMMSCAN_BIN  = "/home/welcome3/anaconda3/envs/oneash_dtu/bin/hmmscan"
HMMPRESS_BIN = "/home/welcome3/anaconda3/envs/oneash_dtu/bin/hmmpress"
PFAM_CPU     = 4   # mirrors layer1_1/config.R's N_CORES

# ---------------------------------------------------------------------------
# TransDecoder (conda installs scripts under opt/transdecoder/util/)
# ---------------------------------------------------------------------------
_TD_UTIL = Path("/home/welcome3/anaconda3/envs/oneash_dtu/opt/transdecoder/util")
TRANSDECODER_LONGORFS = _TD_UTIL / "TransDecoder.LongOrfs"
TRANSDECODER_PREDICT  = _TD_UTIL / "TransDecoder.Predict"
GFFREAD_BIN           = "/home/welcome3/anaconda3/envs/oneash_dtu/bin/gffread"

# ---------------------------------------------------------------------------
# Next-stage selection (J4)
# ---------------------------------------------------------------------------
# protein_change_type values that mean "no real protein-level change"
NULL_PROTEIN_CHANGE_TYPES = {"identical", "no_sequence"}
