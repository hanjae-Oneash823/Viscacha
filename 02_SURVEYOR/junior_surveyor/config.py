"""JUNIOR_SURVEYOR — central configuration."""

from pathlib import Path

from assistant_surveyor.config import (
    JUNIOR_PASS_BIOTYPES, OT_ENDPOINT, OT_TIMEOUT, OT_RETRIES, OT_BACKOFF,
)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

HITS_CSV  = REPO_ROOT / "outputs/assistant_surveyor/hits_enriched.csv"
ALL_GROUPS_CSV = REPO_ROOT / "outputs/DIU_significant_hits/DIU_significant_hits_all_candidate_groups.csv"
OUT_DIR   = REPO_ROOT / "outputs/junior_surveyor"
CACHE_DIR = OUT_DIR / "cache"

TX_ID_MAP = REPO_ROOT / "outputs/annotation/tx_id_map.csv"
PSEUDOBULK_DIR = REPO_ROOT / "outputs/00_PreAggregation_QC/pseudobulk"

# SQANTI3 classification of the full long-read isoform catalog (known +
# novel), including its own ORF prediction per transcript (ORF_seq column).
# Keyed by "isoform" -- unversioned ENST accession for known transcripts,
# novel long-read ID (e.g. transcript1261.chr1.nnic) for novel ones -- the
# same identifier space as this pipeline's ENST_ID/alt_ENST_ID columns.
SQANTI_CLASSIFICATION_CSV = Path(
    "/node212data/welcome3/Grad_proj_2026/DATA/"
    "isoforms_classification_with_tx_name_and_gene_name.csv"
)

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
# Open Targets GraphQL — drugAndClinicalCandidates (third drug-evidence
# source alongside ChEMBL + DGIdb). Reuses assistant_surveyor's endpoint;
# batch size is smaller than L3's (200) since each target block here nests
# full drug + disease rows, not just a score, so the response is much heavier
# per target. Dropped 50->15 when each drug row grew to also carry
# mechanismsOfAction/indications/drugWarnings (full detail, not just name) --
# a 50-gene batch with several drugs each and dozens of indications per drug
# risked oversized/slow responses.
# ---------------------------------------------------------------------------
OT_DRUG_BATCH = 15

# ---------------------------------------------------------------------------
# Pharos/TCRD GraphQL — target development level (tdl) + ligand counts.
# Broadens recall beyond "has an approved/clinical drug" (ChEMBL/DGIdb/OT) to
# "has any known chemical ligand at all" (Tchem tier, tool/probe compounds
# included) -- verified via live schema introspection (2026-07): query root
# is `target(q: {sym: "..."})`, fields `tdl` and `ligandCounts { name value }`
# (rows named "ligand" = all known ligands, "drug" = the approved-drug subset).
# ---------------------------------------------------------------------------
PHAROS_ENDPOINT = "https://pharos-api.ncats.io/graphql"
PHAROS_BATCH    = 50
PHAROS_TIMEOUT  = 30
PHAROS_RETRIES  = 3
PHAROS_BACKOFF  = 3.0

# ---------------------------------------------------------------------------
# Pfam / HMMER (local scan — replaces name-only UniProt domain guessing)
# ---------------------------------------------------------------------------
PFAM_HMM     = REPO_ROOT / "outputs/reference/pfam/Pfam-A.hmm"
HMMSCAN_BIN  = "/home/welcome3/anaconda3/envs/oneash_dtu/bin/hmmscan"
HMMPRESS_BIN = "/home/welcome3/anaconda3/envs/oneash_dtu/bin/hmmpress"
PFAM_CPU     = 6

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

# ---------------------------------------------------------------------------
# Alt-row quality filter (J1c)
# ---------------------------------------------------------------------------
# trial_failure ranked-alternate rows below this AD usage share are dropped
# (new_target is untouched -- single alt_rank=0 row per hit, already a
# significant DTU hit by construction).
MIN_AD_USAGE_PCT = 0.03
