"""MASTER_SURVEYOR — central configuration."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

HITS_CSV = REPO_ROOT / "outputs/junior_surveyor/hits_deep.csv"
OUT_DIR  = REPO_ROOT / "outputs/master_surveyor/plots"
