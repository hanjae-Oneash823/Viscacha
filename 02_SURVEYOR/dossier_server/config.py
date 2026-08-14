"""DOSSIER_SERVER — central configuration."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

HOST = "0.0.0.0"   # server-hosted, reachable off this box -- single-user
                    # internal tool, no auth layer (see plan doc's Context)
PORT = 5057

BASE_DIR   = REPO_ROOT / "outputs/master_surveyor"
CART_JSON  = BASE_DIR / "cart.json"
JOBS_DIR   = BASE_DIR / "jobs"
DOCKING_DIR = BASE_DIR / "docking"
