"""
Surveyor (Layer 2) — configuration loader.

Loads surveyor_config.yaml, resolves all input/output paths relative to the
repository root, and expands ${ENV_VAR} references in the auth block.

Usage:
    from layer2.config import CONFIG
    CONFIG.dtu_significant_csv        # resolved absolute Path
    CONFIG.thresholds["chembl_potency_uM"]
    CONFIG.conditions["active_control"]  # "Active control"
"""

import os
import re
from pathlib import Path

import yaml

# Repository root = parent of the layer2/ directory holding this file.
REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = Path(__file__).resolve().parent / "surveyor_config.yaml"

_ENV_PATTERN = re.compile(r"\$\{([^}]+)\}")


def _expand_env(value):
    """Recursively expand ${VAR} references using os.environ (missing -> '')."""
    if isinstance(value, str):
        return _ENV_PATTERN.sub(lambda m: os.environ.get(m.group(1), ""), value)
    if isinstance(value, dict):
        return {k: _expand_env(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand_env(v) for v in value]
    return value


class SurveyorConfig:
    """Parsed configuration with resolved paths and convenient accessors."""

    def __init__(self, raw: dict):
        self._raw = raw

        # --- Sub-blocks (verbatim from YAML) ---
        self.version          = raw["version"]
        self.conditions       = raw["conditions"]
        self.apis             = raw["apis"]
        self.auth             = _expand_env(raw["auth"])
        self.thresholds       = raw["thresholds"]
        self.rdkit            = raw["rdkit"]
        self.tier_assignment  = raw["tier_assignment"]
        self.ad_pathways      = raw["ad_pathways"]
        self.ad_genes         = raw["ad_genes"]
        self.opentargets_labels = raw["opentargets_labels"]
        self.opentargets_efo  = raw["opentargets_efo"]
        self.braak_robustness = raw["braak_robustness"]
        self.cache            = raw["cache"]
        self.execution        = raw["execution"]

        # --- Resolved input paths ---
        inp = raw["inputs"]
        self.dtu_significant_csv = self._resolve(inp["dtu_significant_csv"])
        self.h5ad_path           = self._resolve(inp["h5ad_path"])
        self.pseudobulk_dir      = self._resolve(inp["pseudobulk_dir"])
        self.gene_de_csv         = self._resolve(inp["gene_de_csv"])

        # --- Resolved output paths ---
        out = raw["output"]
        self.output_base       = self._resolve(out["base_dir"])
        self.candidates_dir    = self._resolve(out["candidates_dir"])
        self.reports_dir       = self._resolve(out["reports_dir"])
        self.alphafold_dir     = self._resolve(out["alphafold_dir"])
        self.run_summary_path  = self._resolve(out["run_summary"])
        self.audit_log_path    = self._resolve(out["audit_log"])
        self.cache_dir         = self._resolve(raw["cache"]["directory"])

    @staticmethod
    def _resolve(rel: str) -> Path:
        """Resolve a repo-root-relative path string to an absolute Path."""
        p = Path(rel)
        return p if p.is_absolute() else (REPO_ROOT / p)

    def ensure_output_dirs(self):
        """Create all output directories (idempotent)."""
        for d in (self.output_base, self.candidates_dir, self.reports_dir,
                  self.alphafold_dir, self.cache_dir):
            d.mkdir(parents=True, exist_ok=True)

    @property
    def raw(self) -> dict:
        return self._raw


def load_config(path: Path = CONFIG_PATH) -> SurveyorConfig:
    with open(path) as f:
        raw = yaml.safe_load(f)
    return SurveyorConfig(raw)


# Module-level singleton — import this directly.
CONFIG = load_config()
