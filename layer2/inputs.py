"""
Surveyor (Layer 2) — input contract, validation, and candidate assembly.

Implements plan sections 2.1-2.6 and the candidate-table / gene-deduplication
parts of section 4.1 (orchestrator steps 1-3).

This module is pure-local: it reads only files that Layer 0 / Layer 1 already
produce. No external API calls, no biopython, no rdkit. It is the foundation
the rest of Surveyor hangs off, and is independently runnable for verification:

    python -m layer2.inputs        # prints a summary of the assembled inputs
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from layer2.config import CONFIG

# Columns Surveyor requires from dtu_significant_all_celltypes.csv (plan 2.1).
REQUIRED_DTU_COLUMNS = [
    "transcript_id", "gene_id", "padj_gene", "padj_tx",
    "psi_AD", "psi_ctrl", "delta_psi", "psi_active_ctrl",
    "robust_to_braak", "cell_type",
]


class InputValidationError(Exception):
    """Raised on a blocking validation failure (plan 2.6). Aborts the run."""


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SignificantTranscript:
    """One significant transcript in one cell type (one row of the DTU CSV)."""
    transcript_id: str          # GENCODE name, e.g. "CAMK2B-204"
    gene_id: str                # gene symbol, e.g. "CAMK2B"
    cell_type: str              # e.g. "Excitatory_neuron"
    enst_id: str                # Ensembl transcript ID (from h5ad var)
    ensg_id: str                # Ensembl gene ID (from h5ad var)
    delta_psi: float
    psi_AD: float
    psi_ctrl: float
    psi_active_ctrl: float
    padj_gene: float
    padj_tx: float
    robust_to_braak: bool

    @property
    def role(self) -> str:
        """ad_enriched (ΔPSI > 0) or control_enriched (ΔPSI < 0). Plan M01."""
        return "ad_enriched" if self.delta_psi > 0 else "control_enriched"


@dataclass
class GeneCandidate:
    """Gene-level unit: all significant transcripts for one gene, across every
    cell type. Drives the gene-level modules (M01, M02, M03, M04, M05, M06),
    which run once per gene regardless of cell-type multiplicity (plan 4.1)."""
    gene_id: str
    ensg_id: str
    transcripts: list[SignificantTranscript] = field(default_factory=list)

    @property
    def cell_types(self) -> list[str]:
        return sorted({t.cell_type for t in self.transcripts})

    @property
    def unique_transcript_ids(self) -> list[str]:
        # de-duplicated, preserving a stable order
        seen, out = set(), []
        for t in self.transcripts:
            if t.transcript_id not in seen:
                seen.add(t.transcript_id)
                out.append(t.transcript_id)
        return out


@dataclass
class CellTypeCandidate:
    """(gene, cell_type) unit: the significant transcripts for one gene in one
    cell type. Drives the cell-type-level modules (M07, M_VIS, M08)."""
    gene_id: str
    cell_type: str
    transcripts: list[SignificantTranscript] = field(default_factory=list)


@dataclass
class SurveyorInputs:
    """The fully validated, assembled input bundle for one Surveyor run."""
    gene_work_list: list[GeneCandidate]
    cell_type_work_list: list[CellTypeCandidate]
    braak_condition_correlation: float | None
    input_file_checksums: dict[str, str]
    enst_map: "pd.DataFrame"          # transcript_name index -> ENST_ID, ENSG_ID
    warnings: list[str] = field(default_factory=list)

    @property
    def n_genes(self) -> int:
        return len(self.gene_work_list)

    @property
    def n_candidates(self) -> int:
        return len(self.cell_type_work_list)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _md5(path: Path) -> str:
    if not path.exists():
        return "absent"
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_enst_map(h5ad_path: Path) -> pd.DataFrame:
    """Read only adata.var (transcript_name index -> ENST_ID, ENSG_ID).
    Opened backed='r' so the expression matrix is never loaded (plan 2.2)."""
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        var = adata.var.copy()
    finally:
        if adata.isbacked:
            adata.file.close()
    if "ENST_ID" not in var.columns:
        raise InputValidationError(
            f"h5ad var has no 'ENST_ID' column (found: {list(var.columns)})"
        )
    return var


def _compute_braak_correlation(pseudobulk_dir: Path, cell_types: list[str],
                               cond_ad: str, cond_ctrl: str) -> tuple[float | None, str | None]:
    """Spearman r between numeric Braak stage and binary condition (AD=1,
    Control=0), over the donors actually used in the DTU test — i.e. the union
    of the PRIMARY pseudobulk metadata across the analysed cell types,
    deduplicated by donor (plan 4.1, gap-2 fix). Active Control donors live only
    in the sensitivity files and are correctly excluded.

    Returns (r, note). r is None if it cannot be computed."""
    frames = []
    for ct in cell_types:
        meta_path = pseudobulk_dir / f"metadata_{ct}.csv"
        if meta_path.exists():
            frames.append(pd.read_csv(meta_path, index_col=0))
    if not frames:
        return None, "no primary metadata files found"

    meta = pd.concat(frames)
    meta = meta[~meta.index.duplicated(keep="first")]
    meta = meta[meta["condition"].isin([cond_ad, cond_ctrl])]
    meta = meta.dropna(subset=["braak_stage", "condition"])
    if meta["condition"].nunique() < 2 or len(meta) < 3:
        return None, "insufficient donors with both Braak and condition"

    cond_binary = (meta["condition"] == cond_ad).astype(int).to_numpy()
    braak = meta["braak_stage"].astype(float).to_numpy()
    if np.std(braak) == 0:
        return None, "Braak stage has zero variance"

    r, _ = spearmanr(braak, cond_binary)
    return float(r), None


# ---------------------------------------------------------------------------
# Validation (plan 2.6)
# ---------------------------------------------------------------------------

def _validate_and_load_dtu(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise InputValidationError(f"DTU significant CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    missing = [c for c in REQUIRED_DTU_COLUMNS if c not in df.columns]
    if missing:
        raise InputValidationError(
            f"DTU CSV missing required columns: {missing}\n"
            f"  found: {list(df.columns)}"
        )
    if len(df) == 0:
        raise InputValidationError("DTU CSV has zero rows — nothing to process.")
    return df


# ---------------------------------------------------------------------------
# Assembly (orchestrator steps 1-3)
# ---------------------------------------------------------------------------

def load_inputs(config=CONFIG) -> SurveyorInputs:
    """Validate inputs, build the gene/cell-type work lists, compute the Braak
    correlation, and record input file checksums. Raises InputValidationError
    on any blocking failure (plan 2.6)."""
    warnings: list[str] = []

    # Step 1 — validation + load DTU table
    df = _validate_and_load_dtu(config.dtu_significant_csv)

    # ENST mapping (blocking if h5ad missing/malformed)
    if not config.h5ad_path.exists():
        raise InputValidationError(f"h5ad not found: {config.h5ad_path}")
    enst_map = _load_enst_map(config.h5ad_path)

    cond = config.conditions

    # Build SignificantTranscript objects, attaching ENST/ENSG from the h5ad var
    transcripts: list[SignificantTranscript] = []
    for _, row in df.iterrows():
        tx = row["transcript_id"]
        if tx in enst_map.index:
            enst = str(enst_map.at[tx, "ENST_ID"])
            ensg = str(enst_map.at[tx, "ENSG_ID"]) if "ENSG_ID" in enst_map.columns else ""
        else:
            enst, ensg = "", ""
            warnings.append(f"transcript '{tx}' not found in h5ad var — no ENST ID")

        transcripts.append(SignificantTranscript(
            transcript_id=tx,
            gene_id=str(row["gene_id"]),
            cell_type=str(row["cell_type"]),
            enst_id=enst,
            ensg_id=ensg,
            delta_psi=float(row["delta_psi"]),
            psi_AD=float(row["psi_AD"]),
            psi_ctrl=float(row["psi_ctrl"]),
            psi_active_ctrl=(float(row["psi_active_ctrl"])
                             if pd.notna(row["psi_active_ctrl"]) else float("nan")),
            padj_gene=float(row["padj_gene"]),
            padj_tx=float(row["padj_tx"]),
            robust_to_braak=bool(row["robust_to_braak"]),
        ))

    # Step 3 — gene deduplication into two work lists
    gene_map: dict[str, GeneCandidate] = {}
    ct_map: dict[tuple[str, str], CellTypeCandidate] = {}
    for t in transcripts:
        g = gene_map.get(t.gene_id)
        if g is None:
            g = GeneCandidate(gene_id=t.gene_id, ensg_id=t.ensg_id)
            gene_map[t.gene_id] = g
        g.transcripts.append(t)
        if not g.ensg_id and t.ensg_id:
            g.ensg_id = t.ensg_id

        key = (t.gene_id, t.cell_type)
        c = ct_map.get(key)
        if c is None:
            c = CellTypeCandidate(gene_id=t.gene_id, cell_type=t.cell_type)
            ct_map[key] = c
        c.transcripts.append(t)

    gene_work_list = sorted(gene_map.values(), key=lambda g: g.gene_id)
    cell_type_work_list = sorted(ct_map.values(), key=lambda c: (c.gene_id, c.cell_type))

    # Step 2 — Braak-condition correlation (over analysed cell types)
    analysed_cts = sorted({t.cell_type for t in transcripts})
    braak_r, braak_note = _compute_braak_correlation(
        config.pseudobulk_dir, analysed_cts, cond["ad"], cond["control"]
    )
    if braak_note:
        warnings.append(f"Braak correlation: {braak_note}")

    # Non-blocking input warnings (plan 2.6)
    if not config.gene_de_csv.exists():
        warnings.append("gene_level_de_results.csv absent — expression pattern "
                        "classification will be 'unavailable' (Layer 1 step15 not run)")

    checksums = {
        "dtu_csv_md5":     _md5(config.dtu_significant_csv),
        "h5ad_md5":        _md5(config.h5ad_path),
        "gene_de_csv_md5": _md5(config.gene_de_csv),
    }

    return SurveyorInputs(
        gene_work_list=gene_work_list,
        cell_type_work_list=cell_type_work_list,
        braak_condition_correlation=braak_r,
        input_file_checksums=checksums,
        enst_map=enst_map,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Standalone verification
# ---------------------------------------------------------------------------

def _summary(inp: SurveyorInputs) -> str:
    lines = [
        "=" * 64,
        "Surveyor input assembly — summary",
        "=" * 64,
        f"Genes (gene_work_list):            {inp.n_genes}",
        f"Candidates (cell_type_work_list):  {inp.n_candidates}",
        f"Braak-condition Spearman r:        "
        f"{inp.braak_condition_correlation:.4f}" if inp.braak_condition_correlation is not None
        else "Braak-condition Spearman r:        unavailable",
        "",
        "Gene work list:",
    ]
    for g in inp.gene_work_list:
        cts = ", ".join(g.cell_types)
        lines.append(f"  {g.gene_id:10} ({g.ensg_id or 'no-ENSG':16})  "
                     f"{len(g.unique_transcript_ids)} tx  | cell types: {cts}")
    lines.append("")
    lines.append("Candidate (gene x cell_type) work list:")
    for c in inp.cell_type_work_list:
        txs = ", ".join(f"{t.transcript_id}[{t.role[:3]}]" for t in c.transcripts)
        lines.append(f"  {c.gene_id:10} | {c.cell_type:18} | {txs}")
    if inp.warnings:
        lines.append("")
        lines.append("Warnings:")
        for w in inp.warnings:
            lines.append(f"  - {w}")
    lines.append("")
    lines.append("Input checksums:")
    for k, v in inp.input_file_checksums.items():
        lines.append(f"  {k:18} {v}")
    lines.append("=" * 64)
    return "\n".join(lines)


if __name__ == "__main__":
    inp = load_inputs()
    print(_summary(inp))
