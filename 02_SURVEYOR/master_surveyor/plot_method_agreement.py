"""Visualize agreement between ColabFold and ESMFold for A-D gate hits.

For each alt isoform that has both models, writes a three-panel figure:
aligned C-alpha overlay, per-residue post-alignment displacement, and pLDDT
traces.  It also writes cohort summaries comparing agreement inside and
outside each spliced/altered span.  This is a method-agreement diagnostic;
it does not measure canonical-versus-alternate structural change.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import Superimposer

from master_surveyor import m2_structures
from master_surveyor.align_utils import canonical_span_to_alt, residue_correspondence
from master_surveyor.gate_matrix import GATES
from master_surveyor.m2b_structure_qc import _ca_atoms_by_residue, _plddt_by_residue
from master_surveyor.config import OUT_DIR, STRUCTURE_CACHE_DIR

GATE_PAIRS = OUT_DIR / "07_gate_pairs.csv"
OUT_SUBDIR = OUT_DIR / "11_method_agreement"
SUMMARY_CSV = OUT_SUBDIR / "method_agreement_summary.csv"

CF_COLOR = "#3977a8"
ESM_COLOR = "#e67e22"
ALTERED = "#c43c55"


def _models(sequence: str) -> tuple[Path | None, Path | None]:
    folder = STRUCTURE_CACHE_DIR / m2_structures.seq_hash(sequence)
    cf = sorted((folder / "alt_colabfold").glob("*_unrelaxed_rank_001_*.pdb"))
    esm = folder / "esmfold.pdb"
    return (cf[0] if cf else None, esm if esm.exists() else None)


def _aligned_arrays(cf_path: Path, esm_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return residue ids, aligned C-alpha coordinates, and pLDDT arrays."""
    cf_atoms = _ca_atoms_by_residue(cf_path)
    esm_atoms = _ca_atoms_by_residue(esm_path)
    ids = np.array(sorted(set(cf_atoms) & set(esm_atoms)), dtype=int)
    if len(ids) < 3:
        raise ValueError("fewer than three matched C-alpha residues")
    cf_list, esm_list = [cf_atoms[i] for i in ids], [esm_atoms[i] for i in ids]
    sup = Superimposer()
    sup.set_atoms(cf_list, esm_list)
    rot, tran = sup.rotran
    cf_xyz = np.array([a.get_coord() for a in cf_list])
    esm_xyz = np.array([a.get_coord() for a in esm_list]) @ rot + tran
    cf_p = _plddt_by_residue(cf_path)
    esm_p = _plddt_by_residue(esm_path)
    return ids, cf_xyz, esm_xyz, np.array([cf_p[i] for i in ids]), np.array([esm_p[i] for i in ids])


def _plot_hit(row: pd.Series, cf: Path, esm: Path, alt_start: int, alt_end: int) -> dict:
    ids, cf_xyz, esm_xyz, cf_p, esm_p = _aligned_arrays(cf, esm)
    dist = np.linalg.norm(cf_xyz - esm_xyz, axis=1)
    altered = (ids >= alt_start) & (ids <= alt_end)
    label = f"{row.gene_name} — {row.alt_transcript_name} ({row.cell_type})"
    fig = plt.figure(figsize=(12, 3.8), constrained_layout=True)
    ax0 = fig.add_subplot(1, 3, 1, projection="3d")
    stride = max(1, len(ids) // 1600)
    ax0.plot(*cf_xyz[::stride].T, lw=0.75, color=CF_COLOR, label="ColabFold rank 1")
    ax0.plot(*esm_xyz[::stride].T, lw=0.75, color=ESM_COLOR, alpha=.8, label="ESMFold")
    if altered.any():
        ax0.scatter(*cf_xyz[altered].T, s=5, color=ALTERED, label="altered span")
    ax0.set_title("Aligned Cα overlay", fontsize=10)
    ax0.set_axis_off(); ax0.legend(loc="upper left", fontsize=7)

    ax1 = fig.add_subplot(1, 3, 2)
    ax1.plot(ids, dist, color="#3f454c", lw=.8)
    ax1.axvspan(alt_start, alt_end, color=ALTERED, alpha=.16, label="altered span")
    ax1.set(title="Post-alignment displacement", xlabel="alternate residue", ylabel="Cα distance (Å)")
    ax1.legend(fontsize=7); ax1.spines[["top", "right"]].set_visible(False)

    ax2 = fig.add_subplot(1, 3, 3)
    ax2.plot(ids, cf_p, color=CF_COLOR, lw=.8, label="ColabFold")
    ax2.plot(ids, esm_p, color=ESM_COLOR, lw=.8, label="ESMFold")
    ax2.axvspan(alt_start, alt_end, color=ALTERED, alpha=.16)
    ax2.set(title="Per-residue confidence", xlabel="alternate residue", ylabel="pLDDT", ylim=(0, 100))
    ax2.legend(fontsize=7); ax2.spines[["top", "right"]].set_visible(False)
    fig.suptitle(label, fontsize=11, y=1.03)
    fig.savefig(OUT_SUBDIR / f"{row.gene_name}_{row.alt_transcript_name}_method_agreement.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    return {
        "gene_name": row.gene_name, "alt_transcript_name": row.alt_transcript_name,
        "cell_type": row.cell_type, "alt_span_start": alt_start, "alt_span_end": alt_end,
        "n_aligned_residues": len(ids), "global_rmsd_angstrom": float(np.sqrt(np.mean(dist**2))),
        "altered_mean_displacement_angstrom": float(np.mean(dist[altered])) if altered.any() else np.nan,
        "outside_mean_displacement_angstrom": float(np.mean(dist[~altered])) if (~altered).any() else np.nan,
        "altered_colabfold_plddt": float(np.mean(cf_p[altered])) if altered.any() else np.nan,
        "altered_esmfold_plddt": float(np.mean(esm_p[altered])) if altered.any() else np.nan,
        "status": "plotted",
    }


def _summaries(table: pd.DataFrame) -> None:
    plotted = table[table.status.eq("plotted")].copy()
    if plotted.empty:
        return
    fig, ax = plt.subplots(figsize=(5.2, 4.2), constrained_layout=True)
    size = 25 + 2 * plotted.n_aligned_residues.clip(upper=700)
    ax.scatter(plotted.altered_colabfold_plddt, plotted.altered_mean_displacement_angstrom,
               s=size, alpha=.75, color=CF_COLOR, edgecolor="white", linewidth=.6)
    for _, r in plotted.iterrows():
        ax.annotate(f"{r.gene_name}-{r.alt_transcript_name.split('-')[-1]}",
                    (r.altered_colabfold_plddt, r.altered_mean_displacement_angstrom), fontsize=7)
    ax.set(xlabel="ColabFold pLDDT in altered span", ylabel="CF–ESM mean displacement (Å)",
           title="Method agreement in altered regions")
    ax.axvline(70, ls="--", lw=.8, color="#888"); ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT_SUBDIR / "11a_method_agreement_confidence_scatter.png", dpi=220)
    plt.close(fig)

    values = [plotted.altered_mean_displacement_angstrom.dropna(), plotted.outside_mean_displacement_angstrom.dropna()]
    fig, ax = plt.subplots(figsize=(4.7, 4.2), constrained_layout=True)
    parts = ax.violinplot(values, showmedians=True)
    for body, color in zip(parts["bodies"], [ALTERED, "#6d7885"]): body.set_facecolor(color); body.set_alpha(.7)
    ax.set(xticks=[1, 2], xticklabels=["Altered span", "Outside span"], ylabel="CF–ESM mean displacement (Å)",
           title="Method disagreement by sequence region")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT_SUBDIR / "11b_method_agreement_altered_vs_outside.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUT_SUBDIR.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(GATE_PAIRS, low_memory=False)
    pairs = pairs[pairs[GATES].all(axis=1)].copy()
    records = []
    for _, row in pairs.iterrows():
        cf, esm = _models(row.alt_protein_seq)
        base = {"gene_name": row.gene_name, "alt_transcript_name": row.alt_transcript_name, "cell_type": row.cell_type}
        if cf is None or esm is None:
            records.append({**base, "status": "missing ESMFold" if esm is None else "missing ColabFold"})
            continue
        span = canonical_span_to_alt(residue_correspondence(row.canonical_protein_seq, row.alt_protein_seq), int(row.changed_aa_start), int(row.changed_aa_end))
        if span is None:
            records.append({**base, "status": "unmappable altered span"}); continue
        try:
            records.append(_plot_hit(row, cf, esm, *span))
        except Exception as exc:
            records.append({**base, "status": f"error: {exc}"})
    table = pd.DataFrame(records)
    table.to_csv(SUMMARY_CSV, index=False)
    _summaries(table)
    print(table[["gene_name", "alt_transcript_name", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
