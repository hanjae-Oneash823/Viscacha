"""
M07 — Expression Context (plan 5, M07).

Primary data: Layer 0 pseudobulk CSVs (per-donor counts). Supplementary: GTEx.

Per (gene, cell_type) candidate this module:
  1. Computes per-donor PSI for each significant transcript and its comparator,
     from the primary pseudobulk (AD/Control) and the sensitivity pseudobulk
     (Active control); summarises mean +/- SD per condition.
  2. Flags whether the comparator transcript is present in the pseudobulk
     (a MANE-Select comparator may be below the Layer 0 prevalence threshold).
  3. Classifies the expression pattern (pure switch / combined DE+DTU / complex)
     when Layer 1 gene-level DE is available; otherwise 'unavailable'.
  4. Computes the GTEx brain-specificity index for off-target context.

Condition labels come from config (Active control has a lowercase 'c', matching
Layer 0). Active control is non-fatal: absent sensitivity file -> two conditions.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from layer2.config import CONFIG
from layer2.inputs import CellTypeCandidate
from layer2.m01_sequence import M01Result
from layer2.utils.api_client import APIClient, APIError
from layer2.utils.audit import AuditLog


@dataclass
class ConditionPSI:
    condition: str
    mean: float
    sd: float
    n: int


@dataclass
class TranscriptPSI:
    transcript_name: str
    role: str                          # ad_enriched | control_enriched | comparator
    per_condition: list[ConditionPSI]
    donor_points: list[tuple] = field(default_factory=list)  # (condition, psi) per donor


@dataclass
class GTExProfile:
    gencode_id: str
    brain_mean_tpm: float
    all_mean_tpm: float
    bsi: float | None
    bsi_label: str
    cardiac_max_tpm: float
    liver_tpm: float


@dataclass
class M07Result:
    gene_id: str
    cell_type: str
    transcript_psi: list[TranscriptPSI] = field(default_factory=list)
    comparator_in_pseudobulk: bool = True
    expression_pattern: str = "unavailable"
    gtex: GTExProfile | None = None
    cell_type_breadth: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Per-donor PSI from pseudobulk
# ---------------------------------------------------------------------------

def _load_pb(pb_dir, cell_type: str, suffix: str = ""):
    cpath = pb_dir / f"counts_{cell_type}{suffix}.csv"
    mpath = pb_dir / f"metadata_{cell_type}{suffix}.csv"
    if not cpath.exists() or not mpath.exists():
        return None, None
    counts = pd.read_csv(cpath, index_col=0)     # donors x transcripts
    meta = pd.read_csv(mpath, index_col=0)
    return counts, meta


def _donor_psi(counts: pd.DataFrame, gene_id: str, tx: str) -> pd.Series:
    """PSI of transcript tx per donor = count(tx) / sum(counts of gene's tx)."""
    gene_cols = [c for c in counts.columns if c.startswith(gene_id + "-")]
    if tx not in counts.columns or not gene_cols:
        return pd.Series(dtype=float)
    gene_total = counts[gene_cols].sum(axis=1)
    return (counts[tx] / gene_total.replace(0, np.nan))


def _summarise(psi: pd.Series, meta: pd.DataFrame, condition: str) -> ConditionPSI:
    donors = [d for d in meta.index[meta["condition"] == condition] if d in psi.index]
    vals = psi.loc[donors].dropna()
    n = len(vals)
    mean = float(vals.mean()) if n else float("nan")
    sd = float(vals.std(ddof=1)) if n > 1 else 0.0
    return ConditionPSI(condition=condition, mean=mean, sd=sd, n=n)


# ---------------------------------------------------------------------------
# GTEx
# ---------------------------------------------------------------------------

def _gtex_profile(gene_symbol: str, client: APIClient, audit: AuditLog,
                  config) -> GTExProfile | None:
    base = config.apis["gtex"].rstrip("/")
    try:
        ref = client.get_json("gtex", f"{base}/reference/gene",
                              params={"geneId": gene_symbol})
        rows = ref.get("data", [])
        if not rows:
            return None
        gencode_id = rows[0]["gencodeId"]
        expr = client.get_json("gtex", f"{base}/expression/medianGeneExpression",
                               params={"gencodeId": gencode_id,
                                       "datasetId": "gtex_v8", "itemsPerPage": 100})
    except (APIError, KeyError) as e:
        audit.warning("gtex", f"{gene_symbol}: {e}")
        return None

    data = expr.get("data", [])
    if not data:
        return None
    brain, cardiac, liver, allt = [], [], [], []
    for x in data:
        tid = str(x.get("tissueSiteDetailId", ""))
        med = float(x.get("median", 0) or 0)
        allt.append(med)
        if tid.startswith("Brain"):
            brain.append(med)
        elif tid.startswith("Heart"):
            cardiac.append(med)
        elif tid == "Liver":
            liver.append(med)

    brain_mean = float(np.mean(brain)) if brain else 0.0
    all_mean = float(np.mean(allt)) if allt else 0.0
    bsi = round(brain_mean / all_mean, 2) if all_mean > 0 else None
    if bsi is None:
        label = "unavailable"
    elif bsi > config.thresholds["brain_specificity_specific"]:
        label = "brain-specific"
    elif bsi >= config.thresholds["brain_specificity_enriched"]:
        label = "brain-enriched"
    elif bsi >= 1.0:
        label = "ubiquitous"
    else:
        label = "peripherally-enriched"
    return GTExProfile(
        gencode_id=gencode_id, brain_mean_tpm=round(brain_mean, 2),
        all_mean_tpm=round(all_mean, 2), bsi=bsi, bsi_label=label,
        cardiac_max_tpm=round(max(cardiac), 2) if cardiac else 0.0,
        liver_tpm=round(liver[0], 2) if liver else 0.0)


# ---------------------------------------------------------------------------
# Expression pattern (needs Layer 1 step15 gene-level DE)
# ---------------------------------------------------------------------------

def _expression_pattern(gene_id: str, cell_type: str, config) -> str:
    """Three-way (plan M07):
      pure_switch     — no gene-level DE (|log2FC| < thr AND padj > thr): the
                        cleanest signal, total gene output constant.
      combined_de_dtu — significant gene-level DE (|log2FC| >= thr AND padj <= thr).
      de_uncertain    — large fold-change estimate but not significant (small n,
                        underpowered), or significant with small effect. Reported
                        honestly rather than collapsed into pure_switch.
    """
    path = config.gene_de_csv
    if not path.exists():
        return "unavailable"
    de = pd.read_csv(path)
    row = de[(de["gene_id"] == gene_id) & (de["cell_type"] == cell_type)]
    if row.empty:
        return "unavailable"
    lfc = float(row.iloc[0]["log2FC"])
    padj = row.iloc[0]["padj"]
    padj = float(padj) if pd.notna(padj) else 1.0
    big = abs(lfc) >= config.thresholds["gene_de_log2fc_threshold"]
    sig = padj <= config.thresholds["gene_de_padj_threshold"]
    if sig and big:
        return "combined_de_dtu"
    if not sig and not big:
        return "pure_switch"
    return "de_uncertain"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(cand: CellTypeCandidate, m01: M01Result, client: APIClient,
        audit: AuditLog, config=CONFIG) -> M07Result:
    ct = cand.cell_type
    cond = config.conditions
    res = M07Result(gene_id=cand.gene_id, cell_type=ct)

    counts, meta = _load_pb(config.pseudobulk_dir, ct)
    counts_s, meta_s = _load_pb(config.pseudobulk_dir, ct, "_sensitivity")
    if counts is None:
        res.warnings.append(f"no primary pseudobulk for {ct}")
        return res

    # Comparator transcript name for this gene (from M01)
    comparator_name = None
    for s in m01.significant:
        if comparator_name is None:
            comparator_name = s.comparator.transcript_name

    # transcripts to profile: significant ones + comparator
    sig_names = [t.transcript_id for t in cand.transcripts]
    roles = {t.transcript_id: t.role for t in cand.transcripts}
    profile_list = [(n, roles.get(n, "ad_enriched")) for n in sig_names]
    if comparator_name and comparator_name not in sig_names:
        profile_list.append((comparator_name, "comparator"))

    res.comparator_in_pseudobulk = (
        comparator_name in counts.columns if comparator_name else False)
    if comparator_name and not res.comparator_in_pseudobulk:
        res.warnings.append(
            f"comparator {comparator_name} below prevalence threshold — "
            f"not in pseudobulk; visualizations will degrade")

    for tname, role in profile_list:
        psi = _donor_psi(counts, cand.gene_id, tname)
        if psi.empty:
            continue
        per_cond = [_summarise(psi, meta, cond["control"]),
                    _summarise(psi, meta, cond["ad"])]
        points = [(meta.at[d, "condition"], float(psi[d]))
                  for d in psi.index if d in meta.index and not pd.isna(psi[d])]
        if counts_s is not None and meta_s is not None:
            psi_s = _donor_psi(counts_s, cand.gene_id, tname)
            if not psi_s.empty:
                per_cond.append(_summarise(psi_s, meta_s, cond["active_control"]))
                points += [(meta_s.at[d, "condition"], float(psi_s[d]))
                           for d in psi_s.index if d in meta_s.index and not pd.isna(psi_s[d])]
        res.transcript_psi.append(TranscriptPSI(tname, role, per_cond, points))

    res.expression_pattern = _expression_pattern(cand.gene_id, ct, config)
    res.gtex = _gtex_profile(cand.gene_id, client, audit, config)

    bsi = res.gtex.bsi if res.gtex else None
    audit.info(f"M07 {cand.gene_id}/{ct}: {len(res.transcript_psi)} tx profiled, "
               f"comparator_in_pb={res.comparator_in_pseudobulk}, "
               f"pattern={res.expression_pattern}, BSI={bsi}")
    return res


if __name__ == "__main__":
    from layer2.inputs import load_inputs
    from layer2 import m01_sequence as m01_mod
    audit = AuditLog(CONFIG.audit_log_path, verbose=False)
    client = APIClient(audit)
    inp = load_inputs()
    m01_cache = {}
    print(f"{'gene':9} {'cell_type':18} {'comp_in_pb':10} {'pattern':14} "
          f"{'BSI':>6} {'label':20} AD-enriched PSI (ctrl->AD)")
    print("-" * 120)
    for cand in inp.cell_type_work_list:
        gene = next(g for g in inp.gene_work_list if g.gene_id == cand.gene_id)
        if gene.gene_id not in m01_cache:
            m01_cache[gene.gene_id] = m01_mod.run(gene, client, audit, inp.enst_map)
        r = run(cand, m01_cache[gene.gene_id], client, audit)
        g = r.gtex
        ad_tx = next((t for t in r.transcript_psi if t.role == "ad_enriched"),
                     r.transcript_psi[0] if r.transcript_psi else None)
        psi_str = "-"
        if ad_tx:
            cm = {c.condition: c.mean for c in ad_tx.per_condition}
            psi_str = f"{ad_tx.transcript_name}: {cm.get('Control',float('nan')):.2f}->{cm.get('AD',float('nan')):.2f}"
        print(f"{r.gene_id:9} {r.cell_type:18} {str(r.comparator_in_pseudobulk):10} "
              f"{r.expression_pattern:14} {str(g.bsi if g else '-'):>6} "
              f"{(g.bsi_label if g else '-'):20} {psi_str}")
