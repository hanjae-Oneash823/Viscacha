"""
M08 — Report Assembly (plan 5, M08; v2.2).

Assembles, for one (gene, cell_type) candidate, the machine-readable JSON
dossier and the self-contained HTML report from every module's output. Also
computes the docking-readiness verdict (the v2.2 priority table, Type N first)
and writes the AlphaFold submission package for docking-relevant verdicts.

Inputs are bundled in a CandidateContext (assembled by the orchestrator) so this
module has no API access of its own.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path

from layer2.config import CONFIG
from layer2.inputs import CellTypeCandidate
from layer2.m01_sequence import M01Result
from layer2.m02_domain import M02Result
from layer2.m03_splice import M03Result, SpliceClassification
from layer2.m04_drug import M04Result
from layer2.m05_disease import DiseaseAssociation
from layer2.m06_pathway import PathwayResult
from layer2.m07_expression import M07Result

# Strength order for picking a candidate-level type from its transcripts.
_TYPE_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "N": 4}


@dataclass
class CandidateContext:
    """Everything needed to assemble one (gene, cell_type) dossier."""
    candidate: CellTypeCandidate
    m01: M01Result
    m02: M02Result
    m03: M03Result
    m04: M04Result
    m05: DiseaseAssociation
    m06: PathwayResult
    m07: M07Result
    braak_correlation: float | None
    input_checksums: dict
    database_versions: dict
    viz_svgs: dict = field(default_factory=dict)   # name -> inline SVG (from M_VIS)


# ---------------------------------------------------------------------------
# Verdict (plan v2.2 priority table)
# ---------------------------------------------------------------------------

def _candidate_type(classifications: list[SpliceClassification]) -> str:
    if not classifications:
        return "N"
    return min((c.splice_type for c in classifications), key=lambda t: _TYPE_ORDER[t])


def compute_verdict(ctx: CandidateContext) -> tuple[str, str, str]:
    """Return (candidate_type, verdict, reason)."""
    cls = [c for c in ctx.m03.classifications
           if c.transcript_name in {t.transcript_id for t in ctx.candidate.transcripts}]
    ctype = _candidate_type(cls)
    coding = any(s.transcript.is_protein_coding for s in ctx.m01.significant)
    has_drug = ctx.m04.drug_target_status == "known"
    bsi = ctx.m07.gtex.bsi if ctx.m07.gtex else None

    # 1. Type N — checked first
    if ctype == "N":
        return ctype, "no-go", "no protein-level change — regulatory candidate (no docking)"
    # 2. M02 unavailable
    if not ctx.m02.available:
        return ctype, "conditional", "domain annotation unavailable — manual splice classification required"
    # 3-6. Type A/B
    if ctype in ("A", "B"):
        if coding and has_drug:
            return ctype, "proceed", "pocket-altering splice with a known drug"
        if coding and (bsi is None or bsi >= 1.0):
            return ctype, "conditional", "pocket-altering splice; no known drug (novel target)"
        if coding and bsi is not None and bsi < 1.0:
            return ctype, "conditional", "pocket-altering splice; peripherally enriched (off-target risk)"
        return ctype, "conditional", "confirm protein-coding status"
    # 7. Type C
    if ctype == "C":
        return ctype, "conditional", "domain-level structural change — indirect/allosteric docking signal"
    # 8-10. Type D
    if bsi is not None and bsi < 1.0:
        return ctype, "no-go", "no docking hypothesis; flag for PPI surface analysis"
    return ctype, "conditional", "splice in disordered/unannotated region — PPI surface analysis recommended"


# ---------------------------------------------------------------------------
# JSON dossier
# ---------------------------------------------------------------------------

def build_dossier(ctx: CandidateContext, ctype: str, verdict: str, reason: str) -> dict:
    cand = ctx.candidate
    cls_by_tx = {c.transcript_name: c for c in ctx.m03.classifications}

    sig = []
    for t in cand.transcripts:
        c = cls_by_tx.get(t.transcript_id)
        sig.append({
            "enst_id": t.enst_id, "transcript_name": t.transcript_id,
            "delta_psi": t.delta_psi, "role": t.role,
            "splice_type": c.splice_type if c else None,
            "tier": c.tier if c else None,
            "narrative": c.narrative if c else None,
        })

    # comparator (from M01, first significant unit)
    comp = None
    if ctx.m01.significant:
        s0 = ctx.m01.significant[0]
        comp = {"enst_id": s0.comparator.enst_id,
                "transcript_name": s0.comparator.transcript_name,
                "source": s0.comparator_source}

    g = ctx.m07.gtex
    return {
        "gene": cand.gene_id,
        "cell_type": cand.cell_type,
        "significant_transcripts": sig,
        "comparator_transcript": comp,
        "candidate_type": ctype,
        "tier": next((c.tier for c in ctx.m03.classifications
                      if c.tier is not None), None),
        "splice_classification": ctype,
        "docking_readiness": verdict,
        "docking_readiness_reason": reason,
        "uniprot": ctx.m02.uniprot,
        "uniprot_reviewed": ctx.m02.reviewed,
        "pdb_structures": [asdict(p) for p in ctx.m02.pdb],
        "drug_target_status": ctx.m04.drug_target_status,
        "ad_relevance": ctx.m04.ad_relevance,
        "n_drugs": len(ctx.m04.drugs),
        "n_docking_compounds": len(ctx.m04.docking_candidates),
        "opentargets_score": ctx.m05.ad_score,
        "opentargets_label": ctx.m05.ad_label,
        "ad_specificity_ratio": ctx.m05.ad_specificity_ratio,
        "has_genetic_evidence": ctx.m05.has_genetic_evidence,
        "brain_specificity_index": g.bsi if g else None,
        "brain_specificity_label": g.bsi_label if g else "unavailable",
        "expression_pattern": ctx.m07.expression_pattern,
        "comparator_in_pseudobulk": ctx.m07.comparator_in_pseudobulk,
        "robust_to_braak": any(t.robust_to_braak for t in cand.transcripts),
        "braak_condition_correlation": ctx.braak_correlation,
        "pathway_ad_flag": len(ctx.m06.ad_pathways) > 0,
        "interaction_ad_genes": [i.partner for i in ctx.m06.ad_gene_interactions],
        "database_versions": ctx.database_versions,
        "parameters": CONFIG.thresholds,
        "input_file_checksums": ctx.input_checksums,
        "generated": datetime.now(timezone.utc).isoformat(),
    }


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;margin:2em auto;
color:#1a1a1a;line-height:1.5;padding:0 1em}
h1{font-size:1.6em;margin-bottom:0}h2{font-size:1.2em;border-bottom:2px solid #eee;padding-bottom:.2em;margin-top:1.6em}
.sub{color:#666;font-size:.9em}.verdict{display:inline-block;padding:.3em .8em;border-radius:4px;font-weight:bold;color:#fff}
.proceed{background:#2ca02c}.conditional{background:#ff7f0e}.no-go{background:#d62728}
table{border-collapse:collapse;width:100%;font-size:.9em;margin:.5em 0}
th,td{border:1px solid #ddd;padding:.3em .5em;text-align:left}th{background:#f5f5f5}
.chk{font-family:monospace;white-space:pre-wrap;background:#fafafa;border:1px solid #eee;padding:.8em;font-size:.85em}
.tag{display:inline-block;background:#eef;border-radius:3px;padding:.1em .5em;margin:.1em;font-size:.85em}
.muted{color:#999}details{margin:.4em 0}summary{cursor:pointer;font-weight:bold}
"""


def _checklist(ctx: CandidateContext, ctype: str, verdict: str, reason: str) -> str:
    m01, m02, m07 = ctx.m01, ctx.m02, ctx.m07
    s0 = m01.significant[0] if m01.significant else None
    lines = []
    if s0:
        lines.append(f"[ok]  Protein sequence: {s0.transcript.protein_length} aa "
                     f"(query), {s0.comparator.protein_length} aa (comparator)")
        bt = s0.transcript.biotype
        lines.append((f"[ok]  Biotype: protein-coding" if s0.transcript.is_protein_coding
                      else f"[!]   Biotype: {bt} — potential non-coding"))
        if ctype == "N":
            lines.append("[N]   Protein-level change: NONE — protein identical/non-coding")
        else:
            lines.append(f"[ok]  Protein-level change: {s0.protein_length_delta:+d} aa vs comparator")
    cls = next((c for c in ctx.m03.classifications), None)
    if cls:
        lines.append(f"[{ctype}]   Splice classification: Type {ctype} — {cls.narrative}")
    usable_pdb = [p for p in m02.pdb if p.resolution_ok]
    if usable_pdb:
        p = usable_pdb[0]
        lines.append(f"[ok]  Canonical structure: PDB {p.pdb_id} ({p.method}, {p.resolution_a} A)")
    else:
        lines.append("[!]   Canonical structure: none usable — AlphaFold required")
    lines.append(f"[{'ok' if ctx.m04.drug_target_status=='known' else '!'}]   "
                 f"Drugs: {ctx.m04.n_direct} AD-direct / {len(ctx.m04.drugs)} total "
                 f"({ctx.m04.drug_target_status})")
    lines.append(f"[ok]  OpenTargets: {ctx.m05.ad_score} ({ctx.m05.ad_label})")
    if m07.gtex:
        lines.append(f"[ok]  Brain specificity: BSI={m07.gtex.bsi} ({m07.gtex.bsi_label})")
    if any(t.robust_to_braak for t in ctx.candidate.transcripts) is False:
        r = ctx.braak_correlation
        lines.append(f"[!]   robust_to_braak = FALSE — Braak/condition collinearity "
                     f"(r={r:.3f}); see note. Does not affect verdict." if r is not None
                     else "[!]   robust_to_braak = FALSE — see note")
    lines.append(f"\n=> VERDICT: {verdict.upper()}\n   Reason: {reason}")
    return "\n".join(lines)


def render_html(ctx: CandidateContext, dossier: dict, ctype: str,
                verdict: str, reason: str) -> str:
    cand = ctx.candidate
    parts = [f"<!doctype html><html><head><meta charset='utf-8'>"
             f"<title>Surveyor — {cand.gene_id} ({cand.cell_type})</title>"
             f"<style>{_CSS}</style></head><body>"]

    # Header
    parts.append(f"<h1>{cand.gene_id} <span class='sub'>· {cand.cell_type.replace('_',' ')}</span></h1>")
    parts.append(f"<p class='sub'>UniProt {ctx.m02.uniprot or '-'} · "
                 f"Type {ctype} · tier {dossier['tier'] or '-'} · "
                 f"<span class='verdict {verdict}'>{verdict.upper()}</span></p>")

    # Executive summary (auto-generated)
    tx_names = ", ".join(t.transcript_id for t in cand.transcripts)
    parts.append(f"<p><b>Summary.</b> {cand.gene_id} shows differential transcript "
                 f"usage ({tx_names}) in {cand.cell_type.replace('_',' ')}. "
                 f"Classified <b>Type {ctype}</b>: {reason}. OpenTargets AD score "
                 f"{ctx.m05.ad_score} ({ctx.m05.ad_label}); expression pattern "
                 f"{ctx.m07.expression_pattern}.</p>")

    # Checklist
    parts.append("<h2>Docking readiness</h2>")
    parts.append(f"<div class='chk'>{_checklist(ctx, ctype, verdict, reason)}</div>")

    # Section 1 — isoform architecture
    parts.append("<h2>1 · Isoform architecture</h2><table>"
                 "<tr><th>transcript</th><th>role</th><th>ΔPSI</th><th>type</th>"
                 "<th>narrative</th></tr>")
    cls_by = {c.transcript_name: c for c in ctx.m03.classifications}
    for t in cand.transcripts:
        c = cls_by.get(t.transcript_id)
        parts.append(f"<tr><td>{t.transcript_id}</td><td>{t.role}</td>"
                     f"<td>{t.delta_psi:+.3f}</td><td>{c.splice_type if c else '-'}</td>"
                     f"<td>{c.narrative if c else '-'}</td></tr>")
    parts.append("</table>")
    if ctx.m02.pdb:
        parts.append("<details><summary>PDB structures</summary><table>"
                     "<tr><th>PDB</th><th>method</th><th>res (Å)</th><th>chains</th></tr>")
        for p in ctx.m02.pdb[:8]:
            parts.append(f"<tr><td>{p.pdb_id}</td><td>{p.method}</td>"
                         f"<td>{p.resolution_a}</td><td>{p.chains}</td></tr>")
        parts.append("</table></details>")

    # Section 2 — drugs
    parts.append("<h2>2 · Drug target profile</h2>")
    if ctx.m04.drugs:
        parts.append("<table><tr><th>drug</th><th>phase</th><th>action</th>"
                     "<th>AD-direct</th></tr>")
        for d in ctx.m04.drugs[:15]:
            parts.append(f"<tr><td>{d.name}</td><td>{d.max_phase}</td>"
                         f"<td>{d.action_type}</td><td>{'yes' if d.ad_direct else '-'}</td></tr>")
        parts.append("</table>")
    else:
        parts.append(f"<p class='muted'>No ChEMBL drugs — {ctx.m04.drug_target_status} target. "
                     f"{len(ctx.m04.docking_candidates)} sub-µM bioactive compounds.</p>")

    # Section 3 — disease association
    parts.append("<h2>3 · Disease association</h2>")
    parts.append(f"<p>OpenTargets AD: <b>{ctx.m05.ad_score}</b> ({ctx.m05.ad_label}); "
                 f"AD-specificity ratio {ctx.m05.ad_specificity_ratio}; "
                 f"genetic evidence: {'yes' if ctx.m05.has_genetic_evidence else 'no'}.</p>")
    if ctx.m05.ad_datatypes:
        parts.append("<p>" + " ".join(f"<span class='tag'>{k}:{v}</span>"
                                       for k, v in ctx.m05.ad_datatypes.items()) + "</p>")

    # Section 4 — biological context
    parts.append("<h2>4 · Biological context</h2>")
    if ctx.m06.ad_gene_interactions:
        parts.append("<p><b>AD-gene interactions:</b> " +
                     ", ".join(f"{i.partner} ({i.score})" for i in ctx.m06.ad_gene_interactions)
                     + "</p>")
    if ctx.m06.interactions:
        parts.append("<p><b>Top STRING partners:</b> " +
                     ", ".join(i.partner for i in ctx.m06.interactions[:10]) + "</p>")
    if ctx.m06.ad_pathways:
        parts.append("<p><b>AD-relevant pathways:</b> " +
                     "; ".join(p.name for p in ctx.m06.ad_pathways) + "</p>")

    # Section 5 — expression
    parts.append("<h2>5 · Expression profile</h2>")
    parts.append(f"<p>Pattern: <b>{ctx.m07.expression_pattern}</b>. "
                 f"Brain specificity: {dossier['brain_specificity_index']} "
                 f"({dossier['brain_specificity_label']}).</p>")
    parts.append("<table><tr><th>transcript</th><th>role</th>" +
                 "".join(f"<th>{c}</th>" for c in ("Control", "AD", "Active control"))
                 + "</tr>")
    for tp in ctx.m07.transcript_psi:
        cm = {c.condition: c for c in tp.per_condition}
        cells = "".join(
            f"<td>{cm[c].mean:.2f}±{cm[c].sd:.2f} (n={cm[c].n})</td>" if c in cm else "<td>-</td>"
            for c in ("Control", "AD", "Active control"))
        parts.append(f"<tr><td>{tp.transcript_name}</td><td>{tp.role}</td>{cells}</tr>")
    parts.append("</table>")
    for name, svg in ctx.viz_svgs.items():
        parts.append(f"<div>{svg}</div>")

    # Footer
    parts.append(f"<h2 class='sub'>Provenance</h2><p class='sub'>Generated "
                 f"{dossier['generated']} · DB versions: {ctx.database_versions} · "
                 f"input md5: {ctx.input_checksums}</p>")
    parts.append("</body></html>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(ctx: CandidateContext, config=CONFIG) -> dict:
    """Assemble dossier + HTML for one candidate. Writes both files and the
    AlphaFold package for docking-relevant verdicts. Returns the run-summary row."""
    config.ensure_output_dirs()
    ctype, verdict, reason = compute_verdict(ctx)
    dossier = build_dossier(ctx, ctype, verdict, reason)

    stem = f"{ctx.candidate.gene_id}_{ctx.candidate.cell_type}"
    json_path = config.candidates_dir / f"{stem}_dossier.json"
    json_path.write_text(json.dumps(dossier, indent=2, default=str))
    html_path = config.reports_dir / f"{stem}_surveyor_report.html"
    html_path.write_text(render_html(ctx, dossier, ctype, verdict, reason))

    # AlphaFold package for docking-relevant verdicts (Type A/B/C, not D/N)
    if ctype in ("A", "B", "C") and verdict in ("proceed", "conditional"):
        _write_alphafold_package(ctx, config)

    return {"gene": ctx.candidate.gene_id, "cell_type": ctx.candidate.cell_type,
            "type": ctype, "tier": dossier["tier"], "verdict": verdict, "reason": reason}


def _write_alphafold_package(ctx: CandidateContext, config):
    stem = f"{ctx.candidate.gene_id}_{ctx.candidate.cell_type}"
    af_dir = config.alphafold_dir / stem
    af_dir.mkdir(parents=True, exist_ok=True)
    usable_pdb = any(p.resolution_ok for p in ctx.m02.pdb)
    for s in ctx.m01.significant:
        if s.role == "ad_enriched" and s.transcript.protein_seq:
            (af_dir / f"{ctx.candidate.gene_id}_{s.transcript.enst_id}_ad_enriched.fasta").write_text(
                f">{s.transcript.transcript_name}|{s.transcript.enst_id}\n{s.transcript.protein_seq}\n")
            if not usable_pdb and s.comparator.protein_seq:
                (af_dir / f"{ctx.candidate.gene_id}_{s.comparator.enst_id}_comparator.fasta").write_text(
                    f">{s.comparator.transcript_name}|{s.comparator.enst_id}\n{s.comparator.protein_seq}\n")
    params = {"recommended_model": "alphafold2_multimer_v3" if False else "alphafold2_ptm",
              "msa_mode": "mmseqs2_uniref_env",
              "use_templates": usable_pdb,
              "note": "comparator FASTA omitted — usable PDB exists" if usable_pdb else ""}
    (af_dir / "submission_params.json").write_text(json.dumps(params, indent=2))
