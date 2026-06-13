"""
M05 — Disease Association (plan 5, M05).

Database: OpenTargets Platform GraphQL.

Per gene, against Alzheimer's disease (MONDO_0004975) plus three other
neurodegenerative diseases (PD, FTD, ALS) for the AD-specificity ratio. Uses the
disease-side query `disease(efoId){ associatedTargets(Bs:[ensg]) }`, which
returns the overall association score and the datatype breakdown filtered to one
target. All four diseases are requested in a single aliased query.

Scores are reported, never thresholded as pass/fail (plan): understudied genes
score low not because evidence contradicts AD relevance but because they are
under-investigated — exactly the candidates a long-read isoform study surfaces.

Not yet implemented (noted, not silently skipped): GWAS lead-SNP extraction and
splice-site proximity. The genetic-evidence datatype score is captured as a
presence signal; per-variant extraction needs the heavier evidences() query and
is deferred to a later pass.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from layer2.config import CONFIG
from layer2.inputs import GeneCandidate
from layer2.utils.api_client import APIClient, APIError
from layer2.utils.audit import AuditLog


@dataclass
class DiseaseAssociation:
    gene_id: str
    ensg_id: str
    ad_score: float
    ad_label: str                      # supported | emerging | novel
    ad_datatypes: dict[str, float]
    has_genetic_evidence: bool
    scores_other: dict[str, float]     # parkinson / frontotemporal / als
    ad_specificity_ratio: float | None
    available: bool = True
    warnings: list[str] = field(default_factory=list)


def _label(score: float, labels: dict) -> str:
    if score >= labels["supported"]:
        return "supported"
    if score >= labels["emerging"]:
        return "emerging"
    return "novel"


def _build_query(ensg: str, efo: dict) -> str:
    # one aliased query: AD returns datatype breakdown; others only the score
    def block(alias: str, efo_id: str, with_dt: bool) -> str:
        dt = " datatypeScores{ id score }" if with_dt else ""
        return (f'{alias}: disease(efoId: "{efo_id}") {{ '
                f'associatedTargets(Bs: ["{ensg}"]) {{ rows {{ score{dt} }} }} }}')
    return "{ " + " ".join([
        block("ad", efo["alzheimer"], True),
        block("pd", efo["parkinson"], False),
        block("ftd", efo["frontotemporal"], False),
        block("als", efo["als"], False),
    ]) + " }"


def _score_of(node: dict) -> tuple[float, list[dict]]:
    rows = ((node or {}).get("associatedTargets") or {}).get("rows") or []
    if not rows:
        return 0.0, []
    return float(rows[0].get("score") or 0.0), rows[0].get("datatypeScores") or []


def run(gene: GeneCandidate, client: APIClient, audit: AuditLog,
        config=CONFIG) -> DiseaseAssociation:
    efo = config.opentargets_efo
    labels = config.opentargets_labels
    url = config.apis["opentargets"]

    if not gene.ensg_id:
        return DiseaseAssociation(gene.gene_id, "", 0.0, "novel", {}, False, {},
                                  None, available=False,
                                  warnings=["no Ensembl gene ID"])

    try:
        data = client.post_json("opentargets", url,
                                {"query": _build_query(gene.ensg_id, efo)})
    except APIError as e:
        audit.error("opentargets", f"{gene.gene_id}: {e}")
        return DiseaseAssociation(gene.gene_id, gene.ensg_id, 0.0, "novel", {},
                                  False, {}, None, available=False,
                                  warnings=[f"OpenTargets unavailable: {e}"])

    d = data.get("data", {}) or {}
    ad_score, ad_dt_rows = _score_of(d.get("ad"))
    ad_datatypes = {r["id"]: round(float(r["score"]), 4) for r in ad_dt_rows}
    pd_score, _ = _score_of(d.get("pd"))
    ftd_score, _ = _score_of(d.get("ftd"))
    als_score, _ = _score_of(d.get("als"))

    others = {"parkinson": round(pd_score, 4), "frontotemporal": round(ftd_score, 4),
              "als": round(als_score, 4)}
    mean_other = (pd_score + ftd_score + als_score) / 3.0
    if mean_other > 0:
        ratio = round(ad_score / mean_other, 3)
    else:
        ratio = None  # AD-only or no ND association

    has_genetic = ad_datatypes.get("genetic_association", 0.0) > 0
    label = _label(ad_score, labels)

    audit.info(f"M05 {gene.gene_id}: AD={ad_score:.3f} ({label}), "
               f"specificity={ratio}, genetic={'Y' if has_genetic else 'n'}")
    return DiseaseAssociation(
        gene_id=gene.gene_id, ensg_id=gene.ensg_id, ad_score=round(ad_score, 4),
        ad_label=label, ad_datatypes=ad_datatypes, has_genetic_evidence=has_genetic,
        scores_other=others, ad_specificity_ratio=ratio,
    )


if __name__ == "__main__":
    from layer2.inputs import load_inputs
    audit = AuditLog(CONFIG.audit_log_path, verbose=False)
    client = APIClient(audit)
    inp = load_inputs()
    print(f"{'gene':9} {'AD':>6} {'label':10} {'spec':>6} {'genetic':8} other(PD/FTD/ALS)")
    print("-" * 80)
    for gene in inp.gene_work_list:
        r = run(gene, client, audit)
        o = r.scores_other
        print(f"{r.gene_id:9} {r.ad_score:6.3f} {r.ad_label:10} "
              f"{str(r.ad_specificity_ratio or '-'):>6} "
              f"{'yes' if r.has_genetic_evidence else 'no':8} "
              f"{o.get('parkinson',0):.2f}/{o.get('frontotemporal',0):.2f}/{o.get('als',0):.2f}")
