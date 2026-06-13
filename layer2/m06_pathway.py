"""
M06 — Pathway and Interaction Analysis (plan 5, M06).

Databases: STRING (protein-protein interactions) and Reactome (pathway
membership).

STRING: high-confidence interaction partners (combined score >= 0.7); the top
partners are cross-referenced against the curated AD-gene list, and any direct
high-confidence interaction with an AD gene is flagged (a candidate that binds
MAPT or APP is mechanistically positioned for AD regardless of its own score).

Reactome: pathway membership for the canonical UniProt accession; pathways whose
names match the AD-relevant pathway list are flagged.

Both databases are non-blocking (plan 7): on failure the section is marked
unavailable and the run continues.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from layer2.config import CONFIG
from layer2.inputs import GeneCandidate
from layer2.utils.api_client import APIClient, APIError
from layer2.utils.audit import AuditLog


@dataclass
class Interaction:
    partner: str
    score: float
    is_ad_gene: bool


@dataclass
class Pathway:
    st_id: str
    name: str
    ad_relevant: bool


@dataclass
class PathwayResult:
    gene_id: str
    interactions: list[Interaction] = field(default_factory=list)
    ad_gene_interactions: list[Interaction] = field(default_factory=list)
    pathways: list[Pathway] = field(default_factory=list)
    ad_pathways: list[Pathway] = field(default_factory=list)
    string_available: bool = True
    reactome_available: bool = True
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# STRING
# ---------------------------------------------------------------------------

def _string_partners(gene: str, client: APIClient, audit: AuditLog,
                     ad_genes: set[str], min_score: float, limit: int = 20
                     ) -> tuple[list[Interaction], bool]:
    url = f"{CONFIG.apis['string'].rstrip('/')}/json/interaction_partners"
    params = {"identifiers": gene, "species": 9606,
              "required_score": int(min_score * 1000), "limit": limit}
    try:
        data = client.get_json("string", url, params=params)
    except APIError as e:
        audit.warning("string", f"{gene}: {e}")
        return [], False
    out = []
    for x in data:
        partner = x.get("preferredName_B", "")
        out.append(Interaction(partner=partner, score=round(float(x.get("score", 0)), 3),
                               is_ad_gene=partner in ad_genes))
    out.sort(key=lambda i: i.score, reverse=True)
    return out, True


# ---------------------------------------------------------------------------
# Reactome
# ---------------------------------------------------------------------------

def _ad_relevant(name: str, ad_pathways: list[str]) -> bool:
    low = name.lower()
    for phrase in ad_pathways:
        # match if all significant words of the phrase appear in the pathway name
        words = [w for w in phrase.split() if len(w) > 3]
        if words and all(w in low for w in words):
            return True
    # common single-token hits
    return any(tok in low for tok in
               ("amyloid", "tau ", "mitophag", "autophag", "axonal transport",
                "long-term potentiation", "synaptic"))


def _reactome_pathways(uniprot: str | None, client: APIClient, audit: AuditLog,
                       ad_pathways: list[str]) -> tuple[list[Pathway], bool]:
    if not uniprot:
        return [], True
    url = f"{CONFIG.apis['reactome'].rstrip('/')}/data/mapping/UniProt/{uniprot}/pathways"
    try:
        data = client.get_json("reactome", url, params={"species": 9606})
    except APIError as e:
        audit.warning("reactome", f"{uniprot}: {e}")
        return [], False
    if not isinstance(data, list):
        return [], True
    out = []
    for p in data:
        name = p.get("displayName", "")
        out.append(Pathway(st_id=p.get("stId", ""), name=name,
                           ad_relevant=_ad_relevant(name, ad_pathways)))
    return out, True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(gene: GeneCandidate, uniprot: str | None, client: APIClient,
        audit: AuditLog, config=CONFIG) -> PathwayResult:
    ad_genes = set(config.ad_genes)
    min_score = config.thresholds["string_confidence_high"]

    interactions, string_ok = _string_partners(
        gene.gene_id, client, audit, ad_genes, min_score)
    pathways, reactome_ok = _reactome_pathways(
        uniprot, client, audit, config.ad_pathways)

    res = PathwayResult(
        gene_id=gene.gene_id,
        interactions=interactions,
        ad_gene_interactions=[i for i in interactions if i.is_ad_gene],
        pathways=pathways,
        ad_pathways=[p for p in pathways if p.ad_relevant],
        string_available=string_ok, reactome_available=reactome_ok,
    )
    audit.info(f"M06 {gene.gene_id}: {len(interactions)} STRING partners "
               f"({len(res.ad_gene_interactions)} AD-genes), "
               f"{len(pathways)} Reactome pathways ({len(res.ad_pathways)} AD-relevant)")
    return res


if __name__ == "__main__":
    from layer2.inputs import load_inputs
    from layer2 import m01_sequence as m01_mod
    audit = AuditLog(CONFIG.audit_log_path, verbose=False)
    client = APIClient(audit)
    inp = load_inputs()
    print(f"{'gene':9} {'partners':>8} {'ad_int':>6} {'paths':>6} {'ad_path':>7}  "
          f"AD-gene interactions / AD pathways")
    print("-" * 110)
    for gene in inp.gene_work_list:
        m01r = m01_mod.run(gene, client, audit, inp.enst_map)
        r = run(gene, m01r.canonical_uniprot, client, audit)
        adi = ", ".join(f"{i.partner}({i.score})" for i in r.ad_gene_interactions) or "-"
        adp = "; ".join(p.name for p in r.ad_pathways[:2]) or "-"
        print(f"{gene.gene_id:9} {len(r.interactions):8} {len(r.ad_gene_interactions):6} "
              f"{len(r.pathways):6} {len(r.ad_pathways):7}  {adi[:40]} | {adp[:40]}")
