"""
M03 — Splice Event Classification and Tier Assignment (plan 5, M03; v2.2).

Internal computation only — consumes M01 (protein-identity gate, comparator)
and M02 (features, splice-affected residues in canonical coordinates, PDB).

Order of decisions per significant transcript:
  1. Type N gate (v2.2): if the isoform encodes a protein identical to its
     comparator, or is non-coding, there is no docking hypothesis. Classify as
     Type N, no tier, route to the regulatory path. Checked first.
  2. Otherwise classify A / B / C / D from where the splice-affected residues
     fall relative to the functional features:
       A — overlaps a pocket feature, or removes a structural-domain entirely
       B — within 15 residues of a pocket feature, or inside a domain that
           contains a pocket
       C — overlaps a structural-domain feature, distal from any pocket
       D — no functional feature in or near the affected region
  3. Tier: A/B -> 1, C -> 2, D/N -> none.

Feature vocabulary is the broadened v2.2 set (pocket vs structural-domain) from
m02_domain, so Repeat / Zinc finger / Region-class cores (KLC1, MGRN1) are not
mis-called Type D.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from layer2.config import CONFIG
from layer2.m01_sequence import M01Result
from layer2.m02_domain import Feature, M02Result, TranscriptDomainMap


@dataclass
class SpliceClassification:
    transcript_name: str
    role: str                       # ad_enriched | control_enriched
    splice_type: str                # A | B | C | D | N
    tier: int | None
    distance_to_pocket: float | None   # residues; None if protein has no pocket
    affected_features: list[Feature]   # features the splice event touches
    confidence: str                    # experimental | by_similarity | predicted | unknown
    narrative: str
    docking_path: bool                 # True for A/B/C/D; False for N

    @property
    def is_type_n(self) -> bool:
        return self.splice_type == "N"


@dataclass
class M03Result:
    gene_id: str
    classifications: list[SpliceClassification] = field(default_factory=list)
    m02_available: bool = True
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Geometry helpers (canonical residue coordinates)
# ---------------------------------------------------------------------------

def _range_distance(affected: set[int], start: int, end: int) -> float:
    """Min residue distance from the affected set to the feature [start,end];
    0 if any affected residue falls inside the feature."""
    if not affected:
        return math.inf
    best = math.inf
    for p in affected:
        if start <= p <= end:
            return 0.0
        best = min(best, abs(p - start), abs(p - end))
    return best


def _fully_within(feat: Feature, residues: set[int]) -> bool:
    return all(c in residues for c in range(feat.start, feat.end + 1))


def _contains(outer: Feature, inner: Feature) -> bool:
    return outer.start <= inner.start and inner.end <= outer.end


def _is_structural(f: Feature) -> bool:
    """A structural-domain feature for classification. 'Region' is heterogeneous:
    a 'Disordered'/'Compositional' region is the Type D case (plan), while a
    named functional region (e.g. 'Transactivation domain') is structural."""
    if not f.is_struct:
        return False
    desc = f.description.lower()
    if f.type == "Region" and ("disorder" in desc or "compositional" in desc):
        return False
    return True


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def _classify_one(tm: TranscriptDomainMap, features: list[Feature],
                  proximity: int) -> tuple[str, float | None, list[Feature], str, str]:
    """Return (type, distance_to_pocket, affected_features, confidence, narrative)
    for a structurally-different transcript."""
    affected = set(tm.affected_canonical)
    deleted = set(tm.deleted_canonical)
    pockets = [f for f in features if f.is_pocket]
    structs = [f for f in features if _is_structural(f)]

    # Distances to pocket features
    if pockets:
        dist_pocket = min(_range_distance(affected, f.start, f.end) for f in pockets)
    else:
        dist_pocket = None

    overlaps_pocket = dist_pocket == 0
    affected_structs = [f for f in structs
                        if _range_distance(affected, f.start, f.end) == 0]
    removed_structs = [f for f in structs if _fully_within(f, deleted)]
    struct_contains_pocket = any(
        any(_contains(sf, pf) for pf in pockets) for sf in affected_structs)

    # --- decide type ---
    if overlaps_pocket or removed_structs:
        stype = "A"
        gov = (next((f for f in pockets
                     if _range_distance(affected, f.start, f.end) == 0), None)
               or (removed_structs[0] if removed_structs else None))
        if removed_structs and not overlaps_pocket:
            narrative = (f"Splice event removes the {gov.type} "
                         f"'{gov.description or gov.type}' ({gov.start}-{gov.end}) "
                         f"entirely — binding architecture is altered.")
        else:
            narrative = (f"Splice-affected residues directly overlap a "
                         f"{gov.type} ({gov.start}-{gov.end}) — the drug binding "
                         f"pocket differs between isoforms.")
        affected_feats = ([f for f in pockets
                           if _range_distance(affected, f.start, f.end) == 0]
                          + removed_structs)
    elif (dist_pocket is not None and dist_pocket <= proximity) or struct_contains_pocket:
        stype = "B"
        gov = min(pockets, key=lambda f: _range_distance(affected, f.start, f.end))
        narrative = (f"Splice event falls {int(dist_pocket)} residue(s) from a "
                     f"{gov.type} ({gov.start}-{gov.end}); pocket geometry may be "
                     f"indirectly affected.") if dist_pocket is not None else (
                     f"Splice event lies within a domain containing the binding pocket.")
        affected_feats = affected_structs + [gov]
    elif affected_structs:
        stype = "C"
        gov = affected_structs[0]
        d = f"{int(dist_pocket)} residues" if dist_pocket is not None else "n/a (no annotated pocket)"
        narrative = (f"Splice event alters the {gov.type} "
                     f"'{gov.description or gov.type}' ({gov.start}-{gov.end}); "
                     f"nearest pocket {d} away — pocket residues conserved.")
        affected_feats = affected_structs
    else:
        stype = "D"
        narrative = ("Splice event falls in an unannotated / disordered region — "
                     "no functional feature in or near the affected residues.")
        affected_feats = []

    confidence = (affected_feats[0].evidence if affected_feats else "unknown")
    return stype, dist_pocket, affected_feats, confidence, narrative


def run(m01: M01Result, m02: M02Result, config=CONFIG) -> M03Result:
    proximity = config.thresholds["splice_domain_proximity_residues"]
    tier1 = set(config.tier_assignment["tier1"])
    tier2 = set(config.tier_assignment["tier2"])

    res = M03Result(gene_id=m01.gene_id, m02_available=m02.available)

    for s in m01.significant:
        tname = s.transcript.transcript_name

        # 1. Type N gate (checked first)
        if not s.has_structural_difference:
            reason = ("non-coding transcript — no protein product"
                      if not s.transcript.is_protein_coding
                      else "protein identical to comparator")
            res.classifications.append(SpliceClassification(
                transcript_name=tname, role=s.role, splice_type="N", tier=None,
                distance_to_pocket=None, affected_features=[], confidence="n/a",
                narrative=f"No protein-level change vs comparator "
                          f"({reason}) — regulatory candidate, no docking hypothesis.",
                docking_path=False))
            continue

        # M02 unavailable -> cannot classify (caller forces CONDITIONAL)
        if not m02.available or tname not in m02.transcript_maps:
            res.classifications.append(SpliceClassification(
                transcript_name=tname, role=s.role, splice_type="D", tier=None,
                distance_to_pocket=None, affected_features=[], confidence="unknown",
                narrative="Domain annotation unavailable — manual splice "
                          "classification required.",
                docking_path=True))
            res.warnings.append(f"{tname}: classified provisionally (no M02 map)")
            continue

        tm = m02.transcript_maps[tname]
        stype, dist, feats, conf, narrative = _classify_one(tm, m02.features, proximity)
        tier = 1 if stype in tier1 else (2 if stype in tier2 else None)
        res.classifications.append(SpliceClassification(
            transcript_name=tname, role=s.role, splice_type=stype, tier=tier,
            distance_to_pocket=dist, affected_features=feats, confidence=conf,
            narrative=narrative, docking_path=True))

    return res


# ---------------------------------------------------------------------------
# Standalone verification
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from layer2.inputs import load_inputs
    from layer2.utils.api_client import APIClient
    from layer2.utils.audit import AuditLog
    from layer2 import m01_sequence as m01_mod
    from layer2 import m02_domain as m02_mod

    audit = AuditLog(CONFIG.audit_log_path, verbose=False)
    client = APIClient(audit)
    inp = load_inputs()

    print(f"{'gene':9} {'tx':14} {'role':16} {'type':5} {'tier':>4} "
          f"{'dist':>5} {'conf':13} narrative")
    print("-" * 130)
    for gene in inp.gene_work_list:
        m01r = m01_mod.run(gene, client, audit, inp.enst_map)
        m02r = m02_mod.run(m01r, client, audit)
        m03r = run(m01r, m02r)
        for c in m03r.classifications:
            d = f"{c.distance_to_pocket:.0f}" if c.distance_to_pocket is not None else "-"
            print(f"{gene.gene_id:9} {c.transcript_name:14} {c.role:16} "
                  f"{c.splice_type:5} {str(c.tier or '-'):>4} {d:>5} {c.confidence:13} "
                  f"{c.narrative[:60]}")
