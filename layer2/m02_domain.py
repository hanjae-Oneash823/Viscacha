"""
M02 — Domain Annotation, PDB Lookup, and Residue Remapping (plan 5, M02).

Database: UniProt REST (https://rest.uniprot.org/uniprotkb/{acc}.json).

For each gene's canonical UniProt accession (resolved by M01) this module:
  1. Retrieves all sequence features (Domain, Binding site, Active site, Repeat,
     Zinc finger, Region, Motif, Modified residue, ...) with evidence type, plus
     the reviewed/unreviewed status and the canonical sequence.
  2. Retrieves PDB cross-references (id, method, resolution, chain ranges) and
     flags those usable as docking templates (resolution <= threshold).
  3. For each significant transcript that has a real protein, aligns the
     canonical UniProt sequence to the isoform (BLOSUM62 global) and:
       - remaps every feature into isoform coordinates (present/truncated/absent)
       - identifies the splice-affected canonical residues (the alignment
         divergence), used by M03 for classification and for PDB coverage
       - flags PTM sites within 5 residues of the affected region

The affected-residue set and the feature table are the inputs M03 consumes.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from Bio.Align import PairwiseAligner, substitution_matrices

from layer2.config import CONFIG
from layer2.m01_sequence import M01Result, SignificantTranscriptM01
from layer2.utils.api_client import APIClient, APIError
from layer2.utils.audit import AuditLog

# Feature vocabulary (broadened in plan v2.2 after empirical validation).
POCKET_TYPES = {"Binding site", "Active site", "Site"}
STRUCT_TYPES = {"Domain", "Repeat", "Zinc finger", "DNA binding",
                "Coiled coil", "Motif", "Region"}
PTM_TYPES = {"Modified residue", "Lipidation", "Glycosylation", "Cross-link"}

# ECO evidence code -> coarse confidence tier
_EXPERIMENTAL_ECO = {"ECO:0000269", "ECO:0000314", "ECO:0000353", "ECO:0000305"}
_BYSIM_ECO = {"ECO:0000250", "ECO:0000266"}
_PREDICTED_ECO = {"ECO:0000255", "ECO:0000256", "ECO:0000259", "ECO:0000213"}

PTM_PROXIMITY_RESIDUES = 5


class M02Error(Exception):
    """Blocking M02 failure (UniProt unavailable) — plan 7."""


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Feature:
    type: str
    start: int                 # canonical (UniProt) coordinates, 1-based
    end: int
    description: str
    evidence: str              # experimental | by_similarity | predicted | unknown

    @property
    def is_pocket(self) -> bool:
        return self.type in POCKET_TYPES

    @property
    def is_struct(self) -> bool:
        return self.type in STRUCT_TYPES


@dataclass
class RemappedFeature:
    feature: Feature
    isoform_start: int | None
    isoform_end: int | None
    status: str                # present | truncated | absent


@dataclass
class PDBStructure:
    pdb_id: str
    method: str
    resolution_a: float | None
    chains: str
    resolution_ok: bool


@dataclass
class TranscriptDomainMap:
    transcript_name: str
    remapped: list[RemappedFeature]
    affected_canonical: list[int]            # sorted residues that diverge
    affected_ranges: list[tuple[int, int]]   # contiguous ranges of the above
    deleted_canonical: set[int]              # subset: deleted in isoform
    ptm_proximity: list[Feature]             # PTM within 5 res of affected
    pdb_covers_splice: dict[str, bool]       # pdb_id -> covers affected region


@dataclass
class M02Result:
    gene_id: str
    uniprot: str | None
    reviewed: bool
    canonical_length: int
    features: list[Feature] = field(default_factory=list)
    pdb: list[PDBStructure] = field(default_factory=list)
    transcript_maps: dict[str, TranscriptDomainMap] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    available: bool = True                    # False if UniProt had no entry


# ---------------------------------------------------------------------------
# UniProt parsing
# ---------------------------------------------------------------------------

def _evidence_tier(evidences: list[dict]) -> str:
    codes = {e.get("evidenceCode") for e in (evidences or [])}
    if codes & _EXPERIMENTAL_ECO:
        return "experimental"
    if codes & _BYSIM_ECO:
        return "by_similarity"
    if codes & _PREDICTED_ECO:
        return "predicted"
    return "unknown"


def _parse_features(entry: dict) -> list[Feature]:
    out = []
    for f in entry.get("features", []):
        loc = f.get("location", {})
        s = (loc.get("start") or {}).get("value")
        e = (loc.get("end") or {}).get("value")
        if s is None or e is None:
            continue
        out.append(Feature(
            type=f.get("type", ""),
            start=int(s), end=int(e),
            description=f.get("description", "") or "",
            evidence=_evidence_tier(f.get("evidences")),
        ))
    return out


def _parse_resolution(value: str | None) -> float | None:
    if not value:
        return None
    m = re.search(r"([\d.]+)", value)
    return float(m.group(1)) if m else None


def _parse_pdb(entry: dict, max_res: float) -> list[PDBStructure]:
    out = []
    for x in entry.get("uniProtKBCrossReferences", []):
        if x.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in x.get("properties", [])}
        res = _parse_resolution(props.get("Resolution"))
        out.append(PDBStructure(
            pdb_id=x.get("id", ""),
            method=props.get("Method", ""),
            resolution_a=res,
            chains=props.get("Chains", ""),
            resolution_ok=(res is not None and res <= max_res),
        ))
    return out


def _pdb_chain_ranges(chains: str) -> list[tuple[int, int]]:
    """Parse a PDB 'Chains' string like 'A/B=411-592, A=177-592' into residue
    ranges (UniProt canonical coordinates)."""
    ranges = []
    for part in chains.split(","):
        m = re.search(r"=(\d+)-(\d+)", part)
        if m:
            ranges.append((int(m.group(1)), int(m.group(2))))
    return ranges


# ---------------------------------------------------------------------------
# Alignment & remapping
# ---------------------------------------------------------------------------

_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALIGNER.mode = "global"
_ALIGNER.open_gap_score = -10
_ALIGNER.extend_gap_score = -0.5


def _align_and_map(canonical: str, isoform: str
                   ) -> tuple[dict[int, int], set[int], set[int]]:
    """Global-align canonical->isoform. Returns:
      pos_map: canonical residue (1-based) -> isoform residue (1-based), for
               positions present in both
      affected: canonical residues that diverge (deleted or substituted, or
                adjacent to an isoform insertion)
      deleted:  subset of affected that are absent from the isoform
    """
    aln = _ALIGNER.align(canonical, isoform)[0]
    c_blocks, i_blocks = aln.aligned    # arrays of [start,end) 0-based

    pos_map: dict[int, int] = {}
    covered_c: set[int] = set()
    affected: set[int] = set()

    for (c0, c1), (i0, i1) in zip(c_blocks, i_blocks):
        for off in range(c1 - c0):
            c_res = c0 + off + 1            # 1-based canonical
            i_res = i0 + off + 1            # 1-based isoform
            pos_map[c_res] = i_res
            covered_c.add(c_res)
            if canonical[c0 + off] != isoform[i0 + off]:
                affected.add(c_res)         # substitution

    # Deletions: canonical residues not covered by any aligned block
    deleted = {c + 1 for c in range(len(canonical)) if (c + 1) not in covered_c}
    affected |= deleted

    # Insertions in the isoform: mark the flanking canonical residue as affected
    prev_c1 = 0
    for (c0, c1), (i0, i1) in zip(c_blocks, i_blocks):
        if i0 > 0 and c0 == prev_c1 and c0 >= 1:
            affected.add(c0)               # canonical residue just before insert
        prev_c1 = c1

    return pos_map, affected, deleted


def _contiguous_ranges(positions: set[int]) -> list[tuple[int, int]]:
    if not positions:
        return []
    s = sorted(positions)
    ranges, start, prev = [], s[0], s[0]
    for p in s[1:]:
        if p == prev + 1:
            prev = p
        else:
            ranges.append((start, prev))
            start = prev = p
    ranges.append((start, prev))
    return ranges


def _remap_feature(feat: Feature, pos_map: dict[int, int],
                   deleted: set[int]) -> RemappedFeature:
    span = range(feat.start, feat.end + 1)
    mapped = [pos_map[c] for c in span if c in pos_map]
    n_deleted = sum(1 for c in span if c in deleted)
    if not mapped:
        status = "absent"
        return RemappedFeature(feat, None, None, status)
    status = "truncated" if n_deleted > 0 else "present"
    return RemappedFeature(feat, min(mapped), max(mapped), status)


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def _fetch_uniprot(client: APIClient, accession: str) -> dict:
    url = f"{CONFIG.apis['uniprot'].rstrip('/')}/{accession}.json"
    return client.get_json("uniprot", url)


def run(m01: M01Result, client: APIClient, audit: AuditLog,
        config=CONFIG) -> M02Result:
    """Run M02 for one gene. Raises M02Error on a blocking failure (plan 7:
    UniProt unavailable -> M03 cannot classify; caller forces CONDITIONAL)."""
    acc = m01.canonical_uniprot
    max_res = config.thresholds["pdb_resolution_max_angstrom"]

    if not acc:
        return M02Result(gene_id=m01.gene_id, uniprot=None, reviewed=False,
                         canonical_length=0, available=False,
                         warnings=["no canonical UniProt accession from M01"])

    try:
        entry = _fetch_uniprot(client, acc)
    except APIError as e:
        raise M02Error(f"{m01.gene_id}: UniProt fetch failed for {acc}: {e}")

    canonical_seq = (entry.get("sequence") or {}).get("value", "")
    reviewed = "reviewed" in (entry.get("entryType", "").lower())
    features = _parse_features(entry)
    pdb = _parse_pdb(entry, max_res)
    audit.db_version("uniprot", entry.get("entryAudit", {}).get("lastSequenceUpdateDate", "?"))

    result = M02Result(
        gene_id=m01.gene_id, uniprot=acc, reviewed=reviewed,
        canonical_length=len(canonical_seq), features=features, pdb=pdb,
    )

    pdb_ranges = {p.pdb_id: _pdb_chain_ranges(p.chains) for p in pdb}

    # Cache alignments canonical->seq by sequence (comparators repeat across pairs)
    align_cache: dict[str, tuple[dict[int, int], set[int], set[int]]] = {}

    def aligned(seq: str):
        if seq not in align_cache:
            align_cache[seq] = _align_and_map(canonical_seq, seq)
        return align_cache[seq]

    # Per significant transcript: the splice-affected region is what differs
    # between the AD-enriched isoform and ITS COMPARATOR, expressed in canonical
    # coordinates (so it lines up with the features). Computed as the symmetric
    # difference of each one's divergence from the shared canonical frame; when
    # the comparator is the canonical/MANE sequence this reduces to the isoform's
    # own divergence.
    for s in m01.significant:
        iso_seq = s.transcript.protein_seq
        comp_seq = s.comparator.protein_seq
        tname = s.transcript.transcript_name
        if not iso_seq or not canonical_seq:
            continue  # non-coding (Type N) — nothing to remap

        pos_map_q, affected_q, deleted_q = aligned(iso_seq)
        if comp_seq:
            _, affected_c, deleted_c = aligned(comp_seq)
        else:
            affected_c, deleted_c = set(), set()

        splice_affected = affected_q ^ affected_c      # differs between the pair
        splice_deleted = deleted_q - deleted_c         # lost in isoform vs comparator
        affected_ranges = _contiguous_ranges(splice_affected)
        remapped = [_remap_feature(f, pos_map_q, deleted_q) for f in features]

        # PTM proximity: PTM feature within N residues of the splice region
        ptm_prox = []
        for f in features:
            if f.type in PTM_TYPES and splice_affected:
                if min(abs((f.start + f.end) // 2 - a)
                       for a in splice_affected) <= PTM_PROXIMITY_RESIDUES:
                    ptm_prox.append(f)

        # PDB coverage of the splice region, in canonical coords
        covers = {}
        for pid, ranges in pdb_ranges.items():
            covers[pid] = any(
                not (ar[1] < pr[0] or ar[0] > pr[1])
                for ar in affected_ranges for pr in ranges
            )

        result.transcript_maps[tname] = TranscriptDomainMap(
            transcript_name=tname, remapped=remapped,
            affected_canonical=sorted(splice_affected), affected_ranges=affected_ranges,
            deleted_canonical=splice_deleted, ptm_proximity=ptm_prox,
            pdb_covers_splice=covers,
        )

    n_pocket = sum(1 for f in features if f.is_pocket)
    n_struct = sum(1 for f in features if f.is_struct)
    audit.info(f"M02 {m01.gene_id} ({acc}, {'reviewed' if reviewed else 'unreviewed'}): "
               f"{len(features)} features ({n_pocket} pocket, {n_struct} struct), "
               f"{len(pdb)} PDB ({sum(p.resolution_ok for p in pdb)} usable)")
    return result


# ---------------------------------------------------------------------------
# Standalone verification
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from layer2.inputs import load_inputs
    from layer2 import m01_sequence as m01_mod

    audit = AuditLog(CONFIG.audit_log_path, verbose=False)
    client = APIClient(audit)
    inp = load_inputs()

    print(f"{'gene':9} {'uniprot':9} {'rev':4} {'feat':>4} {'pkt':>3} {'str':>3} "
          f"{'pdb':>3} {'usable':>6}  affected(canonical) per sig tx")
    print("-" * 100)
    for gene in inp.gene_work_list:
        m01r = m01_mod.run(gene, client, audit, inp.enst_map)
        try:
            m02r = run(m01r, client, audit)
        except M02Error as e:
            print(f"{gene.gene_id:9} M02Error: {e}")
            continue
        npk = sum(f.is_pocket for f in m02r.features)
        nst = sum(f.is_struct for f in m02r.features)
        usable = sum(p.resolution_ok for p in m02r.pdb)
        tmap_str = "; ".join(
            f"{t}:{len(tm.affected_canonical)}res/{len(tm.affected_ranges)}rng"
            for t, tm in m02r.transcript_maps.items()) or "(none - Type N)"
        print(f"{m02r.gene_id:9} {str(m02r.uniprot):9} "
              f"{'Y' if m02r.reviewed else 'n':4} {len(m02r.features):4} "
              f"{npk:3} {nst:3} {len(m02r.pdb):3} {usable:6}  {tmap_str}")
