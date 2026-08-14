"""Canonical<->alt residue position correspondence, from the SAME edlib
alignment convention junior_surveyor/j2_protein_diff.py already uses to
compute hits_deep.csv's changed_aa_start/changed_aa_end (1-indexed,
CANONICAL-sequence coordinates, query=canonical in edlib's NW alignment).

m2b_structure_qc.py needs more than j2_protein_diff.py already exposes: the
ALT-local position for the affected span (structures are numbered 1..len(alt
sequence), not canonical coordinates) and a full column-by-column
correspondence table to use as rigid-body superposition anchors. Rather than
importing junior_surveyor's private `_align_proteins` (a different package,
and it also carries the change-type classifier we don't need here) or
duplicating its full logic, this module re-derives just the position
bookkeeping with the identical edlib call so the two stay consistent with
each other by construction, not by convention alone.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import edlib

_CIGAR_RE = re.compile(r"(\d+)([=XID])")


@dataclass
class AlignedColumn:
    canon_pos: int | None   # 1-indexed canonical position, None if alt-only insertion
    alt_pos: int | None     # 1-indexed alt position, None if canonical-only deletion
    op: str                 # "=" match, "X" mismatch, "I" canonical-only, "D" alt-only


def residue_correspondence(canonical: str, alt: str) -> list[AlignedColumn]:
    """Full per-op correspondence table, same CIGAR semantics as
    j2_protein_diff.py's _align_proteins (query=canonical, mode="NW"):
        = : match            -- both advance
        X : mismatch         -- both advance
        I : canonical has it, alt doesn't -- canonical advances only
        D : alt has it, canonical doesn't -- alt advances only
    """
    result = edlib.align(canonical, alt, mode="NW", task="path")
    ops = [(int(n), op) for n, op in _CIGAR_RE.findall(result["cigar"])]

    canon_pos, alt_pos = 0, 0
    columns: list[AlignedColumn] = []
    for count, op in ops:
        for _ in range(count):
            if op == "=" or op == "X":
                canon_pos += 1
                alt_pos += 1
                columns.append(AlignedColumn(canon_pos, alt_pos, op))
            elif op == "I":
                canon_pos += 1
                columns.append(AlignedColumn(canon_pos, None, op))
            elif op == "D":
                alt_pos += 1
                columns.append(AlignedColumn(None, alt_pos, op))
    return columns


def canonical_span_to_alt(columns: list[AlignedColumn], canon_start: int, canon_end: int) -> tuple[int, int] | None:
    """Map a [canon_start, canon_end] (1-indexed, inclusive) span to the
    corresponding alt-local span, for pulling pLDDT/PAE out of the ALT
    structure (numbered 1..len(alt_seq), not canonical coordinates). Returns
    None if no aligned column falls in the requested canonical range (can
    happen for a pure canonical-only deletion span with no alt anchor at
    all, though in practice the D-op's canonical anchor -- see
    j2_protein_diff.py's "clamped >= 1" comment -- means this is rare).
    """
    alt_positions = [
        c.alt_pos for c in columns
        if c.canon_pos is not None and canon_start <= c.canon_pos <= canon_end and c.alt_pos is not None
    ]
    if not alt_positions:
        return None
    return min(alt_positions), max(alt_positions)


def match_anchors(columns: list[AlignedColumn]) -> list[tuple[int, int]]:
    """(canon_pos, alt_pos) pairs for every exact-match ("=") column -- the
    only positions guaranteed structurally equivalent by sequence identity,
    used as rigid-body superposition anchors in m2b_structure_qc.py rather
    than trusting a raw 1:1 index correspondence (which breaks the moment
    an indel shifts the two sequences out of register).
    """
    return [(c.canon_pos, c.alt_pos) for c in columns if c.op == "="]
