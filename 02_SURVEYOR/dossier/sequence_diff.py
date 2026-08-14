"""DOSSIER — canonical-vs-alt protein alignment for the sequence comparison bar.

Reuses junior_surveyor's J2 classifier (_align_proteins) for the
protein_change_type / summary stats so the dossier never disagrees with
hits_deep.csv. align_segments() is a SEPARATE alignment, deliberately not
edlib: edlib's unit-cost edit distance has no preference for consolidating
gaps, so a single 675-aa deletion can come back as ~150 scattered 1-2 residue
fragments (any distribution of the same total gap length costs the same).
That's fine for J2's summary stats (position of first/last change, counts)
but renders as a nonsense fan of ribbons. Biopython's PairwiseAligner with
an affine gap penalty (steep open, shallow extend) prefers few large gaps
over many small ones -- a biologically sane shape for the ribbon diagram.
"""

from __future__ import annotations

from Bio.Align import PairwiseAligner

from junior_surveyor.j2_protein_diff import _align_proteins

summarize = _align_proteins  # public alias

_ALIGNER = PairwiseAligner()
_ALIGNER.mode = "global"
_ALIGNER.match_score = 2
_ALIGNER.mismatch_score = -1
_ALIGNER.open_gap_score = -10
_ALIGNER.extend_gap_score = -0.5


def align_segments(canonical: str, alt: str) -> list[dict]:
    """One entry per contiguous same-op run, in alignment order -- lets a
    ribbon connect the right canonical range to the right alt range,
    bending correctly around an indel instead of assuming the two
    sequences stay in lockstep.

    Each entry (1-based, inclusive coordinates):
        op: 'match' | 'mismatch' | 'deleted' (canonical-only, gap in alt)
            | 'inserted' (alt-only, gap in canonical)
        canon_start, canon_end: None for an 'inserted' segment
        alt_start, alt_end: None for a 'deleted' segment
    """
    if not canonical or not alt:
        return []
    if canonical == alt:
        return [{"op": "match", "canon_start": 1, "canon_end": len(canonical),
                  "alt_start": 1, "alt_end": len(alt)}]

    alignment = _ALIGNER.align(canonical, alt)[0]
    canon_blocks, alt_blocks = alignment.aligned  # 0-based, end-exclusive, paired by index

    segments: list[dict] = []
    c_pos, a_pos = 0, 0
    for (c_start, c_end), (a_start, a_end) in zip(canon_blocks, alt_blocks):
        if c_start > c_pos:
            segments.append({"op": "deleted", "canon_start": c_pos + 1, "canon_end": c_start,
                              "alt_start": None, "alt_end": None})
        if a_start > a_pos:
            segments.append({"op": "inserted", "canon_start": None, "canon_end": None,
                              "alt_start": a_pos + 1, "alt_end": a_start})

        # An ungapped aligned block can still contain mismatches -- split it
        # into match/mismatch sub-runs residue by residue.
        run_op, run_c0, run_a0 = None, None, None
        for i in range(c_end - c_start):
            ci, ai = c_start + i, a_start + i
            op = "match" if canonical[ci] == alt[ai] else "mismatch"
            if op != run_op:
                if run_op is not None:
                    segments.append({"op": run_op, "canon_start": run_c0 + 1, "canon_end": ci,
                                      "alt_start": run_a0 + 1, "alt_end": ai})
                run_op, run_c0, run_a0 = op, ci, ai
        if run_op is not None:
            segments.append({"op": run_op, "canon_start": run_c0 + 1, "canon_end": c_end,
                              "alt_start": run_a0 + 1, "alt_end": a_end})

        c_pos, a_pos = c_end, a_end

    if c_pos < len(canonical):
        segments.append({"op": "deleted", "canon_start": c_pos + 1, "canon_end": len(canonical),
                          "alt_start": None, "alt_end": None})
    if a_pos < len(alt):
        segments.append({"op": "inserted", "canon_start": None, "canon_end": None,
                          "alt_start": a_pos + 1, "alt_end": len(alt)})

    return _merge_adjacent(segments)


def _merge_adjacent(segments: list[dict]) -> list[dict]:
    """Safety net against zero-length or artificially split adjacent runs."""
    merged: list[dict] = []
    for s in segments:
        prev = merged[-1] if merged else None
        contiguous = prev is not None and prev["op"] == s["op"] and (
            (s["canon_start"] is None or prev["canon_end"] + 1 == s["canon_start"])
            and (s["alt_start"] is None or prev["alt_end"] + 1 == s["alt_start"])
        )
        if contiguous:
            if s["canon_end"] is not None:
                prev["canon_end"] = s["canon_end"]
            if s["alt_end"] is not None:
                prev["alt_end"] = s["alt_end"]
        else:
            merged.append(dict(s))
    return merged
