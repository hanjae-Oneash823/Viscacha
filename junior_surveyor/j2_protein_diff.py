"""J2 — Protein-level consequence from pairwise alignment.

Uses edlib for exact edit-distance alignment (returns CIGAR string) instead of
Biopython's scoring-based PairwiseAligner. This eliminates the micro-gap
artefact where a zero-net-length domain swap was misclassified as internal_indel
due to the aligner opening 1-2 cheap gap pairs instead of recording mismatches.

With unit-cost edit distance:
  - A k-aa same-length domain swap → k mismatches (cost k)
  - vs k insertions + k deletions   → cost 2k  (always worse)
  So length_diff==0 + any divergence → always pure X ops → always "substitution".

New change types vs original:
    N_extension       — extra sequence prepended to N-terminus  (split from extension)
    C_extension       — extra sequence appended to C-terminus   (split from extension)
    internal_insertion — extra sequence inserted mid-protein    (new)
    substitution      — same-length, residue-level swap         (new)

New columns:
    mismatch_aa_count     — residue–residue mismatches (CIGAR X ops)
    indel_aa_count        — gap positions (CIGAR I + D ops)
    changed_aa_fraction   — unique canonical positions affected / canonical length
"""

from __future__ import annotations

import re
import sys

import edlib
import pandas as pd

from junior_surveyor import j2_pfam

_CIGAR_RE = re.compile(r"(\d+)([=XID])")

# Positions from N-terminus where truncation/extension is called "N-terminal"
_N_TERM_TRUNC = 5
# Positions from N-terminus where an extension is called "N_extension"
_N_EXT_ABS    = 10


def _log(msg: str) -> None:
    print(f"  [j2] {msg}", file=sys.stderr, flush=True)


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    return [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]


def _align_proteins(canonical: str, alt: str) -> dict:
    """Align alt to canonical with edlib NW and classify the protein change.

    CIGAR conventions (edlib, query=canonical):
        = : match            — both advance
        X : mismatch         — both advance; canonical position recorded
        I : insertion in query (canonical has it, alt doesn't)
              → deletion from alt's perspective; canonical position recorded
        D : deletion in query (alt has it, canonical doesn't)
              → insertion in alt; anchor = current canonical pos (clamped ≥ 1)
    """
    empty = {
        "protein_change_type": "no_sequence",
        "protein_length_diff": 0,
        "pct_identity":        0.0,
        "changed_aa_start":    0,
        "changed_aa_end":      0,
        "mismatch_aa_count":   0,
        "indel_aa_count":      0,
        "changed_aa_fraction": 0.0,
        "premature_stop":      False,
    }
    if not canonical or not alt:
        return empty

    premature_stop = "*" in alt[:-1]

    if canonical == alt:
        return {**empty, "protein_change_type": "identical",
                "pct_identity": 1.0, "premature_stop": False}

    length_diff = len(alt) - len(canonical)
    canon_len   = len(canonical)

    result = edlib.align(canonical, alt, mode="NW", task="path")
    ops = _parse_cigar(result["cigar"])

    can_pos        = 0
    match_count    = 0
    mismatch_count = 0
    indel_count    = 0
    changed_positions: list[int] = []

    for count, op in ops:
        if op == "=":
            can_pos     += count
            match_count += count
        elif op == "X":
            for _ in range(count):
                can_pos += 1
                mismatch_count += 1
                changed_positions.append(can_pos)
        elif op == "I":
            # Canonical residues absent from alt (deletion in alt)
            for _ in range(count):
                can_pos += 1
                indel_count += 1
                changed_positions.append(can_pos)
        elif op == "D":
            # Alt residues absent from canonical (insertion in alt)
            for _ in range(count):
                indel_count += 1
                # Anchor to current canonical position (clamped to ≥ 1 so
                # an N-terminal insertion before the first residue is not
                # recorded as position 0, which would spuriously trigger
                # the N_truncation threshold check of ≤ 5).
                changed_positions.append(max(can_pos, 1))

    unique_positions    = sorted(set(changed_positions))
    changed_aa_start    = unique_positions[0]  if unique_positions else 0
    changed_aa_end      = unique_positions[-1] if unique_positions else 0
    changed_aa_fraction = round(len(unique_positions) / canon_len, 4) if canon_len else 0.0

    aligned_pairs = match_count + mismatch_count   # residue–residue pairs only
    pct_identity  = round(match_count / aligned_pairs, 4) if aligned_pairs else 0.0

    # ── Relative threshold for N/C extension boundary (5% of protein, min 10) ──
    n_ext_thresh = max(_N_EXT_ABS, int(canon_len * 0.05))
    c_ext_thresh = n_ext_thresh

    # ── Classifier ────────────────────────────────────────────────────────────
    if premature_stop:
        change_type = "frameshift_stop"

    elif length_diff < 0:
        if not unique_positions:
            change_type = "C_truncation"
        elif changed_aa_start <= _N_TERM_TRUNC:
            change_type = "N_truncation"
        elif changed_aa_end >= canon_len - _N_TERM_TRUNC:
            change_type = "C_truncation"
        else:
            change_type = "internal_indel"

    elif length_diff > 0:
        if not unique_positions:
            # No divergence in shared region → pure C-terminal addition
            change_type = "C_extension"
        elif changed_aa_start <= n_ext_thresh:
            change_type = "N_extension"
        elif changed_aa_end >= canon_len - c_ext_thresh:
            change_type = "C_extension"
        else:
            change_type = "internal_insertion"

    else:
        # length_diff == 0: with edlib unit-cost, same-length swaps always
        # produce X ops (not I/D pairs), so indel_count is 0 here in practice.
        change_type = "substitution" if mismatch_count > 0 else "identical"

    return {
        "protein_change_type": change_type,
        "protein_length_diff": length_diff,
        "pct_identity":        pct_identity,
        "changed_aa_start":    changed_aa_start,
        "changed_aa_end":      changed_aa_end,
        "mismatch_aa_count":   mismatch_count,
        "indel_aa_count":      indel_count,
        "changed_aa_fraction": changed_aa_fraction,
        "premature_stop":      premature_stop,
    }


def _find_affected_domain(row: pd.Series, domain_hits: dict[str, list[dict]]) -> str:
    start = row.get("changed_aa_start", 0)
    end   = row.get("changed_aa_end",   0)
    if start == 0:
        return "none"
    canon_key = f"{row.get('canonical_enst', '')}__canonical"
    overlapping = [
        h["name"] for h in domain_hits.get(canon_key, [])
        if h["start"] <= end and h["end"] >= start
    ]
    return ", ".join(dict.fromkeys(overlapping)) if overlapping else "none"


def _canonical_pfam_domains(row: pd.Series, domain_hits: dict[str, list[dict]]) -> tuple[str, int]:
    """Return (comma-separated domain names, count) for the canonical transcript."""
    key   = f"{row.get('canonical_enst', '')}__canonical"
    names = list(dict.fromkeys(h["name"] for h in domain_hits.get(key, [])))
    return ", ".join(names), len(names)


def _domain_gain_loss(row: pd.Series, domain_hits: dict[str, list[dict]]) -> tuple[str, str]:
    if row.get("protein_change_type", "no_sequence") == "no_sequence":
        return "", ""
    canon_key = f"{row.get('canonical_enst', '')}__canonical"
    alt_key   = f"{row.get('ENST_ID', '')}__alt"
    canon_hits = {h["acc"]: h["name"] for h in domain_hits.get(canon_key, [])}
    alt_hits   = {h["acc"]: h["name"] for h in domain_hits.get(alt_key,   [])}
    lost   = list(dict.fromkeys(n for acc, n in canon_hits.items() if acc not in alt_hits))
    gained = list(dict.fromkeys(n for acc, n in alt_hits.items()   if acc not in canon_hits))
    return ", ".join(lost), ", ".join(gained)


def run(df: pd.DataFrame) -> pd.DataFrame:
    _log(f"aligning proteins for {len(df):,} hits …")

    rows = []
    for i, (_, row) in enumerate(df.iterrows()):
        metrics = _align_proteins(row["canonical_protein_seq"], row["alt_protein_seq"])
        rows.append(metrics)
        if (i + 1) % 200 == 0:
            _log(f"  {i+1:,}/{len(df):,} aligned")

    diff_df = pd.DataFrame(rows, index=df.index)
    out = pd.concat([df, diff_df], axis=1)

    _log("scanning Pfam domains (canonical vs alt) …")
    domain_hits = j2_pfam.run(out)

    out["affected_domain"] = out.apply(
        lambda r: _find_affected_domain(r, domain_hits), axis=1)
    gain_loss = out.apply(
        lambda r: _domain_gain_loss(r, domain_hits), axis=1, result_type="expand")
    out["domains_lost"], out["domains_gained"] = gain_loss[0], gain_loss[1]

    pfam_info = out.apply(
        lambda r: _canonical_pfam_domains(r, domain_hits), axis=1, result_type="expand")
    out["pfam_domains"]   = pfam_info[0]
    out["pfam_n_domains"] = pfam_info[1]

    type_counts = out["protein_change_type"].value_counts()
    _log("protein_change_type distribution:")
    for k, v in type_counts.items():
        _log(f"  {k}: {v:,}")

    return out
