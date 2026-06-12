"""
M01 — Sequence Retrieval and Comparator Selection (plan 5, M01).

Database: Ensembl REST.

For each gene candidate this module:
  1. Resolves the MANE Select transcript (mane=1 on the gene lookup).
  2. Selects a structural comparator for each significant transcript, following
     the plan's case rules (paired opposite-sign / single / same-sign / 3+).
  3. Retrieves, per ENST: protein sequence, biotype, exon structure, protein
     length, and the canonical UniProt accession (via the translation xref).
  4. Computes the exon symmetric difference (genomic overlap >= 80% of the
     shorter exon) between each significant transcript and its comparator.
  5. Flags non-coding biotypes and (approximately) frame-shifting exon diffs.

Comparator selection is gene-level: it considers all significant transcripts of
the gene across every cell type.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from layer2.config import CONFIG
from layer2.inputs import GeneCandidate
from layer2.utils.api_client import APIClient, APIError
from layer2.utils.audit import AuditLog

NONCODING_BIOTYPES = {"retained_intron", "nonsense_mediated_decay",
                      "processed_transcript", "non_stop_decay"}


class M01Error(Exception):
    """Blocking M01 failure (no sequence / no comparator) — plan 7."""


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Exon:
    id: str
    start: int          # genomic
    end: int            # genomic
    strand: int
    rank: int           # position within transcript (5'->3')

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class TranscriptRecord:
    enst_id: str
    transcript_name: str
    biotype: str
    protein_seq: str
    protein_length: int
    translation_id: str | None
    uniprot_swissprot: str | None
    uniprot_trembl: str | None
    exons: list[Exon]
    chromosome: str
    strand: int

    @property
    def is_protein_coding(self) -> bool:
        return self.biotype == "protein_coding"

    @property
    def noncoding_flag(self) -> bool:
        return self.biotype in NONCODING_BIOTYPES


@dataclass
class ExonDiff:
    """Symmetric difference of two transcripts' exon sets (by genomic overlap).
    Reports which exons differ, for the M08 isoform-architecture schematic. The
    authoritative frame/coding consequence is taken from the protein-sequence
    comparison on SignificantTranscriptM01, not from exon arithmetic (an
    exon-length proxy false-positives on UTR exons)."""
    unique_to_query: list[Exon]       # exons only in the significant transcript
    unique_to_comparator: list[Exon]
    n_shared: int


@dataclass
class SignificantTranscriptM01:
    """M01 result for one significant transcript and its chosen comparator."""
    transcript: TranscriptRecord
    comparator: TranscriptRecord
    comparator_source: str            # paired_hit | mane_select | highest_ctrl
    exon_diff: ExonDiff
    role: str                         # ad_enriched | control_enriched

    @property
    def protein_identical(self) -> bool:
        """True if the significant transcript and its comparator encode the
        exact same protein. When true there is NO structural docking hypothesis
        — the DTU is regulatory (UTR / non-coding exon usage) only."""
        a, b = self.transcript.protein_seq, self.comparator.protein_seq
        return bool(a) and a == b

    @property
    def protein_length_delta(self) -> int:
        return self.transcript.protein_length - self.comparator.protein_length

    @property
    def has_structural_difference(self) -> bool:
        """Gate for the docking path (consumed by M03/M08): the query must be
        protein-coding, have a sequence, and differ from the comparator."""
        return (self.transcript.is_protein_coding
                and bool(self.transcript.protein_seq)
                and not self.protein_identical)


@dataclass
class M01Result:
    gene_id: str
    ensg_id: str
    mane_select_enst: str | None
    canonical_uniprot: str | None     # gene-level; from MANE Select translation
    significant: list[SignificantTranscriptM01] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Ensembl helpers
# ---------------------------------------------------------------------------

def _ensembl(client: APIClient, path: str, params: dict | None = None):
    url = CONFIG.apis["ensembl"].rstrip("/") + path
    p = {"content-type": "application/json", **(params or {})}
    return client.get_json("ensembl", url, params=p)


def _find_mane_select(client: APIClient, ensg: str) -> tuple[str | None, str | None]:
    """Return (mane_enst, mane_display_name) for a gene, or (None, None)."""
    data = _ensembl(client, f"/lookup/id/{ensg}", {"expand": "1", "mane": "1"})
    for t in data.get("Transcript", []):
        for m in (t.get("MANE") or []):
            if m.get("type") == "MANE_Select":
                return t.get("id"), t.get("display_name")
    # Fallback: Ensembl canonical transcript
    for t in data.get("Transcript", []):
        if t.get("is_canonical"):
            return t.get("id"), t.get("display_name")
    return None, None


def _highest_control_comparator(gene_id: str, cell_type: str, exclude_tx: str,
                                enst_map: pd.DataFrame, config
                                ) -> tuple[str | None, str | None]:
    """Plan M01 fallback #2: the transcript of this gene with the highest mean
    count across CONTROL donors in the significant cell type, excluding the
    significant transcript itself. Returns (enst_id, transcript_name)."""
    counts_path = config.pseudobulk_dir / f"counts_{cell_type}.csv"
    meta_path = config.pseudobulk_dir / f"metadata_{cell_type}.csv"
    if not counts_path.exists() or not meta_path.exists():
        return None, None
    counts = pd.read_csv(counts_path, index_col=0)   # donors x transcripts
    meta = pd.read_csv(meta_path, index_col=0)
    ctrl = [d for d in meta.index[meta["condition"] == config.conditions["control"]]
            if d in counts.index]
    gene_txs = [c for c in counts.columns
                if c.startswith(gene_id + "-") and c != exclude_tx]
    if not ctrl or not gene_txs:
        return None, None
    means = counts.loc[ctrl, gene_txs].mean(axis=0)
    best = str(means.idxmax())
    enst = str(enst_map.at[best, "ENST_ID"]) if best in enst_map.index else None
    return enst, best


def _pick_uniprot(client: APIClient, translation_id: str | None
                  ) -> tuple[str | None, str | None]:
    """Return (swissprot, trembl) accessions for a translation, or (None, None)."""
    if not translation_id:
        return None, None
    try:
        xrefs = _ensembl(client, f"/xrefs/id/{translation_id}")
    except APIError:
        return None, None
    sp = tr = None
    for x in xrefs:
        db = (x.get("dbname") or "").lower()
        if "swissprot" in db and sp is None:
            sp = x.get("primary_id")
        elif "sptrembl" in db and tr is None:
            tr = x.get("primary_id")
    return sp, tr


def retrieve_transcript(client: APIClient, enst: str,
                        transcript_name: str | None = None) -> TranscriptRecord:
    """Fetch protein sequence, biotype, exon structure, and UniProt accession."""
    info = _ensembl(client, f"/lookup/id/{enst}", {"expand": "1"})
    biotype = info.get("biotype", "")
    translation = info.get("Translation") or {}
    translation_id = translation.get("id")
    protein_length = int(translation.get("length") or 0)
    chrom = str(info.get("seq_region_name", ""))
    strand = int(info.get("strand", 0))
    name = transcript_name or info.get("display_name") or enst

    # Exon structure (rank by 5'->3' order, which depends on strand)
    raw_exons = info.get("Exon", [])
    ordered = sorted(raw_exons, key=lambda e: e["start"], reverse=(strand < 0))
    exons = [
        Exon(id=e.get("id", ""), start=int(e["start"]), end=int(e["end"]),
             strand=int(e.get("strand", strand)), rank=i + 1)
        for i, e in enumerate(ordered)
    ]

    # Protein sequence
    protein_seq = ""
    if translation_id:
        seq = _ensembl(client, f"/sequence/id/{enst}", {"type": "protein"})
        protein_seq = seq.get("seq", "") if isinstance(seq, dict) else ""

    sp, tr = _pick_uniprot(client, translation_id)

    return TranscriptRecord(
        enst_id=enst, transcript_name=name, biotype=biotype,
        protein_seq=protein_seq, protein_length=protein_length or len(protein_seq),
        translation_id=translation_id, uniprot_swissprot=sp, uniprot_trembl=tr,
        exons=exons, chromosome=chrom, strand=strand,
    )


# ---------------------------------------------------------------------------
# Exon symmetric difference
# ---------------------------------------------------------------------------

def _overlap_fraction(a: Exon, b: Exon) -> float:
    inter = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    return inter / min(a.length, b.length) if inter > 0 else 0.0


def compute_exon_diff(query: TranscriptRecord, comparator: TranscriptRecord,
                      min_overlap: float) -> ExonDiff:
    """Symmetric difference by genomic overlap >= min_overlap of shorter exon."""
    matched_comp: set[int] = set()
    unique_q: list[Exon] = []
    for qe in query.exons:
        hit = None
        for j, ce in enumerate(comparator.exons):
            if j in matched_comp:
                continue
            if _overlap_fraction(qe, ce) >= min_overlap:
                hit = j
                break
        if hit is None:
            unique_q.append(qe)
        else:
            matched_comp.add(hit)
    unique_c = [ce for j, ce in enumerate(comparator.exons) if j not in matched_comp]
    n_shared = len(query.exons) - len(unique_q)
    return ExonDiff(unique_to_query=unique_q, unique_to_comparator=unique_c,
                    n_shared=n_shared)


# ---------------------------------------------------------------------------
# Comparator selection + main entry
# ---------------------------------------------------------------------------

def run(gene: GeneCandidate, client: APIClient, audit: AuditLog,
        enst_map: pd.DataFrame, config=CONFIG) -> M01Result:
    """Run M01 for one gene. Raises M01Error on a blocking failure (plan 7)."""
    min_overlap = config.thresholds["exon_overlap_min_fraction"]
    warnings: list[str] = []

    # MANE Select (comparator fallback + canonical UniProt source)
    try:
        mane_enst, mane_name = _find_mane_select(client, gene.ensg_id)
    except APIError as e:
        raise M01Error(f"{gene.gene_id}: Ensembl gene lookup failed: {e}")

    # Unique significant transcripts and their ΔPSI signs
    sig_txs = []
    seen = set()
    for t in gene.transcripts:
        if t.transcript_id not in seen:
            seen.add(t.transcript_id)
            sig_txs.append(t)
    signs = {t.transcript_id: (1 if t.delta_psi > 0 else -1) for t in sig_txs}

    # Decide the case (plan M01)
    paired_opposite = (len(sig_txs) == 2 and signs[sig_txs[0].transcript_id]
                       != signs[sig_txs[1].transcript_id])

    # Cache of retrieved TranscriptRecords by ENST (avoid duplicate fetches)
    records: dict[str, TranscriptRecord] = {}

    def get_record(enst: str, name: str | None) -> TranscriptRecord:
        if enst not in records:
            try:
                records[enst] = retrieve_transcript(client, enst, name)
            except APIError as e:
                raise M01Error(f"{gene.gene_id}: sequence retrieval failed for {enst}: {e}")
        return records[enst]

    # Canonical UniProt from MANE Select translation (gene-level, for M02)
    canonical_uniprot = None
    if mane_enst:
        mane_rec = get_record(mane_enst, mane_name)
        canonical_uniprot = mane_rec.uniprot_swissprot or mane_rec.uniprot_trembl
    else:
        warnings.append(f"{gene.gene_id}: no MANE Select or canonical transcript found")

    results: list[SignificantTranscriptM01] = []
    for t in sig_txs:
        q = get_record(t.enst_id, t.transcript_id)

        # Comparator selection (plan M01)
        if paired_opposite:
            partner = next(s for s in sig_txs if s.transcript_id != t.transcript_id)
            comp = get_record(partner.enst_id, partner.transcript_id)
            comp_source = "paired_hit"
        else:
            if len(sig_txs) >= 3:
                audit.info(f"{gene.gene_id}: {len(sig_txs)} significant transcripts "
                           f"— MANE Select used as common reference")
            elif len(sig_txs) == 2:
                audit.warning(gene.gene_id, "two significant transcripts share ΔPSI "
                              "sign — treating each as single-hit vs MANE Select")

            # Use MANE Select unless it IS the significant transcript (then it is
            # not a comparator) — fall through to highest-expressed control tx.
            use_mane = bool(mane_enst) and mane_enst != t.enst_id
            if use_mane:
                comp = get_record(mane_enst, mane_name)
                comp_source = "mane_select"
            else:
                reason = ("significant transcript is itself the MANE Select"
                          if mane_enst == t.enst_id else "no MANE Select")
                fb_enst, fb_name = _highest_control_comparator(
                    gene.gene_id, t.cell_type, t.transcript_id, enst_map, config)
                if not fb_enst:
                    raise M01Error(f"{gene.gene_id}: no comparator available for "
                                   f"{t.transcript_id} ({reason}; no eligible "
                                   f"control transcript in pseudobulk)")
                audit.warning(gene.gene_id, f"{reason} — using highest-expressed "
                              f"control transcript {fb_name} as comparator")
                comp = get_record(fb_enst, fb_name)
                comp_source = "highest_ctrl"

        exon_diff = compute_exon_diff(q, comp, min_overlap)
        if q.noncoding_flag:
            warnings.append(f"{t.transcript_id}: non-coding biotype '{q.biotype}'")

        results.append(SignificantTranscriptM01(
            transcript=q, comparator=comp, comparator_source=comp_source,
            exon_diff=exon_diff, role=t.role,
        ))

    audit.info(f"M01 {gene.gene_id}: {len(results)} significant tx, "
               f"comparator source(s): {sorted({r.comparator_source for r in results})}, "
               f"canonical UniProt: {canonical_uniprot}")

    return M01Result(
        gene_id=gene.gene_id, ensg_id=gene.ensg_id,
        mane_select_enst=mane_enst, canonical_uniprot=canonical_uniprot,
        significant=results, warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Standalone verification
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from layer2.inputs import load_inputs

    audit = AuditLog(CONFIG.audit_log_path, verbose=False)
    client = APIClient(audit)
    inp = load_inputs()

    print(f"{'gene':9} {'role':16} {'tx':14} {'biotype':16} {'aa':>5} "
          f"{'comp':14} {'src':12} {'uniq_exon':>9} {'dAA':>5} {'docking?':>9}")
    print("-" * 120)
    for gene in inp.gene_work_list:
        try:
            res = run(gene, client, audit, inp.enst_map)
        except M01Error as e:
            print(f"{gene.gene_id:9} M01Error: {e}")
            continue
        for s in res.significant:
            verdict = "candidate" if s.has_structural_difference else "TYPE-N"
            print(f"{res.gene_id:9} {s.role:16} {s.transcript.transcript_name:14} "
                  f"{s.transcript.biotype:16} {s.transcript.protein_length:5} "
                  f"{s.comparator.transcript_name:14} {s.comparator_source:12} "
                  f"{len(s.exon_diff.unique_to_query):9} "
                  f"{s.protein_length_delta:+5} {verdict:>9}")
