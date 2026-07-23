"""J0 — Novel transcript ORF prediction via TransDecoder.

For novel IsoQuant transcripts (*.nic / *.nnic) with no CDS annotation in the
Ensembl CDS FASTA, reconstructs full transcript sequences from the genome + GTF
and predicts the best ORF per transcript using TransDecoder with Pfam support.

Pipeline (one-time; result cached as JSON):
  1. Filter GTF to target transcript IDs
  2. gffread  : exons + genome → spliced transcript FASTA
  3. TransDecoder.LongOrfs  : enumerate all ORFs ≥ MIN_PROTEIN_LEN aa
  4. hmmscan  : scan LongOrfs peptides against Pfam-A (guides ORF selection)
  5. TransDecoder.Predict   : choose best ORF per transcript
  6. Parse .transdecoder.cds → {transcript_name: cds_sequence}

Cache: CACHE_DIR / "j0_novel_orfs.json"
Intermediates: CACHE_DIR / "transdecoder/" (kept for inspection)
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

from junior_surveyor.config import (
    CACHE_DIR, GFFREAD_BIN, HMMSCAN_BIN, PFAM_HMM, REPO_ROOT,
    TRANSDECODER_LONGORFS, TRANSDECODER_PREDICT,
)

_GTF    = REPO_ROOT / "outputs/layer1_1/annotation/extended_annotation.filtered.gtf"
_GENOME = REPO_ROOT / "outputs/layer1_1/reference/GRCh38.primary_assembly.genome.fa"
_CACHE  = CACHE_DIR / "j0_novel_orfs.json"
_WORK   = CACHE_DIR / "transdecoder"

MIN_PROTEIN_LEN = 100   # aa
_PFAM_CPU       = 4
_PFAM_EVALUE    = "1e-5"

_TX_ID_RE  = re.compile(r'transcript_id "([^"]+)"')
_P_SUFFIX  = re.compile(r"\.p\d+$")


def _log(msg: str) -> None:
    print(f"  [j0] {msg}", file=sys.stderr, flush=True)


def _check_tools() -> None:
    tools = {
        "gffread":               Path(GFFREAD_BIN),
        "TransDecoder.LongOrfs": TRANSDECODER_LONGORFS,
        "TransDecoder.Predict":  TRANSDECODER_PREDICT,
    }
    missing = [name for name, path in tools.items() if not path.exists()]
    if missing:
        raise RuntimeError(
            f"[j0] missing tools: {missing}\n"
            "Install with:\n"
            "  conda install -n oneash_dtu -c bioconda gffread transdecoder"
        )


# ---------------------------------------------------------------------------
# Step 1 — filter GTF
# ---------------------------------------------------------------------------

def _filter_gtf(transcript_ids: set[str], out_gtf: Path) -> None:
    written = 0
    with open(_GTF) as fh, open(out_gtf, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            m = _TX_ID_RE.search(line)
            if m and m.group(1) in transcript_ids:
                out.write(line)
                written += 1
    _log(f"GTF filtered: {written:,} rows → {len(transcript_ids)} transcripts")


# ---------------------------------------------------------------------------
# Step 2 — gffread: splice exons → transcript FASTA
# ---------------------------------------------------------------------------

def _extract_sequences(gtf: Path, out_fa: Path) -> None:
    cmd = [str(GFFREAD_BIN), str(gtf), "-g", str(_GENOME), "-w", str(out_fa)]
    _log("gffread: extracting spliced transcript sequences …")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"gffread failed:\n{r.stderr}")
    n = sum(1 for line in open(out_fa) if line.startswith(">"))
    _log(f"  {n} transcript sequences written → {out_fa.name}")


# ---------------------------------------------------------------------------
# Step 3 — TransDecoder.LongOrfs
# ---------------------------------------------------------------------------

def _run_long_orfs(fasta: Path) -> Path:
    cmd = [
        str(TRANSDECODER_LONGORFS),
        "-t", str(fasta),
        "-m", str(MIN_PROTEIN_LEN),
    ]
    _log(f"TransDecoder.LongOrfs (min {MIN_PROTEIN_LEN} aa) …")
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_WORK))
    if r.returncode != 0:
        raise RuntimeError(f"TransDecoder.LongOrfs failed:\n{r.stderr}")
    pep = _WORK / f"{fasta.name}.transdecoder_dir" / "longest_orfs.pep"
    n = sum(1 for line in open(pep) if line.startswith(">"))
    _log(f"  {n} candidate ORFs found")
    return pep


# ---------------------------------------------------------------------------
# Step 4 — hmmscan: Pfam support for ORF selection
# ---------------------------------------------------------------------------

def _run_pfam_scan(pep_file: Path) -> Path:
    domtblout = _WORK / "pfam.domtblout"
    cmd = [
        str(HMMSCAN_BIN),
        "--domtblout", str(domtblout),
        "--cpu", str(_PFAM_CPU),
        "-E", _PFAM_EVALUE,
        "--noali",
        str(PFAM_HMM),
        str(pep_file),
    ]
    _log("hmmscan: scanning LongOrfs peptides against Pfam-A …")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"hmmscan failed:\n{r.stderr}")
    n = sum(1 for line in open(domtblout) if not line.startswith("#"))
    _log(f"  {n} Pfam domain hits found")
    return domtblout


# ---------------------------------------------------------------------------
# Step 5 — TransDecoder.Predict
# ---------------------------------------------------------------------------

def _run_predict(fasta: Path, domtblout: Path) -> Path:
    cmd = [
        str(TRANSDECODER_PREDICT),
        "-t", str(fasta),
        "--retain_pfam_hits", str(domtblout),
        "--single_best_only",
    ]
    _log("TransDecoder.Predict: selecting best ORF per transcript …")
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_WORK))
    if r.returncode != 0:
        raise RuntimeError(f"TransDecoder.Predict failed:\n{r.stderr}")
    cds_file = _WORK / f"{fasta.name}.transdecoder.cds"
    return cds_file


# ---------------------------------------------------------------------------
# Step 6 — parse .transdecoder.cds
# ---------------------------------------------------------------------------

def _parse_cds(cds_file: Path) -> dict[str, str]:
    """Return {transcript_name: cds_sequence}, keeping only the top prediction.

    Header example:
      >transcript39631.chr6.nic.p1  type:complete len:150 (+) ...
    The .p1 suffix marks the top-ranked ORF; we keep only that one per transcript.
    """
    result: dict[str, str] = {}
    cur_tx: str | None = None
    parts: list[str] = []

    def _commit() -> None:
        if cur_tx is not None and parts and cur_tx not in result:
            result[cur_tx] = "".join(parts)

    with open(cds_file) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                _commit()
                parts = []
                raw = line[1:].split()[0]           # e.g. transcript39631.chr6.nic.p1
                cur_tx = _P_SUFFIX.sub("", raw)     # → transcript39631.chr6.nic
            else:
                parts.append(line)
    _commit()

    _log(f"parsed {len(result)} CDS predictions")
    return result


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run(transcript_ids: set[str]) -> dict[str, str]:
    """Return {transcript_name: cds_sequence} for *transcript_ids*.

    Loads from cache if available. On first run, executes the full
    gffread → TransDecoder → Pfam pipeline and caches the result.

    Transcripts with no predicted ORF are absent from the returned dict
    (they will remain no_sequence in J2).
    """
    if _CACHE.exists():
        data: dict[str, str] = json.loads(_CACHE.read_text())
        matched = {k: v for k, v in data.items() if k in transcript_ids}
        _log(f"cache hit: {len(matched)}/{len(transcript_ids)} novel CDS sequences loaded")
        return matched

    _check_tools()
    _WORK.mkdir(parents=True, exist_ok=True)
    _log(f"running full pipeline for {len(transcript_ids)} novel transcripts")

    filtered_gtf = _WORK / "novel_transcripts.gtf"
    _filter_gtf(transcript_ids, filtered_gtf)

    tx_fasta = _WORK / "novel_transcripts.fa"
    _extract_sequences(filtered_gtf, tx_fasta)

    pep_file  = _run_long_orfs(tx_fasta)
    domtblout = _run_pfam_scan(pep_file)
    cds_file  = _run_predict(tx_fasta, domtblout)

    orfs = _parse_cds(cds_file)

    _CACHE.write_text(json.dumps(orfs, indent=2))
    _log(f"cached {len(orfs)} predictions → {_CACHE.name}")

    no_orf = transcript_ids - orfs.keys()
    if no_orf:
        _log(f"  {len(no_orf)} transcripts had no ORF ≥ {MIN_PROTEIN_LEN} aa (will stay no_sequence)")

    return {k: v for k, v in orfs.items() if k in transcript_ids}
