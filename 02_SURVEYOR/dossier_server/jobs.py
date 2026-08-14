"""Background job runner for export requests.

ThreadPoolExecutor(max_workers=1), not Celery/Redis: single user, one GPU,
so GPU-bound work (m2_structures' ColabFold/ESMFold calls) must be
serialized anyway -- a distributed queue would add operational weight with
no throughput benefit here. Status is written to a JSON file per job
(outputs/master_surveyor/jobs/<job_id>.json) so it survives a server
restart and the frontend can poll it with a plain GET.
"""

from __future__ import annotations

import json
import sys
import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from dossier_server.config import DOCKING_DIR, JOBS_DIR

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dossier.data import get_hit_rows
from master_surveyor import m1_ligands, m2_structures, m2b_structure_qc, m3_export

_executor = ThreadPoolExecutor(max_workers=1)
_lock = threading.Lock()


def _job_path(job_id: str) -> Path:
    return JOBS_DIR / f"{job_id}.json"


def _write_status(job_id: str, status: dict) -> None:
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    with _lock:
        _job_path(job_id).write_text(json.dumps(status, indent=2, default=str))


def get_job(job_id: str) -> dict | None:
    path = _job_path(job_id)
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _update_item(job_id: str, status: dict, item_idx: int, **fields) -> None:
    status["items"][item_idx].update(fields)
    status["updated_at"] = time.time()
    _write_status(job_id, status)


def _run_item(hit_row: dict, thresholds_used: dict) -> dict:
    gene, cell_type = hit_row["gene_name"], hit_row["cell_type"]
    transcript = hit_row.get("hit_transcript_name") or hit_row["hit_ENST_ID"]
    drug_names = hit_row.get("selected_drugs") or []

    canonical = m2_structures.get_canonical_structure(
        hit_row.get("uniprot_acc", ""), hit_row["canonical_protein_seq"],
    )
    alt = m2_structures.get_alt_structure(hit_row["alt_protein_seq"])
    confidence = m2b_structure_qc.assess(
        alt, canonical.path,
        hit_row["canonical_protein_seq"], hit_row["alt_protein_seq"],
        int(hit_row.get("changed_aa_start") or 0), int(hit_row.get("changed_aa_end") or 0),
    )

    ligands_sdf = DOCKING_DIR / m3_export.hit_folder_name(gene, cell_type, transcript) / "ligands.sdf"
    ligand_statuses = m1_ligands.build_ligands_sdf(drug_names, ligands_sdf) if drug_names else []

    folder = m3_export.export_hit(
        hit_row, canonical.path, canonical.source, alt.seed_models[0],
        confidence, ligand_statuses, ligands_sdf, thresholds_used,
    )
    return {"folder": str(folder), "structure_confidence": confidence.get("verdict", "unknown")}


def _run_job(job_id: str, cart_items: list[dict], thresholds_used: dict) -> None:
    status = get_job(job_id)
    status["status"] = "running"
    _write_status(job_id, status)

    any_failed = False
    for idx, cart_item in enumerate(cart_items):
        _update_item(job_id, status, idx, stage="loading hit", item_status="running")
        try:
            rows = get_hit_rows(cart_item["gene"], cart_item["cell_type"], cart_item["hit_enst"])
            hit_row = rows.iloc[0].to_dict()
            hit_row["selected_drugs"] = cart_item.get("selected_drugs", [])

            _update_item(job_id, status, idx, stage="ligands + structures + QC + export")
            result = _run_item(hit_row, thresholds_used)

            _update_item(job_id, status, idx, stage="done", item_status="done", result=result)
        except Exception as exc:  # noqa: BLE001 -- one item's failure must not abort the batch
            any_failed = True
            _update_item(job_id, status, idx, stage="failed", item_status="error",
                         error=f"{exc}\n{traceback.format_exc()[-2000:]}")

    status = get_job(job_id)
    status["status"] = "done_with_errors" if any_failed else "done"
    status["finished_at"] = time.time()
    _write_status(job_id, status)


def submit_export_job(cart_items: list[dict], thresholds_used: dict) -> str:
    job_id = uuid.uuid4().hex[:12]
    status = {
        "job_id": job_id,
        "status": "queued",
        "created_at": time.time(),
        "thresholds_used": thresholds_used,
        "items": [
            {"gene": c["gene"], "cell_type": c["cell_type"], "hit_enst": c["hit_enst"],
             "stage": "queued", "item_status": "queued"}
            for c in cart_items
        ],
    }
    _write_status(job_id, status)
    _executor.submit(_run_job, job_id, cart_items, thresholds_used)
    return job_id
