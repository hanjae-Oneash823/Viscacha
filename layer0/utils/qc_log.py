import json
from collections import defaultdict
from pathlib import Path


class QCLogger:
    def __init__(self):
        self._log = {
            "step_01": {"id_harmonization": {}, "missing_fields": {}, "flags": []},
            "step_02": {"unmatched_barcodes": None, "flags": []},
            "step_03": {"dropout_by_ct_donor": {}, "flags": []},
            "step_04": {"prevalence_filter": {}, "flags": []},
            "step_05": {"missing_covariates": {}, "rin_status": "absent_confirmed", "flags": []},
            "step_06": {"excluded_donor_ct": [], "flags": []},
        }

    # --- Step 0.1 ---
    def log_id_harmonization(self, n_smc: int, n_po: int, n_total: int):
        self._log["step_01"]["id_harmonization"] = {
            "n_smc": n_smc, "n_po": n_po, "n_total": n_total,
        }

    def log_missing_metadata(self, donor_id: str, fields: list):
        self._log["step_01"]["missing_fields"][donor_id] = fields
        self.flag("step_01", f"donor {donor_id} missing: {fields}")

    # --- Step 0.2 ---
    def log_unmatched_barcodes(self, n_unmatched: int, n_total: int):
        self._log["step_02"]["unmatched_barcodes"] = {
            "n_unmatched": n_unmatched,
            "n_total": n_total,
            "pct_unmatched": round(n_unmatched / n_total * 100, 2) if n_total else 0,
        }

    # --- Step 0.3 ---
    def log_barcode_drop(self, cell_type: str, donor: str,
                         n_before: int, n_after: int, reason: str,
                         dropout_flag_pct: float = 0.20):
        key = f"{cell_type}|{donor}"
        entry = {
            "n_before": int(n_before),
            "n_after": int(n_after),
            "n_dropped": int(n_before - n_after),
            "pct_dropped": round((1 - n_after / n_before) * 100, 1) if n_before else 0,
            "reason": reason,
        }
        self._log["step_03"]["dropout_by_ct_donor"][key] = entry
        if n_before > 0 and (1 - n_after / n_before) > dropout_flag_pct:
            self.flag("step_03",
                      f"HIGH DROPOUT {cell_type} × {donor}: "
                      f"{entry['pct_dropped']}% dropped ({n_before}→{n_after})")

    # --- Step 0.4 ---
    def log_prevalence_filter(self, cell_type: str, n_before: int, n_after: int):
        self._log["step_04"]["prevalence_filter"][cell_type] = {
            "n_before": n_before,
            "n_after": n_after,
            "n_dropped": n_before - n_after,
            "pct_kept": round(n_after / n_before * 100, 1) if n_before else 0,
        }

    # --- Step 0.5 ---
    def log_covariate_missing(self, donors: list, covariate: str):
        self._log["step_05"]["missing_covariates"][covariate] = donors
        self.flag("step_05", f"covariate '{covariate}' missing for donors: {donors}")

    def log_rin_absent(self):
        msg = (
            "RIN scores absent from all sources (xlsx, pptx, AnnData obs). "
            "pct_counts_mt (from adata_sr.obs, transferred in Step 0.2) is used as "
            "RNA quality proxy. median_pct_mt per donor×cell_type is computed in "
            "Step 0.6 and used as the RNA quality covariate in DRIMSeq."
        )
        self._log["step_05"]["rin_status"] = msg
        self.flag("step_05", msg)

    # --- Step 0.6 ---
    def log_min_cell_exclusion(self, cell_type: str, donor: str, n_cells: int):
        self._log["step_06"]["excluded_donor_ct"].append({
            "cell_type": cell_type,
            "donor": donor,
            "n_cells": int(n_cells),
            "reason": f"below min_cells threshold ({n_cells} cells)",
        })

    # --- General flag ---
    def flag(self, step: str, message: str):
        step_key = step if step in self._log else None
        if step_key:
            self._log[step_key]["flags"].append(message)
        print(f"[QC FLAG] [{step}] {message}")

    # --- Write to disk ---
    def write(self, out_dir: Path):
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        json_path = out_dir / "layer0_qc.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(self._log, f, indent=2, ensure_ascii=False)

        txt_path = out_dir / "layer0_qc_report.txt"
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write("=== Viscacha Layer 0 QC Report ===\n\n")
            for step, data in self._log.items():
                f.write(f"--- {step.upper()} ---\n")
                flags = data.get("flags", [])
                if flags:
                    for msg in flags:
                        f.write(f"  FLAG: {msg}\n")
                else:
                    f.write("  No flags.\n")
                f.write("\n")

        print(f"QC log written to {out_dir}")
        return json_path, txt_path
