#!/usr/bin/env python3
"""Summarize masitinib replicate docking in the experimental c-KIT 1T46 pocket."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean, median, stdev


def main() -> None:
    campaign = Path("outputs/docking_campaign/systems/KIT_masitinib")
    rows = []
    for path in sorted((campaign / "runs").glob("masitinib_1T46_*/result.json")):
        result = json.loads(path.read_text())
        top = result["results"][0]
        rows.append({"run": path.parent.name, "seed": result["seed"], "exhaustiveness": result["exhaustiveness"], "vina_affinity_kcal_mol": top["vina_affinity_kcal_mol"]})
    if not rows:
        raise SystemExit("No experimental-template masitinib results found")
    analysis = campaign / "analysis"
    analysis.mkdir(exist_ok=True)
    with (analysis / "masitinib_1T46_replicates.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    scores = [row["vina_affinity_kcal_mol"] for row in rows]
    summary = {
        "pocket_template": "experimental c-KIT kinase domain PDB 1T46 (imatinib-bound)",
        "ligand": "masitinib",
        "n_independent_vina_runs": len(rows),
        "score_mean_kcal_mol": round(mean(scores), 3),
        "score_median_kcal_mol": round(median(scores), 3),
        "score_sd_kcal_mol": round(stdev(scores), 3) if len(scores) > 1 else None,
        "score_range_kcal_mol": [round(min(scores), 3), round(max(scores), 3)],
        "scope": "Canonical c-KIT only; this is not a canonical-versus-KIT-223 affinity comparison.",
        "important_limitation": "The KIT-223 deletion site (canonical residue 715) lies in a kinase-insert segment unresolved in PDB 1T46. A comparative result requires a separately rebuilt and validated loop model.",
    }
    (analysis / "masitinib_1T46_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (analysis / "masitinib_1T46_preliminary.md").write_text(
        "# c-KIT–masitinib preliminary docking\n\n"
        f"Masitinib was docked {len(rows)} times into the validated experimental c-KIT 1T46 ATP pocket. "
        f"The Vina score distribution was {summary['score_mean_kcal_mol']:.2f} ± {summary['score_sd_kcal_mol']:.2f} kcal/mol (mean ± SD); "
        f"median {summary['score_median_kcal_mol']:.2f} kcal/mol.\n\n"
        "This supports a reproducible canonical c-KIT docking setup only. It must not be interpreted as an effect of KIT-223, because residue 715 is in an unresolved kinase-insert segment of 1T46.\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
