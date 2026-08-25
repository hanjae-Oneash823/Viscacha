#!/usr/bin/env python3
"""Download versioned public inputs for the expanded docking campaign.

The candidate list mixes direct redocking controls, exploratory cross-docking,
and structural exclusions.  This script only stages immutable experimental
structures and reviewed protein sequences; receptor/ligand preparation and
docking are separate, auditable steps.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
SYSTEMS = CAMPAIGN / "systems"

STRUCTURES = {
    "BACE1_verubecestat": {
        "5HU1.pdb": "https://files.rcsb.org/download/5HU1.pdb",
    },
    "CHRNA7_encenicline": {
        "7EKP.pdb": "https://files.rcsb.org/download/7EKP.pdb",
        "9QTO.cif": "https://files.rcsb.org/download/9QTO.cif",
    },
    "GABRA2_AZD7325": {
        "9CXC.cif": "https://files.rcsb.org/download/9CXC.cif",
        "9CSB.cif": "https://files.rcsb.org/download/9CSB.cif",
        "6X3X.pdb": "https://files.rcsb.org/download/6X3X.pdb",
    },
    "CACNA1D_isradipine": {
        "6JP5.pdb": "https://files.rcsb.org/download/6JP5.pdb",
    },
    "PDE9A_BI409306": {
        "4GH6.pdb": "https://files.rcsb.org/download/4GH6.pdb",
    },
}

SEQUENCES = {
    "BACE1_verubecestat": {
        "BACE1_501_P56817-1.fasta": "https://rest.uniprot.org/uniprotkb/P56817-1.fasta",
        "BACE1_476_P56817-2.fasta": "https://rest.uniprot.org/uniprotkb/P56817-2.fasta",
        "BACE1_457_P56817-3.fasta": "https://rest.uniprot.org/uniprotkb/P56817-3.fasta",
    },
    "PDE9A_BI409306": {
        "PDE9A_O76083.fasta": "https://rest.uniprot.org/uniprotkb/O76083.fasta",
    },
}


def download(url: str, destination: Path) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "Viscacha-docking/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    if len(payload) < 100:
        raise RuntimeError(f"Unexpectedly short response from {url}: {len(payload)} bytes")
    destination.write_bytes(payload)
    return {
        "url": url,
        "path": str(destination.relative_to(ROOT)),
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def main() -> None:
    records: list[dict[str, object]] = []
    for candidate, files in STRUCTURES.items():
        for filename, url in files.items():
            record = download(url, SYSTEMS / candidate / "inputs" / filename)
            record.update({"candidate": candidate, "kind": "experimental_structure"})
            records.append(record)
    for candidate, files in SEQUENCES.items():
        for filename, url in files.items():
            record = download(url, SYSTEMS / candidate / "inputs" / filename)
            record.update({"candidate": candidate, "kind": "reviewed_sequence"})
            records.append(record)

    manifest = {
        "purpose": "Public experimental templates and reviewed protein sequences for all-candidate docking triage",
        "records": records,
    }
    output = CAMPAIGN / "analysis" / "metadata" / "expanded_inputs_manifest.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
