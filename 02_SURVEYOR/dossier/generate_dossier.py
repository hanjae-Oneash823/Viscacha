"""DOSSIER — CLI entry point.

Usage:
    python generate_dossier.py --gene STAT3 --cell-type Astrocyte

A (gene, cell_type) pair occasionally carries two distinct hits (e.g. a
trial_failure hit and a separate drug_repurposing hit on different
transcripts) -- if so this exits with an error listing the hit_ENST_ID
options; pass --hit-enst to pick one.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dossier.config import OUT_DIR
from dossier.render import render


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--cell-type", required=True)
    ap.add_argument("--hit-enst", default=None, help="disambiguate when the (gene, cell_type) pair has >1 hit")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    html_out = render(args.gene, args.cell_type, args.hit_enst)

    suffix = f"_{args.hit_enst}" if args.hit_enst else ""
    out_path = OUT_DIR / f"{args.gene}_{args.cell_type}{suffix}.html"
    out_path.write_text(html_out)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
