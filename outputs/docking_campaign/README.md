# Docking campaign

This directory contains the complete AutoDock Vina candidate campaign. All CPU-bound docking jobs were run serially with a maximum allocation of 16 CPU cores.

## Start here

- Scientific conclusions: [`docs/KEY_FINDINGS.md`](docs/KEY_FINDINGS.md)
- Detailed workflow: [`docs/ANALYSIS_PROCESS.md`](docs/ANALYSIS_PROCESS.md)
- GNINA comparison: [`docs/GNINA_COMPARISON.md`](docs/GNINA_COMPARISON.md)
- Aggregate results: [`analysis/aggregate/expanded_summary.json`](analysis/aggregate/expanded_summary.json)
- GNINA aggregate: [`analysis/gnina/summary.json`](analysis/gnina/summary.json)
- Candidate disposition table: [`analysis/aggregate/candidate_status.csv`](analysis/aggregate/candidate_status.csv)
- Presentation figures: [`figures/expanded_campaign/README.md`](figures/expanded_campaign/README.md)
- GNINA presentation figures: [`figures/gnina_comparison/README.md`](figures/gnina_comparison/README.md)

## Directory layout

```text
docking_campaign/
├── README.md                 navigation and folder conventions
├── docs/
│   ├── KEY_FINDINGS.md       current scientific interpretation
│   ├── ANALYSIS_PROCESS.md   current detailed workflow
│   ├── GNINA_COMPARISON.md   CNN rescoring methods and findings
│   └── archive/              superseded preliminary reports
├── analysis/
│   ├── aggregate/            campaign-level CSV and JSON summaries
│   ├── gnina/                GNINA pose-level and run-level comparisons
│   └── metadata/             shared input and interface provenance
├── figures/
│   ├── expanded_campaign/    current presentation-ready plots and 3D renders
│   ├── gnina_comparison/     Vina/GNINA consensus and disagreement figures
│   ├── preliminary/          earlier standalone figures
│   └── presentation/         presentation-layout assets and exports
├── systems/                    one self-contained directory per protein–drug system
└── tools/                      external executables not tracked by Git
```

## Per-system convention

Each directory under `systems/` retains the original audit structure where applicable:

- `inputs/`: downloaded experimental structures and reviewed sequences;
- `models/`: predicted alternate structures and model-quality outputs;
- `prepared/`: docking-ready receptors, ligands, and staging metadata;
- `runs/`: raw docked poses plus one `result.json` per run;
- `gnina_rescoring/`: pose-preserving CNN rescoring inputs, outputs, logs, and QC;
- `logs/`: execution logs;
- `analysis/` and `figures/`: system-specific summaries and legacy figures.

Do not compare absolute Vina scores between different protein–drug systems. The matched comparisons and their limitations are documented in `docs/KEY_FINDINGS.md`.
