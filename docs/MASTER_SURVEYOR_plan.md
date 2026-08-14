# MASTER_SURVEYOR — final filtering + BIOVIA Discovery Studio export

## Context

junior_surveyor's J4 gate tags every passing hit with `master_group` (`trial_failure_candidate` / `drug_repurposing_candidate` / `novel_target_candidate`) in `outputs/junior_surveyor/hits_deep.csv`. master_surveyor is the next and final pipeline stage: it must (1) apply a last, group-specific filter, then (2) produce actual 3D protein structures and ligand files, packaged per-hit so they can be opened directly in BIOVIA Discovery Studio for docking.

Decisions already made with the user:
- **Scope**: only `trial_failure_candidate` (43 hits) and `drug_repurposing_candidate` (175 hits) go through this stage — 218 hits total. `novel_target_candidate` (209 hits) is excluded entirely for now (stays in `hits_deep.csv` for a later effort).
- **Filtering rule**: per-group cutoffs are deliberately deferred ("think about it later") — build the mechanism (config-driven, per-group, easy to tighten) but leave the actual thresholds permissive/off so nothing is dropped yet.
- **Protein structures**: predict locally via **ColabFold** (not blind reliance on AlphaFold DB for everything).
- **Ligands**: only for trial_failure/drug_repurposing hits (matches scope above); skip external-library screening for novel_target — not needed since that group is excluded.
- **Deliverable**: per-hit folder containing `protein.pdb` + `ligands.sdf` + a manifest, ready for manual Discovery Studio import.

Infra findings (verified 2026-07-14, read-only checks against the actual box):
- GPU0 has ~22GB free VRAM (GPU1 is busy with another job — must pin `CUDA_VISIBLE_DEVICES=0`).
- Root disk (`/`) has only 39GB free — full AlphaFold2 genetic databases (multi-TB) won't fit. `/home` has 1.5TB free.
- No AlphaFold2/ColabFold/JAX is installed anywhere on the box today.
- `localcolabfold`'s installer script is reachable and installs a self-contained conda env + AF2 params (~a few GB) using ColabFold's remote MMseqs2 API for the MSA step (no local genetic-database download needed) — this is the structure-prediction path that fits the available disk/compute budget.
- `rdkit` (2026.03.3) and `biopython` (1.87) are already installed in the `oneash_dtu` env — no new install needed for ligand conformer generation.
- AlphaFold DB has a live REST API (`alphafold.ebi.ac.uk/api/prediction/<UniProt acc>`) that returns a direct PDB download link when a canonical-sequence structure already exists — free, instant, no GPU. Since the *canonical* protein sequence is (by construction) the plain UniProt/MANE sequence, this covers most canonical structures without spending any local GPU time; only the *alt* isoform sequences (which are novel by definition) actually need local ColabFold prediction. Canonical structures that AFDB doesn't have (no match, obsolete entry, etc.) fall back to local ColabFold too.
- `hits_deep.csv` already carries everything needed as raw material — `alt_protein_seq`, `canonical_protein_seq`, `uniprot_acc`, `chembl_target_id`, `drug_names`, `ot_drug_names`, `affected_domain`, `protein_change_type`, etc. No new sequence-fetching module is needed; only SMILES/ligand-structure fetching is new.
- A known pre-existing data quirk: gene `TBC1D14`/`Oligodendrocyte` has two distinct DTU-significant transcripts both classified `new_target_candidate` → `drug_repurposing_candidate`, i.e. one (gene, cell_type) hit can have 2 rows. Dedup must key on `(gene_name, cell_type, alt_ENST_ID)`, not just `(gene_name, cell_type)`.

## Implementation

New package `02_SURVEYOR/master_surveyor/`, following the same `m0_/m1_/...` step-file convention as `assistant_surveyor` (`l1_`...) and `junior_surveyor` (`j1_`...), each a `run(df) -> df`-style module wired by a `run_master_surveyor.py` orchestrator (mirrors `junior_surveyor/run_junior_surveyor.py`).

### `m0_select.py` — final filter + collapse to one row per hit
- Load `hits_deep.csv`, restrict to `master_group.isin(["trial_failure_candidate", "drug_repurposing_candidate"])`.
- Collapse trial_failure's long-format multi-`alt_rank` rows to one representative row per hit — reuse the existing `_tf_representative_row` logic currently living in `plot_results.py:142` (pick the gate-driving alt with the largest `alt_usage_delta`); move it into this module so both `plot_results.py` and the new export pipeline share one implementation instead of duplicating it.
- Dedup key is `(gene_name, cell_type, alt_ENST_ID)` to correctly keep TBC1D14's two distinct transcripts as two separate hits rather than colliding them.
- Apply group-specific filter predicates from a config block in `master_surveyor/config.py`, e.g.:
  ```python
  # Deliberately permissive placeholders — tighten later, no code changes needed.
  TF_MIN_ABS_DELTA_USAGE = 0.0       # trial_failure: |delta_usage| floor
  TF_REQUIRE_DOMAIN_OVERLAP = False  # trial_failure: affected_domain != "none"
  DR_MIN_CHEMBL_OR_OT_PHASE = 0      # drug_repurposing: clinical-phase floor
  DR_REQUIRE_STRUCTURAL_CHANGE = False  # drug_repurposing: exclude "substitution"-only changes
  ```
  Each flag is applied as a simple boolean mask in `run()`, logged (`kept X/Y after <rule>`) the same way `j4_gate.py` logs branch counts — so raising a threshold later is a one-line config edit, not new code.
- Output: `outputs/master_surveyor/shortlist.csv` (one row per hit, ~218 rows today).

### `m1_ligands.py` — SMILES + 3D conformers for drug candidates
- For each hit's gene, gather candidate drug names already resolved by J3 (`drug_names`, `ot_drug_names`) plus their ChEMBL target id (`chembl_target_id`).
- Extend the ChEMBL molecule fetch (same endpoint/pattern as `junior_surveyor/j3_drug_targets.py:_chembl_drugs_for_target`, adding `canonical_smiles` to the `fields=` query param) to pull SMILES for the named molecules; cache to `outputs/master_surveyor/cache/m1_chembl_smiles.json` following the existing `CACHE_DIR / "j3_*.json"` cache convention.
- For hits with drug *names* but where that pass didn't resolve a SMILES (e.g. Pharos-only or DGIdb-only evidence with no ChEMBL mechanism), fall back to a ChEMBL molecule name search (`/molecule?pref_name__iexact=...`) as a second-tier lookup.
- Generate one 3D conformer per ligand via RDKit ETKDG (`AllChem.EmbedMolecule` + `MMFFOptimizeMolecule`), matching the approach used previously in the retired `layer2/m04_drug.py`.
- Output: per-gene ligand records cached; consumed by `m3_export.py`.

### `m2_structures.py` — protein structures via AlphaFold DB + local ColabFold
- One-time setup (documented as a manual setup note, run once, not by the pipeline script): install `localcolabfold` under `/home/welcome3/tools/localcolabfold` via its official installer, pin `CUDA_VISIBLE_DEVICES=0`.
- Canonical structure: query the AlphaFold DB API (`alphafold.ebi.ac.uk/api/prediction/<uniprot_acc>`) by the hit's `uniprot_acc`; download the returned PDB URL if present and cache it.
- Alt structure (always) and canonical structure (only on an AFDB miss): write a FASTA per unique sequence to `outputs/master_surveyor/cache/fasta/`, dedup across hits that share an identical sequence (e.g. same alt transcript hit in multiple cell types), then batch-invoke `colabfold_batch` (MMseqs2-API MSA mode, no local database) over the FASTA directory. ColabFold already skips sequences whose output exists, so reruns are incremental/cheap.
- Output: `outputs/master_surveyor/cache/structures/<seq_hash>.pdb`, keyed by sequence hash so identical sequences are computed once and reused across every hit that needs them.

### `m3_export.py` — per-hit BIOVIA-ready folders
- For each row in `shortlist.csv`, create `outputs/master_surveyor/docking/<gene>_<cell_type>_<transcript>/` containing:
  - `protein_alt.pdb` and `protein_canonical.pdb` (copied from `m2`'s hash-keyed cache)
  - `ligands.sdf` (multi-molecule SD file, one entry per candidate drug from `m1`, `_Name` property set to the drug name, SMILES/source/evidence-tier as SDF tags)
  - `manifest.json` — gene, cell_type, transcript, `master_group`, `protein_change_type`, `affected_domain`, `delta_usage`, `chi_padj`, drug names + evidence tier/source per ligand, `uniprot_acc`, `chembl_target_id`
- Print a final summary (hits exported, hits missing a structure/ligand and why) the same way other stages log pass/drop counts.

### `run_master_surveyor.py`
Wires `m0 -> m1 -> m2 -> m3` in sequence (mirrors `junior_surveyor/run_junior_surveyor.py`'s structure), then keeps the existing `plot_results.py` visualization suite as a separate, still-independent step (unchanged — it already reads `hits_deep.csv` + `selected_for_next_stage` directly and doesn't need `shortlist.csv`).

### Config additions (`master_surveyor/config.py`)
Add: `COLABFOLD_BIN` path, `CUDA_VISIBLE_DEVICES` pin, `AFDB_API`, `outputs/master_surveyor/{cache,docking}` dirs, the `m0` filter-threshold constants above, RDKit ETKDG params (seed, optimize iterations).

### `rerun_surveyor_pipeline.sh`
Append the new `run_master_surveyor.py` call after the existing `junior_surveyor` step (before the current `plot_results.py` line, which stays).

## Verification
- Run `m0_select.py` alone first and confirm `shortlist.csv` has 218 rows (43 TF + 175 DR) with no TBC1D14-style duplicate collapse regression.
- Run `m1_ligands.py` on a small subset (e.g. 5 genes) and manually inspect a generated `.sdf` in RDKit (`Chem.SDMolSupplier`) to confirm valid 3D coordinates.
- Run `m2_structures.py` on 2-3 sequences first (one AFDB hit, one AFDB miss) to confirm both the AFDB-download path and the local ColabFold path produce a valid PDB before batching the full 218-hit run (which is the long pole — budget real wall-clock time here and checkpoint via ColabFold's own resume behavior).
- Run `m3_export.py` and manually open one resulting folder's `protein_alt.pdb` + `ligands.sdf` in Discovery Studio (or PyMOL/RDKit as a sanity substitute if DS isn't available in this environment) to confirm the files actually load.
- Confirm `rerun_surveyor_pipeline.sh` still runs end-to-end after the new step is inserted.
