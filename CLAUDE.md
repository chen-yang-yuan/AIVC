# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

AIVC is an AI-based virtual tumor cell model trained on subcellular-resolution Xenium 5K spatial transcriptomics.
Given a tumor cell's gene expression and/or H&E patch (plus spatial context), it will predict environmental stress
scores and the cell's subcellular response (pathway activity in nucleus, cytoplasm, and RNA granules), and support
counterfactual "virtual stress" simulation. The code currently covers data processing and outcome definition only;
model implementation has not started (`code/3_model/` holds foundation-model trial runs).

## Documentation rule

- `CLAUDE.md` holds only the essential project description and operating guides. Keep it short.
- `plans/` holds one plan document per project stage. Current documents:
  - `plans/AIVC_plan.md` — background and rationale, data information, model design and architecture.
  - `plans/data_processing.md` — the data processing and outcome-definition pipeline.
  - `plans/AIVC Model Design Tracker.md` — the model-design brainstorm: current state of each design component,
    open questions, dated decision log.
  - `plans/fig1_exploration.md` — the Fig. 1 exploratory analysis in `code/2_exploration/`: nuclear vs cytoplasmic
    expression of tumor cells in the six Xenium 5K samples, artifact controls, association with the model inputs,
    candidate target gene set.
  - Future stages (model implementation, simulation studies, real data analysis, ...) get their own document.
- Maintenance rule: a plan document always reflects the latest status and decisions. When updating, rewrite the
  affected sections and delete stale content; append a dated entry to the document's **Log** section for locked
  decisions and notable updates. Read the relevant plan document before starting non-trivial work.
- `plans/` also holds source documents (e.g. the original ideas doc) for reference.

## Environment and commands

```bash
make env-create          # conda env AIVC-env from code/utils/env.yaml
make env-update          # update --prune, or create if missing
make push                # DEFAULT GOAL: git add -A, commit "Auto-commit: <timestamp>", push origin main
make data-push-Xenium    # rsync data/Xenium* -> HGCC cluster (~/hulab/projects/AIVC/data)
make out-pull            # rsync HGCC output/ -> local output/
make pull_granule_data   # pull processed_data/granule_adata.h5ad files from HGCC
```

Add `-dry` to rsync targets to preview. Run `make` from the repo root.

Heavy Python steps run on the HGCC SLURM cluster: each `.py` has a sibling `.sh` wrapper that activates
`preprocessing-env` and `cd`s into its own directory; submit with `sbatch <name>.sh` from that directory.
Notebooks run locally in Jupyter. R scripts use `here::i_am()` and resolve paths from the repo root.

There is no build, test suite, or linter.

## Layout and path rules

- Python scripts and notebooks use relative paths `../../data/{dataset}/` and `../../output/{dataset}/`, so the
  working directory must be the script's own folder. Do not run them from the repo root.
- `data/` and `output/` are gitignored except `data/_utils/`. `code/old/` and `code/test/` are ignored scratch.
- `output/2_exploration/<dataset>/` holds the Fig. 1 tables and figures; `<dataset>_quick/` holds development runs.
- Per dataset: `data/<dataset>/raw_data/` (10x outputs) → `intermediate_data/` (cleaned AnnData) →
  `processed_data/` (compartment matrices, granules, pathway scores). `data/merged_data/` holds all datasets combined.

## Pipeline at a glance

Run in numeric order; `data/_utils/` first (shared gene panel, cell-type dictionaries, pathway GMT).

- `code/1_preprocessing/`: clean per dataset → merge datasets → embed → plot.
- `code/2_outcome/`: split nuclear/cytoplasmic expression → ssGSEA → SG detection (mcDETECT) → SG ssGSEA →
  stress scores on the merged tumor object.
- `code/3_model/`: scGPT and UNI trial notebooks only.
- `code/2_exploration/`: Fig. 1 exploration. Light steps in `fig1_exploration.ipynb` (local, `AIVC-env`); heavy
  steps (clustering, DE, association) in `1_heavy.py` + `1_heavy.sh` on HGCC; outputs in `output/2_exploration/`.
  Execution order in its README.

## Must-know conventions

- Tumor cells are always `adata.obs["cell_type_merged"] == "Malignant cell"`.
- Transcript compartments: `in_nucleus == 1` is nuclear; `overlaps_nucleus == 0` is cytoplasmic.
- Every script shares the same six-dataset `settings` dict, but the loop is often narrowed to one dataset while
  iterating. Check which datasets a notebook actually runs before assuming outputs exist for all.
- The compartment matrices in `processed_data/` contain all decoded transcripts (no QV >= 20 filter); per-cell depth
  must be computed from them, not from obs `total_counts` (which includes control probes).
