# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

AIVC is an AI-based virtual tumor cell model trained on subcellular-resolution Xenium 5K spatial transcriptomics.
Given a tumor cell's gene expression and/or H&E patch (plus spatial context), it will predict environmental stress
scores and the cell's subcellular response (pathway activity in nucleus, cytoplasm, and RNA granules), and support
counterfactual "virtual stress" simulation. The design is being reframed as a predictive, interpretable model of the
cytoplasmic activity of a data-selected target gene set (state and open decisions in the Tracker). The code covers
data processing, outcome definition and the completed exploratory stage; model implementation has not started
(`code/3_model/` holds foundation-model trial runs).

## Documentation rule

- `CLAUDE.md` holds only the essential project description and operating guides. Keep it short.
- `plans/` holds one plan document per project stage. Current documents:
  - `plans/AIVC_plan.md` — background and rationale, data information, model design and architecture.
  - `plans/data_processing.md` — the data processing and outcome-definition pipeline.
  - `plans/AIVC Model Design Tracker.md` — the model-design brainstorm: current state of each design component,
    open questions, dated decision log.
  - `plans/fig1_exploration.md` — the Fig. 1 exploratory analysis in `code/2_exploration/`: nuclear vs cytoplasmic
    expression of tumor cells in the six Xenium 5K samples, artifact controls, association with the model inputs,
    first per-gene candidate list.
  - `plans/subtype_exploration.md` — the tumor-cell subtype exploration in `code/1_exploration/1_cytoplasmic/`: the 634-gene target
    set (cytoplasm-enriched highly variable genes) and its sensitivity analysis, per-sample,
    uncorrected pooled and sample-corrected pooled (cytoplasmic and nuclear) subtyping, whether the substates are
    genuine, the sample effect, and the implications for the outcome variable.
  - `plans/stress_substate_exploration.md` — the stress vs target-gene-substate analysis: stress axes and groups
    for the five samples, correspondence of the cytoplasmic and nuclear substates with stress, and whether
    cytoplasmic expression of the target genes carries more stress information than nuclear expression.
  - `plans/ncd_ratio_exploration.md` — the current plan of the nucleocytoplasmic discordance (NCD) analysis:
    NCD as the cytoplasm-to-nucleus ratio per gene and cell, score log2((c+1)/(n+1)), cells = tumor cells with at least one nuclear transcript, protein-coding genes only, gene
    filter by tumor expression, target genes selected from the score, substates and their checks. Steps 0 (counts check), 1 (gene
    filter, N = 3000, 2967 genes), 2 (selection statistics, three candidate lists) and 3 (audit of the lists)
    are run and audited (2026-11-19). Primary target list: V-abs, the 500 filtered genes with the largest variance of the score (checks in step 4:
    the first 1000, and the union with the V-null list). Procedure
    in its sections 3 and 4, decisions and rationale in 8, results in 10. Step 4 (pooled embedding and clustering for
    three gene lists and two per-sample corrections; procedure in its section 5) is implemented; the full run on
    HGCC and the audit of list and correction are pending.
  - `plans/ncd_exploration_stale_residual_model.md` — STALE (2026-11-18). The first NCD analysis (binomial
    residual with a per-cell shift and scale; negative for 634 and 545 genes). The advisor judged the model
    wrong; kept for reference. Its scripts are in `code/1_exploration/old/ncd_residual/`, its outputs and logs
    in `output/1_exploration/old/ncd_residual/`.
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
make out-pull            # rsync HGCC output/ -> local output/ (never overwrites newer local files)
make pull_granule_data   # pull processed_data/granule_adata.h5ad files from HGCC
make ncd-push            # NCD ratio step 4: gene lists -> HGCC;  make ncd-pull: results of 4_embedding -> local
```

Add `-dry` to rsync targets to preview. Run `make` from the repo root.

Heavy Python steps run on the HGCC SLURM cluster: each `.py` has a sibling `.sh` wrapper that activates
`preprocessing-env` and `cd`s into its own directory; submit with `sbatch <name>.sh` from that directory.
Notebooks run locally in Jupyter. R scripts use `here::i_am()` and resolve paths from the repo root.

There is no build, test suite, or linter.

## Layout and path rules

- Python scripts and notebooks use relative paths `../../data/{dataset}/` and `../../output/{dataset}/`
  (`../../../` in the subfolders of `code/1_exploration/`, which sit one level deeper), so the
  working directory must be the script's own folder. Do not run them from the repo root.
- `data/` and `output/` are gitignored except `data/_utils/`. `code/old/` and `code/test/` are ignored scratch;
  every archive folder `code/**/old/` is gitignored too (local and Dropbox only). `code/1_exploration/old/` is
  the archive of the superseded first exploration scripts and, in `ncd_residual/`, of the first NCD analysis
  (do not run).
- `output/2_exploration/<dataset>/` holds the Fig. 1 tables and figures; `<dataset>_quick/` holds development runs.
- `output/1_exploration/` holds the gene set (`enrichment/`), the subtype exploration per sample (`<dataset>/`,
  summary in `cross_sample/`) and pooled (`pooled*/`), the gene-set sensitivity (`gene_set_sensitivity/`), the
  stress analysis (`stress/`), the NCD ratio analysis (`ncd_ratio/`, one numbered folder per step) and the
  archive of the first NCD analysis (`old/ncd_residual/`). A re-run empties the analysis subfolders it writes; the embedded `.h5ad` objects from
  HGCC, `enrichment_results.csv` and the cached stress table are kept.
- Per dataset: `data/<dataset>/raw_data/` (10x outputs) → `intermediate_data/` (cleaned AnnData) →
  `processed_data/` (compartment matrices, granules, pathway scores). `data/merged_data/` holds all datasets combined.

## Pipeline at a glance

Run each folder in the numeric order of its files; `data/_utils/` first (shared gene panel, cell-type
dictionaries, pathway GMT), then `1_preprocessing` → `2_outcome` → the two exploration folders.

- `code/1_preprocessing/`: clean per dataset → merge datasets → embed → plot.
- `code/2_outcome/`: split nuclear/cytoplasmic expression → ssGSEA → SG detection (mcDETECT) → SG ssGSEA →
  stress scores on the merged tumor object.
- `code/1_exploration/`: two analyses of the five carcinoma samples, one subfolder each (map in its `README.md`).
  `subtype_utils.py` at the top is shared; the subfolders add `..` to the import path.
  - `1_cytoplasmic/`: target gene selection and tumor-cell subtype exploration (`1`-`5`; embeddings `2`/`4` as
    py+sh on HGCC, notebooks local) and the stress vs substate analysis (`6_stress_substates.ipynb` with
    `stress_utils.py` / `stress_figures.py`, local, about 1 h).
  - `2_ncd_ratio/`: the NCD analysis with the ratio score (`plans/ncd_ratio_exploration.md`), scripts numbered
    by step: `gene_biotype.py` (once, needs the network), `0_counts.py` (local, about 40 s), `1_filter.py` (seconds), `2_selection.py` (4 min), `3_audit.py` (4.5 min), with
    `ncd_ratio_utils.py`; `check_upstream.py` (4 min) audits the input matrices and is not a step. Step 4 is heavy:
    `4_embedding.py` + `4_embedding.sh` (SLURM array of 24 tasks on HGCC; `make ncd-push` before, `make ncd-pull`
    after; steps in `code/1_exploration/README.md`), then `4_review.py` locally. Later steps are implemented
    one at a time after review.
  - Outputs in `output/1_exploration/`. The scripts of the first NCD analysis (`7` to `13`, `ncd_utils.py`) are
    archived in `old/ncd_residual/`.
- `code/3_model/`: scGPT and UNI trial notebooks only.
- `code/2_exploration/`: Fig. 1 exploration. Light steps in `fig1_exploration.ipynb` (local, `AIVC-env`); heavy
  steps (clustering, DE, association) in `1_heavy.py` + `1_heavy.sh` on HGCC; outputs in `output/2_exploration/`.
  Execution order in its README.

## Must-know conventions

- Tumor cells are always `adata.obs["cell_type_merged"] == "Malignant cell"`.
- Transcript compartments: `in_nucleus == 1` is nuclear; `overlaps_nucleus == 0` is cytoplasmic.
- The preprocessing and outcome scripts share the same six-dataset `settings` dict, but the loop is often narrowed
  to one dataset while iterating. Check which datasets a notebook actually runs before assuming outputs exist for
  all. `code/1_exploration/` uses a fixed list of the five carcinoma samples (Skin, a melanoma, is excluded);
  `code/2_exploration/` uses all six (`config.py`).
- The compartment matrices in `processed_data/` contain all decoded transcripts (no QV >= 20 filter); per-cell depth
  must be computed from them, not from obs `total_counts` (which includes control probes).
