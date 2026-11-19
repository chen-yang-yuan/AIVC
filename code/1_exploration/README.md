# code/1_exploration

Exploratory analyses of the tumor cells of the five carcinoma samples (BC, OC, CC, LC, Prostate). One subfolder
per analysis; outputs in `output/1_exploration/`.

| Folder | Analysis | Plan document |
|---|---|---|
| `1_cytoplasmic/` | Cytoplasmic expression: target gene set, tumor-cell substates, stress | `plans/subtype_exploration.md`, `plans/stress_substate_exploration.md` |
| `2_ncd_ratio/` | Nucleocytoplasmic discordance (NCD) as the ratio score | `plans/ncd_ratio_exploration.md` |
| `old/` | Archive (first exploration design, first NCD analysis); gitignored, do not run | |

`subtype_utils.py` is shared by the two analyses (gene programs, substate labels and colours, GMT reader). The
scripts and notebooks of the subfolders add `..` to the import path to load it.

## 1_cytoplasmic

Run order: `1_enrichment.ipynb` -> `2_embedding.py` -> `3_clustering.ipynb` -> `4_pooled_embedding.py` (three
variants) -> `4_pooled.ipynb`, `4_pooled_sample_corrected.ipynb`, `4_pooled_sample_corrected_nuclear.ipynb` ->
`5_target_gene_sensitivity.ipynb` (independent) -> `6_stress_substates.ipynb` (with `stress_utils.py` and
`stress_figures.py`). Notebooks run locally in `AIVC-env`; `.py` + `.sh` on HGCC.

## 2_ncd_ratio

Scripts are numbered by the step of the plan document and write to the output folder of the same number:
`0_counts.py` -> `output/1_exploration/ncd_ratio/0_counts/`, `1_filter.py` -> `1_filter/`, `2_selection.py` -> `2_selection/`, `3_audit.py` -> `3_audit/`. Each step reads the
output of the step before, so run them in order. `ncd_ratio_utils.py` is shared by the steps (constants, score, null, selection statistics, ranking rules).
`check_upstream.py` is not a step: it audits the compartment matrices the steps read against the raw 10x transcripts
(`upstream_check/`).
`gene_biotype.py` is run once before the steps (it needs the network) and writes `data/_utils/gene_biotype.csv`.
The steps load the counts through `ncd_ratio_utils.load_counts`, which keeps the tumor cells with at least one nuclear
transcript and the protein-coding genes. Steps 0 to 3 fix the genes: 2967 filtered genes; primary target list = the
first 500 genes of the V-abs ranking (`ncd_ratio_utils.load_target_genes`).

Step 4 (pooled embedding and clustering) is the only heavy step and runs on HGCC as a SLURM array of 24 tasks
(three gene lists x two corrections x data, null and two transcript halves):

| # | Where | Command |
|---|---|---|
| 1 | Local, repo root | `make push` |
| 2 | HGCC | `cd ~/hulab/projects/AIVC && git pull` |
| 3 | Local, repo root | `make ncd-push-dry`, then `make ncd-push` (the gene lists of steps 1 and 3) |
| 4 | HGCC | `cd ~/hulab/projects/AIVC/code/1_exploration/2_ncd_ratio && sbatch 4_embedding.sh` |
| 5 | HGCC | `squeue -u $USER`; `python3 4_embedding.py --status` (in `preprocessing-env`) |
| 6 | HGCC, if a task failed | `sbatch --array=<ids> 4_embedding.sh` (complete tasks are skipped) |
| 7 | Local, repo root | `make ncd-pull-dry`, then `make ncd-pull` (about 5 GB into `output/1_exploration/ncd_ratio/4_embedding/`) |
| 8 | Local, this folder | `python 4_review.py` (comparison tables and diagnostic figures in `4_embedding/review/`) |

`python 4_embedding.py --all --quick --quick-cells 4000` and `python 4_review.py --quick` are the local
development check (subsample; folder `4_embedding_quick/`, delete afterwards).

## Rules

- Run every script and notebook from its own folder: the paths are relative (`../../../data/`,
  `../../../output/`).
- Logs go to `logs/` inside each subfolder (gitignored).
