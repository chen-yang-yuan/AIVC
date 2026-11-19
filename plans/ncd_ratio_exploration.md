# Nucleocytoplasmic discordance (NCD) as a ratio score: plan

Status 2026-11-19: steps 0 to 3 are run, reviewed and audited (last audit 2026-11-19, no error in any
result), and the input matrices are audited against the raw data (§10.8). Two restrictions were added after the audits, and steps 0 to 3 rerun each time: the tumor
cells without nuclear transcripts are removed, and only protein-coding genes are used. All results below are
from the last rerun. Current state of the genes: gene filter with N = 3000 (2967 genes); target list V-abs
(filtered genes ranked by the variance of the score across cells), primary list = the first 500 genes, the
first 1000 as a nested check and the union with the V-null first 500 (630 genes) as a second check. Step 4,
the pooled embedding and clustering, is implemented and tested on a subsample (§5); the full run on HGCC is
submitted by the user, and its results are not in yet. The exact procedure of the gene selection is in §3
and §4, the decisions and their reasons in §8, the results in §10. This document replaces
`plans/ncd_exploration_stale_residual_model.md` as the plan of the NCD analysis. It follows the maintenance
rule of `CLAUDE.md`: sections describe the latest state, the Log is append-only.

## 1. Purpose and why the analysis restarts

The first NCD analysis (2026-11-14 to 2026-11-18) scored each gene in each cell by a binomial residual after
removing a per-cell shift (alpha) and a per-cell scale (kappa) of the cytoplasmic fraction. It was negative: no
reproducible substates and no stress link. The advisor's judgement (reported by the user on 2026-11-18):

- The cell-level terms, kappa in particular, may be biology and not a confounder. Removing them leaves a residual
  that is uninformative by construction. The negative result is therefore a result about that model, not about
  NCD.
- Two alternatives: (1) the plain cytoplasm-to-nucleus ratio as the NCD indicator; (2) the regression model
  without kappa.
- Alternative (1) comes first. The ratio is the most direct indicator of the cytoplasm-nucleus balance, and it is
  easier to diagnose whether NCD is feasible on a ratio than on a residual. The ratio drops the abundance of the
  gene, should still add information that absolute cytoplasmic expression does not show, and does not depend on
  cell size.

Alternative (2) is kept as a later option and is not planned here.

Question of this analysis: does the ratio score, on genes selected from the score itself, organise tumor cells
into substates that have a spatial pattern and a biological meaning? The substates need not reproduce the
substates of cytoplasmic expression (`plans/subtype_exploration.md`).

Scope: tumor cells (`cell_type_merged == "Malignant cell"`) with at least one nuclear transcript (§2) of the five
carcinoma samples (BC, OC, CC, LC, Prostate), compartment matrices of `processed_data/` (§2), sample-corrected
pooled analysis.

## 2. Score

For gene g in cell i, c = cytoplasmic count and n = nuclear count.

**Working score (decided 2026-11-18):**

    r = log2((c + 1) / (n + 1))

| Case | c | n | r | Ratio it stands for |
|---|---|---|---|---|
| No transcript | 0 | 0 | 0 | 1 (balanced, neutral) |
| Cytoplasm only | 1, 2, 3 | 0 | 1, 1.58, 2 | above 1, larger with more transcripts |
| Nucleus only | 0 | 1, 2, 3 | -1, -1.58, -2 | below 1 |
| Both | 2 | 2 | 0 | 1 |
| Both | 4 | 1 | 1.32 | 2.5 |

- One rule gives the three special cases that were to be hard coded. An empty pair gets the balanced value, a
  cytoplasm-only pair a positive value, a nucleus-only pair a negative value.
- The score is symmetric in the two compartments, and it grows with the evidence: three transcripts against none
  score higher than one against none.
- Positive r: the gene is more cytoplasmic in this cell. Negative r: more nuclear. Unit: log2 of a count ratio.

The raw ratio c / n with hard-coded values for the special cases is not used and not run.

**What c and n are** (audited against the raw data, §10.8):

- The counts come from the compartment matrices `data/<sample>/processed_data/{cytoplasmic,nuclear}_expression_matrix.npz`,
  built by `code/1_preprocessing/1_clean.ipynb` and `code/2_outcome/1_split_expression.ipynb` from the 10x
  transcript table. n = transcripts of the cell with the 10x flag `overlaps_nucleus` = 1, c = the other
  transcripts of the cell. c + n is every gene transcript assigned to the cell. (In the processed transcript
  table `1_clean.ipynb` renames the 10x flag to `in_nucleus`; the convention "nuclear: `in_nucleus == 1`,
  cytoplasmic: `overlaps_nucleus == 0`" of `CLAUDE.md` names the columns of that table and is the same split.)
- The nucleus is the 2D outline of the nuclear segmentation; depth is ignored. Cytoplasm above and below the
  nucleus counts as nuclear. The nuclear share of a cell's counts therefore follows the nucleus share of the
  cell area (correlation 0.81 to 0.89 by sample).
- No quality filter: 10 to 17 % of the transcripts are below Q20, equally in both compartments. Q20 is not
  applied (§8).

**Cells of the analysis** (decided 2026-11-18): the tumor cells with at least one nuclear transcript, summed
over all 5001 panel genes. A tumor cell with zero nuclear transcripts has n = 0 for every gene, so every
detected gene would score positive; this is an artefact of segmentation, not a state of the cell. The removed
cells are the cells without a segmented nucleus plus a few with a small nucleus that holds no transcript. The
matrices are not rebuilt: `ncd_ratio_utils.load_counts` drops these rows.

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Tumor cells | 102,180 | 160,250 | 221,355 | 44,624 | 95,429 |
| Removed: no segmented nucleus | 4686 | 4397 | 8804 | 1568 | 7857 |
| Removed: nucleus without any transcript | 14 | 1 | 72 | 2 | 48 |
| Share removed | 4.6 % | 2.7 % | 4.0 % | 3.5 % | 8.3 % |
| Cells of the analysis | 97,480 | 155,852 | 212,479 | 43,054 | 87,524 |

**Genes of the analysis** (decided 2026-11-18): the protein-coding genes of the panel, 4972 of 5001. The
biotype is the Ensembl biotype (release 116, looked up through mygene.info by `gene_biotype.py`, table
`data/_utils/gene_biotype.csv`). Dropped: 26 long non-coding RNAs (BCAR4, CASC8, CCAT2, CRNDE, CYTOR, DLEU1,
DLEU2, FENDRR, H19, HOTAIR, MEG3, MIAT, NORAD, PART1, PCA3, PCAT1, PCGEM1, PRNCR1, PVT1, SNHG1, SNHG8, SNHG11,
SNHG15, SNHG17, SOX2-OT, UCA1), the immune-receptor constant segments IGHE and TRDC, and CBSL (annotated as an
artifact). `ncd_ratio_utils.load_counts` drops these columns, so every rule below is stated among the coding
genes. The cell rule above counts the transcripts of all panel genes.

**Facts the score has to face** (step 0, §10; the cells and genes of the analysis, the 2967 genes of N = 3000):

| Quantity | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Non-empty (gene, cell) pairs | 8.2 % | 11.4 % | 5.0 % | 25.4 % | 7.7 % |
| Of the non-empty pairs: cytoplasm only | 40 % | 40 % | 51 % | 43 % | 47 % |
| Of the non-empty pairs: nucleus only | 51 % | 49 % | 42 % | 35 % | 44 % |
| Of the non-empty pairs: both compartments | 9 % | 11 % | 7 % | 22 % | 9 % |

- A nucleus-empty pair is not rare: 40 to 51 % of the non-empty pairs have transcripts only in the cytoplasm.
  The reason is counting: 78 to 82 % of the non-empty pairs hold a single transcript (59 % in LC).
- Most scores are therefore 0, +1 or -1. A score differs from 0 only where the gene is detected, so the score
  matrix carries detection, that is abundance, with a sign. Substates of the score can then repeat the
  expression substates for that reason alone. The stale document set pseudo-count log ratios aside for this
  reason (its §3.3). Here it is not a reason to drop the score; it is the main thing the checks of §5 have to
  separate.

## 3. Gene filter (steps 0 and 1)

Purpose: remove genes that are too lowly expressed in tumor cells to carry a score. The kept genes are the
candidates of the target list. The gene universe is the 4972 protein-coding panel genes (§2); the restriction
to coding genes comes first, so that "top N" means the N most expressed coding genes of a sample.

Procedure (`0_counts.py` computes the table, `1_filter.py` applies the rule), on the cells of the analysis (§2):

1. Expression of a gene in a sample = (cytoplasmic count + nuclear count, summed over the cells) / number
   of cells. The two counts together are all transcripts assigned to the cell (§2, measurement).
2. Within each sample the 4972 genes are ranked by this expression (rank 1 = highest; equal values in
   panel order).
3. A gene is kept if its rank is N or better in at least 3 of the 5 samples. Equivalently: its third-best rank
   over the samples (`entry_N`) is N or better.

| Constant (`ncd_ratio_utils.py`) | Value | Meaning |
|---|---|---|
| `N_FILTER` | 3000 | N of the rule |
| `MIN_SAMPLES` | 3 | Samples in which a gene must qualify, here and in every later rule |
| `MIN_CELLS_GE2` | 200 | Cells with two or more transcripts of the gene that a sample needs for the variation statistic of §4.2 |

Result: 2967 genes (`output/1_exploration/ncd_ratio/1_filter/filtered_genes.npy`, panel order).

Two terms used below:

- Detected: a gene is detected in a cell if c + n is at least 1.
- Measurable gene: at least 200 cells with two or more transcripts in at least 3 of 5 samples. 2436 of the 4972
  genes are measurable, and all of them are among the 2967 kept genes.

Leakage is not considered (decided 2026-11-18). All transcripts assigned to a tumor cell are treated as its
own. No gene is flagged or removed because neighbouring cell types express it. Genes of other cell types leave
the analysis only through low expression in tumor cells.

Why coding genes only, why N = 3000 and why 3 of 5 samples: §8; the evidence is in §10.4 and §10.5.

## 4. Target gene list (steps 2 and 3)

The 634 target genes of the earlier stages are not used. Three candidate lists were built from the score on
the 2967 filtered genes (step 2, `2_selection.py`): V-null and V-abs by variation of the score across cells, M
by the mean score. An audit compared them (step 3, `3_audit.py`). **The target list is V-abs; its first 500 genes are the
primary list.** M is not carried further; V-null is kept as the documented alternative (§10.7).

### 4.1 Rules common to the three lists

- **Within sample.** Every statistic is computed within one sample. Nothing is computed on the pooled cells.
- **Detected cells only.** Every statistic of a gene uses only the tumor cells in which the gene is detected.
  An empty pair scores 0 whatever the gene does; over all cells every statistic would mostly measure the
  detection rate.
- **Aggregation over samples.** In each sample the eligible genes are ranked by the statistic of the list
  (rank 1 = largest). The rank of a gene over the samples is its third-best within-sample rank (`entry_K`). The
  list of size n is the n genes with the smallest `entry_K`. This is the same as: a gene qualifies in a sample
  if it is in that sample's top K, it is selected if it qualifies in at least 3 of 5 samples, and K is set so
  that the list has n genes. Genes with the same `entry_K` are ordered by the median of the statistic over the
  samples (larger first).
- **Direction.** No list is restricted to one direction; cytoplasm-leaning and nucleus-leaning genes both
  enter. The direction of a gene is an annotation: cytoplasmic if the median over samples of its centred mean
  score (§4.2) is above 0, nuclear otherwise.

Why selection within samples and not on the pooled cells:

- Selection on the pooled raw score would pick genes whose balance differs between samples (tissue, section,
  segmentation). The per-sample standardisation of step 4 removes exactly that difference, so such genes would
  enter the embedding with nothing left.
- Selection on the pooled score after per-sample standardisation avoids that, but it is a cell-weighted average
  of the within-sample variances. CC (212,000 cells) and OC (156,000) outweigh LC (43,000), and a gene
  can enter on the strength of one sample.

### 4.2 Statistics per sample and gene

Computed by `ncd_ratio_utils.pair_stats`, written to `2_selection/gene_stats.csv`. For one gene in one sample,
over the D cells in which it is detected, with r the score of the gene in each of these cells:

| Column | Definition |
|---|---|
| `cells_detected`, `detection_rate` | D, and D / number of tumor cells |
| `cells_ge2` | Cells with two or more transcripts of the gene (c + n of at least 2) |
| `cytoplasmic_fraction` | p = sum of c / sum of (c + n) over the cells |
| `mean_r` | Mean of r |
| `mean_r_centred` | `mean_r` minus the median of `mean_r` over the 2967 genes in the sample |
| `var_r` | Variance of r (mean of r squared minus squared mean) |
| `null_var_mean`, `null_var_sd` | Mean and SD of the variance of r over 20 draws of the null (below) |
| `excess` | `var_r` / `null_var_mean` |
| `z` | (`var_r` - `null_var_mean`) / `null_var_sd` |
| `excess_t2` | Variance of r on the cells with exactly two transcripts of the gene, divided by the exact binomial variance for two transcripts with the gene's cytoplasmic fraction on these cells. Reported only with at least 200 such cells. A check, not used for any ranking |

**The count-preserving null** (`ncd_ratio_utils.null_counts`). Every (cell, gene) pair keeps its total
t = c + n. The t transcripts are reassigned to the compartments at random: c* is drawn from Binomial(t, p) with
p the cytoplasmic fraction of the gene in the sample, and n* = t - c*. The null has the same sparsity and the
same abundance as the data and lets cells differ in nothing. 20 draws per sample (`N_NULL`), random seed 0
(`SEED`).

Two properties that follow:

- Cells with one transcript of the gene carry no information on variation. Their score is +1 or -1, the
  variance of a two-valued score is fixed by its mean, and it is the same in the data and in the null.
  Differences between cells show only where a gene has two or more transcripts in a cell (18 to 22 % of the
  non-empty pairs, 41 % in LC; §10.3). This is why `cells_ge2` is required for V-null.
- A gene also exceeds the null when the whole cell leans to one compartment. This follows the decision that
  cell-level differences count as biology (§1).

### 4.3 Definition of the three lists

| | V-null | V-abs (the target list) | M |
|---|---|---|---|
| Idea | Variation above what counting alone gives | Plain variation: the most variable genes | Constant lean to one compartment |
| Statistic in a sample | `excess` | `var_r` | Absolute value of `mean_r_centred` |
| Genes ranked in a sample | Only genes with `cells_ge2` of at least 200 and `z` of at least 3 in that sample | All 2967 | All 2967 |
| Extra condition | None | None | The 3 or more samples in which the gene qualifies must agree on the sign of `mean_r_centred` |
| Aggregation | §4.1 | §4.1 | §4.1 |
| Genes that can be ranked | 2350 | 2967 | 2967 |

Constants: `Z_MIN` = 3, `MIN_CELLS_GE2` = 200, `N_NULL` = 20. The rankings are computed by
`ncd_ratio_utils.entry_ranks` and written to `2_selection/gene_summary.csv` (`entry_K_*`, `rank_*`); the first
1000 genes of each are `list_v_null.csv`, `list_v_abs.csv` and `list_m.csv`.

Notes on each list:

- **V-abs.** A cell with one transcript scores +1 or -1, so the variance of a gene is highest when its
  transcripts split half and half, whether or not cells differ. The variance therefore follows the balance and
  the expression level of the gene (§10.6). The audit weighed this against its gene content (§10.7).
- **V-null.** Qualification is on evidence (z), the rank is on the size of the excess. The 531 filtered genes
  that are not measurable cannot enter.
- **M.** The mean is centred within the sample because the overall balance differs by sample (median mean
  score BC -0.18, OC -0.13, CC +0.06, LC +0.11, Prostate -0.01). "Absolute" is measured from 0, so without
  centring the top of the BC and OC rankings would be filled from the nuclear tail and the top of the CC and LC
  rankings from the cytoplasmic tail. No minimum of detected cells is needed: every filtered gene is detected in
  at least 734 cells in at least 3 samples.

### 4.4 Size of the list

| List | Size | Constant | Use |
|---|---|---|---|
| Primary | First 500 genes of the ranking | `LIST_SIZE` | The target genes of the analysis and of the later prediction |
| Larger nested list | First 1000 | `LIST_SIZE_CHECK` | Run in step 4 as a check: do the substates change with more genes? |
| Union with V-null | V-abs first 500 plus V-null first 500: 630 genes | | Run in step 4 as a check: does the choice between V-abs and V-null matter for the substates? |
| Smaller nested list | First 250 | not set | A possible further check: can the target shrink? Not planned |

The ranked files hold the first 1000 genes (`LIST_SIZE_POOL`). Decided by the user on 2026-11-18 (first 1000
with a check at 500, then revisited to 500 with a check at 1000). Reasons in §8.

### 4.5 Audit of the lists (step 3)

Run after step 2 and before any embedding (`3_audit.py`, output `3_audit/`), at the size of the primary list. Question: which list is the most
plausible biologically? Agreed in the review of step 2: M is a reference only, so the audit decides between
V-null and V-abs; the judgement is on gene content first, not on a single number. The two lists share 370 of
500 genes, so three gene sets are compared: shared (370), only V-null (130), only V-abs (130). M and the 634
genes are reference columns.

| Aspect | Exact procedure |
|---|---|
| Overlap | Genes shared at sizes 500 and 1000; median rank of the genes of one list in the other |
| Gene content | Number of genes of each curated program (`subtype_utils.CURATED_PROGRAMS`), hallmark set (`hallmark_pathways_filtered.gmt`) and of the surface and secreted annotation in each gene set; one-sided hypergeometric test against the 2967 filtered genes, not corrected for the number of groups |
| Direction and expression | Direction, detection rate, cytoplasmic fraction per gene set |
| Reproducibility across samples | Spearman over genes of the within-sample statistic between each pair of samples (10 pairs); number of samples in which a listed gene is in the top K |
| Stability across cells | The cells of every sample are split in two random halves. The statistics of §4.2 are recomputed in each half (with 100 instead of 200 cells with two or more transcripts), the lists rebuilt per half over the five samples, and compared with each other and with the full list |
| Reliability of the per-cell score | For every pair with two or more transcripts, each transcript is assigned to half A or half B with probability 0.5 and r is computed in each half. Per gene and sample: Pearson correlation of the two half-scores over these pairs (at least 200 pairs); median over samples. Reference: the same on 5 draws of the null. The reference is slightly below 0, not exactly 0, because the two halves share the total of the pair |
| Dependence on expression | Within each sample the genes are put in 10 bins of the detection rate and the statistic is replaced by its percentile within the bin; the lists are rebuilt from these percentiles by the rule of §4.1. Reported: how many genes of the list remain |

Surface and secreted genes are annotated by gene-name patterns written in the script (`SURFACE`, `SECRETED`;
133 and 63 of the filtered genes; the annotation of every gene is in `3_audit/gene_audit.csv`). This is a
proxy: no atlas annotation is in the repo, and a few matches are not surface or secreted proteins in the strict
sense (for example CD68, SERPINH1).

### 4.6 The chosen list

V-abs (decided by the user on 2026-11-18 on the recommendation of the audit at size 1000; reasons in §8,
evidence and the closer comparison at size 500 in §10.7). Recorded in `ncd_ratio_utils.TARGET_LIST`. Files:
`3_audit/target_genes.csv` (the first 1000 genes of the ranking with their statistics; the flag
`in_primary_list` marks the first 500) and `3_audit/target_genes.npy` (the 1000 gene names in the order of the
ranking); `3_audit/target_genes_union.csv` and `.npy` (the union of the V-abs and V-null first 500, with the
set of each gene: shared, only V-abs, only V-null). `ncd_ratio_utils.load_target_genes("primary")` returns the
500, `"larger"` the 1000 and `"union"` the 630.

## 5. Substates: embedding, clustering and their checks (step 4)

Step 4 embeds and clusters the cells on the score, pooled over the five samples after a correction within each
sample. It is run for three gene lists and two versions of the correction, and an audit then picks one list
and one correction (§5.4). The form follows the sample-corrected pooled analysis of cytoplasmic expression
(`plans/subtype_exploration.md`, `1_cytoplasmic/4_pooled_embedding.py --correction sample`), so the two can be
read side by side.

### 5.1 Procedure

Script `4_embedding.py`; the correction is `ncd_ratio_utils.correct_within_sample`.

**Input.** The cells of the analysis (§2, 596,389 cells; all of them, no further cell filter) and the genes of
one list (§4.4). The score matrix of a sample is cells x genes with r of every pair; an empty pair is 0.

**Correction within each sample, per gene.** Two versions, both run:

| Step | `zscore` | `depth` |
|---|---|---|
| 1. Regress out depth: residual of the least-squares fit of r on an intercept and log10(depth + 1), over the cells of the sample | no | yes |
| 2. z-score: subtract the mean over the cells of the sample, divide by the SD (n - 1 in the denominator) | yes | yes |
| 3. Clip at -10 and +10 | yes | yes |

- Depth of a cell = all its transcripts (c + n) over the 4972 protein-coding genes, not only the genes of the
  list. The regression is the formula of the cytoplasmic analysis (`regress_out_covariate`).
- After step 2 every gene has mean 0 and SD 1 in every sample. This removes the differences between samples
  in the balance of a gene (tissue, section, segmentation) and in its spread, before pooling.
- Step 3 is the convention of the cytoplasmic analysis (there only the upper side is clipped; here both
  sides, because the score has a sign). It touches rarely detected genes: their few non-zero scores lie far
  from the mean.
- Nothing else is removed at cell level. The cell-wide lean to one compartment stays in the data (advisor,
  §1) and is reported.

**Pooled.** The five corrected matrices are stacked. Then, with the parameters of the cytoplasmic analysis:

| Step | Setting | Constant |
|---|---|---|
| PCA | 50 components, arpack solver | `N_COMPS_PCA` |
| Neighbours | 30 nearest neighbours on the first 30 PCs | `N_NEIGHBORS`, `N_PCS` |
| UMAP | scanpy defaults, on the neighbour graph | |
| Clusters | Louvain on the neighbour graph at resolutions 0.5, 0.8, 1.0, 1.2, 1.5 (Leiden if the `louvain` package is missing on HGCC; `summary.txt` says which) | `RES_GRID` |

Random seed 0 everywhere. No resolution is chosen in this step.

### 5.2 Lists, corrections and runs: 24 tasks

| Factor | Values |
|---|---|
| Gene list | `primary` (V-abs first 500), `larger` (V-abs first 1000), `union` (V-abs 500 plus V-null 500, 630 genes) |
| Correction | `zscore`, `depth` |
| Run | `data`, `null`, `half_a`, `half_b` |

| Run | Counts | Question |
|---|---|---|
| `data` | Observed | The result |
| `null` | One draw of the count-preserving null of §4.2: every pair keeps its total, the transcripts are reassigned with the gene's cytoplasmic fraction in the sample | Is the structure more than abundance and counting? |
| `half_a`, `half_b` | Every transcript of a pair goes to half A with probability 0.5, otherwise to half B; the two halves sum to the data | Do the clusters repeat on independent transcripts of the same cells? No UMAP |

Every run goes through the whole procedure of §5.1 by itself (own correction, own PCA, own clusters). The
cell covariates (depth, nuclear share) are always those of the observed cell.

### 5.3 Execution and output

Heavy part on HGCC: SLURM array `4_embedding.sh`, 24 tasks, 16 CPUs and 64 GB each. The exact commands are
in the header of `4_embedding.sh`:

1. Local: `make push`. HGCC: `git pull` in `~/hulab/projects/AIVC`.
2. Local: `make ncd-push` (sends the gene lists of steps 1 and 3).
3. HGCC, in `code/1_exploration/2_ncd_ratio/`: `sbatch 4_embedding.sh`. `python3 4_embedding.py --status`
   lists the complete tasks; `sbatch --array=<ids> 4_embedding.sh` reruns failed ones.
4. Local: `make ncd-pull` (pulls `output/1_exploration/ncd_ratio/4_embedding/`, about 5 GB).

Output per task, `4_embedding/<list>/<correction>/<run>/`:

| File | Content |
|---|---|
| `cells.csv.gz` | Per cell: sample, id, coordinates, areas, depth, nuclear share of counts, detected genes of the list, UMAP, cluster at each resolution |
| `pca.npy`, `pca_variance.csv` | The PCs (order of `cells.csv.gz`), variance explained |
| `gene_scaling.csv` | Per sample and gene: mean and SD of r, slope on depth, detection rate |
| `cluster_gene_means.csv.gz` | Per resolution, cluster and gene: mean corrected score, mean r, detection rate |
| `cluster_sizes.csv` | Clusters by sample, with the median covariates |
| `checks.csv`, `summary.txt` | Correlation of PC1 to PC10 with the cell covariates, sample effect, check of the correction |
| `connectivities.npz` | The neighbour graph (data runs) |

No count matrix is written; the cells are those of `load_counts` in the same order.

Review, local and light: `4_review.py` reads these files and writes `4_embedding/review/`
(`comparison.csv`, `by_resolution.csv`, `list_agreement.csv`, `cluster_markers.csv`, `summary.txt`,
`overview.png`, one UMAP figure per list and correction with data and null by sample, cluster, depth and
nuclear share). These are diagnostic figures. The annotated figures (tissue maps, marker heat maps, next to the
cytoplasmic substates) are step 5 and are made for the chosen list only.

### 5.4 Audit: which list, which correction

Done on the output of `4_review.py`, before step 5. Read at resolution 1.0, with the other resolutions in
`by_resolution.csv`.

| Criterion | Measure | Decides mainly |
|---|---|---|
| Sample effect | Share of a cell's neighbours from its own sample, divided by the sample's share of the cells (1 = perfect mixing); Cramer's V of cluster and sample | Correction |
| Depth and geometry | Largest Spearman of PC1 to PC10 with depth and with the nuclear share of counts; share of the variance of depth and of the nuclear share that lies between clusters | Correction |
| More than counting | Variance explained by the first 30 PCs, data against null; number of leading PCs above the null; clusters of the null | List |
| Repeatable | Adjusted Rand index of the clusters of half A and half B | List, correction |
| Agreement of the lists | Adjusted Rand index of the clusters between the three lists | Whether the choice of the list matters |
| Spatial pattern | Share of a cell's 10 nearest tumor cells in the tissue that are in its cluster, divided by the share expected from the cluster sizes; data against null | List |
| Gene content | Marker genes of a cluster = the 30 genes with the highest mean corrected score; clusters with a program or hallmark set enriched among the markers (hypergeometric against the genes of the list, p below 0.001); surface and secreted genes among the markers | List |

No single number decides. Claude gives a recommendation with the evidence; the choice is the user's.

### 5.5 Expectation and the checks of a substate

Annotation of the clusters by Claude (initial run, audit, refined run) and the figures of step 5 follow the
choice of §5.4.

Expectation (advisor): not all cytoplasmic substates are recovered. The ratio substates should show a spatial
pattern and some should have a biological meaning. One specific expectation: surface markers have an unbalanced
distribution between nucleus and cytoplasm, enter the target list, and mark a class of tumor cells at the tumor
boundary.

Checks that decide whether a cluster is a ratio substate:

| Check | Question |
|---|---|
| Count-preserving null (same pipeline on the reassigned counts of §4.2) | Is the cluster more than abundance and counting? |
| Cytoplasmic substates (agreement of labels, marker overlap) | New information, or expression seen again? Reported, not a pass criterion |
| Split halves of the transcripts | Does the cluster exist in both halves? |
| Samples | Is the cluster present in several samples? |
| Cell geometry (cell area, nucleus share of the area) per cluster | Reported, not removed |
| Space (distance to the tumor boundary, spatial autocorrelation, tissue maps) | Is there a spatial pattern? |

The surface-marker class is looked for explicitly: surface-marker genes in the target list, a cluster marked by
them, and its distance to the tumor boundary. Two reading notes. Transcripts of neighbouring cells assigned to a
boundary tumor cell would give a similar picture; no correction is applied (§3). mRNAs of surface proteins are
translated at the ER, so their imbalance is first a property of the gene; whether it differs between cells is
the question.

## 6. Fallback if single-cell scores are too sparse

Ratio of summed counts over small groups of cells (neighbours in expression space or in the tissue), which
replaces single transcripts by tens of transcripts per gene. An option only; not planned.

## 7. Planned steps

| Step | Content | Output to review |
|---|---|---|
| 0 (`0_counts.py`, done) | Counts check on the most expressed genes: share of empty, cytoplasm-only, nucleus-only and two-sided pairs per sample; pairs and cells per gene with two or more transcripts; distribution of r; detection rates by expression rank | One table and one figure; sets N |
| 1 (`1_filter.py`, done) | Gene filter (§3) | Gene list with detection rates |
| 2 (`2_selection.py`, done) | Score, null, and the three candidate lists V-null, V-abs and M (§4.1 to §4.4) | Per-gene table, the three rankings; sets the list size |
| 3 (`3_audit.py`, done) | Audit of V-null against V-abs, M as reference (§4.5); V-abs chosen (§4.6) | Comparison with recommendation; the target list |
| 4 (`4_embedding.py` on HGCC, `4_review.py` local; implemented, not yet run on the full data) | Pooled embedding and clustering for three gene lists and two corrections, each on the data, the null and two transcript halves (§5.1 to §5.3); audit of list and correction (§5.4) | Diagnostic figures and comparison table; the list and the correction to carry on |
| 5 | Figures and annotation of the clusters of the chosen list | Figures next to the cytoplasmic analysis |
| 6 | Link to stress (environmental and expression stress of `plans/stress_substate_exploration.md`) | Stress tables and figures |

Rules:

- One step at a time. A step is not implemented before the results of the previous step are reviewed with the
  user and both agree to go on. This holds strictly for steps 0 to 3, which fix the genes.
- Where a heavy step runs (local or HGCC) is asked before the run.
- Every script has a `--quick` mode where it reads the counts (20,000 random cells per sample, output in
  `<step>_quick/`). A quick run is a development check, not a result, and its folder is deleted afterwards.
- Scripts are in `code/1_exploration/2_ncd_ratio/`, numbered by step (`0_counts.py`, ...; `ncd_ratio_utils.py`
  shared by the steps; `subtype_utils.py` one level up is shared with the cytoplasmic analysis), outputs in
  `output/1_exploration/ncd_ratio/`, one numbered folder per step (map in the `README.md` there).

## 8. Decisions and their rationale

All decided on 2026-11-18 by the user unless noted; the Log has the order of events.

| Decision | Value | Why | Alternatives considered | Evidence |
|---|---|---|---|---|
| NCD indicator | Cytoplasm-to-nucleus ratio per gene and cell | Advisor: the most direct indicator of the balance; the residual model removed what may be biology | Model without kappa (later option) | §1 |
| Score | r = log2((c + 1) / (n + 1)) | One rule covers the empty, cytoplasm-only and nucleus-only pairs; symmetric; grows with the evidence | Raw ratio c / n with hard-coded special cases: dropped, not run | §2 |
| Leakage | Not considered | All transcripts of a tumor cell are treated as its own | Flagging genes of neighbouring cell types | §3 |
| Cells of the analysis | Tumor cells with at least one nuclear transcript | A cell with zero nuclear transcripts has n = 0 for every gene and scores positive on every detected gene: a segmentation artefact that would add variance to every gene and sit in the embedding as cells with only positive scores. 2.7 to 8.3 % of the tumor cells per sample | Rule by segmentation (nucleus area missing): differs by 137 cells that have a small nucleus without any transcript; the user chose the count rule, which leaves no cell with an empty nucleus | §2, §10.8 |
| Genes of the analysis | Protein-coding panel genes only (4972 of 5001), applied before the expression filter | The project studies NCD as a response to tumor stress. For a coding RNA the cytoplasmic share is what is available for translation; long non-coding RNAs are mostly held in the nucleus as a fixed property of the gene, so their balance is extreme but hardly differs between cells. Applied first, so that the universe is defined once and "top N" means the N most expressed coding genes | Keeping all genes; restricting after the expression filter (same genes within about a dozen at the margin, but the filter would have two parts to explain). The effect is small: V-abs had only NORAD and H19 among its first 1000 | §2, §10.8 |
| Quality filter | None (all decoded transcripts) | The low-quality share is the same in both compartments, so it does not shift the balance; the gain is small (variance of the score 3 to 4 % lower, gene ranking nearly unchanged), the cost is 10 to 17 % fewer transcripts in already sparse scores, and the counts would differ from the earlier stages | Q20: kept as a possible sensitivity check | §10.8 |
| Gene filter, measure | Mean count per cell (cytoplasmic + nuclear), ranked within sample | The filter is about whether tumor cells express the gene at all | | §3 |
| Gene filter, N | 3000 (2967 genes) | The filter only has to remove genes with too little data; the variance ranking then does the selecting. 3000 keeps every one of the 2436 measurable genes. The embedding uses the target list, not the filtered genes, so the sparser extra genes do not reach it. Above 3000 nothing is gained. Re-examined when the list size went to 500: a stricter N only deletes genes from the list (37 of the 500 at N = 2000, 15 at N = 2500), and the deleted genes are subset-expressed state genes (POSTN, CXCL9, SERPINE1, IFIT2, PLAU, HMMR) | 2000 (nearly every kept gene measurable, but 488 measurable genes lost) and 2500 (the value of the rule fixed before the run; 80 measurable genes lost). All three are reasonable | §10.2 to §10.4 |
| Samples required | 3 of 5, in the filter and in every list | The analysis looks for substates shared by the samples; a gene expressed or variable in one or two tumors cannot support one. 3 of 5 tolerates two samples in which a gene is weak | 5 of 5 would drop 697 of the 2967 genes; pooled selection is dominated by the large samples (§4.1) | §10.5 |
| Cells used for the statistics | Cells in which the gene is detected | Empty pairs score 0 for every gene; with them every statistic measures the detection rate | All cells | §4.1 |
| Data required for V-null | 200 cells with two or more transcripts, in the sample | Single-transcript cells carry no information on variation, so the statistic needs cells with two or more transcripts. 200 is a round minimum fixed before step 0; it was not tuned. Also the definition of a measurable gene | | §4.2, §10.1 |
| Null | Count-preserving, 20 draws, evidence z of at least 3 | Same sparsity and abundance as the data. 20 draws are enough to estimate the mean and SD of the null variance; z of 3 is a conventional evidence threshold, fixed before the run. The threshold matters little: 95.0 to 99.5 % of the genes with enough cells pass it, by sample | | §10.6 |
| Candidate lists | V-null, V-abs and M, same size, same aggregation | The user asked whether absolute variation or variation above a null is the better rule and to try both; M was the alternative by mean score | | §4.3 |
| List size | Primary list 500; the first 1000 as nested check; 250 a possible further check | The selected genes are the later prediction target (their NCD, or their cytoplasmic expression if NCD is not feasible), and more genes, above all weakly informative ones, make the prediction harder. The gene statistics give no natural cut: quality falls gradually along the ranking. The first 500 are better detected (median 15 % against 9 % for the second 500), the samples agree on them more, they hold 63 % of the non-zero scores of the 1000 and 37 of its 55 surface genes. The lists are nested, so step 4 can test whether 1000 changes the substates | 1000 as primary (the first decision, when sparsity of the score was the main argument); 250, 750 | §10.7 (rank bands) |
| M | Reference only, not a candidate | M ranks genes that sit constantly in one compartment, not genes that differ between cells, and substates need the latter. 412 of its first 500 genes are nuclear-leaning, a third of the first 500 are not measurable, and it holds about half as many program genes as the variation lists | | §10.6 |
| Target list | V-abs | Decided on the audit at size 1000, on biological plausibility of the gene content: there the genes only in V-abs were surface receptors, junction and signalling genes and those only in V-null mostly nuclear-retained RNAs, and V-abs loses little by the null criterion (nearly all of its own genes exceed the null). Both lists are equally stable under a split of the cells | V-null: less tied to the expression level. At size 500 the comparison is closer (370 of 500 genes shared; V-abs ahead on surface and secreted genes, V-null on interferon, proliferation and TNFa genes and on independence from expression: 363 against 275 genes kept when expression is held fixed). The user kept V-abs and added the union of the two lists (630 genes) as a check of step 4, to see whether the choice matters for the substates | §10.7 |
| Correction before pooling | Within each sample, per gene: z-score, in two versions, without and with a regression on depth first; both run | Per-sample z-scoring removes the sample differences in the balance and spread of a gene, as in the cytoplasmic analysis. Whether depth should also be regressed out is not clear for a ratio: deeper cells have more non-zero scores, so depth can drive the PCs; on the other hand the regression may be unnecessary. The user chose to run both and let the audit decide, for example on the sample effect | Only one version | §5.1, §5.4 |
| Cell-level lean | Not removed, reported | Advisor: the cell-level terms may be biology (§1). Depth is the amount of RNA, not the balance, so the depth version does not remove the lean | Regressing out the nuclear share of counts | §5.1 |
| Embedding parameters | 50 PCs computed, 30 used, 30 neighbours, Louvain at five resolutions | The values of the cytoplasmic analysis, so the two analyses can be read side by side. Not tuned | | §5.1 |
| Runs per list | Data, count-preserving null, two transcript halves | The null shows what abundance and counting alone give; the halves show whether a cluster repeats on independent transcripts. All heavy work of the audit is done in one HGCC round | Data only; data and null | §5.2 |
| Where step 4 runs | HGCC, SLURM array of 24 tasks; the user submits and pulls | The user's choice | Local | §5.3 |
| Choice of list and correction | Audit right after step 4, before the annotation of step 5 | The annotation effort then goes into one list. The diagnostic figures of step 4 are enough for the choice | After step 5 | §5.4 |
| Order of work | One step at a time, each after review | The gene selection fixes everything downstream | | §7 |

Open:

- The model without kappa (advisor's alternative 2), after the ratio analysis.
- 250 genes as a further nested check of the list size.
- Q20 as a sensitivity check, if a result looks fragile.

## 9. Where the earlier analysis is

| Item | Location |
|---|---|
| Plan document of the residual model (stale) | `plans/ncd_exploration_stale_residual_model.md` |
| Scripts `7` to `13` and `ncd_utils.py` | `code/1_exploration/old/ncd_residual/` (do not run from there) |
| Outputs and run logs | `output/1_exploration/old/ncd_residual/` |

## 10. Results

All results of §10.1 to §10.7 are from the last rerun of 2026-11-18, on the cells of the analysis (tumor cells
with at least one nuclear transcript) and the genes of the analysis (protein-coding panel genes; §2). What the
two reruns changed is in §10.8.

### 10.1 Step 0: counts check (run 2026-11-18)

What was run:

- Script `code/1_exploration/2_ncd_ratio/0_counts.py` with `ncd_ratio_utils.py`, local, about 40 s. Log:
  `code/1_exploration/2_ncd_ratio/logs/0_counts.log`.
- Input: the cytoplasmic and nuclear count matrices of the five samples, the 4972 protein-coding genes. Cells: BC 97,480,
  OC 155,852, CC 212,479, LC 43,054, Prostate 87,524.
- Output: `output/1_exploration/ncd_ratio/0_counts/`. Start with `summary.txt` and `counts_overview.png`.

| File | Content |
|---|---|
| `gene_table.csv` | Per sample and gene: expression, rank, detection, cells with two or more transcripts, pair types, mean and variance of r |
| `gene_consensus.csv` | Per gene: the N at which it enters the set, measurable or not, medians over samples |
| `by_n.csv` | Per sample and N: sparsity, pair types, the genes at the margin |
| `pair_composition.csv` | Non-empty pairs by number of transcripts and type |
| `retention.csv` | Known gene groups kept by N |
| `measurable_curve.csv` | Share of measurable genes along the order of entry |

Definitions:

- Gene set of N: genes in the top N by mean count per cell in at least 3 of 5 samples (§3). The set is
  slightly smaller than N (2000: 1950 genes; 2500: 2487; 3000: 2967).
- Measurable gene: at least 200 cells with two or more transcripts in at least 3 of 5 samples. This is the data
  the variation statistic needs (§4.2). 2436 of the 4972 genes are measurable.

Check: the summed counts of the pairs equal the sums of the matrices in every sample. With all 4972 genes the
share of non-empty pairs is BC 5.3 %, OC 7.2 %, CC 3.2 %, LC 16.1 %, Prostate 5.0 %.

### 10.2 Sparsity by N

Range over BC, OC, CC and Prostate; LC separately, because its depth is about three times higher.

| N | Genes | Non-empty pairs | LC | Pairs with two or more transcripts | LC | Non-zero scores per cell (median) | LC |
|---|---|---|---|---|---|---|---|
| 1500 | 1395 | 8.2 to 18.3 % | 39 % | 1.6 to 4.7 % | 18.1 % | 98 to 225 | 500 |
| 2000 | 1950 | 6.7 to 15.2 % | 33 % | 1.3 to 3.6 % | 14.4 % | 112 to 261 | 596 |
| 2500 | 2487 | 5.7 to 12.9 % | 29 % | 1.0 to 2.9 % | 12.0 % | 122 to 283 | 659 |
| 3000 | 2967 | 5.0 to 11.4 % | 25 % | 0.9 to 2.5 % | 10.4 % | 128 to 297 | 691 |
| 4972 | 4972 | 3.2 to 7.2 % | 16 % | 0.6 to 1.6 % | 6.4 % | 135 to 315 | 733 |

The share of non-empty pairs falls slowly with N. From 2000 to 3000 the set grows by about 1000 genes, the
share of non-empty pairs falls by 2 to 4 points (8 in LC), and a cell gains 16 to 36 non-zero scores (95 in LC).

### 10.3 Types of non-empty pairs (N = 3000)

| Transcripts in the pair | Share of non-empty pairs | Cytoplasm only | Nucleus only | Both |
|---|---|---|---|---|
| 1 | 78 to 82 % (LC 59 %) | 44 to 55 % | 45 to 56 % | 0 |
| 2 | 14 to 15 % (LC 21 %) | 27 to 39 % | 27 to 37 % | 31 to 38 % |
| 3 to 4 | 4 to 5 % (LC 13 %) | 15 to 26 % | 14 to 22 % | 56 to 66 % |
| 5 and more | 1 to 2 % (LC 7 %) | 5 to 12 % | 3 to 7 % | 84 to 91 % |

The cytoplasm, nucleus and both columns are shares within the row, ranges over the five samples.

- The advisor expected cytoplasm-only pairs to be rare. They are 40 to 51 % of the non-empty pairs. The cause is
  counting: most pairs hold one transcript, which is in one compartment by necessity. Where a gene has five or
  more transcripts in a cell, 5 to 12 % of the pairs are still cytoplasm-only and 3 to 7 % nucleus-only.
- The balance of cytoplasm-only against nucleus-only differs by sample: nucleus-only is ahead in BC and OC,
  cytoplasm-only in CC, LC and Prostate. The per-sample standardisation has to absorb this sample effect.
- For the variation statistic only the pairs with two or more transcripts count (§4.2): 18 to 22 % of the
  non-empty pairs, 41 % in LC.

### 10.4 Choice of N

| N | Genes | Measurable in the set (of 2436) | Not measurable in the set | Measurable genes left out | Measurable among the last 250 genes entered |
|---|---|---|---|---|---|
| 1500 | 1395 | 1395 | 0 | 1041 | 100 % |
| 2000 | 1950 | 1948 | 2 | 488 | 99 % |
| 2500 | 2487 | 2356 | 131 | 80 | 66 % |
| 3000 | 2967 | 2436 | 531 | 0 | 7 % |
| 3500 | 3469 | 2436 | 1033 | 0 | 0 % |

Known gene groups kept (curated programs and hallmark sets of the earlier stages; 1680 coding panel genes in these
groups, 973 of them measurable):

| | N = 2000 | N = 2500 | N = 3000 |
|---|---|---|---|
| Measurable program and hallmark genes kept (of 973) | 766 (79 %) | 936 (96 %) | 973 (100 %) |
| The 634 target genes kept (633 are coding) | 509 (80 %) | 574 (91 %) | 613 (97 %) |
| Interferon (curated), measurable genes kept | 72 % | 92 % | 100 % |
| Complement / acute phase (curated) | 50 % | 90 % | 100 % |
| Hallmark hypoxia | 75 % | 99 % | 100 % |
| Hallmark TNFA signalling via NFKB | 65 % | 93 % | 100 % |
| Hallmark EMT | 58 % | 87 % | 100 % |

**Conclusion: 2000, 2500 and 3000 are all reasonable choices. The analysis proceeds with N = 3000 for now.**

| N | For | Against |
|---|---|---|
| 2000 | Nearly every gene of the set is measurable; the densest scores | 488 measurable genes (20 %) are not candidates, among them a quarter or more of the interferon, hypoxia and TNFA genes |
| 2500 | The value of the rule fixed before the run (largest N at which at least half of the last 250 genes entered are measurable: 2598, rounded down; it keeps 90 % of the measurable program genes from N = 2304 on). Balance of the two concerns | 80 measurable genes are left out |
| 3000 | All measurable genes are candidates | 531 genes of the set are not measurable; scores slightly sparser over the whole set |

Rationale for 3000 (user's decision, 2026-11-18; Claude agrees):

- The purpose of the filter is to define the candidates for the target list. More measurable candidates is the
  gain that matters, and 3000 keeps nearly all of them.
- The cost is small. The substate analysis embeds the target list (500 genes), not the 2967 filtered
  genes, so the zeros of the extra genes do not reach it. The 531 non-measurable genes cannot enter list
  V-null. They can enter lists V-abs and M.
- The rule fixed before the run counted the non-measurable genes as a cost. Given the point above that cost is
  low, so the rule was stricter than needed.
- Above 3000 nothing is gained: no measurable gene has an entry rank above 3000.

What the filter drops at N = 3000:

- Genes of other cell types, as intended: Myeloid border (curated) keeps 9 of 24 genes, hallmark Allograft
  rejection 41 %, Inflammatory response 50 %, EMT 56 % of their panel genes.
- Surface genes checked by hand: EPCAM, F11R, CD47, CDH1, CD44, ERBB3, ERBB2, EGFR, MUC1, ITGAV, CDH3, ITGB4,
  ITGA6 and MET are in the set. CD274, PROM1, FOLR1 and MSLN are out by the 3-of-5 rule: each is in the top
  3000 of one or two samples only (FOLR1 detected in 77 % of the LC and 10 % of the OC cells, PROM1 in 21 % of
  LC, MSLN in 11 % of OC, CD274 in at most 2.5 %) and in under 1 % of the cells of the other samples.

**N and the target list** (asked again on 2026-11-18, when the list size went to 500 and the coding
restriction was added). Tightening N only deletes genes from the ranking; the order of the remaining genes does
not change. So the effect on the list is exact:

| | N = 2000 | N = 2500 | N = 3000 (used) |
|---|---|---|---|
| Genes kept by the filter | 1950 | 2487 | 2967 |
| Primary list of 500: genes shared with the list at N = 3000 | 463 | 485 | 500 |
| List of 1000: genes shared | 889 | 961 | 1000 |

- The threshold hardly matters for the list of 500: 463 of its genes are already in the top 2000 by expression.
- The genes a stricter N would delete are low on average because only some cells express them, and they vary
  strongly where they are expressed. At N = 2000: FCGBP, POSTN, CXCL9, CP, SFRP4, SERPINE1, COL5A1, IFIT2,
  UBD, ASPM, HMMR, PLAU, TGFB2, PLAT, PLK4, IL18 among 37 genes. At N = 2500: 15 genes, among them SFRP4,
  HMMR, NES, LEF1, INHBA, IL18.
- Removing the 29 non-coding genes from the 5001 does not loosen the threshold: "top 3000" moves by about a
  dozen ranks.
- The user kept N = 3000.

Consequences for step 2:

- The list by variation above the null is selected among the measurable genes; the other lists among all 2967
  genes.
- The comparison with the null has to be made within sample (the pair types differ by sample, §10.3).

### 10.5 Step 1: gene filter (run 2026-11-18)

Script `code/1_exploration/2_ncd_ratio/1_filter.py` (a few seconds; reads the per-gene table of step 0). Output
`output/1_exploration/ncd_ratio/1_filter/`: `filtered_genes.npy` (the list, input of step 2; read by
`ncd_ratio_utils.load_filtered_genes`), `gene_filter.csv` (all 4972 genes with ranks and counts per sample),
`dropped_sample_specific.csv`, `kept_per_sample.csv`, `top_n_overlap.csv`, `groups_kept.csv`, `summary.txt`,
`filter_overview.png`.

Result: **2967 of the 4972 protein-coding genes are kept.**

| Genes | Number |
|---|---|
| Kept, in the top 3000 of all 5 samples | 2270 |
| Kept, in the top 3000 of 4 samples | 386 |
| Kept, in the top 3000 of 3 samples | 311 |
| Dropped, in the top 3000 of 2 samples | 325 |
| Dropped, in the top 3000 of 1 sample | 523 |
| Dropped, in the top 3000 of no sample | 1157 |

Of the kept genes 2436 are measurable (the data V-null needs) and 531 are not; no measurable gene is dropped.

The kept genes per sample:

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Genes of the sample's own top 3000 that are kept | 2762 | 2797 | 2772 | 2733 | 2763 |
| Detection rate, median over genes | 5.4 % | 7.6 % | 2.9 % | 20.4 % | 4.9 % |
| Detection rate, 10th to 90th percentile | 1.2 to 18 % | 1.8 to 26 % | 0.7 to 11 % | 3.8 to 55 % | 1.0 to 18 % |
| Cells with two or more transcripts, median over genes | 563 | 1160 | 578 | 1933 | 423 |
| Genes with fewer than 200 such cells in this sample | 747 | 416 | 712 | 355 | 979 |
| Cytoplasmic fraction, median over genes | 0.42 | 0.44 | 0.53 | 0.55 | 0.49 |
| Genes with a positive mean r on detected cells | 25 % | 20 % | 61 % | 65 % | 47 % |

- The samples agree well on what is expressed: any two per-sample top-3000 lists share 2582 to 2665 genes, and
  91 to 93 % of each sample's own top 3000 are kept.
- Direction: by the median over samples, 1274 kept genes are more cytoplasmic (mean r above 0) and 1693 more
  nuclear, so both directions are well represented among the candidates. The sign agrees in all five samples
  for only 1095 genes, because the overall balance differs by sample (BC and OC lean nuclear; CC and LC lean
  cytoplasmic; Prostate is balanced). This is the sample effect of §10.3 seen per gene; mean scores must be
  compared within sample (relevant for list M).
- The cytoplasmic fraction of a gene ranges from 0.10 to 0.17 up to 0.56 to 0.71 (1st and 99th percentile over
  genes, range over the samples),
  with a sharp upper edge in every sample. No gene is close to fully cytoplasmic, while some are close to fully
  nuclear. The likely reason is the measurement (not tested directly): the nucleus is a 2D outline, so
  cytoplasmic transcripts above and below the nucleus count as nuclear (§2, §10.8).

What is dropped:

- Genes of other cell types, as in §10.4 (Myeloid border 15 of 24 dropped; hallmark Allograft rejection,
  Inflammatory response, Kras signalling and EMT lose 44 to 62 % of their panel genes).
- Sample-specific genes. 848 dropped genes are in the top 3000 of one or two samples; 119 of them have a rank
  of 1000 or better there and 63 a rank of 500 or better (Prostate 19, CC 14, LC 11, BC 10, OC 9). These are
  lineage and tissue genes of the tumor: NKX3-1, TRPM8, FOLH1, HOXB13, KLK4 (prostate); GATA3, PDZK1,
  ANKRD30A (breast); NKX2-1 (lung); DSG3, SERPINB4, A2ML1 (squamous, CC); MUC16, SOX17, FOLR1 (ovary). The full
  list is `dropped_sample_specific.csv`.
- These genes are dropped by design: the analysis looks for substates shared by the samples, and a gene
  expressed in one or two tumors cannot support a shared substate. If they are wanted later (for example for a
  per-sample analysis), the rule can be relaxed; the table gives what is needed.
- Curated programs are nearly complete (Interferon 25 of 27, Complement 11 of 12, Stress response 16 of 17,
  Lysosome and Cholesterol synthesis complete); cell-cycle and growth hallmark sets are complete or nearly so
  (E2F targets 111 of 111, G2M checkpoint 114 of 116, MYC targets 79 of 79). 613 of the 634 target genes are
  kept.

### 10.6 Step 2: selection statistics and candidate lists (run 2026-11-18)

Script `code/1_exploration/2_ncd_ratio/2_selection.py` (local, 4 min; 20 null draws per sample). Output
`output/1_exploration/ncd_ratio/2_selection/`: `gene_stats.csv` (per sample and gene), `gene_summary.csv` (per
gene: medians, rank under the three rules, direction), `list_v_null.csv`, `list_v_abs.csv`, `list_m.csv` (first
1000 genes of each ranking), `list_overlap.csv`, `groups_in_lists.csv`, `summary.txt`,
`selection_overview.png`. Check: mean and variance of r per gene reproduce the table of step 0.

**Almost every gene with enough data varies more than counting, but by little.**

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Genes with at least 200 cells holding two or more transcripts | 2220 | 2551 | 2255 | 2612 | 1988 |
| Of these, above the null (z of at least 3) | 2110 | 2431 | 2243 | 2591 | 1951 |
| Excess (observed / null variance), median | 1.09 | 1.05 | 1.07 | 1.11 | 1.07 |
| Excess, 90th percentile | 1.15 | 1.09 | 1.12 | 1.27 | 1.11 |
| Excess on pairs with exactly two transcripts, median | 1.46 | 1.30 | 1.49 | 1.26 | 1.37 |

- 2350 genes are above the null in at least 3 of 5 samples (1492 in all five). The null does not thin the
  candidates; it only orders them.
- The excess is small: the variance of r is typically 5 to 11 % above counting. Only 19 genes have a median
  excess above 1.2 and 491 above 1.1. EEF1G stands alone at 2.1.
- On pairs with exactly two transcripts the variance is 26 to 49 % above the binomial expectation: the two
  transcripts of a gene in a cell sit in the same compartment more often than chance. Over all detected cells
  this is diluted by the single-transcript cells, which is why the overall excess is small.

**What each ranking follows.**

| Ranking | Spearman of the gene's rank with its detection rate | Reading |
|---|---|---|
| V-null | -0.64 | Favours highly expressed genes: they have more cells with several transcripts, where the excess shows |
| V-abs | -0.56 | Favours highly expressed and balanced genes (variance against detection rate: 0.42 to 0.77 by sample) |
| M | -0.03 | Independent of expression level |

The excess at exactly two transcripts runs the other way (Spearman with the detection rate -0.69 to -0.83):
at equal count, rarely detected genes are the more over-dispersed ones. Its correlation with the overall excess
is 0.22. So the V-null ranking is partly a ranking by expression level.

**The three lists.**

| Size | List | Top K needed | Cytoplasmic / nuclear | Median detection rate | Measurable genes |
|---|---|---|---|---|---|
| 250 | V-null | 337 | 166 / 84 | 14.9 % | 100 % |
| 250 | V-abs | 350 | 213 / 37 | 20.8 % | 100 % |
| 250 | M | 324 | 6 / 244 | 3.4 % | 60 % |
| 500 | V-null | 597 | 338 / 162 | 12.7 % | 100 % |
| 500 | V-abs | 622 | 413 / 87 | 14.8 % | 99 % |
| 500 | M | 632 | 88 / 412 | 4.3 % | 66 % |
| 1000 | V-null | 1066 | 641 / 359 | 10.2 % | 100 % |
| 1000 | V-abs | 1075 | 768 / 232 | 11.3 % | 98 % |
| 1000 | M | 1161 | 377 / 623 | 5.7 % | 76 % |

"Top K needed" is the within-sample K that gives the size. It is close to the size itself, so the samples
agree well on all three rankings.

Along the V-null ranking the median excess of the gene at the rank falls gradually (single genes, so the
values are noisy): 1.14 to 1.16 at ranks 50 and 100, 1.10 at 250, 500 and 750, 1.07 at 1000 and 1500, 1.04 at
2000. There is no clear elbow.

Genes shared (size 1000; size 500 in brackets):

| | V-null | V-abs | M | 634 target genes (613 filtered) |
|---|---|---|---|---|
| V-null | | 770 (370) | 368 (109) | 364 (197) |
| V-abs | | | 240 (65) | 387 (217) |
| M | | | | 210 (49) |

- V-null and V-abs are largely the same list (770 of 1000 and 370 of 500 shared). M is a different list.
- V-null, first genes: EEF1G, POSTN, NOTCH2NLA, NR4A1, YWHAZ, LDHA, H3F3B, VPS13A, ZC3HAV1L, NDRG1, SOST,
  EIF3C, TUBB, CXCL10, PABPC1L, CDH1, GRN, CSNK1A1, TRIM25, PKM. A mix of abundant cytoplasmic genes,
  nuclear-retained transcripts and known state genes (NR4A1, NDRG1, CXCL10).
- V-abs, first genes: EEF1G, H3F3B, YWHAZ, EIF3C, PKM, CSNK1A1, NDRG1, LDHA, RAC1, TUBB, NR4A1, HNRNPD, EIF4G1,
  ERRFI1, CDH1, HSPA8, EPCAM. Mostly abundant, balanced, cytoplasm-leaning genes.
- M, first genes: PABPC1L, LENG8, NOTCH2NLA, U2AF1, SOST, TOMM40, YJEFN3, CYTH2, KCNIP2, RSPH9. The top of M is
  the nuclear tail (412 of the first 500 genes are nuclear), because the nuclear tail is longer: the
  cytoplasmic fraction of a gene goes down to 0.1 but not above about 0.7 (§10.5). Many of these genes are
  rarely detected and have little excess over the null; a third of the first 500 are not measurable.
- Known programs (genes in the first 1000 of V-null / V-abs / M): Interferon 16 / 15 / 9 of 25, Stress response
  10 / 12 / 7 of 16, hallmark Hypoxia 36 / 38 / 16 of 75, E2F targets 57 / 56 / 29 of 111, TNFA signalling
  41 / 41 / 18 of 79. The variation lists hold about twice as many program genes as M.
- Surface genes: CDH1, EPCAM and CD47 are in the first 100 of both variation lists; ERBB3 in the first 250 of
  both; CD44 and F11R in the first 1000 of both, ITGB4 in the first 1000 of V-abs. MUC1, EGFR (V-null), ERBB2
  and ITGAV rank lower. In M only EPCAM and CDH1 are in the first 500.

Reviewed 2026-11-18: M as a reference only; list size 500 with 1000 as the nested check (§4.4); audit in §10.7.

### 10.7 Step 3: audit of V-null against V-abs (run 2026-11-18)

Script `code/1_exploration/2_ncd_ratio/3_audit.py` (local, 4.5 min; procedure in §4.5). The audit is read at
the size of the primary list, 500; the numbers at 1000 are in `summary.txt` and `set_profile.csv`.
Settings: random seed 0; cell split with 100 cells with two or more transcripts required per half; reliability
on pairs with two or more transcripts, at least 200 per gene and sample, 5 null draws; 10 detection-rate bins.
Output `output/1_exploration/ncd_ratio/3_audit/`:

| File | Content |
|---|---|
| `summary.txt`, `audit_overview.png` | The comparison in text and in one figure |
| `set_profile.csv` | One row per gene set (shared, only V-null, only V-abs, the lists, the references) |
| `overlap.csv`, `groups.csv` | Overlaps; programs, surface and secreted genes with hypergeometric p |
| `reproducibility.csv`, `stability.csv`, `expression_dependence.csv` | The robustness tables |
| `gene_audit.csv` | Per gene: set, annotation, reliability, found in both halves, kept under expression matching |
| `half_stats.csv` | Statistics of the cell split per sample, half and gene |
| `genes_shared.csv`, `genes_only_v_null.csv`, `genes_only_v_abs.csv` | The three gene sets, ranked |
| `target_genes.csv`, `target_genes.npy` | The V-abs ranking, first 1000 genes; `in_primary_list` marks the first 500 (§4.6) |
| `target_genes_union.csv`, `target_genes_union.npy` | Union of the V-abs and V-null first 500 (630 genes; column `set`: shared, only V-abs, only V-null), the second check of step 4 |

**The gene sets at 500.** 370 genes are shared; each list has 130 genes of its own. (At 1000: 770 shared, 230
each.)

| | Shared | Only V-null | Only V-abs | M (reference) |
|---|---|---|---|---|
| Genes | 370 | 130 | 130 | 500 |
| Cytoplasmic / nuclear lean | 304 / 66 | 34 / 96 | 109 / 21 | 88 / 412 |
| Median detection rate | 15.0 % | 7.1 % | 14.6 % | 4.3 % |
| Median cytoplasmic fraction | 0.54 | 0.40 | 0.52 | 0.31 |
| Median excess over the null | 1.118 | 1.105 | 1.091 | 1.067 |
| Above the null in all 5 samples | 86 % | 64 % | 76 % | 33 % |
| In the 634 target genes | 173 | 24 | 44 | 49 |
| Surface genes (of 133) | 29 | 4 | 8 | 14 |
| Secreted genes (of 63) | 16 | 4 | 7 | 15 |

- The genes only in V-abs are close to V-null: 126 of the 130 exceed the null in 3 samples, and their median
  rank in V-null is 661. The genes only in V-null are further from V-abs (median rank 989).
- Only V-null is a nuclear-leaning set. Its first genes are rarely detected, nuclear-retained transcripts
  (ZC3HAV1L, SOST, PABPC1L, KAT2A, OGT, KCNC3, LENG8). It also holds state genes: proliferation (E2F1, MYBL2,
  TOP2A, EZH2, SKP2, CHAF1A), inflammation and stress (CXCL2, LIF, PLAUR, GADD45A, IFIT1).
- Only V-abs is a cytoplasm-leaning, well-detected set. It holds surface and signalling genes (F11R, BCAM,
  ALCAM, DDR1, ITGA2, NECTIN2, TGFBR1, IL6ST, CTNNB1, SOX9, ESR1, AHR, TGFB2, PDGFA), next to abundant genes of
  general cell function (GANAB, NONO, HNRNPF, CANX, YWHAQ, KPNB1, PSMC5).

**Gene content of the whole lists at 500.** Mixed: V-abs is ahead on surface and secreted genes, V-null on
several programs.

| Group (genes among the filtered) | V-null | V-abs | M |
|---|---|---|---|
| Surface (133) | 33 (p = 0.011) | 37 (p = 0.0008) | 14 |
| Secreted (63) | 20 (p = 0.002) | 23 (p = 0.0001) | 15 |
| Hallmark apical junction (67) | 17 (p = 0.048) | 16 (p = 0.086) | 7 |
| Hallmark epithelial mesenchymal transition (56) | 22 (p = 0.00005) | 20 (p = 0.0005) | 8 |
| Hallmark hypoxia (75) | 23 (p = 0.002) | 23 (p = 0.002) | 8 |
| Hallmark TNFa signalling via NFkB (79) | 27 (p = 0.0001) | 23 (p = 0.004) | 7 |
| Hallmark E2F targets (111) | 33 (p = 0.0004) | 31 (p = 0.002) | 12 |
| Hallmark G2M checkpoint (114) | 34 (p = 0.0003) | 30 (p = 0.006) | 12 |
| Curated interferon (25) | 12 (p = 0.0003) | 9 (p = 0.016) | 2 |
| Curated lysosome (22) | 5 (p = 0.31) | 7 (p = 0.063) | 5 |

p: hypergeometric, against the 2967 filtered genes, not corrected for the number of groups.

**Robustness at 500.**

| | V-null | V-abs | M |
|---|---|---|---|
| Spearman of the statistic between samples (mean of 10 pairs) | 0.60 | 0.56 | 0.49 |
| Spearman of the statistic between cell halves (mean of 5 samples) | 0.92 | 0.96 | 0.96 |
| List found again in both cell halves (of 500) | 421 | 437 | 460 |
| Spearman of the statistic with the detection rate | 0.36 | 0.59 | -0.02 |
| List kept when expression is held fixed (of 500) | 363 | 275 | 399 |
| Median reliability of the per-cell score | 0.222 | 0.211 | 0.147 |

- Both lists are stable under a split of the cells. Of the genes of their own, each list finds about 70 %
  again in both halves.
- V-abs depends more on the expression level, and more so at 500 than at 1000: 275 of its 500 genes (55 %)
  stay when the statistic is ranked within detection-rate deciles, against 363 (73 %) for V-null. Of the 130
  genes only in V-abs, 33 % stay. After this matching the two lists share 312 genes.
- Reliability does not separate the lists. The correlation of the two half-scores is about 0.2 for every gene
  set, also for all filtered genes (0.202), against -0.01 under the null. The per-cell score carries a signal
  that is shared by nearly all genes. It is a cell-level lean to one compartment (what the archived model
  called the per-cell scale), and it largely follows cell geometry (§10.8). Step 4 has to show whether there
  is structure beyond it.

**The choice at size 500.** V-abs was chosen by the user on 2026-11-18 on the audit at
size 1000, where the genes only in V-abs were the more plausible set and V-abs was ahead on the surface,
junction and EMT groups. At size 500 the comparison is closer:

- For V-abs: more surface and secreted genes (37 and 23 against 33 and 20); its own genes are well detected
  and nearly all exceed the null; slightly more stable under a split of the cells.
- For V-null: more interferon, proliferation and TNFa genes; much less tied to the expression level (73 %
  against 55 % kept under matching); its own genes include known state genes.
- Against V-abs: part of its own genes are abundant genes of general cell function.
- Against V-null: its own genes are led by rarely detected, nuclear-retained transcripts.

The two lists share 370 of 500 genes.

**Decision (user, 2026-11-18): V-abs first 500 stays the primary list.** Step 4 has two checks: the V-abs
first 1000 (nested), and the union of the V-abs and V-null first 500. The union has 630 genes: the 370 shared,
the 130 only in V-abs and the 130 only in V-null (447 cytoplasm-leaning and 183 nucleus-leaning genes, median
detection rate 13.2 %, 241 of the 634 genes, 41 surface and 27 secreted genes). If the substates of the union match those of
the primary list, the choice between the two rules does not matter for the result.

The primary list (V-abs, first 500): 413 cytoplasm-leaning and 87 nucleus-leaning genes, median detection rate
14.8 %, 217 of the 634 target genes of the cytoplasmic analysis, 37 surface and 23 secreted genes, 84 % above
the null in all five samples. The larger list (first 1000): 768 / 232, median detection rate 11.3 %, 387 of
the 634, 55 surface and 31 secreted genes.

**The V-abs ranking in bands** (basis of the list size, §8):

| Rank band | Median variance of r | Median detection rate | In known programs | In the 634 genes | In the top 1000 of all 5 samples | Median reliability |
|---|---|---|---|---|---|---|
| 1 to 100 | 1.20 | 29 % | 57 % | 42 % | 77 % | 0.20 |
| 101 to 250 | 1.16 | 17 % | 61 % | 42 % | 59 % | 0.21 |
| 251 to 500 | 1.13 | 11 % | 46 % | 45 % | 30 % | 0.22 |
| 501 to 750 | 1.11 | 9.5 % | 48 % | 36 % | 14 % | 0.22 |
| 751 to 1000 | 1.10 | 8.1 % | 38 % | 32 % | 1 % | 0.22 |
| 1001 to 1250 | 1.09 | 7.7 % | 39 % | 30 % | 0 % | 0.22 |

The quality of the genes falls gradually with the rank and there is no break. The second 500 are less often
detected, the samples agree on them less (the last-but-one column falls toward rank 1000 partly by
construction, since the rank is the third-best sample), and they hold fewer program genes. The per-cell
reliability does not fall with the rank.

### 10.8 Audit of the compartment count matrices (run 2026-11-18)

All steps read the compartment matrices, not the raw data. Script
`code/1_exploration/2_ncd_ratio/check_upstream.py` (local, 4 min; not a step of the analysis, it changes no
result). Output `output/1_exploration/ncd_ratio/upstream_check/`: `summary.txt`, `upstream_overview.png`,
`rebuild_check.csv`, `qv_by_compartment.csv`, `qv_by_gene.csv`, `variant_summary.csv`, `variant_by_sample.csv`,
`variant_genes.csv`, `geometry.csv`, `geometry_by_gene.csv`.

**A. The matrices are built correctly.** The code of the two notebooks was read and the matrices were rebuilt
independently from `raw_data/transcripts.parquet` (nuclear = 10x flag `overlaps_nucleus` = 1, cytoplasmic = 0,
for the cell ids and genes of the matrices).

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Transcripts of the tumor cells (million) | 35.6 | 81.8 | 43.5 | 71.4 | 32.1 |
| Entries that differ from the stored nuclear matrix | 0 | 0 | 0 | 0 | 0 |
| Entries that differ from the stored cytoplasmic matrix | 0 | 0 | 0 | 0 | 0 |
| With Q20: entries of c + n that differ from the 10x count matrix | 0 | 0 | 0 | 0 | 0 |

Rows are the malignant cells of the current `adata.h5ad` in the same order, ids unique; columns are the panel
genes in the order of the AnnData. The last row also confirms that cells and counts are aligned in the AnnData.

**B. Three properties of the measurement** (not errors of the code):

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Tumor cells without a segmented nucleus | 4686 (4.6 %) | 4397 (2.7 %) | 8804 (4.0 %) | 1568 (3.5 %) | 7857 (8.2 %) |
| Of these, with any nuclear transcript | 0 | 0 | 0 | 0 | 0 |
| Cells with a nucleus but zero nuclear transcripts | 14 | 1 | 72 | 2 | 48 |
| Share of the transcripts in cells with zero nuclear transcripts | 2.3 % | 1.6 % | 2.0 % | 1.3 % | 4.9 % |
| Transcripts below Q20, nuclear | 12.2 % | 12.1 % | 11.7 % | 17.3 % | 10.6 % |
| Transcripts below Q20, cytoplasmic | 12.7 % | 13.1 % | 12.1 % | 17.5 % | 10.3 % |
| Correlation over cells: nuclear share of counts with nucleus share of the cell area | 0.84 | 0.81 | 0.85 | 0.89 | 0.81 |
| Variance of r of a target gene explained by the cell's area ratio (median over genes) | 7 % | 4 % | 6 % | 14 % | 7 % |

- Cells with zero nuclear transcripts: n = 0 for every gene, so every detected gene scores positive. Every
  cell without a segmented nucleus is such a cell. 137 further cells have a small nucleus (median 8 to 12 square
  micrometres by sample, against 35 to 56 over all tumor cells with a nucleus) with no transcript in it; they
  are low-count cells (median 35 to 126 transcripts against 173 to 1432).
- Q20: the share of low-quality transcripts is the same in the two compartments, so it does not shift the
  balance (cytoplasmic share of the transcripts changes by at most 0.3 points). It adds noise.
- Geometry: the cell-wide lean of the score follows the area ratio of nucleus and cell. For a single gene the
  area ratio explains 4 to 14 % of the variance of r. This is the cell-level signal seen in the reliability
  check of §10.7. The correlation is over the cells of the analysis with at least 50 transcripts; the variance
  explained is per gene of the primary list, over its detected cells.

**C. Effect on the gene selection.** The gene filter and the V-abs ranking were computed for four variants,
all on the protein-coding genes. The reference is all tumor cells with all transcripts.

| | All tumor cells | Cells without nuclear transcripts removed (the analysis) | Q20 | Both |
|---|---|---|---|---|
| Genes kept by the filter | 2965 | 2967 | 2961 | 2962 |
| Of these, in the filter of the reference | 2965 | 2964 | 2951 | 2951 |
| Primary list of 500: genes shared with the reference | 500 | 461 | 435 | 425 |
| List of 1000: genes shared | 1000 | 946 | 888 | 869 |
| Spearman of the V-abs rank with the reference | 1 | 0.985 | 0.968 | 0.952 |
| Cytoplasm-leaning genes in the list of 500 | 390 | 413 | 388 | 414 |

- Per gene the statistics hardly move: the Spearman correlation with the mean of r of the reference is 0.997
  or more in every sample and variant, with its variance 0.94 to 0.997. The lists change at their margin,
  where many genes have nearly the same variance.
- Removing the cells lowers the mean score of every gene slightly (median change -0.02 to -0.05) and matters
  most in Prostate, which has the most such cells.
- Q20 lowers the variance of r by 3 to 4 % (median).
- The script checks that its removal variant reproduces the gene filter and the target ranking of the analysis
  exactly.

**Decisions (user, 2026-11-18).** The cells with zero nuclear transcripts are removed; Q20 is not applied
(reasons in §8). Steps 0 to 3 were rerun. Afterwards the genes were restricted to protein-coding genes and the
primary list size set to 500, and steps 0 to 3 rerun again.

**What the reruns changed.**

| | Run 1: all tumor cells, all panel genes | Run 2: cells without nuclear transcripts removed | Run 3 (current): also protein-coding genes only |
|---|---|---|---|
| Cells | 623,838 | 596,389 | 596,389 |
| Genes in the universe | 5001 | 5001 | 4972 |
| Genes kept by the filter | 2965 | 2966 | 2967 |
| Measurable genes among the kept | 2457 | 2445 | 2436 |
| Genes above the null in 3 of 5 samples | 2387 | 2361 | 2350 |
| V-null and V-abs: genes shared at 1000 | 775 | 766 | 770 |
| V-abs first 500: genes in common with the run before | | 461 | 498 |
| V-abs first 1000: genes in common with the run before | | 947 | 997 |

- Run 2 against run 1: the filter keeps 2963 genes in common; the list moves toward cytoplasm-leaning genes,
  because the removed cells had pushed every gene's score upward.
- Run 3 against run 2: the filter loses its 15 non-coding genes (CRNDE, CYTOR, DLEU1, DLEU2, H19, NORAD, PCAT1,
  PRNCR1, PVT1, SNHG1, SNHG15, SNHG17, SOX2-OT, UCA1 and CBSL) and gains 16 coding genes at the margin. In the
  primary list H19 and NORAD leave and ENPP1 and INHBA enter. V-abs had kept nearly all non-coding RNAs out
  already: they are nuclear-retained and vary little.
- N = 3000 stands (§10.4). The audit of step 3 is now read at size 500, where it is closer than at 1000
  (§10.7).

## Log

- 2026-11-18: document created. The advisor judged the residual model of the first NCD analysis wrong (the
  per-cell scale may be biology) and chose the ratio as the indicator. The earlier scripts, outputs and logs were
  moved to the archive folders of §9 and the earlier plan document was renamed and marked stale. Decisions of
  §8 recorded (user): working score log2((c + 1) / (n + 1)), no leakage handling, selection by variation with
  selection by mean score as alternative. No code, no results.
- 2026-11-18 (review by the user): open items settled. N of the filter follows the detection rates of step 0
  (start 2000, then 3000). Target lists of 500 to 1000 genes, set from the variation results, equal size for
  both rules. The raw ratio with hard-coded special cases (empty pair 1, cytoplasm empty 0, nucleus empty capped
  at 5 or 10), until now a sensitivity version, is dropped and will not be run. Selection by mean score is run
  in step 2 next to selection by variation; a new audit step (step 3) compares the two lists and only one goes
  on. Selection is done within samples and aggregated (3 of 5), not on pooled cells (§4.1). Steps are
  implemented strictly one at a time after review.
- 2026-11-18 (variation statistic): the variation statistic of list V is computed on the cells where the gene is
  detected, not on all cells (empty pairs are identical in data and null). Noted that single-transcript cells
  carry no information on variation between cells, so the statistic rests on cells with two or more transcripts
  of the gene; a minimum number of such cells is required. Qualification on evidence, ranking on size.
- 2026-11-18 (step 0 run): `7_ncd_ratio_counts.py` and `ncd_ratio_utils.py` written and run on the five samples.
  Results in §10. Cytoplasm-only pairs are 41 to 52 % of the non-empty pairs (single transcripts). 2459 genes
  are measurable for the variation statistic. The rule fixed before the run gave N = 2500 (2488 genes). The
  estimates from the archived tables that stood in §2 and §3 were replaced by the measured numbers.
- 2026-11-18 (N, user's review of step 0): 2000, 2500 and 3000 are all reasonable; the analysis proceeds with
  N = 3000 for now (2965 genes, 2457 measurable), so that nearly all measurable genes are candidates for the
  target list. Recorded in `ncd_ratio_utils.N_FILTER`. A minimum of detected cells was added to list M (§4.3).
- 2026-11-18 (folder layout): `code/1_exploration/` was split into one subfolder per analysis. The NCD ratio scripts
  are in `code/1_exploration/2_ncd_ratio/`, numbered by step: `7_ncd_ratio_counts.py` became `0_counts.py`. Step 0
  was rerun from the new folder with an identical `summary.txt`.
- 2026-11-18 (step 1 run): `1_filter.py` written and run. 2965 genes kept (rule of §3 with N = 3000), 2457 of
  them measurable; results in §10.5. Dropped by design: genes of other cell types and sample-specific lineage
  genes (65 with a rank of 500 or better in one or two samples). The constants of the rule are now in
  `ncd_ratio_utils` (`N_FILTER`, `MIN_SAMPLES`, `MIN_CELLS_GE2`); step 0 reran with an identical summary. Awaits
  the user's review before step 2.
- 2026-11-18 (review of steps 0 and 1): outputs and documentation of both steps checked; the user agreed to go
  on to step 2. The filter rule stays as decided (sample-specific lineage genes remain dropped). The left-over
  quick run of step 0 (`0_counts_quick/`) was deleted.
- 2026-11-18 (step 2 design, user): variation is selected both ways, above the null (V-null) and by absolute
  variance (V-abs); the audit decides between V-null, V-abs and M on biological plausibility. List M uses the
  mean score centred within sample and no additional threshold (the minimum of detected cells planned earlier
  is dropped: every filtered gene has at least 760 detected cells in 3 samples). Aggregation for all lists: in
  the sample's top K in at least 3 of 5 samples.
- 2026-11-18 (step 2 run): `2_selection.py` written and run. 2387 genes exceed the null in 3 of 5 samples, with
  a small excess (median 5 to 11 % by sample). V-null and V-abs share 775 of their first 1000 genes; M shares
  384 and 233 and is dominated by nuclear-leaning, often rarely detected genes. Results in §10.6. Awaits the
  user's review before step 3.
- 2026-11-18 (review of step 2, user): M is a reference only in the audit, which decides between V-null and
  V-abs. List size 1000, with the first 500 genes of the ranking as a nested sensitivity check in step 4
  (`ncd_ratio_utils.LIST_SIZE`, `LIST_SIZE_CHECK`). The user approved step 3.
- 2026-11-18 (step 3 run): `3_audit.py` written and run. The statistics and the ranking rule of step 2 were
  moved to `ncd_ratio_utils` (`pair_stats`, `entry_ranks`, `ranks_desc`, `third_best`, `top_genes`); a quick
  run of step 2 gave byte-identical tables before and after. Results in §10.7: 775 genes shared; only V-null is
  a nuclear-leaning set led by non-coding RNAs, only V-abs holds surface and signalling genes; V-null is less
  tied to the expression level (837 against 700 genes kept under matching). Claude recommends V-abs. Awaits the
  user's choice before step 4.
- 2026-11-18 (choice of the list, user): V-abs is the target list (1000 genes, first 500 as nested check).
  V-null and M are not carried further. Recorded in `ncd_ratio_utils.TARGET_LIST`; `3_audit.py` writes the list
  as `3_audit/target_genes.csv` and `.npy`. Step 4 is not started.
- 2026-11-18 (audit of steps 0 to 3, on the user's request): scripts, outputs and documents checked against
  each other. No error in the statistics of any step. Fixed: (1) in `3_audit.py` the expression-matched rebuild
  of list M did not follow the rule of step 2 (direction agreement applied after ranking); corrected, the number
  of M genes kept changes from 817 to 882, nothing else changes; (2) the reproducibility row of M now takes the
  direction from the rule of step 2 (same numbers); (3) `0_counts.py` no longer reads the 545-gene file of the
  archived analysis (one row less in `retention.csv`); (4) formatting of the audit summary and figure. Steps 0
  to 3 were rerun: every table is byte-identical to the previous run except the three items above. Document:
  §3 and §4 rewritten as the exact procedure, §8 as a table of decisions with their rationale; corrected "mean
  whole-cell count" (the filter ranks by cytoplasmic plus nuclear count) and the old name "list V" (V-null).
- 2026-11-18 (audit of the compartment matrices, on the user's request): `check_upstream.py` written and run.
  The stored matrices equal an independent rebuild from the raw 10x transcripts in every entry, and with Q20
  their sum equals the 10x count matrix. Three properties of the measurement recorded in §2 and §10.8: 2.7 to
  8.2 % of the tumor cells have no segmented nucleus, 10 to 17 % of the transcripts are below Q20, and the
  nucleus is a 2D outline so that the cell-wide lean follows cell geometry (correlation 0.81 to 0.89).
  Corrected: the sentence in §3 and §10.5 that transcripts overlapping the nucleus belong to neither
  compartment was wrong; c + n is every transcript of the cell. Open for the user: removal of the cells
  without a nucleus and Q20.
- 2026-11-18 (cells of the analysis, user): tumor cells with zero nuclear transcripts are removed (the cells
  without a segmented nucleus plus 137 with an empty small nucleus; the user chose the count rule over the
  nucleus-area rule). Q20 is not applied. `ncd_ratio_utils.load_counts` drops the rows; the matrices are not
  rebuilt. Steps 0 to 3 and `check_upstream.py` were rerun and §2, §3, §4, §8 and §10 rewritten from the new
  outputs (596,389 cells, 2966 filtered genes, 947 of the 1000 target genes unchanged; table in §10.8). No
  conclusion changed. The reason for the list size was reworded (the excess along the V-null ranking falls
  gradually; no plateau up to 1000).
- 2026-11-18 (genes and list size, user): (1) protein-coding genes only, applied before the expression filter
  (4972 of the 5001 panel genes; Ensembl biotype, release 116, through mygene.info because the Ensembl REST
  service was not reachable; `gene_biotype.py`, `data/_utils/gene_biotype.csv`; `ncd_ratio_utils.load_counts`
  drops the other columns). (2) Primary list size 500, the first 1000 as the nested check, 250 documented as
  a possible further check; reason: the selected genes are the later prediction target, and weakly informative
  genes make the prediction harder. (3) N stays 3000 after a check that a stricter N only deletes
  subset-expressed state genes from the list (§10.4). Steps 0 to 3 and `check_upstream.py` rerun; §2, §3, §4,
  §8 and §10 rewritten from the outputs. The primary list changed by 2 genes against the run before. The audit
  of step 3, now read at 500, is closer between V-abs and V-null than at 1000; V-abs was kept as decided and
  the point is put to the user.
- 2026-11-18 (target list at size 500, user): V-abs first 500 stays the primary list. Step 4 gets two checks:
  the V-abs first 1000 (nested) and the union of the V-abs and V-null first 500 (630 genes). `3_audit.py`
  writes the union (`target_genes_union.csv`, `.npy`); `ncd_ratio_utils.load_target_genes` takes "primary",
  "larger" or "union". Step 3 rerun with all existing tables unchanged.
- 2026-11-19 (second audit of steps 0 to 3 and of the upstream check, on the user's request): all scripts
  read again, every number of this document compared with the output tables, all figures looked at. No error
  in any statistic, filter or list. Corrected in this document: the reason why CD274, PROM1, FOLR1 and MSLN
  are dropped (§10.4: expressed in one or two samples only, not "under 1 % of the cells"); the nucleus area
  of the 137 cells with an empty nucleus (§10.8: median 8 to 12 square micrometres); "63 %" for "about two
  thirds" and "95.0 to 99.5 %" for "95 to 99 %" (§8); the percentiles of the cytoplasmic fraction (§10.5);
  the naming of the nuclear flag (§1, §2); the union files added to the file table of §10.7. Code, no rule
  or constant changed: `check_upstream.py` calls its reference variant "reference" (columns
  `..._with_reference`, before `..._with_current`) and computes the geometry part on the cells of the
  analysis (before: all tumor cells with a nucleus; same numbers at the precision given here); the figure of
  step 2 shows the lists of 500; docstrings, labels and run times corrected. Steps 0 to 3 and
  `check_upstream.py` were rerun: every csv and npy of steps 0 to 3 is byte-identical to the run before.
- 2026-11-19 (step 4 implemented, user's go-ahead): `4_embedding.py` with the SLURM array `4_embedding.sh`
  (24 tasks: three gene lists, two corrections, four runs) and the local review `4_review.py`; Makefile targets
  `ncd-push` and `ncd-pull`. Decisions of the user: the score is corrected within each sample before pooling;
  both corrections (z-score only, and depth regression then z-score) are run and audited; data, null and two
  transcript halves per list; the heavy part runs on HGCC, submitted and pulled by the user; the audit of list
  and correction comes right after step 4. Procedure in §5. Tested on a subsample only (all 24 tasks and the review on 4000 cells per sample, one task on 20,000
  cells per sample: 5 min; the test output is deleted). No result yet.

