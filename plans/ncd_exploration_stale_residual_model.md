# Nucleocytoplasmic discordance (NCD) of the target genes: substates and stress (residual model, STALE)

> **Stale since 2026-11-18. Superseded by `plans/ncd_ratio_exploration.md`.** The advisor judged the model of
> this analysis wrong: the per-cell scale kappa (§3.1) may be biology, and removing it leaves an uninformative
> residual. The NCD analysis restarts with the cytoplasm-to-nucleus ratio. The numbers below are correct for the
> residual model; its conclusions (§6 to §8) are not the current position of the project. The text below is
> unchanged and still uses the old locations. Scripts are now in `code/1_exploration/old/ncd_residual/` (do not
> run from there), outputs and run logs in `output/1_exploration/old/ncd_residual/`. This file was
> `plans/ncd_exploration.md`.

Final as of 2026-11-15; the stage is closed. All steps and the two closing tests are run for two gene sets: the
634 target genes and the 545 NCD-variable genes selected on the full panel. **The result is negative for both:
the NCD does not organise tumor cells into reproducible substates and carries almost no stress information.**
At population level a few dozen genes shift their compartment split with the environment by 5 to 27 %, but the
shifts follow cell geometry and the expression of neighbouring cells and form no program (§5.7). §5 holds the
results, §6 the outcome, §7 the discussion of why, §8 the consequence for the model output. §5.8 holds the
figures that show the result next to the cytoplasmic analysis, and one qualification they brought: the
full-depth NCD embedding contains a small cell-cycle structure (proliferating cells split by where their
cell-cycle transcripts sit), which is real, is not a stress state and does not change the outcome. The analysis was
audited on 2026-11-15 (code, outputs and this document) and again on 2026-11-18 (this document, the outputs and
the descriptions in the code against each other). Follows the maintenance rule of `CLAUDE.md`: sections
describe the latest state, the Log is append-only.

## 1. Purpose and scope

Source: the advisor's sketch `plans/AIVC_ideas_2026-11-12.docx` ("Nucleocytoplasmic discordance reveals spatially
organized stress responses in tumors"). Its claims are that NCD is a structured, stress-specific signal, that it
identifies cell states that nuclear or whole-cell clustering misses, that each stress axis has a distinct NCD
signature, that NCD explains stress better than nuclear expression, and that a model can predict NCD from nuclear
expression and be applied to snRNA-seq atlases.

The exploratory stage established two results for the cytoplasmic expression of the 634 target genes
(`plans/subtype_exploration.md` §5.1, `plans/stress_substate_exploration.md` §7): it reveals substates shared by
five tumors, and the substates are linked to stress. This stage asks whether the same holds when the cytoplasmic
expression is replaced by an NCD score per gene and cell:

1. Do NCD substates exist: reproducible groups of cells in the NCD space, shared by the samples?
2. Are they, or the NCD itself, linked to tumor stress, primarily environmental stress?

Expectation, clarified by the user on 2026-11-14: the NCD substates need not be the cytoplasmic substates cell by
cell. Similar categories that are linked to stress are enough. Agreement with the cytoplasmic labels is reported,
never required.

Scope, fixed by the user: the sample-corrected pooled analysis only (no per-sample clustering, no uncorrected
pool). Tumor cells of BC, OC, CC, LC and Prostate (623,838), the stress table of the stress stage, and two gene
sets:

| Gene set (`--gene-set`) | Genes | Origin |
|---|---|---|
| `target634` | 634 (142 leak-prone) | Cytoplasm-enriched highly variable genes, the target set of the earlier stages |
| `ncd545` | 545 (none leak-prone) | Genes with reliable NCD in >= 3 samples on the full panel (§5.5), leak-prone genes removed. Primary new set (user, 2026-11-15) |
| `ncd786` | 786 (241 leak-prone) | The same selection with the leak-prone genes. Defined, not run |

## 2. What the earlier stages said about NCD

- **Fig. 1 exploration** (`plans/fig1_exploration.md` §5.3): per-cell localisation is mostly counting noise
  (reliable for 15 to 26 genes per sample, more in LC and Skin); the environment adds almost nothing to the
  per-gene compartment split and shows in cytoplasmic counts and at pathway level instead.
- **Tracker target form B** (compartment split given the total) is the NCD in model form. The naive model failed
  on it for the SG genes.
- **Nuclear and cytoplasmic expression of the 634 genes disagree cell by cell** (ARI 0.05), which is the raw
  material of a discordance, but depth and biology are not separated in that number.
- **The 634 genes are cytoplasm-enriched by selection**, so the nuclear side of their NCD is the sparse side.

Consequence for the design: every result needs a counting-noise null and a noise ceiling, and "NCD has structure"
must be separated from "NCD inherits the structure of abundance".

## 3. Definition of the NCD score

### 3.1 Model

Cell c of sample s, gene g. k = cytoplasmic count (`overlaps_nucleus == 0`), u = nuclear count
(`in_nucleus == 1`), n = k + u. M_c and Q_c are the cytoplasmic and nuclear depths over the full panel; M'_c and
Q'_c the same over the genes outside the gene set.

    k_cg | n_cg ~ Binomial(n_cg, p_cg),    logit p_cg = alpha_c + kappa_c beta_gs + eta_cg

- All three terms are unknown and estimated from the counts. The model is bilinear (kappa_c multiplies
  beta_gs): a logistic regression across genes for each cell given beta, and across cells for each gene given
  the cell terms; the fit alternates between the two.
- beta_gs: baseline localisation of gene g in sample s (centred; the sample correction).
- alpha_c (**cell shift**) and kappa_c (**cell scale**, mean 1): the cell's own balance over the genes of the set,
  and its contrast between cytoplasm-heavy and nucleus-heavy genes. Fitted by maximum likelihood together with
  beta within the sample (`ncd_utils.fit_cell_shift`), with weak ridge priors (SD 2 and 1) that keep cells with
  very few transcripts in a compartment finite. A blurred compartment assignment (cytoplasm above or below the
  nucleus counted as nuclear) lowers alpha and compresses kappa together, which is the pattern in the data.
- **eta_cg is the gene-level NCD**: is this gene more cytoplasmic in this cell than the cell's other genes
  predict, given its baseline. 0 = as expected, +0.7 = the cytoplasm-to-nucleus ratio is doubled. It is not
  fitted (one value per pair, mostly on a single transcript). The model is fitted with eta = 0, and the score
  S1 (§3.2) is the standardised test statistic of eta_cg = 0: about eta times the square root of the pair's
  information, a z-like value and not an effect size.
- The residuals of a cell are orthogonal to 1 and to beta, so neither a cell-wide shift nor a cell-wide
  compression can enter a gene or program score. Both are kept as cell-level scores of their own (S0).
- Leak-prone genes (§4) are outside the fit and scored as alpha_c + beta_gs, without the scale.

The model conditions on n, the gene's whole-cell count. The score is therefore free of the gene's abundance in
the cell, which is the property that makes it a discordance and not a second view of expression. The identity
log[(k/M) / (u/Q)] = logit(k/n) - logit(M/(M+Q)) shows that the ratio of normalised expression, the difference
of log-normalised expression and the offset-corrected cytoplasmic fraction are the same quantity when counts are
large; they differ only in how they treat zeros and counting noise.

How the model was reached is in §5.1 and §5.2 (a non-target offset, then a leave-one-out offset, then the
per-cell shift, then the scale).

### 3.2 Scores

| Score | Level | Definition | Role |
|---|---|---|---|
| S0 cell shift and cell scale | cell | alpha_c - o_c with o_c = log((M'_c + 0.5) / (Q'_c + 0.5)), centred at the sample median, and kappa_c; standard errors from the inverse of the cell's information matrix I_c = sum_g w x x', x = (1, beta_g), w = n p0 (1 - p0) | Cell-level NCD; candidates for a technical origin |
| S1 residual NCD | cell x gene | r_cg = (k - n p0) / sqrt(w (1 - w x' I_c^-1 x)), p0 = expit(alpha_c + kappa_c beta_gs); 0 when n = 0; clipped at +-10. The factor is the leverage of the fitted cell parameters | **Primary.** Input of the embedding |
| S3 program NCD | cell x program | R_cP = sum_{g in P} (k - n p0) / sqrt(W_P - l_P' I_c^-1 l_P), W_P = sum w, l_P = sum w x over the program; NA below 5 program transcripts in the cell | **Primary readout** for the stress link. Always next to matched random gene sets |

S1 is zero where there is no information and has unit variance under "partitioned as expected", so the (cell,
gene) pairs without a transcript (88 % over the scored cells of the 634 set, 70 % in LC) neither vote nor add
variance.

Programs for S3: the curated programs of `subtype_utils.py` restricted to the gene set, the top 30 marker genes
of each pooled cytoplasmic substate (from the 634-gene DE tables, so restricted to shared genes), and for the new
sets the SG markers (`data/_utils/SG_markers.xlsx`, fraction in SGs > 0.4). Programs below 5 genes are dropped.

### 3.3 Scores not used

- **S2, the difference score** d_cg = log1p(1e4 k / M_c) - log1p(1e4 u / Q_c): the literal "cytoplasmic minus
  nuclear". Not run (user decision after S1 was read). It is not abundance-free: it is, up to noise, the
  cytoplasmic matrix minus a sparse nuclear matrix, so a substate result on S2 can be the cytoplasmic result seen
  through nuclear noise. §7.4 says what it could and could not show.
- Raw ratios or fractions without the cell terms (confounded by geometry and depth), and pseudo-count log ratios
  per gene and cell (dominated by the pseudo-count).

## 4. Pipeline and controls

Sample-corrected pooled, the mirror of `4_pooled_embedding.py --correction sample`.

| Step | What | Note |
|---|---|---|
| 1. Cells | Tumor cells of the five samples. Not scored: no nuclear or no cytoplasmic transcript, M'_c < 10 or Q'_c < 10, fewer than 20 gene-set transcripts in total, fewer than 3 transcripts of the fit genes in either compartment | 81 to 97 % of the tumor cells are scored (634 genes) |
| 2. Leak-prone genes | Mean depth-normalised whole-cell expression in a stromal or immune cell type (>= 200 cells in the sample; non-malignant epithelium not used) at least 4 times the mean in malignant cells, in >= 3 of 5 samples. Outside the fit and the PCA | Uses the cell-type annotation only. 142 of the 634 genes (85 stromal, 57 immune); 1,524 of the 5,001 panel genes |
| 3. Fit | alpha_c, kappa_c and beta_gs jointly within each sample on the non-leak genes | |
| 4. Scores | S0, S1 (all genes of the set), S3 (programs; the stromal and immune leak-prone genes as two flagged leak scores when the set has such genes) | |
| 5. Scale | z-score per gene within sample, clip 10 | As the cytoplasmic run |
| 6. Check | Spearman of PC1 to PC5 with both log depths, o_c, informative genes, cell area, nucleus-to-cell area ratio, S0 and the leak scores; above 0.3 the covariate is regressed out within sample | Never triggered |
| 7. Embed | PCA(50) on the non-leak genes, kNN(30 neighbours, 30 PCs), UMAP on the pooled cells | Parameters of `2_embedding.py` |
| 8. Null | The same pipeline on k simulated from Binomial(n, p0) of the fitted model | Floor for every S1 number |

Controls that were run:

| Control | Question |
|---|---|
| Counting-noise null (step 8), through gate, embedding, clustering and stress tests | What arises from the pattern of n alone |
| Split halves by binomial thinning of both compartments (Spearman-Brown) | Reliability of a score; existence of clusters without selection on shared counts |
| Matched random gene sets: 20 per program, same size and abundance deciles | Is an S3 result specific to the program. A signal shared by all genes of a cell makes any gene set reliable |
| Geometry and depth covariates against PCs and S0; S0 also after adjustment for them | Is a discordance a segmentation property |
| Leak scores and the border substates; every contact result with and without border cells | Do neighbours' transcripts drive a result |
| Spatial block bootstrap (200 um tiles, drawn within sample) and block cross-validation | No cell-level tests |
| Cytoplasmic and nuclear scores and PCs on the same cells with the same statistic | Benchmark |
| Closing tests (§5.7): null PCs as the floor of added features; sign-flip of block scores; whole-cell abundance of the same gene with the same statistic; the gene's expression ratio in neighbouring cell types | Chance level and confounds of the population-level test |

Planned and not run: the S2 arm with its permuted-nuclear control; comparisons by decile of informative genes
per cell (list of everything not run in §9).

## 5. Results

Code in `code/1_exploration/` (each script with a `.sh` wrapper for HGCC; all runs reported here were local,
`AIVC-env`). All outputs are in `output/1_exploration/ncd/`, one folder per gene set (`target634/`, `ncd545/`)
with one numbered folder per step; `README.md` there maps the folder. The paths are defined in one place,
`ncd_utils.gene_set_config`. Every step folder has a `summary.txt` (`2_embedding/`: `embedding_summary.txt` and
`embedding_summary_null.txt`; `6_figures/`: `annotation/summary.txt`).

| Step | Script | Output folder (in `ncd/<set>/`) |
|---|---|---|
| 0, 1. Gate under three models | `7_ncd_gate.py` | `1_gate/` |
| 1. Pooled embedding, real and null | `8_ncd_embedding.py` (`--null`) | `2_embedding/` |
| 2. Stress without clusters | `9_ncd_stress.py` | `3_stress/` |
| 3. Substates | `10_ncd_substates.py` | `4_substates/` |
| 4a. Gene selection on the full panel | `11_ncd_gene_selection.py` | `ncd/gene_selection/` (full panel, outside the set folders) |
| Closing test 1 (part C of step 2) | `9_ncd_stress.py` | `3_stress/` |
| Closing test 2, genes at population level | `12_ncd_gene_stress.py` | `5_gene_stress/` |
| Figures and annotation of the NCD clusters | `13_ncd_figures.py` | `6_figures/` |

§5.1 to §5.4 report the 634 genes; §5.6 reports the 545 genes next to them; §5.7 and §5.8 report both.

### 5.1 Step 0: is the NCD measurable? (gate, 634 genes)

Per sample on the scored cells: coverage, baseline, S1, per-gene overdispersion (mean squared residual; 1 =
counting noise), split-half reliability of genes, programs and PCs, the counting-noise null (5 simulations),
eigenvalues of the gene correlation matrix of S1 against the largest null eigenvalue, matched random gene sets.
Run first under two simple models: `nontarget` (logit p = o_c + beta_g) and `target_loo` (the cell's balance over
its other genes of the set as a plug-in offset).

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Tumor cells scored | 90.1 % | 97.0 % | 86.2 % | 95.7 % | 81.7 % |
| (cell, gene) pairs with n >= 1 / n >= 2 / both compartments | 9.9 / 2.7 / 1.2 % | 13.8 / 3.6 / 2.0 % | 7.6 / 1.8 / 0.8 % | 29.9 / 13.5 / 7.4 % | 10.2 / 2.6 / 1.3 % |
| Informative genes per cell, median | 58 | 81 | 43 | 191 | 59 |
| Median per-gene overdispersion | 1.12 | 1.06 | 1.09 | 1.12 | 1.08 |
| Slope of the non-target offset (1 = transfers) | 0.82 | 0.83 | 0.87 | 0.83 | 0.88 |
| PCs above the null, nontarget / target_loo | 6 / 6 | 7 / 5 | 6 / 5 | 8 / 8 | 4 / 3 |
| Reliability of "all genes of the set", nontarget / target_loo | 0.60 / 0.53 | 0.66 / 0.52 | 0.47 / 0.61 | 0.86 / 0.76 | 0.57 / 0.63 |
| Random gene sets, mean reliability, nontarget / target_loo | 0.25 / 0.18 | 0.18 / 0.11 | 0.22 / 0.16 | 0.32 / 0.14 | 0.20 / 0.15 |

- **Coverage.** 7 to 14 % of the pairs carry any transcript outside LC, 0.8 to 2 % carry both compartments.
- **There is structure beyond counting noise, and most of it is a cell-wide shift.** The set of all 634 genes has
  reliability 0.47 to 0.86, random gene sets reach 0.18 to 0.32, and the median excess of a program over its
  random sets is +0.02. PC1 is this shift. The non-target balance does not transfer one to one to the genes of
  the set (slope 0.82 to 0.88), and the plug-in offset removes the shift only in part.
- **Leakage axes.** A stromal axis (POSTN, COL5A1, SULF1, MMP11, CTHRC1, THY1) in all five samples and a myeloid
  axis (CD68, CYBB, MSR1, ITGB2) in BC, OC, CC and LC: neighbours' transcripts in the expansion-based cytoplasm.
- The gate written before the run (>= 1 PC above the null in >= 3 samples; >= 2 programs with reliability >= 0.1
  in >= 3 samples) passes under every model and reading, but it was too easy, because it did not anticipate the
  shift. The user accepted the gate with three design changes: per-cell terms fitted jointly with the gene
  baselines, leak-prone genes out of the fit, matched random gene sets in every program result.

### 5.2 Step 1: cell model and pooled embedding (634 genes)

With the shift alone the first residual axis was still the shift, loading on the genes in proportion to beta
(Spearman -0.86 in LC). Gene-specific loadings on the shift and gene-specific slopes on o_c and depth did not
remove it; the cell scale kappa_c does. 492 genes enter the fit and the PCA.

| Model `cell_shift` | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Cells scored (share of tumor cells) | 91,149 (89.2 %) | 155,273 (96.9 %) | 188,872 (85.3 %) | 42,611 (95.5 %) | 77,298 (81.0 %) |
| Null calibration: median overdispersion / program reliability | 1.006 / 0.01 | 1.005 / 0.00 | 1.009 / 0.01 | 1.000 / 0.00 | 1.004 / 0.01 |
| Median per-gene overdispersion, real | 1.09 | 1.04 | 1.08 | 1.06 | 1.06 |
| Cell shift: SD / median SE / split-half reliability | 0.42 / 0.26 / 0.38 | 0.33 / 0.20 / 0.43 | 0.42 / 0.31 / 0.18 | 0.34 / 0.12 / 0.69 | 0.40 / 0.24 / 0.34 |
| Cell scale: SD / median SE / split-half reliability | 0.55 / 0.68 / 0.27 | 0.55 / 0.71 / 0.23 | 0.51 / 0.82 / 0.18 | 0.64 / 0.45 / 0.64 | 0.54 / 0.74 / 0.26 |
| PCs above the null (also split-half reliable) | 3 (3) | 1 (1) | 2 (2) | 5 (5) | 1 (0) |
| PC1: eigenvalue / null maximum / reliability | 1.32 / 1.19 / 0.25 | 1.39 / 1.28 / 0.18 | 1.16 / 1.13 / 0.13 | 1.91 / 1.27 / 0.54 | 1.26 / 1.26 / 0.09 |
| Count-weighted gene reliability / random gene sets | 0.13 / 0.19 | 0.08 / 0.10 | 0.11 / 0.18 | 0.11 / 0.12 | 0.11 / 0.15 |

- **The cell model works.** The top PCs are uncorrelated with shift, scale, o_c, depth, informative genes and the
  leak scores (|rho| <= 0.11); the null is calibrated; no stromal or myeloid axis remains.
- **The cell shift is real and partly geometric.** Reliability 0.18 to 0.69; it correlates with o_c (-0.22 to
  -0.41), nuclear depth (up to 0.45), cell area (0.10 to 0.25) and the nucleus-to-cell ratio (up to 0.20).
- **What remains is small outside LC.** 1 to 3 PCs just above the null in four samples, 5 in LC. The axes: cell
  cycle (BIRC5, TOP2A, MYBL2, TK1, TPX2), EEF1G against the rest (EEF1G is 4 to 16 % of the fit-gene transcripts
  and has its own reliable NCD), ER / secretory (PDIA4, PDIA3, CDH1, HYOU1), interferon in LC (an interferon
  axis also exists in OC, below the null).
- **Programs above their matched random sets in >= 3 samples:** S-phase and the ER stress marker set (5 samples),
  the S-phase and G2/M marker sets, G2/M, hallmark interferon gamma. Not above: Hypoxia, Stress response,
  TNFA/NFKB, p53, Apoptosis, Complement. Single genes have reliable NCD of their own, and a random set inherits
  it, so random sets are read against the per-gene reliability, not against 0.

Pooled embedding: 555,203 cells x 634 genes.

| | Real | Null |
|---|---|---|
| PCA variance of PC1 to PC5 (1 = one gene) | 1.183, 1.142, 1.108, 1.068, 1.059 | 1.159, 1.058, 1.057, 1.057, 1.055 |
| Variance explained by 30 PCs | 6.72 % | 6.69 % |
| Largest \|rho\| of PC1 to PC5 with a covariate | 0.19 (cell scale) | 0.11 (cell scale) |
| R2 of sample on PC1 to PC10 | 0.00 | 0.00 |

Samples mix completely and no covariate is regressed out. **The pooled NCD space is close to the counting-noise
floor**: about four PCs exceed the null by 0.01 to 0.08 gene-units. The null has one structural axis of its own
(the residuals of a cell sum to zero, so the most abundant gene is anticorrelated with the rest); real PC1
(EEF1G) is mostly that axis.

### 5.3 Step 2: NCD against stress, without clusters (634 genes)

Scored cells joined to the cached stress table (all found). Scores and axes are rank-standardised within sample.

- Part A: cell shift, cell scale, leak scores and S3 of every program against the 28 axes (the 20
  expression-free ones are the fair evidence), the environmental and the expression archetype and the
  cytoplasmic substate. Effect: r for axes, standardised mean difference d for groups. **Reproducible** = |r| >=
  0.03 or |d| >= 0.05, block-bootstrap CI excluding 0, same sign in >= 3 samples and, for programs, above the 95th
  percentile of the matched random gene sets. The thresholds are low on purpose, because S3 is mostly noise.
- Two gene arms: leak-prone genes out (primary) and all genes in.
- Part C: ridge regression of each axis on a feature set, 5 spatial block folds.

| Part C, features | Median block-CV R2, expression-free axes (20) | Signature axes (5) |
|---|---|---|
| 50 NCD PCs | -0.0005 | 0.0002 |
| 50 null-embedding PCs | -0.0007 | -0.0002 |
| Cell shift + cell scale | -0.0005 | 0.0004 |
| Geometry and depth | 0.006 | 0.012 |
| 50 nuclear PCs | 0.069 | 0.229 |
| 50 cytoplasmic PCs | 0.091 | 0.234 |
| Cytoplasmic + NCD PCs | 0.094 | 0.245 |

The NCD PCs carry no stress information beyond the null (+0.0002, CI 0.0002 to 0.0003). Added to the cytoplasmic
PCs they give +0.002.

| Part A: program, stress variable | NCD, leak genes out | NCD, all genes | Cytoplasmic score | Nuclear score |
|---|---|---|---|---|
| Hypoxia, hypoxia signature (r) | 0.003 | 0.005 | 0.58 | 0.53 |
| Hypoxia, vessel distance (r) | 0.018 | 0.015 | 0.15 | 0.13 |
| Hypoxia, Core isolation archetype (d) | 0.012 | 0.009 | 0.27 | 0.25 |
| Interferon, IFNG signature (r) | 0.003 | 0.004 | 0.30 | 0.28 |
| Interferon, Immune attack archetype (d) | -0.006 | 0.006 | 0.10 | 0.09 |
| Stress response, heat shock (r) | -0.016 | -0.016 | 0.21 | 0.17 |
| TNFA/NFKB, TNFA/NFKB signature (r) | -0.017 | -0.014 | 0.44 | 0.41 |

- Of 1,160 program tests on axes and archetypes in the primary arm, 10 are reproducible (all-genes arm 16 of
  1,280); the cytoplasmic and nuclear scores of the same programs reach 422 and 375 of 920.
- No stress program has a reproducible NCD association with its own axis or archetype, in either arm. The ten
  of the primary arm: four with the nucleus-to-cell ratio (S-phase, G2/M, Autophagy, r -0.03 to -0.04), one with
  a signature axis, and five with archetypes (Interferon lower in Stromal contact, hallmark interferon gamma
  higher in Core isolation, hallmark Complement and the Myeloid border markers lower in Crowding, Complement
  lower in the interferon expression archetype; |d| 0.05 to 0.10, intervals reaching down to 0.01 to 0.03). None
  of the five is found again in the 545 set: four are tested there and are not reproducible (|d| at most 0.04),
  and the Myeloid border marker set does not exist in that set.
- **Selection artefact, not evidence:** program NCD is higher in the cytoplasmic substate of the same name
  (d 0.10 to 0.49). The substate labels are defined on the cytoplasmic counts of the same genes, so cells selected
  for high cytoplasmic counts have a positive residual through shared counting noise.
- **The leak-score control failed as designed:** the leak scores do not follow stromal or myeloid contact
  (largest |r| 0.02). Leakage adds transcripts, so it shows in a gene's count and not in its split given the
  count. At the level of the cell the leak-prone genes are therefore not the danger for a conditional score
  that step 0 suggested, and the two gene arms agree. At population level the leakage is visible in the split
  (§5.7).
- **Cell shift:** geometry first (cell area r 0.15, nucleus-to-cell ratio 0.11, crowding -0.06, all five
  samples). After adjustment for geometry and depth small associations remain (hypoxia signature -0.05, heat
  shock -0.04, Crowding archetype +0.08); border cells have a lower adjusted shift (d -0.25). None registers in
  part C.

### 5.4 Step 3: NCD substates (634 genes)

**Existence, by split halves.** Both compartments are thinned into halves; half 1 is embedded and clustered
(Louvain 1.0, the setting of the cytoplasmic run), and each cluster's gene profile (mean scaled S1) is compared
between half 1 and the independent half 2. Labels shuffled within sample give the floor. A cluster is a
**reproducible NCD substate** if all four hold (criterion agreed with the user on 2026-11-15, replacing a
correlation-only criterion that every real component passes, however small):

1. >= 1 % of the cells in >= 3 samples;
2. profile correlation between the halves >= 0.3 and above the shuffled-label maximum;
3. not a single-gene axis: the correlation stays >= 0.3 without the 5 genes of the largest half-1 effect;
4. multi-gene effect: >= 10 genes with a half-2 effect >= 0.05 SD and the half-1 sign.

A cluster has a program-level profile (used for naming) if a curated program's half-2 mean is >= 0.05 SD in
absolute value and above the 95th percentile of its matched random gene sets.

Result for the 634 genes: 30 half-1 clusters, 24 of them with >= 1 % of the cells in >= 3 samples. **None passes.**

| Part of the criterion | 634 genes, the 24 clusters of sufficient size |
|---|---|
| Profile correlation between the halves | 0.29 to 0.85; 23 pass the earlier correlation-only criterion |
| The same without the 5 largest genes | 0.18 to 0.51, median 0.32 |
| Genes replicated in half 2 (10 required) | 0 to 3, median 1 |
| Largest program mean in half 2 (0.05 SD required) | 0.003 to 0.012 SD; no cluster has a program-level profile |

Nearly every profile is one gene against the rest. The leading genes include EEF1G, LDHA, ANXA2, HSPA8, TUBB,
MYH9, CD47, CSNK1A1, APP, KIF5B, and NDRG1 against SLC2A1 in both directions (full list in `4_substates/summary.txt`).
These are abundant genes with a reliable NCD of their own. The one profile with several genes is an ER cluster
(PDIA4, GRN, PDIA3 lower; 3 replicated genes, 10 required). One cluster of 0.4 % of the cells replicates 17
genes; it fails the size requirement and its correlation without the largest genes is 0.15.

**Full-depth clusters against stress** (Louvain on the graphs of step 1; null clusters = floor, cytoplasmic
substates on the same cells = benchmark, which reproduces its published values).

| | NCD clusters | Null clusters | Cytoplasmic substates |
|---|---|---|---|
| Clusters at resolution 1.0 (0.5) | 11 (2) | 9 (2) | 10 labels |
| Environmental archetype: Cramer's V | 0.021 | 0.020 | 0.151 |
| Environmental archetype: AMI | 0.0006 | 0.0006 | 0.031 |
| Reproducible environmental enrichments | 0 | 0 | 16 |
| Expression archetype: V / AMI | 0.035 / 0.0013 | 0.040 / 0.0017 | 0.237 / 0.049 |
| Reproducible expression enrichments | 1 | 2 | 18 |
| R2 of the expression-free axes beyond sample, median | 0.0010 | 0.0011 | 0.042 |
| ARI / AMI with the cytoplasmic labels (reported only) | 0.005 / 0.013 | 0.002 / 0.002 | |

The NCD clusters equal the null clusters on every stress measure. The one reproducible enrichment of an NCD
cluster (Hypoxia expression archetype, log2 0.65) is matched by two in the null clusters: which genes of a cell
carry transcripts depends on expression, and the null keeps that pattern.

The overall indices hide one small structure: two of the 11 full-depth clusters are a cell-cycle pair with no
counterpart in the null (annotation in §5.8). It did not form a cluster at half depth, so it is absent from the
split-half test above.

### 5.5 Step 4a: NCD-variable genes on the full panel

The 634 genes were selected for cytoplasmic enrichment and whole-cell variability, not for variable
localisation. Per sample, all 5,001 panel genes: cell model on the non-leak panel genes, per-gene split-half
reliability of S1 (two thinnings) and the counting-noise null. Rule fixed beforehand: reliability >= 0.1 with
>= 200 informative cells in >= 3 of 5 samples.

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Cells scored | 96,283 | 155,705 | 208,605 | 42,923 | 85,129 |
| Genes with a reliability | 4,863 | 4,961 | 4,698 | 4,489 | 4,433 |
| Median reliability, real / null | 0.089 / 0.00 | 0.048 / 0.00 | 0.077 / 0.00 | 0.072 / 0.00 | 0.067 / 0.00 |
| Median overdispersion, real / null | 1.058 / 1.001 | 1.032 / 1.001 | 1.055 / 1.003 | 1.053 / 1.000 | 1.048 / 1.001 |
| Genes with reliability >= 0.1 | 1,998 | 219 | 1,419 | 968 | 814 |

- **786 genes pass** (441 in three samples, 271 in four, 74 in five); 545 are not leak-prone.
- **252 of the 634 target genes pass**, 2.4 times the expectation.
- Enriched among the selected genes (fold): TNFA/NFKB 2.3, EMT/stroma 2.4, Myeloid border 3.1, SG markers 1.8,
  Hypoxia 1.8, S-phase 1.7, G2/M 1.7, interferon gamma 1.6, mTORC1 1.7, ER stress 1.5. Leak-prone genes are not
  enriched (1.08 and 0.94).
- Caveats: reliability follows how often a gene is observed (Spearman 0.42); OC contributes few genes.
- Check: per-gene reliabilities agree with the 634-gene gate (Spearman 0.90 to 0.97 per sample).

### 5.6 Step 4b: the 545 NCD-variable genes through steps 0 to 3

Differences from the 634 run, all forced by the gene set: no leak-prone genes, so all 545 genes are in the fit
and the PCA, there are no leak scores and no second gene arm. Complement / acute phase, Stress response, Lysosome
and Myeloid border fall below 5 genes and are dropped, as are the Stroma border and Myeloid border marker sets;
SG markers are added (27 programs). The cytoplasmic and nuclear benchmarks are those of the 634 genes.

| Question | Number | 545 genes | 634 genes |
|---|---|---|---|
| Cells scored | | 544,789 | 555,203 |
| Per-gene reliability, median per sample | gate, cell model | 0.08 to 0.15 | 0.06 to 0.12 |
| PCs above the null per sample (BC, OC, CC, LC, Prostate) | gate | 3, 1, 1, 7, 1 | 3, 1, 2, 5, 1 |
| Structure beyond counting noise, pooled | variance of 30 PCs, real / null | 6.36 % / 6.34 % | 6.72 % / 6.69 % |
| Stress information without clusters | block-CV R2 of the NCD PCs on the expression-free axes, real / null | 0.0002 / -0.0002 | -0.0005 / -0.0007 |
| | cytoplasmic / nuclear PCs, same cells | 0.090 / 0.068 | 0.091 / 0.069 |
| | cytoplasmic + NCD over cytoplasmic | +0.0011 | +0.002 |
| Program associations | reproducible, with an environmental or expression archetype | 0 of 324 | 5 of 348, \|d\| 0.05 to 0.10 |
| | reproducible, with an expression-free axis | 5 of 540, all \|r\| <= 0.04 | 4 of 580, all \|r\| <= 0.04 |
| Substates exist | half-1 clusters passing the four-part criterion | 0 of 29 | 0 of 30 |
| | passing the earlier correlation-only criterion | 27 of 29 | 23 of 30 |
| Substates linked to stress | V with the environmental archetype, NCD / null clusters / cytoplasmic substates | 0.022 / 0.020 / 0.150 | 0.021 / 0.020 / 0.151 |
| | reproducible environmental enrichments, NCD / null / cytoplasmic | 1 / 0 / 16 | 0 / 0 / 16 |

- **Gate.** Null calibrated in all samples. Program reliabilities are higher than in the 634 set (median per
  sample 0.12 to 0.22 against 0.10 to 0.18) but equal to their matched random sets (0.11 to 0.21): selecting genes for reliable NCD makes every gene set reliable. Specific in >= 3 samples:
  the S-phase and ER stress marker sets only.
- **Embedding.** Samples mix, no covariate reaches the threshold (largest |rho| 0.11). Pooled eigenvalues 1.174,
  1.133, 1.102, 1.079 against 1.156, 1.101, 1.056, 1.055 in the null.
- **Stress without clusters.** The NCD PCs are above the null embedding with a CI that excludes zero (+0.0005),
  at 0.5 % of what the cytoplasmic PCs carry. The five program associations with expression-free axes:
  Interferon with myeloid and immune proximity (r 0.03 to 0.04, 4 samples), G2/M with the nucleus ratio (-0.03),
  one marker set with cell area. No stress program follows its own axis (Hypoxia with the hypoxia signature
  0.004 against 0.58 for the cytoplasmic score; TNFA/NFKB -0.01 against 0.44).
- **SG markers.** Their NCD follows the SG transcript axis (r 0.12, all five samples, random sets 0.03) and the
  granule count (0.06). Four other programs have small associations with the two granule axes (|r| 0.03 to
  0.04). The SG-marker association is the largest NCD association with any axis in either gene set. It is partly
  mechanical: granules are detected from cytoplasmic transcripts, partly of the same genes, so a cell with more
  SG-marker transcripts in granules has a cytoplasmic shift of these genes by construction. It shows that the
  score detects a real relocalisation; it is not a link to environmental stress.
- **Substates.** No cluster passes. Without the five largest genes the profile correlation is 0.06 to 0.50
  (median 0.33); replicated genes per cluster 0 to 5 (median 1; at most 3 in clusters of sufficient size); best program mean 0.003 to 0.013 SD. The
  clusters are single genes with their own reliable NCD. The one enrichment (a cluster of 0.3 % of the cells,
  Immune attack archetype, log2 0.63, five samples) is below the size required of a substate.
- Check: the per-gene reliabilities of this gate agree with the selection run (Spearman 0.84 to 0.92).

### 5.7 Closing tests (both gene sets)

Two tests proposed by the discussion and approved by the user. Their rules were fixed before the run.

**Test 1. Is the cytoplasmic advantage over nuclear the partition?** Part C of step 2 with four more feature
sets: nuclear + cytoplasmic PCs, nuclear + NCD PCs, and nuclear or cytoplasmic + null PCs as the floor of 50
added features. Rule: the partition explains the advantage if the NCD PCs add to the nuclear PCs at least 10 %
of what the cytoplasmic PCs add, beyond the floor.

| Median block-CV R2, expression-free axes | 634 genes | 545 genes |
|---|---|---|
| Nuclear PCs | 0.069 | 0.068 |
| Cytoplasmic PCs | 0.091 | 0.090 |
| Nuclear + cytoplasmic PCs | 0.114 | 0.112 |
| Nuclear + NCD PCs | 0.073 | 0.071 |
| Nuclear + null PCs | 0.069 | 0.068 |
| Gain of cytoplasmic over nuclear alone (CI) | +0.047 (0.046 to 0.050) | +0.047 (0.046 to 0.050) |
| Gain of NCD over the floor, on nuclear (CI) | +0.0038 (0.0034 to 0.0040) | +0.0026 (0.0024 to 0.0029) |
| Share of the cytoplasmic gain | 8 % | 6 % |
| Gain of NCD over the floor, on cytoplasmic (CI) | +0.0023 (0.0021 to 0.0026) | +0.0011 (0.0010 to 0.0013) |

- **Rule not met.** The NCD PCs add a real but small amount to the nuclear PCs, 6 to 8 % of what the cytoplasmic
  PCs add.
- **The two compartments add to each other.** Nuclear + cytoplasmic carries 0.114, more than either alone, and
  the nuclear PCs add +0.027 to the cytoplasmic ones. This is the behaviour of two noisy measurements of the same
  quantity, abundance. The cytoplasmic advantage is therefore mostly a second, deeper read of abundance.
- Limit: the NCD PCs are near the noise floor, so this test cannot exclude partition information that is lost in
  counting noise at the level of the cell. Test 2 covers that side.

**Test 2. Does the NCD of single genes follow stress at population level?** `12_ncd_gene_stress.py`. The cell
model is fitted as before. For gene g and stress variable a (a 0/1 archetype indicator or an axis in SD units)
the one-step estimate of delta in logit p = logit p0 + delta a is sum (k - n p0)(a - abar) / sum w (a - abar)^2,
with numerators and denominators summed over the samples. For an archetype it is the difference of the gene's
logit between the archetype and the other cells (0.1 = cytoplasm-to-nucleus ratio about 10 % higher, relative to
the cell's other genes); for an axis the change per SD. Counting noise averages out over cells and the model
conditions on n, so abundance cannot bias the estimate.

- Inference from the spatial blocks: cluster-robust standard error, z, Benjamini-Hochberg q within the primary
  family (environmental archetypes and the 20 expression-free axes) and the secondary family (expression
  archetypes, signature and granule axes). **Reproducible** = q < 0.05, CI excluding 0, same sign in >= 3
  samples, and |effect| >= 0.10 (archetype) or 0.05 per SD (axis).
- Nulls: counts simulated from the fitted model and refitted (calibration); random signs on the block scores,
  100 draws (reproducible pairs expected by chance).
- Context: the gene's whole-cell abundance with the same statistic; programs pooled over their genes against
  matched random gene sets; the gene's expression ratio in stromal and immune cells over malignant cells; all
  tests again without the border substates.
- Rule: positive if the reproducible pairs of the primary family exceed the chance level and are not explained
  by the neighbour control.
- Two changes from the approved plan: the cluster-robust standard error replaces block resampling, and the
  sign-flip null replaces a permutation of stress variables between blocks (not well defined for cell-level
  variables in blocks of different size).

| Genes that are not leak-prone | 634 genes (492) | 545 genes |
|---|---|---|
| SD of z, real / simulated null (1 = counting noise) | 1.63 / 1.05 | 1.86 / 1.05 |
| Reproducible pairs, primary family | 51 (26 genes) | 193 (73 genes) |
| of these: archetypes / cell area and nucleus ratio / other axes | 10 / 10 / 31 | 16 / 44 / 133 |
| Expected by chance (sign-flip null, 95th percentile; simulated null) | 0; 0 | 0; 1 |
| Also reproducible without the border substates | 37 | 149 |
| Size of the reproducible effects, median (largest) | 0.065 (0.16) | 0.076 (0.27) |
| Share of reproducible pairs by the gene's neighbour ratio below 1 / 1 to 2 / 2 to 4 | 0.3 / 0.3 / 0.7 % | 0.5 / 1.0 / 2.4 % |
| Abundance of the same genes: reproducible pairs, primary family | 2,597 | 2,501 |
| Programs above their matched random gene sets, primary family | 0 | 0 |

- **The count is far above chance, so gene-level shifts of the compartment split exist.** The first half of the
  rule is met. The effects are small: a 5 to 27 % change of a gene's cytoplasm-to-nucleus ratio.
- **The second half of the rule is not met: most of the pairs have a technical explanation.**
  - Geometry: a fifth of the pairs are with cell area or the nucleus-to-cell ratio (IRAK1, BIRC5, CDKN3 and
    other cell-cycle genes, more nuclear in cells with a large nucleus).
  - Neighbours' transcripts: the rate of reproducible pairs rises 2 to 5 times with the gene's expression in
    stromal and immune cells, below the 4-fold leak rule. The strongest genes are of this kind: UBD, CXCL10 and
    CXCL11 with immune and myeloid proximity. The leak-prone genes of the 634 set show it openly: 355 of their
    3,512 primary tests are reproducible, and 75 to 78 % of the stromal and immune genes are more cytoplasmic
    with fibroblast and immune proximity. So leakage does reach the split at population level, which corrects
    the reading of the leak-score control in §5.3 for this level.
  - Coupling with abundance: across genes the partition effect follows the abundance effect on the
    environmental axes (Spearman 0.13 to 0.31; not on the two geometry axes), and where the abundance effect of
    a reproducible pair is itself reproducible (31 of 51 pairs, 634; 114 of 193, 545) the signs agree in 87 and
    83 %. A gene is more cytoplasmic where it is more expressed. This is what
    transcripts of neighbouring tumor cells of the same state produce in an expansion-based cytoplasm, and the
    opposite of the induction lag. A negative coupling, as the lag predicts, appears on the hypoxia, TNFA/NFKB,
    interferon and SG signatures (-0.03 to -0.20) and on the Unstressed expression archetype (-0.09).
- **No program.** The pooled NCD of a stress program on its own variable is at most 0.03 (Hypoxia with vascular
  deprivation 0.02 to 0.03, not above its random sets; TNFA/NFKB in its expression archetype -0.03). SG markers
  with the SG transcript axis 0.08 is the one program above its random sets, and it is partly mechanical (§5.6).
- **What remains as candidates.** Outside geometry, 10 genes (634) and 19 genes (545) have a reproducible pair
  and a neighbour ratio below 2, among them P4HA1 (more cytoplasmic far from vessels, where it is also
  induced), MRAP2 with crowding, CCDC113 and NEK2. They cannot be separated from tumor-to-tumor
  spill-over with this segmentation. Lists in `reproducible_pairs.csv` of each output folder.
- Checks: the simulated null is calibrated (5.9 to 6.1 % of |z| above 1.96, at most one call). On the 201 genes
  shared by the two sets the effects agree (Spearman 0.975), although the cell model is fitted on different
  genes. The abundance benchmark finds the known response (2,500 pairs).

**Reading.** Test 2 is positive by count and negative by attribution. Gene-level shifts with the environment are
real, a few percent to 27 %, gene-specific, and mostly explained by geometry and by transcripts of neighbouring
cells. Nothing in them supports a stress-specific relocalisation program.

### 5.8 Figures and annotation of the NCD clusters (both gene sets)

`13_ncd_figures.py` draws the figure types of the cytoplasmic exploration for the NCD, with the cytoplasmic
reference and the counting-noise null in every comparison and on the colour scale of the cytoplasmic figure. It
reads stored objects only (no fit, embedding or clustering). Added 2026-11-15 at the user's request, after the
stage was closed.

**Annotation of the full-depth clusters** (Louvain 1.0; initial run, audit, refined run, as for the cytoplasmic
clusters). Markers = genes with a mean scaled NCD at least 0.10 SD above the other cells; label from the
hypergeometric overlap with the programs and from the program NCD scores (thresholds of `subtype_utils`). The
same procedure runs on the null clusters. A program name is kept only with at least 3 program markers, a
stronger overlap than any null cluster, a counterpart in the split-half test and at least 1 % of the cells in
3 samples. Two tests say whether a cluster is one of the cytoplasmic substates: same cells (fold of a
cytoplasmic substate above the largest fold of a null cluster, pooled and in 3 samples) and stress
(reproducible environmental enrichment). Table: `6_figures/annotation/ncd_cluster_annotation.csv`, reviewed
labels with reasons in `REVIEWED` of the script.

| | 634 genes | 545 genes |
|---|---|---|
| Real / null clusters | 11 / 9 | 12 / 10 |
| Clusters named by a single gene (EEF1G, LDHA, HSPA8, NDRG1 vs SLC2A1, CD47, VEGFA, SOST, ...) | 8 | 9 |
| Null clusters named by a single gene | 9 | 10 |
| Largest program overlap of a null cluster | 4 genes | 5 genes |
| Clusters with an automatic program name that meets all conditions | 0 | 0 |
| Clusters with a reviewed descriptive label | 3 | 3 |

Reviewed clusters:

| Cluster | Cells | What it is | Substates of its cells (fold: cytoplasmic / nuclear label) |
|---|---|---|---|
| 634 N08, G2/M genes cytoplasmic | 2.0 %, 5 samples | 22 G2/M markers (BIRC5, KPNA2, CCNB1, TOP2A at +0.7 to +0.9 SD), program z 1.0 (null at most 0.27) | G2/M 6.1 / 1.1; 59 % of the cluster |
| 634 N07, G2/M genes nuclear | 5.1 % | the same genes at -0.5 to -0.7 SD | G2/M 1.3 / 4.6 |
| 634 N04, ER genes nuclear | 13.4 % | APP, PDIA3, PDIA4, CDH1, GRN at -0.3 to -0.6 SD | ER stress 0.5 / 2.5 |
| 545 N11, S-phase genes cytoplasmic | 0.3 %, below 1 % in every sample | 20 S-phase markers (TYMS, MYBL2, BIRC5, DHFR) | S-phase 3.0 / 1.4 |
| 545 N05, cell-cycle genes nuclear | 7.5 % | TUBB, TYMS, BIRC5, DHFR, MYBL2 at -0.3 to -0.9 SD | S-phase 1.2 / 2.0 |
| 545 N09, ER genes cytoplasmic | 5.0 % | 7 ER stress markers and membrane genes (ATP2A2, DSG2, CDH1, TFRC) | ER stress 1.3 / 0.6 |

**Reading of the annotation.**

- Most NCD clusters are single-gene clusters, as in the null. Their names would carry no meaning.
- The cell-cycle pair is real: many coherent genes, absent from the null, present in all samples. Both
  clusters of the 634 set are proliferating cells (abundance of the S-phase and G2/M genes 0.6 to 0.8 SD above
  the other cells; in the 545 set 0.6 to 0.8 SD for N11 and 0.2 to 0.3 SD for N05). The NCD splits them by where these transcripts sit. This is the cell-cycle component already seen as
  a weak PC and as a specific program in §5.1 and §5.2.
- It is not the cytoplasmic G2/M substate. N08 holds 12 % of that substate. The cytoplasmic-heavy half carries
  the cytoplasmic label and the nuclear-heavy half the nuclear label, which is the count-sharing effect of
  §5.3 seen from the NCD side. A state found by NCD would be enriched under both labels.
- No name was kept as a program name: N08, the only cluster with an automatic program name, has no
  counterpart at half depth. N04 and several single-gene clusters do have one (N04: the ER cluster of §5.4).
- Same-cells test: passed by N08 (fold 6.1) and, in the 634 set, by four single-gene clusters and N07 with
  folds of 1.3 to 1.7 against a null maximum of 1.28; in the 545 set by N11 only. The small folds reflect that a gene needs
  transcripts to have an NCD (the null clusters reach 1.28 the same way), so only N08 and N11 count.
- Stress test: no cluster of the 634 set. In the 545 set N11 is enriched in Immune attack (log2 0.63, 300
  cells, 5 samples); it is 0.3 % of the cells.
- No cluster has a hypoxia, interferon or stress-response gene profile in either set. N09 of the 634 set
  (NDRG1) is enriched in the Hypoxia expression archetype (log2 0.65, the enrichment of §5.4) and holds the
  cytoplasmic Hypoxia substate at fold 1.7; it is a single-gene cluster, and two null clusters have an
  expression enrichment of the same kind.

**Consequence.** The statement "no reproducible NCD substates" (§5.4, §5.6) stands for the split-half test and
needs one qualification: at full depth the NCD shows a small cell-cycle structure (2 % and 5 % of the cells in
the 634 set). It is a cell-cycle effect on transcript location, plausibly the timing of synthesis and export or
mitosis, and it has no environmental link. The stress conclusion and outcome D of §6 are unchanged.

**Figures** (per gene set, in `output/1_exploration/ncd/<set>/6_figures/`).

| File | What it shows | What to see |
|---|---|---|
| `overview.png` | Rows cytoplasmic, NCD, null; UMAP by own clusters, by cytoplasmic substate, by environmental archetype; cluster x archetype enrichment | The first figure to open. The cytoplasmic UMAP has islands and a coloured enrichment map; the NCD row looks like the null row |
| `clustering/umap_clusters.png`, `umap_cyto_substates.png`, `umap_cyto_substates_highlight.png` | The three embeddings by their clusters and by the cytoplasmic substates, one substate at a time | Each substate is a compact region on the cytoplasmic UMAP and is spread over the whole NCD and null UMAP |
| `clustering/umap_sample.png`, `umap_cell_scores.png`, `umap_stress.png` | NCD UMAP by sample, by cell parameters, depth and geometry, by stress archetypes and six axes | The main gradient of the NCD UMAP follows single genes and depth, not stress |
| `clustering/cluster_composition.png`, `sankey_cyto_vs_ncd.png`, `pc_variance.png` | Sample shares, flow of the cytoplasmic substates into NCD and null clusters, variance of the PCs | Every substate spreads over all NCD clusters, except G2/M into N08 in the 634 set; the NCD PC spectrum lies on the null |
| `annotation/cluster_vs_cyto_substate.png` | Cytoplasmic and nuclear substates inside the NCD and null clusters | The cell-cycle pair and the ER cluster with opposite sign under the two labels; null panels blank |
| `annotation/heatmap_program_scores.png`, `DE/heatmap_cluster_genes.png` | Program scores and top genes per group, cytoplasmic substates against NCD and null clusters | Cytoplasmic substates have blocks of marker genes; NCD and null clusters have one gene each, except the cell-cycle pair (634: N07 and N08) |
| `annotation/split_half_profiles.png`, `split_half_summary.png` | Gene profile of every half-1 cluster in the two halves; the four parts of the criterion | One dot off the origin per cluster; no cluster passes |
| `spatial/spatial_labels_<sample>.png` | Tissue maps: cytoplasmic substates, NCD clusters, null clusters, environmental archetype | Substates form territories; NCD and null clusters are salt and pepper |
| `spatial/spatial_<program>_vs_<axis>.png` | Stress axis, cytoplasmic program score and program NCD per cell and as 200 um block mean, null block mean. Three pairs; the third is Stress response vs crowding in the 634 set and ER stress vs crowding in the 545 set (Stress response has too few genes there) | The cytoplasmic score has tissue patterns that follow the axis; the block mean of the program NCD is flat like the null |
| `spatial/spatial_cell_shift.png` | Cell shift and nucleus-to-cell ratio | The cell-wide shift follows geometry |
| `stress/profile_heatmap.png`, `stress_subtype_enrichment.png` (and `_null`), `correspondence_summary.png` | Stress profiles and archetype enrichments of the three labellings; V, AMI, R2, reproducible pairs | 16 environmental and 18 expression enrichments for the cytoplasmic substates; 0 and 1 for the NCD clusters of the 634 set, 1 and 0 for the 545 set |
| `stress/sankey_cyto_vs_stress.png`, `sankey_ncd_vs_stress.png`, `boxplot_axes_by_*.png` | Flows between archetypes and labels; raw axes per label | Opaque ribbons only in the cytoplasmic figure |
| `stress/program_axis_heatmap.png` | Programs x stress axes for the cytoplasmic score, the nuclear score and the program NCD | Two coloured panels and an empty one |
| `stress/information_closing.png` | R2 per feature set; gain over the nuclear PCs (closing test 1) | NCD PCs add a sliver, cytoplasmic PCs a bar |
| `genes/volcano.png`, `partition_vs_abundance.png`, `pairs_by_class_and_neighbour_ratio.png`, `forest_top_pairs.png`, `gene_examples.png` | Closing test 2: gene effects for partition, simulated noise and abundance; coupling with abundance; hit rate by neighbour expression; strongest pairs; examples without the model | Partition effects are few and small next to abundance, rise with expression in neighbouring cell types and follow the abundance effect |

Differences from the plan of the figures: the dot plot of clusters by axes was dropped (same content as
`profile_heatmap.png`); the top-gene figure is a heat map, because the dot size of an expression dot plot has no
NCD analogue.

## 6. Reading rules and outcome

| Outcome | Evidence | Conclusion |
|---|---|---|
| A. Stress-linked NCD substates | Reproducible NCD substates (four-part criterion of §5.4), at least one with a reproducible environmental enrichment, and an overall association with the environmental archetype above the null clusters | The advisor's hypothesis holds for the gene set |
| B. NCD substates without a stress link | Reproducible NCD substates, no reproducible environmental enrichment | NCD has cell states that are not stress states |
| C. Inherited structure | Structure only for a score that carries abundance, unchanged under a permuted nuclear term | The result is abundance |
| D. No reproducible structure | Clusters do not reproduce between halves, or no more than the null clusters | No substates at this depth for these genes |

Without clusters, "linked" means reproducible associations of S3 or S0 with the expression-free axes or the
environmental archetype, and NCD PCs with a block-CV R2 above the null embedding, read against the cytoplasmic
and nuclear PCs on the same cells.

**Outcome: D for the 634 genes and D for the 545 genes.** None of the criteria uses agreement with the
cytoplasmic labels, so the outcome holds under the loose expectation of §1. The closing tests (§5.7) do not
change it: they concern single genes at population level and the source of the cytoplasmic advantage, not
substates. The cell-cycle structure of the full-depth clusters (§5.8) does not change it either: it does not
pass the split-half criterion, and with it the nearest outcome would be B (structure without a stress link),
not A.

## 7. Discussion: why the result is negative, and what follows

Written from the results of §5.1 to §5.6. The two tests it proposed are run (§5.7); §7.6 says what they changed.

### 7.1 Is the result negative under the loose criterion?

Yes, for both gene sets. The loose criterion asks for reproducible categories of cells in NCD space that are
linked to stress, whatever their relation to the cytoplasmic substates. The evidence does not use the
cytoplasmic labels:

- Clusters of the NCD embedding are as associated with the environmental archetype as clusters of pure counting
  noise (V 0.021 and 0.022 against 0.020; cytoplasmic substates 0.15).
- No cluster of substate size keeps a multi-gene profile in independent reads (at most 3 replicated genes; 10
  required).
- The NCD PCs predict the expression-free stress axes barely above the null embedding: +0.0002 (634) and
  +0.0005 (545), with intervals that exclude zero, both below 1 % of what the cytoplasmic PCs carry.
- No stress program's NCD follows its own stress axis or archetype in either set. Program associations with
  any archetype: 5 of 348 for the 634 genes (|d| at most 0.10), none of them found again among the 0 of 324 for
  the 545 genes.

What exists, for completeness: a cell-cycle and an ER / secretory component that are specific but too weak to
form reproducible clusters at half depth; at full depth they appear as a small cell-cycle pair and an ER
cluster without an environmental link (§5.8; they are not stress categories); interferon NCD with myeloid proximity at r 0.03 to 0.04;
SG-marker NCD with the SG axis at r 0.12 (partly mechanical); one cluster of 0.3 % of the cells.

### 7.2 Why "the same prediction target" does not mean "the same substates"

The equivalence of §8 holds **given the input**: cytoplasmic = nuclear x partition, so a model that receives the
nuclear (or whole-cell) expression and predicts the cytoplasmic expression is predicting the partition.

The cytoplasmic substates were found **without conditioning**, in the cytoplasmic expression itself, which is
abundance x partition. Abundance, how much of a gene the cell holds, is shared by the two compartments and is
what the substates are made of: S-phase, hypoxia or interferon cells are cells that hold more transcripts of
those genes. Evidence from a separate comparison on the same cells (§5.3, part C): the nuclear counts are a
physically separate set of transcripts of the same genes, and their PCs predict the stress axes nearly as well
as the cytoplasmic PCs (0.069 against 0.091). If the stress signal sat in the partition, the nuclear
expression would not carry it and the NCD would; the reverse is observed. (The comparison shows a similar
amount of stress information, not identical information; the substates themselves were found on the
cytoplasmic expression only.) NCD, as defined here, removes exactly the shared part and keeps the partition.

So both statements are true. For a model with nuclear input, NCD and cytoplasmic expression are the same target.
For unsupervised structure they are opposite ends: cytoplasmic expression is mostly abundance, NCD is what is
left when abundance is taken out. The negative result measures how little of the cytoplasmic signal is
partition.

### 7.3 Three reasons, with their evidence

1. **Measurement.** Per gene and cell the split is a binomial on very few transcripts: about three quarters of
   the pairs with any transcript hold exactly one (55 % in LC). Variance beyond
   counting noise is 4 to 10 % of the residual variance (median overdispersion 1.04 to 1.10 under the cell
   model, over samples and both sets). A program score reaches a reliability of about 0.2, a single gene about 0.1.
2. **Structure.** The reliable part varies gene by gene, not as programs. Program reliability equals that of
   matched random gene sets, selecting NCD-variable genes raised all gene sets alike, and clusters replicate as
   single genes (EEF1G, LDHA, CD47, NDRG1 against SLC2A1). Clustering needs genes that move together across
   cells. Abundance has that property, because transcription programs co-regulate genes. Export and stability
   are regulated largely transcript by transcript, and this data shows no program-wide relocalisation.
3. **Size of the stress link.** This is not only a power problem. Corrected for the unreliability of the score,
   the correlation of a program's NCD with an expression-free axis is at most 0.11 (634 genes) and 0.09 (545 genes), with a 95th
   percentile of 0.05, where the cytoplasmic score of a program reaches 0.27 without any correction. A stress link of NCD is
   therefore not excluded, but it is bounded at a small value, and it is not specific to stress programs (the
   same programs do not exceed their matched random gene sets). This is a bound, not a proof of absence. The
   population-level test (§5.7) confirms the bound from the other side: with counting noise averaged out, the
   pooled NCD of a stress program moves by at most 0.03 logit with its own stress variable.

Technical floor under all three: the cell-wide shift, the largest reliable NCD signal, follows cell area and the
nucleus-to-cell ratio. In a 5 um section with 2D segmentation the compartment assignment itself varies between
cells, and the cell model has to spend the cell-level signal on removing that.

### 7.4 Can the definition be improved?

| Candidate | What it could change | What it cannot |
|---|---|---|
| Difference score S2 (cytoplasmic minus nuclear log expression) | It would very likely show clusters and a stress link | It carries abundance: the clusters would be the cytoplasmic substates seen through nuclear noise (outcome C). It answers a different question and would not show discordance |
| Cytoplasmic expression conditional on the observed nuclear expression (the residual of a nuclear-to-cytoplasmic prediction) | Closest to what the model does. Would also show structure | Nuclear counts are a noisy measure of abundance, so the residual still contains abundance that the nuclear counts missed. Test 1 of §5.7 shows this is the main part: nuclear + NCD recovers 6 to 8 % of what nuclear + cytoplasmic adds |
| A coarser unit than the cell (neighbourhood, pseudobulk) | Removes counting noise (reason 1) | Not reason 2 or 3. The archetype comparisons of §5.3 are already group means over tens of thousands of cells with tight intervals, and they are near zero |
| Gene classes with a known localisation mechanism (SG markers, nuclear-retained transcripts) | The one sizeable signal is of this kind | Narrow: a handful of genes, partly mechanical for SG markers, and the panel lacks the classic retained transcripts (MALAT1, NEAT1) |
| Better data (3D segmentation, deeper panels, nascent or unspliced reads) | Reasons 1 and the technical floor | Not available in this project |

No redefinition that keeps NCD free of abundance is likely to produce substates in this data. Definitions that
produce substates do so by letting abundance back in.

### 7.5 Is the NCD direction wrong?

The advisor's sketch contains three separable claims.

- **NCD as a basis for tumor cell states.** Not supported. Two gene sets, one of them selected for NCD
  variability, give no reproducible categories. I do not see a variant of the analysis that would change this
  with the present data.
- **NCD as a stress-specific signal.** Bounded small at program and cell level. At the level of single genes
  across the population, shifts with the environment exist but are small and mostly follow geometry and
  neighbours' transcripts (§5.7). No stress-specific relocalisation is shown.
- **NCD as the model's estimand.** Unchanged in principle (§8), but the result changes the expectation. The
  partition is the only part of the cytoplasmic expression that nuclear expression does not already determine,
  and it is small and unstructured. A model evaluated on NCD should be expected to gain little over the
  constant-partition predictor, and a model evaluated on cytoplasmic expression will owe most of its accuracy to
  abundance.

So the direction is not wrong as a definition, and it is the honest form of the prediction target. It is wrong as
a source of substates in this data, and it lowers what the model can be expected to show.

### 7.6 What the closing tests added

Both tests of §5.7 were proposed here and run on 2026-11-15.

1. **The cytoplasmic advantage is mostly abundance.** Nuclear and cytoplasmic PCs add to each other like two
   reads of one quantity, and the NCD PCs recover 6 to 8 % of what the cytoplasmic PCs add to the nuclear ones.
   This supports §7.2 with a number.
2. **At population level single genes do shift, but not as a stress response.** With full power the test finds
   reproducible gene-level effects far above chance, of 5 to 27 %, concentrated in geometry axes and in genes
   expressed by neighbouring cells, and coupled to the gene's own abundance in the direction that spill-over
   predicts. No program moves. The candidates that survive the controls are listed, not claimed.
3. A third lesson is technical. In this data the compartment split of a gene is sensitive to transcripts of
   neighbouring cells at the 5 to 27 % level, which is the size of the effects one would hope to find. Any later
   use of the partition as a target has to control for it (3D or membrane-stained segmentation, or a restriction
   to transcripts well inside the cell).

The stage closes with these conclusions.

## 8. NCD or absolute cytoplasmic expression as the model output

The question raised by the advisor: NCD is more meaningful biologically, and predicting cytoplasmic expression
from whole-cell expression looks circular while predicting NCD looks independent of it. The user's view: both
amount to predicting the nucleus-cytoplasm partition. Position of this document:

1. **Given the input, the two targets are the same prediction problem.** With whole-cell expression T as input,
   C = T x p and NCD = logit p. With nuclear expression N as input, C = N x exp(NCD). The unknown is one number
   per gene and cell, the partition, and a model for one target converts into a model for the other.
2. **The circularity is in the evaluation, not in the target.** "Predicted against observed cytoplasmic
   expression" is high for a model that knows nothing but each gene's average share. NCD removes that trivial
   component, so its headline metric is honest. The same honesty is obtained for the cytoplasmic target by
   reporting the increment over the constant-partition null (M0 in the tracker).
3. **Only the conditional form is free of abundance.** Observed difference or ratio scores share reads with the
   whole-cell count (§3.3). Biologically the ratio cancels the transcription rate only at steady state; during
   induction nuclear transcripts lead cytoplasmic ones, so NCD mixes recent transcriptional change with export
   and stability.
4. **With nuclear input the absolute cytoplasmic target is the cleaner one technically.** Input and target reads
   are disjoint. An observed NCD target contains the input with a negative sign, so noise in the nuclear counts
   produces a spurious negative relation unless the nuclear reads are split between input and target.
5. **What this stage adds.** The partition component is measured: a few percent of the residual variance per
   informative pair, gene-specific, with a stress link bounded near zero at program and cell level (§7.3). The
   honest increment over the constant-partition null will therefore be small whichever target is used. With
   nuclear input a second source of apparent gain exists that is not partition: nuclear counts measure abundance
   with noise, and context can inform the model about the abundance the nuclear counts missed. Test 1 of §5.7
   puts this source at more than 90 % of what the cytoplasmic expression adds to the nuclear expression for
   stress.
6. **Conditional means of the partition are learnable but confounded.** A model learns population-level
   effects, and those exist (§5.7). In this data they are mostly geometry and neighbours' transcripts, so a
   model would learn them as if they were biology unless the evaluation controls for both.

Recommendation: keep NCD (the partition given the total) as the estimand that defines what is new, train on the
counts, and report every result as the increment over the constant-partition null. Do not present NCD as a
source of cell states or as a stress-specific signal. The partition alone is not a sufficient basis for the
model in this data: what a model could learn of it beyond the constant-partition null is small and mostly
technical.

## 9. Implementation and status

| Step | Code (`code/1_exploration/`) | Status |
|---|---|---|
| 0. Gate (§5.1) | `7_ncd_gate.py`, `ncd_utils.py` | run for `target634` and `ncd545` |
| 1. Cell model and pooled embedding (§5.2) | `ncd_utils.py`, `7_ncd_gate.py` (model `cell_shift`), `8_ncd_embedding.py` | run for both sets, real and null |
| 2. Stress without clusters (§5.3) | `9_ncd_stress.py` with `stress_utils.py` | run for both sets |
| 3. Substates (§5.4) | `10_ncd_substates.py` | run for both sets |
| 4a. Gene selection (§5.5) | `11_ncd_gene_selection.py` | run |
| Closing test 1 (§5.7) | `9_ncd_stress.py` (part C) | run for both sets |
| Closing test 2 (§5.7) | `12_ncd_gene_stress.py` | run for both sets |
| Figures and cluster annotation (§5.8) | `13_ncd_figures.py` with `stress_figures.py` | run for both sets |
| Audit | code, outputs and document (scripts 7 to 11) | done 2026-11-15 |
| Audit | document, README, tracker and docstrings against the outputs (scripts 7 to 13) | done 2026-11-18 |

The stage is closed.

Usage: every script runs from its own folder, `python <script>.py [--quick] [--gene-set <name>]`; order 7, 8, 8
`--null`, 9, 10, 12 with the same gene set. Script 11 works on the full panel, takes no `--gene-set` and runs
once, before a new set is used; 12 needs only the leak-gene table of 7 or 8. `--quick` uses 20,000 cells per
sample (script 13: 50,000 cells in total) and writes to `<step>_quick/` next to the step folder. A full run
empties the step folder it writes, with two exceptions: script 8 keeps an existing embedded object unless
`--force` is given, and script 13 with `--only <parts>` redraws the named parts and leaves the rest. Runtime of
a full chain locally: about 1.5 hours (634 genes) and 1.2 hours (545 genes); script 12 adds 7 to 9 minutes.
Script 13 runs last (about 2 minutes per set) and writes only to its own `6_figures/` folder.

Deliberately not run: the `ncd786` chain (its 241 additional genes are leak-prone; they do not change a
cell-level score, §5.3, and at population level they carry neighbours' transcripts, §5.7); the S2 arm with
its permuted-nuclear control; comparisons by decile of informative genes per cell; per-sample or uncorrected
analyses; the predictability analysis and a whole-cell reference arm (user decision of
2026-11-14).

## 10. Decisions and open items

Decided by the user:

- 2026-11-14: scope (sample-corrected pool only); S1 primary; the gate accepted with the three design changes;
  one step at a time, each reviewed; NCD substates need not match the cytoplasmic substates.
- 2026-11-15: `ncd545` primary new set, `ncd786` sensitivity; four-part existence criterion; local runs; audit
  and discussion from existing results, with new diagnostics as proposals only; then both proposed tests
  approved and run, and the stage closed. The `ncd786` chain was then left out (reason in §9).

Open:

- The consequences for the model design (target, evaluation against the constant-partition null) are decisions
  of the tracker, not of this stage.
- The candidate genes of §5.7 are not followed up. A follow-up needs a segmentation that excludes spill-over.
- EEF1G dominates the fit-gene counts of the 634 set (4 to 16 %) and one axis per sample. It does not change any
  conclusion; noted in case the set is reused.
- `stress_utils.block_resample` takes the sample from the text before the first underscore of the block name,
  which is "Xenium" for every sample, so the bootstrap of `6_stress_substates.ipynb` drew blocks across the pool
  and not within sample. Scripts 9 and 10 use their own stratified resampling and are not affected. A fix means
  re-running notebook 6.
- From the advisor's sketch, for later stages: reproducibility across sections and platforms, gene identity of
  high-NCD genes, clinical relevance.

## Log

- 2026-11-14: document created from the advisor's sketch of 2026-11-12 and the user's instructions. Scope fixed
  to the sample-corrected pooled analysis. NCD defined as the offset-corrected, sample-centred log
  cytoplasm-to-nucleus ratio (eta); binomial residual S1 primary, log difference S2 as comparison arm, program
  score S3 as readout. Position on the output metric recorded in §7. No code, no results.
- 2026-11-14 (step 0): user decisions recorded (S1 primary; no predictability analysis, no whole-cell arm; same
  634 genes). `ncd_utils.py` and `7_ncd_gate.py` written and run on the five samples. A second offset
  (`target_loo`) and matched random gene sets were added during the run. The gate passes under every reading,
  but most of the reliable signal is a cell-wide shift of the target genes (random gene sets reach reliability
  0.18 to 0.32; median program excess +0.02), followed by stromal and myeloid leakage axes and a cell-cycle axis;
  program-specific excess is small and absent for Hypoxia and Stress response. Three design changes proposed
  for the downstream steps; nothing downstream implemented.
- 2026-11-14 (after step 0): the user accepted the gate and the three design changes. Model revised to
  logit p = alpha_c + beta_gs + eta_cg with alpha_c fitted jointly (the `target_loo` plug-in offset is confounded
  by the cell's gene composition); cell shift S0 = alpha_c - o_c kept as a cell-level score; leak-prone genes
  defined from the cell-type annotation and excluded from the fit and the embedding input; matched random gene
  sets required for every S3 result. Step 1 (revised score and pooled embedding) specified in §5.2 with its
  criteria; sections 5.2 to 5.4 renumbered to 5.3 to 5.5.
- 2026-11-14 (step 1): `fit_cell_shift`, leak-gene table, matched random sets, gate model `cell_shift`,
  `8_ncd_embedding.py` and the two SLURM wrappers written; full runs local. Deviations from the approved step:
  cell scale kappa_c added next to the shift (first residual axis was the shift with loadings proportional to
  beta), random-set criterion re-read against per-gene reliability, weak ridge prior on alpha_c. 142 leak-prone
  genes. The cell model removes the shift from the residual axes in all samples and the null is calibrated.
  What remains: 1 to 3 weak PCs outside LC, 5 in LC; cell cycle, EEF1G, ER / secretory, interferon; specific
  programs S-phase, G2/M, ER stress markers, interferon gamma; none for Hypoxia or Stress response. Pooled
  embedding: samples mix, no covariate regressed, but 30 PCs explain 6.72 % against 6.69 % in the null. User
  instruction: heavy steps go to HGCC from now on.
- 2026-11-14 (step 2): on the user's choice the stress tests that need no clusters were run before the substate
  step (`9_ncd_stress.py`). NCD PCs carry no stress information beyond the null embedding (median block-CV R2 on
  the expression-free axes -0.0005 vs -0.0007; cytoplasmic 0.091, nuclear 0.069; cytoplasmic + NCD +0.002). No
  stress program has a reproducible NCD association with its own axis or archetype, with or without the
  leak-prone genes (r <= 0.02 where the cytoplasmic score reaches 0.15 to 0.58). The cell shift follows cell
  area and nucleus ratio (r 0.15, 0.11); after adjustment small associations remain (hypoxia and heat-shock
  signatures -0.05). The leak-score control failed: leakage is in the count, not in the split. Program NCD by
  same-name cytoplasmic substate (d 0.1 to 0.5) is flagged as selection on shared counts.
- 2026-11-14 (steps 3 and 4a): the user clarified that NCD substates need not match the cytoplasmic ones cell by
  cell; §6 rewritten. Step 3 (`10_ncd_substates.py`): half-1 clusters pass the pre-set profile correlation (23 of
  30) but replicate only 4 to 8 % of their profile, each as one abundant gene; full-depth NCD clusters equal the
  null clusters on every stress measure (environmental archetype V 0.021 vs 0.020; cytoplasmic substates 0.151).
  Outcome D for the 634 genes. Step 4a (`11_ncd_gene_selection.py`): 786 of 5,001 panel genes have reliable NCD in
  >= 3 samples (545 not leak-prone; 252 of the 634), enriched for TNFA/NFKB, hypoxia, cell cycle, interferon and
  SG markers. Re-running on the new set waits for the user.
- 2026-11-15 (step 4b, 545 genes): user decisions recorded (545 primary, 786 sensitivity, tighter existence
  criterion, local runs). `--gene-set` added to scripts 7 to 10; four-part existence criterion in
  `10_ncd_substates.py`. Full chain on the 545 non-leak NCD-variable genes: null calibrated; 30 PCs explain 6.36 %
  against 6.34 % in the null; NCD PCs R2 on the expression-free axes 0.0002 against -0.0002 (cytoplasmic 0.090);
  no program associated with a stress archetype; 0 of 29 half-1 clusters pass the new criterion (27 pass the old
  one); NCD clusters V 0.022 with the environmental archetype against 0.020 for the null clusters. Outcome D
  again. One sizeable association, SG-marker NCD with the SG transcript axis (r 0.12), is partly mechanical. The
  786 chain and the 634 substate rerun were not run.
- 2026-11-15 (audit and discussion): code, outputs and this document audited at the user's request. Code: no
  error that affects a result; cleanup only (one shared sample loader in place of three copies, unused names
  removed, gene-set-neutral names and docstrings, summary columns only for existing arms and programs). Quick
  runs of both gene sets gave identical tables before and after the cleanup. Both full chains were re-run with
  the final code: every table is identical to the earlier run; bootstrap intervals of the 634 substate step
  differ in the third decimal because the random draws of the new criterion precede them. The 634 substate step
  is now under the four-part criterion (0 of 30; 23 pass the earlier one). Document: restructured to the current
  state; sections renumbered (stress without clusters is §5.3, substates §5.4, gene selection §5.5, the 545 set
  §5.6, discussion §7, model output §8, status §9, decisions §10), so section numbers in the entries above refer
  to the earlier layout. One number corrected: the cytoplasmic and nuclear program scores reach 422 and 375
  reproducible associations of 920 (322 and 292 were written before). Added the five small archetype
  associations of the 634 set, which the earlier text left out. No script retired; development outputs and
  superseded logs removed. Discussion added as §7: the result is negative under the loose criterion for both
  sets; the equivalence of the two prediction targets holds given the input and does not carry over to
  unsupervised structure, because the substates are made of abundance, which NCD removes; two tests proposed,
  not run.
- 2026-11-15 (closing tests, stage closed): the two tests of §7.6 run for both gene sets, rules fixed before the
  run. Test 1 (part C of `9_ncd_stress.py` extended; all earlier numbers of the step unchanged): on the nuclear
  PCs the NCD PCs add +0.0038 (634) and +0.0026 (545) over a null-PC floor, 8 and 6 % of the +0.047 that the
  cytoplasmic PCs add; rule (10 %) not met. Test 2 (new `12_ncd_gene_stress.py`): 51 reproducible gene x stress
  pairs in 26 genes (634) and 193 in 73 genes (545) in the primary family against 0 expected by chance, effects
  of 0.05 to 0.27 logit; a fifth are geometry axes, the rate rises with the gene's expression in neighbouring
  cell types, leak-prone genes are reproducible in 355 of 3,512 tests, and the partition effect is positively
  coupled to the abundance effect; no program above its random gene sets. Read as positive by count and negative
  by attribution. Deviations from the approved plan: cluster-robust standard errors in place of block
  resampling, sign-flip null in place of a permutation between blocks. §5.7 added, §7.6 rewritten, §8 extended
  (points 5 and 6, recommendation). The NCD stage is closed.
- 2026-11-15 (figures and cluster annotation): new `13_ncd_figures.py` (39 figures per gene set, stored objects
  only) with the figure types of the cytoplasmic exploration, each next to the cytoplasmic reference and the
  null. NCD clusters annotated with the procedure of the cytoplasmic clusters (initial run, audit, refined
  run), with the null clusters as the floor. Result: 8 of 11 (634) and 9 of 12 (545) clusters are single-gene
  clusters as in the null; a cell-cycle pair exists at full depth (634: G2/M genes cytoplasmic 2.0 % and
  nuclear 5.1 % of the cells; 545: S-phase genes cytoplasmic 0.3 %, cell-cycle genes nuclear 7.5 %) and a
  weaker ER / secretory cluster; no hypoxia, interferon or stress-response cluster; no environmental
  enrichment in the 634 set, one in the 545 set (the 0.3 % cluster with Immune attack). §5.8 added; §5.4 and §6
  qualified: the cell-cycle structure is real and was missed by the split-half test; outcome D and the stress
  conclusion unchanged.
- 2026-11-15 (folders reorganized): all NCD outputs moved into `output/1_exploration/ncd/<set>/` with one
  numbered folder per step (`1_gate`, `2_embedding`, `3_stress`, `4_substates`, `5_gene_stress`, `6_figures`),
  for `target634` and `ncd545`; the earlier locations (`pooled_sample_corrected_ncd*/`, `stress_ncd*/`,
  `ncd/gate*/`, `ncd/gene_stress*/`, named in the entries above) no longer exist. Paths are defined in
  `ncd_utils.gene_set_config`; scripts 7 to 13 and their wrappers updated. Nothing was recomputed; script 13
  (`--only overview`) reproduced identical files from the new paths. Logs renamed to
  `ncd_<set>_<script>.log`. `README.md` added in the output folder.
- 2026-11-18 (audit of documents and outputs): at the user's request every number of §3 to §8 was compared with
  the stored tables, and the README of the output folder, the tracker, `CLAUDE.md` and the docstrings were
  checked against the folders and the code. All tables of this document match the outputs; no conclusion
  changes; nothing was recomputed. Corrected statements: EEF1G is 4 to 16 % of the fit-gene transcripts (4 to
  14 % was its share of all 634 genes); the myeloid leakage axis of step 0 is above the null in four samples,
  not in CC only; the interferon axis under the cell model is above the null in LC only; 27 programs in the
  545 set (28 counted the all-genes reference); program reliabilities of the 545 gate given as per-sample
  medians (0.12 to 0.22; the earlier 0.18 to 0.26 could not be reproduced); agreement with the selection run
  0.84 to 0.92; sign agreement of partition and abundance effects restated with its definition (87 and 83 %);
  the same-cells test is passed by four single-gene clusters and N07, not five single-gene clusters; N04 and
  several single-gene clusters have a half-depth counterpart, N08 has none; N09 (NDRG1) carries the Hypoxia
  expression enrichment of §5.4. Reworded: share of empty pairs and of single-transcript pairs (§3.2, §7.3),
  the nuclear-expression argument of §7.2, the NCD PCs against the null in §7.1 (same wording for both sets),
  the weak components of §7.1 (stale since §5.8), set-specific remarks in the figure table. Added to §3.1: all
  three model terms are estimated, eta is not fitted and S1 is its test statistic. §9: usage corrected (script
  11 without `--gene-set`, quick size of script 13, what a run empties), one list of what was not run. Log
  files are named `ncd_<set>_<nn>_<step>.log` (the entry above gives a shorter pattern). Also updated: README of
  the output folder, exploratory-stage sections of the tracker, `CLAUDE.md`, docstrings and wrapper comments of
  scripts 7 to 13 (no logic changed).
- 2026-11-18 (stale): the advisor judged the residual model wrong (the per-cell scale may be biology). Document
  renamed from `ncd_exploration.md` and marked stale; scripts moved to `code/1_exploration/old/ncd_residual/`,
  outputs and logs to `output/1_exploration/old/ncd_residual/`. The analysis restarts with the ratio score in
  `plans/ncd_ratio_exploration.md`. No number of this document was changed.
