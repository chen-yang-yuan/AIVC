# AIVC Model Design Tracker

Current as of 2026-11-14.

## Purpose

This tracker records the AIVC model-design brainstorm: the current state of each design component, live candidates, and the decisions locked. It supplements the repo documents `plans/AIVC_plan.md` (background, data, prior model design), and `plans/data_processing.md` (pipeline status), which keep the original design until decisions here are ready to be folded back in, and is informed by the three exploration documents (`plans/fig1_exploration.md`, `plans/subtype_exploration.md`, `plans/stress_substate_exploration.md`).

Conventions: every section except the Decision log describes the latest state only and is rewritten as thinking changes, with stale material removed rather than annotated; the Decision log is append-only and dated. A decision is only locked once it appears in the log.

## The model in one paragraph

AIVC (the only name used in the manuscript) is a predictive, interpretable model of cytoplasmic activity in tumor cells. It takes a cell's identity embedding (own expression, optionally tumor subtype), its neighborhood embedding, and, when H&E is available, its morphology embedding, and predicts the cytoplasmic expression of a selected gene set or, if per-gene signal is too sparse, the cytoplasmic activity of selected pathways. It is trained on subcellular spatial transcriptomics (Xenium 5K, CosMx 6K), which assign reads to nucleus and cytoplasm directly. Its distinguishing features are the compartment-level target and the attribution of each gene's or pathway's cytoplasmic activity to the input modalities. No perturbation data are involved. The target gene set is not restricted to stress-response genes; it will be chosen from the data, with biologically relevant genes (stress response first) preferred among those whose cytoplasmic expression is sensitive to subtype, morphology or neighbors.

## Framing

| Question | Current answer | Status |
| --- | --- | --- |
| Q1. Training data | Xenium 5K and CosMx 6K (expression, compartment-assigned reads, neighbors, H&E), FFPE samples only for now; fresh-frozen samples are deferred because gene localization patterns differ between preparations. More than 20 such datasets are available; the six processed Xenium 5K samples (BC, OC, CC, LC, Prostate, Skin) are all FFPE (verified 2026-09-26 from `experiment.xenium`). The subtype and stress explorations use the five carcinoma samples; Skin (melanoma) is excluded there. | Settled |
| Q2a. Inputs | Identity and neighborhood embeddings required; morphology optional. No perturbation embedding. Whether identity is built from total or nuclear expression is open (see Input design). | Settled except identity source |
| Q2b. Outputs | Cytoplasmic expression of a data-selected gene set; pathway-level cytoplasmic activity as the fallback if per-gene signal is too sparse. Working gene set: the 634 cytoplasm-enriched highly variable genes of the subtype exploration (recommended, not locked). | Form settled; gene set recommended, not locked |
| Q3. Framing | A predictive model of cytoplasmic activity with modality attribution. Counterfactuals (swapping the neighborhood) are a use of the trained model, not its definition. | Settled |
| Q4. Application data | Deferred until the model performs well on the training datasets. Candidates: subcellular ST with smaller panels (standard Xenium, MERSCOPE), scRNA-seq with a constructed neighborhood, and possibly Visium. H&E alone is out of scope because identity and neighborhood are required. | Deferred |

A principle that applies throughout: a target that is a deterministic function of an input is trivial to predict. Total expression of a gene contains most of its cytoplasmic count, so the treatment of the target genes on the input side (see Output design) decides what the model actually learns.

## Input design

| Embedding | Required? | Content | Source | Notes |
| --- | --- | --- | --- | --- |
| Identity | Yes | Own gene expression, optionally tumor subtype / cancer type. Source compartment open: total counts (then target genes must be excluded, or the target defined as a split) or nuclear counts (then input and target reads are disjoint by construction; see below) | scGPT or a trained encoder on expression; subtype as categorical | Cancer type may be treated as a confounder to regularize against. |
| Neighborhood | Yes | Cell-type composition and distances within a radius; aggregated expression of neighboring cells; tumor-cell crowding; boundary vs core | Neighbor graph; scGPT or Nicheformer on neighbor expression | Use the full neighborhood by default; the immune-only kernel tested in the naive model carried no signal for the SG target. |
| Morphology | Optional (data availability) | Cell-centered H&E patch | UNI2-h | Dropped at training time with some probability so the model tolerates its absence. |

Nuclear expression as the identity source (proposed). Nuclear reads approximate nascent transcription, cytoplasmic reads the mature, functional mRNA pool; building identity from the nucleus and the target from the cytoplasm separates the two by read assignment rather than by gene exclusion, so no panel gene has to be withheld from the input. It also puts identity on the better-segmented compartment (the DAPI-based nucleus) and the target on the noisier one (the expansion-based cytoplasm). Two things to check before adopting it: nuclear counts are a minority of reads, so a nuclear identity embedding is sparser than a total one, and nuclear-only clustering must recover the same tumor subtypes as total expression (a Fig. 1 analysis). One thing to keep in mind: the nuclear and cytoplasmic counts of the same gene remain strongly correlated through abundance, so a nuclear-identity-only model (M0) still learns per-gene nuclear-to-cytoplasmic ratios; that ratio is a biologically meaningful quantity (export and stability), and the neighborhood and morphology increments over M0 remain the environmental evidence. A hybrid is also possible: nuclear counts of all genes plus total counts of non-target genes.

Fusion: the three embeddings are combined by attention so that, for each output gene or pathway, the attention weights give a modality attribution (see Interpretability). Detailed architecture beyond this (frozen foundation models vs trained encoders, graph vs pooled neighborhood) remains deferred.

Existing assets: scGPT and UNI2-h weights in `data/_utils/`; trial embedding notebooks in `code/3_model/`; the exploration code in `code/1_exploration/` and `code/2_exploration/` (compartment statistics, controls, subtyping, stress axes).

## Output design

Primary target: cytoplasmic expression of a selected gene set. Fallback: cytoplasmic activity of selected pathways, computed by aggregating cytoplasmic counts over each pathway's genes (count aggregation rather than rank-based ssGSEA, which is unstable at per-cell sparsity). Pathway scores can also be derived post hoc from predicted per-gene counts if the per-gene target is retained.

Gene-set selection: the aim is genes whose cytoplasmic expression varies across cells beyond counting noise and is associated with tumor subtype, morphology or neighborhood, with adequate counts. The exploration produced the 634-gene set (highly variable and cytoplasm-enriched across samples; see Exploratory stage), which is the working target set.

Target form (open). Three definitions of "cytoplasmic activity" are compatible with the framing and differ in what the model learns and what the identity input may contain:

- A. Absolute cytoplasmic count from total-expression identity, k\_g \~ NB. The target gene's total count must be excluded from the identity embedding, otherwise the prediction is trivial. An imputation task with a compartment-level target and a neighborhood input.
- B. Compartment split given total, k\_g \~ Binomial(n\_g, p\_g) with the per-cell global offset from the naive model. Total counts may stay in the identity embedding. The model learns where a gene's transcripts sit, not how many there are. The naive model showed this signal is small for the 129 SG genes, and the Fig. 1 analysis found that the environment adds almost nothing to the per-gene split on the whole panel. The first NCD analysis tested this form directly with a residual model (`plans/ncd_exploration_stale_residual_model.md`): no substates and a stress link at the counting-noise floor for two gene sets. The advisor judged that model wrong on 2026-11-18 (the per-cell scale it removes may be biology); the test is repeated with the ratio score (`plans/ncd_ratio_exploration.md`, design only).
- C. Absolute cytoplasmic count from nuclear-expression identity, k\_g \~ NB with the nuclear counts of all genes as input. No gene exclusion needed; input and target reads are disjoint. The model learns the cytoplasmic pool given the nascent state, and the environmental question becomes whether the neighborhood and morphology change that mapping.

A and B are complementary (predicting both n\_g and p\_g gives the cytoplasmic count as their product with abundance and localization separated); C replaces the gene-exclusion device of A with read-level separation. Evidence so far (Exploratory stage and Open questions): B does not carry enough environmental signal to stand alone, and nuclear expression is a defensible identity source where nuclear depth is adequate. The subtype exploration adds a fourth element that applies to A and C alike: a per-gene, per-sample intercept, so that the predicted quantity is the sample-relative log cytoplasmic activity.

Evaluation against the imputation null. Under either form, identity-only prediction is the natural null model. Three fits are reported: M0 identity only, M1 identity + neighborhood, M2 identity + neighborhood + morphology (and neighborhood + morphology without identity as a check). The increments M1 − M0 and M2 − M1 are the evidence that the environment and morphology carry information about cytoplasmic activity beyond own expression, and they are the quantitative counterpart of the attention attributions.

### Evidence from the naive model (retained lessons)

The naive model (a predictability screen of cytoplasmic SG-gene localization from own expression and the immune
neighborhood) is no longer part of the project and its code and write-up are not kept in the repo. What it
established for the SG gene set, retained as design lessons:

- Per-cell cytoplasmic localization of the 129 SG genes is 65–100% counting noise (median 7–29 SG transcripts per cell), reliability 0.00–0.34; beyond that, own expression carries almost no information in five of six samples (LC is the exception and does not transfer), and the immune-only neighborhood carries none.
- An expression-matched non-SG control set is far more predictable than the SG set, so residual skill is generic compartment or segmentation structure, not SG biology.
- Retained: the offset-binomial target form and the evaluation discipline (reliability ceiling, matched-control gene set, leave-one-sample-out and within-sample references, synthetic-signal plumbing check). These are re-implemented in `code/2_exploration/`.

Design lessons: any per-cell compartment target must pool enough transcripts, by choosing well-expressed genes or aggregating to pathways; a gene list from one study is fragile and the target set should be validated on these data; and any new target should be screened with the same control suite (swap the target, keep the controls) before architecture work.

## Interpretability

Interpretability is a stated requirement of the model, not an afterthought. Components:

- Attention-based fusion of the identity, neighborhood and morphology embeddings, with per-output attention weights read out as a modality attribution: which modality drives the cytoplasmic activity of each gene or pathway, in each cell and averaged across cells, samples and cancer types.
- Ablation increments (M0, M1, M2 above) as the model-agnostic check that attention attributions agree with what removing a modality actually costs.
- Within-modality attribution where useful: which neighborhood features (cell types, distances) or which identity genes drive a given output, by gradient or attention over the neighborhood graph.
- Counterfactual readouts as a downstream use: swapping the neighborhood embedding and reading the change in cytoplasmic activity, with pseudo dose-response along a neighborhood axis as the display.

## Benchmarking

Benchmark on the Xenium 5K and CosMx 6K datasets against methods that solve adjacent tasks; no existing method predicts compartment-level expression from these inputs, so comparisons adapt each method to the target.

- Gene imputation from own expression (SpaGE, Tangram, gimVI, scGPT imputation): the identity-only null, adapted to predict cytoplasmic counts.
- Spatial-context expression prediction (niche-conditioned models such as Nicheformer-style or GNN-based spatial expression predictors): the identity + neighborhood comparison.
- Simple statistical baselines: a GLM on the same features, and a per-gene offset-binomial GLM.

Metrics: likelihood-based (NB / binomial deviance explained) and R² on log counts, always relative to each sample's reliability ceiling; leave-one-sample-out and leave-one-cancer-type-out as the generalization estimates.

## Scope and data

### Scope statement

AIVC is a methodology paper, not a foundation model. Its contribution is the model design, the compartment-level target and evaluation pipeline, and the findings on the cancer types represented in the training data. Pan-cancer generalization is presented as potential, supported by leave-one-cancer-type-out evidence where it exists, and released weights are described as covering the trained cancer types and platforms only.

### Data strategy

A nested protocol:

1. Partition the available datasets once into development sets and 2–4 discovery sets (chosen for biological interest, from cancer types with at least two datasets). Discovery sets are excluded from every training, selection and tuning step.
2. Credibility: leave-one-dataset-out, leave-one-cancer-type-out and leave-one-platform-out cross-validation across development sets only. Design decisions rest on this.
3. Final model: trained on all development sets.
4. Discovery: apply the final model to the discovery sets. Observed data validates the factual predictions; the findings are the modality attributions, neighborhood counterfactuals and per-sample heterogeneity, which the observed data cannot show directly.
5. Release: final model, code and the evaluation harness. Code and the target-definition pipeline are the primary asset.
6. Platform: Xenium 5K and CosMx 6K share only part of their panels; the shared gene set, or a panel-agnostic encoder, decides whether one model spans both.

### Application datasets (deferred)

Work on these starts only after the model performs well on the training datasets.

- Subcellular ST with smaller panels (standard Xenium, MERSCOPE): the identity encoder must tolerate a different and smaller gene set. Options: train on the panel intersection; gene-masking augmentation during training; a gene-token encoder that is panel-agnostic by construction (a pretrained foundation model is the extreme case, to be treated as a candidate initialization rather than assumed superior).
- scRNA-seq: no neighborhood exists, so one must be constructed (sample-level cell-type composition as a pseudo-neighborhood, or the model run with the neighborhood treated as missing), and there is no compartment ground truth, so the output is a prediction to be validated indirectly.
- Visium: spot-level, no compartments, no single cells; a possible application of the neighborhood component only. To be discussed.
- Out of scope: H&E alone, because identity and neighborhood are required inputs.

## Manuscript plan (six figures)

| Figure | Content |
| --- | --- |
| 1 | Observations on cytoplasmic vs nuclear expression in the Xenium 5K datasets that motivate the work: which genes and pathways show differential compartment patterns, and how those patterns relate to subtype, morphology and neighborhood. |
| 2 | Overall workflow: inputs, target definition, model, attribution. |
| 3–4 | Benchmarking and cross-validation on Xenium 5K and CosMx 6K. |
| 5–6 | Applications to standard Xenium and/or scRNA-seq datasets. |

## Exploratory stage (completed 2026-11)

Four exploratory analyses of the processed Xenium 5K samples are complete; each has its own plan document with
design, results and caveats. Their conclusions enter this tracker as recommendations; nothing below is locked until
it appears in the Decision log.

| Analysis | Samples | Document | Code / outputs |
| --- | --- | --- | --- |
| Fig. 1: nuclear vs cytoplasmic expression, artifact controls, association with the model inputs | six (BC, OC, CC, LC, Prostate, Skin) | `plans/fig1_exploration.md` | `code/2_exploration/`, `output/2_exploration/` |
| Target gene set and tumor-cell substates: per sample, pooled, sample-corrected pooled (cytoplasmic and nuclear), gene-set sensitivity | five carcinomas (Skin, a melanoma, excluded) | `plans/subtype_exploration.md` | `code/1_exploration/1_cytoplasmic/` (`1`-`5`), `output/1_exploration/` |
| Stress vs substates: 28 stress axes, two stress labellings, cytoplasmic vs nuclear | the same five | `plans/stress_substate_exploration.md` | `code/1_exploration/1_cytoplasmic/6_stress_substates.ipynb`, `output/1_exploration/stress/` |
| Nucleocytoplasmic discordance (NCD), residual model: substates and stress link of the compartment split given the total, for the 634 target genes and 545 NCD-variable genes (closed 2026-11-15; superseded 2026-11-18, the advisor judged the model wrong) | the same five | `plans/ncd_exploration_stale_residual_model.md` | archived: `code/1_exploration/old/ncd_residual/`, `output/1_exploration/old/ncd_residual/` |
| NCD as the cytoplasm-to-nucleus ratio (restart 2026-11-18; steps 0 to 3, the counts check, the gene filter, the candidate gene lists and their audit, are run; tumor cells with nuclear transcripts, protein-coding genes; primary target list V-abs with 500 genes; step 4, the pooled embedding for three gene lists and two corrections, is implemented and waits for the HGCC run; in progress) | the same five | `plans/ncd_ratio_exploration.md` | `code/1_exploration/2_ncd_ratio/`, `output/1_exploration/ncd_ratio/` |

The design rules of the Fig. 1 stage (depth matching, reliability ceiling against counting noise, positive controls
for compartment assignment, matched-control gene sets, tumor cells only) are specified in `plans/fig1_exploration.md`
section 3 and remain the evaluation discipline for later stages.

### Exploration results (recommendations, not locked)

- Target gene set: 634 genes, highly variable in >= 2 and cytoplasm-enriched in >= 3 of the five carcinoma samples
  (`output/1_exploration/enrichment/overlap_genes_3plus_samples.npy`). The set is stable to the HVG basis; a
  detection-rate criterion is not a workable substitute. It replaces the ~40-gene first list of the Fig. 1 stage as
  the working target set (16 of those 43 genes are in it).
- Structure of the target: the cytoplasmic expression of these genes organises tumor cells into substates. S-phase,
  G2/M, Interferon and Hypoxia recur across the five tumors once the per-sample mean is removed; ER stress and
  Complement / acute phase partly; the two border labels (Stroma, Myeloid) are neighbours' transcripts captured by
  the cytoplasm segmentation, i.e. leakage that a neighbourhood input predicts trivially. "Stress response" lumps
  different programs and has to be split before it is used as a target.
- Sample effect: a per-gene mean shift concentrated in lineage genes; without a per-sample correction the pooled
  analysis clusters tumors, with it the samples mix. Recommended outcome: the sample-relative log cytoplasmic
  activity eta, with a count link `log mu = log d + a_gs + eta` (cytoplasmic depth offset, per-gene per-sample
  intercept) and a program-level primary readout with a per-gene track (`plans/subtype_exploration.md` section 7).
  Between-sample level differences are conditioned on, not modelled.
- Identity source: the nuclear expression of the 634 genes does not reproduce their cytoplasmic substates cell by
  cell (ARI 0.05), so a nuclear identity input (target form C) is non-redundant with the outcome; nuclear depth
  limits BC, CC and Prostate. The identity-only null needs a split-half reliability ceiling.
- Environment: the cytoplasmic substates are connected to tumor stress coherently but weakly (Cramer's V <= 0.40,
  AMI <= 0.06), and cytoplasmic expression is more closely related to the expression-free stress axes than nuclear
  expression (19 of 20 axes). The advantage is concentrated in cells with few nuclear transcripts; the
  thinned-cytoplasm depth control has not been run.
- NCD, residual model (superseded 2026-11-18: the advisor judged the model wrong because the per-cell scale it
  removes may be biology; the analysis restarts with the ratio score, `plans/ncd_ratio_exploration.md`, and the
  result below holds for the residual model only). The compartment split of a gene given its whole-cell count, target form B: negative for the 634 target
  genes and for 545 genes selected for NCD variability on the full panel. No reproducible NCD substates; clusters
  of the NCD embedding are as associated with the environmental archetype as clusters of counting noise (V 0.02
  against 0.15 for the cytoplasmic substates); the NCD PCs carry under 1 % of the stress information of the
  cytoplasmic PCs. The only structure at full depth is a small cell-cycle pair (proliferating cells split by
  where their cell-cycle transcripts sit) without a stress link. Added to the nuclear PCs, the NCD PCs recover
  6 to 8 % of what the cytoplasmic PCs add, so the cytoplasmic advantage is mostly a second read of abundance.
  At population level a few dozen genes shift their split with the environment by 5 to 27 %, mostly explained
  by cell geometry and neighbours' transcripts. Consequence proposed in
  `plans/ncd_exploration_stale_residual_model.md` section 8:
  the partition alone is not a sufficient basis for the model; report increments over the constant-partition
  null whichever target is used.
- Open checks before locking: split-half reliability of the outcome per sample, the per-gene test of whether the
  sample intercept is needed, the thinned-cytoplasm control, the split of the Stress response label, and whether the
  set and substates hold in Skin.

## Ideas under discussion

- Decomposition framing: present the model as partitioning cytoplasmic activity into identity-driven, neighborhood-driven and morphology-driven components, so that a small environmental increment is a quantified finding rather than a failure.
- Cancer type as confounder: hold out by sample; add a pan-cancer regularizer or adversarial cancer-type head.
- Per-sample reporting of compartment signal against each sample's reliability ceiling, so heterogeneity across specimens becomes a finding.
- Stress-invariant identity: if identity and target genes co-express strongly, train the identity encoder to be uninformative of the target residual (adversarial term), or restrict identity genes to lineage and housekeeping sets.

## Open questions

- [ ] Target form A (absolute count, target genes out of identity), B (compartment split given total) or C (absolute count from nuclear identity), or a combination. Fig. 1 evidence (`plans/fig1_exploration.md` 5.3): the environment adds almost nothing to the per-gene split (B) but shows in cytoplasmic counts and at pathway level; recommendation: count-based primary target (A or C) with B as a secondary readout. Subtype exploration (`plans/subtype_exploration.md` 7): count link with cytoplasmic depth offset and a per-gene, per-sample intercept; the prediction is the sample-relative log cytoplasmic activity. NCD is under re-examination with the ratio score (`plans/ncd_ratio_exploration.md`, design only, 2026-11-18), after the advisor judged the residual model wrong; the following holds for the residual model only. First NCD analysis (`plans/ncd_exploration_stale_residual_model.md` 5 to 8, closed 2026-11-15): form B was tested directly as the NCD of two gene sets and is negative as a source of substates and of stress information; the cytoplasmic advantage over nuclear expression is mostly a second read of abundance, and the partition adds 6 to 8 % of it. Recommendation of that document (section 8, not locked): keep the partition given the total as the estimand, train on the counts, and report every result as the increment over the constant-partition null.
- [ ] Identity source: total or nuclear expression. Fig. 1 evidence: nuclear-only clustering recovers the total partition as well as cytoplasmic clustering does, but clusterability is depth-limited; nuclear identity is defensible in LC, Skin, OC and marginal in BC, CC, Prostate. Subtype exploration: on the 634 target genes the nuclear and cytoplasmic substates of a cell agree at chance level (ARI 0.05), so nuclear identity is non-redundant with the outcome.
- [ ] Target gene set and whether the primary target is per gene or per pathway. Fig. 1 evidence: reliable per-gene localisation for 15-209 genes per sample; ~40 recurrent candidate / cyto-only-DE genes; environment signal mainly at pathway level (KRAS up, EMT, ER stress, TNFa). Recommendation: pathway-level primary target with a per-gene track. Subtype exploration: working target set of 634 genes; program-level readout on the recurrent substate programs (S-phase, G2/M, Interferon, Hypoxia, ER stress, Complement / acute phase), which are defined on the data rather than taken from the GMT; the Fig. 1 per-gene list overlaps the 634 only partly (16 of 43).
- [x] FFPE status of the six processed samples: all FFPE (2026-09-26). The remaining datasets' status is still to be checked as they are ingested; fresh-frozen samples are excluded from training for now.
- [ ] Tumor-cell definition across datasets: BC/OC/CC use 10x supervised labels, LC/Prostate/Skin manual cluster
  annotation without a non-malignant epithelial class. Step 0b of the Fig. 1 exploration (run 2026-09-26) found the supervised
  tumor sets clean and benign-like or low-signal clusters in 20-28 % of the LC, Prostate and Skin tumor cells;
  harmonisation (marker QC, CNV inference) is still to be decided.
- [ ] Neighborhood definition and radius; which features enter the embedding.
- [ ] Attention design that yields clean per-output modality attributions (per-gene attention heads vs shared fusion with per-gene readout).
- [ ] Which datasets become discovery sets, and whether one model spans Xenium and CosMx.
- [ ] Uncertainty estimation: likelihood-based heads vs a separate method.

## Decision log

| Date | Decision |
| --- | --- |
| 2026-09-20 | Tracker created. Training data: Xenium 5K and CosMx 6K, because all modalities are present. Cohort-level application is limited to about one figure and does not drive the design. |
| 2026-09-20 | Individual stress granules are not the smallest unit of the outcome because of sparsity; cytoplasm-level quantities replace them. |
| 2026-09-20 | Detailed architecture (frozen foundation models vs GNN, attention vs concatenation) deferred until the biological question and input/output structure are fixed. |
| 2026-09-22 | Naive model complete with a pre-registered FAIL: per-cell cytoplasmic localization of the 129 SG genes is not predictable from own expression or immune neighborhood across cancer types (`plans/naive_model.md`). Offset-binomial target form and the control suite are retained for future targets. |
| 2026-09-23 | Data inventory updated: more than 20 Xenium 5K / CosMx 6K datasets available, six processed. Nested data strategy (development sets with cross-validation; reserved discovery sets) proposed. |
| 2026-09-23 | Scope fixed: AIVC is a methodology paper answering questions within the cancer types and platforms in the training data; pan-cancer applicability is presented as potential, not a claim. |
| 2026-09-26 | Perturbation data dropped entirely; no perturbation embedding or module. Inputs are identity and neighborhood (required) and morphology (optional) only. |
| 2026-09-26 | Model reframed as a predictive, interpretable model of cytoplasmic activity: cytoplasmic expression of a data-selected gene set, or pathway-level cytoplasmic activity if per-gene signal is too sparse. The target gene set is not restricted to stress-response genes. Attention-based modality attribution is a required feature; benchmarking against adjacent methods is planned. |
| 2026-09-26 | Application datasets (standard Xenium, MERSCOPE, scRNA-seq, possibly Visium) deferred until the model performs well on Xenium 5K; H&E-only application is out of scope. Six-figure manuscript plan adopted. |
| 2026-09-26 | Next step: exploratory analysis of the processed Xenium 5K datasets on nuclear vs cytoplasmic differential patterns, to produce Fig. 1 and select the target gene set. |
| 2026-09-26 | Training data restricted to FFPE samples for now; fresh-frozen samples deferred because gene localization patterns differ between preparation methods. Nuclear expression added as a candidate identity source (target form C). Fig. 1 plan expanded with basic compartment comparisons (clustering, differential expression) ahead of the target-selection analyses. |
| 2026-09-26 | All six processed Xenium 5K samples verified FFPE from `experiment.xenium`. The naive model is dropped from the project (code and write-up not kept); its retained lessons stay in this tracker. Fig. 1 exploration implemented in `code/2_exploration/` with light steps in a notebook and heavy steps on HGCC (`plans/fig1_exploration.md`). Positive controls replaced by panel lncRNAs (MEG3, MIAT, PVT1, CRNDE, HOTAIR nuclear; NORAD cytoplasmic). Morphology axis for this stage = segmentation-derived features only. The compartment matrices keep all decoded transcripts (no QV filter); recorded as a caveat. |
| 2026-09-27 | Fig. 1 exploration run on all six datasets; results and recommendations in `plans/fig1_exploration.md` section 5. The pre-registered "cytoplasm-only spatially coherent cluster" panel does not exist in these data (OC shows a continuous residual axis only); cytoplasm-specific DE across nuclear subtypes and pathway-level environment associations are the positive findings. Decisions on target form, identity source and gene set remain open for the user. |
| 2026-11-14 | Exploratory stage closed: target-gene selection (634 genes) with its sensitivity analysis, per-sample and pooled subtyping (cytoplasmic and nuclear), and the stress vs substate analysis are complete and audited (`plans/subtype_exploration.md`, `plans/stress_substate_exploration.md`). Their conclusions are recorded in "Exploratory stage" as recommendations. No decision is locked by this entry: target form, identity source, the gene set and the outcome definition (sample-relative log cytoplasmic activity) remain open for the user. |
| 2026-11-14 | The advisor's sketch of 2026-11-12 (`plans/AIVC_ideas_2026-11-12.docx`) proposes nucleocytoplasmic discordance (NCD), the cytoplasmic-nuclear expression difference, as the outcome in place of absolute cytoplasmic expression, with prediction from nuclear expression and application to snRNA-seq. An exploration is planned in `plans/ncd_exploration.md` (design only): NCD as the offset-corrected, sample-centred log cytoplasm-to-nucleus ratio, i.e. the linear predictor of target form B; sample-corrected pooled substates and stress link with counting-noise and abundance controls. Its section 7 records the position that NCD and absolute cytoplasmic expression are the same prediction problem given the input. Nothing is locked by this entry. |
| 2026-11-15 | NCD exploration run and audited (`plans/ncd_exploration.md`). Negative for the 634 target genes and for the 545 NCD-variable genes selected on the full panel: the NCD (binomial residual given the gene's whole-cell count, with a per-cell shift and scale) shows no reproducible substates, and its link to stress is at the counting-noise floor (cytoplasmic substates V 0.15 with the environmental archetype, NCD clusters 0.02, null clusters 0.02). The discussion (its section 7) attributes the cytoplasmic substates to abundance, which is shared by the compartments and removed by NCD; the partition component is small and gene-specific. Consequence for the design (its section 8): NCD stays the estimand that defines what is new, results are to be reported as increments over the constant-partition null M0, and that increment is expected to be small. |
| 2026-11-15 | NCD stage closed after two closing tests on both gene sets (`plans/ncd_exploration.md` section 5.7). (1) Added to the nuclear PCs, the NCD PCs recover 6 to 8 % of the stress information that the cytoplasmic PCs add; nuclear and cytoplasmic PCs add to each other like two reads of abundance. (2) At population level a few dozen genes shift their compartment split with the environment by 5 to 27 %, far above chance, but the shifts concentrate in geometry axes and in genes expressed by neighbouring cells and are coupled to the gene's own abundance as spill-over predicts; no stress program moves. Conclusions for the design (its section 8): the partition alone is not a sufficient basis for the model in this data; NCD remains the estimand that says what is new; every result is to be reported as the increment over the constant-partition null M0; a model of the partition must control for geometry and neighbours' transcripts. Target form, identity source and gene set remain open for the user. |
| 2026-11-18 | The advisor judged the residual model of the NCD exploration wrong: the per-cell scale (kappa) may be biology, and removing it makes the residual uninformative. Two alternatives were named: the plain cytoplasm-to-nucleus ratio, and the model without kappa. The ratio comes first. The first analysis is archived (`plans/ncd_exploration_stale_residual_model.md`, renamed from `plans/ncd_exploration.md` in the rows above; `code/1_exploration/old/ncd_residual/`; `output/1_exploration/old/ncd_residual/`) and its conclusions are not the current position. New plan: `plans/ncd_ratio_exploration.md` (design only): score log2((c+1)/(n+1)) per gene and cell, genes filtered by tumor expression, leakage not considered, target genes selected from the score by variation across cells in both directions. Nothing in the model design is locked by this entry. |
