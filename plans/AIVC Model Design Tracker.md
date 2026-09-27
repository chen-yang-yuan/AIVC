# AIVC Model Design Tracker

Current as of 2026-09-26 (evening).

## Purpose

This tracker records the AIVC model-design brainstorm: the current state of each design component, live candidates, and the decisions locked. It supplements the repo documents `plans/AIVC_plan.md` (background, data, prior model design), `plans/data_processing.md` (pipeline status) and `plans/naive_model.md` (the SG-localization screen), which stay as they are until decisions here are ready to be folded back in.

Conventions: every section except the Decision log describes the latest state only and is rewritten as thinking changes, with stale material removed rather than annotated; the Decision log is append-only and dated. A decision is only locked once it appears in the log.

## The model in one paragraph

AIVC (the only name used in the manuscript) is a predictive, interpretable model of cytoplasmic activity in tumor cells. It takes a cell's identity embedding (own expression, optionally tumor subtype), its neighborhood embedding, and, when H&E is available, its morphology embedding, and predicts the cytoplasmic expression of a selected gene set or, if per-gene signal is too sparse, the cytoplasmic activity of selected pathways. It is trained on subcellular spatial transcriptomics (Xenium 5K, CosMx 6K), which assign reads to nucleus and cytoplasm directly. Its distinguishing features are the compartment-level target and the attribution of each gene's or pathway's cytoplasmic activity to the input modalities. No perturbation data are involved. The target gene set is not restricted to stress-response genes; it will be chosen from the data, with biologically relevant genes (stress response first) preferred among those whose cytoplasmic expression is sensitive to subtype, morphology or neighbors.

## Framing

| Question | Current answer | Status |
| --- | --- | --- |
| Q1. Training data | Xenium 5K and CosMx 6K (expression, compartment-assigned reads, neighbors, H&E), FFPE samples only for now; fresh-frozen samples are deferred because gene localization patterns differ between preparations. More than 20 such datasets are available; the six processed Xenium 5K samples (BC, OC, CC, LC, Prostate, Skin) are all FFPE (verified 2026-09-26 from `experiment.xenium`). | Settled |
| Q2a. Inputs | Identity and neighborhood embeddings required; morphology optional. No perturbation embedding. Whether identity is built from total or nuclear expression is open (see Input design). | Settled except identity source |
| Q2b. Outputs | Cytoplasmic expression of a data-selected gene set; pathway-level cytoplasmic activity as the fallback if per-gene signal is too sparse. Gene set to be determined by the Fig. 1 exploratory analysis. | Form settled; gene set open |
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

Existing assets: scGPT and UNI2-h weights in `data/_utils/`; trial embedding notebooks in `code/3_model/`; naive-model feature pipeline (per-fold SVD, neighborhood geometry cache) in `naive_model/`.

## Output design

Primary target: cytoplasmic expression of a selected gene set. Fallback: cytoplasmic activity of selected pathways, computed by aggregating cytoplasmic counts over each pathway's genes (count aggregation rather than rank-based ssGSEA, which is unstable at per-cell sparsity). Pathway scores can also be derived post hoc from predicted per-gene counts if the per-gene target is retained.

Gene-set selection (to be done in the Fig. 1 analysis): genes whose cytoplasmic expression varies across cells beyond counting noise and is associated with tumor subtype, morphology or neighborhood, filtered for adequate counts, with biologically relevant genes (stress response first) preferred when several qualify.

Target form (open). Three definitions of "cytoplasmic activity" are compatible with the framing and differ in what the model learns and what the identity input may contain:

- A. Absolute cytoplasmic count from total-expression identity, k\_g \~ NB. The target gene's total count must be excluded from the identity embedding, otherwise the prediction is trivial. An imputation task with a compartment-level target and a neighborhood input.
- B. Compartment split given total, k\_g \~ Binomial(n\_g, p\_g) with the per-cell global offset from the naive model. Total counts may stay in the identity embedding. The model learns where a gene's transcripts sit, not how many there are. The naive model showed this signal is small for the 129 SG genes; whether it is larger for a data-selected gene set is what the Fig. 1 analysis will show.
- C. Absolute cytoplasmic count from nuclear-expression identity, k\_g \~ NB with the nuclear counts of all genes as input. No gene exclusion needed; input and target reads are disjoint. The model learns the cytoplasmic pool given the nascent state, and the environmental question becomes whether the neighborhood and morphology change that mapping.

A and B are complementary (predicting both n\_g and p\_g gives the cytoplasmic count as their product with abundance and localization separated); C replaces the gene-exclusion device of A with read-level separation. The Fig. 1 results decide: whether B carries enough signal to stand alone, and whether nuclear expression is a sufficient identity source for C.

Evaluation against the imputation null. Under either form, identity-only prediction is the natural null model. Three fits are reported: M0 identity only, M1 identity + neighborhood, M2 identity + neighborhood + morphology (and neighborhood + morphology without identity as a check). The increments M1 − M0 and M2 − M1 are the evidence that the environment and morphology carry information about cytoplasmic activity beyond own expression, and they are the quantitative counterpart of the attention attributions.

### Evidence from the naive model (retained lessons)

The naive model (a predictability screen of cytoplasmic SG-gene localization from own expression and the immune
neighborhood) is no longer part of the project and its code and write-up are not kept in the repo. What it
established for the SG gene set, retained as design lessons:

- Per-cell cytoplasmic localization of the 129 SG genes is 65–100% counting noise (median 7–29 SG transcripts per cell), reliability 0.00–0.34; beyond that, own expression carries almost no information in five of six samples (LC is the exception and does not transfer), and the immune-only neighborhood carries none.
- An expression-matched non-SG control set is far more predictable than the SG set, so residual skill is generic compartment or segmentation structure, not SG biology.
- Retained: the offset-binomial target form and the evaluation discipline (reliability ceiling, matched-control gene set, leave-one-sample-out and within-sample references, synthetic-signal plumbing check). These are re-implemented in `code/2_exploration/`.

Design lessons: any per-cell compartment target must pool enough transcripts, by choosing well-expressed genes or aggregating to pathways; a gene list from one study is fragile and the target set should be validated on these data; and any new target should be screened with the naive pipeline (swap the target, keep the controls) before architecture work.

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
- Simple statistical baselines: the naive-model GLM on the same features, and a per-gene offset-binomial GLM.

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

## Next step: exploratory analysis for Fig. 1

Goal: show, in an intuitive way, that cytoplasmic expression in tumor cells carries information that nuclear (or total) expression does not, characterize the differential patterns between the two compartments in the six processed Xenium 5K samples (FFPE samples only), and from them select the target gene set. The narrative runs from basic to specific: the compartments are different (Part A), the difference is biology rather than artifact (Part B), and the difference relates to the model's inputs (Part C).

Design rules that apply to every analysis: (i) depth matching, so any nuclear-vs-cytoplasmic comparison downsamples both compartments to the same per-cell depth (rarefaction) or works on fractions with a per-cell leave-one-gene-out offset (measured 2026-09-26: nuclear reads are 37–54% of in-cell reads, more than cytoplasmic in BC and OC, so matching costs little); (ii) a reliability ceiling against binomial counting noise, so that a difference is only reported where it exceeds it; (iii) positive controls for compartment assignment using transcripts with known localization. MALAT1, NEAT1 and XIST are not on the 5K panel; the panel lncRNAs MEG3, MIAT, PVT1, CRNDE and HOTAIR are nuclear-retained in the data and NORAD is cytoplasmic, and they serve as the controls, together with segmentation-method and area-ratio strata checks; (iv) tumor cells only, per sample first, then across samples.

Implementation: `code/2_exploration/` (design, status and results in `plans/fig1_exploration.md`; execution order in the folder README). Morphology in this stage is segmentation-derived (nucleus/cell area ratio, cell area, segmentation method) because H&E exists only for OC; the UNI morphology axis is deferred.

Part A: the two compartments are different.

1. Global description: per sample, the fraction of reads in nucleus vs cytoplasm; per gene, the distribution of cytoplasmic fraction across cells relative to the global offset; the positive controls; and the gene classes that emerge (nuclear-retained, cytoplasm-enriched, balanced).
2. Per-gene compartment concordance: for each gene, the correlation of nuclear and cytoplasmic counts across tumor cells after depth matching. Genes with high abundance but low concordance are where the compartments diverge, and they are the first candidates for the target set.
3. Compartment-specific clustering: cluster tumor cells three times on nuclear, cytoplasmic and total expression with an identical pipeline (same normalization, HVG selection and depth matching). Compare partitions by ARI/NMI; identify clusters that appear only in the cytoplasmic partition; and show whether those clusters are spatially coherent (do they align with tumor boundary, immune-rich regions or morphology clusters). A cytoplasm-only cluster that is spatially coherent is the single most intuitive Fig. 1 panel. This analysis also answers whether nuclear-only clustering recovers the subtypes seen in total expression, which decides the nuclear-identity option.
4. Compartment-specific differential expression: define subtypes from the nuclear partition (transcriptional identity), then run differential expression across those subtypes on nuclear counts and on cytoplasmic counts separately. Compare the top-gene lists: genes differential in the cytoplasm but not the nucleus reflect post-transcriptional regulation (export, stability, localization) rather than transcription, and are the second source of target candidates. Repeat with the cytoplasmic partition as the grouping to show the converse.

Part B: the difference is biology rather than artifact.

5. Matched-control check: for every gene set proposed from Parts A and C, compare against an expression-matched control set drawn from the rest of the panel, so that segmentation, cell-size and depth effects are not mistaken for biology. Genes whose compartment pattern survives this check are kept.
6. Cross-sample reproducibility: which per-gene compartment patterns and which cytoplasm-specific clusters recur across the six samples and cancer types, and which are sample-specific. Recurrent patterns anchor the main figure; sample-specific ones are reported as heterogeneity.

Part C: the difference relates to the model's inputs.

7. Association with the three input axes: for reliable genes, how cytoplasmic fraction and cytoplasmic count vary with tumor subtype (nuclear partition), morphology (UNI embedding clusters) and neighborhood features (cell-type composition, crowding, boundary vs core), each compared with the same association computed on nuclear counts. Output: genes whose cytoplasmic activity is sensitive to each axis beyond what the nucleus shows.
8. Pathway level: analyses 2, 4 and 7 on count-aggregated pathway scores (stress-response pathways first, then the broader GMT), to establish whether pathway targets carry signal where per-gene targets do not.
9. Biological relevance filter: intersect the candidate list with stress-response and other annotated programs; report which programs are represented, and prefer annotated genes when several candidates qualify.

Output of this step: the candidate target gene set (or pathway set), the choice among target forms A, B and C, the identity source (total vs nuclear), and the Fig. 1 panels. Likely panel order: compartment read fractions and positive controls (1); concordance scatter with divergent genes highlighted (2); three-way clustering with a spatial map of the cytoplasm-only cluster (3); DE comparison as paired volcano or overlap plots (4); association with neighborhood or morphology for the top candidates (7); pathway summary (8).

## Ideas under discussion

- Decomposition framing: present the model as partitioning cytoplasmic activity into identity-driven, neighborhood-driven and morphology-driven components, so that a small environmental increment is a quantified finding rather than a failure.
- Cancer type as confounder: hold out by sample; add a pan-cancer regularizer or adversarial cancer-type head.
- Per-sample reporting of compartment signal against each sample's reliability ceiling, so heterogeneity across specimens becomes a finding.
- Stress-invariant identity: if identity and target genes co-express strongly, train the identity encoder to be uninformative of the target residual (adversarial term), or restrict identity genes to lineage and housekeeping sets.

## Open questions

- [ ] Target form A (absolute count, target genes out of identity), B (compartment split given total) or C (absolute count from nuclear identity), or a combination; decided by the Fig. 1 results.
- [ ] Identity source: total or nuclear expression; decided by whether nuclear-only clustering recovers the subtypes (Fig. 1 analysis 3).
- [ ] Target gene set and whether the primary target is per gene or per pathway.
- [x] FFPE status of the six processed samples: all FFPE (2026-09-26). The remaining datasets' status is still to be checked as they are ingested; fresh-frozen samples are excluded from training for now.
- [ ] Tumor-cell definition across datasets: BC/OC/CC use 10x supervised labels, LC/Prostate/Skin manual cluster
  annotation without a non-malignant epithelial class. Step 0b of the Fig. 1 exploration reports the consistency;
  harmonisation (marker QC, CNV inference) to be decided from that report.
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
