# Tumor-cell subtype exploration (`code/1_exploration/`)

Current as of 2026-11-14. Summarises the per-sample, uncorrected pooled and sample-corrected pooled (cytoplasmic and
nuclear) subtyping of the tumor cells of five Xenium 5K samples, answers the two extension questions (are the
per-sample subtypes genuine, and does nuclear expression carry them), and draws the implications for the outcome
variable of the AIVC model. Follows the maintenance rule of `CLAUDE.md`: sections describe the latest state, the Log
is append-only.

## 1. Purpose and scope

The model predicts cytoplasmic activity of tumor cells. Before an outcome can be defined, we need to know whether
the cytoplasmic expression of a cytoplasm-enriched gene set organises tumor cells into substates, whether those
substates recur across tumors of different cancer types, how strong the sample effect is, and whether the nuclear
expression of the same genes carries the same information. Five FFPE Xenium 5K samples are used (BC, OC, CC, LC,
Prostate): the exploration is restricted to the carcinoma samples, and Skin (melanoma, a different tumor lineage) is
excluded from the gene selection onward. Tumor cells are
`cell_type_merged == "Malignant cell"`: BC 102,180, OC 160,250, CC 221,355, LC 44,624, Prostate 95,429 (623,838
pooled).

## 2. Shared procedure

| Step | Where | What |
|---|---|---|
| Gene set | `1_enrichment.ipynb` | 634 genes that are highly variable (seurat_v3) in >= 2 samples and cytoplasm-enriched (log2FC cyto vs nuc > 0, FDR < 0.05) in >= 3 of 5 samples (`output/1_exploration/enrichment/overlap_genes_3plus_samples.npy`). |
| Embedding | `2_embedding.py` (per sample), `4_pooled_embedding.py` (pooled), HGCC | Compartment matrix (all decoded transcripts), `normalize_total(1e4)` on the full panel with the compartment total as size factor, subset to the gene set, log1p, regress out log10(compartment depth + 1), z-score per gene (clip 10), PCA(50), kNN(30 neighbours, 30 PCs), UMAP, t-SNE. |
| Clustering and annotation | `subtype_utils.py`, used by `3_clustering.ipynb` and the three `4_pooled*.ipynb` | Louvain sweep (0.5 to 1.5), Wilcoxon DE, program scores on `layers["lognorm"]` (hallmark/GOBP + curated Interferon, Complement / acute phase, Stress response, Lysosome, Cholesterol synthesis, Myeloid border). A cluster gets a program label when its markers overlap the program (hypergeometric p < 1e-3) or its mean program-score z > 0.5; overrides: Low depth (median < 20 gene-set transcripts), Stroma border (EMT/stroma cluster within 25 um of fibroblasts). Claude reviews every run and fixes the resolution and manual labels in the notebook config with dated reasoning; no cluster is hand-labelled by the user. |
| Pooled comparison | `comparison/` of each pooled run | Label presence pooled vs per sample, per-sample label agreement (share of a per-sample label's cells that carry the same pooled label; recovered = the pooled clustering has the label and the median agreement over the samples that have it is >= 0.5), ARI/NMI, batch mixing (same-sample kNN fraction over its expected value, neighbourhood entropy, sample R2 of PCs), sample-specific clusters (top-sample share / expected share > 3; CC, 35 % of the cells, has a maximum ratio of 2.8 and can never be flagged). Cells without a stored kNN neighbour (identical vectors, i.e. zero counts) are left out of the kNN statistics. |

Depth of the two compartments on the gene set (median transcripts per tumor cell):

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| Cytoplasmic: panel / gene set | 130 / 39 | 200 / 67 | 90 / 30 | 704 / 215 | 140 / 43 |
| Nuclear: panel / gene set | 165 / 35 | 246 / 63 | 74 / 20 | 568 / 138 | 127 / 30 |
| Cells with 0 gene-set transcripts, cyto / nuc | 0.1 % / 4.9 % | 0.0 % / 2.8 % | 0.1 % / 4.4 % | 0.1 % / 3.6 % | 0.2 % / 8.9 % |

LC is a deep sample; the gene set is cytoplasm-enriched by construction, so the nuclear matrix is sparser on it.

### 2.1 Sensitivity of the target gene set to the selection criteria (`5_target_gene_sensitivity.ipynb`)

Requested by the advisor (2026-11-07): how many target genes other search criteria give and how they overlap
with the 634. The notebook runs locally (AIVC-env), reads the five tumor-cell objects and cytoplasmic matrices,
and writes to `output/1_exploration/gene_set_sensitivity/` (`gene_sample_criteria.csv` per gene and sample,
`scenario_counts.csv` per scenario and threshold pair, `hvg_basis_agreement.csv` / `.png`, `detection_vs_logfc.csv`,
`trajectory_counts.png`, `trajectory_overlap.png`, `summary.txt`). The 634-gene set and every analysis built on it
are unchanged. Varied: (1) the enrichment criterion, logFC/FDR (current) or cytoplasmic detection (count > 0 in
>= p % of the sample's tumor cells, p = 5, 10, 20, 30); (2) the HVG basis, whole-cell counts (current) or the
cytoplasmic matrix; (3) the thresholds, HVG in >= N samples (N = 0..5, 0 = no HVG filter) and enriched or detected
in >= M samples (M = 1..5). Fixed: top 2500 HVGs, and the per-sample bottom-10 % total-count drop as the gene
universe of both enrichment criteria. The current rule reproduces the saved 634 genes exactly.

At the current thresholds (N >= 2, M >= 3):

| Enrichment criterion | HVG basis | Genes | Overlap with 634 | Lost | New |
|---|---|---|---|---|---|
| logFC/FDR | whole cell (current) | 634 | 634 | 0 | 0 |
| logFC/FDR | cytoplasm | 675 | 592 | 42 | 83 |
| detect >= 5 % | whole cell | 378 | 269 | 365 | 109 |
| detect >= 5 % | cytoplasm | 429 | 249 | 385 | 180 |
| detect >= 10 % | whole cell | 115 | 96 | 538 | 19 |
| detect >= 10 % | cytoplasm | 130 | 85 | 549 | 45 |
| detect >= 20 % | whole cell | 24 | 24 | 610 | 0 |
| detect >= 30 % | whole cell | 7 | 7 | 627 | 0 |

- **HVG basis**: the two top-2500 lists share 1,949-2,007 genes per sample (Jaccard 0.64-0.67); 592 of
  the 634 are cytoplasmic HVGs in >= 2 samples. The gene set is robust to this choice (Jaccard 0.83 with the 634).
- **Detection criterion**: it is an abundance filter, not a compartment criterion, and the per-cell
  cytoplasmic detection of the 634 genes is low in four of five samples (median 2.7-5.5 %; LC 17 %), so it mostly
  removes genes. At 10 % in >= 3 samples only 115 genes remain (96 of the 634); at 5 % 378 genes. The genes it
  adds are nucleus-enriched by the pseudobulk logFC (101 of 109 new genes at 5 %, 16 of 19 at 10 %), i.e. abundant
  genes whose transcripts sit mainly in the nucleus. The per-sample detection counts of the 634 are dominated by
  LC (at 10 %: BC 89, OC 171, CC 57, LC 434, Prostate 97), so the detection route also tilts the set to the deep
  sample.
- **Thresholds**: with the current criteria the count runs from 2,294 (no HVG filter, M = 1) to 60 (N = 5, M = 5);
  along the current row (N = 2) 1,463 / 910 / 634 / 447 / 235 genes for M = 1..5. The M >= 1 and M >= 2 sets
  contain all 634 by construction; M >= 4 loses 187 of them. Detection scenarios fall by an order of
  magnitude between M = 1 and M = 5 (`trajectory_counts.png`).

Conclusion: the 634 set is stable to the HVG basis and is defined mainly by the logFC criterion and the M >= 3
cut; a detection criterion at 10 % is not a workable substitute on this panel because per-gene detection is sparse
and sample-dependent, and at any threshold it trades compartment specificity for abundance.

## 3. Per-sample analysis (cytoplasmic, `3_clustering.ipynb`)

Resolutions after review: BC 1.0, OC 1.5, CC 1.2, LC 1.5, Prostate 1.5; 13 to 17 clusters per sample. Outputs in
`output/1_exploration/<sample>/` and `cross_sample/`.

Label share (% of tumor cells) per sample:

| Label | BC | OC | CC | LC | Prostate | Samples |
|---|---|---|---|---|---|---|
| S-phase | 9.5 | 13.8 | 10.4 | 6.6 | 2.3 | 5 |
| G2/M | 7.9 | 10.9 | 10.0 | 0 | 0 | 3 |
| Interferon | 5.3 | 3.1 | 6.5 | 12.8 | 0 | 4 |
| Hypoxia | 0 | 14.4 | 13.2 | 20.4 | 5.7 | 4 |
| Stress response | 11.2 | 1.7 | 14.9 | 0 | 15.7 | 4 |
| ER stress | 9.9 | 0 | 0 | 7.8 | 6.5 | 3 |
| Complement / acute phase | 6.2 | 0 | 0 | 12.2 | 4.5 | 3 |
| MYC targets | 0 | 0 | 1.7 | 0 | 18.9 | 2 |
| Lysosome | 2.8 | 0 | 0 | 9.8 | 0 | 2 |
| EMT/stroma | 0 | 0 | 11.0 | 0 | 0 | 1 |
| Cholesterol synthesis | 0 | 0 | 0 | 0 | 1.5 | 1 |
| Baseline | 36.9 | 49.3 | 22.4 | 21.4 | 35.8 | 5 |
| Stroma border | 3.6 | 3.5 | 2.9 | 4.8 | 0 | 4 |
| Myeloid border | 0 | 3.4 | 2.5 | 4.1 | 0 | 3 |
| Low depth | 6.7 | 0 | 4.6 | 0.1 | 9.1 | 4 |

Findings:

- Every sample yields interpretable substates with the same program vocabulary. Cell cycle (S-phase in all five,
  G2/M in the three samples with enough proliferating cells), Interferon, Hypoxia and the two border labels are the
  recurrent ones; their marker lists are near-identical across samples (PCNA/MYBL2/MCM4; TOP2A/BIRC5/TPX2;
  IFIT1/IFIT3/STAT1/CXCL10; NDRG1/LDHA/SLC2A1/P4HA1; COL5A1/POSTN/MMP14; ITGB2/SIRPA/CD68/LIPA).
- The labels mostly have the expected characteristics within each sample. Interferon has the highest interferon
  ssGSEA score in all four samples that have it. Hypoxia has the highest hypoxia ssGSEA score in OC and CC, and is
  second in LC (after ER stress) and third in Prostate (after Stress response and Complement / acute phase; the
  Prostate Hypoxia cluster was a reviewed mTORC1 cluster). The cell-cycle labels are the deepest in BC, OC and CC, not in LC (ER stress) or
  Prostate (Complement / acute phase). Stroma border is the label closest to fibroblasts in all four samples that
  have it (median 15 to 21 um); Myeloid border is the label closest to myeloid cells in LC only. Low depth has a
  median cytoplasmic depth of 23 to 35 transcripts in BC, CC and Prostate; in LC it is 33 cells without any.
- Several labels are sample-specific at this stage: EMT/stroma (CC), Cholesterol synthesis (Prostate), MYC targets
  (Prostate, CC), Lysosome (BC, LC). "Stress response" is a reviewed umbrella label that covers different programs
  in different samples (heat shock HSPB8/HSPA8 in BC, p53/TNFA MALL/HOPX/CDKN1A in CC, TNFA/NFKB ZFAND5/GADD45A/
  NFE2L2 in Prostate, GADD45A/HSPA8/SQSTM1 in a small OC cluster, 1.7 %).
- Baseline (21 to 49 %) holds clusters whose markers are lineage or housekeeping genes (FOXA1, NKX3-1, XBP1, EEF1G);
  they are not substates in the program sense.
- Whether the per-sample labels are the same state in different tumors cannot be decided per sample; this is what
  the pooled runs test.

## 4. Uncorrected pooled analysis (cytoplasmic, `4_pooled.ipynb`, `output/1_exploration/pooled/`)

Depth regression and scaling on the pooled matrix, no correction for sample (the reference for the raw sample
effect). Resolution 1.0, 16 clusters, 9 labels.

- The samples separate completely. Same-sample kNN fraction is 1.00 in every sample (ratio to expectation at its
  maximum: BC 6.1, OC 3.9, CC 2.8, LC 14.0, Prostate 6.5), neighbourhood entropy 0.00, and sample explains 50 to 72 %
  of the variance of PC1 to PC4. UMAP and t-SNE show five islands.
- The islands are driven by lineage genes, the top markers of the one-sample clusters: H19 (OC), XBP1/TSPAN13/
  FOXA1 (BC), SORD/NKX3-1/FOXA1 (Prostate), SLC40A1/PDIA4/EPDR1 (LC). The three largest clusters and the eighth
  are one-sample "Baseline" clusters (OC, BC, Prostate, LC; 98 to 99.6 % one sample). 8 of 16 clusters are flagged
  sample-specific (ratio > 3); three more are 97 to 99 % CC, which the ratio criterion cannot flag.
- No per-sample label is recovered: agreement is above 0.5 only inside CC (Hypoxia 0.84, Stress response 0.70,
  Interferon 0.59, S-phase 0.58) and for OC S-phase (0.73); the BC, LC and Prostate program labels collapse into
  that sample's Baseline island (agreement <= 0.03, except Complement 0.35 in BC and 0.15 in Prostate; the Prostate
  Stress response cells go to the Prostate Stroma border cluster, 0.47). ARI per sample 0.00 to 0.22. G2/M, ER
  stress, Lysosome, Cholesterol synthesis, EMT/stroma, MYC targets have no pooled cluster.
- Conclusion: pooling without a sample correction clusters tumors, not substates. The sample effect is a per-gene
  mean shift concentrated in lineage genes, large relative to the within-sample substate signal.

## 5. Sample-corrected pooled analysis (`--correction sample`)

Depth regression and z-scoring run within each sample (by construction the pooled X is the stacked per-sample X),
then one PCA/kNN/UMAP on all cells. Every gene has mean 0 and unit variance within every sample, so sample-level
mean differences and sample-specific depth slopes are removed before pooling; differences in within-sample
covariance are not.

### 5.1 Cytoplasmic (`4_pooled_sample_corrected.ipynb`, `output/1_exploration/pooled_sample_corrected/`)

Resolution 1.0, 16 clusters, 11 labels (reviewed: cluster 4 G2/M, 3 ER stress, 5 Stress response, 9 Myeloid
border, 10 Low depth, 11 Hypoxia, 0/12 Baseline).

Batch mixing: same-sample ratio BC 1.47, OC 1.21, CC 1.55, LC 1.45, Prostate 2.25 (1 = perfectly mixed);
neighbourhood entropy 0.67 to 0.84; sample R2 of PC1 to PC10 = 0.00. UMAP and t-SNE coloured by sample show one
cloud. 9 of 16 clusters, 92 % of the cells, contain every sample at close to its expected share (top-sample ratio
1.0 to 1.9); the exceptions are one LC-dominated Hypoxia cluster (7,001 cells, 82 % LC), a 355-cell LC Baseline
cluster, three tiny Low depth clusters (33 to 91 cells) and two Prostate-enriched clusters (Complement / acute
phase 45 % Prostate, ratio 3.0; Low depth 38 %, ratio 2.5). The Stress response cluster is 68 % CC (ratio 1.9).

Agreement of the per-sample labels with the pooled labels (share of cells; recovered = median over the samples
that have the label >= 0.5):

| Per-sample label | BC | OC | CC | LC | Prostate | Verdict |
|---|---|---|---|---|---|---|
| S-phase | 0.84 | 0.74 | 0.76 | 0.47 | 0.42 | recovered (below 0.5 in LC/Prostate, split with G2/M and Baseline) |
| G2/M | 0.86 | 0.92 | 0.82 | . | . | recovered |
| Interferon | 0.72 | 0.86 | 0.60 | 0.48 | . | recovered (LC just below 0.5) |
| Hypoxia | . | 0.79 | 0.79 | 0.67 | 0.24 | recovered (Prostate's Hypoxia cluster was a reviewed mTORC1 cluster) |
| Stroma border | 0.73 | 0.84 | 0.74 | 0.90 | . | recovered |
| Myeloid border | . | 0.55 | 0.71 | 0.68 | . | recovered |
| ER stress | 0.46 | . | . | 0.36 | 0.66 | partial (median 0.46) |
| Complement / acute phase | 0.32 | . | . | 0.03 | 0.48 | partial (LC's cluster was a medium-confidence Apoptosis overlap) |
| Stress response | 0.15 | 0.06 | 0.80 | . | 0.06 | CC only; the BC/OC/Prostate clusters scatter (different programs under one label) |
| Lysosome | 0.00 | . | . | 0.00 | . | not reproduced (BC -> Myeloid border 0.47, LC -> ER stress 0.55) |
| EMT/stroma | . | . | 0.00 | . | . | not reproduced (-> Stroma border 0.49) |
| MYC targets | . | . | 0.00 | . | 0.00 | not reproduced (-> S-phase / Baseline) |
| Cholesterol synthesis | . | . | . | . | 0.00 | not reproduced (-> ER stress 0.35) |

ARI of per-sample vs pooled labels within a sample: BC 0.25, OC 0.27, CC 0.36, LC 0.22, Prostate 0.15 (NMI 0.19 to
0.40), against 0.00 to 0.22 uncorrected. The pooled labels keep their biology: the pooled Hypoxia label has the
highest hypoxia ssGSEA score (0.038 vs -0.007 Baseline), Interferon the highest interferon ssGSEA score (0.073 vs
0.029), G2/M the highest depth (median 166; S-phase, 139, is at the level of Baseline), Stress response the
largest cells (97 um2, nucleus ratio 0.44) and the sparsest neighbourhoods, Stroma border the lowest nucleus ratio
after it. Per-sample label shares are similar
across samples (S-phase 7 to 17 %, Hypoxia 10 to 22 %, ER stress 8 to 21 %, Baseline 12 to 27 %) because the
correction removes between-sample differences in average activity by construction.

### 5.2 Nuclear (`4_pooled_sample_corrected_nuclear.ipynb`, `output/1_exploration/pooled_sample_corrected_nuclear/`)

The same procedure on the nuclear matrix of the same 634 genes, compared cell by cell with the cytoplasmic
sample-corrected labels of 5.1 (section 11 of the notebook, `comparison/*_vs_pooled_cyto*`). Resolution 1.0, 21
clusters, 12 labels (reviewed: 0/14 Baseline, 4 G2/M, 11/19 Stress response, 15 Myeloid border, 16 Hypoxia).

- Samples mix about as in the cytoplasmic run (same-sample ratio BC 1.47, OC 1.07, CC 1.75, LC 1.93, Prostate
  1.80; entropy 0.63 to 0.81; sample R2 = 0; the 28,694 cells without kNN neighbours are excluded); the
  same LC-only Hypoxia cluster appears (5,260 cells, 83 % LC).
- The nuclear clustering finds the same program clusters at the cluster level: S-phase (9.6 %), G2/M (9.6 %),
  Interferon (8.3 %), Hypoxia (11.6 %), ER stress (4.7 %), Complement (2.0 %), with largely the same top markers
  as the cytoplasmic clusters (the nuclear Interferon cluster also carries MALL/CDKN1A/SAMD9/HOPX, markers of the
  cytoplasmic Stress response, and is the best nuclear match of that label, 0.26). It adds two labels the cytoplasmic pool does not have, Cholesterol synthesis (4.5 %, DHCR24/
  SQLE/HMGCR/FDFT1/EBP, 2.8 to 5.8 % of every sample) and EMT/stroma (8.0 %, CDKN1A/COL5A1/MMP14/ITGA6, not at the
  fibroblast border), and loses Stress response (2.0 % vs 8.5 %) and Stroma border (absent). 29,197 cells (4.7 %)
  have zero nuclear gene-set transcripts (27,449 of them no nuclear transcript at all); they form one zero-count
  Low depth cluster per sample (27,495 cells); nuclear Low depth is 8.5 % of cells against 2.8 % cytoplasmic.
- Cell by cell, the two partitions do not match. ARI 0.05, NMI 0.07 over all cells (0.03 to 0.11 per sample).
  Agreement of each cytoplasmic label with the same nuclear label: Baseline 0.45, G2/M 0.43, Hypoxia 0.34, S-phase
  0.27, Interferon 0.27, Low depth 0.13, Complement 0.09, ER stress 0.08, Myeloid border 0.03, Stress response 0.02,
  Stroma border 0.00. Within single samples only OC G2/M (0.57) and Baseline in OC and LC (0.56, 0.54) pass 0.5;
  the next best are LC Hypoxia 0.48 and OC Hypoxia 0.47.
- Per-cell program scores correlate weakly between compartments: median Spearman 0.23 over 23 programs; S-phase
  0.52, Estrogen late 0.43, Hypoxia 0.43, Estrogen early 0.41, G2/M 0.37, ER stress 0.33; Stress response 0.14,
  Cholesterol synthesis 0.09, Autophagy 0.09. Within a single sample the correlations are lower (S-phase 0.23 to
  0.49, Hypoxia 0.08 to 0.33).
- Neighbourhood purity of the cytoplasmic labels on the nuclear kNN graph is a third of the purity on the
  cytoplasmic graph (median ratio 0.34; S-phase 0.22 vs 0.67, G2/M 0.25 vs 0.71, Hypoxia 0.24 vs 0.67, Interferon
  0.14 vs 0.64, Myeloid border 0.04 vs 0.55; cells without neighbours in a graph are excluded from that graph).
  Against random expectation the nuclear graph still carries some
  structure for cell cycle and Hypoxia (0.22 to 0.25 vs 0.09 to 0.14), nearly none for the border labels.

## 6. Comparison and answers to the two extension questions

### (a) Are the per-sample subtypes genuine, and does the sample effect diminish?

Genuine, for the recurrent program substates. S-phase, G2/M, Interferon and Hypoxia, plus the two border labels,
reappear in the sample-corrected pool as clusters that contain every sample at its expected share, carry the same
markers, keep the same characteristics (ssGSEA scores, depth, proximity), and absorb 60 to 92 % of the cells that
carried the same label per sample in BC, OC and CC (OC Myeloid border 55 %), and 42 to 90 % in LC and Prostate,
whose cell-cycle clusters were smaller (Prostate Hypoxia, a reviewed mTORC1 cluster, 24 %). ER stress and
Complement / acute phase exist in the pool (12 % and 4 %) but match the per-sample clusters only partly (0.32 to
0.66; LC Complement 0.03), so their per-sample boundaries are less stable. Four
per-sample labels are not reproduced (Lysosome, EMT/stroma, MYC targets, Cholesterol synthesis): their cells join
neighbouring pooled states, so as cytoplasmic clusters they are sample-specific partitions rather than shared
substates. "Stress response" is a shared label but not a shared state: only the CC program (MALL/HOPX/CDKN1A)
forms a pooled cluster, and it stays 68 % CC.

The sample effect diminishes with the per-sample correction, and the UMAP/t-SNE mixing is real but partly
guaranteed by construction. Uncorrected: same-sample fraction 1.00, entropy 0, sample R2 of PC1 to PC4 0.50 to 0.72,
five islands. Corrected: same-sample ratio 1.2 to 1.6 (Prostate 2.25), entropy 0.67 to 0.84, sample R2 0.00, one
cloud. Because the correction forces every gene to mean 0 and variance 1 within each sample, the removal of the
sample mean is not a finding; the findings are (i) that a per-gene, per-sample mean shift is sufficient (no
residual sample axis in the PCs, 9 of 16 clusters with 92 % of the cells mixed), and (ii) what remains after it:
Prostate cells still sit 2.25 times more often next to Prostate cells, LC keeps a private Hypoxia cluster (deep
sample, larger within-sample variance), and the Stress response cluster is CC-dominated. These residuals are
within-sample covariance differences that a mean/variance correction cannot remove.

Caveat: the correction cannot distinguish a batch shift from a genuine difference in average program activity
between tumors (a tumor that is hypoxic throughout loses that information). Pooled label frequencies are therefore
relative within each sample.

### (b) Does nuclear expression reveal the same subtype pattern?

At the level of programs, partly: the nuclear pool has S-phase, G2/M, Interferon, Hypoxia, ER stress and Complement
clusters with the same markers, so these programs are transcribed and detectable in the nucleus. At the level of
cells, no: the nuclear and cytoplasmic assignments of the same cell agree at chance level (ARI 0.05), no cytoplasmic
program label is recovered by its nuclear counterpart (agreement over all cells 0.27 to 0.43 for S-phase, G2/M,
Interferon and Hypoxia, 0 to 0.09 for the other program and border labels; only OC G2/M passes 0.5 within a
sample), per-cell program scores correlate 0.1 to 0.5 between compartments, and cytoplasmic labels are three
times less pure in the nuclear kNN graph. The labels that vanish in the nucleus are those that are cytoplasmic in
origin or in mechanism: Stroma border and Myeloid border (transcripts of neighbouring cells captured by the
expansion-based cytoplasm segmentation), Stress response (mature-mRNA programs), and ER stress largely. The nucleus
instead resolves a Cholesterol synthesis state (SREBP targets) in every sample and an EMT/stroma state away from the
border, which the cytoplasm of the same genes does not.

Two reasons combine, and the data do not separate them: (i) depth, since the gene set is cytoplasm-enriched by
selection, the nuclear matrix has 20 to 138 gene-set transcripts per cell (cytoplasm 30 to 215) and 4.7 % of the
cells have none, so per-cell nuclear assignments are noisier; (ii) biology, since nuclear reads approximate nascent
transcription and cytoplasmic reads the mature pool, and the two differ by export, stability and localisation
(the Fig. 1 exploration found disattenuated nuclear-cytoplasmic correlations of 0.46 to 0.66 per gene at tile
level). The answer is specific to this cytoplasm-enriched gene set; on the full panel the Fig. 1 exploration found
nuclear-only clustering to recover the total partition as well as cytoplasmic clustering does.

## 7. Implications for the outcome variable of the model

What the exploration establishes about a cytoplasmic outcome:

1. The within-sample transform used here (compartment size factor on the full panel, log1p, regress out log
   cytoplasmic depth within sample, z-score per gene within sample) produces a quantity whose structure is shared
   across five tumors of five cancer types, with program states that keep their biology and their characteristics.
   The same quantity without the per-sample step is dominated by the tumor of origin. A sample-relative definition
   of cytoplasmic activity is therefore the one that can be learned from several samples at once.
2. The sample effect is a per-gene mean shift (and a per-sample depth slope), concentrated in lineage genes. It is
   removed by per-gene, per-sample intercepts; nothing more elaborate (embedding-level batch integration) was
   needed to mix the samples.
3. The reproducible signal is at program level. With a median of 30 to 70 gene-set transcripts per cell (LC 215),
   per-gene per-cell values are mostly zeros and the substates only emerge from the covariance of many genes; the
   per-sample and pooled clusters are defined by programs of 5 to 40 genes. This agrees with the Fig. 1
   recommendation of a pathway-level primary target with a per-gene track.
4. The nuclear expression of these genes is not a proxy for their cytoplasmic state cell by cell, so a nuclear
   identity input (target form C of the tracker) is non-redundant with the cytoplasmic outcome: the model would
   have to predict the cytoplasmic substate from nuclear identity plus neighbourhood and morphology, which is the
   intended question. The identity-only null will be weak for this gene set, so the environmental increments need
   a noise ceiling (split-half reliability of the cytoplasmic outcome) to be interpretable, and nuclear depth
   limits BC, CC and Prostate.

### 7.1 Observed quantities

- Unit: tumor cell c of sample s, gene g of the target set, cytoplasmic count k_cg from the compartment matrix
  (`overlaps_nucleus == 0`), depth d_c = total cytoplasmic transcripts of the full panel (not of the target set,
  and not `total_counts`, which includes control probes). Cells with fewer than about 20 target-set cytoplasmic
  transcripts (0 to 9 % per sample, the Low depth clusters) are excluded from the outcome or down-weighted.
- Program P: K_cP = sum over g in P of k_cg, the count-aggregated cytoplasmic activity.

### 7.2 What the predictive model predicts

The model (identity, neighbourhood and optional morphology embeddings fused by attention, per the tracker) outputs
for every cell c one real number per target gene, eta_cg = f_g(x_c), and one per program, eta_cP. The number is the
**sample-relative log cytoplasmic activity**: the log fold change of the cytoplasmic expression of g in cell c
relative to the average tumor cell of the same sample at the same cytoplasmic depth. eta_cg = 0 means a typical
cell of that tumor, +0.7 means twice the typical cytoplasmic level, -0.7 half. It is dimensionless, free of depth
and free of the tumor's baseline level of g. Neither the count nor the sample baseline is predicted by the network;
they enter through the link below.

### 7.3 Link to the actual expression of gene g in cell c

Count link (primary, negative-binomial likelihood, used for training and for likelihood-based evaluation):

    k_cg ~ NB(mu_cg, theta_g),    log mu_cg = log d_c + a_gs + eta_cg

- log d_c is a fixed offset: the model is not asked to predict how many transcripts a cell has.
- a_gs is a per-gene, per-sample intercept: a learned parameter table (samples x genes), fitted jointly with the
  network and, for a new sample, re-fitted alone with the network frozen from that sample's tumor cells (a mean per
  gene; a few hundred cells suffice). It is the tumor's baseline level of g (lineage expression, panel run, FFPE
  chemistry), the quantity that the per-sample correction of section 5 removed.
- theta_g is a per-gene dispersion.
- The identifiability constraint eta has mean zero within each sample for each gene (a penalty during training, or
  recentering at the end) fixes the split between a_gs and eta, so that a_gs is the baseline and eta the deviation.

Reading the prediction back on the expression scale:

- expected cytoplasmic count:          mu_cg = d_c * exp(a_gs) * exp(eta_cg)
- normalised expression (counts per 10k, the `lognorm` scale before log1p):  1e4 * exp(a_gs + eta_cg)
- fold change against the sample's typical cell:  exp(eta_cg)
- program activity: the same three lines with K_cP, a_Ps and eta_cP.

Continuous link (equivalent Gaussian head, the form validated by the clustering of sections 3 to 5): the model
predicts r_cg, the within-sample z-score of [log1p(1e4 k_cg / d_c) regressed on log10(d_c + 1)], which is exactly
X of the sample-corrected embedded objects. The link back is
log1p(1e4 k_cg / d_c) = m_gs + b_gs * log10(d_c + 1) + sd_gs * r_cg, with the sample parameters (m_gs, b_gs, sd_gs)
stored next to the model. The two heads carry the same information; the count head is preferred because it keeps
the counting noise in the likelihood and does not need the variance scaling, which also removes genuine
differences in within-tumor heterogeneity (LC).

Program-level readout (recommended primary): eta_cP for S-phase, G2/M, Interferon, Hypoxia, ER stress and
Complement / acute phase, with Cholesterol synthesis and EMT/stroma as nuclear-informed candidates. Stress response
is a candidate only after it is split into its component programs (section 8); as labelled now, only its CC-type
program (MALL/HOPX/CDKN1A) forms a pooled cluster. Per-gene track: the recurrent per-gene candidates of the Fig. 1
exploration that are in the target set (NDRG1, LDHA, EEF1G, ...). The Fig. 1 list (43 genes, built on six samples
with other criteria, `plans/fig1_exploration.md` section 5.3) shares only 16 genes with the 634; a candidate
outside the set (e.g. YWHAZ) would have to be added to the target genes explicitly. Program scores can also be derived post hoc from the per-gene predictions through the count link.

### 7.4 What the per-sample intercept does and does not cost

The intercept is a term of the model, not the outcome: the fitted model does know that tumor A expresses g above
tumor B (a_gA > a_gB), it only does not explain that difference from the inputs. Between-sample level
differences are conditioned on, not modelled. Three consequences:

- Not every between-sample difference goes into the intercept. With a shared f_g, what the inputs explain stays in
  eta: a sample with 30 % hypoxic cells against 10 % keeps its higher mean of a hypoxia gene attributed to those
  cells. The caveat is identifiability: the identity embedding identifies the sample through lineage genes, so a
  flexible f_g and the intercept can trade a level shift without changing the likelihood. The honest statement is
  therefore: within-sample contrasts are interpretable, between-sample level differences are not.
- The information that is given up was not recoverable. With five samples, one per cancer type, in two panel runs,
  a between-sample level difference is confounded with cancer type, tissue, FFPE chemistry and run. The
  uncorrected pooled run (section 4) shows what a model without the intercept learns: a tumor classifier on H19,
  NKX3-1, SORD and XBP1, which would contaminate the environment and morphology attributions, the model's claimed
  findings. The intercept is the same device as a batch covariate in differential expression.
- What stays answerable across samples: whether a neighbourhood or morphology feature shifts eta_cP by the same
  amount in BC and in CC (f_g is shared), the counterfactual readouts (swap the neighbourhood, read the change in
  eta; within-sample by nature), and a descriptive comparison of the baselines a_gs with cancer type and panel run.
  Not answerable: why tumor A is more hypoxic overall than tumor B; that is a cohort question that needs tens of
  samples with covariates, outside the within-tumor scope fixed in the tracker.

Options kept open: (i) a per-gene test of whether the intercept is needed, fitting each gene with a shared and with
per-sample intercepts and comparing held-out likelihood; the uncorrected run suggests the sample effect sits in 10
to 20 lineage genes, so most target genes may keep an absolute interpretation; (ii) random intercepts with
sample-level covariates (cancer type, panel run) once the 20-plus available datasets are processed, which turns
the fixed baselines into a decomposed between-sample variance; (iii) a per-sample depth slope b_s next to the
offset if the slopes differ between samples (the exploration allowed them to).

### 7.5 Open checks before locking

(i) split-half reliability of r_cg and K_cP per sample, to set the ceiling for the identity-only null and the
environmental increments; (ii) the per-gene intercept test above; (iii) leakage: the border labels show that
cytoplasmic counts include neighbours' transcripts, so a neighbourhood input predicts them trivially for stroma-
and myeloid-adjacent cells; the matched-control discipline of the Fig. 1 exploration applies.

## 8. Open items

- Decide on the outcome form (section 7) and record it in the tracker's decision log; the recommendation is
  listed there as not locked. The environment-to-outcome link that the recommendation relies on is tested in
  `plans/stress_substate_exploration.md` (stress axes vs the cytoplasmic and nuclear substates; cytoplasmic
  expression carries more stress information on 19 of 20 environmental axes and in both label-based approaches,
  pending the thinned-cytoplasm control).
- The Stress response vocabulary lumps different programs (heat shock, p53/TNFA, NFKB); split it before it is used
  as a target.
- Residual sample effects after the correction (Prostate, LC Hypoxia cluster, panel-run grouping) are not
  removed by mean/variance correction; evaluate by leave-one-sample-out rather than by further correction.
- Skin (melanoma) is excluded throughout; the gene set and the conclusions cover five carcinoma samples. Whether
  the 634 set and the substates hold in melanoma is untested.

## 9. Code, outputs and re-run rules

Run order in `code/1_exploration/` (notebooks locally in `AIVC-env` from that folder; `.py` + `.sh` on HGCC, then
`make out-pull`): `1_enrichment.ipynb` -> `2_embedding.py` -> `3_clustering.ipynb` -> `4_pooled_embedding.py` (three
variants: `--correction none`, `--correction sample`, `--correction sample --compartment nuclear`) -> `4_pooled.ipynb`,
`4_pooled_sample_corrected.ipynb`, `4_pooled_sample_corrected_nuclear.ipynb` (in this order: the nuclear notebook
reads the tables of the sample-corrected cytoplasmic run, and all three read the per-sample tables of
`3_clustering.ipynb`) -> `5_target_gene_sensitivity.ipynb` (independent) -> `6_stress_substates.ipynb` (reads the two
sample-corrected pooled objects). `old/` is the archive of the superseded first design (full panel, six samples,
parameter grid) and is not part of the pipeline.

| Folder in `output/1_exploration/` | Content | Written by | On re-run |
|---|---|---|---|
| `enrichment/` | `overlap_genes_3plus_samples.npy` (the 634 genes), `overlap_genes_counts.png` | `1_enrichment.ipynb` | overwritten |
| `<sample>/` (root) | `enrichment_results.csv` (per-gene compartment statistics; the Wilcoxon column is spelled `wilconxon_stat`), `adata_cytoplasmic_embedded.h5ad`, `adata_nuclear_embedded.h5ad` | `1_enrichment.ipynb`, `2_embedding.py` | kept by the clustering notebook. The per-sample nuclear objects are read by no current notebook (the nuclear analysis is pooled). |
| `<sample>/{clustering,DE,annotation,spatial,characteristics}/` | sweep, cluster statistics, DE tables, annotation, per-cell labels and characteristics, figures | `3_clustering.ipynb` | emptied and rewritten |
| `cross_sample/` | label composition and characteristics across samples, all annotations, all per-cell labels (summaries; not read downstream) | `3_clustering.ipynb` | overwritten, only when all five samples were run |
| `pooled/`, `pooled_sample_corrected/`, `pooled_sample_corrected_nuclear/` | the pooled object (rewritten by the notebook with labels, scores and `uns["subtyping"]`: date, resolution, manual labels, notebook, compartment, correction) and the six subfolders (the five above plus `comparison/`, with `summary.txt`) | `4_pooled_embedding.py`, the three pooled notebooks | subfolders emptied and rewritten; object rewritten |
| `gene_set_sensitivity/` | section 2.1 | `5_target_gene_sensitivity.ipynb` | overwritten |
| `stress/` | `plans/stress_substate_exploration.md` | `6_stress_substates.ipynb` | emptied and rewritten (cached stress table kept) |

Louvain and the embeddings are deterministic here: re-running a notebook on the same object reproduces the cluster
ids and labels exactly (checked on every re-run, last 2026-11-14), so the reviewed configs stay valid.

## Log

- 2026-10-05/06: per-sample clustering reviewed (resolutions BC 1.0, OC 1.5, CC 1.2, LC 1.5, Prostate 1.5).
- 2026-11-03: analysis functions moved to `subtype_utils.py` and the per-sample notebook re-run with the reviewed
  resolutions. Uncorrected pooled run (complete sample separation, lineage genes) and sample-corrected pooled run
  (S-phase, G2/M, Interferon, Hypoxia recovered; samples mix) reviewed. Nuclear counterpart run and reviewed:
  nuclear expression does not recover the cytoplasmic substates cell by cell (ARI 0.05).
- 2026-11-04: `make out-pull` made safe (`--update`), all three pooled notebooks re-executed from the pulled
  embeddings; results identical to the reviewed runs. This document written as the summary of the stage and the
  basis for the outcome-variable decision. Section 7 extended with the model's prediction target (sample-relative
  log cytoplasmic activity eta), the count and continuous links to the observed expression, and the interpretive
  scope of the per-gene, per-sample intercept.
- 2026-11-07: target-gene-set sensitivity analysis added (`5_target_gene_sensitivity.ipynb`, section 2.1):
  the 634 set is robust to the HVG basis (Jaccard 0.83 with cytoplasmic HVGs) and the detection criterion is not a
  substitute for the logFC criterion (10 % in >= 3 samples leaves 115 genes, added genes are nucleus-enriched).
  Gene set and downstream analyses unchanged.
- 2026-11-08: pointer to the stress vs substate analysis (`plans/stress_substate_exploration.md`) added to §8; it
  uses the cytoplasmic and nuclear pooled labellings of this stage unchanged.
- 2026-11-14 (final audit of the exploratory stage): two defects of the pooled notebooks fixed. (1) Cells without
  kNN neighbours (zero counts: 282 uncorrected, 132 sample-corrected, 28,694 nuclear) were counted as
  same-sample fraction / entropy / purity 0; they are now excluded. (2) `uns["subtyping"]` records each notebook's
  own name, compartment and correction. The three pooled notebooks were re-run: every cluster, label, agreement
  and characteristics table is byte-identical to the previous run; only the nuclear batch-mixing and purity values
  changed (section 5.2; purity ratio 0.32 -> 0.34, same-sample ratio of CC / LC / Prostate 1.65 / 1.45 / 1.58 ->
  1.75 / 1.93 / 1.80) and one entropy value of the sample-corrected run in the third decimal. Document checked
  against the output tables and corrected (rule definitions, rounding, cluster and zero-count numbers, claims
  restricted to where the tables support them); section 9 added. Skin is excluded as a melanoma. Stale comments
  and unused code removed from the scripts; `code/1_exploration/old/` marked as archive. No analysis choice and no
  conclusion changed.
