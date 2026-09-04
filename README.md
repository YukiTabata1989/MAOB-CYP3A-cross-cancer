# Recurrent tumor-associated loss of MAOB transcript expression and member-dependent CYP3A changes across 14 human cancer types

This repository contains the analysis code, concordance-audit code, processed result tables, figure-generation code, and reproducibility information associated with the manuscript:

**"Recurrent tumor-associated loss of MAOB transcript expression and member-dependent CYP3A changes across 14 human cancer types."**

The study examines the normal-tissue expression architecture of **MAOB, CYP3A4, CYP3A5, CYP3A7, and CYP3A43** and their paired tumor–adjacent non-tumor transcript differences across 14 TCGA cancer types. The primary analysis uses count-level TCGA RNA-sequencing data and GTEx v8 tissue-median RNA expression. Human Protein Atlas RNA data are used only for a descriptive comparison of normal-tissue rankings.

> **Scope note:** This repository contains transcriptomic analyses only. The study does not measure GGA abundance, protein abundance, enzyme activity, or metabolic flux, and does not establish temporal or causal relationships.

## Repository structure

```text
MAOB_CYP3A4_crosscancer/
├── README.md
├── LICENSE
├── CITATION.cff
├── Supplementary_Tables_S1_S13.xlsx
├── code/
│   ├── 01_round2B_primary_analysis.R
│   ├── 02_round2C_xena_audit.R
│   ├── 03_round2D_figure_reconstruction.R
│   ├── 04_round2E_figure_corrections.R
│   ├── 05_round2F_publication_layout.R
│   ├── 06_round2G_final_figure_corrections.R
│   ├── 07_round2H_supplementary_figures.R
│   ├── 08_round3A_denominator_aligned_recurrence.R
│   ├── 09_round3B_figure_rebuild.R
│   ├── 10_round3C_figure_annotation_fixes.R
│   └── 11_round3D_coexpression_asterisks.R
├── docs/
│   ├── Supplementary_Methods.docx
│   └── Figure_legends.docx
├── figures/
│   ├── Figure1_round2G_GTEx_v8_normal_tissue_atlas.pdf
│   ├── Figure2_round2G_count_level_paired_heatmap.pdf
│   ├── Figure3A_round2F_tissue_identity_distribution.pdf
│   ├── Figure3B_round2F_matched_background_context.pdf
│   ├── Figure4A_round3C_denominator_aligned_recurrence_null.pdf
│   ├── Figure4B_round3C_coverage_sensitivity.pdf
│   ├── Figure4C_round3C_MAOB_organ_collapsed_sensitivity.pdf
│   ├── SupplementaryFigureS1A_round3D_nonLIHC_coexpression_change_heatmap.pdf
│   ├── SupplementaryFigureS1B_round2G_representative_coexpression_rank_scatter.pdf
│   ├── SupplementaryFigureS2_round2H_normal_baseline_paired_difference.pdf
│   ├── SupplementaryFigureS3_round2H_HPA_only_RNA_expression.pdf
│   └── SupplementaryFigureS4_round3B_pseudo_target_calibration_QQ.pdf
├── logs/
│   └── sessioninfo/
│       ├── round2B_sessionInfo.txt
│       ├── round2C_sessionInfo.txt
│       ├── round2F_sessionInfo.txt
│       ├── round2G_sessionInfo.txt
│       ├── round2H_sessionInfo.txt
│       ├── round3A_sessionInfo.txt
│       ├── round3B_sessionInfo.txt
│       ├── round3C_sessionInfo.txt
│       └── round3D_sessionInfo.txt
└── results/
    ├── round2B/
    │   ├── logs/
    │   └── tables/
    ├── round2C/
    │   ├── logs/
    │   └── tables/
    └── round3A/
        ├── logs/
        └── tables/
```

## Analysis overview

### Round 2B — primary count-level analysis

`01_round2B_primary_analysis.R`

This is the source of every effect estimate reported in the manuscript. It:

- retrieves and analyzes harmonized TCGA STAR count data;
- fits paired tumor-versus-adjacent-non-tumor models using `limma-voom`;
- includes LIHC MAOB and CYP3A4 in the paired pan-cancer analysis;
- excludes the entire LIHC project from the MAOB–CYP3A coexpression module;
- constructs GTEx-based tissue-identity reference sets;
- performs expression-property-matched background comparisons;
- performs the same-gene cross-cancer recurrence analysis;
- evaluates the exploratory normal-baseline/tumor-difference relationship;
- performs non-LIHC MAOB–CYP3A coexpression analyses;
- performs the HPA-only descriptive normal-tissue comparison.

The recurrence, matched-background and coexpression inference produced by this script was subsequently superseded by Round 3A; see the version note below.

### Round 2C — transformed-expression concordance audit

`02_round2C_xena_audit.R`

A sensitivity audit using previously cached UCSC Xena Toil transformed-expression data. It does **not** replace the count-level Round 2B analysis.

Round 2C re-ran the audit alone, so that the LIHC MAOB and CYP3A4 cells, which had been masked under an earlier rule, are included. The audit results reported in the manuscript, covering 43 overlapping estimable cells, are those in `results/round2C/`. Round 2C verified by checksum that no Round 2B output was altered.

### Round 2D to 2H — manuscript figures

`03_round2D_figure_reconstruction.R` through `07_round2H_supplementary_figures.R`

These stages reconstruct, correct and lay out the manuscript figures from frozen analysis outputs without refitting any statistical model. Round 2D reconstructs the figures, Round 2E applies targeted corrections, Round 2F renders at the publication width of 174 mm with a common minimum text size, Round 2G finalises selected figures including continuous colour-bar orientation, and Round 2H finalises Supplementary Figures S2 and S3.

### Round 3A — denominator-aligned recurrence, calibration, and revised inference

`08_round3A_denominator_aligned_recurrence.R`

This stage re-executed three inferential components without refitting any differential-expression model.

**Denominator alignment.** In the Round 2B recurrence null, a target could be scored on a complete cross-cancer profile while most eligible background genes were scored on incomplete ones, so the empirical tail probabilities were not on the same footing. Round 3A fixes a common set of cancers for each target and coverage setting and re-scores the target and every background gene on that set.

**Calibration.** The aligned procedure was calibrated by treating each eligible background gene in turn as a pseudo-target and comparing it with the remainder under the identical procedure. Tie-randomized rank values are the primary diagnostic; the discrete +1 values are exported alongside.

**Two-sided matched-background inference.** The matched-background comparison now uses two-sided empirical probabilities, since no directional hypothesis was specified in advance for any cell. The one-sided values are retained in the output.

**Permutation-based coexpression inference.** The coexpression change is now tested by a paired label permutation test, in which the adjacent-non-tumor and tumor labels are swapped jointly within each patient, rather than by the frequency with which a bootstrap distribution crosses zero.

**Organ-collapsed sensitivity.** The recurrence analysis was repeated with the three renal and two lung histologies collapsed to organ units.

Round 3A verified by checksum that the Round 2B differential-expression output was unchanged, so that no effect estimate was recomputed at any stage.

### Round 3B to 3D — figures rebuilt from the Round 3A output

`09_round3B_figure_rebuild.R`, `10_round3C_figure_annotation_fixes.R`, `11_round3D_coexpression_asterisks.R`

Round 3B rebuilds Figure 4 from the denominator-aligned results, adds Figure 4C for the organ-collapsed sensitivity analysis, and produces Supplementary Figure S4 for the pseudo-target calibration. The calibration figure is drawn from the frozen values in `166_round3A_pseudo_target_calibration_values.csv` rather than copied from an earlier render, so that the final-size text audit and font embedding apply to it as to every other figure.

Round 3C corrects annotation placement in Figure 4A, 4B and 4C without altering any plotted value.

Round 3D regenerates Supplementary Figure S1A so that the asterisks mark the eight cells identified by the paired-label permutation test rather than the seven identified by the earlier bootstrap criterion.

None of these stages recomputes a statistic; each reads frozen Round 3A output and applies graphical transformations only.

## Recommended execution order

```text
1. Round 2B
2. Round 2C
3. Round 3A
4. Round 2D
5. Round 2E
6. Round 2F
7. Round 2G
8. Round 2H
9. Round 3B
10. Round 3C
11. Round 3D
```

Round 3A depends only on the Round 2B and Round 2C outputs. Rounds 2D to 2H produced the figures that draw on the Round 2B results, and Rounds 3B to 3D produced or replaced those that draw on the Round 3A results.

## Software environment

The final analysis and figure scripts were developed under:

- **R 4.6.1**

Major packages include:

- TCGAbiolinks
- SummarizedExperiment
- edgeR
- limma
- data.table
- dplyr
- tidyr
- purrr
- readr
- stringr
- ggplot2
- scales
- sandwich
- lmtest
- tibble
- statmod

Stochastic procedures use the fixed base seed:

```r
set.seed(23016802)
```

Exact package versions are recorded in `logs/sessioninfo/`.

## Project directory

Scripts use the environment variable:

```text
MAOB_CYP3A_PROJECT_DIR
```

Example in R:

```r
Sys.setenv(
  MAOB_CYP3A_PROJECT_DIR = "/path/to/MAOB_CYP3A4_crosscancer"
)
```

Public-release scripts should not depend on an author-specific absolute path.

## Public data sources

All external data were retrieved on **29 August 2026**. Because these resources are updated over time, re-running the scripts at a later date may return different underlying data. The result tables in `results/` correspond to the retrieval date given here.

### TCGA / NCI Genomic Data Commons

Harmonized TCGA RNA-sequencing data were obtained from the NCI Genomic Data Commons using the STAR - Counts workflow, retrieved on 29 August 2026. The primary analysis uses patient-matched Primary Tumor and Solid Tissue Normal samples from 14 projects, comprising 661 patients with one matched pair each.

Raw TCGA files are **not redistributed in this repository**. The query parameters used to retrieve them are recorded in `results/round2B/tables/101_round2B_configuration.csv`, and the resulting per-project sample counts and aliquot selections are recorded in `results/round2B/tables/120_round2B_paired_project_qc.csv`.

### GTEx v8

The normal-tissue reference is the official GTEx v8 tissue-level gene-median TPM file, retrieved on 29 August 2026:

```text
GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
```

This file is a fixed v8 release and is not expected to change.

### Human Protein Atlas

The HPA-only RNA tissue file `rna_tissue_hpa.tsv.zip` was retrieved on 29 August 2026 and is used for a descriptive comparison of tissue rankings. The HPA consensus dataset is not used, because it incorporates GTEx information and would therefore not be independent of the atlas under comparison.

### UCSC Xena

Previously cached UCSC Xena Toil transformed-expression data are used only for the Round 2C concordance audit. The cache was created during an earlier exploratory stage and was reused rather than re-downloaded, so that the audit is not affected by drift in the source resource.

## Supplementary tables

`Supplementary_Tables_S1_S13.xlsx` contains the supplementary tables as numbered in the manuscript, one sheet per table, with a Contents sheet mapping each sheet to its source file under `results/` and naming the analysis run it came from. Sheets S5, S6a, S6b, S7, S9a, S9b, S12a, S12b and S13 come from Round 3A; S3a and S3b come from Round 2C; the remainder come from Round 2B.

## Important output tables

Round 2B, effect estimates and context:

```text
121_round2B_target_paired_limma_voom_results.csv
140_GTEx_v8_mapping_qc.csv
142_GTEx_v8_target_tissue_specificity.csv
145_round2B_target_normal_abundance_context.csv
151_round2B_tissue_identity_reference_genes.csv
152_round2B_target_residual_vs_tissue_identity.csv
171_round2B_baseline_disruption_organ_level_count.csv
173_round2B_baseline_disruption_pooled_model_count.csv
190_round2B_HPA_only_RNA_target_tissues.csv
194_round2B_HPA_only_validation_status.csv
```

Round 2C, transformed-expression audit:

```text
123_round2C_Xena_count_level_concordance_cells.csv
124_round2C_Xena_count_level_concordance_summary.csv
199_round2C_integrity_checks.csv
```

Round 3A, the inferential results reported in the manuscript:

```text
160_round3A_matched_background_cell_results.csv
162_round3A_denominator_aligned_same_gene_recurrence_null.csv
163_round3A_denominator_aligned_same_gene_profiles.csv
164_round3A_denominator_aligned_coverage_sensitivity.csv
166_round3A_pseudo_target_calibration_values.csv
167_round3A_pseudo_target_calibration_summary.csv
168_round3A_MAOB_organ_collapsed_recurrence_sensitivity.csv
180_round3A_nonLIHC_coexpression_permutation_comparison.csv
181_round3A_coexpression_permutation_summary.csv
199_round3A_integrity_checks.csv
```

Note on file naming: `194_round2B_HPA_only_validation_status.csv` retains an internal name from an earlier stage. The HPA comparison is descriptive concordance of tissue rankings and is not a validation of the GTEx atlas, as stated in the manuscript.

## Figures

Final manuscript figures are generated from frozen analysis outputs.

Main figures:

- Figure 1, normal-tissue expression architecture
- Figure 2, paired tumor–adjacent non-tumor differences
- Figure 3A, tissue-identity context
- Figure 3B, expression-matched background context
- Figure 4A, denominator-aligned same-gene recurrence null at the 80% setting
- Figure 4B, coverage sensitivity for the recurrent-loss and negative-median endpoints
- Figure 4C, organ-collapsed recurrence sensitivity for MAOB

Supplementary figures:

- Supplementary Figure S1A, MAOB–CYP3A coexpression change in non-LIHC cancers
- Supplementary Figure S1B, representative coexpression cells
- Supplementary Figure S2, exploratory normal-baseline relationship
- Supplementary Figure S3, HPA-only normal-tissue RNA expression
- Supplementary Figure S4, pseudo-target calibration of the recurrence null

All figures were rendered at a nominal width of 174 mm with a minimum effective text size of at least 7 pt, verified by a final-size text audit in each figure stage. Categorical colour scales use the Okabe-Ito palette, and the orientation of every continuous colour scale was checked against its endpoint fills.

## Reproducibility and integrity checks

The Round 2B pipeline includes explicit checks for:

- completeness and uniqueness of the cancer-by-gene target grid;
- full-rank paired design matrices;
- inclusion of LIHC MAOB and CYP3A4 in the paired differential-expression analysis;
- exclusion of LIHC from the coexpression analysis;
- target estimability after count filtering;
- GTEx organ mapping;
- matched-background construction;
- same-gene recurrence-null construction;
- figure dimensions and text-size requirements;
- input manifests and session information.

Round 3A adds checks for:

- the frozen differential-expression table retaining its 70 rows and its checksum;
- identity of the target effect estimates between the differential-expression and matched-background outputs;
- selection of exactly one common cancer set per target and coverage setting;
- every exported background profile sharing the selected denominator;
- completeness of the coverage grid;
- export of both the tie-randomized and the +1 pseudo-target calibration values;
- use of the median rule in the organ-collapsed sensitivity analysis;
- coverage of the 40 primary non-LIHC coexpression cells by 5,000 paired label permutations, with no LIHC row in any output.

Results with insufficient expression, sample size, variation, organ mapping, or background support are recorded as not estimable rather than imputed. Not estimable is distinct from evidence of no difference, and the two are kept separate throughout the output tables and figures.

## Important version note

Effect estimates come from **Round 2B**. The inferential results reported in the manuscript for the cross-cancer recurrence analysis, the matched-background comparison and the coexpression analysis come from **Round 3A**, and the transformed-expression audit comes from **Round 2C**.

Earlier or superseded analyses should not be used to reproduce the manuscript results. In particular:

- Stage 1 used transformed TPM values and was exploratory;
- the earlier Round 2 analysis applied a historical LIHC masking rule that excluded LIHC MAOB and CYP3A4 from the paired fit and therefore does not match the final manuscript;
- the concordance audit reported in the manuscript is the Round 2C version covering 43 overlapping cells, not the earlier 41-cell version;
- the Round 2B recurrence null did not align denominators between the target and the background genes, and its probabilities are superseded by those in `results/round3A/`;
- the Round 2B matched-background output used directional probabilities for inference, and is superseded by the two-sided values in `160_round3A_matched_background_cell_results.csv`;
- the Round 2B coexpression output used bootstrap probabilities for inference, and is superseded by the permutation values in `180_round3A_nonLIHC_coexpression_permutation_comparison.csv`.

Neither Stage 1 nor the earlier Round 2 analysis is part of the authoritative publication pipeline, and neither is included here. The superseded Round 2B outputs listed above are retained so that the effect of each revision can be checked.

## Two-paper boundary

This repository pertains to the cross-cancer tumor-versus-adjacent-non-tumor study.

The following are outside this repository and belong to a separate LIHC-focused project:

- LIHC-specific MAOB–CYP3A4 coexpression or correlation;
- MAOB-low/high or CYP3A4-based tumor subgrouping;
- survival analyses;
- Cox models and clinical adjustment;
- HNF4A analyses;
- HCCDB18 or ICGC-LIRI-JP validation;
- pathway enrichment;
- immune infiltration;
- tumor-purity analyses.

No LIHC MAOB–CYP3A correlation is generated by the present pipeline.

## Interpretation boundary

The results do not establish that:

- MAOB loss precedes tumor development;
- MAOB loss causes tumor development;
- MAOB or CYP3A transcript abundance determines GGA abundance;
- transcript abundance directly reflects protein abundance or enzyme activity.

GGA-related metabolism provides the biological rationale for target selection but is not a measured endpoint.

## Citation

A manuscript citation will be added after preprint deposition or journal publication.

```text
[Author list]. Recurrent tumor-associated loss of MAOB transcript expression
and member-dependent CYP3A changes across 14 human cancer types.
[Journal / preprint server]. [Year]. [DOI]
```

## License

Code in this repository is released under the MIT License. See `LICENSE`.

Data originating from TCGA/GDC, GTEx, HPA, and UCSC Xena remain subject to the terms and licenses of the respective source resources.
