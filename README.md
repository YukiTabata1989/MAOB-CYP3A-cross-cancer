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
├── Supplementary_Tables_S1_S11.xlsx
├── code/
│   ├── 01_round2B_primary_analysis.R
│   ├── 02_round2C_xena_audit.R
│   ├── 03_round2D_figure_reconstruction.R
│   ├── 04_round2E_figure_corrections.R
│   ├── 05_round2F_publication_layout.R
│   ├── 06_round2G_final_figure_corrections.R
│   └── 07_round2H_supplementary_figures.R
├── docs/
│   ├── Supplementary_Methods.docx
│   └── Figure_legends.docx
├── figures/
│   ├── Figure1_round2G_GTEx_v8_normal_tissue_atlas.pdf
│   ├── Figure2_round2G_count_level_paired_heatmap.pdf
│   ├── Figure3A_round2F_tissue_identity_distribution.pdf
│   ├── Figure3B_round2F_matched_background_context.pdf
│   ├── Figure4A_round2G_same_gene_recurrence_null.pdf
│   ├── Figure4B_round2F_coverage_sensitivity.pdf
│   ├── SupplementaryFigureS1A_round2G_nonLIHC_coexpression_change_heatmap.pdf
│   ├── SupplementaryFigureS1B_round2G_representative_coexpression_rank_scatter.pdf
│   ├── SupplementaryFigureS2_round2H_normal_baseline_paired_difference.pdf
│   └── SupplementaryFigureS3_round2H_HPA_only_RNA_expression.pdf
├── logs/
│   └── sessioninfo/
│       ├── round2B_sessionInfo.txt
│       ├── round2C_sessionInfo.txt
│       ├── round2F_sessionInfo.txt
│       ├── round2G_sessionInfo.txt
│       └── round2H_sessionInfo.txt
└── results/
    ├── round2B/
    │   ├── logs/
    │   └── tables/
    └── round2C/
        ├── logs/
        └── tables/
```

## Analysis overview

### Round 2B — primary count-level analysis

`01_round2B_primary_analysis.R`

This is the authoritative analysis used for the manuscript results. It:

- retrieves and analyzes harmonized TCGA STAR count data;
- fits paired tumor-versus-adjacent-non-tumor models using `limma-voom`;
- includes LIHC MAOB and CYP3A4 in the paired pan-cancer analysis;
- excludes the entire LIHC project from the MAOB–CYP3A coexpression module;
- constructs GTEx-based tissue-identity reference sets;
- performs expression-property-matched background comparisons;
- performs the covariance-preserving same-gene cross-cancer recurrence analysis;
- evaluates the exploratory normal-baseline/tumor-difference relationship;
- performs non-LIHC MAOB–CYP3A coexpression analyses;
- performs the HPA-only descriptive normal-tissue comparison.

### Round 2C — transformed-expression concordance audit

`02_round2C_xena_audit.R`

This script performs a sensitivity audit using previously cached UCSC Xena Toil transformed-expression data. It does **not** replace the count-level Round 2B analysis.

Round 2C re-ran the audit alone, so that the LIHC MAOB and CYP3A4 cells, which had been masked under an earlier rule, are included. The audit results reported in the manuscript, covering 43 overlapping estimable cells, are those in `results/round2C/`. Round 2C verified by checksum that no Round 2B output was altered.

### Round 2D — initial manuscript figure reconstruction

`03_round2D_figure_reconstruction.R`

Reconstructs manuscript figures from frozen Round 2B and Round 2C outputs without refitting the statistical analyses.

### Round 2E — targeted figure corrections

`04_round2E_figure_corrections.R`

Applies targeted corrections to Figure 1, Figure 3B, and Supplementary Figure S1B while retaining frozen analytical values.

### Round 2F — publication-size figure layout

`05_round2F_publication_layout.R`

Renders the manuscript figures at the target publication width of 174 mm and the common minimum text-size requirement.

### Round 2G — final corrections for selected figures

`06_round2G_final_figure_corrections.R`

Produces the final versions of selected figures, including corrections to continuous color bars and figure annotations, without recomputing biological or statistical results.

### Round 2H — final Supplementary Figures S2 and S3

`07_round2H_supplementary_figures.R`

Produces the final versions of Supplementary Figure S2 and Supplementary Figure S3 from frozen Round 2B outputs.

## Recommended execution order

```text
1. Round 2B
2. Round 2C
3. Round 2D
4. Round 2E
5. Round 2F
6. Round 2G
7. Round 2H
```

Round 2D to 2H are figure-generation and revision stages and do not replace the primary statistical analysis.

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

`Supplementary_Tables_S1_S11.xlsx` contains the supplementary tables as numbered in the manuscript, one sheet per table, with a Contents sheet mapping each sheet to its source file under `results/`. All sheets derive from Round 2B except S3a and S3b, which come from the Round 2C audit.

## Important output tables

Important Round 2B outputs include:

```text
121_round2B_target_paired_limma_voom_results.csv
140_GTEx_v8_mapping_qc.csv
142_GTEx_v8_target_tissue_specificity.csv
145_round2B_target_normal_abundance_context.csv
151_round2B_tissue_identity_reference_genes.csv
152_round2B_target_residual_vs_tissue_identity.csv
160_round2B_matched_background_cell_results.csv
162_round2B_same_gene_recurrence_null.csv
163_round2B_same_gene_null_gene_profiles.csv
164_round2B_same_gene_recurrence_sensitivity.csv
171_round2B_baseline_disruption_organ_level_count.csv
173_round2B_baseline_disruption_pooled_model_count.csv
180_round2B_nonLIHC_coexpression_rewiring.csv
181_round2B_representative_coexpression_cells.csv
182_round2B_representative_coexpression_points.csv
190_round2B_HPA_only_RNA_target_tissues.csv
194_round2B_HPA_only_validation_status.csv
```

Important Round 2C outputs include:

```text
123_round2C_Xena_count_level_concordance_cells.csv
124_round2C_Xena_count_level_concordance_summary.csv
199_round2C_integrity_checks.csv
```

Note on file naming: `194_round2B_HPA_only_validation_status.csv` retains an internal name from an earlier stage. The HPA comparison is descriptive concordance of tissue rankings and is not a validation of the GTEx atlas, as stated in the manuscript.

## Figures

Final manuscript figures are generated from frozen analysis outputs.

Main figures:

- Figure 1, normal-tissue expression architecture
- Figure 2, paired tumor–adjacent non-tumor differences
- Figure 3A, tissue-identity context
- Figure 3B, expression-matched background context
- Figure 4A, same-gene recurrence null
- Figure 4B, coverage sensitivity of the recurrence test

Supplementary figures:

- Supplementary Figure S1A, MAOB–CYP3A coexpression change in non-LIHC cancers
- Supplementary Figure S1B, representative coexpression cells
- Supplementary Figure S2, exploratory normal-baseline relationship
- Supplementary Figure S3, HPA-only normal-tissue RNA expression

All figures were rendered at a nominal width of 174 mm with a minimum effective text size of at least 7 pt. Categorical color scales use the Okabe-Ito palette.

## Reproducibility and integrity checks

The pipeline includes explicit checks for:

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

Results with insufficient expression, sample size, variation, organ mapping, or background support are recorded as not estimable rather than imputed. Not estimable is distinct from evidence of no difference, and the two are kept separate throughout the output tables and figures.

## Important version note

The manuscript results correspond to **Round 2B and later**.

Earlier exploratory or superseded analyses should not be used to reproduce the manuscript results. In particular:

- Stage 1 used transformed TPM values and was exploratory;
- the earlier Round 2 analysis applied a historical LIHC masking rule that excluded LIHC MAOB and CYP3A4 from the paired fit and therefore does not match the final manuscript;
- the concordance audit reported in the manuscript is the Round 2C version covering 43 overlapping cells, not the earlier 41-cell version.

Neither Stage 1 nor the earlier Round 2 analysis is part of the authoritative publication pipeline, and neither is included here.

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
