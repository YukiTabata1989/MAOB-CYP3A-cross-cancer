#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 2C: Xena transformed-expression concordance audit only
# =============================================================================
#
# PURPOSE
#   Recompute only the UCSC Xena Toil transformed-expression audit after
#   returning LIHC-MAOB and LIHC-CYP3A4 to Paper 1. The Stage 1 raw Xena caches
#   are reused before the historical LIHC mask was applied.
#
# CHANGES FROM THE PREVIOUS AUDIT
#   1) LIHC-MAOB and LIHC-CYP3A4 are included.
#   2) All 70 Xena cancer-by-gene cells are retested.
#   3) Benjamini-Hochberg correction is recalculated across all 70 finite Xena
#      P values; no previous Stage 1 q value is reused.
#   4) Concordance with the authoritative Round 2B count-level estimates is
#      evaluated across the 43 cells estimable in both analyses.
#   5) Xena and GDC patient sets are paired within source but are not forced to
#      be identical. Source-specific pair counts and case-set differences are
#      recorded explicitly for every cell.
#
# OUT OF SCOPE / NOT EXECUTED
#   This script does not fit or refit count-level limma-voom models and does not
#   run tissue-identity, matched-background, same-gene recurrence, baseline,
#   coexpression, HPA, survival, clinical, subgroup, immune, pathway, or tumor-
#   purity analyses. In particular, no LIHC coexpression vector is accessed.
#
# INPUTS EXPECTED UNDER MAOB_CYP3A_PROJECT_DIR
#   data/raw/xena_toil_<GENE>.rds, for all five target genes
#   data/raw/TcgaTargetGTEX_phenotype.txt.gz
#   data/processed/round2_project_cache/TCGA-*/paired_raw_input.rds
#   results/round2B/tables/121_round2B_target_paired_limma_voom_results.csv
#
# OUTPUTS
#   results/round2C/tables/123_round2C_Xena_count_level_concordance_cells.csv
#   results/round2C/tables/124_round2C_Xena_count_level_concordance_summary.csv
#   results/round2C/tables/199_round2C_integrity_checks.csv
#   results/round2C/logs/round2C_progress.log
#   results/round2C/logs/round2C_sessionInfo.txt
#
# DATA ACCESS NOTE
#   UCSCXenaTools/UCSCXenaShiny is intentionally not used here. Reusing the
#   exact Stage 1 raw cache avoids source drift and guarantees that the only
#   analytical change is removal of the historical two-cell LIHC mask.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration --------------------------------------------------------

analysis_id <- "MAOB_CYP3A_Xena_audit_round2C"
required_r_version <- "4.6.1"

target_genes <- c("MAOB", "CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
target_grid_n <- 14L * length(target_genes)
expected_xena_bh_n <- 70L
expected_count_level_overlap <- 43L
expected_xena_pair_total <- 644L
effect_shift_threshold <- log2(1.5)
minimum_pairs_primary <- 10L

project_dir <- normalizePath(
  Sys.getenv("MAOB_CYP3A_PROJECT_DIR", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)

raw_dir <- file.path(project_dir, "data", "raw")
cache_dir <- file.path(
  project_dir, "data", "processed", "round2_project_cache"
)
round2b_result_dir <- file.path(project_dir, "results", "round2B")
round2b_table_dir <- file.path(round2b_result_dir, "tables")
round2c_result_dir <- file.path(project_dir, "results", "round2C")
table_dir <- file.path(round2c_result_dir, "tables")
log_dir <- file.path(round2c_result_dir, "logs")

invisible(lapply(
  c(round2c_result_dir, table_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

round2b_target_file <- file.path(
  round2b_table_dir,
  "121_round2B_target_paired_limma_voom_results.csv"
)
phenotype_file <- file.path(raw_dir, "TcgaTargetGTEX_phenotype.txt.gz")
xena_cache_files <- setNames(
  file.path(raw_dir, paste0("xena_toil_", target_genes, ".rds")),
  target_genes
)

cells_output_file <- file.path(
  table_dir,
  "123_round2C_Xena_count_level_concordance_cells.csv"
)
summary_output_file <- file.path(
  table_dir,
  "124_round2C_Xena_count_level_concordance_summary.csv"
)
integrity_output_file <- file.path(
  table_dir,
  "199_round2C_integrity_checks.csv"
)
progress_log_file <- file.path(log_dir, "round2C_progress.log")
session_info_file <- file.path(log_dir, "round2C_sessionInfo.txt")

required_packages <- c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "tibble"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 2C:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 2C requires R ", required_r_version,
    "; current version is ", as.character(getRversion()), ".",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
})

# ---- 1. Utilities ------------------------------------------------------------

timestamp_utc <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

replace_with_temp_file <- function(tmp, path) {
  if (!file.exists(tmp)) {
    stop("Temporary output was not created: ", tmp, call. = FALSE)
  }
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not remove the previous output: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Could not move completed temporary output to: ", path, call. = FALSE)
  }
  invisible(path)
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(x, tmp, na = "")
  replace_with_temp_file(tmp, path)
}

write_lines_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(x, tmp, useBytes = TRUE)
  replace_with_temp_file(tmp, path)
}

write_lines_atomic(
  paste0("[", timestamp_utc(), "] Starting ", analysis_id),
  progress_log_file
)

log_message <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", timestamp_utc(), "] ", msg)
  message(msg)
  cat(line, "\n", file = progress_log_file, append = TRUE, sep = "")
  invisible(msg)
}

assert_true <- function(condition, message_text) {
  if (!isTRUE(condition)) stop(message_text, call. = FALSE)
  invisible(TRUE)
}

safe_correlation <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(unname(stats::cor(x[ok], y[ok], method = method)))
}

md5_snapshot <- function(directory) {
  paths <- sort(list.files(
    directory,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  ))
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  if (length(paths) == 0L) {
    stop("No files were found for checksum snapshot: ", directory, call. = FALSE)
  }
  values <- unname(tools::md5sum(paths))
  names(values) <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  values
}

required_input_files <- c(
  round2b_target_file,
  phenotype_file,
  unname(xena_cache_files)
)
missing_input_files <- required_input_files[!file.exists(required_input_files)]
if (length(missing_input_files) > 0L) {
  stop(
    "Required Round 2C input files are missing:\n  ",
    paste(missing_input_files, collapse = "\n  "),
    call. = FALSE
  )
}

round2b_md5_before <- md5_snapshot(round2b_result_dir)
log_message(
  "Input paths validated; checksum snapshot contains ",
  length(round2b_md5_before), " Round 2B files."
)

# ---- 2. Fixed cancer mapping and authoritative Round 2B table ---------------

cancer_key <- tibble::tribble(
  ~tcga_code, ~detailed_category,
  "BLCA", "Bladder Urothelial Carcinoma",
  "BRCA", "Breast Invasive Carcinoma",
  "COAD", "Colon Adenocarcinoma",
  "ESCA", "Esophageal Carcinoma",
  "HNSC", "Head & Neck Squamous Cell Carcinoma",
  "KICH", "Kidney Chromophobe",
  "KIRC", "Kidney Clear Cell Carcinoma",
  "KIRP", "Kidney Papillary Cell Carcinoma",
  "LIHC", "Liver Hepatocellular Carcinoma",
  "LUAD", "Lung Adenocarcinoma",
  "LUSC", "Lung Squamous Cell Carcinoma",
  "PRAD", "Prostate Adenocarcinoma",
  "STAD", "Stomach Adenocarcinoma",
  "THCA", "Thyroid Carcinoma"
)

tcga_codes <- cancer_key$tcga_code
full_grid <- tidyr::expand_grid(
  tcga_code = factor(tcga_codes, levels = tcga_codes),
  gene = factor(target_genes, levels = target_genes)
) %>%
  mutate(
    tcga_code = as.character(tcga_code),
    gene = as.character(gene)
  )

round2b_target <- readr::read_csv(
  round2b_target_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_round2b_columns <- c(
  "tcga_code", "gene", "n_pairs", "logFC", "q_value_global",
  "analysis_status", "robust_count_level_change",
  "count_level_interpretation"
)
assert_true(
  all(required_round2b_columns %in% names(round2b_target)),
  paste(
    "Round 2B target table lacks required columns:",
    paste(setdiff(required_round2b_columns, names(round2b_target)), collapse = ", ")
  )
)
assert_true(nrow(round2b_target) == target_grid_n, "Round 2B target grid is not 70 rows.")
assert_true(
  !anyDuplicated(round2b_target[c("tcga_code", "gene")]),
  "Round 2B target table contains duplicate cancer-gene cells."
)
assert_true(
  setequal(round2b_target$tcga_code, tcga_codes) &&
    setequal(round2b_target$gene, target_genes),
  "Round 2B target table does not contain the expected 14 cancers and five genes."
)

round2b_target <- full_grid %>%
  left_join(round2b_target, by = c("tcga_code", "gene"))

log_message(
  "Loaded authoritative Round 2B table: ", nrow(round2b_target),
  " cells; ", sum(is.finite(round2b_target$logFC)),
  " finite count-level effects."
)

# ---- 3. Round 2B patient sets from cached paired raw inputs -----------------

load_round2b_cases <- function(tcga_code) {
  cache_file <- file.path(
    cache_dir,
    paste0("TCGA-", tcga_code),
    "paired_raw_input.rds"
  )
  if (!file.exists(cache_file)) {
    stop("Missing Round 2B paired-input cache: ", cache_file, call. = FALSE)
  }
  obj <- readRDS(cache_file)
  if (is.null(obj$sample_meta)) {
    stop("sample_meta is absent from: ", cache_file, call. = FALSE)
  }
  sample_meta <- tibble::as_tibble(obj$sample_meta)
  if (!all(c("case_id", "condition") %in% names(sample_meta))) {
    stop("case_id/condition is absent from: ", cache_file, call. = FALSE)
  }
  sample_meta <- sample_meta %>%
    transmute(
      tcga_code = tcga_code,
      case_id = toupper(as.character(case_id)),
      condition = as.character(condition)
    ) %>%
    filter(
      !is.na(case_id), nzchar(case_id),
      condition %in% c("Normal", "Tumor")
    )

  patient_qc <- sample_meta %>%
    count(case_id, condition, name = "n_samples") %>%
    tidyr::complete(
      case_id,
      condition = c("Normal", "Tumor"),
      fill = list(n_samples = 0L)
    ) %>%
    group_by(case_id) %>%
    summarise(
      n_normal = n_samples[condition == "Normal"],
      n_tumor = n_samples[condition == "Tumor"],
      .groups = "drop"
    )
  if (any(patient_qc$n_normal != 1L | patient_qc$n_tumor != 1L)) {
    stop(
      "Round 2B cache is not exactly paired for project ", tcga_code, ".",
      call. = FALSE
    )
  }
  tibble::tibble(
    tcga_code = tcga_code,
    case_id = sort(unique(patient_qc$case_id))
  )
}

round2b_cases <- bind_rows(lapply(tcga_codes, load_round2b_cases))
round2b_case_counts <- round2b_cases %>%
  count(tcga_code, name = "round2B_patient_cases")

round2b_declared_counts <- round2b_target %>%
  distinct(tcga_code, n_pairs) %>%
  rename(round2B_declared_pairs = n_pairs)
assert_true(
  nrow(round2b_declared_counts) == length(tcga_codes),
  "Round 2B n_pairs is not constant within every project."
)

round2b_case_count_check <- round2b_case_counts %>%
  left_join(round2b_declared_counts, by = "tcga_code")
assert_true(
  all(
    round2b_case_count_check$round2B_patient_cases ==
      round2b_case_count_check$round2B_declared_pairs
  ),
  "Round 2B patient-cache counts do not match the 121 table."
)
assert_true(
  sum(round2b_case_counts$round2B_patient_cases) == 661L,
  "Round 2B patient caches do not contain the expected 661 pairs."
)
log_message("Validated Round 2B patient caches: 661 paired cases.")

# ---- 4. Load Stage 1 raw Xena values and phenotype --------------------------

load_xena_gene_cache <- function(gene) {
  cache_file <- xena_cache_files[[gene]]
  x <- tibble::as_tibble(readRDS(cache_file))
  required <- c("sample", "xena_log2_tpm_p001")
  if (!all(required %in% names(x))) {
    stop(
      "Raw Xena cache lacks sample/xena_log2_tpm_p001: ",
      cache_file,
      call. = FALSE
    )
  }
  if ("gene" %in% names(x)) {
    observed_genes <- unique(as.character(x$gene[!is.na(x$gene)]))
    if (length(observed_genes) > 0L && !identical(observed_genes, gene)) {
      stop("Unexpected gene label in raw Xena cache: ", cache_file, call. = FALSE)
    }
  }
  out <- x %>%
    transmute(
      sample = as.character(sample),
      gene = gene,
      xena_log2_tpm_p001 = as.numeric(xena_log2_tpm_p001),
      tpm_recomputed = pmax(2^xena_log2_tpm_p001 - 0.001, 0),
      log2_tpm1_recomputed = log2(tpm_recomputed + 1)
    )
  if (anyDuplicated(out[c("sample", "gene")])) {
    stop("Duplicate sample-gene rows in raw Xena cache: ", cache_file, call. = FALSE)
  }
  out
}

xena_expression <- bind_rows(lapply(target_genes, load_xena_gene_cache))
log_message(
  "Loaded Stage 1 raw Xena caches: ", nrow(xena_expression),
  " sample-gene values."
)

phenotype <- data.table::fread(
  phenotype_file,
  sep = "\t",
  encoding = "Latin-1",
  data.table = FALSE,
  showProgress = FALSE
)
expected_phenotype_columns <- c(
  "sample", "detailed_category", "primary_tissue", "primary_site",
  "sample_type", "sex", "study"
)
assert_true(
  ncol(phenotype) == length(expected_phenotype_columns),
  paste0(
    "Unexpected Xena phenotype structure: expected ",
    length(expected_phenotype_columns), " columns but found ", ncol(phenotype), "."
  )
)
names(phenotype) <- expected_phenotype_columns
phenotype <- tibble::as_tibble(phenotype) %>%
  mutate(sample = as.character(sample))
assert_true(!anyDuplicated(phenotype$sample), "Xena phenotype contains duplicate sample IDs.")

xena_joined <- xena_expression %>%
  inner_join(phenotype, by = "sample")
matched_fraction <- nrow(xena_joined) / nrow(xena_expression)
assert_true(
  is.finite(matched_fraction) && matched_fraction >= 0.99,
  paste0(
    "Less than 99% of Xena expression rows matched phenotype metadata: ",
    signif(matched_fraction, 4)
  )
)

xena_tcga <- xena_joined %>%
  filter(
    study == "TCGA",
    sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
  ) %>%
  inner_join(cancer_key, by = "detailed_category") %>%
  mutate(case_id = toupper(substr(sample, 1L, 12L))) %>%
  filter(tcga_code %in% tcga_codes)

assert_true(nrow(xena_tcga) > 0L, "No eligible TCGA rows were found in Xena caches.")
log_message(
  "Mapped Xena expression to the 14 audit cancers: ", nrow(xena_tcga),
  " sample-gene rows."
)

# ---- 5. Recompute all 70 Xena paired tests ----------------------------------

# This reproduces the Stage 1 rule: multiple technical/aliquot values within a
# patient, gene, and sample type are collapsed by their median before pairing.
xena_case_values <- xena_tcga %>%
  group_by(tcga_code, gene, case_id, sample_type) %>%
  summarise(
    log2_tpm1 = median(log2_tpm1_recomputed, na.rm = TRUE),
    n_aliquots = n(),
    .groups = "drop"
  )

xena_paired <- xena_case_values %>%
  select(tcga_code, gene, case_id, sample_type, log2_tpm1) %>%
  pivot_wider(names_from = sample_type, values_from = log2_tpm1) %>%
  rename(
    tumor_log2_tpm1 = `Primary Tumor`,
    normal_log2_tpm1 = `Solid Tissue Normal`
  ) %>%
  filter(
    is.finite(tumor_log2_tpm1),
    is.finite(normal_log2_tpm1)
  ) %>%
  mutate(paired_shift = tumor_log2_tpm1 - normal_log2_tpm1)

assert_true(
  !anyDuplicated(xena_paired[c("tcga_code", "gene", "case_id")]),
  "Xena paired table contains duplicate cancer-gene-patient rows."
)

xena_pair_counts <- xena_paired %>%
  count(tcga_code, gene, name = "xena_pairs")
assert_true(
  nrow(xena_pair_counts) == target_grid_n,
  "Xena paired data do not contain all 70 cancer-gene cells."
)
xena_project_pair_counts <- xena_pair_counts %>%
  group_by(tcga_code) %>%
  summarise(
    n_gene_cells = n(),
    n_distinct_pair_counts = n_distinct(xena_pairs),
    xena_pairs = first(xena_pairs),
    .groups = "drop"
  )
assert_true(
  all(
    xena_project_pair_counts$n_gene_cells == length(target_genes) &
      xena_project_pair_counts$n_distinct_pair_counts == 1L
  ),
  "Xena patient-pair counts are not identical across the five genes within a project."
)
assert_true(
  sum(xena_project_pair_counts$xena_pairs) == expected_xena_pair_total,
  paste0(
    "Expected ", expected_xena_pair_total,
    " source-specific Xena pairs across 14 projects, but found ",
    sum(xena_project_pair_counts$xena_pairs), "."
  )
)

xena_test_results <- xena_paired %>%
  group_by(tcga_code, gene) %>%
  group_modify(~ {
    shift <- .x$paired_shift
    shift <- shift[is.finite(shift)]
    n_pairs <- length(shift)
    if (n_pairs < minimum_pairs_primary) {
      return(tibble::tibble(
        stage1_n_pairs = n_pairs,
        stage1_shift = if (n_pairs > 0L) median(shift) else NA_real_,
        stage1_mean_shift = if (n_pairs > 0L) mean(shift) else NA_real_,
        stage1_p_value = NA_real_,
        stage1_status = "insufficient_pairs_lt10"
      ))
    }
    wt <- tryCatch(
      suppressWarnings(stats::wilcox.test(
        shift,
        mu = 0,
        alternative = "two.sided",
        exact = FALSE
      )),
      error = function(e) e
    )
    if (inherits(wt, "error") || !is.finite(unname(wt$p.value))) {
      return(tibble::tibble(
        stage1_n_pairs = n_pairs,
        stage1_shift = median(shift),
        stage1_mean_shift = mean(shift),
        stage1_p_value = NA_real_,
        stage1_status = "wilcoxon_not_estimable"
      ))
    }
    tibble::tibble(
      stage1_n_pairs = n_pairs,
      stage1_shift = median(shift),
      stage1_mean_shift = mean(shift),
      stage1_p_value = unname(wt$p.value),
      stage1_status = "estimable_xena"
    )
  }) %>%
  ungroup()

xena_test_results <- full_grid %>%
  left_join(xena_test_results, by = c("tcga_code", "gene")) %>%
  mutate(
    stage1_status = if_else(
      is.na(stage1_status),
      "missing_xena_paired_cell",
      stage1_status
    )
  )

xena_bh_eligible <- xena_test_results$stage1_status == "estimable_xena" &
  is.finite(xena_test_results$stage1_p_value)
xena_bh_n <- sum(xena_bh_eligible)
assert_true(
  xena_bh_n == expected_xena_bh_n,
  paste0(
    "Expected 70 finite Xena tests for BH correction, but found ",
    xena_bh_n, ". No q values were written."
  )
)

xena_test_results$stage1_q_value_bh70 <- NA_real_
xena_test_results$stage1_q_value_bh70[xena_bh_eligible] <- stats::p.adjust(
  xena_test_results$stage1_p_value[xena_bh_eligible],
  method = "BH"
)
xena_test_results <- xena_test_results %>%
  mutate(
    stage1_passes_statistical_screen =
      is.finite(stage1_q_value_bh70) & stage1_q_value_bh70 < 0.05,
    stage1_passes_effect_screen =
      is.finite(stage1_shift) & abs(stage1_shift) >= effect_shift_threshold,
    stage1_robust =
      stage1_passes_statistical_screen & stage1_passes_effect_screen,
    stage1_interpretation = case_when(
      stage1_status != "estimable_xena" ~ "not_estimable_xena",
      stage1_robust & stage1_shift < 0 ~ "loss",
      stage1_robust & stage1_shift > 0 ~ "gain",
      TRUE ~ "not_robust"
    )
  )

log_message(
  "Recomputed ", xena_bh_n,
  " Xena Wilcoxon tests and BH-adjusted all 70 P values."
)

# ---- 6. Compare source-specific patient sets -------------------------------

xena_cases_by_cell <- xena_paired %>%
  distinct(tcga_code, gene, case_id)

compare_case_sets <- function(tcga_code, gene) {
  xena_set <- sort(unique(xena_cases_by_cell$case_id[
    xena_cases_by_cell$tcga_code == tcga_code &
      xena_cases_by_cell$gene == gene
  ]))
  round2b_set <- sort(unique(round2b_cases$case_id[
    round2b_cases$tcga_code == tcga_code
  ]))
  n_common <- length(intersect(xena_set, round2b_set))
  n_xena_only <- length(setdiff(xena_set, round2b_set))
  n_round2b_only <- length(setdiff(round2b_set, xena_set))
  tibble::tibble(
    tcga_code = tcga_code,
    gene = gene,
    xena_patient_cases = length(xena_set),
    round2B_patient_cases = length(round2b_set),
    n_common_patient_cases = n_common,
    n_xena_only_cases = n_xena_only,
    n_round2B_only_cases = n_round2b_only,
    exact_patient_set_match = setequal(xena_set, round2b_set),
    xena_patient_set_key = paste(xena_set, collapse = ";"),
    round2B_patient_set_key = paste(round2b_set, collapse = ";"),
    patient_set_relation = if_else(
      setequal(xena_set, round2b_set),
      "exact_match",
      "source_specific_mismatch_recorded"
    )
  )
}

patient_set_comparison <- bind_rows(lapply(seq_len(nrow(full_grid)), function(i) {
  compare_case_sets(full_grid$tcga_code[[i]], full_grid$gene[[i]])
}))

mismatched_projects <- patient_set_comparison %>%
  filter(!exact_patient_set_match) %>%
  distinct(tcga_code) %>%
  pull(tcga_code)
log_message(
  "Patient-set audit complete. Exact Xena/GDC sets differ in: ",
  ifelse(length(mismatched_projects) == 0L,
         "none", paste(mismatched_projects, collapse = ", ")),
  ". Differences are recorded and are not treated as execution errors."
)

# ---- 7. Cell-level concordance and leverage-aware summary ------------------

round2b_for_join <- round2b_target %>%
  transmute(
    tcga_code,
    gene,
    round2B_n_pairs = n_pairs,
    round2B_analysis_status = analysis_status,
    round2B_logFC = logFC,
    round2B_q = q_value_global,
    round2B_robust = as.logical(robust_count_level_change),
    round2B_interpretation = count_level_interpretation
  )

concordance_cells <- full_grid %>%
  left_join(round2b_for_join, by = c("tcga_code", "gene")) %>%
  left_join(xena_test_results, by = c("tcga_code", "gene")) %>%
  left_join(patient_set_comparison, by = c("tcga_code", "gene")) %>%
  mutate(
    finite_overlap = is.finite(stage1_shift) & is.finite(round2B_logFC),
    direction_concordant = case_when(
      finite_overlap ~ sign(stage1_shift) == sign(round2B_logFC),
      TRUE ~ NA
    ),
    robust_call_agreement = case_when(
      finite_overlap ~ stage1_robust == round2B_robust,
      TRUE ~ NA
    ),
    status = case_when(
      stage1_status != "estimable_xena" ~ paste0("xena_", stage1_status),
      !is.finite(round2B_logFC) ~ paste0(
        "count_level_not_estimable:", round2B_analysis_status
      ),
      TRUE ~ "overlap_estimable"
    )
  ) %>%
  select(
    # Existing Round 2B concordance columns are retained first.
    tcga_code, gene,
    round2B_logFC, round2B_q, round2B_robust, round2B_interpretation,
    stage1_shift, stage1_robust, stage1_interpretation,
    finite_overlap, direction_concordant,
    # Round 2C additions.
    robust_call_agreement,
    round2B_n_pairs, stage1_n_pairs,
    stage1_mean_shift, stage1_p_value, stage1_q_value_bh70,
    stage1_passes_statistical_screen, stage1_passes_effect_screen,
    round2B_analysis_status, stage1_status,
    xena_patient_cases, round2B_patient_cases,
    n_common_patient_cases, n_xena_only_cases, n_round2B_only_cases,
    exact_patient_set_match, patient_set_relation,
    status
  )

assert_true(nrow(concordance_cells) == target_grid_n, "Concordance output is not 70 rows.")
assert_true(
  !anyDuplicated(concordance_cells[c("tcga_code", "gene")]),
  "Concordance output contains duplicate cells."
)

overlap <- concordance_cells %>% filter(finite_overlap)
assert_true(
  nrow(overlap) == expected_count_level_overlap,
  paste0(
    "Expected 43 Xena/count-level overlapping cells, but found ",
    nrow(overlap), "."
  )
)

lihc_required <- overlap %>%
  filter(tcga_code == "LIHC", gene %in% c("MAOB", "CYP3A4"))
assert_true(
  nrow(lihc_required) == 2L &&
    all(is.finite(lihc_required$stage1_p_value)) &&
    all(is.finite(lihc_required$stage1_q_value_bh70)),
  "LIHC-MAOB and LIHC-CYP3A4 were not both restored to the Xena audit."
)

overlap_without_lihc_cyp3a4 <- overlap %>%
  filter(!(tcga_code == "LIHC" & gene == "CYP3A4"))

concordance_summary <- tibble::tibble(
  n_overlap = nrow(overlap),
  direction_concordance_fraction = mean(
    overlap$direction_concordant,
    na.rm = TRUE
  ),
  robust_call_agreement = mean(
    overlap$robust_call_agreement,
    na.rm = TRUE
  ),
  pearson_r = safe_correlation(
    overlap$stage1_shift,
    overlap$round2B_logFC,
    method = "pearson"
  ),
  spearman_rho = safe_correlation(
    overlap$stage1_shift,
    overlap$round2B_logFC,
    method = "spearman"
  ),
  n_overlap_excluding_LIHC_CYP3A4 = nrow(overlap_without_lihc_cyp3a4),
  pearson_r_excluding_LIHC_CYP3A4 = safe_correlation(
    overlap_without_lihc_cyp3a4$stage1_shift,
    overlap_without_lihc_cyp3a4$round2B_logFC,
    method = "pearson"
  ),
  spearman_rho_excluding_LIHC_CYP3A4 = safe_correlation(
    overlap_without_lihc_cyp3a4$stage1_shift,
    overlap_without_lihc_cyp3a4$round2B_logFC,
    method = "spearman"
  ),
  xena_BH_tests = xena_bh_n,
  interpretation = paste(
    "Source-specific paired concordance audit between transformed Xena Toil",
    "expression and authoritative GDC count-level Round 2B estimates;",
    "LIHC-MAOB and LIHC-CYP3A4 included; patient sets were paired within",
    "source but were not required to be identical; not an independent",
    "validation dataset"
  )
)

log_message(
  "Concordance summary: ", nrow(overlap),
  " overlapping cells; direction agreement=",
  signif(concordance_summary$direction_concordance_fraction, 4),
  "; robust-call agreement=",
  signif(concordance_summary$robust_call_agreement, 4), "."
)

# ---- 8. Write outputs and verify the audit-only boundary --------------------

write_csv_atomic(concordance_cells, cells_output_file)
write_csv_atomic(concordance_summary, summary_output_file)

round2b_md5_after <- md5_snapshot(round2b_result_dir)
round2b_files_unchanged <- identical(round2b_md5_before, round2b_md5_after)

project_patient_checks <- patient_set_comparison %>%
  group_by(tcga_code) %>%
  summarise(
    gene_cells = n(),
    all_gene_cells_have_same_xena_case_set =
      n_distinct(xena_patient_set_key) == 1L &&
      n_distinct(round2B_patient_set_key) == 1L &&
      n_distinct(xena_patient_cases) == 1L &&
      n_distinct(n_common_patient_cases) == 1L &&
      n_distinct(n_xena_only_cases) == 1L &&
      n_distinct(n_round2B_only_cases) == 1L,
    exact_match_all_gene_cells = all(exact_patient_set_match),
    xena_pairs = first(xena_patient_cases),
    round2B_pairs = first(round2B_patient_cases),
    common_pairs = first(n_common_patient_cases),
    xena_only = first(n_xena_only_cases),
    round2B_only = first(n_round2B_only_cases),
    .groups = "drop"
  ) %>%
  transmute(
    check = paste0("patient_set_", tcga_code),
    result = if_else(
      exact_match_all_gene_cells,
      "exact_match",
      "source_specific_mismatch_recorded"
    ),
    expected = "comparison_recorded_for_all_5_genes",
    passed = gene_cells == 5L & all_gene_cells_have_same_xena_case_set,
    details = paste0(
      "Xena=", xena_pairs,
      "; Round2B=", round2B_pairs,
      "; common=", common_pairs,
      "; Xena_only=", xena_only,
      "; Round2B_only=", round2B_only
    )
  )

integrity_checks <- bind_rows(
  tibble::tribble(
    ~check, ~result, ~expected, ~passed, ~details,
    "R_version",
    as.character(getRversion()), required_r_version,
    as.character(getRversion()) == required_r_version,
    "Exact requested runtime version",
    "round2B_target_grid_rows",
    as.character(nrow(round2b_target)), "70", nrow(round2b_target) == 70L,
    "Authoritative 121 table was read only",
    "round2B_patient_cache_total",
    as.character(sum(round2b_case_counts$round2B_patient_cases)), "661",
    sum(round2b_case_counts$round2B_patient_cases) == 661L,
    "Counts agree with the 121 table",
    "Xena_source_pair_total",
    as.character(sum(xena_project_pair_counts$xena_pairs)),
    as.character(expected_xena_pair_total),
    sum(xena_project_pair_counts$xena_pairs) == expected_xena_pair_total,
    "Stage 1 source-specific pairs across the 14 audit projects",
    "round2C_output_grid_rows",
    as.character(nrow(concordance_cells)), "70",
    nrow(concordance_cells) == 70L,
    "Full 14-cancer by 5-gene grid",
    "round2C_unique_cells",
    as.character(!anyDuplicated(concordance_cells[c("tcga_code", "gene")])),
    "TRUE",
    !anyDuplicated(concordance_cells[c("tcga_code", "gene")]),
    "No duplicated cancer-gene cells",
    "Xena_BH_test_count",
    as.character(xena_bh_n), "70", xena_bh_n == 70L,
    "BH recalculated from all 70 new Xena P values",
    "count_level_overlap_cells",
    as.character(nrow(overlap)), "43", nrow(overlap) == 43L,
    "Finite Xena shift and finite Round 2B logFC",
    "LIHC_MAOB_included",
    as.character(any(
      overlap$tcga_code == "LIHC" & overlap$gene == "MAOB"
    )),
    "TRUE",
    any(overlap$tcga_code == "LIHC" & overlap$gene == "MAOB"),
    "Paired-expression audit only",
    "LIHC_CYP3A4_included",
    as.character(any(
      overlap$tcga_code == "LIHC" & overlap$gene == "CYP3A4"
    )),
    "TRUE",
    any(overlap$tcga_code == "LIHC" & overlap$gene == "CYP3A4"),
    "Paired-expression audit only",
    "patient_set_comparisons_recorded",
    as.character(nrow(patient_set_comparison)), "70",
    nrow(patient_set_comparison) == 70L,
    "Source-specific patient differences are retained, not imputed",
    "round2B_files_unchanged",
    as.character(round2b_files_unchanged), "TRUE",
    round2b_files_unchanged,
    paste0(length(round2b_md5_before), " Round 2B files checked by MD5"),
    "main_analysis_modules_recomputed",
    "FALSE", "FALSE", TRUE,
    paste(
      "No limma-voom, tissue identity, matched background, recurrence,",
      "baseline, coexpression, or HPA module is present"
    ),
    "forbidden_LIHC_or_clinical_modules_run",
    "FALSE", "FALSE", TRUE,
    paste(
      "No LIHC coexpression, survival, clinical, subgroup, immune, pathway,",
      "or purity analysis is present"
    )
  ),
  project_patient_checks
)

write_csv_atomic(integrity_checks, integrity_output_file)

critical_integrity <- integrity_checks %>%
  filter(!grepl("^patient_set_", check))
if (any(!critical_integrity$passed)) {
  failed <- critical_integrity$check[!critical_integrity$passed]
  stop(
    "Round 2C integrity checks failed: ",
    paste(failed, collapse = ", "),
    call. = FALSE
  )
}
if (any(!project_patient_checks$passed)) {
  failed <- project_patient_checks$check[!project_patient_checks$passed]
  stop(
    "Patient-set comparison was inconsistent across genes in: ",
    paste(failed, collapse = ", "),
    call. = FALSE
  )
}

write_lines_atomic(capture.output(sessionInfo()), session_info_file)

log_message("Round 2C complete.")
log_message("Cells:     ", cells_output_file)
log_message("Summary:   ", summary_output_file)
log_message("Integrity: ", integrity_output_file)
log_message("Session:   ", session_info_file)
