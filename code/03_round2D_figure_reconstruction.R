#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 2D: manuscript figure reconstruction only
# =============================================================================
#
# PURPOSE
#   Reconstruct manuscript figures from frozen Round 2B / Round 2C outputs.
#   No differential-expression model, empirical-null model, tissue-identity
#   definition, matched-background selection, baseline model, coexpression
#   statistic, or HPA analysis is fitted or recomputed here.
#
# IMPORTANT INPUT EXCEPTION FOR FIGURE 1A
#   The 5-gene x 54-tissue values are not present in Tables 141/142. Figure 1A
#   therefore reads the existing Round 2B cache
#     data/processed/round2_GTEx_v8_gene_metrics.rds
#   which contains the official GTEx v8 tissue-median TPM matrix used by Round
#   2B. The cache is read only. If it is absent, this script stops rather than
#   downloading or recreating it.
#
# FIGURES WRITTEN
#   Figure1_round2D_GTEx_v8_normal_tissue_atlas.{pdf,png}
#   Figure2_round2D_count_level_paired_heatmap.{pdf,png}
#   Figure3A_round2D_tissue_identity_distribution.{pdf,png}
#   Figure3B_round2D_matched_background_context.{pdf,png}
#   Figure4A_round2D_same_gene_recurrence_null.{pdf,png}
#   Figure4B_round2D_coverage_sensitivity.{pdf,png}
#   SupplementaryFigureS1A_round2D_nonLIHC_coexpression_change_heatmap.{pdf,png}
#   SupplementaryFigureS1B_round2D_representative_coexpression_rank_scatter.{pdf,png}
#
# EXPLICITLY OUT OF SCOPE
#   - Any refit or recalculation of the Round 2B / Round 2C analyses
#   - LIHC coexpression
#   - Survival, clinical, subgroup, HNF4A, immune, purity, or pathway analyses
#   - Missing-value imputation
#   - Reproduction of unchanged Supplementary Figures S2/S3
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration --------------------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round2D_figures"
required_r_version <- "4.6.1"

target_genes <- c("MAOB", "CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
cyp_partners <- c("CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
tcga_order <- c(
  "BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC",
  "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA"
)
mapped_tcga_order <- c(
  "BRCA", "COAD", "ESCA", "KICH", "KIRC", "KIRP",
  "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA"
)
coexpression_tcga_order <- setdiff(tcga_order, "LIHC")

gene_palette <- c(
  MAOB = "#222222",
  CYP3A4 = "#D55E00",
  CYP3A5 = "#0072B2",
  CYP3A7 = "#009E73",
  CYP3A43 = "#CC79A7"
)

project_dir <- normalizePath(
  Sys.getenv("MAOB_CYP3A_PROJECT_DIR", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)

processed_dir <- file.path(project_dir, "data", "processed")
round2b_dir <- file.path(project_dir, "results", "round2B")
round2b_table_dir <- file.path(round2b_dir, "tables")
round2c_table_dir <- file.path(project_dir, "results", "round2C", "tables")
round2d_dir <- file.path(project_dir, "results", "round2D")
figure_dir <- file.path(round2d_dir, "figures")
log_dir <- file.path(round2d_dir, "logs")

invisible(lapply(
  c(round2d_dir, figure_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

input_files <- c(
  target_results = file.path(
    round2b_table_dir,
    "121_round2B_target_paired_limma_voom_results.csv"
  ),
  gtex_mapping = file.path(
    round2b_table_dir,
    "140_GTEx_v8_mapping_qc.csv"
  ),
  gtex_specificity = file.path(
    round2b_table_dir,
    "142_GTEx_v8_target_tissue_specificity.csv"
  ),
  identity_reference = file.path(
    round2b_table_dir,
    "151_round2B_tissue_identity_reference_genes.csv"
  ),
  identity_targets = file.path(
    round2b_table_dir,
    "152_round2B_target_residual_vs_tissue_identity.csv"
  ),
  matched_background = file.path(
    round2b_table_dir,
    "160_round2B_matched_background_cell_results.csv"
  ),
  recurrence_null = file.path(
    round2b_table_dir,
    "162_round2B_same_gene_recurrence_null.csv"
  ),
  recurrence_profiles = file.path(
    round2b_table_dir,
    "163_round2B_same_gene_null_gene_profiles.csv"
  ),
  recurrence_sensitivity = file.path(
    round2b_table_dir,
    "164_round2B_same_gene_recurrence_sensitivity.csv"
  ),
  coexpression_cells = file.path(
    round2b_table_dir,
    "180_round2B_nonLIHC_coexpression_rewiring.csv"
  ),
  coexpression_representatives = file.path(
    round2b_table_dir,
    "181_round2B_representative_coexpression_cells.csv"
  ),
  coexpression_points = file.path(
    round2b_table_dir,
    "182_round2B_representative_coexpression_points.csv"
  ),
  round2c_concordance = file.path(
    round2c_table_dir,
    "123_round2C_Xena_count_level_concordance_cells.csv"
  ),
  gtex_reference_cache = file.path(
    processed_dir,
    "round2_GTEx_v8_gene_metrics.rds"
  )
)

progress_log_file <- file.path(log_dir, "round2D_progress.log")
input_manifest_file <- file.path(log_dir, "round2D_input_manifest.csv")
session_info_file <- file.path(log_dir, "round2D_sessionInfo.txt")
completion_report_file <- file.path(log_dir, "round2D_completion_report.txt")

required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "tibble", "ggplot2", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 2D:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 2D requires R ", required_r_version,
    "; current version is ", as.character(getRversion()), ".",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

# ---- 1. Utilities and read-only input manifest ------------------------------

timestamp_utc <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

replace_with_temp_file <- function(tmp, path) {
  if (!file.exists(tmp)) {
    stop("Temporary output was not created: ", tmp, call. = FALSE)
  }
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not remove previous output: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Could not move completed temporary output to: ", path, call. = FALSE)
  }
  invisible(path)
}

write_csv_atomic <- function(x, path) {
  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path), fileext = ".csv"
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(x, tmp, na = "")
  replace_with_temp_file(tmp, path)
}

write_lines_atomic <- function(x, path) {
  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path), fileext = ".txt"
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

assert_has_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(
      label, " lacks required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

safe_md5 <- function(path) {
  unname(tools::md5sum(path)[[1]])
}

missing_inputs <- unname(input_files)[!file.exists(unname(input_files))]
if (length(missing_inputs) > 0L) {
  stop(
    "Required frozen Round 2B / 2C inputs are missing:\n  ",
    paste(missing_inputs, collapse = "\n  "),
    call. = FALSE
  )
}

input_md5_before <- vapply(input_files, safe_md5, character(1))
input_manifest_records <- list()

record_input <- function(label, path, rows, details = "") {
  input_manifest_records[[length(input_manifest_records) + 1L]] <<-
    tibble::tibble(
      input = label,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      rows_or_primary_dimension = as.integer(rows),
      md5 = safe_md5(path),
      details = details
    )
  log_message(
    "Input ", label, ": ", rows, " rows/dimension; ", details,
    ifelse(nzchar(details), "", "")
  )
  invisible(NULL)
}

read_frozen_csv <- function(label, required_columns, expected_rows = NULL) {
  path <- input_files[[label]]
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  assert_has_columns(x, required_columns, label)
  if (!is.null(expected_rows)) {
    assert_true(
      nrow(x) == expected_rows,
      paste0(label, " must contain ", expected_rows, " rows; found ", nrow(x), ".")
    )
  }
  record_input(label, path, nrow(x), paste0(ncol(x), " columns"))
  x
}

theme_manuscript <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey25"),
      plot.tag = element_text(face = "bold", size = rel(1.25)),
      panel.grid.minor = element_blank(),
      legend.key.height = grid::unit(0.45, "cm")
    )
}

draw_single_plot <- function(plot) {
  print(plot)
}

draw_plot_grid <- function(plots, nrow, ncol, widths = NULL, heights = NULL) {
  if (is.null(widths)) widths <- rep(1, ncol)
  if (is.null(heights)) heights <- rep(1, nrow)
  assert_true(
    length(plots) == nrow * ncol,
    "Plot-grid dimensions do not match the number of plots."
  )
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = nrow, ncol = ncol,
    widths = grid::unit(widths, "null"),
    heights = grid::unit(heights, "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(plots)) {
    row_i <- ((i - 1L) %/% ncol) + 1L
    col_i <- ((i - 1L) %% ncol) + 1L
    print(
      plots[[i]],
      newpage = FALSE,
      vp = grid::viewport(layout.pos.row = row_i, layout.pos.col = col_i)
    )
  }
  grid::popViewport()
  invisible(NULL)
}

render_atomic <- function(path, type, width, height, dpi, draw_function) {
  ext <- if (type == "pdf") ".pdf" else ".png"
  tmp <- tempfile(
    pattern = paste0(".", tools::file_path_sans_ext(basename(path)), "."),
    tmpdir = dirname(path), fileext = ext
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  device_open <- FALSE
  tryCatch({
    if (type == "pdf") {
      grDevices::pdf(
        file = tmp, width = width, height = height,
        onefile = FALSE, useDingbats = FALSE, family = "Helvetica"
      )
    } else {
      grDevices::png(
        filename = tmp,
        width = round(width * dpi),
        height = round(height * dpi),
        res = dpi,
        bg = "white"
      )
    }
    device_open <- TRUE
    draw_function()
    grDevices::dev.off()
    device_open <- FALSE
  }, finally = {
    if (device_open && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  })
  replace_with_temp_file(tmp, path)
}

save_plot_pair <- function(plot, stem, width, height, dpi = 320) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  render_atomic(
    pdf_path, "pdf", width, height, dpi,
    function() draw_single_plot(plot)
  )
  render_atomic(
    png_path, "png", width, height, dpi,
    function() draw_single_plot(plot)
  )
  log_message("Wrote figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

save_plot_grid_pair <- function(
    plots, stem, nrow, ncol, widths, heights, width, height, dpi = 320
) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  draw_fun <- function() {
    draw_plot_grid(
      plots = plots, nrow = nrow, ncol = ncol,
      widths = widths, heights = heights
    )
  }
  render_atomic(pdf_path, "pdf", width, height, dpi, draw_fun)
  render_atomic(png_path, "png", width, height, dpi, draw_fun)
  log_message("Wrote composite figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

# ---- 2. Read frozen tables and validate core contracts ----------------------

target_results <- read_frozen_csv(
  "target_results",
  c(
    "tcga_code", "gene", "n_pairs", "logFC", "analysis_status",
    "q_value_global", "robust_count_level_change",
    "count_level_interpretation"
  ),
  expected_rows = 70L
)

gtex_mapping <- read_frozen_csv(
  "gtex_mapping",
  c(
    "tcga_code", "mapping_valid_for_baseline_model",
    "matched_gtex_tissues"
  ),
  expected_rows = 14L
)

gtex_specificity <- read_frozen_csv(
  "gtex_specificity",
  c(
    "gene_id", "gene_symbol", "max_tpm", "cross_tissue_median_tpm",
    "tau", "breadth_0_1tpm", "breadth_1tpm"
  ),
  expected_rows = 5L
)

identity_reference <- read_frozen_csv(
  "identity_reference",
  c("tcga_code", "gene_id", "gene_symbol", "logFC"),
  expected_rows = 2400L
)

identity_targets <- read_frozen_csv(
  "identity_targets",
  c(
    "tcga_code", "gene", "logFC", "analysis_status",
    "identity_median_logFC", "n_identity_genes", "identity_status",
    "more_negative_than_identity"
  ),
  expected_rows = 60L
)

matched_background <- read_frozen_csv(
  "matched_background",
  c(
    "tcga_code", "gene", "target_logFC", "background_median_logFC",
    "background_percentile", "matching_status",
    "q_empirical_directional", "selective_disruption_candidate"
  ),
  expected_rows = 60L
)

recurrence_null <- read_frozen_csv(
  "recurrence_null",
  c(
    "gene", "minimum_coverage_fraction", "n_cancers",
    "observed_loss_fraction", "n_candidate_genes", "status",
    "q_recurrent_loss", "q_median_negative"
  ),
  expected_rows = 5L
)

recurrence_profiles <- read_frozen_csv(
  "recurrence_profiles",
  c(
    "target_gene", "loss_fraction", "coverage_fraction",
    "eligible_primary_null"
  )
)

recurrence_sensitivity <- read_frozen_csv(
  "recurrence_sensitivity",
  c(
    "gene", "minimum_coverage_fraction", "n_candidate_genes", "status",
    "q_recurrent_loss", "q_median_negative"
  ),
  expected_rows = 25L
)

coexpression_cells <- read_frozen_csv(
  "coexpression_cells",
  c(
    "tcga_code", "partner", "delta_rho", "status", "q_value_primary",
    "passes_coexpression_effect", "robust_rewiring"
  ),
  expected_rows = 52L
)

coexpression_representatives <- read_frozen_csv(
  "coexpression_representatives",
  c("tcga_code", "partner", "rho_normal", "rho_tumor", "delta_rho"),
  expected_rows = 3L
)

coexpression_points <- read_frozen_csv(
  "coexpression_points",
  c(
    "tcga_code", "partner", "case_id", "condition",
    "MAOB_rank_fraction", "partner_rank_fraction", "panel"
  ),
  expected_rows = 434L
)

round2c_concordance <- read_frozen_csv(
  "round2c_concordance",
  c(
    "tcga_code", "gene", "round2B_logFC", "stage1_shift",
    "finite_overlap", "status"
  ),
  expected_rows = 70L
)

assert_true(
  !anyDuplicated(target_results[c("tcga_code", "gene")]),
  "Target table contains duplicate cancer-gene cells."
)
assert_true(
  setequal(target_results$tcga_code, tcga_order) &&
    setequal(target_results$gene, target_genes),
  "Target table does not contain the expected 14-cancer x 5-gene grid."
)
assert_true(
  !any(coexpression_cells$tcga_code == "LIHC") &&
    !any(coexpression_points$tcga_code == "LIHC"),
  "Paper-boundary failure: LIHC coexpression data are present."
)

lihc_count_cells <- target_results %>%
  filter(tcga_code == "LIHC", gene %in% c("MAOB", "CYP3A4"))
lihc_round2c_cells <- round2c_concordance %>%
  filter(tcga_code == "LIHC", gene %in% c("MAOB", "CYP3A4"))
assert_true(
  nrow(lihc_count_cells) == 2L &&
    all(is.finite(lihc_count_cells$logFC)) &&
    !any(lihc_count_cells$analysis_status == "reserved_not_recomputed"),
  "Figure 2 contract failed: LIHC-MAOB/CYP3A4 are not both finite."
)
assert_true(
  nrow(lihc_round2c_cells) == 2L &&
    all(lihc_round2c_cells$finite_overlap),
  "Round 2C contract failed: the two LIHC audit cells are not restored."
)

# ---- 3. Figure 1: official GTEx v8 atlas and specificity --------------------

gtex_reference <- readRDS(input_files[["gtex_reference_cache"]])
assert_true(
  is.list(gtex_reference) &&
    all(c("gene_metrics", "tpm_matrix", "tissue_columns") %in%
          names(gtex_reference)),
  "GTEx reference cache has an unexpected structure."
)

gtex_metrics_cache <- tibble::as_tibble(gtex_reference$gene_metrics)
gtex_tpm_matrix <- gtex_reference$tpm_matrix
gtex_tissues <- as.character(gtex_reference$tissue_columns)
assert_true(
  is.matrix(gtex_tpm_matrix) && ncol(gtex_tpm_matrix) == 54L &&
    length(gtex_tissues) == 54L,
  "Figure 1A requires exactly 54 official GTEx v8 tissue categories."
)
assert_true(
  identical(colnames(gtex_tpm_matrix), gtex_tissues),
  "GTEx cache tissue-column labels are inconsistent."
)
record_input(
  "gtex_reference_cache",
  input_files[["gtex_reference_cache"]],
  nrow(gtex_tpm_matrix),
  paste0(ncol(gtex_tpm_matrix), " tissue categories in cached TPM matrix")
)

gtex_target_metrics_cache <- gtex_metrics_cache %>%
  filter(gene_symbol %in% target_genes) %>%
  arrange(match(gene_symbol, target_genes))
assert_true(
  nrow(gtex_target_metrics_cache) == 5L &&
    identical(gtex_target_metrics_cache$gene_symbol, target_genes),
  "The GTEx cache does not contain one row for each target gene."
)

metric_check <- gtex_specificity %>%
  arrange(match(gene_symbol, target_genes)) %>%
  select(
    gene_symbol, max_tpm, cross_tissue_median_tpm,
    tau, breadth_0_1tpm, breadth_1tpm
  ) %>%
  left_join(
    gtex_target_metrics_cache %>%
      select(
        gene_symbol,
        max_tpm_cache = max_tpm,
        cross_tissue_median_tpm_cache = cross_tissue_median_tpm,
        tau_cache = tau,
        breadth_0_1tpm_cache = breadth_0_1tpm,
        breadth_1tpm_cache = breadth_1tpm
      ),
    by = "gene_symbol"
  )
metric_differences <- c(
  abs(metric_check$max_tpm - metric_check$max_tpm_cache),
  abs(
    metric_check$cross_tissue_median_tpm -
      metric_check$cross_tissue_median_tpm_cache
  ),
  abs(metric_check$tau - metric_check$tau_cache),
  abs(metric_check$breadth_0_1tpm - metric_check$breadth_0_1tpm_cache),
  abs(metric_check$breadth_1tpm - metric_check$breadth_1tpm_cache)
)
assert_true(
  all(is.finite(metric_differences)) && max(metric_differences) < 1e-10,
  "GTEx cache metrics do not exactly agree with frozen Table 142."
)

target_matrix_rows <- match(
  gtex_target_metrics_cache$gene_id,
  rownames(gtex_tpm_matrix)
)
assert_true(
  all(is.finite(target_matrix_rows)),
  "One or more target Ensembl IDs are absent from the GTEx TPM matrix."
)

gtex_target_tpm <- gtex_tpm_matrix[
  target_matrix_rows, gtex_tissues, drop = FALSE
]
rownames(gtex_target_tpm) <- gtex_target_metrics_cache$gene_symbol

atlas_long <- as.data.frame(gtex_target_tpm, check.names = FALSE) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "tissue",
    values_to = "tpm"
  ) %>%
  mutate(log2_tpm1 = log2(tpm + 1))

mapped_gtex_tissues <- gtex_mapping %>%
  filter(mapping_valid_for_baseline_model, !is.na(matched_gtex_tissues)) %>%
  select(matched_gtex_tissues) %>%
  separate_rows(matched_gtex_tissues, sep = ";") %>%
  transmute(tissue = str_trim(matched_gtex_tissues)) %>%
  distinct() %>%
  pull(tissue)
assert_true(
  length(mapped_gtex_tissues) == 12L &&
    all(mapped_gtex_tissues %in% gtex_tissues),
  paste0(
    "Expected 12 GTEx tissue categories used by valid mappings; found ",
    length(mapped_gtex_tissues), "."
  )
)

gastrointestinal_anchor <- c(
  "Liver",
  "Small Intestine - Terminal Ileum",
  "Colon - Transverse",
  "Colon - Sigmoid",
  "Stomach",
  "Esophagus - Gastroesophageal Junction",
  "Esophagus - Mucosa",
  "Esophagus - Muscularis",
  "Pancreas"
)

tissue_order_table <- tibble::tibble(tissue = gtex_tissues) %>%
  mutate(
    anchor_rank = match(tissue, gastrointestinal_anchor),
    tissue_group = case_when(
      !is.na(anchor_rank) ~ 1L,
      str_detect(tissue, "^(Kidney|Bladder|Adrenal)") ~ 2L,
      str_detect(tissue, "^(Lung|Minor Salivary)") ~ 3L,
      str_detect(tissue, "^(Breast|Adipose|Skin)") ~ 4L,
      str_detect(
        tissue,
        "^(Prostate|Testis|Ovary|Uterus|Vagina|Cervix|Fallopian)"
      ) ~ 5L,
      str_detect(tissue, "^(Thyroid|Pituitary)") ~ 6L,
      str_detect(tissue, "^(Heart|Artery|Muscle|Nerve)") ~ 7L,
      str_detect(tissue, "^(Spleen|Whole Blood|Cells)") ~ 8L,
      str_detect(tissue, "^Brain") ~ 9L,
      TRUE ~ 10L
    ),
    within_group_rank = if_else(
      tissue_group == 1L,
      anchor_rank,
      as.integer(rank(tissue, ties.method = "first"))
    )
  ) %>%
  arrange(tissue_group, within_group_rank, tissue)
tissue_order <- tissue_order_table$tissue

atlas_long <- atlas_long %>%
  mutate(
    tissue = factor(tissue, levels = tissue_order),
    gene = factor(gene, levels = rev(target_genes))
  )
tissue_axis_labels <- setNames(
  ifelse(
    tissue_order %in% mapped_gtex_tissues,
    paste0(tissue_order, " *"),
    tissue_order
  ),
  tissue_order
)

plot_atlas <- ggplot(atlas_long, aes(x = tissue, y = gene, fill = log2_tpm1)) +
  geom_tile(color = "white", linewidth = 0.20) +
  scale_fill_viridis_c(
    option = "C", begin = 0.05, end = 0.95,
    name = "log2(TPM + 1)"
  ) +
  scale_x_discrete(labels = tissue_axis_labels, drop = FALSE) +
  labs(
    x = NULL, y = NULL,
    title = "Normal-tissue expression architecture",
    subtitle = paste(
      "Official GTEx v8 tissue-median TPM;",
      "* denotes a tissue category used in a valid TCGA-organ mapping"
    ),
    tag = "A"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(face = "italic"),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

specificity_plot_data <- gtex_specificity %>%
  mutate(
    gene_symbol = factor(gene_symbol, levels = target_genes),
    label_x = tau,
    label_y = case_when(
      gene_symbol == "MAOB" ~ breadth_1tpm - 0.065,
      gene_symbol == "CYP3A5" ~ breadth_1tpm + 0.060,
      gene_symbol == "CYP3A4" ~ breadth_1tpm + 0.060,
      gene_symbol == "CYP3A7" ~ breadth_1tpm - 0.035,
      gene_symbol == "CYP3A43" ~ breadth_1tpm + 0.065,
      TRUE ~ breadth_1tpm
    )
  )

plot_specificity <- ggplot(
  specificity_plot_data,
  aes(x = tau, y = breadth_1tpm, fill = gene_symbol)
) +
  geom_point(shape = 21, color = "black", size = 3.4, stroke = 0.5) +
  geom_text(
    aes(x = label_x, y = label_y, label = gene_symbol),
    size = 3.1, show.legend = FALSE
  ) +
  scale_fill_manual(values = gene_palette, guide = "none") +
  scale_x_continuous(
    limits = c(0.72, 1.01),
    breaks = c(0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
  ) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Tau tissue-specificity statistic",
    y = "Tissue breadth at >=1 TPM",
    title = "Specificity and expression breadth",
    subtitle = "Five prespecified target genes",
    tag = "B"
  ) +
  theme_manuscript(base_size = 9) +
  theme(aspect.ratio = 1)

save_plot_grid_pair(
  plots = list(plot_atlas, plot_specificity),
  stem = "Figure1_round2D_GTEx_v8_normal_tissue_atlas",
  nrow = 1L, ncol = 2L,
  widths = c(4.3, 1.35), heights = 1,
  width = 18.5, height = 8.0
)

# ---- 4. Figure 2: paired count-level heatmap -------------------------------

pair_labels <- target_results %>%
  distinct(tcga_code, n_pairs) %>%
  mutate(
    tcga_code = factor(tcga_code, levels = tcga_order),
    x_label = paste0(as.character(tcga_code), "\nn=", n_pairs)
  ) %>%
  arrange(tcga_code)
assert_true(nrow(pair_labels) == 14L, "Pair counts are not unique by cancer type.")

figure2_data <- target_results %>%
  left_join(
    pair_labels %>% select(tcga_code, x_label) %>%
      mutate(tcga_code = as.character(tcga_code)),
    by = "tcga_code"
  ) %>%
  mutate(
    label = case_when(
      analysis_status == "low_expression_not_estimable" ~ "L",
      robust_count_level_change ~ "*",
      TRUE ~ ""
    ),
    gene = factor(gene, levels = rev(target_genes)),
    x_label = factor(x_label, levels = pair_labels$x_label)
  )
assert_true(
  !any(figure2_data$label == "R"),
  "Figure 2 still contains a reserved-cell label."
)

plot_figure2 <- ggplot(
  figure2_data,
  aes(x = x_label, y = gene, fill = logFC)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 3.1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    na.value = "grey85", name = "limma logFC\n(tumor - normal)"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Paired tumor-adjacent non-tumor expression differences",
    subtitle = paste0(
      "* global target FDR<0.05 and |logFC|>=log2(1.5); ",
      "L=filtered for low expression"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

save_plot_pair(
  plot_figure2,
  "Figure2_round2D_count_level_paired_heatmap",
  width = 11, height = 5
)

# ---- 5. Figure 3A: tissue-identity distributions and targets ---------------

identity_reference_counts <- identity_reference %>%
  count(tcga_code, name = "n_identity")
assert_true(
  nrow(identity_reference_counts) == 12L &&
    setequal(identity_reference_counts$tcga_code, mapped_tcga_order) &&
    all(identity_reference_counts$n_identity == 200L),
  "Tissue-identity reference must contain 200 genes for each of 12 cancers."
)
assert_true(
  setequal(identity_targets$tcga_code, mapped_tcga_order) &&
    all(identity_targets$n_identity_genes == 200L),
  "Target tissue-identity table does not match the 12 x 200 design."
)

identity_target_estimable <- identity_targets %>%
  filter(is.finite(logFC))
n_identity_target_estimable <- nrow(identity_target_estimable)
n_more_negative <- sum(
  identity_target_estimable$more_negative_than_identity,
  na.rm = TRUE
)
assert_true(
  n_identity_target_estimable == 39L && n_more_negative == 12L,
  paste0(
    "Frozen identity table must contain 39 estimable target cells and 12 ",
    "below-median cells; found ", n_identity_target_estimable, " and ",
    n_more_negative, "."
  )
)

identity_ne_counts <- identity_targets %>%
  group_by(tcga_code) %>%
  summarise(
    n_estimable_targets = sum(is.finite(logFC)),
    n_not_estimable_targets = sum(!is.finite(logFC)),
    .groups = "drop"
  ) %>%
  mutate(
    x_label = paste0(tcga_code, "\nNE=", n_not_estimable_targets)
  )

identity_code_labels <- identity_ne_counts$x_label[
  match(mapped_tcga_order, identity_ne_counts$tcga_code)
]

identity_reference_plot <- identity_reference %>%
  left_join(
    identity_ne_counts %>% select(tcga_code, x_label),
    by = "tcga_code"
  ) %>%
  mutate(x_label = factor(x_label, levels = identity_code_labels))

identity_target_plot <- identity_target_estimable %>%
  left_join(
    identity_ne_counts %>% select(tcga_code, x_label),
    by = "tcga_code"
  ) %>%
  mutate(
    x_label = factor(x_label, levels = identity_code_labels),
    gene = factor(gene, levels = target_genes),
    relative_to_identity = if_else(
      more_negative_than_identity,
      "Below identity median",
      "At or above identity median"
    )
  )

plot_identity <- ggplot(
  identity_reference_plot,
  aes(x = x_label, y = logFC)
) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_violin(
    fill = "grey90", color = "grey55", scale = "width", trim = FALSE,
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.13, outlier.shape = NA, fill = "white",
    color = "grey25", linewidth = 0.45
  ) +
  geom_point(
    data = identity_target_plot,
    aes(
      x = x_label, y = logFC, fill = gene,
      shape = relative_to_identity, group = gene
    ),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.68),
    color = "black", size = 2.5, stroke = 0.45
  ) +
  scale_fill_manual(values = gene_palette, name = "Target gene") +
  scale_shape_manual(
    values = c(
      "Below identity median" = 25,
      "At or above identity median" = 21
    ),
    name = "Target position"
  ) +
  labs(
    x = NULL,
    y = "Paired limma logFC (tumor - adjacent non-tumor)",
    title = "Target differences within tissue-identity gene distributions",
    subtitle = paste0(
      "Grey distributions: 200 GTEx-defined identity genes per cancer; ",
      "39 target cells estimable, 12 below the identity median. ",
      "NE counts are shown beneath each cancer code."
    )
  ) +
  theme_manuscript(base_size = 9.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2)
  )

save_plot_pair(
  plot_identity,
  "Figure3A_round2D_tissue_identity_distribution",
  width = 12.5, height = 7.0
)

# ---- 6. Figure 3B: matched background with directional FDR -----------------

matched_background <- matched_background %>%
  mutate(
    estimable = str_detect(matching_status, "^estimable") &
      is.finite(target_logFC),
    gene = factor(gene, levels = target_genes),
    tcga_code = factor(tcga_code, levels = rev(mapped_tcga_order))
  )
matched_estimable <- matched_background %>% filter(estimable)
matched_not_estimable <- matched_background %>% filter(!estimable)

assert_true(
  nrow(matched_estimable) == 39L && nrow(matched_not_estimable) == 21L,
  "Matched-background figure requires 39 estimable and 21 non-estimable cells."
)
minimum_directional_fdr <- min(
  matched_estimable$q_empirical_directional,
  na.rm = TRUE
)
assert_true(
  is.finite(minimum_directional_fdr) &&
    abs(minimum_directional_fdr - 0.0862354892205638) < 1e-10 &&
    !any(matched_estimable$q_empirical_directional < 0.05) &&
    !any(matched_estimable$selective_disruption_candidate),
  "Matched-background FDR contract does not match the frozen Results."
)

matched_estimable <- matched_estimable %>%
  mutate(
    directional_fdr_band = case_when(
      q_empirical_directional < 0.10 ~ "q < 0.10",
      q_empirical_directional < 0.25 ~ "0.10 <= q < 0.25",
      TRUE ~ "q >= 0.25"
    ),
    directional_fdr_band = factor(
      directional_fdr_band,
      levels = c("q < 0.10", "0.10 <= q < 0.25", "q >= 0.25")
    )
  )
minimum_fdr_cell <- matched_estimable %>%
  slice_min(q_empirical_directional, n = 1L, with_ties = FALSE)

plot_matched <- ggplot(
  matched_background,
  aes(x = target_logFC, y = tcga_code)
) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_segment(
    data = matched_estimable,
    aes(
      x = background_median_logFC, xend = target_logFC,
      y = tcga_code, yend = tcga_code
    ),
    inherit.aes = FALSE,
    color = "grey65", linewidth = 0.45
  ) +
  geom_point(
    data = matched_estimable,
    aes(
      x = target_logFC, y = tcga_code,
      fill = background_percentile, shape = directional_fdr_band
    ),
    inherit.aes = FALSE,
    color = "black", size = 2.7, stroke = 0.45
  ) +
  geom_label(
    data = matched_not_estimable,
    aes(x = 0, y = tcga_code, label = "NE"),
    inherit.aes = FALSE,
    fill = "grey90", color = "grey35", linewidth = 0.2,
    label.padding = grid::unit(0.08, "lines"), size = 2.25
  ) +
  geom_text(
    data = minimum_fdr_cell,
    aes(
      x = target_logFC, y = tcga_code,
      label = paste0(
        "min q=", formatC(q_empirical_directional, digits = 4, format = "f")
      )
    ),
    inherit.aes = FALSE,
    nudge_x = 0.35, hjust = 0, size = 2.6
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "grey80", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    name = "Matched-background\npercentile"
  ) +
  scale_shape_manual(
    values = c(
      "q < 0.10" = 24,
      "0.10 <= q < 0.25" = 22,
      "q >= 0.25" = 21
    ),
    drop = FALSE,
    name = "Directional FDR"
  ) +
  facet_wrap(~ gene, ncol = 3, drop = FALSE) +
  labs(
    x = "Paired limma logFC (tumor - adjacent non-tumor)", y = NULL,
    title = "Cancer-level target differences within expression-matched backgrounds",
    subtitle = paste0(
      "Segments begin at the matched-background median; minimum directional ",
      "FDR=0.0862 and no cell had FDR<0.05. NE is a status label without a ",
      "numeric x-position. BLCA and HNSC are absent because no valid GTEx ",
      "mapping was specified."
    )
  ) +
  theme_manuscript(base_size = 8.8) +
  theme(legend.position = "bottom") +
  guides(
    fill = guide_colorbar(order = 1, barwidth = grid::unit(4.0, "cm")),
    shape = guide_legend(order = 2, nrow = 1)
  )

save_plot_pair(
  plot_matched,
  "Figure3B_round2D_matched_background_context",
  width = 13.0, height = 8.0
)

# ---- 7. Figure 4A: same-gene recurrence null -------------------------------

recurrence_evaluable_genes <- recurrence_null %>%
  filter(status == "estimable_same_gene_null") %>%
  pull(gene)
assert_true(
  setequal(recurrence_evaluable_genes, c("MAOB", "CYP3A4", "CYP3A7")),
  "Unexpected set of evaluable genes in primary recurrence null."
)

recurrence_null_plot_data <- recurrence_profiles %>%
  filter(
    eligible_primary_null,
    target_gene %in% recurrence_evaluable_genes
  ) %>%
  mutate(gene = factor(target_gene, levels = target_genes))

profile_counts <- recurrence_null_plot_data %>%
  count(gene, name = "profiles_in_plot") %>%
  mutate(gene = as.character(gene))
expected_profile_counts <- recurrence_null %>%
  filter(status == "estimable_same_gene_null") %>%
  transmute(gene, expected = as.integer(n_candidate_genes))
profile_count_check <- expected_profile_counts %>%
  left_join(profile_counts, by = "gene")
assert_true(
  all(profile_count_check$expected == profile_count_check$profiles_in_plot) &&
    identical(
      profile_count_check$expected[match(
        c("MAOB", "CYP3A4", "CYP3A7"), profile_count_check$gene
      )],
      c(66L, 92L, 313L)
    ),
  "Eligible primary recurrence-profile counts do not equal 66/92/313."
)

recurrence_observed_plot <- recurrence_null %>%
  filter(status == "estimable_same_gene_null") %>%
  mutate(gene = factor(gene, levels = target_genes))

candidate_labels <- recurrence_observed_plot %>%
  transmute(
    gene,
    y = -0.115,
    label = paste0("eligible profiles: ", as.integer(n_candidate_genes))
  )

not_evaluable_annotations <- recurrence_null %>%
  filter(status != "estimable_same_gene_null") %>%
  mutate(
    gene = factor(gene, levels = target_genes),
    y = 0.55,
    label = case_when(
      status == "insufficient_same_gene_candidates" ~ paste0(
        "Not evaluable\n", as.integer(n_candidate_genes),
        " candidates (<20)"
      ),
      status == "insufficient_cancers" ~ paste0(
        "Not evaluable\n", as.integer(n_cancers), " cancers"
      ),
      TRUE ~ "Not evaluable"
    )
  )

plot_recurrence <- ggplot(
  recurrence_null_plot_data,
  aes(x = gene, y = loss_fraction)
) +
  geom_violin(
    fill = "grey88", color = "grey45", scale = "width", trim = TRUE
  ) +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white") +
  geom_point(
    data = recurrence_observed_plot,
    aes(x = gene, y = observed_loss_fraction),
    inherit.aes = FALSE,
    shape = 23, size = 3.4, fill = "#B2182B", color = "black"
  ) +
  geom_text(
    data = recurrence_observed_plot,
    aes(
      x = gene,
      y = pmin(1.08, observed_loss_fraction + 0.08),
      label = paste0("q=", formatC(q_recurrent_loss, digits = 2, format = "g"))
    ),
    inherit.aes = FALSE, size = 3
  ) +
  geom_text(
    data = candidate_labels,
    aes(x = gene, y = y, label = label),
    inherit.aes = FALSE, size = 2.8, color = "grey25"
  ) +
  geom_label(
    data = not_evaluable_annotations,
    aes(x = gene, y = y, label = label),
    inherit.aes = FALSE,
    fill = "grey93", color = "grey30", linewidth = 0.25, size = 3
  ) +
  scale_x_discrete(drop = FALSE, limits = target_genes) +
  scale_y_continuous(
    limits = c(-0.18, 1.12),
    breaks = seq(0, 1, 0.25),
    labels = function(x) ifelse(x < 0, "", formatC(x, digits = 2, format = "f"))
  ) +
  labs(
    x = NULL,
    y = "Fraction of cancers with loss",
    title = "Cross-cancer recurrence against the same-gene empirical null",
    subtitle = paste0(
      "Grey distributions: eligible background-gene profiles at 80% ",
      "coverage; diamonds: observed target genes. Not evaluable is distinct ",
      "from evidence of no change."
    )
  ) +
  theme_manuscript(base_size = 10) +
  theme(axis.text.x = element_text(face = "italic"))

save_plot_pair(
  plot_recurrence,
  "Figure4A_round2D_same_gene_recurrence_null",
  width = 10.0, height = 5.8
)

# ---- 8. Figure 4B: coverage sensitivity ------------------------------------

assert_true(
  setequal(recurrence_sensitivity$gene, target_genes) &&
    setequal(
      recurrence_sensitivity$minimum_coverage_fraction,
      c(0.5, 0.6, 0.7, 0.8, 0.9)
    ) &&
    !anyDuplicated(
      recurrence_sensitivity[c("gene", "minimum_coverage_fraction")]
    ),
  "Coverage-sensitivity table is not the expected 5-gene x 5-coverage grid."
)

expected_not_evaluable <- tibble::tribble(
  ~gene, ~minimum_coverage_fraction, ~status,
  "MAOB", 0.9, "insufficient_same_gene_candidates",
  "CYP3A4", 0.9, "insufficient_same_gene_candidates",
  "CYP3A5", 0.8, "insufficient_same_gene_candidates",
  "CYP3A5", 0.9, "insufficient_same_gene_candidates",
  "CYP3A7", 0.9, "insufficient_same_gene_candidates",
  "CYP3A43", 0.5, "insufficient_cancers",
  "CYP3A43", 0.6, "insufficient_cancers",
  "CYP3A43", 0.7, "insufficient_cancers",
  "CYP3A43", 0.8, "insufficient_cancers",
  "CYP3A43", 0.9, "insufficient_cancers"
)
observed_not_evaluable <- recurrence_sensitivity %>%
  filter(status != "estimable_same_gene_null") %>%
  select(gene, minimum_coverage_fraction, status) %>%
  arrange(gene, minimum_coverage_fraction)
expected_not_evaluable <- expected_not_evaluable %>%
  arrange(gene, minimum_coverage_fraction)
assert_true(
  identical(observed_not_evaluable, expected_not_evaluable),
  "Coverage-sensitivity not-evaluable pattern differs from the frozen table."
)

sensitivity_evaluable <- recurrence_sensitivity %>%
  filter(
    status == "estimable_same_gene_null",
    is.finite(q_recurrent_loss)
  ) %>%
  mutate(gene = factor(gene, levels = target_genes))

not_evaluable_y <- c(
  MAOB = 1.55,
  CYP3A4 = 1.40,
  CYP3A5 = 1.27,
  CYP3A7 = 1.15,
  CYP3A43 = 1.03
)
sensitivity_not_evaluable <- recurrence_sensitivity %>%
  filter(status != "estimable_same_gene_null") %>%
  mutate(
    gene = factor(gene, levels = target_genes),
    display_y = unname(not_evaluable_y[as.character(gene)]),
    status_short = case_when(
      status == "insufficient_same_gene_candidates" ~
        "Insufficient candidate profiles",
      status == "insufficient_cancers" ~ "Insufficient cancers",
      TRUE ~ "Not evaluable"
    )
  )

maob_q <- sensitivity_evaluable %>%
  filter(gene == "MAOB") %>%
  arrange(minimum_coverage_fraction)
cyp3a4_q <- sensitivity_evaluable %>%
  filter(gene == "CYP3A4") %>%
  arrange(minimum_coverage_fraction)
assert_true(
  all(maob_q$minimum_coverage_fraction == c(0.5, 0.6, 0.7, 0.8)) &&
    all(maob_q$q_recurrent_loss < 0.05) &&
    sum(cyp3a4_q$q_recurrent_loss < 0.05) == 1L &&
    cyp3a4_q$minimum_coverage_fraction[
      which(cyp3a4_q$q_recurrent_loss < 0.05)
    ] == 0.8,
  "MAOB/CYP3A4 coverage-sensitivity pattern differs from the frozen Results."
)

plot_sensitivity <- ggplot(
  sensitivity_evaluable,
  aes(
    x = minimum_coverage_fraction,
    y = q_recurrent_loss,
    color = gene, group = gene
  )
) +
  annotate(
    "rect", xmin = 0.485, xmax = 0.915, ymin = 0.95, ymax = 1.72,
    fill = "grey95", color = NA
  ) +
  geom_hline(
    yintercept = 0.05, linetype = "dashed",
    color = "#B2182B", linewidth = 0.65
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.5) +
  geom_point(
    data = sensitivity_not_evaluable,
    aes(
      x = minimum_coverage_fraction,
      y = display_y,
      color = gene,
      shape = status_short
    ),
    inherit.aes = FALSE,
    size = 2.7, stroke = 0.9
  ) +
  annotate(
    "text", x = 0.70, y = 1.63,
    label = "Not evaluable status strip (marker height is not a q value)",
    size = 3.0, color = "grey25"
  ) +
  scale_color_manual(values = gene_palette, drop = FALSE, name = "Gene") +
  scale_shape_manual(
    values = c(
      "Insufficient candidate profiles" = 4,
      "Insufficient cancers" = 8
    ),
    name = "Not-evaluable reason"
  ) +
  scale_x_continuous(
    limits = c(0.485, 0.915),
    breaks = c(0.5, 0.6, 0.7, 0.8, 0.9),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_y_log10(
    limits = c(0.015, 1.72),
    breaks = c(0.02, 0.05, 0.10, 0.25, 0.50, 1.00),
    labels = c("0.02", "0.05", "0.10", "0.25", "0.50", "1.00")
  ) +
  labs(
    x = "Minimum coverage fraction",
    y = "Recurrent-loss FDR (log scale)",
    title = "Coverage sensitivity of the same-gene recurrence test",
    subtitle = paste0(
      "MAOB remained below FDR 0.05 from 50%-80% coverage; CYP3A4 was ",
      "below 0.05 only at 80%. All 25 gene-coverage cells are represented."
    )
  ) +
  theme_manuscript(base_size = 10) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(
    color = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2)
  )

save_plot_pair(
  plot_sensitivity,
  "Figure4B_round2D_coverage_sensitivity",
  width = 10.5, height = 6.3
)

# ---- 9. Supplementary Figure S1A: coexpression-change heatmap --------------

assert_true(
  setequal(coexpression_cells$tcga_code, coexpression_tcga_order) &&
    setequal(coexpression_cells$partner, cyp_partners) &&
    sum(coexpression_cells$status == "primary_n_ge_30") == 40L &&
    sum(coexpression_cells$status == "exploratory_n_20_29") == 4L &&
    sum(coexpression_cells$status == "insufficient_n_lt20") == 8L &&
    sum(coexpression_cells$robust_rewiring) == 7L,
  "Coexpression table does not match the frozen 52-cell status structure."
)

coexpression_plot_data <- coexpression_cells %>%
  mutate(
    partner = factor(partner, levels = rev(cyp_partners)),
    tcga_code = factor(tcga_code, levels = coexpression_tcga_order),
    label = case_when(
      robust_rewiring ~ "*",
      status == "exploratory_n_20_29" ~ "E",
      TRUE ~ ""
    )
  )

plot_coexpression_heatmap <- ggplot(
  coexpression_plot_data,
  aes(x = tcga_code, y = partner, fill = delta_rho)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(
    data = coexpression_plot_data %>%
      filter(status == "exploratory_n_20_29"),
    fill = "white", alpha = 0.48, color = "white", linewidth = 0.3
  ) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#762A83", mid = "white", high = "#1B7837", midpoint = 0,
    limits = c(-1, 1), oob = scales::squish, na.value = "grey85",
    name = "delta rho\n(tumor - normal)"
  ) +
  labs(
    x = NULL, y = "MAOB partner",
    title = "Within-patient MAOB-CYP3A coexpression change",
    subtitle = paste0(
      "Non-LIHC cancers only; *=primary FDR<0.05 and |delta rho|>=0.20; ",
      "E=exploratory n=20-29 (desaturated); grey=not estimable"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

save_plot_pair(
  plot_coexpression_heatmap,
  "SupplementaryFigureS1A_round2D_nonLIHC_coexpression_change_heatmap",
  width = 10.5, height = 4.7
)

# ---- 10. Supplementary Figure S1B: representative rank plots ---------------

assert_true(
  setequal(
    paste(
      coexpression_representatives$tcga_code,
      coexpression_representatives$partner,
      sep = "-"
    ),
    paste(
      unique(coexpression_points[c("tcga_code", "partner")])$tcga_code,
      unique(coexpression_points[c("tcga_code", "partner")])$partner,
      sep = "-"
    )
  ),
  "Representative coexpression cells and point table are inconsistent."
)

coexpression_points_plot <- coexpression_points %>%
  mutate(
    condition = factor(condition, levels = c("Normal", "Tumor")),
    panel = factor(panel, levels = unique(panel))
  )

plot_coexpression_representative <- ggplot(
  coexpression_points_plot,
  aes(
    x = MAOB_rank_fraction,
    y = partner_rank_fraction,
    color = condition
  )
) +
  geom_point(alpha = 0.65, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_grid(condition ~ panel) +
  scale_color_manual(values = c(Normal = "#2166AC", Tumor = "#B2182B")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Within-condition MAOB expression rank",
    y = "Within-condition CYP3A expression rank",
    title = "Representative non-LIHC coexpression changes",
    subtitle = "Rank-scale display corresponding to the frozen Spearman summaries"
  ) +
  theme_manuscript(base_size = 9) +
  theme(legend.position = "none")

save_plot_pair(
  plot_coexpression_representative,
  "SupplementaryFigureS1B_round2D_representative_coexpression_rank_scatter",
  width = 10.5, height = 6.8
)

# ---- 11. Final read-only checks, manifest, and session information ----------

input_manifest <- bind_rows(input_manifest_records)
assert_true(
  nrow(input_manifest) == length(input_files),
  paste0(
    "Input manifest should contain ", length(input_files),
    " records; found ", nrow(input_manifest), "."
  )
)
write_csv_atomic(input_manifest, input_manifest_file)

input_md5_after <- vapply(input_files, safe_md5, character(1))
assert_true(
  identical(input_md5_before, input_md5_after),
  "One or more frozen Round 2B / 2C inputs changed during Round 2D."
)

expected_figure_stems <- c(
  "Figure1_round2D_GTEx_v8_normal_tissue_atlas",
  "Figure2_round2D_count_level_paired_heatmap",
  "Figure3A_round2D_tissue_identity_distribution",
  "Figure3B_round2D_matched_background_context",
  "Figure4A_round2D_same_gene_recurrence_null",
  "Figure4B_round2D_coverage_sensitivity",
  "SupplementaryFigureS1A_round2D_nonLIHC_coexpression_change_heatmap",
  "SupplementaryFigureS1B_round2D_representative_coexpression_rank_scatter"
)
expected_figure_files <- unlist(lapply(
  expected_figure_stems,
  function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".png")))
))
assert_true(
  all(file.exists(expected_figure_files)) &&
    all(file.info(expected_figure_files)$size > 1000),
  "One or more Round 2D figure files are missing or unexpectedly small."
)

write_lines_atomic(capture.output(sessionInfo()), session_info_file)

completion_lines <- c(
  paste0("Analysis: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Project directory: ", project_dir),
  "Scope: figure reconstruction from frozen Round 2B / Round 2C values only",
  "Official GTEx v8 cache: read only; 54 tissue categories validated",
  "Figure 2: LIHC-MAOB and LIHC-CYP3A4 finite; no reserved legend",
  "Figure 3A: 2400 identity-gene rows; 39 target cells; 12 below median",
  paste0(
    "Figure 3B: minimum directional FDR=",
    formatC(minimum_directional_fdr, digits = 6, format = "f"),
    "; no FDR<0.05"
  ),
  "Figure 4A: eligible profiles MAOB=66, CYP3A4=92, CYP3A7=313",
  "Figure 4B: all 25 coverage-gene statuses represented",
  "Supplementary Figure S1: 52 non-LIHC cells; LIHC coexpression absent",
  paste0("Input files unchanged by MD5: ", length(input_files)),
  paste0("Figure files written: ", length(expected_figure_files)),
  paste0("Figures: ", figure_dir),
  paste0("Logs: ", log_dir)
)
write_lines_atomic(completion_lines, completion_report_file)

log_message("Round 2D complete.")
log_message("Figures: ", figure_dir)
log_message("Manifest: ", input_manifest_file)
log_message("Session: ", session_info_file)
log_message("Completion report: ", completion_report_file)
