#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 2E: revision of three manuscript figures only
# =============================================================================
#
# PURPOSE
#   Re-render only the following three figures from frozen Round 2B outputs:
#     1. Figure 1: revised GTEx mapping annotation
#     2. Figure 3B: two-sided matched-background FDR display
#     3. Supplementary Figure S1B: explicit representative-cell rule
#
# ANALYTICAL BOUNDARY
#   - No differential-expression model is fitted.
#   - No tissue-identity, matched-background selection, recurrence-null,
#     baseline, coexpression, HPA, survival, clinical, subgroup, immune,
#     purity, pathway, or LIHC-coexpression analysis is run.
#   - Frozen values are read unchanged. The only permitted derived statistic is
#     Benjamini-Hochberg adjustment of empirical_p_two_sided across the 39
#     estimable matched-background cells, as specified in the manuscript.
#   - No missing value is imputed.
#
# OUTPUTS
#   results/round2E/figures/
#     Figure1_round2E_GTEx_v8_normal_tissue_atlas.{pdf,png}
#     Figure3B_round2E_matched_background_context.{pdf,png}
#     SupplementaryFigureS1B_round2E_representative_coexpression_rank_scatter.{pdf,png}
#
#   results/round2E/logs/
#     round2E_progress.log
#     round2E_input_manifest.csv
#     round2E_integrity_checks.csv
#     round2E_completion_report.txt
#     round2E_sessionInfo.txt
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration --------------------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round2E_figures"
required_r_version <- "4.6.1"

target_genes <- c("MAOB", "CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
mapped_tcga_order <- c(
  "BRCA", "COAD", "ESCA", "KICH", "KIRC", "KIRP",
  "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA"
)

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
  mustWork = TRUE
)

processed_dir <- file.path(project_dir, "data", "processed")
round2b_table_dir <- file.path(project_dir, "results", "round2B", "tables")
round2e_dir <- file.path(project_dir, "results", "round2E")
figure_dir <- file.path(round2e_dir, "figures")
log_dir <- file.path(round2e_dir, "logs")

invisible(lapply(
  c(round2e_dir, figure_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

input_files <- c(
  gtex_mapping = file.path(
    round2b_table_dir,
    "140_GTEx_v8_mapping_qc.csv"
  ),
  gtex_specificity = file.path(
    round2b_table_dir,
    "142_GTEx_v8_target_tissue_specificity.csv"
  ),
  matched_background = file.path(
    round2b_table_dir,
    "160_round2B_matched_background_cell_results.csv"
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
  gtex_reference_cache = file.path(
    processed_dir,
    "round2_GTEx_v8_gene_metrics.rds"
  )
)

progress_log_file <- file.path(log_dir, "round2E_progress.log")
input_manifest_file <- file.path(log_dir, "round2E_input_manifest.csv")
integrity_file <- file.path(log_dir, "round2E_integrity_checks.csv")
session_info_file <- file.path(log_dir, "round2E_sessionInfo.txt")
completion_report_file <- file.path(log_dir, "round2E_completion_report.txt")

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr",
  "tibble", "ggplot2", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 2E:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 2E requires R ", required_r_version,
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

# ---- 1. Utilities ------------------------------------------------------------

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
    tmpdir = dirname(path),
    fileext = ".csv"
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(x, tmp, na = "")
  replace_with_temp_file(tmp, path)
}

write_lines_atomic <- function(x, path) {
  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".txt"
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

safe_md5 <- function(path) unname(tools::md5sum(path)[[1]])

missing_inputs <- unname(input_files)[!file.exists(unname(input_files))]
if (length(missing_inputs) > 0L) {
  stop(
    "Required frozen Round 2B inputs are missing:\n  ",
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
    "Input ", label, ": ", rows,
    " rows/dimension; ", details
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
      paste0(
        label, " must contain ", expected_rows,
        " rows; found ", nrow(x), "."
      )
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

draw_single_plot <- function(plot) print(plot)

draw_plot_grid <- function(plots, nrow, ncol, widths = NULL, heights = NULL) {
  if (is.null(widths)) widths <- rep(1, ncol)
  if (is.null(heights)) heights <- rep(1, nrow)
  assert_true(
    length(plots) == nrow * ncol,
    "Plot-grid dimensions do not match the number of plots."
  )
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = nrow,
    ncol = ncol,
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
    tmpdir = dirname(path),
    fileext = ext
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  device_open <- FALSE
  tryCatch({
    if (type == "pdf") {
      grDevices::pdf(
        file = tmp,
        width = width,
        height = height,
        onefile = FALSE,
        useDingbats = FALSE,
        family = "Helvetica"
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
      plots = plots,
      nrow = nrow,
      ncol = ncol,
      widths = widths,
      heights = heights
    )
  }
  render_atomic(pdf_path, "pdf", width, height, dpi, draw_fun)
  render_atomic(png_path, "png", width, height, dpi, draw_fun)
  log_message("Wrote composite figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

integrity_checks <- list()
record_check <- function(check, result, expected, passed, details = "") {
  integrity_checks[[length(integrity_checks) + 1L]] <<- tibble::tibble(
    check = check,
    result = as.character(result),
    expected = as.character(expected),
    passed = isTRUE(passed),
    details = details
  )
  assert_true(isTRUE(passed), paste0("Integrity check failed: ", check))
  invisible(NULL)
}

# ---- 2. Read frozen inputs ---------------------------------------------------

gtex_mapping <- read_frozen_csv(
  "gtex_mapping",
  c(
    "tcga_code", "organ", "mapping_valid_for_baseline_model",
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

matched_background <- read_frozen_csv(
  "matched_background",
  c(
    "tcga_code", "gene", "target_logFC", "background_median_logFC",
    "background_percentile", "empirical_p_two_sided", "matching_status"
  ),
  expected_rows = 60L
)

coexpression_cells <- read_frozen_csv(
  "coexpression_cells",
  c(
    "tcga_code", "partner", "delta_rho", "status",
    "q_value_primary", "passes_coexpression_effect", "robust_rewiring"
  ),
  expected_rows = 52L
)

coexpression_representatives <- read_frozen_csv(
  "coexpression_representatives",
  c(
    "tcga_code", "partner", "rho_normal", "rho_tumor", "delta_rho"
  ),
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

# ---- 3. Figure 1: GTEx v8 normal-tissue atlas -------------------------------

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
  is.matrix(gtex_tpm_matrix) &&
    ncol(gtex_tpm_matrix) == 54L &&
    length(gtex_tissues) == 54L &&
    identical(colnames(gtex_tpm_matrix), gtex_tissues),
  "Figure 1 requires exactly 54 consistently labelled GTEx v8 tissues."
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

valid_mapping <- gtex_mapping %>%
  filter(mapping_valid_for_baseline_model, !is.na(matched_gtex_tissues))

mapped_gtex_tissues <- valid_mapping %>%
  select(matched_gtex_tissues) %>%
  separate_rows(matched_gtex_tissues, sep = ";") %>%
  transmute(tissue = str_trim(matched_gtex_tissues)) %>%
  distinct() %>%
  pull(tissue)

n_mapped_organ_baselines <- valid_mapping %>%
  distinct(organ) %>%
  nrow()

assert_true(
  length(mapped_gtex_tissues) == 12L &&
    n_mapped_organ_baselines == 9L &&
    all(mapped_gtex_tissues %in% gtex_tissues),
  paste0(
    "Expected 12 GTEx categories contributing to nine mapped organ baselines; ",
    "found ", length(mapped_gtex_tissues), " categories and ",
    n_mapped_organ_baselines, " organ baselines."
  )
)

record_check(
  "figure1_mapping_annotation",
  paste0(length(mapped_gtex_tissues), " categories / ",
         n_mapped_organ_baselines, " organ baselines"),
  "12 categories / 9 organ baselines",
  length(mapped_gtex_tissues) == 12L && n_mapped_organ_baselines == 9L,
  "Asterisks denote tissue categories contributing to mapped baselines."
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
    option = "C",
    begin = 0.05,
    end = 0.95,
    name = "log2(TPM + 1)"
  ) +
  scale_x_discrete(labels = tissue_axis_labels, drop = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    title = "Normal-tissue expression architecture",
    subtitle = paste(
      "Official GTEx v8 tissue-median TPM;",
      paste0(
        "* denotes a tissue category contributing to one of the ",
        n_mapped_organ_baselines,
        " TCGA-mapped organ baselines"
      )
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
    size = 3.1,
    show.legend = FALSE
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
  stem = "Figure1_round2E_GTEx_v8_normal_tissue_atlas",
  nrow = 1L,
  ncol = 2L,
  widths = c(4.3, 1.35),
  heights = 1,
  width = 18.5,
  height = 8.0
)

# ---- 4. Figure 3B: matched background with two-sided FDR --------------------

matched_background <- matched_background %>%
  mutate(
    estimable = str_detect(matching_status, "^estimable") &
      is.finite(target_logFC) &
      is.finite(empirical_p_two_sided),
    gene = factor(gene, levels = target_genes),
    tcga_code = factor(tcga_code, levels = rev(mapped_tcga_order))
  )

matched_estimable <- matched_background %>% filter(estimable)
matched_not_estimable <- matched_background %>% filter(!estimable)

assert_true(
  nrow(matched_estimable) == 39L && nrow(matched_not_estimable) == 21L,
  "Matched-background figure requires 39 estimable and 21 non-estimable cells."
)

# This is the only permitted statistical derivation in Round 2E.
matched_estimable <- matched_estimable %>%
  mutate(
    q_empirical_two_sided = p.adjust(
      empirical_p_two_sided,
      method = "BH"
    )
  )

minimum_two_sided_fdr <- min(
  matched_estimable$q_empirical_two_sided,
  na.rm = TRUE
)

raw_minimum_cell <- matched_estimable %>%
  arrange(empirical_p_two_sided, gene, tcga_code) %>%
  slice(1L)

expected_minimum_q <- 0.1724709784411277
expected_minimum_p <- 0.005994005994005994
n_minimum_q_ties <- sum(
  abs(matched_estimable$q_empirical_two_sided - minimum_two_sided_fdr) < 1e-12
)

assert_true(
  abs(minimum_two_sided_fdr - expected_minimum_q) < 1e-10 &&
    !any(matched_estimable$q_empirical_two_sided < 0.05),
  "Two-sided matched-background FDR does not match the frozen manuscript result."
)
assert_true(
  as.character(raw_minimum_cell$tcga_code) == "KICH" &&
    as.character(raw_minimum_cell$gene) == "MAOB" &&
    abs(raw_minimum_cell$empirical_p_two_sided - expected_minimum_p) < 1e-12,
  "KICH-MAOB is not the expected smallest raw two-sided empirical P value."
)
assert_true(
  n_minimum_q_ties > 1L,
  "Expected multiple cells to share the minimum BH-adjusted two-sided FDR."
)

matched_estimable <- matched_estimable %>%
  mutate(
    two_sided_fdr_band = case_when(
      q_empirical_two_sided < 0.25 ~ "q < 0.25",
      q_empirical_two_sided < 0.50 ~ "0.25 <= q < 0.50",
      TRUE ~ "q >= 0.50"
    ),
    two_sided_fdr_band = factor(
      two_sided_fdr_band,
      levels = c("q < 0.25", "0.25 <= q < 0.50", "q >= 0.50")
    )
  )

fdr_band_counts <- matched_estimable %>%
  count(two_sided_fdr_band, .drop = FALSE)

record_check(
  "figure3B_two_sided_BH_family",
  nrow(matched_estimable),
  39L,
  nrow(matched_estimable) == 39L,
  "BH correction was applied only across the 39 estimable cells."
)
record_check(
  "figure3B_minimum_two_sided_FDR",
  formatC(minimum_two_sided_fdr, digits = 10, format = "f"),
  formatC(expected_minimum_q, digits = 10, format = "f"),
  abs(minimum_two_sided_fdr - expected_minimum_q) < 1e-10,
  paste0(
    "Smallest raw P: KICH-MAOB, p=", 
    formatC(expected_minimum_p, digits = 6, format = "f"),
    "; minimum q is shared by ", n_minimum_q_ties, " cells."
  )
)
record_check(
  "figure3B_two_sided_FDR_bands_nonempty",
  paste(fdr_band_counts$n, collapse = "/"),
  "all three bands nonempty",
  all(fdr_band_counts$n > 0L),
  paste0(
    "Bands: ",
    paste(
      paste0(as.character(fdr_band_counts$two_sided_fdr_band), "=",
             fdr_band_counts$n),
      collapse = "; "
    )
  )
)

plot_matched <- ggplot(
  matched_background,
  aes(x = target_logFC, y = tcga_code)
) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_segment(
    data = matched_estimable,
    aes(
      x = background_median_logFC,
      xend = target_logFC,
      y = tcga_code,
      yend = tcga_code
    ),
    inherit.aes = FALSE,
    color = "grey65",
    linewidth = 0.45
  ) +
  geom_point(
    data = matched_estimable,
    aes(
      x = target_logFC,
      y = tcga_code,
      fill = background_percentile,
      shape = two_sided_fdr_band
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 2.7,
    stroke = 0.45
  ) +
  geom_label(
    data = matched_not_estimable,
    aes(x = 0, y = tcga_code, label = "NE"),
    inherit.aes = FALSE,
    fill = "grey90",
    color = "grey35",
    linewidth = 0.2,
    label.padding = grid::unit(0.08, "lines"),
    size = 2.25
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "grey80",
    high = "#B2182B",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "Matched-background\npercentile"
  ) +
  scale_shape_manual(
    values = c(
      "q < 0.25" = 24,
      "0.25 <= q < 0.50" = 22,
      "q >= 0.50" = 21
    ),
    drop = FALSE,
    name = "Two-sided FDR"
  ) +
  facet_wrap(~ gene, ncol = 3, drop = FALSE) +
  labs(
    x = "Paired limma logFC (tumor - adjacent non-tumor)",
    y = NULL,
    title = "Cancer-level target differences within expression-matched backgrounds",
    subtitle = paste0(
      "Segments begin at the matched-background median; minimum two-sided ",
      "FDR = 0.1725; no cell had FDR < 0.05.\n",
      "NE is a status label without a numeric x-position; BLCA and HNSC are ",
      "absent because no valid GTEx mapping was specified."
    )
  ) +
  theme_manuscript(base_size = 8.8) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(lineheight = 1.05)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barwidth = grid::unit(4.0, "cm")),
    shape = guide_legend(order = 2, nrow = 1)
  )

save_plot_pair(
  plot_matched,
  "Figure3B_round2E_matched_background_context",
  width = 13.5,
  height = 8.2
)

# ---- 5. Supplementary Figure S1B: explicit representative-cell rule --------

robust_cells <- coexpression_cells %>%
  filter(robust_rewiring, is.finite(delta_rho))

robust_cancers <- robust_cells %>%
  distinct(tcga_code) %>%
  arrange(tcga_code)

largest_abs_delta_by_cancer <- robust_cells %>%
  group_by(tcga_code) %>%
  arrange(desc(abs(delta_rho)), partner, .by_group = TRUE) %>%
  slice_head(n = 1L) %>%
  ungroup() %>%
  arrange(tcga_code) %>%
  select(tcga_code, partner, delta_rho)

representative_key <- coexpression_representatives %>%
  arrange(tcga_code) %>%
  select(tcga_code, partner, delta_rho)

selection_check <- largest_abs_delta_by_cancer %>%
  rename(expected_partner = partner, expected_delta_rho = delta_rho) %>%
  left_join(
    representative_key %>%
      rename(observed_partner = partner, observed_delta_rho = delta_rho),
    by = "tcga_code"
  )

assert_true(
  nrow(robust_cells) == 7L &&
    nrow(robust_cancers) == 3L &&
    setequal(robust_cancers$tcga_code, c("BRCA", "KIRC", "KIRP")),
  "Expected seven robust cells localized to BRCA, KIRC, and KIRP."
)
assert_true(
  nrow(selection_check) == 3L &&
    all(selection_check$expected_partner == selection_check$observed_partner) &&
    all(
      abs(
        selection_check$expected_delta_rho -
          selection_check$observed_delta_rho
      ) < 1e-12
    ),
  paste0(
    "Frozen representative cells do not follow the stated rule: largest ",
    "absolute delta rho in each cancer containing robust changes."
  )
)

point_keys <- coexpression_points %>%
  distinct(tcga_code, partner) %>%
  arrange(tcga_code, partner)
representative_keys <- coexpression_representatives %>%
  distinct(tcga_code, partner) %>%
  arrange(tcga_code, partner)

assert_true(
  identical(point_keys, representative_keys),
  "Representative coexpression cells and point table are inconsistent."
)
assert_true(
  !any(coexpression_cells$tcga_code == "LIHC") &&
    !any(coexpression_points$tcga_code == "LIHC"),
  "LIHC data must be absent from the coexpression inputs."
)

record_check(
  "supplementary_S1B_selection_rule",
  paste0(
    largest_abs_delta_by_cancer$tcga_code, "-",
    largest_abs_delta_by_cancer$partner,
    collapse = "; "
  ),
  "largest |delta rho| in each robust-change cancer",
  TRUE,
  paste0(
    "Seven robust cells occurred in three cancers; one representative was ",
    "selected per cancer by maximum absolute delta rho."
  )
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
    subtitle = paste0(
      "Rank-scale display; the cell with the largest |delta rho| in each of ",
      "the three cancers containing robust changes"
    )
  ) +
  theme_manuscript(base_size = 9) +
  theme(legend.position = "none")

save_plot_pair(
  plot_coexpression_representative,
  "SupplementaryFigureS1B_round2E_representative_coexpression_rank_scatter",
  width = 10.8,
  height = 6.8
)

# ---- 6. Final integrity, manifest, and session information ------------------

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
record_check(
  "frozen_inputs_unchanged",
  identical(input_md5_before, input_md5_after),
  TRUE,
  identical(input_md5_before, input_md5_after),
  paste0(length(input_files), " frozen input files checked by MD5.")
)

record_check(
  "main_analysis_modules_recomputed",
  FALSE,
  FALSE,
  TRUE,
  paste0(
    "Only three figures were rendered; BH adjustment of the existing ",
    "two-sided empirical P values was the sole permitted derivation."
  )
)

record_check(
  "forbidden_LIHC_or_clinical_modules_run",
  FALSE,
  FALSE,
  TRUE,
  paste0(
    "No LIHC coexpression, survival, clinical, subgroup, immune, purity, ",
    "pathway, or HNF4A analysis was run."
  )
)

integrity_table <- bind_rows(integrity_checks)
write_csv_atomic(integrity_table, integrity_file)

expected_figure_stems <- c(
  "Figure1_round2E_GTEx_v8_normal_tissue_atlas",
  "Figure3B_round2E_matched_background_context",
  "SupplementaryFigureS1B_round2E_representative_coexpression_rank_scatter"
)
expected_figure_files <- unlist(lapply(
  expected_figure_stems,
  function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".png")))
))

assert_true(
  length(expected_figure_files) == 6L &&
    all(file.exists(expected_figure_files)) &&
    all(file.info(expected_figure_files)$size > 1000),
  "One or more Round 2E figure files are missing or unexpectedly small."
)

write_lines_atomic(capture.output(sessionInfo()), session_info_file)

completion_lines <- c(
  paste0("Analysis: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Project directory: ", project_dir),
  "Scope: three figure revisions from frozen Round 2B inputs only",
  paste0(
    "Figure 1: ", length(mapped_gtex_tissues),
    " marked GTEx categories contributing to ",
    n_mapped_organ_baselines, " TCGA-mapped organ baselines"
  ),
  paste0(
    "Figure 3B: BH across 39 two-sided empirical P values; minimum FDR=",
    formatC(minimum_two_sided_fdr, digits = 6, format = "f"),
    "; no FDR<0.05; minimum shared by ", n_minimum_q_ties, " cells"
  ),
  paste0(
    "Supplementary Figure S1B: largest |delta rho| cell selected in each ",
    "of the three cancers containing seven robust changes"
  ),
  paste0("Input files unchanged by MD5: ", length(input_files)),
  paste0("Figure files written: ", length(expected_figure_files)),
  paste0("Figures: ", figure_dir),
  paste0("Logs: ", log_dir)
)
write_lines_atomic(completion_lines, completion_report_file)

log_message("Round 2E complete.")
log_message("Figures: ", figure_dir)
log_message("Manifest: ", input_manifest_file)
log_message("Integrity: ", integrity_file)
log_message("Session: ", session_info_file)
log_message("Completion report: ", completion_report_file)

