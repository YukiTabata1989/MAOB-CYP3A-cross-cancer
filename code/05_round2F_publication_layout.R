#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 2F: publication-size typography and layout revision only
# =============================================================================
#
# PURPOSE
#   Re-render all eight manuscript figures at their final 174-mm publication
#   width. Scientific content and frozen values are unchanged from Round 2D,
#   with the three Round 2E corrections retained for Figure 1, Figure 3B, and
#   Supplementary Figure S1B. Figure 1A alone is transposed for legibility.
#   No analytical model or biological statistic is fitted or recomputed here.
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
#   Figure1_round2F_GTEx_v8_normal_tissue_atlas.{pdf,png}
#   Figure2_round2F_count_level_paired_heatmap.{pdf,png}
#   Figure3A_round2F_tissue_identity_distribution.{pdf,png}
#   Figure3B_round2F_matched_background_context.{pdf,png}
#   Figure4A_round2F_same_gene_recurrence_null.{pdf,png}
#   Figure4B_round2F_coverage_sensitivity.{pdf,png}
#   SupplementaryFigureS1A_round2F_nonLIHC_coexpression_change_heatmap.{pdf,png}
#   SupplementaryFigureS1B_round2F_representative_coexpression_rank_scatter.{pdf,png}
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

analysis_id <- "MAOB_CYP3A_tissue_disruption_round2F_figures"
required_r_version <- "4.6.1"
publication_width_mm <- 174
minimum_text_pt <- 7
publication_font_family <- Sys.getenv(
  "MAOB_CYP3A_FIGURE_FONT",
  unset = "Arial"
)

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

# The legacy palette is retained where Round 2F requested no colour change.
gene_palette <- c(
  MAOB = "#222222",
  CYP3A4 = "#D55E00",
  CYP3A5 = "#0072B2",
  CYP3A7 = "#009E73",
  CYP3A43 = "#CC79A7"
)

# Okabe-Ito categorical palette, used only for Figures 3A, 4A, and 4B.
gene_palette_okabe_ito <- c(
  MAOB = "#000000",
  CYP3A4 = "#E69F00",
  CYP3A5 = "#56B4E9",
  CYP3A7 = "#009E73",
  CYP3A43 = "#CC79A7"
)

project_dir <- normalizePath(
  Sys.getenv("MAOB_CYP3A_PROJECT_DIR", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

processed_dir <- file.path(project_dir, "data", "processed")
round2b_dir <- file.path(project_dir, "results", "round2B")
round2b_table_dir <- file.path(round2b_dir, "tables")
round2c_table_dir <- file.path(project_dir, "results", "round2C", "tables")
round2f_dir <- file.path(project_dir, "results", "round2F")
figure_dir <- file.path(round2f_dir, "figures")
log_dir <- file.path(round2f_dir, "logs")

invisible(lapply(
  c(round2f_dir, figure_dir, log_dir),
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

progress_log_file <- file.path(log_dir, "round2F_progress.log")
input_manifest_file <- file.path(log_dir, "round2F_input_manifest.csv")
text_audit_file <- file.path(log_dir, "round2F_text_size_audit.csv")
text_summary_file <- file.path(log_dir, "round2F_text_size_summary.csv")
figure_spec_file <- file.path(log_dir, "round2F_figure_dimensions.csv")
pdf_device_file <- file.path(log_dir, "round2F_pdf_device_manifest.csv")
integrity_file <- file.path(log_dir, "round2F_integrity_checks.csv")
session_info_file <- file.path(log_dir, "round2F_sessionInfo.txt")
completion_report_file <- file.path(log_dir, "round2F_completion_report.txt")

required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "tibble", "ggplot2", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 2F:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 2F requires R ", required_r_version,
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

wrap_plot_text <- function(x, width = 92L) {
  stringr::str_wrap(x, width = width)
}

# Text sizes below are points. Geom text sizes elsewhere are millimetres, as
# defined by ggplot2; these are also audited after the figures are written.
theme_manuscript <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = publication_font_family) +
    theme(
      text = element_text(size = base_size, family = publication_font_family),
      axis.text = element_text(size = 8.0),
      axis.title = element_text(size = 9.0),
      strip.text = element_text(size = 8.0),
      legend.text = element_text(size = 7.5),
      legend.title = element_text(size = 8.0),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 8.0, color = "grey25"),
      plot.caption = element_text(size = 7.5),
      plot.tag = element_text(face = "bold", size = 11.0),
      panel.grid.minor = element_blank(),
      legend.key.height = grid::unit(0.45, "cm")
    )
}

theme_heatmap <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = publication_font_family) +
    theme(
      text = element_text(size = base_size, family = publication_font_family),
      axis.text = element_text(size = 8.0),
      axis.title = element_text(size = 9.0),
      strip.text = element_text(size = 8.0),
      legend.text = element_text(size = 7.5),
      legend.title = element_text(size = 8.0),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 8.0, color = "grey25"),
      plot.caption = element_text(size = 7.5),
      plot.tag = element_text(face = "bold", size = 11.0),
      panel.grid = element_blank(),
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

pdf_device_records <- list()
text_audit_records <- list()
figure_output_records <- list()

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

attempt_graphics_device <- function(label, open_function) {
  device_before <- as.integer(grDevices::dev.cur())
  warning_messages <- character()
  error_message <- NULL

  opened <- tryCatch(
    withCallingHandlers({
      open_function()
      if (as.integer(grDevices::dev.cur()) == device_before) {
        stop(label, " returned without opening a graphics device.")
      }
      TRUE
    }, warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error_message <<- conditionMessage(e)
      FALSE
    }
  )

  # Close a partially opened device after a failed attempt.
  if (!opened && as.integer(grDevices::dev.cur()) != device_before) {
    try(grDevices::dev.off(), silent = TRUE)
  }

  list(
    opened = isTRUE(opened),
    error = error_message,
    warnings = warning_messages
  )
}

device_attempt_message <- function(x, default) {
  messages <- c(x$error, x$warnings)
  messages <- unique(messages[!is.na(messages) & nzchar(messages)])
  if (length(messages) == 0L) default else paste(messages, collapse = " | ")
}

open_embedded_pdf_device <- function(filename, width_in, height_in) {
  is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
  quartz_attempt <- list(opened = FALSE, error = NULL, warnings = character())
  cairo_attempt <- list(opened = FALSE, error = NULL, warnings = character())

  # On macOS, Quartz is native and avoids the failed-cairo-DLL warning seen in
  # some CRAN/RStudio installations. It is therefore attempted first.
  if (is_macos) {
    quartz_attempt <- attempt_graphics_device(
      "quartz_pdf",
      function() grDevices::quartz(
        type = "pdf",
        file = filename,
        width = width_in,
        height = height_in,
        family = publication_font_family,
        bg = "white"
      )
    )
    if (quartz_attempt$opened) return("quartz_pdf")
    if (file.exists(filename)) unlink(filename)
  }

  cairo_attempt <- attempt_graphics_device(
    "cairo_pdf",
    function() grDevices::cairo_pdf(
      filename = filename,
      width = width_in,
      height = height_in,
      # Each manuscript figure is a single-page PDF. onefile=FALSE can insert
      # a page number into a tempfile name and break the atomic exact-path test.
      onefile = TRUE,
      family = publication_font_family,
      bg = "white"
    )
  )
  if (cairo_attempt$opened) return("cairo_pdf")
  if (file.exists(filename)) unlink(filename)

  stop(
    paste0(
      "No vector PDF device with embedded-font support could be opened. ",
      "Quartz: ", device_attempt_message(quartz_attempt, "not applicable"),
      ". Cairo: ", device_attempt_message(cairo_attempt, "not available"),
      "."
    ),
    call. = FALSE
  )
}

render_atomic <- function(
    path, type, width_mm, height_mm, dpi, draw_function, figure_stem
) {
  ext <- if (type == "pdf") ".pdf" else ".png"
  tmp <- tempfile(
    pattern = paste0(".", tools::file_path_sans_ext(basename(path)), "."),
    tmpdir = dirname(path), fileext = ext
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  device_open <- FALSE
  pdf_method <- NA_character_
  tryCatch({
    if (type == "pdf") {
      pdf_method <- open_embedded_pdf_device(tmp, width_in, height_in)
    } else {
      grDevices::png(
        filename = tmp,
        width = round(width_in * dpi),
        height = round(height_in * dpi),
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
  if (type == "pdf") {
    pdf_device_records[[length(pdf_device_records) + 1L]] <<-
      tibble::tibble(
        figure = figure_stem,
        pdf_device = pdf_method,
        font_family = publication_font_family,
        vector_pdf = TRUE,
        embedded_font_device = TRUE
      )
    log_message(
      "PDF device for ", figure_stem, ": ", pdf_method,
      " (vector, embedded-font device; family=", publication_font_family, ")"
    )
  }
  invisible(path)
}

label_as_text <- function(x) {
  if (is.null(x)) return("")
  paste(vapply(
    as.list(x),
    function(z) paste(deparse(z, width.cutoff = 500L), collapse = ""),
    character(1)
  ), collapse = " | ")
}

collect_text_grobs <- function(
    grob, figure, panel, path = "root", inherited_fontsize = 12
) {
  local_fontsize <- inherited_fontsize
  if (!is.null(grob$gp) && !is.null(grob$gp$fontsize)) {
    candidate <- suppressWarnings(as.numeric(grob$gp$fontsize))
    if (length(candidate) > 0L && any(is.finite(candidate))) {
      local_fontsize <- min(candidate[is.finite(candidate)])
    }
  }

  rows <- list()
  if (inherits(grob, "text")) {
    label <- label_as_text(grob$label)
    visible <- nzchar(gsub("[[:space:]]", "", label))
    if (visible) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        figure = figure,
        panel = panel,
        grob_path = path,
        label = label,
        fontsize_pt = as.numeric(local_fontsize)
      )
    }
  }

  child_sets <- list()
  if (!is.null(grob$grobs) && length(grob$grobs) > 0L) {
    child_sets[["grobs"]] <- grob$grobs
  }
  if (!is.null(grob$children) && length(grob$children) > 0L) {
    child_sets[["children"]] <- as.list(grob$children)
  }

  for (set_name in names(child_sets)) {
    children <- child_sets[[set_name]]
    child_names <- names(children)
    if (is.null(child_names)) child_names <- rep("", length(children))
    for (i in seq_along(children)) {
      child_name <- child_names[[i]]
      if (!nzchar(child_name)) child_name <- as.character(i)
      rows[[length(rows) + 1L]] <- collect_text_grobs(
        children[[i]],
        figure = figure,
        panel = panel,
        path = paste(path, set_name, child_name, sep = "/"),
        inherited_fontsize = local_fontsize
      )
    }
  }
  dplyr::bind_rows(rows)
}

audit_plot_text <- function(plot, figure, panel) {
  forced_grob <- grid::grid.force(ggplot2::ggplotGrob(plot))
  audit <- collect_text_grobs(
    forced_grob,
    figure = figure,
    panel = panel
  )
  assert_true(
    nrow(audit) > 0L,
    paste0("No visible text grobs were found for ", figure, " / ", panel, ".")
  )
  assert_true(
    all(is.finite(audit$fontsize_pt)),
    paste0("Unresolved text size in ", figure, " / ", panel, ".")
  )
  audit
}

audit_figure_text <- function(plots, stem) {
  if (inherits(plots, "ggplot")) {
    plots <- list(panel_1 = plots)
  } else if (!is.list(plots)) {
    plots <- list(panel_1 = plots)
  }
  panel_names <- names(plots)
  if (is.null(panel_names) || any(!nzchar(panel_names))) {
    panel_names <- paste0("panel_", seq_along(plots))
  }
  audit <- dplyr::bind_rows(lapply(
    seq_along(plots),
    function(i) audit_plot_text(plots[[i]], stem, panel_names[[i]])
  ))
  minimum_observed <- min(audit$fontsize_pt)
  assert_true(
    minimum_observed + 1e-8 >= minimum_text_pt,
    paste0(
      stem, " contains text smaller than ", minimum_text_pt,
      " pt; observed minimum=", formatC(minimum_observed, digits = 3, format = "f"),
      " pt."
    )
  )
  text_audit_records[[length(text_audit_records) + 1L]] <<- audit
  log_message(
    "Final-size text audit PASS: ", stem,
    "; minimum=", formatC(minimum_observed, digits = 3, format = "f"),
    " pt; visible text grobs=", nrow(audit), "; threshold>=",
    formatC(minimum_text_pt, digits = 1, format = "f"), " pt"
  )
  invisible(audit)
}

save_plot_pair <- function(
    plot, stem, width_mm = publication_width_mm, height_mm, dpi = 300
) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  render_atomic(
    pdf_path, "pdf", width_mm, height_mm, dpi,
    function() draw_single_plot(plot), stem
  )
  render_atomic(
    png_path, "png", width_mm, height_mm, dpi,
    function() draw_single_plot(plot), stem
  )
  audit_figure_text(plot, stem)
  figure_output_records[[length(figure_output_records) + 1L]] <<-
    tibble::tibble(
      figure = stem,
      width_mm = width_mm,
      height_mm = height_mm,
      png_dpi = dpi,
      layout = "single ggplot"
    )
  log_message("Wrote figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

save_plot_grid_pair <- function(
    plots, stem, nrow, ncol, widths, heights,
    width_mm = publication_width_mm, height_mm, dpi = 300
) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  draw_fun <- function() {
    draw_plot_grid(
      plots = plots, nrow = nrow, ncol = ncol,
      widths = widths, heights = heights
    )
  }
  render_atomic(pdf_path, "pdf", width_mm, height_mm, dpi, draw_fun, stem)
  render_atomic(png_path, "png", width_mm, height_mm, dpi, draw_fun, stem)
  audit_figure_text(plots, stem)
  figure_output_records[[length(figure_output_records) + 1L]] <<-
    tibble::tibble(
      figure = stem,
      width_mm = width_mm,
      height_mm = height_mm,
      png_dpi = dpi,
      layout = paste0(nrow, "x", ncol, " composite")
    )
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
    "background_percentile", "matching_status", "empirical_p_two_sided",
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
    tissue = factor(tissue, levels = rev(tissue_order)),
    gene = factor(gene, levels = target_genes)
  )
tissue_axis_labels <- setNames(
  ifelse(
    tissue_order %in% mapped_gtex_tissues,
    paste0(tissue_order, " *"),
    tissue_order
  ),
  tissue_order
)

plot_atlas <- ggplot(atlas_long, aes(x = gene, y = tissue, fill = log2_tpm1)) +
  geom_tile(color = "white", linewidth = 0.20) +
  scale_fill_viridis_c(
    option = "C", begin = 0.05, end = 0.95,
    name = "log2(TPM + 1)"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(labels = tissue_axis_labels, drop = FALSE) +
  labs(
    x = NULL, y = NULL,
    title = "Normal-tissue expression architecture",
    subtitle = wrap_plot_text(paste(
      "Official GTEx v8 tissue-median TPM;",
      paste0(
        "* denotes a tissue category contributing to one of the ",
        n_mapped_organ_baselines,
        " TCGA-mapped organ baselines"
      )
    ), width = 96L),
    tag = "A"
  ) +
  theme_heatmap(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", size = 8.0),
    axis.text.y = element_text(
      angle = 0, hjust = 1, vjust = 0.5, size = 7.2
    ),
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
  theme_manuscript(base_size = 10) +
  theme(aspect.ratio = 0.42)

save_plot_grid_pair(
  plots = list(A = plot_atlas, B = plot_specificity),
  stem = "Figure1_round2F_GTEx_v8_normal_tissue_atlas",
  nrow = 2L, ncol = 1L,
  widths = 1, heights = c(2.25, 1),
  width_mm = publication_width_mm, height_mm = 265
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
    subtitle = wrap_plot_text(paste0(
      "* global target FDR<0.05 and |logFC|>=log2(1.5); ",
      "L=filtered for low expression"
    ), width = 105L)
  ) +
  theme_heatmap(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.0),
    axis.text.y = element_text(face = "italic", size = 8.0)
  )

save_plot_pair(
  plot_figure2,
  "Figure2_round2F_count_level_paired_heatmap",
  width_mm = publication_width_mm, height_mm = 92
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
  scale_fill_manual(values = gene_palette_okabe_ito, name = "Target gene") +
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
    subtitle = wrap_plot_text(paste0(
      "Grey distributions: 200 GTEx-defined identity genes per cancer; ",
      "39 target cells estimable, 12 below the identity median. ",
      "NE counts are shown beneath each cancer code."
    ), width = 105L)
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.0),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2)
  )

save_plot_pair(
  plot_identity,
  "Figure3A_round2F_tissue_identity_distribution",
  width_mm = publication_width_mm, height_mm = 132
)

# ---- 6. Figure 3B: matched background with two-sided FDR -------------------

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

# This reproduces the already-defined Round 2E display quantity. It does not
# fit or change the matched-background analysis.
matched_estimable <- matched_estimable %>%
  mutate(
    q_empirical_two_sided = p.adjust(empirical_p_two_sided, method = "BH")
  )

minimum_two_sided_fdr <- min(
  matched_estimable$q_empirical_two_sided,
  na.rm = TRUE
)
expected_minimum_q <- 0.1724709784411277
expected_minimum_p <- 0.005994005994005994
raw_minimum_cell <- matched_estimable %>%
  arrange(empirical_p_two_sided, gene, tcga_code) %>%
  slice(1L)
n_minimum_q_ties <- sum(
  abs(matched_estimable$q_empirical_two_sided - minimum_two_sided_fdr) < 1e-12
)

assert_true(
  abs(minimum_two_sided_fdr - expected_minimum_q) < 1e-10 &&
    !any(matched_estimable$q_empirical_two_sided < 0.05),
  "Two-sided matched-background FDR differs from the frozen Round 2E result."
)
assert_true(
  as.character(raw_minimum_cell$tcga_code) == "KICH" &&
    as.character(raw_minimum_cell$gene) == "MAOB" &&
    abs(raw_minimum_cell$empirical_p_two_sided - expected_minimum_p) < 1e-12 &&
    n_minimum_q_ties > 1L,
  "The minimum two-sided matched-background P/FDR contract is not preserved."
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
assert_true(
  all(fdr_band_counts$n > 0L),
  "All three two-sided FDR shape bands must be represented."
)

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
      fill = background_percentile, shape = two_sided_fdr_band
    ),
    inherit.aes = FALSE,
    color = "black", size = 2.7, stroke = 0.45
  ) +
  geom_label(
    data = matched_not_estimable,
    aes(x = 0, y = tcga_code, label = "NE"),
    inherit.aes = FALSE,
    fill = "grey90", color = "grey35", linewidth = 0.2,
    label.padding = grid::unit(0.08, "lines"), size = 2.65
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "grey80", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
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
    x = "Paired limma logFC (tumor - adjacent non-tumor)", y = NULL,
    title = "Cancer-level target differences within expression-matched backgrounds",
    subtitle = wrap_plot_text(paste0(
      "Segments begin at the matched-background median; minimum two-sided ",
      "FDR = 0.1725; no cell had FDR < 0.05. ",
      "NE is a status label without a numeric x-position; BLCA and HNSC are ",
      "absent because no valid GTEx mapping was specified."
    ), width = 100L)
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.subtitle = element_text(size = 8.0, lineheight = 1.08)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barwidth = grid::unit(4.0, "cm")),
    shape = guide_legend(order = 2, nrow = 1)
  )

save_plot_pair(
  plot_matched,
  "Figure3B_round2F_matched_background_context",
  width_mm = publication_width_mm, height_mm = 160
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

recurrence_annotation_rows <- recurrence_null %>%
  mutate(gene = factor(gene, levels = target_genes)) %>%
  transmute(
    gene,
    q_label = if_else(
      status == "estimable_same_gene_null",
      paste0("q=", formatC(q_recurrent_loss, digits = 3, format = "f")),
      ""
    ),
    profile_label = if_else(
      status == "estimable_same_gene_null",
      paste0("eligible profiles: ", as.integer(n_candidate_genes)),
      ""
    )
  )

observed_q_labels <- recurrence_annotation_rows %>%
  filter(nzchar(q_label)) %>%
  arrange(match(as.character(gene), target_genes)) %>%
  pull(q_label)
assert_true(
  identical(observed_q_labels, c("q=0.032", "q=0.032", "q=0.376")),
  "Figure 4A three-decimal q labels differ from the frozen values."
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
    aes(x = gene, y = observed_loss_fraction, fill = gene),
    inherit.aes = FALSE,
    shape = 23, size = 3.4, color = "black", stroke = 0.55
  ) +
  geom_text(
    data = recurrence_annotation_rows,
    aes(x = gene, y = 1.145, label = q_label),
    inherit.aes = FALSE,
    size = 2.85,
    fontface = "bold"
  ) +
  geom_text(
    data = recurrence_annotation_rows,
    aes(x = gene, y = 1.065, label = profile_label),
    inherit.aes = FALSE,
    size = 2.65,
    color = "grey25"
  ) +
  geom_label(
    data = not_evaluable_annotations,
    aes(x = gene, y = y, label = label),
    inherit.aes = FALSE,
    fill = "grey93", color = "grey30", linewidth = 0.25, size = 3
  ) +
  scale_fill_manual(
    values = gene_palette_okabe_ito,
    guide = "none",
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE, limits = target_genes) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = function(x) formatC(x, digits = 2, format = "f"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(
    x = NULL,
    y = "Fraction of cancers with loss",
    title = "Cross-cancer recurrence against the same-gene empirical null",
    subtitle = wrap_plot_text(paste0(
      "Grey distributions: eligible background-gene profiles at 80% ",
      "coverage; diamonds: observed target genes. Not evaluable is distinct ",
      "from evidence of no change."
    ), width = 105L)
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", size = 8.0),
    plot.subtitle = element_text(
      size = 8.0,
      color = "grey25",
      margin = margin(b = 32)
    ),
    plot.margin = margin(t = 6, r = 8, b = 6, l = 6)
  )

save_plot_pair(
  plot_recurrence,
  "Figure4A_round2F_same_gene_recurrence_null",
  width_mm = publication_width_mm, height_mm = 124
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
  scale_color_manual(
    values = gene_palette_okabe_ito,
    drop = FALSE,
    name = "Gene"
  ) +
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
    subtitle = wrap_plot_text(paste0(
      "MAOB remained below FDR 0.05 from 50%-80% coverage; CYP3A4 was ",
      "below 0.05 only at 80%. All 25 gene-coverage cells are represented."
    ), width = 105L)
  ) +
  theme_manuscript(base_size = 10) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(
    color = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2)
  )

save_plot_pair(
  plot_sensitivity,
  "Figure4B_round2F_coverage_sensitivity",
  width_mm = publication_width_mm, height_mm = 124
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
    subtitle = wrap_plot_text(paste0(
      "Non-LIHC cancers only; *=primary FDR<0.05 and |delta rho|>=0.20; ",
      "E=exploratory n=20-29 (desaturated); grey=not estimable"
    ), width = 105L)
  ) +
  theme_heatmap(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.0),
    axis.text.y = element_text(face = "italic", size = 8.0)
  )

save_plot_pair(
  plot_coexpression_heatmap,
  "SupplementaryFigureS1A_round2F_nonLIHC_coexpression_change_heatmap",
  width_mm = publication_width_mm, height_mm = 88
)

# ---- 10. Supplementary Figure S1B: representative rank plots ---------------

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
    all(abs(
      selection_check$expected_delta_rho -
        selection_check$observed_delta_rho
    ) < 1e-12),
  paste0(
    "Frozen representative cells do not follow the stated rule: largest ",
    "absolute delta rho in each cancer containing robust changes."
  )
)

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
    subtitle = wrap_plot_text(paste0(
      "Rank-scale display; the cell with the largest |delta rho| in each of ",
      "the three cancers containing robust changes"
    ), width = 100L)
  ) +
  theme_manuscript(base_size = 10) +
  theme(legend.position = "none")

save_plot_pair(
  plot_coexpression_representative,
  "SupplementaryFigureS1B_round2F_representative_coexpression_rank_scatter",
  width_mm = publication_width_mm, height_mm = 132
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
  "One or more frozen Round 2B / 2C inputs changed during Round 2F."
)

text_audit <- bind_rows(text_audit_records)
text_summary <- text_audit %>%
  group_by(figure) %>%
  summarise(
    n_visible_text_grobs = n(),
    minimum_fontsize_pt = min(fontsize_pt),
    maximum_fontsize_pt = max(fontsize_pt),
    threshold_pt = minimum_text_pt,
    all_text_at_least_7pt = all(fontsize_pt >= minimum_text_pt),
    .groups = "drop"
  )
figure_specs <- bind_rows(figure_output_records)
pdf_device_manifest <- bind_rows(pdf_device_records)

assert_true(
  n_distinct(text_summary$figure) == 8L &&
    all(text_summary$all_text_at_least_7pt),
  "The final-size text audit did not pass for all eight figures."
)
assert_true(
  nrow(figure_specs) == 8L &&
    all(figure_specs$width_mm == publication_width_mm) &&
    all(figure_specs$png_dpi == 300),
  "All eight figures must be written at 174 mm width and PNG 300 dpi."
)
assert_true(
  nrow(pdf_device_manifest) == 8L &&
    all(pdf_device_manifest$vector_pdf) &&
    all(pdf_device_manifest$embedded_font_device),
  "All eight PDFs must use a vector embedded-font device."
)

write_csv_atomic(text_audit, text_audit_file)
write_csv_atomic(text_summary, text_summary_file)
write_csv_atomic(figure_specs, figure_spec_file)
write_csv_atomic(pdf_device_manifest, pdf_device_file)

expected_figure_stems <- c(
  "Figure1_round2F_GTEx_v8_normal_tissue_atlas",
  "Figure2_round2F_count_level_paired_heatmap",
  "Figure3A_round2F_tissue_identity_distribution",
  "Figure3B_round2F_matched_background_context",
  "Figure4A_round2F_same_gene_recurrence_null",
  "Figure4B_round2F_coverage_sensitivity",
  "SupplementaryFigureS1A_round2F_nonLIHC_coexpression_change_heatmap",
  "SupplementaryFigureS1B_round2F_representative_coexpression_rank_scatter"
)
expected_figure_files <- unlist(lapply(
  expected_figure_stems,
  function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".png")))
))
assert_true(
  all(file.exists(expected_figure_files)) &&
    all(file.info(expected_figure_files)$size > 1000),
  "One or more Round 2F figure files are missing or unexpectedly small."
)

integrity_checks <- tibble::tribble(
  ~check, ~result, ~expected, ~passed, ~details,
  "frozen_inputs_unchanged",
  as.character(sum(input_md5_before == input_md5_after)),
  as.character(length(input_files)),
  identical(input_md5_before, input_md5_after),
  "MD5 before and after figure generation",
  "eight_figure_pairs_written",
  as.character(length(expected_figure_files)),
  "16",
  length(expected_figure_files) == 16L && all(file.exists(expected_figure_files)),
  "Eight PDFs and eight PNGs",
  "final_publication_width",
  paste(unique(figure_specs$width_mm), collapse = ","),
  "174 mm",
  all(figure_specs$width_mm == 174),
  "No post-hoc down-scaling assumed",
  "png_resolution",
  paste(unique(figure_specs$png_dpi), collapse = ","),
  "300 dpi",
  all(figure_specs$png_dpi == 300),
  "All PNG outputs",
  "minimum_text_size",
  formatC(min(text_audit$fontsize_pt), digits = 3, format = "f"),
  ">=7.000 pt",
  min(text_audit$fontsize_pt) >= 7,
  "All visible text grobs inspected after output at final dimensions",
  "vector_embedded_font_pdf_device",
  paste(unique(pdf_device_manifest$pdf_device), collapse = ","),
  "8 vector embedded-font PDF devices",
  nrow(pdf_device_manifest) == 8L &&
    all(pdf_device_manifest$embedded_font_device),
  paste0("Font family: ", publication_font_family),
  "figure1_transposed_atlas",
  paste0(n_distinct(atlas_long$tissue), " tissues x ",
         n_distinct(atlas_long$gene), " genes"),
  "54 tissues x 5 genes",
  n_distinct(atlas_long$tissue) == 54L && n_distinct(atlas_long$gene) == 5L,
  "Tissues on y axis; genes on x axis",
  "figure1_mapped_baselines",
  as.character(n_mapped_organ_baselines),
  "9",
  n_mapped_organ_baselines == 9L,
  "Twelve starred GTEx categories contribute to nine organ baselines",
  "figure3B_two_sided_FDR",
  formatC(minimum_two_sided_fdr, digits = 4, format = "f"),
  "0.1725 and zero q<0.05",
  abs(minimum_two_sided_fdr - expected_minimum_q) < 1e-10 &&
    !any(matched_estimable$q_empirical_two_sided < 0.05),
  "Round 2E two-sided display retained",
  "figure4A_q_labels",
  paste(observed_q_labels, collapse = "; "),
  "q=0.032; q=0.032; q=0.376",
  identical(observed_q_labels, c("q=0.032", "q=0.032", "q=0.376")),
  "Independent q row above the plotting panel",
  "paper2_firewall",
  as.character(any(coexpression_cells$tcga_code == "LIHC")),
  "FALSE",
  !any(coexpression_cells$tcga_code == "LIHC"),
  paste0(
    "No LIHC coexpression, survival, clinical, subgroup, immune, or ",
    "purity analysis"
  )
)
assert_true(
  all(integrity_checks$passed),
  "One or more final Round 2F integrity checks failed."
)
write_csv_atomic(integrity_checks, integrity_file)

completion_lines <- c(
  paste0("Analysis: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Project directory: ", project_dir),
  "Scope: final-size figure formatting from frozen Round 2B / Round 2C values only",
  paste0("Final width: ", publication_width_mm, " mm for all eight figures"),
  "PNG resolution: 300 dpi",
  paste0(
    "Final-size text audit: PASS for 8/8 figures; minimum=",
    formatC(min(text_audit$fontsize_pt), digits = 3, format = "f"),
    " pt (required >=7 pt)"
  ),
  paste0(
    "PDF output: vector embedded-font device; family=",
    publication_font_family
  ),
  "Official GTEx v8 cache: read only; 54 tissue categories validated",
  paste0(
    "Figure 1: atlas transposed to 54 tissue rows x 5 gene columns; ",
    "9 mapped organ baselines"
  ),
  "Figure 2: LIHC-MAOB and LIHC-CYP3A4 finite; no reserved legend",
  "Figure 3A: 2400 identity-gene rows; 39 target cells; 12 below median",
  paste0(
    "Figure 3B: minimum two-sided FDR=",
    formatC(minimum_two_sided_fdr, digits = 6, format = "f"),
    "; no FDR<0.05"
  ),
  paste0(
    "Figure 4A: q labels in external annotation row: ",
    paste(observed_q_labels, collapse = "; "),
    "; eligible profiles MAOB=66, CYP3A4=92, CYP3A7=313"
  ),
  "Figure 4B: all 25 coverage-gene statuses represented",
  paste0(
    "Supplementary Figure S1: 52 non-LIHC cells; LIHC coexpression absent; ",
    "representatives follow the largest-|delta rho| rule"
  ),
  paste0("Input files unchanged by MD5: ", length(input_files)),
  paste0("Figure files written: ", length(expected_figure_files)),
  paste0("Figures: ", figure_dir),
  paste0("Logs: ", log_dir)
)
write_lines_atomic(completion_lines, completion_report_file)

log_message("Round 2F complete.")
log_message("Figures: ", figure_dir)
log_message("Manifest: ", input_manifest_file)
log_message("Text audit: ", text_audit_file)
log_message("Integrity checks: ", integrity_file)
log_message("Completion report: ", completion_report_file)
log_message("Session information will be written last: ", session_info_file)
write_lines_atomic(capture.output(sessionInfo()), session_info_file)
