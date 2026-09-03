#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 2H: Supplementary Figures S2 and S3 only
# =============================================================================
#
# PURPOSE
#   Re-render only the following two supplementary figures from frozen Round 2B
#   outputs:
#     S2: normal baseline versus paired tumor-associated difference
#     S3: HPA-only normal-tissue RNA expression
#
# ANALYTICAL BOUNDARY
#   - No differential-expression, baseline-association, permutation,
#     recurrence, coexpression, survival, clinical, immune, purity, or pathway
#     analysis is rerun.
#   - The S2 effect estimate and P value are read from frozen Table 173.
#   - S2 point coordinates are reconstructed from frozen Table 171 exactly as
#     in Round 2B. Within-gene centering and the plotted lm smooth are display
#     operations only; they do not replace or update the frozen inference.
#   - S3 values are read from frozen Table 190; log2(value + 1) is the same
#     display transformation used in Round 2B.
#   - No missing values are imputed.
#
# OUTPUTS
#   results/round2H/figures/
#     SupplementaryFigureS2_round2H_normal_baseline_paired_difference.{pdf,png}
#     SupplementaryFigureS3_round2H_HPA_only_RNA_expression.{pdf,png}
#
#   results/round2H/logs/
#     round2H_progress.log
#     round2H_input_manifest.csv
#     round2H_text_size_audit.csv
#     round2H_figure_dimensions.csv
#     round2H_pdf_device_manifest.csv
#     round2H_colourbar_endpoint_audit.csv
#     round2H_integrity_checks.csv
#     round2H_completion_report.txt
#     round2H_sessionInfo.txt
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Dependencies and paths ----------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round2H_figures"
required_r_version <- "4.6.1"
publication_width_mm <- 174
minimum_text_pt <- 7
publication_font_family <- Sys.getenv(
  "MAOB_CYP3A_FIGURE_FONT",
  unset = "Arial"
)

target_genes <- c("MAOB", "CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
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

round2b_table_dir <- file.path(project_dir, "results", "round2B", "tables")
round2h_dir <- file.path(project_dir, "results", "round2H")
figure_dir <- file.path(round2h_dir, "figures")
log_dir <- file.path(round2h_dir, "logs")

invisible(lapply(
  c(round2h_dir, figure_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

input_files <- c(
  baseline_organ = file.path(
    round2b_table_dir,
    "171_round2B_baseline_disruption_organ_level_count.csv"
  ),
  baseline_model = file.path(
    round2b_table_dir,
    "173_round2B_baseline_disruption_pooled_model_count.csv"
  ),
  hpa_targets = file.path(
    round2b_table_dir,
    "190_round2B_HPA_only_RNA_target_tissues.csv"
  ),
  hpa_status = file.path(
    round2b_table_dir,
    "194_round2B_HPA_only_validation_status.csv"
  )
)

progress_log_file <- file.path(log_dir, "round2H_progress.log")
input_manifest_file <- file.path(log_dir, "round2H_input_manifest.csv")
text_audit_file <- file.path(log_dir, "round2H_text_size_audit.csv")
figure_spec_file <- file.path(log_dir, "round2H_figure_dimensions.csv")
pdf_device_file <- file.path(log_dir, "round2H_pdf_device_manifest.csv")
colourbar_audit_file <- file.path(
  log_dir,
  "round2H_colourbar_endpoint_audit.csv"
)
integrity_file <- file.path(log_dir, "round2H_integrity_checks.csv")
completion_report_file <- file.path(log_dir, "round2H_completion_report.txt")
session_info_file <- file.path(log_dir, "round2H_sessionInfo.txt")

required_packages <- c(
  "dplyr", "readr", "stringr", "tibble", "ggplot2", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 2H:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 2H requires R ", required_r_version,
    "; current version is ", as.character(getRversion()), ".",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

# ---- 1. General utilities ----------------------------------------------------

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
      rows = as.integer(rows),
      md5 = safe_md5(path),
      details = details
    )
  log_message("Input ", label, ": ", rows, " rows; ", details)
  invisible(NULL)
}

read_frozen_csv <- function(label, required_columns, expected_rows) {
  path <- input_files[[label]]
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  assert_has_columns(x, required_columns, label)
  assert_true(
    nrow(x) == expected_rows,
    paste0(
      label, " must contain ", expected_rows,
      " rows; found ", nrow(x), "."
    )
  )
  record_input(label, path, nrow(x), paste0(ncol(x), " columns"))
  x
}

wrap_plot_text <- function(x, width = 96L) {
  stringr::str_wrap(x, width = width)
}

theme_manuscript <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = publication_font_family) +
    theme(
      text = element_text(size = base_size, family = publication_font_family),
      axis.text = element_text(size = 8.0),
      axis.title = element_text(size = 9.0),
      legend.text = element_text(size = 7.5),
      legend.title = element_text(size = 8.0),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 8.0, color = "grey25"),
      plot.caption = element_text(size = 7.5),
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
      legend.text = element_text(size = 7.5),
      legend.title = element_text(size = 8.0),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 8.0, color = "grey25"),
      plot.caption = element_text(size = 7.5),
      panel.grid = element_blank(),
      legend.key.height = grid::unit(0.45, "cm")
    )
}

# Draw the continuous guide as vector rectangles. This prevents a Quartz PDF
# raster strip from being vertically inverted while tick labels stay fixed.
colourbar_formals <- names(formals(ggplot2::guide_colourbar))
colourbar_display_mode <- if ("display" %in% colourbar_formals) {
  "rectangles"
} else if ("raster" %in% colourbar_formals) {
  "raster=FALSE"
} else {
  stop(
    "The installed ggplot2 guide_colourbar API cannot request a non-raster guide.",
    call. = FALSE
  )
}

fixed_vertical_colourbar <- function() {
  args <- list(
    reverse = FALSE,
    direction = "vertical",
    nbin = 256
  )
  if ("display" %in% colourbar_formals) {
    args$display <- "rectangles"
  } else {
    args$raster <- FALSE
  }
  do.call(ggplot2::guide_colourbar, args)
}

# ---- 2. Atomic final-size figure output ------------------------------------

pdf_device_records <- list()
text_audit_records <- list()
figure_output_records <- list()
colourbar_audit_records <- list()

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
  if (!opened && as.integer(grDevices::dev.cur()) != device_before) {
    try(grDevices::dev.off(), silent = TRUE)
  }
  list(opened = isTRUE(opened), error = error_message, warnings = warning_messages)
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
    tmpdir = dirname(path),
    fileext = ext
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
    grob, figure, path = "root", inherited_fontsize = 12
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
    if (nzchar(gsub("[[:space:]]", "", label))) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        figure = figure,
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
        path = paste(path, set_name, child_name, sep = "/"),
        inherited_fontsize = local_fontsize
      )
    }
  }
  dplyr::bind_rows(rows)
}

audit_figure_text <- function(plot, stem) {
  forced_grob <- grid::grid.force(ggplot2::ggplotGrob(plot))
  audit <- collect_text_grobs(forced_grob, figure = stem)
  assert_true(nrow(audit) > 0L, paste0("No text grobs found for ", stem, "."))
  assert_true(
    all(is.finite(audit$fontsize_pt)),
    paste0("Unresolved text size in ", stem, ".")
  )
  minimum_observed <- min(audit$fontsize_pt)
  assert_true(
    minimum_observed + 1e-8 >= minimum_text_pt,
    paste0(
      stem, " contains text smaller than 7 pt; observed minimum=",
      formatC(minimum_observed, digits = 3, format = "f"), " pt."
    )
  )
  text_audit_records[[length(text_audit_records) + 1L]] <<- audit
  log_message(
    "Final-size text audit PASS: ", stem,
    "; minimum=", formatC(minimum_observed, digits = 3, format = "f"),
    " pt; threshold>=7.0 pt"
  )
  invisible(audit)
}

save_plot_pair <- function(
    plot, stem, height_mm, width_mm = publication_width_mm, dpi = 300
) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  render_atomic(
    pdf_path, "pdf", width_mm, height_mm, dpi,
    function() print(plot), stem
  )
  selected_pdf_device <-
    pdf_device_records[[length(pdf_device_records)]]$pdf_device[[1]]
  render_atomic(
    png_path, "png", width_mm, height_mm, dpi,
    function() print(plot), stem
  )

  audit_figure_text(plot, stem)
  figure_output_records[[length(figure_output_records) + 1L]] <<-
    tibble::tibble(
      figure = stem,
      width_mm = width_mm,
      height_mm = height_mm,
      png_dpi = dpi
    )
  log_message(
    "Wrote figure pair: ", stem,
    "; width=", width_mm, " mm; height=", height_mm,
    " mm; PNG=", dpi, " dpi; PDF device=", selected_pdf_device
  )
  invisible(c(pdf_path, png_path))
}

# ---- 3. Continuous-colourbar endpoint audit --------------------------------

format_scale_label <- function(x) {
  vapply(
    as.list(x),
    function(z) {
      if (is.expression(z) || is.language(z)) {
        paste(deparse(z, width.cutoff = 500L), collapse = "")
      } else {
        as.character(z)
      }
    },
    character(1)
  )
}

audit_continuous_colourbar <- function(
    plot, figure, low_meaning, high_meaning
) {
  built <- ggplot2::ggplot_build(plot)
  fill_scale <- built$plot$scales$get_scales("fill")
  assert_true(
    !is.null(fill_scale) && inherits(fill_scale, "ScaleContinuous"),
    paste0(figure, " does not contain the expected continuous fill scale.")
  )

  limits <- as.numeric(fill_scale$get_limits())
  assert_true(
    length(limits) == 2L && all(is.finite(limits)) && limits[[1]] < limits[[2]],
    paste0(figure, " has invalid continuous fill limits.")
  )
  breaks_all <- as.numeric(fill_scale$get_breaks())
  keep <- is.finite(breaks_all) &
    breaks_all >= limits[[1]] - 1e-12 &
    breaks_all <= limits[[2]] + 1e-12
  breaks <- breaks_all[keep]
  labels_all <- format_scale_label(fill_scale$get_labels(breaks_all))
  labels <- labels_all[keep]
  assert_true(
    length(breaks) >= 2L && length(labels) == length(breaks),
    paste0(figure, " requires at least two finite colourbar breaks/labels.")
  )

  extreme_indices <- c(which.min(breaks), which.max(breaks))
  result <- tibble::tibble(
    figure = figure,
    colourbar_end = c("bottom/minimum", "top/maximum"),
    scale_limit_value = c(limits[[1]], limits[[2]]),
    scale_limit_fill_colour = as.character(fill_scale$map(limits)),
    extreme_guide_break = breaks[extreme_indices],
    extreme_guide_label = labels[extreme_indices],
    fill_colour_at_guide_break = as.character(
      fill_scale$map(breaks[extreme_indices])
    ),
    biological_direction = c(low_meaning, high_meaning),
    guide_reverse = FALSE,
    guide_direction = "vertical",
    guide_display = colourbar_display_mode,
    endpoint_and_label_order_agree = TRUE
  )
  assert_true(
    result$extreme_guide_break[[1]] < result$extreme_guide_break[[2]] &&
      all(!is.na(result$scale_limit_fill_colour)) &&
      all(!is.na(result$fill_colour_at_guide_break)),
    paste0(figure, " colourbar endpoint audit failed.")
  )
  colourbar_audit_records[[length(colourbar_audit_records) + 1L]] <<- result

  for (i in seq_len(nrow(result))) {
    log_message(
      "Colourbar PASS [", figure, "; ", result$colourbar_end[[i]], "]: ",
      "scale endpoint=", formatC(
        result$scale_limit_value[[i]], digits = 4, format = "fg"
      ),
      ", endpoint fill=", result$scale_limit_fill_colour[[i]],
      "; corresponding extreme guide label=",
      result$extreme_guide_label[[i]],
      " at break=", formatC(
        result$extreme_guide_break[[i]], digits = 4, format = "fg"
      ),
      " with fill=", result$fill_colour_at_guide_break[[i]],
      "; meaning=", result$biological_direction[[i]],
      "; reverse=FALSE; display=", colourbar_display_mode
    )
  }
  invisible(result)
}

# ---- 4. Read frozen inputs ---------------------------------------------------

baseline_organ <- read_frozen_csv(
  "baseline_organ",
  c(
    "organ", "gene", "count_level_logFC", "baseline_z",
    "baseline_tpm", "baseline_log2_tpm1"
  ),
  expected_rows = 29L
)

baseline_model <- read_frozen_csv(
  "baseline_model",
  c(
    "n_cells", "n_organs", "n_genes", "beta_baseline_z",
    "robust_se", "p_value", "status"
  ),
  expected_rows = 1L
)

hpa_targets <- read_frozen_csv(
  "hpa_targets",
  c("gene", "tissue", "n_hpa_rows", "hpa_expression"),
  expected_rows = 200L
)

hpa_status <- read_frozen_csv(
  "hpa_status",
  c(
    "status", "source", "expression_column", "n_target_rows", "n_tissues"
  ),
  expected_rows = 1L
)

assert_true(
  !anyDuplicated(baseline_organ[c("organ", "gene")]) &&
    n_distinct(baseline_organ$organ) == 9L &&
    setequal(baseline_organ$gene, target_genes) &&
    all(is.finite(baseline_organ$baseline_z)) &&
    all(is.finite(baseline_organ$count_level_logFC)),
  "Frozen S2 organ-level table does not match the expected 29-cell design."
)
assert_true(
  baseline_model$n_cells[[1]] == 29L &&
    baseline_model$n_organs[[1]] == 9L &&
    baseline_model$n_genes[[1]] == 5L &&
    baseline_model$status[[1]] == "estimable_cluster_HC2" &&
    abs(baseline_model$beta_baseline_z[[1]] - (-0.6437335477423287)) < 1e-12 &&
    abs(baseline_model$p_value[[1]] - 0.0845671476214282) < 1e-12,
  "Frozen S2 model summary differs from the fact-locked values."
)
assert_true(
  !anyDuplicated(hpa_targets[c("gene", "tissue")]) &&
    n_distinct(hpa_targets$gene) == 5L &&
    n_distinct(hpa_targets$tissue) == 40L &&
    setequal(hpa_targets$gene, target_genes) &&
    all(is.finite(hpa_targets$hpa_expression)) &&
    all(hpa_targets$hpa_expression >= 0),
  "Frozen S3 table must contain a complete 5-gene x 40-tissue grid."
)
assert_true(
  hpa_status$status[[1]] == "completed" &&
    tolower(hpa_status$expression_column[[1]]) == "ntpm" &&
    hpa_status$n_target_rows[[1]] == 200L &&
    hpa_status$n_tissues[[1]] == 40L,
  "HPA status does not confirm 200 rows, 40 tissues, and the nTPM column."
)

# ---- 5. Supplementary Figure S2 --------------------------------------------

# For a model containing only an intercept and factor(gene), the fitted value
# for each observation is its gene mean. The following display coordinate is
# therefore exactly the residual used by the original Round 2B figure.
plot_baseline_data <- baseline_organ %>%
  mutate(gene = factor(gene, levels = target_genes)) %>%
  group_by(gene) %>%
  mutate(
    within_gene_logFC_residual =
      count_level_logFC - mean(count_level_logFC)
  ) %>%
  ungroup()

assert_true(
  all(is.finite(plot_baseline_data$within_gene_logFC_residual)) &&
    max(abs(
      plot_baseline_data %>%
        group_by(gene) %>%
        summarise(m = mean(within_gene_logFC_residual), .groups = "drop") %>%
        pull(m)
    )) < 1e-12,
  "S2 within-gene display residuals were not reconstructed correctly."
)

s2_title <- paste(
  "Exploratory relationship between normal baseline and paired",
  "tumor-associated difference"
)
s2_subtitle <- paste0(
  "Organ-level secondary analysis including LIHC paired effects; ",
  "cluster-HC2 beta=-0.644, p=0.0846"
)
assert_true(
  !str_detect(s2_title, regex("disruption", ignore_case = TRUE)) &&
    identical(
      s2_subtitle,
      paste0(
        "Organ-level secondary analysis including LIHC paired effects; ",
        "cluster-HC2 beta=-0.644, p=0.0846"
      )
    ),
  "S2 title/subtitle wording contract failed."
)

plot_s2 <- ggplot(
  plot_baseline_data,
  aes(
    x = baseline_z,
    y = within_gene_logFC_residual,
    color = gene
  )
) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_point(size = 2.3, alpha = 0.9) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey75",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = gene_palette_okabe_ito,
    breaks = target_genes,
    drop = FALSE,
    name = "Gene"
  ) +
  labs(
    x = "Within-gene standardized GTEx normal baseline",
    y = "Paired logFC residual after removing gene means",
    title = wrap_plot_text(s2_title, width = 82L),
    subtitle = wrap_plot_text(s2_subtitle, width = 100L)
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 8.0, lineheight = 1.08),
    plot.margin = margin(t = 7, r = 8, b = 7, l = 7)
  )

save_plot_pair(
  plot_s2,
  "SupplementaryFigureS2_round2H_normal_baseline_paired_difference",
  width_mm = publication_width_mm,
  height_mm = 125
)

# ---- 6. Supplementary Figure S3 --------------------------------------------

hpa_tissue_order <- hpa_targets %>%
  group_by(tissue) %>%
  summarise(max_expression = max(hpa_expression), .groups = "drop") %>%
  arrange(max_expression, tissue) %>%
  pull(tissue)

plot_hpa_data <- hpa_targets %>%
  mutate(
    gene = factor(gene, levels = target_genes),
    # With tissue on the y axis, this preserves the Round 2B low-to-high order:
    # low-expression tissues at the bottom and high-expression tissues at top.
    tissue = factor(tissue, levels = hpa_tissue_order),
    hpa_log2_expression = log2(hpa_expression + 1)
  )

s3_title <- "HPA-only normal-tissue RNA expression"
s3_subtitle <- paste(
  "Human Protein Atlas RNA data; nTPM column; descriptive",
  "comparison of tissue ranking"
)
assert_true(
  !str_detect(
    paste(s3_title, s3_subtitle),
    regex("validation|independent", ignore_case = TRUE)
  ),
  "S3 figure text must not contain 'validation' or 'independent'."
)

plot_s3 <- ggplot(
  plot_hpa_data,
  aes(x = gene, y = tissue, fill = hpa_log2_expression)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    option = "C",
    name = "log2(HPA RNA+1)",
    guide = fixed_vertical_colourbar()
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    title = s3_title,
    subtitle = wrap_plot_text(s3_subtitle, width = 92L)
  ) +
  theme_heatmap(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", size = 8.0),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 7.2
    ),
    legend.position = "right",
    plot.subtitle = element_text(size = 8.0, lineheight = 1.08),
    plot.margin = margin(t = 7, r = 8, b = 7, l = 7)
  )

audit_continuous_colourbar(
  plot_s3,
  "SupplementaryFigureS3_round2H_HPA_only_RNA_expression",
  low_meaning = "low HPA RNA expression / dark purple",
  high_meaning = "high HPA RNA expression / yellow"
)
save_plot_pair(
  plot_s3,
  "SupplementaryFigureS3_round2H_HPA_only_RNA_expression",
  width_mm = publication_width_mm,
  height_mm = 190
)

# ---- 7. Final checks and logs -----------------------------------------------

input_md5_after <- vapply(input_files, safe_md5, character(1))
assert_true(
  identical(input_md5_before, input_md5_after),
  "One or more frozen Round 2B inputs changed during Round 2H."
)

input_manifest <- bind_rows(input_manifest_records)
text_audit <- bind_rows(text_audit_records)
figure_specs <- bind_rows(figure_output_records)
pdf_device_manifest <- bind_rows(pdf_device_records)
colourbar_audit <- bind_rows(colourbar_audit_records)

assert_true(
  nrow(input_manifest) == 4L,
  "The Round 2H input manifest must contain four frozen inputs."
)
assert_true(
  n_distinct(text_audit$figure) == 2L &&
    all(text_audit$fontsize_pt >= minimum_text_pt),
  "The 7-pt text audit did not pass for both Round 2H figures."
)
assert_true(
  nrow(figure_specs) == 2L &&
    all(figure_specs$width_mm == publication_width_mm) &&
    all(figure_specs$png_dpi == 300),
  "Both Round 2H figures must be 174 mm wide with 300-dpi PNG output."
)
assert_true(
  nrow(pdf_device_manifest) == 2L &&
    all(pdf_device_manifest$vector_pdf) &&
    all(pdf_device_manifest$embedded_font_device),
  "Both Round 2H PDFs must use a vector embedded-font device."
)
assert_true(
  nrow(colourbar_audit) == 2L &&
    n_distinct(colourbar_audit$figure) == 1L &&
    all(colourbar_audit$endpoint_and_label_order_agree) &&
    all(!colourbar_audit$guide_reverse) &&
    all(colourbar_audit$guide_direction == "vertical"),
  "The S3 continuous colourbar did not pass its endpoint audit."
)

write_csv_atomic(input_manifest, input_manifest_file)
write_csv_atomic(text_audit, text_audit_file)
write_csv_atomic(figure_specs, figure_spec_file)
write_csv_atomic(pdf_device_manifest, pdf_device_file)
write_csv_atomic(colourbar_audit, colourbar_audit_file)

expected_stems <- c(
  "SupplementaryFigureS2_round2H_normal_baseline_paired_difference",
  "SupplementaryFigureS3_round2H_HPA_only_RNA_expression"
)
expected_files <- unlist(lapply(
  expected_stems,
  function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".png")))
))
assert_true(
  length(expected_files) == 4L &&
    all(file.exists(expected_files)) &&
    all(file.info(expected_files)$size > 1000),
  "One or more Round 2H figure files are missing or unexpectedly small."
)

integrity_checks <- tibble::tribble(
  ~check, ~result, ~expected, ~passed, ~details,
  "frozen_inputs_unchanged",
  as.character(sum(input_md5_before == input_md5_after)),
  "4",
  identical(input_md5_before, input_md5_after),
  "Input MD5 before and after rendering",
  "only_two_figure_pairs",
  as.character(length(expected_files)),
  "4 files",
  length(expected_files) == 4L && all(file.exists(expected_files)),
  "Only S2 and S3 are generated",
  "final_publication_width",
  paste(unique(figure_specs$width_mm), collapse = ","),
  "174 mm",
  all(figure_specs$width_mm == 174),
  "Direct device output at 174 mm; no post-hoc down-scaling assumed",
  "minimum_text_size",
  formatC(min(text_audit$fontsize_pt), digits = 3, format = "f"),
  ">=7.000 pt",
  min(text_audit$fontsize_pt) >= 7,
  "Visible text grobs inspected at final dimensions",
  "S2_fact_locked_model_label",
  s2_subtitle,
  "cluster-HC2 beta=-0.644, p=0.0846",
  str_detect(s2_subtitle, fixed("cluster-HC2 beta=-0.644, p=0.0846")),
  "Effect estimate and P value read from Table 173",
  "S2_title_terminology",
  s2_title,
  "No internal 'disruption' wording",
  !str_detect(s2_title, regex("disruption", ignore_case = TRUE)),
  "Paired tumor-associated difference terminology used",
  "S2_okabe_ito_palette",
  paste(names(gene_palette_okabe_ito), gene_palette_okabe_ito, collapse = "; "),
  "Five target genes",
  identical(names(gene_palette_okabe_ito), target_genes),
  "Categorical palette changed; data values unchanged",
  "S3_transposed_grid",
  paste0(n_distinct(plot_hpa_data$tissue), " tissues x ",
         n_distinct(plot_hpa_data$gene), " genes"),
  "40 tissues on y x 5 genes on x",
  n_distinct(plot_hpa_data$tissue) == 40L &&
    n_distinct(plot_hpa_data$gene) == 5L,
  "Horizontal tissue labels; complete 200-cell grid",
  "S3_wording_boundary",
  paste(s3_title, s3_subtitle),
  "No validation/independent wording",
  !str_detect(
    paste(s3_title, s3_subtitle),
    regex("validation|independent", ignore_case = TRUE)
  ),
  "Descriptive comparison wording used",
  "S3_colourbar_contract",
  paste0("reverse=FALSE; direction=vertical; display=", colourbar_display_mode),
  "minimum/purple at bottom; maximum/yellow at top",
  nrow(colourbar_audit) == 2L &&
    all(colourbar_audit$endpoint_and_label_order_agree) &&
    identical(
      colourbar_audit$biological_direction,
      c(
        "low HPA RNA expression / dark purple",
        "high HPA RNA expression / yellow"
      )
    ),
  paste0(
    "Actual endpoint colours and corresponding guide labels recorded in ",
    "round2H_colourbar_endpoint_audit.csv"
  )
)
assert_true(
  all(integrity_checks$passed),
  "One or more Round 2H integrity checks failed."
)
write_csv_atomic(integrity_checks, integrity_file)

completion_lines <- c(
  paste0("Analysis: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Project directory: ", project_dir),
  "Scope: formatting-only re-render of Supplementary Figures S2 and S3",
  "S2: terminology aligned to paired tumor-associated difference",
  "S2: fact-locked subtitle retained: cluster-HC2 beta=-0.644, p=0.0846",
  "S2: five-gene categorical colours changed to Okabe-Ito",
  "S3: tissues transposed to the y axis; genes placed on the x axis",
  paste0("S3 colourbar: reverse=FALSE; display=", colourbar_display_mode),
  paste0("S3 colourbar endpoint audit: PASS; ", colourbar_audit_file),
  paste0(
    "Final-size text audit: PASS; minimum=",
    formatC(min(text_audit$fontsize_pt), digits = 3, format = "f"),
    " pt"
  ),
  "Both figures: width=174 mm; PNG=300 dpi; vector PDF",
  paste0("Input files unchanged by MD5: ", length(input_files)),
  paste0("Figures: ", figure_dir),
  paste0("Logs: ", log_dir)
)
write_lines_atomic(completion_lines, completion_report_file)

log_message("Round 2H complete.")
log_message("Figures: ", figure_dir)
log_message("Integrity checks: ", integrity_file)
log_message("Session information will be written last: ", session_info_file)
write_lines_atomic(capture.output(sessionInfo()), session_info_file)
