#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 2G: targeted figure corrections only
# =============================================================================
#
# PURPOSE
#   Re-render five figures from frozen Round 2B values/caches:
#     1. Figure 1: correct the continuous-colourbar orientation only.
#     2. Figure 2: correct the continuous-colourbar orientation only.
#     3. Figure 4A: separate subtitle, q/profile annotation rows, and caption.
#     4. Supplementary Figure S1A: correct the continuous-colourbar
#        orientation only.
#     5. Supplementary Figure S1B: short facet strips and in-panel statistics.
#
# ANALYTICAL BOUNDARY
#   - No model, correlation, bootstrap, differential-expression result,
#     recurrence null, or other biological statistic is recomputed.
#   - Values are read from frozen Tables 121, 140, 142, 162, 163, 180, 181,
#     and 182, plus the read-only official-GTEx-v8 cache used in Round 2F.
#   - Formatting to two or three displayed decimal places is display-only.
#   - No LIHC coexpression, survival, clinical, subgroup, immune, purity, or
#     pathway analysis is performed.
#   - Figures 3A, 3B, and 4B are not regenerated or overwritten.
#   - Continuous scales and tile colours are unchanged. The three affected
#     legends use a non-raster vertical guide with reverse=FALSE to prevent
#     device-specific reversal of the colour strip relative to its labels.
#
# OUTPUTS
#   results/round2G/figures/
#     Figure1_round2G_GTEx_v8_normal_tissue_atlas.{pdf,png}
#     Figure2_round2G_count_level_paired_heatmap.{pdf,png}
#     Figure4A_round2G_same_gene_recurrence_null.{pdf,png}
#     SupplementaryFigureS1A_round2G_nonLIHC_coexpression_change_heatmap.{pdf,png}
#     SupplementaryFigureS1B_round2G_representative_coexpression_rank_scatter.{pdf,png}
#
#   results/round2G/logs/
#     round2G_progress.log
#     round2G_input_manifest.csv
#     round2G_text_size_audit.csv
#     round2G_figure_dimensions.csv
#     round2G_pdf_device_manifest.csv
#     round2G_colourbar_endpoint_audit.csv
#     round2G_integrity_checks.csv
#     round2G_completion_report.txt
#     round2G_sessionInfo.txt
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration --------------------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round2G_figures"
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
coexpression_tcga_order <- setdiff(tcga_order, "LIHC")

gene_palette <- c(
  MAOB = "#222222",
  CYP3A4 = "#D55E00",
  CYP3A5 = "#0072B2",
  CYP3A7 = "#009E73",
  CYP3A43 = "#CC79A7"
)

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
round2b_table_dir <- file.path(project_dir, "results", "round2B", "tables")
round2g_dir <- file.path(project_dir, "results", "round2G")
figure_dir <- file.path(round2g_dir, "figures")
log_dir <- file.path(round2g_dir, "logs")

invisible(lapply(
  c(round2g_dir, figure_dir, log_dir),
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
  recurrence_null = file.path(
    round2b_table_dir,
    "162_round2B_same_gene_recurrence_null.csv"
  ),
  recurrence_profiles = file.path(
    round2b_table_dir,
    "163_round2B_same_gene_null_gene_profiles.csv"
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

progress_log_file <- file.path(log_dir, "round2G_progress.log")
input_manifest_file <- file.path(log_dir, "round2G_input_manifest.csv")
text_audit_file <- file.path(log_dir, "round2G_text_size_audit.csv")
figure_spec_file <- file.path(log_dir, "round2G_figure_dimensions.csv")
pdf_device_file <- file.path(log_dir, "round2G_pdf_device_manifest.csv")
colourbar_audit_file <- file.path(
  log_dir,
  "round2G_colourbar_endpoint_audit.csv"
)
integrity_file <- file.path(log_dir, "round2G_integrity_checks.csv")
completion_report_file <- file.path(log_dir, "round2G_completion_report.txt")
session_info_file <- file.path(log_dir, "round2G_sessionInfo.txt")

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "tibble", "ggplot2", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 2G:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 2G requires R ", required_r_version,
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
      rows = as.integer(rows),
      md5 = safe_md5(path),
      details = details
    )
  log_message("Input ", label, ": ", rows, " rows; ", details)
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
      plot.caption = element_text(size = 7.5, color = "grey25"),
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

wrap_plot_text <- function(x, width = 92L) {
  stringr::str_wrap(x, width = width)
}

# A raster colour strip can be vertically inverted by some Quartz-based PDF
# renderers while its tick labels remain in numeric order. Drawing the strip as
# vector rectangles avoids that device-specific mismatch. reverse=FALSE fixes
# the intended vertical convention: minimum/bottom and maximum/top.
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

fixed_vertical_colourbar <- function(order = 0L) {
  args <- list(
    reverse = FALSE,
    direction = "vertical",
    nbin = 256,
    order = order
  )
  if ("display" %in% colourbar_formals) {
    args$display <- "rectangles"
  } else {
    args$raster <- FALSE
  }
  do.call(ggplot2::guide_colourbar, args)
}

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

pdf_device_records <- list()
text_audit_records <- list()
figure_output_records <- list()
colourbar_audit_records <- list()

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
  assert_true(
    length(breaks) >= 2L,
    paste0(figure, " requires at least two finite colourbar breaks.")
  )
  labels_all <- format_scale_label(fill_scale$get_labels(breaks_all))
  labels <- labels_all[keep]

  endpoint_values <- c(limits[[1]], limits[[2]])
  endpoint_colours <- as.character(fill_scale$map(endpoint_values))
  extreme_indices <- c(which.min(breaks), which.max(breaks))
  extreme_breaks <- breaks[extreme_indices]
  extreme_labels <- labels[extreme_indices]
  extreme_colours <- as.character(fill_scale$map(extreme_breaks))

  result <- tibble::tibble(
    figure = figure,
    colourbar_end = c("bottom/minimum", "top/maximum"),
    scale_limit_value = endpoint_values,
    scale_limit_fill_colour = endpoint_colours,
    extreme_guide_break = extreme_breaks,
    extreme_guide_label = extreme_labels,
    fill_colour_at_guide_break = extreme_colours,
    biological_direction = c(low_meaning, high_meaning),
    guide_reverse = FALSE,
    guide_direction = "vertical",
    guide_display = colourbar_display_mode,
    endpoint_and_label_order_agree = c(
      extreme_breaks[[1]] < extreme_breaks[[2]],
      extreme_breaks[[2]] > extreme_breaks[[1]]
    )
  )
  assert_true(
    all(result$endpoint_and_label_order_agree) &&
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
  log_message("Wrote figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

draw_plot_grid <- function(plots, nrow, ncol, widths, heights) {
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

save_plot_grid_pair <- function(
    plots, stem, nrow, ncol, widths, heights,
    height_mm, width_mm = publication_width_mm, dpi = 300
) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  draw_fun <- function() {
    draw_plot_grid(plots, nrow, ncol, widths, heights)
  }
  render_atomic(
    pdf_path, "pdf", width_mm, height_mm, dpi, draw_fun, stem
  )
  render_atomic(
    png_path, "png", width_mm, height_mm, dpi, draw_fun, stem
  )
  invisible(lapply(
    seq_along(plots),
    function(i) audit_figure_text(
      plots[[i]], stem
    )
  ))
  figure_output_records[[length(figure_output_records) + 1L]] <<-
    tibble::tibble(
      figure = stem,
      width_mm = width_mm,
      height_mm = height_mm,
      png_dpi = dpi
    )
  log_message("Wrote composite figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

# ---- 2. Read frozen inputs ---------------------------------------------------

target_results <- read_frozen_csv(
  "target_results",
  c(
    "tcga_code", "gene", "n_pairs", "logFC", "analysis_status",
    "robust_count_level_change"
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

recurrence_null <- read_frozen_csv(
  "recurrence_null",
  c(
    "gene", "minimum_coverage_fraction", "n_cancers",
    "observed_loss_fraction", "n_candidate_genes", "status",
    "q_recurrent_loss"
  ),
  expected_rows = 5L
)

recurrence_profiles <- read_frozen_csv(
  "recurrence_profiles",
  c("target_gene", "loss_fraction", "eligible_primary_null")
)

coexpression_cells <- read_frozen_csv(
  "coexpression_cells",
  c("tcga_code", "partner", "delta_rho", "status", "robust_rewiring"),
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
    "MAOB_rank_fraction", "partner_rank_fraction",
    "rho_normal", "rho_tumor", "delta_rho"
  ),
  expected_rows = 434L
)

assert_true(
  !anyDuplicated(target_results[c("tcga_code", "gene")]) &&
    setequal(target_results$tcga_code, tcga_order) &&
    setequal(target_results$gene, target_genes),
  "Target table does not contain the expected 14-cancer x 5-gene grid."
)

lihc_figure2_cells <- target_results %>%
  filter(tcga_code == "LIHC", gene %in% c("MAOB", "CYP3A4"))
assert_true(
  nrow(lihc_figure2_cells) == 2L &&
    all(is.finite(lihc_figure2_cells$logFC)) &&
    !any(lihc_figure2_cells$analysis_status == "reserved_not_recomputed"),
  "Figure 2 contract failed: LIHC-MAOB/CYP3A4 are not both finite."
)

assert_true(
  !any(coexpression_cells$tcga_code == "LIHC") &&
    !any(coexpression_points$tcga_code == "LIHC"),
  "Paper-boundary failure: LIHC coexpression data are present."
)

# ---- 3. Figure 1: GTEx atlas with corrected colourbar -----------------------

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
    length(gtex_tissues) == 54L &&
    identical(colnames(gtex_tpm_matrix), gtex_tissues),
  "Figure 1A requires the frozen 54-category official GTEx v8 matrix."
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
  "The GTEx cache does not contain one row for every target gene."
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
n_mapped_organ_baselines <- valid_mapping %>% distinct(organ) %>% nrow()
assert_true(
  length(mapped_gtex_tissues) == 12L &&
    n_mapped_organ_baselines == 9L &&
    all(mapped_gtex_tissues %in% gtex_tissues),
  "Expected 12 GTEx categories contributing to nine organ baselines."
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
    name = "log2(TPM + 1)",
    guide = fixed_vertical_colourbar()
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

audit_continuous_colourbar(
  plot_atlas,
  "Figure1_round2G_GTEx_v8_normal_tissue_atlas",
  low_meaning = "low expression / dark purple",
  high_meaning = "high expression / yellow"
)
save_plot_grid_pair(
  plots = list(plot_atlas, plot_specificity),
  stem = "Figure1_round2G_GTEx_v8_normal_tissue_atlas",
  nrow = 2L, ncol = 1L,
  widths = 1, heights = c(2.25, 1),
  width_mm = publication_width_mm, height_mm = 265
)

# ---- 4. Figure 2: paired heatmap with corrected colourbar ------------------

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
    pair_labels %>%
      select(tcga_code, x_label) %>%
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
    na.value = "grey85", name = "limma logFC\n(tumor - normal)",
    guide = fixed_vertical_colourbar()
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

audit_continuous_colourbar(
  plot_figure2,
  "Figure2_round2G_count_level_paired_heatmap",
  low_meaning = "negative tumor-minus-normal logFC / blue",
  high_meaning = "positive tumor-minus-normal logFC / red"
)
save_plot_pair(
  plot_figure2,
  "Figure2_round2G_count_level_paired_heatmap",
  width_mm = publication_width_mm,
  height_mm = 92
)

# ---- 5. Figure 4A: recurrence annotation spacing ----------------------------

recurrence_evaluable_genes <- recurrence_null %>%
  filter(status == "estimable_same_gene_null") %>%
  pull(gene)
assert_true(
  setequal(recurrence_evaluable_genes, c("MAOB", "CYP3A4", "CYP3A7")),
  "Unexpected evaluable genes in the primary recurrence null."
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
      profile_count_check$expected[
        match(c("MAOB", "CYP3A4", "CYP3A7"), profile_count_check$gene)
      ],
      c(66L, 92L, 313L)
    ),
  "Eligible recurrence-profile counts do not equal 66/92/313."
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
  "Figure 4A q labels differ from the frozen values."
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

figure4a_subtitle <- paste0(
  "Grey distributions: eligible background-gene profiles at 80% coverage; ",
  "diamonds: observed target genes."
)
figure4a_caption <- "Not evaluable is distinct from evidence of no change."
assert_true(
  !grepl("\n", figure4a_subtitle, fixed = TRUE),
  "Figure 4A subtitle must remain one line."
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
    aes(x = gene, y = 1.105, label = q_label),
    inherit.aes = FALSE,
    size = 2.85,
    fontface = "bold"
  ) +
  geom_text(
    data = recurrence_annotation_rows,
    aes(x = gene, y = 1.035, label = profile_label),
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
    subtitle = figure4a_subtitle,
    caption = figure4a_caption
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", size = 8.0),
    plot.subtitle = element_text(
      size = 8.0,
      color = "grey25",
      margin = margin(b = 36)
    ),
    plot.caption.position = "plot",
    plot.caption = element_text(
      size = 7.5,
      color = "grey25",
      hjust = 0,
      margin = margin(t = 7)
    ),
    plot.margin = margin(t = 10, r = 8, b = 6, l = 6)
  )

save_plot_pair(
  plot_recurrence,
  "Figure4A_round2G_same_gene_recurrence_null",
  width_mm = publication_width_mm,
  height_mm = 132
)

# ---- 6. Supplementary Figure S1A: corrected colourbar ----------------------

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
    name = "delta rho\n(tumor - normal)",
    guide = fixed_vertical_colourbar()
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

audit_continuous_colourbar(
  plot_coexpression_heatmap,
  "SupplementaryFigureS1A_round2G_nonLIHC_coexpression_change_heatmap",
  low_meaning = "negative delta rho / purple",
  high_meaning = "positive delta rho / green"
)
save_plot_pair(
  plot_coexpression_heatmap,
  "SupplementaryFigureS1A_round2G_nonLIHC_coexpression_change_heatmap",
  width_mm = publication_width_mm,
  height_mm = 88
)

# ---- 7. Supplementary Figure S1B: statistics outside facet strips ----------

robust_cells <- coexpression_cells %>%
  filter(robust_rewiring, is.finite(delta_rho))
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
    setequal(unique(robust_cells$tcga_code), c("BRCA", "KIRC", "KIRP")),
  "Expected seven robust cells in BRCA, KIRC, and KIRP."
)
assert_true(
  nrow(selection_check) == 3L &&
    all(selection_check$expected_partner == selection_check$observed_partner) &&
    all(abs(
      selection_check$expected_delta_rho -
        selection_check$observed_delta_rho
    ) < 1e-12),
  "Representative cells do not follow the frozen largest-|delta rho| rule."
)

point_metric_check <- coexpression_points %>%
  distinct(tcga_code, partner, rho_normal, rho_tumor, delta_rho) %>%
  rename(
    point_rho_normal = rho_normal,
    point_rho_tumor = rho_tumor,
    point_delta_rho = delta_rho
  ) %>%
  left_join(
    coexpression_representatives %>%
      rename(
        representative_rho_normal = rho_normal,
        representative_rho_tumor = rho_tumor,
        representative_delta_rho = delta_rho
      ),
    by = c("tcga_code", "partner")
  )
assert_true(
  nrow(point_metric_check) == 3L &&
    all(abs(
      point_metric_check$point_rho_normal -
        point_metric_check$representative_rho_normal
    ) < 1e-12) &&
    all(abs(
      point_metric_check$point_rho_tumor -
        point_metric_check$representative_rho_tumor
    ) < 1e-12) &&
    all(abs(
      point_metric_check$point_delta_rho -
        point_metric_check$representative_delta_rho
    ) < 1e-12),
  "Tables 181 and 182 do not contain identical rho/delta values."
)

panel_order <- coexpression_representatives %>%
  arrange(match(tcga_code, c("BRCA", "KIRC", "KIRP"))) %>%
  transmute(panel_title = paste0(tcga_code, ": MAOB-", partner)) %>%
  pull(panel_title)

representative_display <- coexpression_representatives %>%
  mutate(
    panel_title = paste0(tcga_code, ": MAOB-", partner),
    panel_title = factor(panel_title, levels = panel_order),
    normal_label = paste0(
      "normal rho=", formatC(rho_normal, digits = 2, format = "f")
    ),
    tumor_label = paste0(
      "tumor rho=", formatC(rho_tumor, digits = 2, format = "f"),
      "\ndelta=", formatC(delta_rho, digits = 2, format = "f")
    )
  )

annotation_data <- bind_rows(
  representative_display %>%
    transmute(
      tcga_code, partner, panel_title,
      condition = "Normal", statistic_label = normal_label,
      label_y = -0.17
    ),
  representative_display %>%
    transmute(
      tcga_code, partner, panel_title,
      condition = "Tumor", statistic_label = tumor_label,
      label_y = -0.235
    )
) %>%
  mutate(condition = factor(condition, levels = c("Normal", "Tumor")))

expected_stat_labels <- c(
  "normal rho=0.36", "normal rho=0.88", "normal rho=0.94",
  "tumor rho=-0.09\ndelta=-0.45",
  "tumor rho=0.15\ndelta=-0.73",
  "tumor rho=0.26\ndelta=-0.68"
)
assert_true(
  identical(annotation_data$statistic_label, expected_stat_labels),
  "Displayed S1B rho/delta labels differ from the frozen values."
)
assert_true(
  length(panel_order) == 3L &&
    all(!grepl("\n", panel_order, fixed = TRUE)),
  "S1B facet strips must contain three single-line cancer/gene labels only."
)

coexpression_points_plot <- coexpression_points %>%
  left_join(
    representative_display %>%
      select(tcga_code, partner, panel_title),
    by = c("tcga_code", "partner")
  ) %>%
  mutate(
    condition = factor(condition, levels = c("Normal", "Tumor")),
    panel_title = factor(panel_title, levels = panel_order)
  )

annotation_bands <- annotation_data %>%
  distinct(condition, panel_title)

plot_coexpression_representative <- ggplot(
  coexpression_points_plot,
  aes(
    x = MAOB_rank_fraction,
    y = partner_rank_fraction,
    color = condition
  )
) +
  geom_rect(
    data = annotation_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.005),
    inherit.aes = FALSE,
    fill = "grey97",
    color = NA
  ) +
  geom_hline(yintercept = 0, color = "grey78", linewidth = 0.35) +
  geom_point(alpha = 0.65, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  geom_text(
    data = annotation_data,
    aes(
      x = 0.98,
      y = label_y,
      label = statistic_label,
      color = condition
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0,
    size = 2.65,
    lineheight = 0.95,
    fontface = "plain"
  ) +
  facet_grid(condition ~ panel_title) +
  scale_color_manual(values = c(Normal = "#2166AC", Tumor = "#B2182B")) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = 0)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(-0.26, 1), ratio = 1, clip = "on") +
  labs(
    x = "Within-condition MAOB expression rank",
    y = "Within-condition CYP3A expression rank",
    title = "Representative non-LIHC coexpression changes",
    subtitle = paste0(
      "Rank-scale display; the cell with the largest |delta rho| in each of ",
      "the three cancers containing robust changes"
    )
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 8.0, face = "bold"),
    strip.text.y = element_text(size = 8.0),
    panel.spacing = grid::unit(2.5, "mm")
  )

save_plot_pair(
  plot_coexpression_representative,
  "SupplementaryFigureS1B_round2G_representative_coexpression_rank_scatter",
  width_mm = publication_width_mm,
  height_mm = 142
)

# ---- 8. Final checks and session information --------------------------------

input_md5_after <- vapply(input_files, safe_md5, character(1))
assert_true(
  identical(input_md5_before, input_md5_after),
  "One or more frozen Round 2B inputs changed during Round 2G."
)

input_manifest <- bind_rows(input_manifest_records)
text_audit <- bind_rows(text_audit_records)
figure_specs <- bind_rows(figure_output_records)
pdf_device_manifest <- bind_rows(pdf_device_records)
colourbar_audit <- bind_rows(colourbar_audit_records)

assert_true(
  n_distinct(text_audit$figure) == 5L &&
    all(text_audit$fontsize_pt >= minimum_text_pt),
  "The 7-pt text audit did not pass for all five Round 2G figures."
)
assert_true(
  nrow(figure_specs) == 5L &&
    all(figure_specs$width_mm == 174) &&
    all(figure_specs$png_dpi == 300),
  "All five figures must be 174 mm wide with 300-dpi PNG output."
)
assert_true(
  nrow(pdf_device_manifest) == 5L &&
    all(pdf_device_manifest$embedded_font_device),
  "All five PDFs must use a vector embedded-font device."
)
assert_true(
  nrow(colourbar_audit) == 6L &&
    n_distinct(colourbar_audit$figure) == 3L &&
    all(colourbar_audit$endpoint_and_label_order_agree) &&
    all(!colourbar_audit$guide_reverse) &&
    all(colourbar_audit$guide_direction == "vertical") &&
    all(colourbar_audit$guide_display == colourbar_display_mode),
  "The three corrected colourbars did not pass their endpoint audit."
)

write_csv_atomic(input_manifest, input_manifest_file)
write_csv_atomic(text_audit, text_audit_file)
write_csv_atomic(figure_specs, figure_spec_file)
write_csv_atomic(pdf_device_manifest, pdf_device_file)
write_csv_atomic(colourbar_audit, colourbar_audit_file)

expected_stems <- c(
  "Figure1_round2G_GTEx_v8_normal_tissue_atlas",
  "Figure2_round2G_count_level_paired_heatmap",
  "Figure4A_round2G_same_gene_recurrence_null",
  "SupplementaryFigureS1A_round2G_nonLIHC_coexpression_change_heatmap",
  "SupplementaryFigureS1B_round2G_representative_coexpression_rank_scatter"
)
expected_files <- unlist(lapply(
  expected_stems,
  function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".png")))
))
assert_true(
  length(expected_files) == 10L &&
    all(file.exists(expected_files)) &&
    all(file.info(expected_files)$size > 1000),
  "One or more Round 2G figure files are missing or unexpectedly small."
)

integrity_checks <- tibble::tribble(
  ~check, ~result, ~expected, ~passed, ~details,
  "frozen_inputs_unchanged",
  as.character(sum(input_md5_before == input_md5_after)),
  as.character(length(input_files)),
  identical(input_md5_before, input_md5_after),
  "Input MD5 before and after rendering",
  "only_five_targeted_figure_pairs",
  as.character(length(expected_files)),
  "10 files",
  length(expected_files) == 10L && all(file.exists(expected_files)),
  "Figures 3A/3B/4B are not regenerated",
  "minimum_text_size",
  formatC(min(text_audit$fontsize_pt), digits = 3, format = "f"),
  ">=7.000 pt",
  min(text_audit$fontsize_pt) >= 7,
  "Visible text grobs inspected at final dimensions",
  "colourbar_non_raster_contract",
  paste0("reverse=FALSE; direction=vertical; display=", colourbar_display_mode),
  "3 figures x 2 ends pass",
  nrow(colourbar_audit) == 6L &&
    n_distinct(colourbar_audit$figure) == 3L &&
    all(colourbar_audit$endpoint_and_label_order_agree) &&
    all(!colourbar_audit$guide_reverse),
  paste0(
    "Endpoint fill colours and corresponding extreme guide labels are in ",
    "round2G_colourbar_endpoint_audit.csv"
  ),
  "figure1_colourbar_semantics",
  paste(colourbar_audit$biological_direction[
    colourbar_audit$figure ==
      "Figure1_round2G_GTEx_v8_normal_tissue_atlas"
  ], collapse = "; "),
  "low/purple; high/yellow",
  identical(
    colourbar_audit$biological_direction[
      colourbar_audit$figure ==
        "Figure1_round2G_GTEx_v8_normal_tissue_atlas"
    ],
    c("low expression / dark purple", "high expression / yellow")
  ),
  "Tile scale unchanged; legend minimum is bottom and maximum is top",
  "figure2_colourbar_semantics",
  paste(colourbar_audit$biological_direction[
    colourbar_audit$figure ==
      "Figure2_round2G_count_level_paired_heatmap"
  ], collapse = "; "),
  "negative/blue; positive/red",
  identical(
    colourbar_audit$biological_direction[
      colourbar_audit$figure ==
        "Figure2_round2G_count_level_paired_heatmap"
    ],
    c(
      "negative tumor-minus-normal logFC / blue",
      "positive tumor-minus-normal logFC / red"
    )
  ),
  "Tile scale unchanged; legend minimum is bottom and maximum is top",
  "S1A_colourbar_semantics",
  paste(colourbar_audit$biological_direction[
    colourbar_audit$figure == paste0(
      "SupplementaryFigureS1A_round2G_",
      "nonLIHC_coexpression_change_heatmap"
    )
  ], collapse = "; "),
  "negative/purple; positive/green",
  identical(
    colourbar_audit$biological_direction[
      colourbar_audit$figure == paste0(
        "SupplementaryFigureS1A_round2G_",
        "nonLIHC_coexpression_change_heatmap"
      )
    ],
    c("negative delta rho / purple", "positive delta rho / green")
  ),
  "Tile scale unchanged; legend minimum is bottom and maximum is top",
  "figure4A_single_line_subtitle",
  as.character(!grepl("\n", figure4a_subtitle, fixed = TRUE)),
  "TRUE",
  !grepl("\n", figure4a_subtitle, fixed = TRUE),
  "Second sentence moved unchanged to the caption",
  "figure4A_annotation_order",
  "q y=1.105; profiles y=1.035",
  "q above profiles; both above panel",
  1.105 > 1.035 && 1.035 > 1,
  "Subtitle bottom margin=36 pt; q row lowered from Round 2F",
  "figure4A_q_labels",
  paste(observed_q_labels, collapse = "; "),
  "q=0.032; q=0.032; q=0.376",
  identical(observed_q_labels, c("q=0.032", "q=0.032", "q=0.376")),
  "Three-decimal display retained",
  "S1B_strip_labels",
  paste(panel_order, collapse = "; "),
  "3 single-line cancer/gene labels",
  length(panel_order) == 3L && all(!grepl("\n", panel_order, fixed = TRUE)),
  "Statistics removed from strips",
  "S1B_statistics_outside_data_range",
  paste(range(annotation_data$label_y), collapse = " to "),
  "all y<0 while rank data are within 0-1",
  all(annotation_data$label_y < 0) &&
    all(coexpression_points_plot$partner_rank_fraction >= 0) &&
    all(coexpression_points_plot$partner_rank_fraction <= 1),
  "Dedicated grey annotation band prevents overlap with data points",
  "S1B_values_unchanged",
  paste(annotation_data$statistic_label, collapse = "; "),
  paste(expected_stat_labels, collapse = "; "),
  identical(annotation_data$statistic_label, expected_stat_labels),
  "Values read from frozen Tables 181/182",
  "paper2_firewall",
  as.character(any(coexpression_cells$tcga_code == "LIHC")),
  "FALSE",
  !any(coexpression_cells$tcga_code == "LIHC"),
  "LIHC coexpression remains excluded"
)
assert_true(
  all(integrity_checks$passed),
  "One or more Round 2G integrity checks failed."
)
write_csv_atomic(integrity_checks, integrity_file)

completion_lines <- c(
  paste0("Analysis: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Project directory: ", project_dir),
  "Scope: five targeted figure corrections from frozen Round 2B values/caches",
  paste0(
    "Colourbars: reverse=FALSE; direction=vertical; display=",
    colourbar_display_mode
  ),
  paste0("Colourbar endpoint audit: PASS; ", colourbar_audit_file),
  "Figure 1: low=dark purple, high=yellow",
  "Figure 2: negative=blue, positive=red",
  "Supplementary Figure S1A: negative=purple, positive=green",
  "Figure 4A: one-line subtitle; second sentence moved to caption",
  paste0(
    "Figure 4A: q/profile rows separated; q labels=",
    paste(observed_q_labels, collapse = "; ")
  ),
  "Supplementary Figure S1B: statistics removed from facet strips",
  "Supplementary Figure S1B: rho/delta shown in a dedicated y<0 annotation band",
  paste0(
    "Final-size text audit: PASS; minimum=",
    formatC(min(text_audit$fontsize_pt), digits = 3, format = "f"),
    " pt"
  ),
  "All five figures: width=174 mm; PNG=300 dpi; vector PDF",
  paste0("Input files unchanged by MD5: ", length(input_files)),
  paste0("Figures: ", figure_dir),
  paste0("Logs: ", log_dir)
)
write_lines_atomic(completion_lines, completion_report_file)

log_message("Round 2G complete.")
log_message("Figures: ", figure_dir)
log_message("Integrity checks: ", integrity_file)
log_message("Session information will be written last: ", session_info_file)
write_lines_atomic(capture.output(sessionInfo()), session_info_file)
