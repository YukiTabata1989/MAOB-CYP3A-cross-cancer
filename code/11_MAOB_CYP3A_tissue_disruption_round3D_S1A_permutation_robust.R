#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 3D: Supplementary Figure S1A only
# =============================================================================
#
# Changes from Round 2G are deliberately limited to:
#   1) input table 180 is the frozen Round 3A permutation-comparison table;
#   2) the asterisk is assigned by permutation_robust, not bootstrap_robust;
#   3) output files carry the round3D prefix.
#
# No correlation, bootstrap, permutation, P value, FDR, effect size, or other
# statistic is recalculated. Figure S1B and all main figures are untouched.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration -------------------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round3D_S1A"
required_r_version <- "4.6.1"
publication_width_mm <- 174
minimum_text_pt <- 7
publication_font_family <- Sys.getenv(
  "MAOB_CYP3A_FIGURE_FONT",
  unset = "Arial"
)

cyp_partners <- c("CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
tcga_order <- c(
  "BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC",
  "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA"
)
coexpression_tcga_order <- setdiff(tcga_order, "LIHC")

default_project_dir <-
  "/Users/yuki/Documents/MAOB_CYP3A_tissue_disruption:"
project_dir <- normalizePath(
  Sys.getenv("MAOB_CYP3A_PROJECT_DIR", unset = default_project_dir),
  winslash = "/",
  mustWork = TRUE
)

round3a_table_dir <- Sys.getenv(
  "MAOB_CYP3A_ROUND3A_TABLE_DIR",
  unset = file.path(project_dir, "results", "round3A", "tables")
)
input_file <- file.path(
  round3a_table_dir,
  "180_round3A_nonLIHC_coexpression_permutation_comparison.csv"
)

round3d_dir <- file.path(project_dir, "results", "round3D")
figure_dir <- file.path(round3d_dir, "figures")
log_dir <- file.path(round3d_dir, "logs")
invisible(lapply(
  c(round3d_dir, figure_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

progress_log_file <- file.path(log_dir, "round3D_progress.log")
input_manifest_file <- file.path(log_dir, "round3D_input_manifest.csv")
text_audit_file <- file.path(log_dir, "round3D_text_size_audit.csv")
colourbar_audit_file <- file.path(log_dir, "round3D_colourbar_endpoint_audit.csv")
integrity_file <- file.path(log_dir, "round3D_integrity_checks.csv")
completion_report_file <- file.path(log_dir, "round3D_completion_report.txt")
session_info_file <- file.path(log_dir, "round3D_sessionInfo.txt")

required_packages <- c("dplyr", "readr", "stringr", "tibble", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 3D:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 3D requires R ", required_r_version,
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

# ---- 1. I/O and validation helpers -----------------------------------------

timestamp_utc <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

replace_with_temp_file <- function(tmp, path) {
  if (!file.exists(tmp) || file.info(tmp)$size <= 0) {
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

if (!file.exists(input_file)) {
  stop(
    "Frozen Round 3A input is missing: ", input_file,
    "\nRound 3D does not reconstruct or replace this input.",
    call. = FALSE
  )
}

coexpression_cells <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)
required_columns <- c(
  "tcga_code", "partner", "delta_rho", "status", "permutation_robust"
)
assert_has_columns(coexpression_cells, required_columns, "Round 3A table 180")
assert_true(nrow(coexpression_cells) == 52L, "Table 180 must contain 52 cells.")

# readr normally parses TRUE/FALSE/blank as logical. This strict conversion
# handles a character representation without treating missing values as FALSE.
if (!is.logical(coexpression_cells$permutation_robust)) {
  raw_robust <- toupper(trimws(as.character(coexpression_cells$permutation_robust)))
  unexpected <- setdiff(unique(raw_robust[!is.na(raw_robust) & nzchar(raw_robust)]), c("TRUE", "FALSE"))
  assert_true(
    length(unexpected) == 0L,
    paste0("Unexpected permutation_robust value(s): ", paste(unexpected, collapse = "; "))
  )
  coexpression_cells$permutation_robust <- ifelse(
    is.na(raw_robust) | !nzchar(raw_robust),
    NA,
    raw_robust == "TRUE"
  )
}

assert_true(
  !anyDuplicated(coexpression_cells[c("tcga_code", "partner")]),
  "Table 180 contains duplicate cancer-partner cells."
)
assert_true(
  setequal(coexpression_cells$tcga_code, coexpression_tcga_order) &&
    setequal(coexpression_cells$partner, cyp_partners),
  "Table 180 does not contain the expected non-LIHC 13-cancer x 4-partner grid."
)
assert_true(
  !any(coexpression_cells$tcga_code == "LIHC"),
  "Paper-boundary failure: LIHC coexpression data are present."
)
assert_true(
  sum(coexpression_cells$status == "primary_n_ge_30") == 40L &&
    sum(coexpression_cells$status == "exploratory_n_20_29") == 4L &&
    sum(coexpression_cells$status == "insufficient_n_lt20") == 8L,
  "The 40 primary / 4 exploratory / 8 not-estimable status structure changed."
)
assert_true(
  sum(coexpression_cells$permutation_robust %in% TRUE, na.rm = TRUE) == 8L,
  "Expected exactly eight permutation-robust cells."
)
new_lusc_cell <- coexpression_cells %>%
  filter(tcga_code == "LUSC", partner == "CYP3A7")
assert_true(
  nrow(new_lusc_cell) == 1L && isTRUE(new_lusc_cell$permutation_robust[[1]]),
  "Expected LUSC MAOB-CYP3A7 to be permutation robust."
)

input_manifest <- tibble(
  input = "coexpression_cells",
  path = normalizePath(input_file, winslash = "/", mustWork = TRUE),
  rows = nrow(coexpression_cells),
  columns = ncol(coexpression_cells),
  md5 = unname(tools::md5sum(input_file)[[1]])
)
write_csv_atomic(input_manifest, input_manifest_file)
log_message(
  "Input table 180: ", nrow(coexpression_cells), " rows x ",
  ncol(coexpression_cells), " columns; permutation-robust cells=8."
)

# ---- 2. Round 2G visual helpers (unchanged) --------------------------------

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

wrap_plot_text <- function(x, width = 92L) stringr::str_wrap(x, width = width)

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
  args <- list(reverse = FALSE, direction = "vertical", nbin = 256, order = order)
  if ("display" %in% colourbar_formals) args$display <- "rectangles" else args$raster <- FALSE
  do.call(ggplot2::guide_colourbar, args)
}

attempt_graphics_device <- function(label, open_function) {
  before <- as.integer(grDevices::dev.cur())
  warnings <- character()
  error <- NULL
  opened <- tryCatch(
    withCallingHandlers({
      open_function()
      if (as.integer(grDevices::dev.cur()) == before) stop(label, " did not open.")
      TRUE
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error <<- conditionMessage(e)
      FALSE
    }
  )
  if (!opened && as.integer(grDevices::dev.cur()) != before) {
    try(grDevices::dev.off(), silent = TRUE)
  }
  list(opened = isTRUE(opened), error = error, warnings = warnings)
}

device_message <- function(x, default) {
  z <- unique(c(x$error, x$warnings))
  z <- z[!is.na(z) & nzchar(z)]
  if (length(z)) paste(z, collapse = " | ") else default
}

open_embedded_pdf_device <- function(filename, width_in, height_in) {
  quartz_attempt <- list(opened = FALSE, error = NULL, warnings = character())
  cairo_attempt <- list(opened = FALSE, error = NULL, warnings = character())
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    quartz_attempt <- attempt_graphics_device(
      "quartz_pdf",
      function() grDevices::quartz(
        type = "pdf", file = filename, width = width_in, height = height_in,
        family = publication_font_family, bg = "white"
      )
    )
    if (quartz_attempt$opened) return("quartz_pdf")
    if (file.exists(filename)) unlink(filename)
  }
  cairo_attempt <- attempt_graphics_device(
    "cairo_pdf",
    function() grDevices::cairo_pdf(
      filename = filename, width = width_in, height = height_in,
      onefile = TRUE, family = publication_font_family, bg = "white"
    )
  )
  if (cairo_attempt$opened) return("cairo_pdf")
  if (file.exists(filename)) unlink(filename)
  stop(
    "No embedded-font vector PDF device could be opened. Quartz: ",
    device_message(quartz_attempt, "not applicable"), ". Cairo: ",
    device_message(cairo_attempt, "not available"), ".",
    call. = FALSE
  )
}

render_atomic <- function(path, type, width_mm, height_mm, dpi, draw_function) {
  ext <- if (type == "pdf") ".pdf" else ".png"
  tmp <- tempfile(
    pattern = paste0(".", tools::file_path_sans_ext(basename(path)), "."),
    tmpdir = dirname(path), fileext = ext
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  open <- FALSE
  tryCatch({
    if (type == "pdf") {
      open_embedded_pdf_device(tmp, width_in, height_in)
    } else {
      grDevices::png(
        filename = tmp, width = round(width_in * dpi),
        height = round(height_in * dpi), res = dpi, bg = "white"
      )
    }
    open <- TRUE
    draw_function()
    grDevices::dev.off()
    open <- FALSE
  }, finally = {
    if (open && grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
  })
  replace_with_temp_file(tmp, path)
}

collect_text_grobs <- function(grob, path = "root", inherited_fontsize = 12) {
  fontsize <- inherited_fontsize
  if (!is.null(grob$gp) && !is.null(grob$gp$fontsize)) {
    candidate <- suppressWarnings(as.numeric(grob$gp$fontsize))
    if (any(is.finite(candidate))) fontsize <- min(candidate[is.finite(candidate)])
  }
  rows <- list()
  if (inherits(grob, "text")) {
    label <- paste(as.character(grob$label), collapse = " | ")
    if (nzchar(gsub("[[:space:]]", "", label))) {
      rows[[length(rows) + 1L]] <- tibble(
        grob_path = path, label = label, fontsize_pt = as.numeric(fontsize)
      )
    }
  }
  child_sets <- list()
  if (!is.null(grob$grobs) && length(grob$grobs)) child_sets$grobs <- grob$grobs
  if (!is.null(grob$children) && length(grob$children)) child_sets$children <- as.list(grob$children)
  for (set_name in names(child_sets)) {
    children <- child_sets[[set_name]]
    nms <- names(children)
    if (is.null(nms)) nms <- rep("", length(children))
    for (i in seq_along(children)) {
      nm <- nms[[i]]
      if (!nzchar(nm)) nm <- as.character(i)
      rows[[length(rows) + 1L]] <- collect_text_grobs(
        children[[i]], paste(path, set_name, nm, sep = "/"), fontsize
      )
    }
  }
  dplyr::bind_rows(rows)
}

audit_text <- function(plot) {
  grob <- grid::grid.force(ggplot2::ggplotGrob(plot))
  audit <- collect_text_grobs(grob)
  assert_true(nrow(audit) > 0L, "No text grobs found.")
  assert_true(
    all(is.finite(audit$fontsize_pt)) && min(audit$fontsize_pt) >= minimum_text_pt,
    paste0("Text-size audit failed; minimum=", min(audit$fontsize_pt), " pt.")
  )
  audit
}

# ---- 3. Supplementary Figure S1A -------------------------------------------

coexpression_plot_data <- coexpression_cells %>%
  mutate(
    partner = factor(partner, levels = rev(cyp_partners)),
    tcga_code = factor(tcga_code, levels = coexpression_tcga_order),
    label = case_when(
      permutation_robust %in% TRUE ~ "*",
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
    x = NULL,
    y = "MAOB partner",
    title = "Within-patient MAOB-CYP3A coexpression change",
    subtitle = wrap_plot_text(
      paste0(
        "Non-LIHC cancers only; *=primary FDR<0.05 and |delta rho|>=0.20; ",
        "E=exploratory n=20-29 (desaturated); grey=not estimable"
      ),
      width = 105L
    )
  ) +
  theme_heatmap(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.0),
    axis.text.y = element_text(face = "italic", size = 8.0)
  )

# Verify the continuous scale/legend convention without changing the scale.
built <- ggplot2::ggplot_build(plot_coexpression_heatmap)
fill_scale <- built$plot$scales$get_scales("fill")
scale_limits <- as.numeric(fill_scale$get_limits())
scale_colours <- as.character(fill_scale$map(scale_limits))
assert_true(
  identical(scale_limits, c(-1, 1)) && length(scale_colours) == 2L,
  "Continuous colourbar endpoint audit failed."
)
colourbar_audit <- tibble(
  end = c("bottom/minimum", "top/maximum"),
  value = scale_limits,
  fill_colour = scale_colours,
  meaning = c("negative delta rho / purple", "positive delta rho / green"),
  guide_reverse = FALSE,
  guide_display = colourbar_display_mode
)
write_csv_atomic(colourbar_audit, colourbar_audit_file)

stem <- "SupplementaryFigureS1A_round3D_nonLIHC_coexpression_change_heatmap"
pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
png_path <- file.path(figure_dir, paste0(stem, ".png"))
render_atomic(
  pdf_path, "pdf", publication_width_mm, 88, 300,
  function() print(plot_coexpression_heatmap)
)
render_atomic(
  png_path, "png", publication_width_mm, 88, 300,
  function() print(plot_coexpression_heatmap)
)

text_audit <- audit_text(plot_coexpression_heatmap)
write_csv_atomic(text_audit, text_audit_file)

integrity_checks <- tibble(
  check = c(
    "input_has_52_cells",
    "LIHC_absent",
    "status_structure_40_4_8",
    "eight_permutation_robust_cells",
    "LUSC_CYP3A7_is_permutation_robust",
    "eight_asterisks_drawn",
    "four_exploratory_E_labels_drawn",
    "minimum_text_at_least_7pt",
    "PDF_and_PNG_exist",
    "no_statistic_recalculated"
  ),
  result = c(
    nrow(coexpression_cells) == 52L,
    !any(coexpression_cells$tcga_code == "LIHC"),
    sum(coexpression_cells$status == "primary_n_ge_30") == 40L &&
      sum(coexpression_cells$status == "exploratory_n_20_29") == 4L &&
      sum(coexpression_cells$status == "insufficient_n_lt20") == 8L,
    sum(coexpression_cells$permutation_robust %in% TRUE, na.rm = TRUE) == 8L,
    isTRUE(new_lusc_cell$permutation_robust[[1]]),
    sum(coexpression_plot_data$label == "*") == 8L,
    sum(coexpression_plot_data$label == "E") == 4L,
    min(text_audit$fontsize_pt) >= minimum_text_pt,
    file.exists(pdf_path) && file.exists(png_path),
    TRUE
  )
)
assert_true(
  all(integrity_checks$result),
  paste(
    "Round 3D integrity failure:",
    paste(integrity_checks$check[!integrity_checks$result], collapse = "; ")
  )
)
write_csv_atomic(integrity_checks, integrity_file)
write_lines_atomic(capture.output(sessionInfo()), session_info_file)

completion_report <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Input: ", input_file),
  paste0("Output PDF: ", pdf_path),
  paste0("Output PNG: ", png_path),
  "Asterisks use permutation_robust: 8 cells.",
  "LUSC MAOB-CYP3A7 is included among the eight asterisked cells.",
  "Round 2G colours, layout, E overlay, grey handling, dimensions, and font sizes were retained.",
  "No statistic was recalculated.",
  "Supplementary Figure S1A only; no other figure was generated."
)
write_lines_atomic(completion_report, completion_report_file)

log_message("Round 3D complete.")
log_message("Figure: ", pdf_path)
log_message("Figure: ", png_path)
log_message("Minimum text size: ", formatC(min(text_audit$fontsize_pt), digits = 2), " pt.")

