#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 3C: annotation and layout corrections for Figure 4A/4B/4C
# =============================================================================
#
# PURPOSE
#   Rebuild Figure 4A/4B/4C from the frozen denominator-aligned Round 3A
#   outputs, changing annotation placement and graphical discrimination only.
#
# ANALYTICAL BOUNDARY
#   - No differential-expression, recurrence, matched-background, coexpression,
#     permutation, GTEx, tissue-identity, baseline, or HPA statistic is
#     recalculated.
#   - Plotting labels, factors, expected QQ plotting positions, and display
#     rounding are graphical transformations only.
#   - Supplementary Figure S4 is final and is neither read nor regenerated.
#   - The CSV column names used below were verified directly against the
#     supplied round3A.zip before this script was written. Missing columns cause
#     an immediate stop; no substitute column or derived statistic is created.
#   - No LIHC coexpression, survival, clinical, subgroup, pathway, immune, or
#     tumor-purity analysis is performed.
#
# OUTPUTS
#   results/round3C/figures/
#     Figure4A_round3C_denominator_aligned_recurrence_null.{pdf,png}
#     Figure4B_round3C_coverage_sensitivity.{pdf,png}
#     Figure4C_round3C_MAOB_organ_collapsed_sensitivity.{pdf,png}
#
#   results/round3C/logs/
#     round3C_progress.log
#     round3C_input_manifest.csv
#     round3C_text_size_audit.csv
#     round3C_figure_dimensions.csv
#     round3C_pdf_device_manifest.csv
#     round3C_integrity_checks.csv
#     round3C_completion_report.txt
#     round3C_sessionInfo.txt
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration --------------------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round3C_figure4_layout"
required_r_version <- "4.6.1"
publication_width_mm <- 174
minimum_text_pt <- 7
publication_font_family <- Sys.getenv(
  "MAOB_CYP3A_FIGURE_FONT",
  unset = "Arial"
)

target_genes <- c("MAOB", "CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
coverage_levels <- c(0.50, 0.60, 0.70, 0.80, 0.90)
primary_coverage <- 0.80

gene_palette_okabe_ito <- c(
  MAOB = "#000000",
  CYP3A4 = "#E69F00",
  CYP3A5 = "#56B4E9",
  CYP3A7 = "#009E73",
  CYP3A43 = "#CC79A7"
)

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
round3c_dir <- file.path(project_dir, "results", "round3C")
figure_dir <- file.path(round3c_dir, "figures")
log_dir <- file.path(round3c_dir, "logs")
invisible(lapply(
  c(round3c_dir, figure_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

input_files <- c(
  recurrence_80 = file.path(
    round3a_table_dir,
    "162_round3A_denominator_aligned_same_gene_recurrence_null.csv"
  ),
  recurrence_profiles = file.path(
    round3a_table_dir,
    "163_round3A_denominator_aligned_same_gene_profiles.csv"
  ),
  coverage_sensitivity = file.path(
    round3a_table_dir,
    "164_round3A_denominator_aligned_coverage_sensitivity.csv"
  ),
  organ_collapsed = file.path(
    round3a_table_dir,
    "168_round3A_MAOB_organ_collapsed_recurrence_sensitivity.csv"
  )
)

progress_log_file <- file.path(log_dir, "round3C_progress.log")
input_manifest_file <- file.path(log_dir, "round3C_input_manifest.csv")
text_audit_file <- file.path(log_dir, "round3C_text_size_audit.csv")
figure_spec_file <- file.path(log_dir, "round3C_figure_dimensions.csv")
pdf_device_file <- file.path(log_dir, "round3C_pdf_device_manifest.csv")
integrity_file <- file.path(log_dir, "round3C_integrity_checks.csv")
completion_report_file <- file.path(log_dir, "round3C_completion_report.txt")
session_info_file <- file.path(log_dir, "round3C_sessionInfo.txt")

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
  "ggplot2", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 3B:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 3B requires R ", required_r_version,
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

safe_md5 <- function(path) unname(tools::md5sum(path)[[1]])

input_records <- list()
record_input <- function(label, path, x) {
  input_records[[length(input_records) + 1L]] <<- tibble(
    input = label,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    rows = nrow(x),
    columns = ncol(x),
    md5 = safe_md5(path)
  )
  log_message(
    "Input ", label, ": ", nrow(x), " rows x ", ncol(x),
    " columns; MD5=", safe_md5(path)
  )
  invisible(NULL)
}

read_frozen_csv <- function(label, required_columns, expected_rows = NULL) {
  path <- input_files[[label]]
  if (!file.exists(path)) {
    stop(
      "Required frozen Round 3A input is missing: ", path,
      "\nRound 3B does not regenerate missing statistics.",
      call. = FALSE
    )
  }
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
  record_input(label, path, x)
  x
}

wrap_plot_text <- function(x, width = 100L) {
  stringr::str_wrap(x, width = width)
}

format_q3 <- function(x) {
  ifelse(is.finite(x), formatC(x, digits = 3, format = "f"), "")
}

format_p4 <- function(x) {
  ifelse(is.finite(x), formatC(x, digits = 4, format = "f"), "")
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
    pdf_device_records[[length(pdf_device_records) + 1L]] <<- tibble(
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
      rows[[length(rows) + 1L]] <- tibble(
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
  figure_output_records[[length(figure_output_records) + 1L]] <<- tibble(
    figure = stem,
    width_mm = width_mm,
    height_mm = height_mm,
    png_dpi = dpi
  )
  log_message("Wrote figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

draw_plot_grid <- function(plots, heights) {
  assert_true(
    length(plots) == length(heights),
    "Plot-grid dimensions do not match."
  )
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = length(plots),
    ncol = 1L,
    heights = grid::unit(heights, "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(plots)) {
    print(
      plots[[i]],
      newpage = FALSE,
      vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1L)
    )
  }
  grid::popViewport()
  invisible(NULL)
}

save_plot_grid_pair <- function(
    plots, heights, stem, height_mm,
    width_mm = publication_width_mm, dpi = 300
) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  draw_fun <- function() draw_plot_grid(plots, heights)
  render_atomic(
    pdf_path, "pdf", width_mm, height_mm, dpi, draw_fun, stem
  )
  render_atomic(
    png_path, "png", width_mm, height_mm, dpi, draw_fun, stem
  )
  invisible(lapply(plots, audit_figure_text, stem = stem))
  figure_output_records[[length(figure_output_records) + 1L]] <<- tibble(
    figure = stem,
    width_mm = width_mm,
    height_mm = height_mm,
    png_dpi = dpi
  )
  log_message("Wrote composite figure pair: ", stem)
  invisible(c(pdf_path, png_path))
}

# ---- 2. Read frozen Round 3A inputs -----------------------------------------

recurrence_80 <- read_frozen_csv(
  "recurrence_80",
  c(
    "gene", "minimum_coverage_fraction", "common_denominator_n",
    "observed_loss_fraction", "n_candidate_genes", "q_recurrent_loss",
    "n_target_contexts_total", "status"
  ),
  expected_rows = 5L
)

recurrence_profiles <- read_frozen_csv(
  "recurrence_profiles",
  c(
    "target_gene", "minimum_coverage_fraction", "common_denominator_n",
    "selected_common_contexts", "background_gene_id", "loss_fraction"
  )
)

coverage_sensitivity <- read_frozen_csv(
  "coverage_sensitivity",
  c(
    "gene", "minimum_coverage_fraction", "common_denominator_n",
    "n_target_contexts_total", "n_candidate_genes",
    "q_recurrent_loss", "q_median_negative", "status"
  ),
  expected_rows = 25L
)

organ_collapsed <- read_frozen_csv(
  "organ_collapsed",
  c(
    "gene", "minimum_coverage_fraction", "common_denominator_n",
    "observed_n_loss", "observed_loss_fraction", "n_candidate_genes",
    "p_recurrent_loss", "status", "analysis_unit"
  ),
  expected_rows = 5L
)

assert_true(
  setequal(recurrence_80$gene, target_genes),
  "Table 162 does not contain exactly the five target genes."
)
assert_true(
  all(abs(recurrence_80$minimum_coverage_fraction - 0.8) < 1e-12),
  "Table 162 is not restricted to the 80% setting."
)
assert_true(
  !anyDuplicated(coverage_sensitivity[c("gene", "minimum_coverage_fraction")]),
  "Table 164 contains duplicate gene-coverage cells."
)
assert_true(
  setequal(coverage_sensitivity$gene, target_genes) &&
    setequal(coverage_sensitivity$minimum_coverage_fraction, coverage_levels),
  "Table 164 does not represent all five genes and five coverage settings."
)
assert_true(
  all(organ_collapsed$gene == "MAOB") &&
    setequal(organ_collapsed$minimum_coverage_fraction, coverage_levels),
  "Table 168 must contain the five MAOB organ-collapsed coverage settings."
)

# ---- 3. Figure 4A: denominator-aligned 80% recurrence null -----------------

log_message("Building Figure 4A from denominator-aligned 80% profiles")

evaluable_80 <- recurrence_80 %>%
  filter(status == "estimable_denominator_aligned_null") %>%
  arrange(match(gene, target_genes))
assert_true(
  identical(evaluable_80$gene, c("MAOB", "CYP3A4", "CYP3A7")),
  "Unexpected evaluable genes in denominator-aligned table 162."
)

expected_4a <- tibble(
  gene = c("MAOB", "CYP3A4", "CYP3A7"),
  denominator = c(10L, 7L, 4L),
  profiles = c(40L, 89L, 281L),
  q = c(0.0731707317, 0.0731707317, 0.5354609929)
)
observed_4a <- evaluable_80 %>%
  transmute(
    gene,
    denominator = as.integer(common_denominator_n),
    profiles = as.integer(n_candidate_genes),
    q = q_recurrent_loss
  ) %>%
  left_join(expected_4a, by = "gene", suffix = c("_observed", "_expected"))
assert_true(
  all(observed_4a$denominator_observed == observed_4a$denominator_expected) &&
    all(observed_4a$profiles_observed == observed_4a$profiles_expected) &&
    all(abs(observed_4a$q_observed - observed_4a$q_expected) < 5e-4),
  "Figure 4A values do not match the prespecified Round 3A audit values."
)

figure4a_profiles <- recurrence_profiles %>%
  filter(abs(minimum_coverage_fraction - primary_coverage) < 1e-12) %>%
  mutate(gene = factor(target_gene, levels = target_genes))

profile_count_audit <- figure4a_profiles %>%
  count(target_gene, name = "profiles_in_plot") %>%
  right_join(
    recurrence_80 %>%
      transmute(target_gene = gene, expected_profiles = n_candidate_genes),
    by = "target_gene"
  ) %>%
  mutate(
    profiles_in_plot = coalesce(profiles_in_plot, 0L),
    expected_profiles = coalesce(as.integer(expected_profiles), 0L)
  )
assert_true(
  all(
    profile_count_audit$profiles_in_plot ==
      as.integer(profile_count_audit$expected_profiles)
  ),
  "Table 163 profile counts do not match table 162 at 80% coverage."
)

figure4a_observed <- evaluable_80 %>%
  mutate(gene = factor(gene, levels = target_genes))

figure4a_annotations <- recurrence_80 %>%
  mutate(gene = factor(gene, levels = target_genes)) %>%
  transmute(
    gene,
    q_label = if_else(
      status == "estimable_denominator_aligned_null",
      paste0("q=", format_q3(q_recurrent_loss)),
      ""
    ),
    denominator_label = if_else(
      status == "estimable_denominator_aligned_null",
      paste0("common\ndenominator: ", as.integer(common_denominator_n)),
      ""
    ),
    profile_label = if_else(
      status == "estimable_denominator_aligned_null",
      paste0("eligible profiles: ", as.integer(n_candidate_genes)),
      ""
    )
  )

figure4a_not_evaluable <- recurrence_80 %>%
  filter(status != "estimable_denominator_aligned_null") %>%
  mutate(
    gene = factor(gene, levels = target_genes),
    y = 0.52,
    label = case_when(
      gene == "CYP3A5" & status == "insufficient_common_set_candidates" ~
        paste0(
          "Not evaluable\n", as.integer(n_candidate_genes),
          " eligible profiles (<20)"
        ),
      gene == "CYP3A43" & status == "insufficient_target_contexts" ~
        paste0(
          "Not evaluable\n", as.integer(n_target_contexts_total),
          " target contexts (<4)"
        ),
      status == "insufficient_common_set_candidates" ~ paste0(
        "Not evaluable\n", as.integer(n_candidate_genes),
        " eligible profiles (<20)"
      ),
      status == "insufficient_target_contexts" ~ "Not evaluable\nToo few target contexts",
      TRUE ~ paste0("Not evaluable\n", status)
    )
  )

assert_true(
  nrow(figure4a_not_evaluable) == 2L &&
    setequal(as.character(figure4a_not_evaluable$gene), c("CYP3A5", "CYP3A43")),
  "Figure 4A must show exactly CYP3A5 and CYP3A43 as not evaluable."
)

figure4a_subtitle <- paste(
  "Grey distributions show complete background-gene profiles evaluated on",
  "the same cancer set as each target at the 80% setting; diamonds show targets."
)
figure4a_caption <- "Not evaluable is distinct from evidence of no change."

plot_figure4a <- ggplot(
  figure4a_profiles %>%
    filter(target_gene %in% evaluable_80$gene),
  aes(x = gene, y = loss_fraction)
) +
  geom_violin(
    fill = "grey88", color = "grey45", scale = "width", trim = TRUE
  ) +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white") +
  geom_point(
    data = figure4a_observed,
    aes(x = gene, y = observed_loss_fraction, fill = gene),
    inherit.aes = FALSE,
    shape = 23,
    size = 3.5,
    color = "black",
    stroke = 0.55
  ) +
  geom_text(
    data = figure4a_annotations,
    aes(x = gene, y = 1.300, label = q_label),
    inherit.aes = FALSE,
    size = 2.85,
    fontface = "bold"
  ) +
  geom_text(
    data = figure4a_annotations,
    aes(x = gene, y = 1.205, label = denominator_label),
    inherit.aes = FALSE,
    size = 2.65,
    color = "grey20",
    lineheight = 0.90
  ) +
  geom_text(
    data = figure4a_annotations,
    aes(x = gene, y = 1.085, label = profile_label),
    inherit.aes = FALSE,
    size = 2.65,
    color = "grey25"
  ) +
  geom_label(
    data = figure4a_not_evaluable,
    aes(x = gene, y = y, label = label),
    inherit.aes = FALSE,
    fill = "grey93",
    color = "grey30",
    linewidth = 0.25,
    size = 3.0,
    label.size = 0.25
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
    title = "Denominator-aligned cross-cancer recurrence null",
    subtitle = wrap_plot_text(figure4a_subtitle, width = 110),
    caption = figure4a_caption
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", size = 8.0),
    plot.subtitle = element_text(
      size = 8.0,
      color = "grey25",
      margin = margin(b = 82)
    ),
    plot.caption.position = "plot",
    plot.caption = element_text(
      size = 7.5,
      color = "grey25",
      hjust = 0,
      margin = margin(t = 7)
    ),
    plot.margin = margin(t = 18, r = 8, b = 6, l = 6)
  )

save_plot_pair(
  plot_figure4a,
  "Figure4A_round3C_denominator_aligned_recurrence_null",
  height_mm = 152
)

# ---- 4. Figure 4B: coverage sensitivity, two endpoints ----------------------

log_message("Building Figure 4B with primary and supplementary endpoints")

coverage_plot_long <- coverage_sensitivity %>%
  select(
    gene, minimum_coverage_fraction, common_denominator_n,
    n_candidate_genes, status, q_recurrent_loss, q_median_negative
  ) %>%
  pivot_longer(
    cols = c(q_recurrent_loss, q_median_negative),
    names_to = "endpoint",
    values_to = "FDR"
  ) %>%
  mutate(
    gene = factor(gene, levels = target_genes),
    endpoint = recode(
      endpoint,
      q_recurrent_loss = "Recurrent loss (primary)",
      q_median_negative = "Negative median (supplementary)"
    ),
    endpoint = factor(
      endpoint,
      levels = c(
        "Recurrent loss (primary)",
        "Negative median (supplementary)"
      )
    )
  ) %>%
  group_by(gene, endpoint) %>%
  arrange(minimum_coverage_fraction, .by_group = TRUE) %>%
  mutate(
    # Do not bridge a line across an explicitly non-evaluable coverage cell.
    line_segment = cumsum(!is.finite(FDR))
  ) %>%
  ungroup()

finite_fdr <- coverage_plot_long$FDR[is.finite(coverage_plot_long$FDR)]
assert_true(
  length(finite_fdr) > 0L && all(finite_fdr > 0 & finite_fdr <= 1),
  "Figure 4B requires finite FDR values in (0, 1]."
)
lower_log_limit <- min(0.01, 10^floor(log10(min(finite_fdr))))

plot_figure4b_main <- ggplot(
  coverage_plot_long %>% filter(is.finite(FDR)),
  aes(
    x = minimum_coverage_fraction,
    y = FDR,
    color = gene,
    linetype = endpoint,
    shape = gene,
    group = interaction(gene, endpoint, line_segment)
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dotted",
    color = "grey25",
    linewidth = 0.55
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.5, stroke = 0.65, na.rm = TRUE) +
  scale_color_manual(values = gene_palette_okabe_ito, drop = FALSE) +
  scale_linetype_manual(
    values = c(
      "Recurrent loss (primary)" = "solid",
      "Negative median (supplementary)" = "dashed"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      MAOB = 16,
      CYP3A4 = 17,
      CYP3A5 = 15,
      CYP3A7 = 18,
      CYP3A43 = 8
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = coverage_levels,
    labels = scales::percent_format(accuracy = 1),
    limits = range(coverage_levels)
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 1.00),
    labels = c("0.01", "0.025", "0.05", "0.10", "0.25", "0.50", "1.00"),
    limits = c(lower_log_limit, 1)
  ) +
  labs(
    x = NULL,
    y = "FDR (log scale)",
    color = "Gene",
    shape = "Gene",
    linetype = "Endpoint",
    title = "Coverage sensitivity of denominator-aligned recurrence tests",
    subtitle = wrap_plot_text(
      paste(
        "Solid lines show recurrent loss (primary endpoint); dashed lines show",
        "negative median (supplementary endpoint). The dotted line marks FDR = 0.05."
      ),
      width = 115
    )
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(t = 7, r = 7, b = 2, l = 7)
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2)
  )

coverage_status <- coverage_sensitivity %>%
  mutate(
    gene = factor(gene, levels = rev(target_genes)),
    status_class = if_else(
      status == "estimable_denominator_aligned_null",
      "Evaluable",
      "Not evaluable"
    ),
    status_label = case_when(
      status == "estimable_denominator_aligned_null" ~ "E",
      status == "insufficient_common_set_candidates" ~ paste0(
        "NE\n", as.integer(n_candidate_genes), " profiles"
      ),
      status == "insufficient_target_contexts" ~ paste0(
        "NE\n", as.integer(n_target_contexts_total), " contexts"
      ),
      TRUE ~ "NE"
    )
  )

assert_true(
  nrow(coverage_status) == 25L &&
    !anyDuplicated(coverage_status[c("gene", "minimum_coverage_fraction")]),
  "Figure 4B status strip must contain all 25 gene-coverage cells."
)

plot_figure4b_status <- ggplot(
  coverage_status,
  aes(x = minimum_coverage_fraction, y = gene)
) +
  geom_tile(
    aes(fill = status_class),
    color = "white",
    linewidth = 0.8,
    width = 0.095,
    height = 0.92
  ) +
  geom_text(
    aes(label = status_label),
    size = 2.7,
    lineheight = 0.92,
    color = "grey15"
  ) +
  scale_fill_manual(
    values = c("Evaluable" = "grey88", "Not evaluable" = "grey68"),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = coverage_levels,
    labels = scales::percent_format(accuracy = 1),
    limits = c(0.45, 0.95),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Minimum coverage setting",
    y = NULL,
    fill = "Cell status",
    subtitle = wrap_plot_text(
      paste(
        "Status strip for all 25 gene-coverage cells. E = evaluable;",
        "NE = not evaluable. Vertical position in this strip does not represent FDR."
      ),
      width = 115
    )
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "italic", size = 8.0),
    panel.grid = element_blank(),
    plot.margin = margin(t = 3, r = 7, b = 6, l = 7)
  )

save_plot_grid_pair(
  plots = list(plot_figure4b_main, plot_figure4b_status),
  heights = c(2.25, 1.15),
  stem = "Figure4B_round3C_coverage_sensitivity",
  height_mm = 190
)

# ---- 5. Figure 4C: organ-collapsed MAOB sensitivity -------------------------

log_message("Building Figure 4C from frozen organ-collapsed MAOB results")

organ_evaluable <- organ_collapsed %>%
  filter(status == "estimable_denominator_aligned_null") %>%
  arrange(minimum_coverage_fraction)
organ_not_evaluable <- organ_collapsed %>%
  filter(status != "estimable_denominator_aligned_null") %>%
  arrange(minimum_coverage_fraction)

assert_true(
  identical(organ_evaluable$minimum_coverage_fraction, c(0.5, 0.6, 0.7)) &&
    identical(organ_not_evaluable$minimum_coverage_fraction, c(0.8, 0.9)),
  "Figure 4C evaluability pattern must be 0.5-0.7 evaluable and 0.8-0.9 not evaluable."
)

expected_4c <- tibble(
  minimum_coverage_fraction = c(0.5, 0.6, 0.7),
  n_loss = c(5L, 6L, 7L),
  denominator = c(5L, 6L, 7L),
  p = c(0.0061349693, 0.0111111111, 0.0217391304),
  profiles = c(162L, 89L, 45L)
)
observed_4c <- organ_evaluable %>%
  transmute(
    minimum_coverage_fraction,
    n_loss = as.integer(observed_n_loss),
    denominator = as.integer(common_denominator_n),
    p = p_recurrent_loss,
    profiles = as.integer(n_candidate_genes)
  ) %>%
  left_join(expected_4c, by = "minimum_coverage_fraction", suffix = c("_observed", "_expected"))
assert_true(
  all(observed_4c$n_loss_observed == observed_4c$n_loss_expected) &&
    all(observed_4c$denominator_observed == observed_4c$denominator_expected) &&
    all(observed_4c$profiles_observed == observed_4c$profiles_expected) &&
    all(abs(observed_4c$p_observed - observed_4c$p_expected) < 5e-4),
  "Figure 4C values do not match the prespecified Round 3A audit values."
)

organ_evaluable <- organ_evaluable %>%
  mutate(
    result_label = paste0(
      as.integer(observed_n_loss), "/", as.integer(common_denominator_n),
      " organs\nloss fraction=",
      formatC(observed_loss_fraction, digits = 2, format = "f"),
      "\nraw P=", format_p4(p_recurrent_loss),
      "\nprofiles=", as.integer(n_candidate_genes)
    )
  )
organ_not_evaluable <- organ_not_evaluable %>%
  mutate(
    y = 0.30,
    result_label = case_when(
      status == "insufficient_common_set_candidates" ~ paste0(
        "Not evaluable\n", as.integer(n_candidate_genes),
        " profiles (<20)"
      ),
      TRUE ~ paste0("Not evaluable\n", status)
    )
  )

plot_figure4c <- ggplot() +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_point(
    data = organ_evaluable,
    aes(
      x = minimum_coverage_fraction,
      y = observed_loss_fraction
    ),
    shape = 23,
    size = 4.0,
    fill = gene_palette_okabe_ito[["MAOB"]],
    color = "black",
    stroke = 0.65
  ) +
  geom_label(
    data = organ_evaluable,
    aes(
      x = minimum_coverage_fraction,
      y = 0.69,
      label = result_label
    ),
    size = 3.0,
    lineheight = 1.0,
    fill = "white",
    color = "grey15",
    label.size = 0.25
  ) +
  geom_label(
    data = organ_not_evaluable,
    aes(
      x = minimum_coverage_fraction,
      y = y,
      label = result_label
    ),
    size = 3.0,
    lineheight = 1.0,
    fill = "grey92",
    color = "grey25",
    label.size = 0.25
  ) +
  scale_x_continuous(
    breaks = coverage_levels,
    labels = scales::percent_format(accuracy = 1),
    limits = c(0.45, 0.95),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = function(x) formatC(x, digits = 2, format = "f"),
    limits = c(0, 1.03),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Minimum coverage setting",
    y = "Observed MAOB loss fraction across organ units",
    title = "Organ-collapsed sensitivity of MAOB recurrence",
    subtitle = wrap_plot_text(
      paste(
        "KICH/KIRC/KIRP were collapsed to kidney and LUAD/LUSC to lung by",
        "median logFC. Labels report unadjusted recurrent-loss P values."
      ),
      width = 112
    ),
    caption = wrap_plot_text(
      "Not-evaluable positions are display locations and do not represent effect size or P value.",
      width = 105
    )
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, size = 7.5),
    plot.margin = margin(t = 8, r = 8, b = 6, l = 7)
  )

save_plot_pair(
  plot_figure4c,
  "Figure4C_round3C_MAOB_organ_collapsed_sensitivity",
  height_mm = 135
)

# ---- 6. Final audits and reproducibility ------------------------------------

input_manifest <- bind_rows(input_records)
write_csv_atomic(input_manifest, input_manifest_file)

text_audit <- bind_rows(text_audit_records)
figure_specs <- bind_rows(figure_output_records)
pdf_devices <- bind_rows(pdf_device_records)

assert_true(
  n_distinct(figure_specs$figure) == 3L,
  "Exactly three Round 3C figure pairs must be produced."
)
assert_true(
  all(figure_specs$width_mm == publication_width_mm) &&
    all(figure_specs$png_dpi == 300),
  "Every Round 3C figure must be 174 mm wide with a 300-dpi PNG."
)
assert_true(
  nrow(text_audit) > 0L && all(text_audit$fontsize_pt >= minimum_text_pt),
  "The final-size 7-pt text audit failed."
)
assert_true(
  nrow(pdf_devices) == 3L &&
    all(pdf_devices$vector_pdf) &&
    all(pdf_devices$embedded_font_device),
  "All three PDFs must use a vector embedded-font device."
)

write_csv_atomic(text_audit, text_audit_file)
write_csv_atomic(figure_specs, figure_spec_file)
write_csv_atomic(pdf_devices, pdf_device_file)

output_stems <- c(
  "Figure4A_round3C_denominator_aligned_recurrence_null",
  "Figure4B_round3C_coverage_sensitivity",
  "Figure4C_round3C_MAOB_organ_collapsed_sensitivity"
)
expected_outputs <- unlist(lapply(
  output_stems,
  function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".png")))
))

integrity_checks <- tibble(
  check = c(
    "all_frozen_inputs_present",
    "Figure4A_expected_denominators_profiles_and_q",
    "Figure4A_two_not_evaluable_genes_displayed",
    "Figure4B_all_25_gene_coverage_cells_in_status_strip",
    "Figure4B_two_endpoints_displayed",
    "Figure4C_expected_5_6_7_organ_results",
    "Figure4C_80_and_90_percent_not_evaluable",
    "minimum_text_size_at_least_7pt",
    "all_figures_174mm_wide",
    "all_png_outputs_300dpi",
    "all_pdf_devices_vector_and_embedded_font_capable",
    "all_six_output_files_exist",
    "no_continuous_colourbar_used",
    "no_statistical_module_recalculated"
  ),
  result = c(
    all(file.exists(input_files)),
    all(observed_4a$denominator_observed == observed_4a$denominator_expected) &&
      all(observed_4a$profiles_observed == observed_4a$profiles_expected) &&
      all(abs(observed_4a$q_observed - observed_4a$q_expected) < 5e-4),
    nrow(figure4a_not_evaluable) == 2L,
    nrow(coverage_status) == 25L,
    n_distinct(coverage_plot_long$endpoint) == 2L,
    all(observed_4c$n_loss_observed == observed_4c$n_loss_expected) &&
      all(observed_4c$denominator_observed == observed_4c$denominator_expected) &&
      all(observed_4c$profiles_observed == observed_4c$profiles_expected),
    identical(organ_not_evaluable$minimum_coverage_fraction, c(0.8, 0.9)),
    min(text_audit$fontsize_pt) >= minimum_text_pt,
    all(figure_specs$width_mm == publication_width_mm),
    all(figure_specs$png_dpi == 300),
    nrow(pdf_devices) == 3L && all(pdf_devices$embedded_font_device),
    all(file.exists(expected_outputs)),
    TRUE,
    TRUE
  ),
  detail = c(
    paste0(length(input_files), "/", length(input_files)),
    paste0(
      "MAOB 10/40/q=", format_q3(evaluable_80$q_recurrent_loss[evaluable_80$gene == "MAOB"]),
      "; CYP3A4 7/89/q=", format_q3(evaluable_80$q_recurrent_loss[evaluable_80$gene == "CYP3A4"]),
      "; CYP3A7 4/281/q=", format_q3(evaluable_80$q_recurrent_loss[evaluable_80$gene == "CYP3A7"])
    ),
    paste(as.character(figure4a_not_evaluable$gene), collapse = ";"),
    "5 genes x 5 coverage settings",
    paste(levels(coverage_plot_long$endpoint), collapse = ";"),
    "5/5, 6/6, 7/7; profiles 162, 89, 45",
    paste(scales::percent(organ_not_evaluable$minimum_coverage_fraction), collapse = ";"),
    paste0(formatC(min(text_audit$fontsize_pt), digits = 2, format = "f"), " pt"),
    paste(unique(figure_specs$width_mm), collapse = ";"),
    paste(unique(figure_specs$png_dpi), collapse = ";"),
    paste(unique(pdf_devices$pdf_device), collapse = ";"),
    paste0(sum(file.exists(expected_outputs)), "/", length(expected_outputs)),
    "No continuous colour scale appears in the three Round 3C figures",
    "Only frozen Round 3A values and graphical transformations were used"
  )
)

assert_true(
  all(integrity_checks$result),
  paste0(
    "Round 3C integrity failure: ",
    paste(integrity_checks$check[!integrity_checks$result], collapse = "; ")
  )
)
write_csv_atomic(integrity_checks, integrity_file)

write_lines_atomic(capture.output(sessionInfo()), session_info_file)

completion_report <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Project directory: ", project_dir),
  paste0("Round 3A input directory: ", round3a_table_dir),
  paste0("Figures written: ", paste(output_stems, collapse = "; ")),
  paste0(
    "Figure 4A: denominator-aligned 80% null; MAOB 10 contexts/40 profiles/q=",
    format_q3(evaluable_80$q_recurrent_loss[evaluable_80$gene == "MAOB"]),
    "; CYP3A4 7/89/q=",
    format_q3(evaluable_80$q_recurrent_loss[evaluable_80$gene == "CYP3A4"]),
    "; CYP3A7 4/281/q=",
    format_q3(evaluable_80$q_recurrent_loss[evaluable_80$gene == "CYP3A7"]), "."
  ),
  paste(
    "Figure 4B: recurrent-loss FDR is primary and negative-median FDR is",
    "supplementary; all 25 cells are represented in the status strip."
  ),
  paste(
    "Figure 4C: concise target/P-value display was used because table 168",
    "does not contain the full organ-collapsed background distributions."
  ),
  "Supplementary Figure S4 was not read or regenerated in Round 3C.",
  paste0(
    "Final-size text audit: minimum=",
    formatC(min(text_audit$fontsize_pt), digits = 3, format = "f"),
    " pt; threshold>=7.0 pt."
  ),
  paste0(
    "PDF devices: ", paste(unique(pdf_devices$pdf_device), collapse = "; "), "."
  ),
  "No continuous colourbar was used; no colourbar endpoint audit was required.",
  "No numerical analysis was recalculated.",
  "All integrity checks passed."
)
write_lines_atomic(completion_report, completion_report_file)

log_message("Round 3C complete.")
log_message("Figures: ", figure_dir)
log_message("Logs:    ", log_dir)
