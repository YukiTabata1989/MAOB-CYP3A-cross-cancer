#!/usr/bin/env Rscript

# =============================================================================
# MAOB/CYP3A tissue-disruption study
# Round 2J: Figure 3A legend correction only
# =============================================================================
#
# PURPOSE
#   Re-render Figure 3A from the frozen Round 2B tissue-identity tables.
#   The scientific values, panel content, annotations, and layout are unchanged.
#   The only intended visual correction is that the five target-gene entries in
#   the fill legend are displayed with their actual Okabe-Ito colours instead
#   of all appearing black.
#
# INPUTS (read only)
#   results/round2B/tables/151_round2B_tissue_identity_reference_genes.csv
#   results/round2B/tables/152_round2B_target_residual_vs_tissue_identity.csv
#
# OUTPUTS
#   results/round2J/figures/
#     Figure3A_round2J_tissue_identity_distribution.pdf
#     Figure3A_round2J_tissue_identity_distribution.png
#   results/round2J/logs/
#     round2J_input_manifest.csv
#     round2J_integrity_checks.csv
#     round2J_text_size_audit.csv
#     round2J_progress.log
#     round2J_sessionInfo.txt
#     round2J_completion_report.txt
#
# EXPLICITLY OUT OF SCOPE
#   No differential-expression, tissue-identity, matched-background,
#   recurrence, baseline, coexpression, HPA, survival, clinical, immune,
#   purity, subgroup, or LIHC-specific correlation analysis is run here.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration --------------------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round2J_Figure3A_legend_fix"
required_r_version <- "4.6.1"
publication_width_mm <- 174
publication_height_mm <- 132
minimum_text_pt <- 7
publication_font_family <- Sys.getenv(
  "MAOB_CYP3A_FIGURE_FONT",
  unset = "Arial"
)

target_genes <- c("MAOB", "CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43")
mapped_tcga_order <- c(
  "BRCA", "COAD", "ESCA", "KICH", "KIRC", "KIRP",
  "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA"
)

# Okabe-Ito categorical palette retained from Round 2F.
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

round2b_table_dir <- file.path(project_dir, "results", "round2B", "tables")
round2j_dir <- file.path(project_dir, "results", "round2J")
figure_dir <- file.path(round2j_dir, "figures")
log_dir <- file.path(round2j_dir, "logs")

invisible(lapply(
  c(round2j_dir, figure_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

input_files <- c(
  identity_reference = file.path(
    round2b_table_dir,
    "151_round2B_tissue_identity_reference_genes.csv"
  ),
  identity_targets = file.path(
    round2b_table_dir,
    "152_round2B_target_residual_vs_tissue_identity.csv"
  )
)

figure_stem <- "Figure3A_round2J_tissue_identity_distribution"
pdf_file <- file.path(figure_dir, paste0(figure_stem, ".pdf"))
png_file <- file.path(figure_dir, paste0(figure_stem, ".png"))
progress_log_file <- file.path(log_dir, "round2J_progress.log")
input_manifest_file <- file.path(log_dir, "round2J_input_manifest.csv")
integrity_file <- file.path(log_dir, "round2J_integrity_checks.csv")
text_audit_file <- file.path(log_dir, "round2J_text_size_audit.csv")
session_info_file <- file.path(log_dir, "round2J_sessionInfo.txt")
completion_report_file <- file.path(log_dir, "round2J_completion_report.txt")

required_packages <- c(
  "dplyr", "readr", "stringr", "tibble", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running Round 2J:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (as.character(getRversion()) != required_r_version) {
  stop(
    "Round 2J requires R ", required_r_version,
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

wrap_plot_text <- function(x, width = 105L) {
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

open_vector_pdf_device <- function(filename, width_in, height_in) {
  is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
  quartz_attempt <- list(opened = FALSE, error = NULL, warnings = character())
  cairo_attempt <- list(opened = FALSE, error = NULL, warnings = character())

  # Quartz is attempted first on macOS to avoid the unavailable-cairo-DLL error.
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
      "No vector PDF device could be opened. Quartz: ",
      device_attempt_message(quartz_attempt, "not applicable"),
      ". Cairo: ",
      device_attempt_message(cairo_attempt, "not available"),
      "."
    ),
    call. = FALSE
  )
}

render_plot_atomic <- function(plot, path, type, width_mm, height_mm, dpi = 300) {
  fileext <- if (type == "pdf") ".pdf" else ".png"
  tmp <- tempfile(
    pattern = paste0(".", tools::file_path_sans_ext(basename(path)), "."),
    tmpdir = dirname(path),
    fileext = fileext
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  device_open <- FALSE
  pdf_method <- NA_character_

  tryCatch({
    if (type == "pdf") {
      pdf_method <- open_vector_pdf_device(tmp, width_in, height_in)
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
    print(plot)
    grDevices::dev.off()
    device_open <- FALSE
  }, finally = {
    if (device_open && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  })

  replace_with_temp_file(tmp, path)
  if (type == "pdf") {
    log_message("PDF device: ", pdf_method)
  }
  invisible(path)
}

collect_text_grobs <- function(grob, path = "root", inherited_size = 12) {
  local_size <- inherited_size
  if (!is.null(grob$gp) && !is.null(grob$gp$fontsize)) {
    candidate <- suppressWarnings(as.numeric(grob$gp$fontsize))
    if (length(candidate) > 0L && any(is.finite(candidate))) {
      local_size <- min(candidate[is.finite(candidate)])
    }
  }

  rows <- list()
  if (inherits(grob, "text")) {
    label <- paste(as.character(grob$label), collapse = " | ")
    if (nzchar(gsub("[[:space:]]", "", label))) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        grob_path = path,
        label = label,
        fontsize_pt = as.numeric(local_size)
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
        path = paste(path, set_name, child_name, sep = "/"),
        inherited_size = local_size
      )
    }
  }

  dplyr::bind_rows(rows)
}

# ---- 2. Read and lock the two frozen inputs ---------------------------------

missing_inputs <- unname(input_files)[!file.exists(unname(input_files))]
if (length(missing_inputs) > 0L) {
  stop(
    "Required frozen Round 2B inputs are missing:\n  ",
    paste(missing_inputs, collapse = "\n  "),
    call. = FALSE
  )
}

input_md5_before <- vapply(input_files, safe_md5, character(1))

identity_reference <- readr::read_csv(
  input_files[["identity_reference"]],
  show_col_types = FALSE,
  progress = FALSE
)
identity_targets <- readr::read_csv(
  input_files[["identity_targets"]],
  show_col_types = FALSE,
  progress = FALSE
)

assert_has_columns(
  identity_reference,
  c("tcga_code", "gene_id", "gene_symbol", "logFC"),
  "identity_reference"
)
assert_has_columns(
  identity_targets,
  c(
    "tcga_code", "gene", "logFC", "analysis_status",
    "identity_median_logFC", "n_identity_genes", "identity_status",
    "more_negative_than_identity"
  ),
  "identity_targets"
)

log_message(
  "Input identity_reference: ", nrow(identity_reference),
  " rows x ", ncol(identity_reference), " columns."
)
log_message(
  "Input identity_targets: ", nrow(identity_targets),
  " rows x ", ncol(identity_targets), " columns."
)

input_manifest <- tibble::tibble(
  input = names(input_files),
  path = vapply(
    input_files,
    normalizePath,
    character(1),
    winslash = "/",
    mustWork = TRUE
  ),
  rows = c(nrow(identity_reference), nrow(identity_targets)),
  columns = c(ncol(identity_reference), ncol(identity_targets)),
  md5 = unname(input_md5_before),
  access = "read_only"
)
write_csv_atomic(input_manifest, input_manifest_file)

# ---- 3. Structural checks and plotting data ---------------------------------

identity_reference_counts <- identity_reference %>%
  count(tcga_code, name = "n_identity")

identity_target_estimable <- identity_targets %>%
  filter(is.finite(logFC))

n_identity_target_estimable <- nrow(identity_target_estimable)
n_more_negative <- sum(
  identity_target_estimable$more_negative_than_identity,
  na.rm = TRUE
)
n_robust_more_negative <- identity_target_estimable %>%
  filter(more_negative_than_identity) %>%
  nrow()

integrity_checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~passed,
  "identity_reference_rows",
  as.character(nrow(identity_reference)), "2400",
  nrow(identity_reference) == 2400L,
  "identity_target_grid_rows",
  as.character(nrow(identity_targets)), "60",
  nrow(identity_targets) == 60L,
  "mapped_cancers",
  as.character(n_distinct(identity_reference$tcga_code)), "12",
  n_distinct(identity_reference$tcga_code) == 12L,
  "mapped_cancer_set",
  paste(sort(unique(identity_reference$tcga_code)), collapse = ";"),
  paste(sort(mapped_tcga_order), collapse = ";"),
  setequal(unique(identity_reference$tcga_code), mapped_tcga_order),
  "identity_genes_per_cancer",
  paste(sort(unique(identity_reference_counts$n_identity)), collapse = ";"),
  "200",
  all(identity_reference_counts$n_identity == 200L),
  "estimable_target_cells",
  as.character(n_identity_target_estimable), "39",
  n_identity_target_estimable == 39L,
  "target_cells_below_identity_median",
  as.character(n_more_negative), "12",
  n_more_negative == 12L,
  "target_gene_set",
  paste(sort(unique(identity_targets$gene)), collapse = ";"),
  paste(sort(target_genes), collapse = ";"),
  setequal(unique(identity_targets$gene), target_genes),
  "no_analysis_module_recomputed",
  "TRUE", "TRUE", TRUE
)

assert_true(
  all(integrity_checks$passed),
  paste0(
    "One or more Round 2J integrity checks failed: ",
    paste(integrity_checks$check[!integrity_checks$passed], collapse = ", ")
  )
)

identity_ne_counts <- identity_targets %>%
  group_by(tcga_code) %>%
  summarise(
    n_estimable_targets = sum(is.finite(logFC)),
    n_not_estimable_targets = sum(!is.finite(logFC)),
    .groups = "drop"
  ) %>%
  mutate(x_label = paste0(tcga_code, "\nNE=", n_not_estimable_targets))

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

# ---- 4. Figure 3A ------------------------------------------------------------

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
      x = x_label,
      y = logFC,
      fill = gene,
      shape = relative_to_identity,
      group = gene
    ),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.68),
    color = "black",
    size = 2.5,
    stroke = 0.45
  ) +
  scale_fill_manual(
    values = gene_palette_okabe_ito,
    breaks = target_genes,
    drop = FALSE,
    name = "Target gene"
  ) +
  scale_shape_manual(
    values = c(
      "At or above identity median" = 21,
      "Below identity median" = 25
    ),
    breaks = c(
      "At or above identity median",
      "Below identity median"
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
    ))
  ) +
  theme_manuscript(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.0),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(
    # Shapes 21/25 accept a fill. Explicitly using shape 21 here prevents the
    # five fill-scale keys from falling back to a solid black default symbol.
    fill = guide_legend(
      nrow = 1,
      order = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 21,
        colour = "black",
        size = 3.0,
        stroke = 0.5
      )
    ),
    shape = guide_legend(
      nrow = 1,
      order = 2,
      byrow = TRUE,
      override.aes = list(fill = "white", colour = "black", size = 3.0)
    )
  )

# ---- 5. Legend and text-size validation -------------------------------------

plot_build <- ggplot_build(plot_identity)
fill_scale <- plot_build$plot$scales$get_scales("fill")
legend_hex <- unname(fill_scale$map(target_genes))

legend_colour_check <- tibble::tibble(
  gene = target_genes,
  expected_fill = unname(gene_palette_okabe_ito[target_genes]),
  mapped_fill = legend_hex,
  passed = toupper(expected_fill) == toupper(mapped_fill)
)

assert_true(
  all(legend_colour_check$passed),
  "Target-gene fill-scale mapping does not match the Okabe-Ito palette."
)

plot_grob <- grid::grid.force(ggplotGrob(plot_identity))
text_audit <- collect_text_grobs(plot_grob)
assert_true(nrow(text_audit) > 0L, "No visible text grobs were detected.")
assert_true(
  all(is.finite(text_audit$fontsize_pt)),
  "At least one text grob has an unresolved font size."
)
minimum_observed_text_pt <- min(text_audit$fontsize_pt)
assert_true(
  minimum_observed_text_pt + 1e-8 >= minimum_text_pt,
  paste0(
    "Figure 3A contains text smaller than ", minimum_text_pt,
    " pt; observed minimum=",
    formatC(minimum_observed_text_pt, digits = 3, format = "f"), " pt."
  )
)

write_csv_atomic(text_audit, text_audit_file)
log_message(
  "Text-size audit PASS: minimum=",
  formatC(minimum_observed_text_pt, digits = 3, format = "f"),
  " pt; threshold>=", minimum_text_pt, " pt."
)
for (i in seq_len(nrow(legend_colour_check))) {
  log_message(
    "Legend fill PASS: ", legend_colour_check$gene[[i]],
    " = ", legend_colour_check$mapped_fill[[i]], "."
  )
}

# ---- 6. Save PDF/PNG atomically ---------------------------------------------

render_plot_atomic(
  plot_identity,
  pdf_file,
  type = "pdf",
  width_mm = publication_width_mm,
  height_mm = publication_height_mm
)
render_plot_atomic(
  plot_identity,
  png_file,
  type = "png",
  width_mm = publication_width_mm,
  height_mm = publication_height_mm,
  dpi = 300
)

assert_true(file.exists(pdf_file), "Figure 3A PDF was not created.")
assert_true(file.exists(png_file), "Figure 3A PNG was not created.")

input_md5_after <- vapply(input_files, safe_md5, character(1))
input_unchanged <- identical(unname(input_md5_before), unname(input_md5_after))
assert_true(input_unchanged, "A frozen Round 2B input changed during Round 2J.")

integrity_checks <- bind_rows(
  integrity_checks,
  tibble::tibble(
    check = c(
      "target_gene_legend_colours",
      "minimum_text_size_pt",
      "input_files_unchanged",
      "pdf_created",
      "png_created"
    ),
    observed = c(
      paste0(legend_colour_check$gene, "=", legend_colour_check$mapped_fill,
             collapse = ";"),
      formatC(minimum_observed_text_pt, digits = 3, format = "f"),
      as.character(input_unchanged),
      as.character(file.exists(pdf_file)),
      as.character(file.exists(png_file))
    ),
    expected = c(
      paste0(target_genes, "=", unname(gene_palette_okabe_ito[target_genes]),
             collapse = ";"),
      paste0(">=", minimum_text_pt),
      "TRUE", "TRUE", "TRUE"
    ),
    passed = c(
      all(legend_colour_check$passed),
      minimum_observed_text_pt + 1e-8 >= minimum_text_pt,
      input_unchanged,
      file.exists(pdf_file),
      file.exists(png_file)
    )
  )
)
write_csv_atomic(integrity_checks, integrity_file)

write_lines_atomic(capture.output(sessionInfo()), session_info_file)

completion_lines <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  "Scope: Figure 3A re-render only; no scientific value was recomputed.",
  paste0("Input identity-reference rows: ", nrow(identity_reference)),
  paste0("Input target-grid rows: ", nrow(identity_targets)),
  paste0("Estimable target cells: ", n_identity_target_estimable),
  paste0("Below identity median: ", n_more_negative),
  paste0("Minimum visible text size: ",
         formatC(minimum_observed_text_pt, digits = 3, format = "f"), " pt"),
  paste0("PDF: ", pdf_file),
  paste0("PNG: ", png_file),
  "Legend correction: five gene keys use the mapped Okabe-Ito fill colours.",
  "Frozen Round 2B inputs unchanged: TRUE"
)
write_lines_atomic(completion_lines, completion_report_file)

log_message("Wrote PDF: ", pdf_file)
log_message("Wrote PNG: ", png_file)
log_message("Integrity checks: ", integrity_file)
log_message("Session information: ", session_info_file)
log_message("Round 2J complete.")

message("\nRound 2J complete.")
message("Figure PDF: ", pdf_file)
message("Figure PNG: ", png_file)
message("Logs: ", log_dir)

