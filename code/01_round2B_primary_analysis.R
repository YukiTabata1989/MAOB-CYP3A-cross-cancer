#!/usr/bin/env Rscript

# =============================================================================
# Tissue-specific disruption of MAOB and CYP3A family expression
# Round 2B: LIHC-complete paired analysis and covariance-preserving recurrence
# =============================================================================
#
# PURPOSE
#   1) Confirm TCGA tumor-versus-adjacent-normal shifts from GDC STAR raw counts
#      using a paired limma-voom model (patient + condition).
#   2) Include LIHC MAOB/CYP3A4 in the paired expression analysis while keeping
#      every LIHC coexpression/correlation analysis outside Paper 1.
#   3) Ask whether target-gene shifts exceed the loss of expression expected for
#      genes matched on normal abundance, tissue specificity, expression breadth,
#      organ specificity, and TCGA adjacent-normal abundance.
#   4) Test cross-cancer recurrence with the same background gene followed across
#      cancers, preserving gene-level cross-cancer covariance.
#   5) Re-estimate the normal-baseline/disruption relationship using count-level
#      TCGA effects and official GTEx v8 tissue-median TPM.
#   6) Estimate MAOB-CYP3A coexpression rewiring in non-LIHC cancers only.
#   7) Validate normal-tissue rankings using HPA-only RNA data, not the HPA
#      consensus dataset (which contains GTEx information).
#
# INTERPRETATION BOUNDARY
#   This is an observational transcriptomic analysis. It cannot establish
#   whether MAOB loss precedes cancer, causes cancer, changes GGA abundance, or
#   changes enzyme activity. The matched-background analysis tests selectivity,
#   not temporal or mechanistic causality.
#
# TWO-PAPER BOUNDARY
#   Paper 1 owns paired tumor-versus-adjacent-non-tumor expression disruption,
#   including LIHC MAOB and CYP3A4 as two cells in the pan-cancer design.
#   Paper 2 owns every LIHC-specific MAOB-CYP3A4 relationship/heterogeneity and
#   every clinical endpoint. Accordingly, LIHC is excluded as an entire cancer
#   from the coexpression module before expression vectors are accessed.
#   No subgroup, survival, Cox, interaction, clinical-adjustment, HNF4A,
#   HCCDB18, ICGC-LIRI-JP, pathway, immune, or purity analysis is present.
#
# RUNNING
#   Set MAOB_CYP3A_PROJECT_DIR to the same project directory used for Rounds 1/2,
#   then run this script. Existing raw-input and non-LIHC Round 2 fit caches are
#   reused. LIHC is refit once into a distinct Round 2B cache. Outputs are written
#   under results/round2B and do not overwrite Round 1 or Round 2 results.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(23016802)

# ---- 0. Configuration --------------------------------------------------------

analysis_id <- "MAOB_CYP3A_tissue_disruption_round2B"

primary_genes   <- c("MAOB", "CYP3A4", "CYP3A5", "CYP3A7")
secondary_genes <- c("CYP3A43")
target_genes    <- c(primary_genes, secondary_genes)
cyp_partners    <- setdiff(target_genes, "MAOB")

lihc_expression_genes_returned_to_paper1 <- c("MAOB", "CYP3A4")
lihc_excluded_from_coexpression <- TRUE

effect_logfc_threshold       <- log2(1.5)
minimum_pairs_primary        <- 10L
minimum_pairs_coexpression   <- 30L
minimum_pairs_exploratory    <- 20L
minimum_detection_fraction   <- 0.20
coexpression_effect_threshold <- 0.20
coexpression_bootstrap_reps  <- 2000L
baseline_permutation_reps    <- 5000L
same_gene_primary_min_coverage <- 0.80
same_gene_sensitivity_coverages <- c(0.50, 0.60, 0.70, 0.80, 0.90)
same_gene_min_candidates      <- 20L
matched_background_max       <- 1000L
matched_background_min       <- 200L
matched_background_caliper   <- 2.0
tissue_identity_top_n        <- 200L
tissue_identity_min_n        <- 50L

project_dir <- normalizePath(
  Sys.getenv("MAOB_CYP3A_PROJECT_DIR", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)

raw_dir       <- file.path(project_dir, "data", "raw")
processed_dir <- file.path(project_dir, "data", "processed")
source_table_dir <- file.path(project_dir, "results", "tables")
round2b_result_dir <- file.path(project_dir, "results", "round2B")
table_dir     <- file.path(round2b_result_dir, "tables")
figure_dir    <- file.path(round2b_result_dir, "figures")
log_dir       <- file.path(round2b_result_dir, "logs")
gdc_dir       <- file.path(raw_dir, "GDCdata_round2")
cache_dir     <- file.path(processed_dir, "round2_project_cache")

invisible(lapply(
  c(raw_dir, processed_dir, table_dir, figure_dir, log_dir, gdc_dir, cache_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "ggplot2", "scales", "sandwich", "lmtest", "tibble", "statmod"
)
bioc_packages <- c("TCGAbiolinks", "SummarizedExperiment", "edgeR", "limma")
required_packages <- c(cran_packages, bioc_packages)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  missing_cran <- intersect(missing_packages, cran_packages)
  missing_bioc <- intersect(missing_packages, bioc_packages)
  install_text <- c(
    "Install the missing packages before running Round 2B:",
    paste(" ", paste(missing_packages, collapse = ", ")),
    if (length(missing_cran) > 0L) {
      paste0(
        "install.packages(c(",
        paste(sprintf("\"%s\"", missing_cran), collapse = ", "),
        "))"
      )
    },
    if (length(missing_bioc) > 0L) {
      paste0(
        "BiocManager::install(c(",
        paste(sprintf("\"%s\"", missing_bioc), collapse = ", "),
        "))"
      )
    }
  )
  stop(paste(install_text, collapse = "\n"), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

message("Project directory: ", project_dir)

# ---- 1. Utilities, scope contract, and inherited mapping --------------------

timestamp_utc <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

message_log_file <- file.path(log_dir, "round2B_progress.log")
log_message <- function(...) {
  x <- paste0(..., collapse = "")
  message(x)
  cat(paste0("[", timestamp_utc(), "] ", x, "\n"),
      file = message_log_file, append = TRUE)
  invisible(x)
}

save_rds_atomic <- function(object, path, compress = TRUE) {
  tmp <- paste0(path, ".tmp")
  saveRDS(object, tmp, compress = compress)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    stop("Could not atomically save cache: ", path, call. = FALSE)
  }
  invisible(path)
}

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  readr::write_csv(x, tmp, na = "")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    stop("Could not atomically save table: ", path, call. = FALSE)
  }
  invisible(path)
}

save_plot_pair <- function(plot, stem, width, height, dpi = 320) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  # Use the standard PDF device. Some macOS R installations report Cairo as
  # available but cannot load its DLL at device-open time.
  ggplot2::ggsave(
    pdf_path, plot = plot, width = width, height = height,
    device = "pdf", useDingbats = FALSE
  )
  ggplot2::ggsave(
    png_path, plot = plot, width = width, height = height, dpi = dpi
  )
  invisible(c(pdf_path, png_path))
}

strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

normalise_names <- function(x) {
  tolower(gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", x)))
}

first_existing_column <- function(x, candidates) {
  hit <- intersect(candidates, names(x))
  if (length(hit) == 0L) NA_character_ else hit[[1]]
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

scope_contract <- tibble::tribble(
  ~analysis_item, ~status, ~implementation,
  "GDC STAR-count paired differential expression", "allowed",
    "patient + condition limma-voom; LIHC MAOB/CYP3A4 included in Paper 1",
  "Matched-background selective-disruption test", "allowed",
    "matched on GTEx baseline, tau, breadth, organ specificity, and TCGA normal abundance",
  "Same-gene cross-cancer recurrence null", "allowed",
    "one background gene profile retained across cancers; coverage sensitivity reported",
  "Tissue-identity reference loss", "allowed",
    "GTEx-defined reference genes; targets excluded from reference set",
  "Normal baseline versus count-level disruption", "secondary",
    "organ-level exploratory summary including paired LIHC target effects",
  "Within-patient coexpression rewiring", "allowed_non_LIHC_only",
    "same paired TCGA patients; entire LIHC project excluded before vector access",
  "GTEx-versus-TCGA direct testing", "not_run",
    "source and disease status are confounded",
  "HPA normal validation", "allowed",
    "HPA-only RNA tissue dataset; consensus dataset is not used",
  "Causal, temporal, GGA, or enzyme-activity inference", "forbidden",
    "not identifiable from bulk cross-sectional expression data",
  "LIHC endpoint or detailed clinical analysis", "forbidden",
    "reserved for the separate LIHC manuscript",
  "LIHC paired MAOB/CYP3A4 expression shifts", "allowed_in_pancancer_context",
    "Paper 1 endpoint; no LIHC-specific standalone clinical interpretation",
  "Any LIHC MAOB-CYP3A coexpression/correlation", "forbidden",
    "entire LIHC project excluded from the coexpression module",
  "HNF4A, external liver cohorts, pathway, immune, or purity analysis", "forbidden",
    "outside this manuscript and/or reserved"
)

write_csv_atomic(scope_contract, file.path(table_dir, "100_round2B_scope_contract.csv"))

config_table <- tibble::tibble(
  parameter = c(
    "analysis_id", "effect_logfc_threshold", "minimum_pairs_primary",
    "minimum_pairs_coexpression", "minimum_pairs_exploratory",
    "minimum_detection_fraction", "coexpression_effect_threshold",
    "coexpression_bootstrap_reps", "baseline_permutation_reps",
    "same_gene_primary_min_coverage",
    "same_gene_sensitivity_coverages", "same_gene_min_candidates",
    "matched_background_max",
    "matched_background_min", "matched_background_caliper",
    "tissue_identity_top_n", "tissue_identity_min_n",
    "lihc_expression_genes_returned_to_paper1",
    "lihc_excluded_from_coexpression"
  ),
  value = as.character(c(
    analysis_id, effect_logfc_threshold, minimum_pairs_primary,
    minimum_pairs_coexpression, minimum_pairs_exploratory,
    minimum_detection_fraction, coexpression_effect_threshold,
    coexpression_bootstrap_reps, baseline_permutation_reps,
    same_gene_primary_min_coverage,
    paste(same_gene_sensitivity_coverages, collapse = ";"),
    same_gene_min_candidates, matched_background_max,
    matched_background_min, matched_background_caliper,
    tissue_identity_top_n, tissue_identity_min_n,
    paste(lihc_expression_genes_returned_to_paper1, collapse = ";"),
    lihc_excluded_from_coexpression
  ))
)
write_csv_atomic(config_table, file.path(table_dir, "101_round2B_configuration.csv"))

mapping_file <- file.path(source_table_dir, "01_cancer_tissue_mapping.csv")
if (!file.exists(mapping_file)) {
  stop(
    "Round 1 mapping table is missing: ", mapping_file,
    "\nRun Round 1 in this project directory before Round 2B.",
    call. = FALSE
  )
}

cancer_key <- readr::read_csv(mapping_file, show_col_types = FALSE)
required_mapping_columns <- c(
  "tcga_code", "tier", "organ", "gtex_pattern",
  "mapping_valid_for_baseline_model"
)
if (!all(required_mapping_columns %in% names(cancer_key))) {
  stop("Round 1 cancer mapping table has an unexpected structure.", call. = FALSE)
}

primary_cancers <- cancer_key %>%
  filter(tier == "Tier1") %>%
  arrange(tcga_code)

requested_codes_text <- Sys.getenv("MAOB_CYP3A_ROUND2_CANCERS", unset = "")
if (nzchar(requested_codes_text)) {
  requested_codes <- trimws(strsplit(requested_codes_text, ",", fixed = TRUE)[[1]])
  invalid_codes <- setdiff(requested_codes, primary_cancers$tcga_code)
  if (length(invalid_codes) > 0L) {
    stop("Unknown/non-Tier1 code(s): ", paste(invalid_codes, collapse = ", "))
  }
  primary_cancers <- primary_cancers %>% filter(tcga_code %in% requested_codes)
  log_message(
    "Subset run requested via MAOB_CYP3A_ROUND2_CANCERS: ",
    paste(primary_cancers$tcga_code, collapse = ", ")
  )
}

if (nrow(primary_cancers) == 0L) stop("No Tier1 cancers selected.", call. = FALSE)

is_excluded_coexpression_cancer <- function(tcga_code) {
  isTRUE(lihc_excluded_from_coexpression) && tcga_code == "LIHC"
}

# ---- 2. GDC query, paired-aliquot selection, and count extraction ------------

normalise_sample_type <- function(x) {
  y <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    y %in% c("primary tumor", "primary tumour", "tp") ~ "Primary Tumor",
    y %in% c("solid tissue normal", "normal", "nt") ~ "Solid Tissue Normal",
    TRUE ~ NA_character_
  )
}

infer_sample_type_from_barcode <- function(barcode) {
  code <- substr(as.character(barcode), 14L, 15L)
  dplyr::case_when(
    code == "01" ~ "Primary Tumor",
    code == "11" ~ "Solid Tissue Normal",
    TRUE ~ NA_character_
  )
}

find_tcga_barcode_vector <- function(col_data, fallback) {
  is_barcode <- function(v) {
    mean(grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}", as.character(v)))
  }

  if (is_barcode(fallback) > 0.8) return(as.character(fallback))

  scores <- vapply(col_data, function(v) {
    tryCatch(is_barcode(v), error = function(e) 0)
  }, numeric(1))
  if (length(scores) > 0L && max(scores, na.rm = TRUE) > 0.8) {
    return(as.character(col_data[[which.max(scores)]]))
  }
  stop("Could not identify TCGA sample barcodes in the prepared object.", call. = FALSE)
}

manifest_sample_metadata <- function(manifest) {
  x <- as.data.frame(manifest)
  original_names <- names(x)
  names(x) <- normalise_names(names(x))

  barcode_col <- first_existing_column(
    x,
    c("cases", "sample_submitter_id", "aliquot_submitter_id", "submitter_id")
  )
  type_col <- first_existing_column(x, c("sample_type", "sample_types"))

  if (is.na(barcode_col)) {
    stop(
      "The GDC manifest lacks a recognizable sample/case barcode column. Columns: ",
      paste(original_names, collapse = ", "),
      call. = FALSE
    )
  }

  barcode <- as.character(x[[barcode_col]])
  barcode <- sub("[,;|].*$", "", barcode)
  sample_type <- if (!is.na(type_col)) {
    normalise_sample_type(x[[type_col]])
  } else {
    infer_sample_type_from_barcode(barcode)
  }
  sample_type[is.na(sample_type)] <- infer_sample_type_from_barcode(barcode[is.na(sample_type)])

  tibble::tibble(
    manifest_row = seq_len(nrow(x)),
    manifest_barcode = barcode,
    case_id = substr(barcode, 1L, 12L),
    sample_type = sample_type
  )
}

choose_unstranded_assay <- function(se) {
  assay_names <- SummarizedExperiment::assayNames(se)
  if ("unstranded" %in% assay_names) return("unstranded")

  candidate <- assay_names[
    grepl("unstrand", assay_names, ignore.case = TRUE) &
      !grepl("tpm|fpkm", assay_names, ignore.case = TRUE)
  ]
  if (length(candidate) == 1L) return(candidate[[1]])

  stop(
    "Unstranded raw-count assay not found. Available assays: ",
    paste(assay_names, collapse = ", "),
    call. = FALSE
  )
}

extract_gene_annotation <- function(se, counts) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  names(rd) <- normalise_names(names(rd))

  gene_id_col <- first_existing_column(rd, c("gene_id", "ensembl_gene_id"))
  symbol_col  <- first_existing_column(rd, c("gene_name", "gene_symbol", "external_gene_name"))
  type_col    <- first_existing_column(rd, c("gene_type", "gene_biotype", "type"))

  gene_id <- if (!is.na(gene_id_col)) rd[[gene_id_col]] else rownames(counts)
  gene_symbol <- if (!is.na(symbol_col)) rd[[symbol_col]] else NA_character_
  gene_type <- if (!is.na(type_col)) rd[[type_col]] else NA_character_

  tibble::tibble(
    source_row = seq_len(nrow(counts)),
    gene_id = strip_ensembl_version(gene_id),
    gene_symbol = as.character(gene_symbol),
    gene_type = as.character(gene_type),
    total_count = rowSums(counts, na.rm = TRUE)
  )
}

collapse_duplicate_gene_ids <- function(counts, annotation) {
  valid <- !is.na(annotation$gene_id) & nzchar(annotation$gene_id)
  counts <- counts[valid, , drop = FALSE]
  annotation <- annotation[valid, , drop = FALSE]

  if (!anyDuplicated(annotation$gene_id)) {
    rownames(counts) <- annotation$gene_id
    return(list(counts = counts, annotation = annotation %>% select(-source_row)))
  }

  collapsed_counts <- rowsum(counts, group = annotation$gene_id, reorder = FALSE)
  representative <- annotation %>%
    arrange(gene_id, desc(total_count), source_row) %>%
    group_by(gene_id) %>%
    slice(1L) %>%
    ungroup() %>%
    select(-source_row, -total_count)
  representative <- representative[match(rownames(collapsed_counts), representative$gene_id), ]

  list(counts = collapsed_counts, annotation = representative)
}

prepare_paired_project_input <- function(tcga_code) {
  project <- paste0("TCGA-", tcga_code)
  project_cache_dir <- file.path(cache_dir, project)
  dir.create(project_cache_dir, recursive = TRUE, showWarnings = FALSE)
  input_cache <- file.path(project_cache_dir, "paired_raw_input.rds")

  if (file.exists(input_cache)) {
    log_message(project, ": loading cached paired raw-count input")
    return(readRDS(input_cache))
  }

  log_message(project, ": querying GDC STAR - Counts metadata")
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )

  manifest <- TCGAbiolinks::getResults(query)
  if (is.null(manifest) || nrow(manifest) == 0L) {
    stop(project, ": no STAR-count files returned by GDCquery.", call. = FALSE)
  }
  manifest_meta <- manifest_sample_metadata(manifest)

  paired_cases <- manifest_meta %>%
    filter(!is.na(case_id), !is.na(sample_type)) %>%
    group_by(case_id) %>%
    summarise(
      has_tumor = any(sample_type == "Primary Tumor"),
      has_normal = any(sample_type == "Solid Tissue Normal"),
      .groups = "drop"
    ) %>%
    filter(has_tumor, has_normal) %>%
    pull(case_id)

  if (length(paired_cases) < 2L) {
    stop(project, ": fewer than two paired cases in the GDC manifest.", call. = FALSE)
  }

  keep_manifest <- manifest_meta$case_id %in% paired_cases
  query_paired <- query
  query_paired$results[[1]] <- as.data.frame(manifest)[keep_manifest, , drop = FALSE]

  manifest_export <- as.data.frame(manifest)[keep_manifest, , drop = FALSE]
  manifest_export[] <- lapply(manifest_export, function(x) {
    if (is.list(x)) vapply(x, paste, character(1), collapse = ";") else x
  })
  write_csv_atomic(
    bind_cols(manifest_export,
              manifest_meta[keep_manifest, c("case_id", "sample_type")]),
    file.path(table_dir, paste0("110_", tcga_code, "_GDC_paired_manifest.csv"))
  )
  save_rds_atomic(query_paired, file.path(project_cache_dir, "paired_query.rds"))

  log_message(
    project, ": downloading ", sum(keep_manifest),
    " files for ", length(paired_cases), " potential paired cases"
  )
  TCGAbiolinks::GDCdownload(
    query = query_paired,
    method = "api",
    files.per.chunk = 20L,
    directory = gdc_dir
  )

  log_message(project, ": preparing SummarizedExperiment")
  se <- TCGAbiolinks::GDCprepare(
    query = query_paired,
    directory = gdc_dir,
    summarizedExperiment = TRUE
  )

  assay_name <- choose_unstranded_assay(se)
  counts <- as.matrix(SummarizedExperiment::assay(se, assay_name))
  storage.mode(counts) <- "integer"

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  names(cd) <- normalise_names(names(cd))
  sample_barcode <- find_tcga_barcode_vector(cd, colnames(counts))
  type_col <- first_existing_column(cd, c("sample_type", "definition", "short_letter_code"))
  sample_type <- if (!is.na(type_col)) {
    normalise_sample_type(cd[[type_col]])
  } else {
    rep(NA_character_, ncol(counts))
  }
  inferred_type <- infer_sample_type_from_barcode(sample_barcode)
  sample_type[is.na(sample_type)] <- inferred_type[is.na(sample_type)]

  sample_meta_all <- tibble::tibble(
    original_index = seq_len(ncol(counts)),
    sample_barcode = sample_barcode,
    case_id = substr(sample_barcode, 1L, 12L),
    sample_type = sample_type,
    library_size = colSums(counts, na.rm = TRUE)
  ) %>%
    filter(
      !is.na(case_id),
      sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
    )

  # If multiple aliquots exist for a case/type, retain the largest library.
  # This rule is deterministic; all candidates are recorded in the QC table.
  aliquot_selection <- sample_meta_all %>%
    group_by(case_id, sample_type) %>%
    arrange(desc(library_size), sample_barcode, .by_group = TRUE) %>%
    mutate(selected_aliquot = row_number() == 1L) %>%
    ungroup()

  selected <- aliquot_selection %>%
    filter(selected_aliquot) %>%
    group_by(case_id) %>%
    filter(
      any(sample_type == "Primary Tumor"),
      any(sample_type == "Solid Tissue Normal")
    ) %>%
    ungroup() %>%
    mutate(
      condition = if_else(sample_type == "Primary Tumor", "Tumor", "Normal"),
      condition = factor(condition, levels = c("Normal", "Tumor"))
    ) %>%
    arrange(case_id, condition)

  pair_counts <- selected %>% count(case_id)
  if (nrow(pair_counts) < 2L || any(pair_counts$n != 2L)) {
    stop(project, ": paired sample selection failed after preparation.", call. = FALSE)
  }

  aliquot_selection <- aliquot_selection %>%
    mutate(
      retained_paired_case = case_id %in% selected$case_id,
      selection_reason = case_when(
        selected_aliquot & retained_paired_case ~ "retained_largest_library",
        !selected_aliquot ~ "duplicate_aliquot_smaller_library",
        TRUE ~ "case_not_complete_after_preparation"
      )
    )
  write_csv_atomic(
    aliquot_selection,
    file.path(table_dir, paste0("111_", tcga_code, "_aliquot_selection.csv"))
  )

  selected_counts <- counts[, selected$original_index, drop = FALSE]
  colnames(selected_counts) <- selected$sample_barcode
  annotation <- extract_gene_annotation(se, counts)
  collapsed <- collapse_duplicate_gene_ids(selected_counts, annotation)

  selected <- selected %>%
    select(case_id, sample_barcode, sample_type, condition, library_size)

  output <- list(
    tcga_code = tcga_code,
    project = project,
    assay_name = assay_name,
    counts = collapsed$counts,
    annotation = collapsed$annotation,
    sample_meta = selected,
    n_pairs = n_distinct(selected$case_id),
    created_utc = timestamp_utc()
  )
  save_rds_atomic(output, input_cache)
  rm(se, counts, selected_counts)
  invisible(gc())
  output
}

# ---- 3. Per-project paired limma-voom analysis ------------------------------

pick_target_rows <- function(annotation, counts, genes = target_genes) {
  row_mean <- rowMeans(counts, na.rm = TRUE)
  purrr::map_dfr(genes, function(g) {
    idx <- which(annotation$gene_symbol == g)
    if (length(idx) == 0L) {
      return(tibble::tibble(gene = g, row_index = NA_integer_, gene_id = NA_character_))
    }
    best <- idx[which.max(row_mean[idx])]
    tibble::tibble(
      gene = g,
      row_index = as.integer(best),
      gene_id = annotation$gene_id[[best]]
    )
  })
}

make_target_expression_matrices <- function(
    counts, annotation, target_index, effective_library_size) {
  n_samples <- ncol(counts)
  raw_target <- matrix(
    NA_real_,
    nrow = length(target_genes),
    ncol = n_samples,
    dimnames = list(target_genes, colnames(counts))
  )

  present <- target_index %>% filter(!is.na(row_index))
  if (nrow(present) > 0L) {
    for (i in seq_len(nrow(present))) {
      raw_target[present$gene[[i]], ] <- counts[present$row_index[[i]], ]
    }
  }

  target_logcpm <- matrix(
    NA_real_,
    nrow = nrow(raw_target),
    ncol = ncol(raw_target),
    dimnames = dimnames(raw_target)
  )
  present_rows <- which(rowSums(is.finite(raw_target)) == n_samples)
  if (length(present_rows) > 0L) {
    target_logcpm[present_rows, ] <- edgeR::cpm(
      raw_target[present_rows, , drop = FALSE],
      lib.size = effective_library_size,
      log = TRUE,
      prior.count = 0.5
    )
  }

  list(raw_counts = raw_target, logcpm = target_logcpm)
}

run_paired_limma_project <- function(project_input) {
  tcga_code <- project_input$tcga_code
  project <- project_input$project
  project_cache_dir <- file.path(cache_dir, project)
  # Non-LIHC projects are analytically unchanged and can reuse the Round 2 fit.
  # LIHC must be refit because Round 2 removed MAOB/CYP3A4 before filtering.
  analysis_cache <- if (tcga_code == "LIHC") {
    file.path(project_cache_dir, "paired_limma_voom_results_round2B_unmasked.rds")
  } else {
    file.path(project_cache_dir, "paired_limma_voom_results.rds")
  }

  if (file.exists(analysis_cache)) {
    log_message(project, ": loading cached paired limma-voom results")
    return(readRDS(analysis_cache))
  }

  counts <- project_input$counts
  annotation <- project_input$annotation
  sample_meta <- project_input$sample_meta
  n_pairs <- project_input$n_pairs

  if (!identical(colnames(counts), sample_meta$sample_barcode)) {
    stop(project, ": count columns and sample metadata are not aligned.", call. = FALSE)
  }

  target_index <- pick_target_rows(annotation, counts)

  # Paper 1 owns paired expression disruption for all five genes in all cancers,
  # including LIHC MAOB/CYP3A4. No target row is removed before the model.
  model_keep_rows <- rep(TRUE, nrow(counts))
  model_counts <- counts[model_keep_rows, , drop = FALSE]
  model_annotation <- annotation[model_keep_rows, , drop = FALSE]

  condition <- factor(as.character(sample_meta$condition), levels = c("Normal", "Tumor"))
  patient <- factor(sample_meta$case_id)
  design <- stats::model.matrix(~ patient + condition)
  if (!"conditionTumor" %in% colnames(design)) {
    stop(project, ": paired design lacks conditionTumor coefficient.", call. = FALSE)
  }
  if (qr(design)$rank != ncol(design)) {
    stop(project, ": paired design matrix is not full rank.", call. = FALSE)
  }

  dge0 <- edgeR::DGEList(counts = model_counts)
  keep_expression <- edgeR::filterByExpr(
    dge0,
    group = condition,
    min.count = 10L,
    min.total.count = 15L
  )
  if (sum(keep_expression) < 1000L) {
    stop(project, ": unexpectedly few genes passed filterByExpr.", call. = FALSE)
  }

  dge <- dge0[keep_expression, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  effective_library_size <- dge$samples$lib.size * dge$samples$norm.factors

  log_message(
    project, ": fitting paired limma-voom model (", n_pairs,
    " pairs; ", sum(keep_expression), " retained genes)"
  )
  voom_object <- limma::voom(dge, design = design, plot = FALSE)
  fit <- limma::lmFit(voom_object, design = design)
  fit <- limma::eBayes(fit, robust = TRUE)
  tt <- limma::topTable(
    fit,
    coef = "conditionTumor",
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )

  retained_annotation <- model_annotation[keep_expression, , drop = FALSE]
  retained_annotation <- retained_annotation[
    match(rownames(tt), retained_annotation$gene_id),
    , drop = FALSE
  ]
  normal_idx <- condition == "Normal"
  tumor_idx <- condition == "Tumor"

  all_gene_de <- tibble::tibble(
    tcga_code = tcga_code,
    gene_id = rownames(tt),
    gene_symbol = retained_annotation$gene_symbol,
    gene_type = retained_annotation$gene_type,
    n_pairs = n_pairs,
    logFC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    p_value = tt$P.Value,
    q_value_project = tt$adj.P.Val,
    B = tt$B,
    mean_normal_logcpm = rowMeans(voom_object$E[, normal_idx, drop = FALSE]),
    mean_tumor_logcpm = rowMeans(voom_object$E[, tumor_idx, drop = FALSE])
  )

  target_expression <- make_target_expression_matrices(
    counts = counts,
    annotation = annotation,
    target_index = target_index,
    effective_library_size = effective_library_size
  )

  target_de <- purrr::map_dfr(target_genes, function(g) {
    index_row <- target_index %>% filter(gene == g)
    de_rows <- all_gene_de %>% filter(gene_symbol == g)
    if (nrow(de_rows) > 1L) de_rows <- de_rows %>% slice_max(AveExpr, n = 1L)

    if (is.na(index_row$row_index[[1]])) {
      return(tibble::tibble(
        tcga_code = tcga_code, gene = g, gene_id = NA_character_, n_pairs = n_pairs,
        logFC = NA_real_, AveExpr = NA_real_, t = NA_real_, p_value = NA_real_,
        q_value_project = NA_real_, B = NA_real_, mean_normal_logcpm = NA_real_,
        mean_tumor_logcpm = NA_real_, analysis_status = "gene_not_found"
      ))
    }

    if (nrow(de_rows) == 0L) {
      return(tibble::tibble(
        tcga_code = tcga_code, gene = g, gene_id = index_row$gene_id[[1]],
        n_pairs = n_pairs, logFC = NA_real_, AveExpr = NA_real_, t = NA_real_,
        p_value = NA_real_, q_value_project = NA_real_, B = NA_real_,
        mean_normal_logcpm = NA_real_, mean_tumor_logcpm = NA_real_,
        analysis_status = "low_expression_not_estimable"
      ))
    }

    de_rows %>%
      transmute(
        tcga_code, gene = g, gene_id, n_pairs, logFC, AveExpr, t, p_value,
        q_value_project, B, mean_normal_logcpm, mean_tumor_logcpm,
        analysis_status = if_else(
          n_pairs >= minimum_pairs_primary,
          "primary_count_level",
          "insufficient_pairs_lt10"
        )
      )
  })

  output <- list(
    tcga_code = tcga_code,
    project = project,
    n_pairs = n_pairs,
    all_gene_de = all_gene_de,
    target_de = target_de,
    target_counts = target_expression$raw_counts,
    target_logcpm = target_expression$logcpm,
    sample_meta = sample_meta,
    filter_summary = tibble::tibble(
      tcga_code = tcga_code,
      n_pairs = n_pairs,
      n_genes_input = nrow(model_counts),
      n_genes_retained = sum(keep_expression),
      retained_fraction = mean(keep_expression),
      design_rank = qr(design)$rank,
      design_columns = ncol(design)
    ),
    created_utc = timestamp_utc()
  )
  save_rds_atomic(output, analysis_cache)
  rm(voom_object, fit, tt, dge, dge0)
  invisible(gc())
  output
}

# Process projects serially to keep peak memory predictable. Each project is
# cached independently, so rerunning resumes at the first incomplete project.
project_analyses <- vector("list", nrow(primary_cancers))
names(project_analyses) <- primary_cancers$tcga_code

for (i in seq_len(nrow(primary_cancers))) {
  code <- primary_cancers$tcga_code[[i]]
  project_input <- prepare_paired_project_input(code)
  project_analyses[[code]] <- run_paired_limma_project(project_input)
  rm(project_input)
  invisible(gc())
}

pair_qc <- purrr::map_dfr(project_analyses, "filter_summary") %>%
  left_join(
    primary_cancers %>% select(tcga_code, organ, tier),
    by = "tcga_code"
  ) %>%
  arrange(tcga_code)
write_csv_atomic(pair_qc, file.path(table_dir, "120_round2B_paired_project_qc.csv"))

target_results <- purrr::map_dfr(project_analyses, "target_de") %>%
  left_join(
    primary_cancers %>%
      select(tcga_code, organ, mapping_valid_for_baseline_model, mapping_note),
    by = "tcga_code"
  ) %>%
  mutate(
    reserved = FALSE,
    analysis_status = case_when(
      n_pairs < minimum_pairs_primary ~ "insufficient_pairs_lt10",
      TRUE ~ analysis_status
    )
  )

eligible_target_p <- target_results$analysis_status == "primary_count_level" &
  is.finite(target_results$p_value)
target_results$q_value_global <- NA_real_
target_results$q_value_global[eligible_target_p] <- p.adjust(
  target_results$p_value[eligible_target_p],
  method = "BH"
)
target_results <- target_results %>%
  mutate(
    passes_statistical_screen = !is.na(q_value_global) & q_value_global < 0.05,
    passes_effect_screen = !is.na(logFC) & abs(logFC) >= effect_logfc_threshold,
    robust_count_level_change = passes_statistical_screen & passes_effect_screen,
    count_level_interpretation = case_when(
      analysis_status != "primary_count_level" ~ "not_estimable_primary",
      robust_count_level_change & logFC < 0 ~ "loss",
      robust_count_level_change & logFC > 0 ~ "gain",
      TRUE ~ "not_robust"
    )
  )

write_csv_atomic(
  target_results,
  file.path(table_dir, "121_round2B_target_paired_limma_voom_results.csv")
)

target_status_matrix <- target_results %>%
  select(
    tcga_code, organ, gene, n_pairs, analysis_status, logFC,
    q_value_global, robust_count_level_change, count_level_interpretation
  ) %>%
  arrange(tcga_code, match(gene, target_genes))
write_csv_atomic(
  target_status_matrix,
  file.path(table_dir, "122_round2B_target_status_matrix.csv")
)

# Optional transformed-expression concordance audit against Round 1. This is a
# robustness audit only; Round 2B count-level estimates remain authoritative.
stage1_paired_file <- file.path(source_table_dir, "21_TCGA_paired_results.csv")
if (file.exists(stage1_paired_file)) {
  stage1_paired <- readr::read_csv(stage1_paired_file, show_col_types = FALSE)
  required_stage1_columns <- c(
    "tcga_code", "gene", "median_paired_shift", "robust_targeted_change",
    "interpretation"
  )
  if (all(required_stage1_columns %in% names(stage1_paired))) {
    stage1_round2b_concordance <- target_results %>%
      select(
        tcga_code, gene, round2B_logFC = logFC,
        round2B_q = q_value_global,
        round2B_robust = robust_count_level_change,
        round2B_interpretation = count_level_interpretation
      ) %>%
      inner_join(
        stage1_paired %>%
          select(
            tcga_code, gene,
            stage1_shift = median_paired_shift,
            stage1_robust = robust_targeted_change,
            stage1_interpretation = interpretation
          ),
        by = c("tcga_code", "gene")
      ) %>%
      mutate(
        finite_overlap = is.finite(stage1_shift) & is.finite(round2B_logFC),
        direction_concordant = finite_overlap &
          sign(stage1_shift) == sign(round2B_logFC)
      )

    finite_concordance <- stage1_round2b_concordance %>%
      filter(finite_overlap)
    concordance_summary <- tibble::tibble(
      n_overlap = nrow(finite_concordance),
      direction_concordance_fraction = mean(
        finite_concordance$direction_concordant, na.rm = TRUE
      ),
      pearson_r = if (nrow(finite_concordance) >= 3L) {
        stats::cor(
          finite_concordance$stage1_shift,
          finite_concordance$round2B_logFC,
          method = "pearson"
        )
      } else {
        NA_real_
      },
      spearman_rho = if (nrow(finite_concordance) >= 3L) {
        stats::cor(
          finite_concordance$stage1_shift,
          finite_concordance$round2B_logFC,
          method = "spearman"
        )
      } else {
        NA_real_
      },
      interpretation = paste(
        "Robustness audit across transformed-expression and count-level",
        "paired analyses; not an independent validation dataset"
      )
    )
    write_csv_atomic(
      stage1_round2b_concordance,
      file.path(table_dir, "123_round1_round2B_concordance_cells.csv")
    )
    write_csv_atomic(
      concordance_summary,
      file.path(table_dir, "124_round1_round2B_concordance_summary.csv")
    )
  }
}

# A compact figure confirms effect direction and shows all non-estimable cells.
code_labels <- pair_qc %>%
  transmute(tcga_code, x_label = paste0(tcga_code, "\nn=", n_pairs))

paired_plot_data <- tidyr::expand_grid(
  tcga_code = primary_cancers$tcga_code,
  gene = target_genes
) %>%
  left_join(target_results, by = c("tcga_code", "gene")) %>%
  left_join(code_labels, by = "tcga_code") %>%
  mutate(
    label = case_when(
      analysis_status == "low_expression_not_estimable" ~ "L",
      robust_count_level_change ~ "*",
      TRUE ~ ""
    ),
    gene = factor(gene, levels = rev(target_genes)),
    x_label = factor(x_label, levels = code_labels$x_label)
  )

plot_paired_count <- ggplot(
  paired_plot_data,
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
    title = "Count-level paired TCGA expression shifts",
    subtitle = paste0(
      "* global target FDR<0.05 and |logFC|>=log2(1.5); ",
      "LIHC paired MAOB/CYP3A4 included; L=filtered for low expression"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )

save_plot_pair(
  plot_paired_count,
  "Figure2_round2B_count_level_paired_heatmap",
  width = 11, height = 5
)

# ---- 4. Official GTEx v8 tissue-median reference ----------------------------

gtex_median_filename <- paste0(
  "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_",
  "gene_median_tpm.gct.gz"
)
gtex_median_urls <- c(
  paste0(
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/",
    gtex_median_filename
  ),
  paste0(
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/",
    gtex_median_filename
  )
)
gtex_median_file <- file.path(
  raw_dir,
  "GTEx_v8_gene_median_tpm.gct.gz"
)

if (!file.exists(gtex_median_file)) {
  log_message("Downloading official GTEx v8 gene-median TPM by tissue")
  downloaded <- FALSE
  last_download_error <- NULL
  for (candidate_url in gtex_median_urls) {
    downloaded <- tryCatch({
      utils::download.file(
        candidate_url, gtex_median_file, mode = "wb", quiet = FALSE
      )
      file.exists(gtex_median_file) && file.info(gtex_median_file)$size > 1000000
    }, error = function(e) {
      last_download_error <<- conditionMessage(e)
      FALSE
    })
    if (isTRUE(downloaded)) break
    if (file.exists(gtex_median_file)) unlink(gtex_median_file)
  }
  if (!isTRUE(downloaded)) {
    stop(
      "Could not download the GTEx v8 tissue-median file. Last error: ",
      last_download_error,
      call. = FALSE
    )
  }
}

gtex_cache_file <- file.path(processed_dir, "round2_GTEx_v8_gene_metrics.rds")
if (file.exists(gtex_cache_file)) {
  gtex_reference <- readRDS(gtex_cache_file)
} else {
  log_message("Reading and summarizing GTEx v8 tissue-median matrix")
  gtex_raw <- data.table::fread(
    gtex_median_file,
    skip = 2L,
    data.table = FALSE,
    check.names = FALSE
  )
  if (ncol(gtex_raw) < 10L) {
    stop("GTEx median GCT has an unexpected structure.", call. = FALSE)
  }

  id_col <- names(gtex_raw)[[1]]
  description_col <- names(gtex_raw)[[2]]
  tissue_columns <- setdiff(names(gtex_raw), c(id_col, description_col))
  gtex_matrix <- as.matrix(gtex_raw[, tissue_columns, drop = FALSE])
  storage.mode(gtex_matrix) <- "numeric"
  rownames(gtex_matrix) <- strip_ensembl_version(gtex_raw[[id_col]])

  max_tpm <- apply(gtex_matrix, 1L, max, na.rm = TRUE)
  cross_tissue_median_tpm <- apply(gtex_matrix, 1L, median, na.rm = TRUE)
  tau <- rep(NA_real_, nrow(gtex_matrix))
  positive_max <- is.finite(max_tpm) & max_tpm > 0
  tau[positive_max] <- (
    ncol(gtex_matrix) -
      rowSums(gtex_matrix[positive_max, , drop = FALSE] /
                max_tpm[positive_max], na.rm = TRUE)
  ) / (ncol(gtex_matrix) - 1L)

  gene_metrics <- tibble::tibble(
    gene_id = rownames(gtex_matrix),
    gene_symbol = as.character(gtex_raw[[description_col]]),
    max_tpm = max_tpm,
    cross_tissue_median_tpm = cross_tissue_median_tpm,
    tau = tau,
    breadth_0_1tpm = rowMeans(gtex_matrix >= 0.1, na.rm = TRUE),
    breadth_1tpm = rowMeans(gtex_matrix >= 1, na.rm = TRUE)
  )

  # Defensive de-duplication after removing Ensembl version suffixes.
  if (anyDuplicated(gene_metrics$gene_id)) {
    retain <- gene_metrics %>%
      mutate(row_index = row_number()) %>%
      arrange(gene_id, desc(max_tpm), row_index) %>%
      group_by(gene_id) %>%
      slice(1L) %>%
      ungroup()
    gtex_matrix <- gtex_matrix[retain$row_index, , drop = FALSE]
    rownames(gtex_matrix) <- retain$gene_id
    gene_metrics <- retain %>% select(-row_index)
  }

  gtex_reference <- list(
    gene_metrics = gene_metrics,
    tpm_matrix = gtex_matrix,
    tissue_columns = tissue_columns,
    source_urls = gtex_median_urls,
    created_utc = timestamp_utc()
  )
  save_rds_atomic(gtex_reference, gtex_cache_file)
  rm(gtex_raw, gtex_matrix, gene_metrics)
  invisible(gc())
}

gtex_metrics <- gtex_reference$gene_metrics
gtex_matrix <- gtex_reference$tpm_matrix
gtex_tissues <- gtex_reference$tissue_columns

build_gtex_cancer_baseline <- function(key_row) {
  pattern <- key_row$gtex_pattern[[1]]
  if (is.na(pattern) || !nzchar(pattern)) return(NULL)
  matched_tissues <- grep(pattern, gtex_tissues, value = TRUE)
  if (length(matched_tissues) == 0L) {
    warning(
      "No GTEx tissue matched pattern '", pattern,
      "' for ", key_row$tcga_code[[1]]
    )
    return(NULL)
  }

  baseline_tpm <- if (length(matched_tissues) == 1L) {
    gtex_matrix[, matched_tissues]
  } else {
    apply(gtex_matrix[, matched_tissues, drop = FALSE], 1L, median, na.rm = TRUE)
  }

  gtex_metrics %>%
    transmute(
      tcga_code = key_row$tcga_code[[1]],
      organ = key_row$organ[[1]],
      mapping_valid_for_baseline_model = as.logical(
        key_row$mapping_valid_for_baseline_model[[1]]
      ),
      gtex_pattern = pattern,
      matched_gtex_tissues = paste(matched_tissues, collapse = ";"),
      gene_id, gene_symbol,
      baseline_tpm = as.numeric(baseline_tpm),
      baseline_log2_tpm1 = log2(as.numeric(baseline_tpm) + 1),
      organ_specificity_log2 = log2(
        (as.numeric(baseline_tpm) + 0.1) /
          (cross_tissue_median_tpm + 0.1)
      ),
      max_tpm, cross_tissue_median_tpm, tau,
      breadth_0_1tpm, breadth_1tpm
    )
}

gtex_baseline <- purrr::map_dfr(
  split(primary_cancers, seq_len(nrow(primary_cancers))),
  build_gtex_cancer_baseline
)
save_rds_atomic(
  gtex_baseline,
  file.path(processed_dir, "round2B_GTEx_v8_cancer_baselines.rds")
)

gtex_mapping_qc <- gtex_baseline %>%
  distinct(
    tcga_code, organ, mapping_valid_for_baseline_model,
    gtex_pattern, matched_gtex_tissues
  ) %>%
  right_join(
    primary_cancers %>%
      select(tcga_code, organ, mapping_valid_for_baseline_model,
             gtex_pattern, mapping_note),
    by = c(
      "tcga_code", "organ", "mapping_valid_for_baseline_model",
      "gtex_pattern"
    )
  ) %>%
  arrange(tcga_code)
write_csv_atomic(gtex_mapping_qc, file.path(table_dir, "140_GTEx_v8_mapping_qc.csv"))

gtex_target_atlas <- gtex_baseline %>%
  filter(gene_symbol %in% target_genes) %>%
  arrange(tcga_code, match(gene_symbol, target_genes))
write_csv_atomic(
  gtex_target_atlas,
  file.path(table_dir, "141_GTEx_v8_target_organ_baselines.csv")
)

gtex_target_metrics <- gtex_metrics %>%
  filter(gene_symbol %in% target_genes) %>%
  arrange(match(gene_symbol, target_genes))
write_csv_atomic(
  gtex_target_metrics,
  file.path(table_dir, "142_GTEx_v8_target_tissue_specificity.csv")
)

target_abundance_context <- target_results %>%
  left_join(
    gtex_baseline %>%
      filter(gene_symbol %in% target_genes) %>%
      group_by(tcga_code, gene_symbol) %>%
      slice_max(baseline_tpm, n = 1L, with_ties = FALSE) %>%
      ungroup() %>%
      select(
        tcga_code, gene = gene_symbol, baseline_tpm,
        baseline_log2_tpm1, organ_specificity_log2
      ),
    by = c("tcga_code", "gene")
  ) %>%
  mutate(
    gtex_normal_abundance_class = case_when(
      !is.finite(baseline_tpm) ~ "unavailable",
      baseline_tpm < 0.1 ~ "below_0.1_TPM",
      baseline_tpm < 1 ~ "0.1_to_below_1_TPM",
      TRUE ~ "at_least_1_TPM"
    ),
    low_normal_abundance_caution = is.finite(baseline_tpm) & baseline_tpm < 1,
    interpretation_note = if_else(
      low_normal_abundance_caution,
      paste(
        "Large fold changes from a low normal-tissue baseline require",
        "abundance-aware interpretation"
      ),
      NA_character_
    )
  )
write_csv_atomic(
  target_abundance_context,
  file.path(table_dir, "145_round2B_target_normal_abundance_context.csv")
)

# Optional concordance audit against the Round 1 target atlas. This audit does
# not enter any Round 2 hypothesis test.
stage1_atlas_file <- file.path(source_table_dir, "10_GTEx_normal_atlas.csv")
if (file.exists(stage1_atlas_file)) {
  stage1_atlas <- readr::read_csv(stage1_atlas_file, show_col_types = FALSE)
  if (all(c("primary_tissue", "gene", "median_tpm") %in% names(stage1_atlas))) {
    round2_tissue_target <- purrr::map_dfr(gtex_tissues, function(tissue) {
      idx <- which(gtex_metrics$gene_symbol %in% target_genes)
      tibble::tibble(
        primary_tissue = tissue,
        gene = gtex_metrics$gene_symbol[idx],
        round2_official_median_tpm = as.numeric(gtex_matrix[idx, tissue])
      )
    })
    stage1_concordance <- stage1_atlas %>%
      select(primary_tissue, gene, stage1_sample_median_tpm = median_tpm) %>%
      inner_join(round2_tissue_target, by = c("primary_tissue", "gene"))
    write_csv_atomic(
      stage1_concordance,
      file.path(table_dir, "143_GTEx_round1_round2B_baseline_audit.csv")
    )
  }
}

# ---- 5. Tissue-identity loss and matched-background empirical null ----------

all_gene_de <- purrr::map_dfr(project_analyses, "all_gene_de")

valid_background_cancers <- primary_cancers %>%
  filter(mapping_valid_for_baseline_model) %>%
  pull(tcga_code)

build_tissue_identity_summary <- function(tcga_code) {
  de <- all_gene_de %>% filter(.data$tcga_code == .env$tcga_code)
  baseline <- gtex_baseline %>%
    filter(
      .data$tcga_code == .env$tcga_code,
      mapping_valid_for_baseline_model
    )

  joined <- de %>%
    inner_join(
      baseline %>%
        select(
          gene_id, baseline_tpm, baseline_log2_tpm1,
          organ_specificity_log2, tau, breadth_1tpm
        ),
      by = "gene_id"
    ) %>%
    filter(
      grepl("^protein[_ -]?coding$", gene_type, ignore.case = TRUE),
      !(gene_symbol %in% target_genes),
      baseline_tpm >= 1,
      organ_specificity_log2 >= 1,
      is.finite(logFC)
    ) %>%
    arrange(desc(organ_specificity_log2), desc(baseline_tpm))

  selected <- joined %>% slice_head(n = tissue_identity_top_n)
  n_selected <- nrow(selected)
  summary <- tibble::tibble(
    tcga_code = tcga_code,
    n_identity_candidates = nrow(joined),
    n_identity_genes = n_selected,
    identity_median_logFC = if (n_selected >= tissue_identity_min_n) {
      median(selected$logFC, na.rm = TRUE)
    } else {
      NA_real_
    },
    identity_mean_logFC = if (n_selected >= tissue_identity_min_n) {
      mean(selected$logFC, na.rm = TRUE)
    } else {
      NA_real_
    },
    identity_fraction_loss_threshold = if (n_selected >= tissue_identity_min_n) {
      mean(selected$logFC <= -effect_logfc_threshold, na.rm = TRUE)
    } else {
      NA_real_
    },
    identity_status = if_else(
      n_selected >= tissue_identity_min_n,
      "estimable",
      "insufficient_identity_genes"
    )
  )

  list(
    summary = summary,
    genes = selected %>%
      transmute(
        tcga_code, gene_id, gene_symbol, baseline_tpm,
        organ_specificity_log2, logFC, p_value, q_value_project
      )
  )
}

identity_objects <- purrr::map(valid_background_cancers, build_tissue_identity_summary)
identity_summary <- purrr::map_dfr(identity_objects, "summary") %>%
  left_join(primary_cancers %>% select(tcga_code, organ), by = "tcga_code")
identity_genes <- purrr::map_dfr(identity_objects, "genes")

write_csv_atomic(
  identity_summary,
  file.path(table_dir, "150_round2B_tissue_identity_loss_summary.csv")
)
write_csv_atomic(
  identity_genes,
  file.path(table_dir, "151_round2B_tissue_identity_reference_genes.csv")
)

target_identity_residual <- target_results %>%
  filter(tcga_code %in% valid_background_cancers) %>%
  left_join(
    identity_summary %>%
      select(tcga_code, identity_median_logFC, n_identity_genes, identity_status),
    by = "tcga_code"
  ) %>%
  mutate(
    residual_vs_identity = logFC - identity_median_logFC,
    more_negative_than_identity = !is.na(residual_vs_identity) &
      residual_vs_identity < 0
  )
write_csv_atomic(
  target_identity_residual,
  file.path(table_dir, "152_round2B_target_residual_vs_tissue_identity.csv")
)

standardized_matching_distance <- function(pool, target_row, feature_names) {
  combined <- bind_rows(
    pool %>% select(all_of(feature_names)),
    target_row %>% select(all_of(feature_names))
  )
  usable <- feature_names[vapply(combined[feature_names], function(x) {
    all(is.finite(x)) && stats::sd(x) > 0
  }, logical(1))]
  if (length(usable) < 3L) return(rep(NA_real_, nrow(pool)))

  scaled <- scale(combined[, usable, drop = FALSE])
  target_scaled <- scaled[nrow(scaled), , drop = FALSE]
  sqrt(rowSums(
    sweep(
      scaled[seq_len(nrow(pool)), , drop = FALSE],
      2L,
      as.numeric(target_scaled),
      "-"
    )^2
  ))
}

match_background_for_target <- function(tcga_code, target_gene) {
  target_de_row <- target_results %>%
    filter(
      .data$tcga_code == .env$tcga_code,
      gene == .env$target_gene,
      analysis_status == "primary_count_level",
      !reserved,
      is.finite(logFC)
    )
  if (nrow(target_de_row) != 1L) {
    return(list(
      summary = tibble::tibble(
        tcga_code = tcga_code, gene = target_gene,
        matching_status = "target_not_estimable"
      ),
      genes = tibble::tibble()
    ))
  }

  baseline <- gtex_baseline %>%
    filter(
      .data$tcga_code == .env$tcga_code,
      mapping_valid_for_baseline_model
    )
  target_gtex <- baseline %>%
    filter(gene_symbol == .env$target_gene) %>%
    arrange(desc(baseline_tpm)) %>%
    slice(1L)
  if (nrow(target_gtex) != 1L) {
    return(list(
      summary = tibble::tibble(
        tcga_code = tcga_code, gene = target_gene,
        matching_status = "target_missing_from_GTEx"
      ),
      genes = tibble::tibble()
    ))
  }

  de <- all_gene_de %>% filter(.data$tcga_code == .env$tcga_code)
  target_tcga <- de %>%
    filter(gene_symbol == .env$target_gene) %>%
    arrange(desc(AveExpr)) %>%
    slice(1L)
  if (nrow(target_tcga) != 1L) {
    return(list(
      summary = tibble::tibble(
        tcga_code = tcga_code, gene = target_gene,
        matching_status = "target_missing_after_filter"
      ),
      genes = tibble::tibble()
    ))
  }

  pool <- de %>%
    inner_join(
      baseline %>%
        select(
          gene_id, baseline_tpm, baseline_log2_tpm1, tau,
          breadth_1tpm, organ_specificity_log2
        ),
      by = "gene_id"
    ) %>%
    filter(
      grepl("^protein[_ -]?coding$", gene_type, ignore.case = TRUE),
      !(gene_symbol %in% target_genes),
      is.finite(logFC), is.finite(p_value),
      is.finite(baseline_log2_tpm1), is.finite(tau),
      is.finite(breadth_1tpm), is.finite(organ_specificity_log2),
      is.finite(mean_normal_logcpm)
    ) %>%
    distinct(gene_id, .keep_all = TRUE)

  target_features <- tibble::tibble(
    baseline_log2_tpm1 = target_gtex$baseline_log2_tpm1,
    tau = target_gtex$tau,
    breadth_1tpm = target_gtex$breadth_1tpm,
    organ_specificity_log2 = target_gtex$organ_specificity_log2,
    mean_normal_logcpm = target_tcga$mean_normal_logcpm
  )
  if (nrow(pool) < 50L || any(!is.finite(unlist(target_features)))) {
    return(list(
      summary = tibble::tibble(
        tcga_code = tcga_code, gene = target_gene,
        matching_status = "insufficient_finite_matching_data"
      ),
      genes = tibble::tibble()
    ))
  }
  feature_names <- names(target_features)
  distance <- standardized_matching_distance(pool, target_features, feature_names)
  pool$matching_distance <- distance
  pool <- pool %>%
    filter(is.finite(matching_distance)) %>%
    arrange(matching_distance)

  within_caliper <- pool %>%
    filter(matching_distance <= matched_background_caliper) %>%
    slice_head(n = matched_background_max)
  if (nrow(within_caliper) < matched_background_min) {
    matched <- pool %>% slice_head(n = matched_background_min)
    caliper_status <- "caliper_relaxed_to_minimum"
  } else {
    matched <- within_caliper
    caliper_status <- "within_caliper"
  }

  if (nrow(matched) < 50L) {
    return(list(
      summary = tibble::tibble(
        tcga_code = tcga_code, gene = target_gene,
        matching_status = "insufficient_background_genes"
      ),
      genes = tibble::tibble()
    ))
  }

  target_logfc <- target_de_row$logFC[[1]]
  lower_p <- (1 + sum(matched$logFC <= target_logfc)) / (nrow(matched) + 1)
  upper_p <- (1 + sum(matched$logFC >= target_logfc)) / (nrow(matched) + 1)
  directional_p <- if (target_logfc < 0) lower_p else upper_p
  two_sided_p <- min(1, 2 * min(lower_p, upper_p))

  summary <- tibble::tibble(
    tcga_code = tcga_code,
    gene = target_gene,
    target_gene_id = target_de_row$gene_id[[1]],
    target_logFC = target_logfc,
    target_q_value_global = target_de_row$q_value_global[[1]],
    target_baseline_tpm = target_gtex$baseline_tpm[[1]],
    target_tau = target_gtex$tau[[1]],
    target_breadth_1tpm = target_gtex$breadth_1tpm[[1]],
    target_organ_specificity_log2 = target_gtex$organ_specificity_log2[[1]],
    target_mean_normal_logcpm = target_tcga$mean_normal_logcpm[[1]],
    n_background = nrow(matched),
    background_median_baseline_log2_tpm1 = median(matched$baseline_log2_tpm1),
    background_median_tau = median(matched$tau),
    background_median_breadth_1tpm = median(matched$breadth_1tpm),
    background_median_organ_specificity_log2 = median(
      matched$organ_specificity_log2
    ),
    background_median_normal_logcpm = median(matched$mean_normal_logcpm),
    background_median_logFC = median(matched$logFC),
    background_mean_logFC = mean(matched$logFC),
    excess_logFC_vs_background = target_logfc - median(matched$logFC),
    background_percentile = mean(matched$logFC <= target_logfc),
    empirical_p_lower = lower_p,
    empirical_p_upper = upper_p,
    empirical_p_directional = directional_p,
    empirical_p_two_sided = two_sided_p,
    median_matching_distance = median(matched$matching_distance),
    maximum_matching_distance = max(matched$matching_distance),
    matching_status = paste("estimable", caliper_status, sep = ":")
  )

  genes <- matched %>%
    transmute(
      tcga_code = .env$tcga_code,
      target_gene = .env$target_gene,
      background_rank = row_number(),
      background_gene_id = gene_id,
      background_gene_symbol = gene_symbol,
      background_logFC = logFC,
      background_p_value = p_value,
      baseline_tpm,
      baseline_log2_tpm1,
      tau,
      breadth_1tpm,
      organ_specificity_log2,
      mean_normal_logcpm,
      matching_distance
    )
  list(summary = summary, genes = genes)
}

matching_grid <- tidyr::expand_grid(
  tcga_code = valid_background_cancers,
  gene = target_genes
)

matching_objects <- purrr::pmap(
  matching_grid,
  function(tcga_code, gene) match_background_for_target(tcga_code, gene)
)
matched_background_summary <- purrr::map_dfr(matching_objects, "summary")
matched_background_genes <- purrr::map_dfr(matching_objects, "genes")

matched_background_summary$q_empirical_directional <- NA_real_
eligible_empirical <- grepl("^estimable", matched_background_summary$matching_status) &
  is.finite(matched_background_summary$empirical_p_directional)
matched_background_summary$q_empirical_directional[eligible_empirical] <- p.adjust(
  matched_background_summary$empirical_p_directional[eligible_empirical],
  method = "BH"
)
matched_background_summary <- matched_background_summary %>%
  mutate(
    selective_disruption_candidate =
      !is.na(q_empirical_directional) & q_empirical_directional < 0.05 &
      abs(target_logFC) >= effect_logfc_threshold,
    selective_direction = case_when(
      selective_disruption_candidate & target_logFC < 0 ~ "selective_loss_candidate",
      selective_disruption_candidate & target_logFC > 0 ~ "selective_gain_candidate",
      grepl("^estimable", matching_status) ~ "not_selective_after_matching",
      TRUE ~ "not_estimable"
    )
  ) %>%
  left_join(primary_cancers %>% select(tcga_code, organ), by = "tcga_code")

write_csv_atomic(
  matched_background_summary,
  file.path(table_dir, "160_round2B_matched_background_cell_results.csv")
)
write_csv_atomic(
  matched_background_genes,
  file.path(table_dir, "161_round2B_matched_background_gene_lists.csv")
)

matching_balance_summary <- matched_background_summary %>%
  filter(grepl("^estimable", matching_status)) %>%
  transmute(
    tcga_code, organ, gene, n_background, matching_status,
    median_matching_distance, maximum_matching_distance,
    target_baseline_log2_tpm1 = log2(target_baseline_tpm + 1),
    background_median_baseline_log2_tpm1,
    difference_baseline_log2_tpm1 = target_baseline_log2_tpm1 -
      background_median_baseline_log2_tpm1,
    target_tau, background_median_tau,
    difference_tau = target_tau - background_median_tau,
    target_breadth_1tpm, background_median_breadth_1tpm,
    difference_breadth_1tpm = target_breadth_1tpm -
      background_median_breadth_1tpm,
    target_organ_specificity_log2,
    background_median_organ_specificity_log2,
    difference_organ_specificity_log2 = target_organ_specificity_log2 -
      background_median_organ_specificity_log2,
    target_mean_normal_logcpm, background_median_normal_logcpm,
    difference_normal_logcpm = target_mean_normal_logcpm -
      background_median_normal_logcpm
  )
write_csv_atomic(
  matching_balance_summary,
  file.path(table_dir, "165_round2B_matching_balance_summary.csv")
)

# Covariance-preserving recurrence null. Each null unit is one background gene
# observed across multiple cancers; its cross-cancer profile is kept intact.
# Because no single matched gene need occur in every cancer, the primary test
# requires >=80% coverage and compares loss/gain fractions plus median logFC.
build_same_gene_null <- function(
    target_gene,
    minimum_coverage_fraction = same_gene_primary_min_coverage,
    return_profiles = FALSE) {
  observed_cells <- matched_background_summary %>%
    filter(
      gene == .env$target_gene,
      grepl("^estimable", matching_status),
      is.finite(target_logFC)
    ) %>%
    arrange(tcga_code)

  empty_summary <- function(status) {
    tibble::tibble(
      gene = target_gene,
      minimum_coverage_fraction = minimum_coverage_fraction,
      minimum_coverage_cancers = if (nrow(observed_cells) > 0L) {
        ceiling(minimum_coverage_fraction * nrow(observed_cells))
      } else {
        NA_integer_
      },
      n_cancers = nrow(observed_cells),
      observed_n_loss = NA_integer_, observed_loss_fraction = NA_real_,
      observed_n_gain = NA_integer_, observed_gain_fraction = NA_real_,
      observed_median_logFC = NA_real_, n_candidate_genes = NA_integer_,
      null_median_loss_fraction = NA_real_,
      null_median_gain_fraction = NA_real_,
      null_median_logFC = NA_real_,
      p_recurrent_loss = NA_real_, p_recurrent_gain = NA_real_,
      p_median_negative = NA_real_, status = status
    )
  }

  if (nrow(observed_cells) < 4L) {
    result <- list(
      summary = empty_summary("insufficient_cancers"),
      profiles = tibble::tibble()
    )
    return(if (return_profiles) result else result$summary)
  }

  minimum_coverage_cancers <- ceiling(
    minimum_coverage_fraction * nrow(observed_cells)
  )
  profiles <- matched_background_genes %>%
    filter(
      .data$target_gene == .env$target_gene,
      tcga_code %in% observed_cells$tcga_code,
      is.finite(background_logFC)
    ) %>%
    distinct(tcga_code, background_gene_id, .keep_all = TRUE) %>%
    group_by(background_gene_id, background_gene_symbol) %>%
    summarise(
      n_coverage = n_distinct(tcga_code),
      covered_cancers = paste(sort(unique(tcga_code)), collapse = ";"),
      n_loss = sum(background_logFC <= -effect_logfc_threshold),
      loss_fraction = n_loss / n_coverage,
      n_gain = sum(background_logFC >= effect_logfc_threshold),
      gain_fraction = n_gain / n_coverage,
      median_logFC = median(background_logFC),
      mean_logFC = mean(background_logFC),
      mean_matching_distance = mean(matching_distance),
      maximum_matching_distance = max(matching_distance),
      .groups = "drop"
    ) %>%
    mutate(
      target_gene = target_gene,
      n_target_cancers = nrow(observed_cells),
      coverage_fraction = n_coverage / n_target_cancers,
      eligible_primary_null = n_coverage >= minimum_coverage_cancers
    ) %>%
    arrange(desc(n_coverage), desc(loss_fraction), median_logFC)

  eligible_profiles <- profiles %>% filter(eligible_primary_null)
  if (nrow(eligible_profiles) < same_gene_min_candidates) {
    result <- list(
      summary = empty_summary("insufficient_same_gene_candidates") %>%
        mutate(n_candidate_genes = nrow(eligible_profiles)),
      profiles = profiles
    )
    return(if (return_profiles) result else result$summary)
  }

  observed_logfc <- observed_cells$target_logFC
  observed_n_loss <- sum(observed_logfc <= -effect_logfc_threshold)
  observed_n_gain <- sum(observed_logfc >= effect_logfc_threshold)
  observed_loss_fraction <- observed_n_loss / length(observed_logfc)
  observed_gain_fraction <- observed_n_gain / length(observed_logfc)
  observed_median <- median(observed_logfc)

  summary <- tibble::tibble(
    gene = target_gene,
    minimum_coverage_fraction = minimum_coverage_fraction,
    minimum_coverage_cancers = minimum_coverage_cancers,
    n_cancers = nrow(observed_cells),
    observed_n_loss = observed_n_loss,
    observed_loss_fraction = observed_loss_fraction,
    observed_n_gain = observed_n_gain,
    observed_gain_fraction = observed_gain_fraction,
    observed_median_logFC = observed_median,
    n_candidate_genes = nrow(eligible_profiles),
    null_median_loss_fraction = median(eligible_profiles$loss_fraction),
    null_median_gain_fraction = median(eligible_profiles$gain_fraction),
    null_median_logFC = median(eligible_profiles$median_logFC),
    p_recurrent_loss = (
      1 + sum(eligible_profiles$loss_fraction >= observed_loss_fraction)
    ) / (nrow(eligible_profiles) + 1),
    p_recurrent_gain = (
      1 + sum(eligible_profiles$gain_fraction >= observed_gain_fraction)
    ) / (nrow(eligible_profiles) + 1),
    p_median_negative = (
      1 + sum(eligible_profiles$median_logFC <= observed_median)
    ) / (nrow(eligible_profiles) + 1),
    status = "estimable_same_gene_null"
  )
  result <- list(summary = summary, profiles = profiles)
  if (return_profiles) result else result$summary
}

same_gene_objects <- purrr::map(
  target_genes,
  ~ build_same_gene_null(.x, return_profiles = TRUE)
)
recurrence_null <- purrr::map_dfr(same_gene_objects, "summary") %>%
  mutate(
    q_recurrent_loss = p.adjust(p_recurrent_loss, method = "BH"),
    q_recurrent_gain = p.adjust(p_recurrent_gain, method = "BH"),
    q_median_negative = p.adjust(p_median_negative, method = "BH")
  )
same_gene_profiles <- purrr::map_dfr(same_gene_objects, "profiles")

recurrence_sensitivity <- tidyr::expand_grid(
  gene = target_genes,
  minimum_coverage_fraction = same_gene_sensitivity_coverages
) %>%
  purrr::pmap_dfr(function(gene, minimum_coverage_fraction) {
    build_same_gene_null(
      gene,
      minimum_coverage_fraction = minimum_coverage_fraction,
      return_profiles = FALSE
    )
  }) %>%
  group_by(minimum_coverage_fraction) %>%
  mutate(
    q_recurrent_loss = p.adjust(p_recurrent_loss, method = "BH"),
    q_recurrent_gain = p.adjust(p_recurrent_gain, method = "BH"),
    q_median_negative = p.adjust(p_median_negative, method = "BH")
  ) %>%
  ungroup()

write_csv_atomic(
  recurrence_null,
  file.path(table_dir, "162_round2B_same_gene_recurrence_null.csv")
)
write_csv_atomic(
  same_gene_profiles,
  file.path(table_dir, "163_round2B_same_gene_null_gene_profiles.csv")
)
write_csv_atomic(
  recurrence_sensitivity,
  file.path(table_dir, "164_round2B_same_gene_recurrence_sensitivity.csv")
)

# Figure 3A: observed recurrent-loss fraction against same-gene null profiles.
recurrence_plot_null <- same_gene_profiles %>%
  filter(eligible_primary_null) %>%
  select(target_gene, loss_fraction) %>%
  rename(gene = target_gene) %>%
  mutate(gene = factor(gene, levels = target_genes))
recurrence_plot_observed <- recurrence_null %>%
  filter(status == "estimable_same_gene_null") %>%
  select(gene, observed_loss_fraction, p_recurrent_loss, q_recurrent_loss) %>%
  mutate(gene = factor(gene, levels = target_genes))

plot_recurrence <- ggplot(
  recurrence_plot_null,
  aes(x = gene, y = loss_fraction)
) +
  geom_violin(fill = "grey88", color = "grey45", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white") +
  geom_point(
    data = recurrence_plot_observed,
    aes(x = gene, y = observed_loss_fraction),
    inherit.aes = FALSE, shape = 23, size = 3.4,
    fill = "#B2182B", color = "black"
  ) +
  geom_text(
    data = recurrence_plot_observed,
    aes(
      x = gene, y = pmin(1.08, observed_loss_fraction + 0.08),
      label = paste0("q=", formatC(q_recurrent_loss, digits = 2, format = "g"))
    ),
    inherit.aes = FALSE, size = 3
  ) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  labs(
    x = NULL, y = "Fraction of cancers with loss",
    title = "Cross-cancer recurrence against a covariance-preserving same-gene null",
    subtitle = paste0(
      "Grey distributions: matched genes observed in >=",
      round(100 * same_gene_primary_min_coverage),
      "% of target-evaluable cancers; diamonds: target genes"
    )
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title.position = "plot")
save_plot_pair(
  plot_recurrence,
  "Figure3A_round2B_same_gene_recurrence_null",
  width = 9.2, height = 5.4
)

# Figure 3B: cancer-level target shifts and their percentile within the matched
# background. This retains effect-size context without claiming cell-wise FDR.
selective_plot_data <- matched_background_summary %>%
  filter(grepl("^estimable", matching_status), is.finite(target_logFC)) %>%
  mutate(
    gene = factor(gene, levels = target_genes),
    tcga_code = factor(tcga_code, levels = rev(primary_cancers$tcga_code)),
    percentile_label = paste0(round(100 * background_percentile), "%")
  )

plot_selective <- ggplot(
  selective_plot_data,
  aes(x = target_logFC, y = tcga_code, color = background_percentile)
) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_segment(
    aes(x = background_median_logFC, xend = target_logFC, yend = tcga_code),
    color = "grey65", linewidth = 0.45
  ) +
  geom_point(size = 2.4) +
  scale_color_gradient2(
    low = "#2166AC", mid = "grey75", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    name = "Matched-background\npercentile"
  ) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(
    x = "Paired limma logFC (tumor - adjacent non-tumor)", y = NULL,
    title = "Cancer-level target shifts within organ-specificity-matched backgrounds",
    subtitle = "Segments begin at the matched-background median; color indicates target percentile"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title.position = "plot"
  )

save_plot_pair(
  plot_selective,
  "Figure3B_round2B_matched_background_context",
  width = 12, height = 7.3
)

identity_plot_data <- target_identity_residual %>%
  filter(identity_status == "estimable") %>%
  mutate(
    gene = factor(gene, levels = rev(target_genes)),
    tcga_code = factor(tcga_code, levels = primary_cancers$tcga_code)
  )
plot_identity_residual <- ggplot(
  identity_plot_data,
  aes(x = tcga_code, y = gene, fill = residual_vs_identity)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-3, 3), oob = scales::squish,
    na.value = "grey85", name = "Target logFC -\nidentity median"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Target shifts relative to GTEx-defined tissue-identity loss"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )
save_plot_pair(
  plot_identity_residual,
  "SupplementaryFigure_round2B_tissue_identity_residual",
  width = 10.5, height = 4.8
)

# ---- 6. Normal baseline versus count-level cancer disruption ----------------

baseline_cancer_level <- target_results %>%
  filter(
    analysis_status == "primary_count_level",
    !reserved,
    mapping_valid_for_baseline_model,
    is.finite(logFC)
  ) %>%
  inner_join(
    gtex_baseline %>%
      filter(
        mapping_valid_for_baseline_model,
        gene_symbol %in% target_genes
      ) %>%
      select(
        tcga_code, gene = gene_symbol, baseline_tpm,
        baseline_log2_tpm1, tau, breadth_1tpm,
        matched_gtex_tissues
      ),
    by = c("tcga_code", "gene")
  ) %>%
  select(
    tcga_code, organ, gene, n_pairs, logFC, p_value, q_value_global,
    robust_count_level_change, baseline_tpm, baseline_log2_tpm1,
    tau, breadth_1tpm, matched_gtex_tissues
  )

# Kidney and lung contain several TCGA histologies but only one GTEx organ
# baseline. Collapse histologies to one organ-gene unit to avoid pseudo-replication.
baseline_organ_level <- baseline_cancer_level %>%
  group_by(organ, gene) %>%
  summarise(
    n_cancer_types = n(),
    cancer_types = paste(sort(tcga_code), collapse = ";"),
    total_pairs = sum(n_pairs),
    baseline_tpm = first(baseline_tpm),
    baseline_log2_tpm1 = first(baseline_log2_tpm1),
    tau = first(tau),
    count_level_logFC = median(logFC),
    minimum_logFC = min(logFC),
    maximum_logFC = max(logFC),
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    baseline_z = if (n() >= 2L && stats::sd(baseline_log2_tpm1) > 0) {
      as.numeric(scale(baseline_log2_tpm1))
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

write_csv_atomic(
  baseline_cancer_level,
  file.path(table_dir, "170_round2B_baseline_disruption_cancer_level_count.csv")
)
write_csv_atomic(
  baseline_organ_level,
  file.path(table_dir, "171_round2B_baseline_disruption_organ_level_count.csv")
)

gene_specific_baseline_tests <- baseline_organ_level %>%
  group_by(gene) %>%
  group_modify(~ {
    d <- .x %>%
      filter(is.finite(baseline_log2_tpm1), is.finite(count_level_logFC))
    if (nrow(d) < 5L || n_distinct(d$baseline_log2_tpm1) < 3L) {
      return(tibble::tibble(
        n_organs = nrow(d), rho = NA_real_, p_value = NA_real_,
        status = "insufficient_organs_or_variation"
      ))
    }
    ct <- suppressWarnings(stats::cor.test(
      d$baseline_log2_tpm1,
      d$count_level_logFC,
      method = "spearman",
      exact = FALSE
    ))
    tibble::tibble(
      n_organs = nrow(d),
      rho = unname(ct$estimate),
      p_value = unname(ct$p.value),
      status = "estimable"
    )
  }) %>%
  ungroup() %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

fit_clustered_baseline_model <- function(data) {
  d <- data %>%
    filter(
      is.finite(baseline_z), is.finite(count_level_logFC),
      !is.na(gene), !is.na(organ)
    ) %>%
    droplevels()
  if (nrow(d) < 12L || n_distinct(d$organ) < 5L || n_distinct(d$gene) < 2L) {
    return(list(
      summary = tibble::tibble(
        n_cells = nrow(d), n_organs = n_distinct(d$organ),
        n_genes = n_distinct(d$gene), beta_baseline_z = NA_real_,
        robust_se = NA_real_, t_value = NA_real_, p_value = NA_real_,
        status = "insufficient_data"
      ),
      fit = NULL
    ))
  }

  fit <- stats::lm(count_level_logFC ~ baseline_z + factor(gene), data = d)
  result <- tryCatch({
    vc <- sandwich::vcovCL(fit, cluster = d$organ, type = "HC2")
    ct <- lmtest::coeftest(fit, vcov. = vc)
    tibble::tibble(
      n_cells = nrow(d),
      n_organs = n_distinct(d$organ),
      n_genes = n_distinct(d$gene),
      beta_baseline_z = unname(ct["baseline_z", "Estimate"]),
      robust_se = unname(ct["baseline_z", "Std. Error"]),
      t_value = unname(ct["baseline_z", "t value"]),
      p_value = unname(ct["baseline_z", "Pr(>|t|)"]),
      status = "estimable_cluster_HC2"
    )
  }, error = function(e) {
    tibble::tibble(
      n_cells = nrow(d), n_organs = n_distinct(d$organ),
      n_genes = n_distinct(d$gene), beta_baseline_z = NA_real_,
      robust_se = NA_real_, t_value = NA_real_, p_value = NA_real_,
      status = paste0("model_error:", conditionMessage(e))
    )
  })
  list(summary = result, fit = fit)
}

pooled_baseline_fit <- fit_clustered_baseline_model(baseline_organ_level)
pooled_baseline_model <- pooled_baseline_fit$summary %>%
  mutate(
    model = "logFC ~ within-gene baseline_z + gene fixed effects",
    cluster = "organ",
    interpretation_boundary = paste(
      "Association between normal baseline and cross-sectional cancer shift;",
      "not temporal or causal"
    )
  )

set.seed(23016802)
permutation_data <- baseline_organ_level %>%
  filter(is.finite(baseline_z), is.finite(count_level_logFC))
observed_beta <- pooled_baseline_model$beta_baseline_z[[1]]
permuted_beta <- rep(NA_real_, baseline_permutation_reps)
if (is.finite(observed_beta)) {
  for (b in seq_len(baseline_permutation_reps)) {
    permuted <- permutation_data %>%
      group_by(gene) %>%
      mutate(permuted_logFC = sample(count_level_logFC, size = n(), replace = FALSE)) %>%
      ungroup()
    perm_fit <- tryCatch(
      stats::lm(permuted_logFC ~ baseline_z + factor(gene), data = permuted),
      error = function(e) NULL
    )
    if (!is.null(perm_fit)) {
      permuted_beta[[b]] <- stats::coef(perm_fit)[["baseline_z"]]
    }
  }
}
valid_permuted_beta <- permuted_beta[is.finite(permuted_beta)]
permutation_summary <- tibble::tibble(
  observed_beta = observed_beta,
  requested_permutations = baseline_permutation_reps,
  valid_permutations = length(valid_permuted_beta),
  permutation_p_two_sided = if (
      is.finite(observed_beta) && length(valid_permuted_beta) > 0L
  ) {
    (1 + sum(abs(valid_permuted_beta) >= abs(observed_beta))) /
      (length(valid_permuted_beta) + 1)
  } else {
    NA_real_
  },
  permutation_p_negative = if (
      is.finite(observed_beta) && length(valid_permuted_beta) > 0L
  ) {
    (1 + sum(valid_permuted_beta <= observed_beta)) /
      (length(valid_permuted_beta) + 1)
  } else {
    NA_real_
  }
)

make_sensitivity_data <- function(data, baseline_min_tpm = -Inf) {
  data %>%
    filter(baseline_tpm >= baseline_min_tpm) %>%
    group_by(gene) %>%
    mutate(
      baseline_z = if (n() >= 2L && stats::sd(baseline_log2_tpm1) > 0) {
        as.numeric(scale(baseline_log2_tpm1))
      } else {
        NA_real_
      }
    ) %>%
    ungroup()
}

baseline_sensitivity <- bind_rows(
  fit_clustered_baseline_model(
    make_sensitivity_data(baseline_organ_level, 0)
  )$summary %>% mutate(sensitivity = "all_detectable"),
  fit_clustered_baseline_model(
    make_sensitivity_data(baseline_organ_level, 1)
  )$summary %>% mutate(sensitivity = "normal_baseline_ge_1_TPM")
)

leave_one_organ_out <- purrr::map_dfr(sort(unique(baseline_organ_level$organ)), function(x) {
  fit_clustered_baseline_model(
    make_sensitivity_data(baseline_organ_level %>% filter(organ != x), 0)
  )$summary %>%
    mutate(omitted_type = "organ", omitted_value = x)
})
leave_one_gene_out <- purrr::map_dfr(target_genes, function(x) {
  fit_clustered_baseline_model(
    make_sensitivity_data(baseline_organ_level %>% filter(gene != x), 0)
  )$summary %>%
    mutate(omitted_type = "gene", omitted_value = x)
})
leave_one_out <- bind_rows(leave_one_organ_out, leave_one_gene_out)

write_csv_atomic(
  gene_specific_baseline_tests,
  file.path(table_dir, "172_round2B_baseline_disruption_gene_tests_count.csv")
)
write_csv_atomic(
  pooled_baseline_model,
  file.path(table_dir, "173_round2B_baseline_disruption_pooled_model_count.csv")
)
write_csv_atomic(
  permutation_summary,
  file.path(table_dir, "174_round2B_baseline_disruption_permutation_count.csv")
)
write_csv_atomic(
  baseline_sensitivity,
  file.path(table_dir, "175_round2B_baseline_disruption_sensitivity.csv")
)
write_csv_atomic(
  leave_one_out,
  file.path(table_dir, "176_round2B_baseline_disruption_leave_one_out.csv")
)

plot_baseline_data <- baseline_organ_level %>%
  filter(
    is.finite(baseline_z), is.finite(count_level_logFC),
    !is.na(gene), !is.na(organ)
  )
if (nrow(plot_baseline_data) >= 3L) {
  gene_only_fit <- stats::lm(
    count_level_logFC ~ factor(gene), data = plot_baseline_data
  )
  plot_baseline_data$within_gene_logFC_residual <- stats::residuals(gene_only_fit)
} else {
  plot_baseline_data$within_gene_logFC_residual <- NA_real_
}

baseline_beta_label <- if (is.finite(pooled_baseline_model$beta_baseline_z[[1]])) {
  paste0(
    "cluster-HC2 beta=",
    formatC(pooled_baseline_model$beta_baseline_z[[1]], digits = 3, format = "f"),
    ", p=",
    formatC(pooled_baseline_model$p_value[[1]], digits = 3, format = "g")
  )
} else {
  "pooled model not estimable"
}

plot_baseline <- plot_baseline_data %>%
  ggplot(aes(
    x = baseline_z,
    y = within_gene_logFC_residual,
    color = gene
  )) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_point(size = 2.3, alpha = 0.9) +
  geom_smooth(
    aes(group = 1), method = "lm", se = TRUE,
    color = "black", fill = "grey75", linewidth = 0.7
  ) +
  labs(
    x = "Within-gene standardized GTEx normal baseline",
    y = "Paired logFC residual after removing gene means",
    title = "Exploratory normal-baseline/disruption relationship",
    subtitle = paste0(
      "Organ-level secondary analysis including LIHC paired effects; ",
      baseline_beta_label
    )
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title.position = "plot")

save_plot_pair(
  plot_baseline,
  "SupplementaryFigure_round2B_baseline_vs_count_disruption",
  width = 8.8, height = 6.2
)

# ---- 7. Within-patient normal-to-tumor coexpression rewiring ----------------

safe_spearman_value <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 4L || n_distinct(x) < 3L || n_distinct(y) < 3L) {
    return(NA_real_)
  }
  suppressWarnings(unname(stats::cor(x, y, method = "spearman")))
}

bootstrap_delta_spearman <- function(
    normal_x, normal_y, tumor_x, tumor_y, n_reps, seed_offset) {
  n <- length(normal_x)
  set.seed(23016802 + seed_offset)
  boot_delta <- rep(NA_real_, n_reps)
  for (b in seq_len(n_reps)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    rho_n <- safe_spearman_value(normal_x[idx], normal_y[idx])
    rho_t <- safe_spearman_value(tumor_x[idx], tumor_y[idx])
    boot_delta[[b]] <- rho_t - rho_n
  }
  boot_delta <- boot_delta[is.finite(boot_delta)]
  if (length(boot_delta) < ceiling(0.8 * n_reps)) {
    return(tibble::tibble(
      bootstrap_valid = length(boot_delta), ci_low = NA_real_, ci_high = NA_real_,
      bootstrap_p_two_sided = NA_real_
    ))
  }
  tibble::tibble(
    bootstrap_valid = length(boot_delta),
    ci_low = unname(stats::quantile(boot_delta, 0.025, type = 6)),
    ci_high = unname(stats::quantile(boot_delta, 0.975, type = 6)),
    bootstrap_p_two_sided = min(
      1,
      2 * min(
        (1 + sum(boot_delta <= 0)) / (length(boot_delta) + 1),
        (1 + sum(boot_delta >= 0)) / (length(boot_delta) + 1)
      )
    )
  )
}

coexpression_one_cell <- function(project_analysis, partner, seed_offset) {
  tcga_code <- project_analysis$tcga_code
  n_pairs_available <- project_analysis$n_pairs

  # Paper boundary: no LIHC expression vector may enter this module.
  if (is_excluded_coexpression_cancer(tcga_code)) {
    stop(
      "Paper-boundary violation: LIHC entered the coexpression module.",
      call. = FALSE
    )
  }

  meta <- project_analysis$sample_meta %>%
    mutate(column_index = seq_len(n())) %>%
    select(case_id, condition, sample_barcode, column_index) %>%
    pivot_wider(
      names_from = condition,
      values_from = c(sample_barcode, column_index),
      names_sep = "_"
    ) %>%
    filter(!is.na(column_index_Normal), !is.na(column_index_Tumor)) %>%
    arrange(case_id)

  n_idx <- meta$column_index_Normal
  t_idx <- meta$column_index_Tumor
  expr <- project_analysis$target_logcpm
  raw <- project_analysis$target_counts

  required_rows <- c("MAOB", partner)
  if (!all(required_rows %in% rownames(expr))) {
    return(tibble::tibble(
      tcga_code = tcga_code, partner = partner, n_pairs = nrow(meta),
      detection_normal_MAOB = NA_real_, detection_tumor_MAOB = NA_real_,
      detection_normal_partner = NA_real_, detection_tumor_partner = NA_real_,
      rho_normal = NA_real_, rho_tumor = NA_real_, delta_rho = NA_real_,
      bootstrap_valid = NA_integer_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_, status = "gene_not_found"
    ))
  }

  normal_maob <- expr["MAOB", n_idx]
  tumor_maob <- expr["MAOB", t_idx]
  normal_partner <- expr[partner, n_idx]
  tumor_partner <- expr[partner, t_idx]
  complete <- is.finite(normal_maob) & is.finite(tumor_maob) &
    is.finite(normal_partner) & is.finite(tumor_partner)

  normal_maob <- normal_maob[complete]
  tumor_maob <- tumor_maob[complete]
  normal_partner <- normal_partner[complete]
  tumor_partner <- tumor_partner[complete]
  n_idx_complete <- n_idx[complete]
  t_idx_complete <- t_idx[complete]
  n_pairs <- sum(complete)

  detection <- c(
    normal_maob = mean(raw["MAOB", n_idx_complete] > 0, na.rm = TRUE),
    tumor_maob = mean(raw["MAOB", t_idx_complete] > 0, na.rm = TRUE),
    normal_partner = mean(raw[partner, n_idx_complete] > 0, na.rm = TRUE),
    tumor_partner = mean(raw[partner, t_idx_complete] > 0, na.rm = TRUE)
  )

  status <- case_when(
    n_pairs < minimum_pairs_exploratory ~ "insufficient_n_lt20",
    any(!is.finite(detection)) ~ "missing_detection_information",
    min(detection) < minimum_detection_fraction ~ "low_detection",
    any(vapply(
      list(normal_maob, tumor_maob, normal_partner, tumor_partner),
      function(x) n_distinct(x) < 3L || stats::sd(x) == 0,
      logical(1)
    )) ~ "low_variability",
    n_pairs < minimum_pairs_coexpression ~ "exploratory_n_20_29",
    TRUE ~ "primary_n_ge_30"
  )

  if (!status %in% c("primary_n_ge_30", "exploratory_n_20_29")) {
    return(tibble::tibble(
      tcga_code = tcga_code, partner = partner, n_pairs = n_pairs,
      detection_normal_MAOB = detection[["normal_maob"]],
      detection_tumor_MAOB = detection[["tumor_maob"]],
      detection_normal_partner = detection[["normal_partner"]],
      detection_tumor_partner = detection[["tumor_partner"]],
      rho_normal = NA_real_, rho_tumor = NA_real_, delta_rho = NA_real_,
      bootstrap_valid = NA_integer_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_, status = status
    ))
  }

  rho_normal <- safe_spearman_value(normal_maob, normal_partner)
  rho_tumor <- safe_spearman_value(tumor_maob, tumor_partner)
  delta_rho <- rho_tumor - rho_normal
  boot <- bootstrap_delta_spearman(
    normal_maob, normal_partner, tumor_maob, tumor_partner,
    n_reps = coexpression_bootstrap_reps,
    seed_offset = seed_offset
  )

  tibble::tibble(
    tcga_code = tcga_code, partner = partner, n_pairs = n_pairs,
    detection_normal_MAOB = detection[["normal_maob"]],
    detection_tumor_MAOB = detection[["tumor_maob"]],
    detection_normal_partner = detection[["normal_partner"]],
    detection_tumor_partner = detection[["tumor_partner"]],
    rho_normal = rho_normal, rho_tumor = rho_tumor, delta_rho = delta_rho,
    bootstrap_valid = as.integer(boot$bootstrap_valid[[1]]),
    ci_low = boot$ci_low[[1]], ci_high = boot$ci_high[[1]],
    p_value = boot$bootstrap_p_two_sided[[1]], status = status
  )
}

coexpression_project_analyses <- project_analyses[
  !vapply(names(project_analyses), is_excluded_coexpression_cancer, logical(1))
]

coexpression_results <- purrr::imap_dfr(coexpression_project_analyses, function(x, code) {
  purrr::map2_dfr(
    cyp_partners,
    seq_along(cyp_partners),
    function(partner, j) {
      code_index <- match(code, names(coexpression_project_analyses))
      coexpression_one_cell(
        x,
        partner,
        seed_offset = 1000L * code_index + j
      )
    }
  )
}) %>%
  left_join(primary_cancers %>% select(tcga_code, organ), by = "tcga_code")

coexpression_results$q_value_primary <- NA_real_
eligible_coexpression <- coexpression_results$status == "primary_n_ge_30" &
  is.finite(coexpression_results$p_value)
coexpression_results$q_value_primary[eligible_coexpression] <- p.adjust(
  coexpression_results$p_value[eligible_coexpression],
  method = "BH"
)
coexpression_results <- coexpression_results %>%
  mutate(
    passes_coexpression_effect = is.finite(delta_rho) &
      abs(delta_rho) >= coexpression_effect_threshold,
    robust_rewiring = status == "primary_n_ge_30" &
      !is.na(q_value_primary) & q_value_primary < 0.05 &
      passes_coexpression_effect,
    rewiring_direction = case_when(
      robust_rewiring & delta_rho > 0 ~ "more_positive_in_tumor",
      robust_rewiring & delta_rho < 0 ~ "more_negative_in_tumor",
      status %in% c("primary_n_ge_30", "exploratory_n_20_29") ~ "not_robust",
      TRUE ~ "not_estimable"
    )
  )

if (any(coexpression_results$tcga_code == "LIHC")) {
  stop("Paper-boundary failure: a LIHC coexpression result was created.", call. = FALSE)
}

write_csv_atomic(
  coexpression_results,
  file.path(table_dir, "180_round2B_nonLIHC_coexpression_rewiring.csv")
)

coexpression_plot_data <- coexpression_results %>%
  mutate(
    partner = factor(partner, levels = rev(cyp_partners)),
    tcga_code = factor(
      tcga_code,
      levels = primary_cancers$tcga_code[
        !vapply(
          primary_cancers$tcga_code,
          is_excluded_coexpression_cancer,
          logical(1)
        )
      ]
    ),
    label = case_when(
      robust_rewiring ~ "*",
      status == "exploratory_n_20_29" ~ "E",
      TRUE ~ ""
    )
  )

plot_coexpression <- ggplot(
  coexpression_plot_data,
  aes(x = tcga_code, y = partner, fill = delta_rho)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#762A83", mid = "white", high = "#1B7837", midpoint = 0,
    limits = c(-1, 1), oob = scales::squish, na.value = "grey85",
    name = "delta rho\n(tumor - normal)"
  ) +
  labs(
    x = NULL, y = "MAOB partner",
    title = "Within-patient MAOB-CYP3A coexpression rewiring",
    subtitle = paste0(
      "Non-LIHC cancers only; *=primary FDR<0.05 and |delta rho|>=0.20; ",
      "E=exploratory n=20-29"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )

save_plot_pair(
  plot_coexpression,
  "Figure4A_round2B_nonLIHC_coexpression_heatmap",
  width = 10.5, height = 4.6
)

# Representative rank-scatter panels: select the largest absolute robust delta
# within BRCA, KIRC, and KIRP. Ranks align the visual with Spearman correlation.
representative_coexpression_cells <- coexpression_results %>%
  filter(
    robust_rewiring,
    tcga_code %in% c("BRCA", "KIRC", "KIRP")
  ) %>%
  group_by(tcga_code) %>%
  slice_max(abs(delta_rho), n = 1L, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(tcga_code)

extract_coexpression_points <- function(tcga_code, partner) {
  # Use distinct scalar names inside tibble(). Otherwise the newly created
  # `partner` column masks the function argument during sequential evaluation,
  # turning a one-row matrix index into an n-by-n extraction.
  code_value <- as.character(tcga_code[[1]])
  partner_value <- as.character(partner[[1]])
  analysis <- coexpression_project_analyses[[code_value]]
  if (is.null(analysis)) return(tibble::tibble())
  meta <- analysis$sample_meta %>%
    mutate(column_index = seq_len(n()))
  expr <- analysis$target_logcpm
  if (!all(c("MAOB", partner_value) %in% rownames(expr))) {
    return(tibble::tibble())
  }

  tibble::tibble(
    tcga_code = code_value,
    partner = partner_value,
    case_id = meta$case_id,
    condition = as.character(meta$condition),
    MAOB_logCPM = as.numeric(expr["MAOB", meta$column_index]),
    partner_logCPM = as.numeric(expr[partner_value, meta$column_index])
  ) %>%
    filter(is.finite(MAOB_logCPM), is.finite(partner_logCPM)) %>%
    group_by(tcga_code, partner, condition) %>%
    mutate(
      MAOB_rank_fraction = rank(MAOB_logCPM, ties.method = "average") / n(),
      partner_rank_fraction = rank(partner_logCPM, ties.method = "average") / n()
    ) %>%
    ungroup()
}

representative_coexpression_points <- if (
    nrow(representative_coexpression_cells) > 0L
) {
  purrr::pmap_dfr(
    representative_coexpression_cells %>% select(tcga_code, partner),
    extract_coexpression_points
  ) %>%
    left_join(
      representative_coexpression_cells %>%
        select(tcga_code, partner, rho_normal, rho_tumor, delta_rho),
      by = c("tcga_code", "partner")
    ) %>%
    mutate(
      panel = paste0(
        tcga_code, ": MAOB-", partner,
        "\nnormal rho=", formatC(rho_normal, digits = 2, format = "f"),
        ", tumor rho=", formatC(rho_tumor, digits = 2, format = "f"),
        ", delta=", formatC(delta_rho, digits = 2, format = "f")
      ),
      rho_label = if_else(
        condition == "Normal",
        paste0("Normal rho=", formatC(rho_normal, digits = 2, format = "f")),
        paste0("Tumor rho=", formatC(rho_tumor, digits = 2, format = "f"))
      )
    )
} else {
  tibble::tibble()
}

write_csv_atomic(
  representative_coexpression_cells,
  file.path(table_dir, "181_round2B_representative_coexpression_cells.csv")
)
write_csv_atomic(
  representative_coexpression_points,
  file.path(table_dir, "182_round2B_representative_coexpression_points.csv")
)

if (nrow(representative_coexpression_points) > 0L) {
  plot_representative_coexpression <- ggplot(
    representative_coexpression_points,
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
      title = "Representative non-LIHC coexpression rewiring",
      subtitle = "Rank-scale visualization corresponding to the Spearman analysis"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title.position = "plot"
    )
  save_plot_pair(
    plot_representative_coexpression,
    "Figure4B_round2B_representative_coexpression_rank_scatter",
    width = 10.5, height = 6.8
  )
}

# ---- 8. HPA-only normal-tissue validation -----------------------------------

hpa_url <- "https://www.proteinatlas.org/download/tsv/rna_tissue_hpa.tsv.zip"
hpa_zip <- file.path(raw_dir, "rna_tissue_hpa.tsv.zip")
hpa_tsv <- file.path(raw_dir, "rna_tissue_hpa.tsv")

hpa_validation_status <- tryCatch({
  if (!file.exists(hpa_zip)) {
    log_message("Downloading HPA-only RNA tissue data")
    utils::download.file(hpa_url, hpa_zip, mode = "wb", quiet = FALSE)
  }
  if (!file.exists(hpa_tsv)) {
    extracted <- utils::unzip(hpa_zip, exdir = raw_dir)
    if (!file.exists(hpa_tsv)) {
      candidate <- extracted[grepl("rna_tissue_hpa.*\\.tsv$", extracted)]
      if (length(candidate) == 1L) hpa_tsv <- candidate[[1]]
    }
  }
  if (!file.exists(hpa_tsv)) stop("HPA-only TSV was not found after unzip.")

  hpa <- data.table::fread(hpa_tsv, data.table = FALSE, check.names = FALSE)
  names(hpa) <- normalise_names(names(hpa))
  gene_col <- first_existing_column(hpa, c("gene_name", "gene_symbol"))
  ensembl_col <- first_existing_column(hpa, c("gene", "gene_id"))
  tissue_col <- first_existing_column(hpa, c("tissue", "tissue_name"))
  expression_col <- first_existing_column(hpa, c("ntpm", "ptpm", "tpm"))

  if (is.na(ensembl_col) || is.na(tissue_col) || is.na(expression_col)) {
    stop(
      "HPA-only file has unexpected columns: ",
      paste(names(hpa), collapse = ", ")
    )
  }

  hpa_gene_id <- strip_ensembl_version(hpa[[ensembl_col]])
  hpa_gene_symbol <- if (!is.na(gene_col)) {
    as.character(hpa[[gene_col]])
  } else {
    gtex_metrics$gene_symbol[match(hpa_gene_id, gtex_metrics$gene_id)]
  }

  hpa_targets <- tibble::tibble(
    gene_id = hpa_gene_id,
    gene = hpa_gene_symbol,
    tissue = as.character(hpa[[tissue_col]]),
    hpa_expression = safe_numeric(hpa[[expression_col]])
  ) %>%
    filter(gene %in% target_genes, is.finite(hpa_expression)) %>%
    group_by(gene, tissue) %>%
    summarise(
      n_hpa_rows = n(),
      hpa_expression = median(hpa_expression),
      .groups = "drop"
    ) %>%
    arrange(gene, desc(hpa_expression))

  if (nrow(hpa_targets) == 0L) stop("No target genes found in HPA-only data.")
  write_csv_atomic(
    hpa_targets,
    file.path(table_dir, "190_round2B_HPA_only_RNA_target_tissues.csv")
  )

  hpa_top_tissues <- hpa_targets %>%
    group_by(gene) %>%
    arrange(desc(hpa_expression), tissue, .by_group = TRUE) %>%
    mutate(tissue_rank = row_number()) %>%
    filter(tissue_rank <= 10L) %>%
    ungroup()
  write_csv_atomic(
    hpa_top_tissues,
    file.path(table_dir, "191_round2B_HPA_only_top10_tissues.csv")
  )

  hpa_organ_patterns <- tibble::tribble(
    ~organ, ~hpa_pattern,
    "adrenal gland", "adrenal",
    "brain", "brain|cerebr|cortex|hippocamp|amygdala|hypothalam|putamen|caudate",
    "breast", "breast",
    "colon", "colon",
    "esophagus", "esophag",
    "kidney", "kidney",
    "liver", "liver",
    "lung", "lung",
    "ovary", "ovary",
    "pancreas", "pancreas",
    "prostate", "prostate",
    "skin", "skin",
    "stomach", "stomach",
    "testis", "testis",
    "thyroid", "thyroid",
    "uterus", "uterus|endometrium"
  )

  hpa_organ_targets <- purrr::pmap_dfr(
    hpa_organ_patterns,
    function(organ, hpa_pattern) {
      hpa_targets %>%
        filter(str_detect(tolower(tissue), hpa_pattern)) %>%
        group_by(gene) %>%
        summarise(
          hpa_n_tissues = n_distinct(tissue),
          hpa_organ_expression = median(hpa_expression),
          matched_hpa_tissues = paste(sort(unique(tissue)), collapse = ";"),
          .groups = "drop"
        ) %>%
        mutate(organ = organ, .before = 1L)
    }
  )

  gtex_organ_targets <- gtex_baseline %>%
    filter(
      mapping_valid_for_baseline_model,
      gene_symbol %in% target_genes
    ) %>%
    group_by(organ, gene = gene_symbol) %>%
    summarise(
      gtex_baseline_tpm = median(baseline_tpm),
      .groups = "drop"
    )

  hpa_gtex_validation <- gtex_organ_targets %>%
    inner_join(hpa_organ_targets, by = c("organ", "gene"))
  write_csv_atomic(
    hpa_gtex_validation,
    file.path(table_dir, "192_round2B_HPA_only_GTEx_organ_concordance_data.csv")
  )

  hpa_gtex_tests <- hpa_gtex_validation %>%
    group_by(gene) %>%
    group_modify(~ {
      d <- .x %>%
        filter(is.finite(gtex_baseline_tpm), is.finite(hpa_organ_expression))
      if (nrow(d) < 5L || n_distinct(d$gtex_baseline_tpm) < 3L ||
          n_distinct(d$hpa_organ_expression) < 3L) {
        return(tibble::tibble(
          n_organs = nrow(d), rho = NA_real_, p_value = NA_real_,
          status = "insufficient_organs_or_variation"
        ))
      }
      ct <- suppressWarnings(stats::cor.test(
        d$gtex_baseline_tpm,
        d$hpa_organ_expression,
        method = "spearman",
        exact = FALSE
      ))
      tibble::tibble(
        n_organs = nrow(d), rho = unname(ct$estimate),
        p_value = unname(ct$p.value), status = "descriptive_validation"
      )
    }) %>%
    ungroup() %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))
  write_csv_atomic(
    hpa_gtex_tests,
    file.path(table_dir, "193_round2B_HPA_only_GTEx_organ_concordance_tests.csv")
  )

  tissue_order <- hpa_targets %>%
    group_by(tissue) %>%
    summarise(max_expression = max(hpa_expression), .groups = "drop") %>%
    arrange(max_expression) %>%
    pull(tissue)
  plot_hpa_data <- hpa_targets %>%
    mutate(
      gene = factor(gene, levels = rev(target_genes)),
      tissue = factor(tissue, levels = tissue_order)
    )
  plot_hpa <- ggplot(
    plot_hpa_data,
    aes(x = tissue, y = gene, fill = log2(hpa_expression + 1))
  ) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "C", name = "log2(HPA RNA+1)") +
    labs(
      x = NULL, y = NULL,
      title = "HPA-only normal-tissue RNA validation",
      subtitle = paste0(
        "Independent HPA tissue dataset; expression column used: ", expression_col
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1),
      plot.title.position = "plot"
    )
  save_plot_pair(
    plot_hpa,
    "SupplementaryFigure_round2B_HPA_only_validation",
    width = 13, height = 4.7
  )

  tibble::tibble(
    status = "completed",
    source = "HPA-only RNA tissue dataset",
    expression_column = expression_col,
    n_target_rows = nrow(hpa_targets),
    n_tissues = n_distinct(hpa_targets$tissue),
    message = NA_character_
  )
}, error = function(e) {
  warning("HPA-only validation skipped: ", conditionMessage(e))
  tibble::tibble(
    status = "skipped_with_error",
    source = "HPA-only RNA tissue dataset",
    expression_column = NA_character_,
    n_target_rows = NA_integer_, n_tissues = NA_integer_,
    message = conditionMessage(e)
  )
})
write_csv_atomic(
  hpa_validation_status,
  file.path(table_dir, "194_round2B_HPA_only_validation_status.csv")
)

# ---- 9. Final integrity checks and reproducibility log -----------------------

expected_target_cells <- nrow(primary_cancers) * length(target_genes)
if (nrow(target_results) != expected_target_cells) {
  stop(
    "Target result grid is incomplete: expected ", expected_target_cells,
    ", found ", nrow(target_results), call. = FALSE
  )
}
if (anyDuplicated(target_results[c("tcga_code", "gene")])) {
  stop("Duplicate cancer-gene cells in target results.", call. = FALSE)
}
if (any(pair_qc$n_pairs < 2L) || any(pair_qc$design_rank != pair_qc$design_columns)) {
  stop("Pair-count or design-rank integrity check failed.", call. = FALSE)
}

lihc_allgene <- all_gene_de %>% filter(tcga_code == "LIHC")
if ("LIHC" %in% primary_cancers$tcga_code) {
  lihc_returned_targets <- target_results %>%
    filter(
      tcga_code == "LIHC",
      gene %in% lihc_expression_genes_returned_to_paper1
    )
  if (nrow(lihc_returned_targets) != length(lihc_expression_genes_returned_to_paper1) ||
      any(lihc_returned_targets$analysis_status == "reserved_not_recomputed") ||
      !all(lihc_expression_genes_returned_to_paper1 %in% lihc_allgene$gene_symbol)) {
    stop(
      "LIHC MAOB/CYP3A4 were not successfully returned to the paired fit.",
      call. = FALSE
    )
  }
}

expected_coexpression_cells <- sum(
  !vapply(
    primary_cancers$tcga_code,
    is_excluded_coexpression_cancer,
    logical(1)
  )
) * length(cyp_partners)
if (nrow(coexpression_results) != expected_coexpression_cells ||
    any(coexpression_results$tcga_code == "LIHC")) {
  stop("Non-LIHC coexpression grid integrity check failed.", call. = FALSE)
}
if (!"organ_specificity_log2" %in% names(matched_background_genes)) {
  stop("Organ specificity was not included in matched-background output.", call. = FALSE)
}

integrity_summary <- tibble::tibble(
  check = c(
    "selected_projects", "target_grid_complete", "unique_target_cells",
    "paired_designs_full_rank", "LIHC_MAOB_CYP3A4_returned_to_paired_fit",
    "entire_LIHC_coexpression_module_excluded",
    "organ_specificity_in_background_matching",
    "same_gene_recurrence_null_created", "GTEx_direct_source_test",
    "HPA_dataset"
  ),
  result = c(
    paste(primary_cancers$tcga_code, collapse = ";"),
    as.character(nrow(target_results) == expected_target_cells),
    as.character(!anyDuplicated(target_results[c("tcga_code", "gene")])),
    as.character(all(pair_qc$design_rank == pair_qc$design_columns)),
    ifelse("LIHC" %in% primary_cancers$tcga_code, "passed", "not_applicable_subset"),
    as.character(!any(coexpression_results$tcga_code == "LIHC")),
    as.character("organ_specificity_log2" %in% names(matched_background_genes)),
    as.character(nrow(recurrence_null) == length(target_genes)),
    "not_run",
    hpa_validation_status$status[[1]]
  )
)
write_csv_atomic(
  integrity_summary,
  file.path(table_dir, "199_round2B_integrity_checks.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "round2B_sessionInfo.txt")
)

subset_run <- nzchar(requested_codes_text)
completion_report <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Completed: ", timestamp_utc()),
  paste0("Project directory: ", project_dir),
  paste0("Run type: ", ifelse(subset_run, "user-requested subset", "full Tier1")),
  paste0("Cancer projects: ", paste(primary_cancers$tcga_code, collapse = ", ")),
  paste0("Total paired cases: ", sum(pair_qc$n_pairs)),
  paste0("Target genes: ", paste(target_genes, collapse = ", ")),
  "Primary DE: GDC STAR unstranded counts, TMM, filterByExpr, paired limma-voom.",
  "Target statistical multiplicity: BH across all estimable target cells, including LIHC MAOB/CYP3A4.",
  paste(
    "Selective disruption: matched on GTEx baseline, tau, breadth, organ specificity,",
    "and TCGA adjacent-normal abundance."
  ),
  paste0(
    "Recurrence null: same background gene retained across cancers; primary minimum coverage ",
    round(100 * same_gene_primary_min_coverage), "% and coverage sensitivity analysis."
  ),
  "Baseline model: count-level logFC, organ aggregation, gene fixed effects, organ-cluster HC2.",
  "Coexpression: non-LIHC cancers only; patient bootstrap within the same paired cases.",
  paste0("HPA-only validation: ", hpa_validation_status$status[[1]]),
  "Paper 1 boundary: LIHC MAOB/CYP3A4 paired expression shifts are included in pan-cancer analyses.",
  "Paper 2 boundary: the entire LIHC coexpression module was excluded before vector access.",
  paste(
    "No LIHC subgroup, survival, Cox, interaction, clinical adjustment, HNF4A,",
    "HCCDB18, ICGC-LIRI-JP, pathway, immune, or purity analysis was run."
  ),
  "No GTEx-versus-TCGA direct hypothesis test was run.",
  "No causal, temporal, GGA-abundance, or enzyme-activity inference was made."
)
writeLines(
  completion_report,
  file.path(log_dir, "round2B_completion_report.txt")
)

log_message("Round 2B complete.")
log_message("Tables:  ", table_dir)
log_message("Figures: ", figure_dir)
log_message("Logs:    ", log_dir)
