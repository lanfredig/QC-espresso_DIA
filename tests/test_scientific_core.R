###############################################################################
# File: tests/test_scientific_core.R
# Purpose:
#   Execute the validated QC-espresso scientific modules in their intended
#   sequence using the article example, then run every current focused suite.
#
# Usage:
#   Rscript --no-init-file tests/test_scientific_core.R
###############################################################################

resolve_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- resolve_root()
setwd(ROOT)

results <- list()
record <- function(name, success, details = "") {
  results[[length(results) + 1L]] <<- data.frame(
    Test = name,
    Success = isTRUE(success),
    Details = details,
    stringsAsFactors = FALSE
  )
  cat(sprintf("%-54s %s\n", name, if (isTRUE(success)) "PASS" else "FAIL"))
  if (!isTRUE(success) && nzchar(details)) cat("  Reason: ", details, "\n", sep = "")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(microseq)
  library(scales)
  library(ggrepel)
  library(RColorBrewer)
})

PROJECT_ROOT <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
params <- list(
  report_path = file.path(PROJECT_ROOT, "example-data", "dataset_4_contaminans.parquet"),
  fasta_path = file.path(PROJECT_ROOT, "example-data", "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta"),
  output_dir = file.path(tempdir(), paste0("qc_espresso_integrated_", Sys.getpid()))
)

dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)

plot_file <- file.path(params$output_dir, "integrated_plots.pdf")
grDevices::pdf(plot_file, width = 12, height = 8)
on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

source(file.path("R", "globals.R"))
source(file.path("R", "io_load.R"))

# Defensive defaults used by plotting modules when not already defined globally.
if (!exists("PANELS_PER_ROW")) PANELS_PER_ROW <- 3L
if (!exists("AR_DENSITY")) AR_DENSITY <- 1
if (!exists("AR_BOXPLOT")) AR_BOXPLOT <- 1
if (!exists("legend_fill", mode = "function")) legend_fill <- function() ggplot2::theme()
if (!exists("legend_color", mode = "function")) legend_color <- function() ggplot2::theme()

modules <- c(
  "qc_identification.R",
  "qc_abundance.R",
  "qc_chromatography.R",
  "qc_chromatography_internalstd.R",
  "qc_completeness.R",
  "qc_coverage.R",
  "qc_diann_quality_metrics.R",
  "qc_intensity.R",
  "qc_ion_mobility.R",
  "qc_mass_error.R",
  "qc_peptides.R",
  "qc_figures_of_merit.R"
)

cat("\n============================================================\n")
cat("INTEGRATED SCIENTIFIC MODULE EXECUTION\n")
cat("============================================================\n")

for (module in modules) {
  module_path <- file.path("R", module)
  outcome <- tryCatch(
    {
      sys.source(module_path, envir = globalenv())
      list(success = TRUE, details = "")
    },
    error = function(e) list(success = FALSE, details = conditionMessage(e))
  )
  record(module, outcome$success, outcome$details)
}

# Integration invariants.
invariants <- list(
  "Filtered report remains available" = function() nrow(report) == 28748L,
  "Sixteen runs remain available" = function() dplyr::n_distinct(report$Run) == 16L,
  "Peptide evidence remains available" = function() nrow(peptides_clean) == 1872L,
  "FASTA reference remains available" = function() nrow(fasta_ref) == 20416L,
  "Completeness global table created" = function() exists("completeness_global_df") && nrow(completeness_global_df) == 5L,
  "Completeness per-run table created" = function() exists("completeness_per_run_df") && nrow(completeness_per_run_df) == 16L,
  "Coverage table created" = function() exists("coverage_results") && nrow(coverage_results) > 0L,
  "Figures of Merit created" = function() exists("figures_of_merit_df") && nrow(figures_of_merit_df) == 16L,
  "Figures of Merit CSV created" = function() exists("figures_of_merit_csv") && file.exists(figures_of_merit_csv)
)

for (name in names(invariants)) {
  outcome <- tryCatch(
    isTRUE(invariants[[name]]()),
    error = function(e) FALSE
  )
  record(name, outcome, if (outcome) "" else "Integration invariant was not satisfied")
}

grDevices::dev.off()

# Run all existing focused suites as child R processes.
focused_suites <- c(
  "run_tests.R",
  "test_multiprotein_abundance.R",
  "test_missed_cleavages.R",
  "test_optional_sections.R",
  "test_completeness.R",
  "test_coverage.R",
  "test_intensity.R"
)

cat("\n============================================================\n")
cat("FOCUSED TEST SUITES\n")
cat("============================================================\n")

rscript <- file.path(R.home("bin"), "Rscript")
for (suite in focused_suites) {
  suite_path <- file.path("tests", suite)
  if (!file.exists(suite_path)) {
    record(paste0("Focused suite: ", suite), FALSE, "Test file was not found")
    next
  }

  log_file <- tempfile(fileext = ".log")
  status <- system2(
    rscript,
    args = c("--no-init-file", shQuote(suite_path)),
    stdout = log_file,
    stderr = log_file
  )
  log_text <- paste(readLines(log_file, warn = FALSE), collapse = "\n")
  unlink(log_file)

  record(
    paste0("Focused suite: ", suite),
    identical(status, 0L),
    if (identical(status, 0L)) "" else tail(strsplit(log_text, "\n", fixed = TRUE)[[1]], 12) |> paste(collapse = "\n")
  )
}

result_table <- do.call(rbind, results)
passed <- sum(result_table$Success)
failed <- sum(!result_table$Success)

output_dir <- file.path("tests", "test-results", "latest")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  result_table,
  file.path(output_dir, "scientific-core-results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

cat("\n============================================================\n")
cat("SCIENTIFIC CORE TEST SUMMARY\n")
cat("============================================================\n")
cat("Tests passed: ", passed, "\n", sep = "")
cat("Tests failed: ", failed, "\n", sep = "")
cat("FINAL STATUS: ", if (failed == 0L) "PASS" else "FAIL", "\n", sep = "")
cat("============================================================\n")

quit(status = if (failed == 0L) 0L else 1L)
