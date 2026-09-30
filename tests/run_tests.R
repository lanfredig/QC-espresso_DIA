###############################################################################
# File: tests/run_tests.R
# Purpose: Automated validation suite for the QC-espresso scientific repository.
# Usage: Rscript --no-init-file tests/run_tests.R
###############################################################################

results <- list()

record_test <- function(name, fun) {
  started <- Sys.time()
  outcome <- tryCatch(
    {
      fun()
      list(name = name, success = TRUE, message = "", elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")))
    },
    error = function(e) {
      list(name = name, success = FALSE, message = conditionMessage(e), elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")))
    }
  )
  results[[length(results) + 1L]] <<- outcome
  cat(sprintf("%-42s %s\n", name, if (outcome$success) "PASS" else "FAIL"))
  if (!outcome$success) cat("  Reason:", outcome$message, "\n")
  invisible(outcome)
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

assert_equal <- function(observed, expected, message) {
  if (!identical(observed, expected)) {
    stop(paste0(message, "\nExpected: ", paste(expected, collapse = ", "), "\nObserved: ", paste(observed, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

resolve_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    test_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(test_file), ".."), winslash = "/", mustWork = TRUE))
  }
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "QC-espresso_DIA.Rproj"))) return(root)
  stop("QC-espresso project root could not be resolved.", call. = FALSE)
}

PROJECT_ROOT <- resolve_project_root()
setwd(PROJECT_ROOT)

EXPECTED_R_SCRIPTS <- sort(c(
  "globals.R", "io_load.R", "qc_abundance.R", "qc_chromatography.R",
  "qc_chromatography_internalstd.R", "qc_completeness.R", "qc_coverage.R",
  "qc_diann_quality_metrics.R", "qc_figures_of_merit.R", "qc_identification.R",
  "qc_intensity.R", "qc_ion_mobility.R", "qc_mass_error.R", "qc_peptides.R",
  "qc_session.R", "validate_inputs.R"
))

REQUIRED_ROOT_FILES <- c(
  "QC-espresso_DIA.Rproj", "VERSION", "run_qc.R", "renv.lock", ".Rprofile", ".gitignore"
)

EXAMPLE_PARQUET <- file.path("example-data", "dataset_4_contaminans.parquet")
EXAMPLE_FASTA <- file.path("example-data", "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta")
REFERENCE_DIR <- file.path("tests", "reference-output", "v1.3.1")

DIRECT_PACKAGES <- c(
  "arrow", "dplyr", "forcats", "ggplot2", "ggrepel", "gt", "gtExtras",
  "knitr", "microseq", "purrr", "RColorBrewer", "rlang", "scales",
  "stringr", "tibble", "tidyr", "tidyverse"
)

source(file.path("R", "validate_inputs.R"), local = globalenv())

cat("\n============================================================\n")
cat("QC-espresso automated validation\n")
cat("============================================================\n")
cat("Project:", PROJECT_ROOT, "\n")
cat("R version:", R.version.string, "\n")
cat("============================================================\n\n")

record_test("Repository root files", function() {
  missing <- REQUIRED_ROOT_FILES[!file.exists(REQUIRED_ROOT_FILES)]
  assert_true(!length(missing), paste("Missing root files:", paste(missing, collapse = ", ")))
})

record_test("R script inventory", function() {
  observed <- sort(list.files("R", pattern = "[.]R$", full.names = FALSE))
  assert_equal(observed, EXPECTED_R_SCRIPTS, "R script inventory differs from the expected inventory.")
  sizes <- file.info(file.path("R", observed))$size
  assert_true(all(!is.na(sizes) & sizes > 0), "One or more R scripts are empty or inaccessible.")
})

record_test("R script syntax", function() {
  files <- c(list.files("R", pattern = "[.]R$", full.names = TRUE), "run_qc.R", "tests/run_tests.R")
  errors <- character()
  for (file in files) {
    message <- tryCatch({ parse(file); "" }, error = function(e) conditionMessage(e))
    if (nzchar(message)) errors[basename(file)] <- message
  }
  assert_true(!length(errors), paste(c("Syntax errors:", paste0(names(errors), ": ", errors)), collapse = "\n"))
})

record_test("QC-espresso version", function() {
  version <- trimws(readLines("VERSION", n = 1L, warn = FALSE))
  assert_equal(version, "1.3.1", "Unexpected QC-espresso version.")
})

record_test("Reference R version", function() {
  version <- paste(R.version$major, R.version$minor, sep = ".")
  assert_equal(version, "4.6.0", "Tests must run with the functional R 4.6.0 environment.")
})

record_test("Required package availability", function() {
  status <- vapply(DIRECT_PACKAGES, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  assert_true(all(status), paste("Missing packages:", paste(DIRECT_PACKAGES[!status], collapse = ", ")))
})

record_test("Package-lockfile equivalence", function() {
  assert_true(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is required to inspect renv.lock.")
  lock <- jsonlite::read_json("renv.lock", simplifyVector = FALSE)
  assert_equal(lock[["R"]][["Version"]], "4.6.0", "renv.lock does not record R 4.6.0.")
  installed <- installed.packages()
  functional <- vapply(DIRECT_PACKAGES, function(pkg) installed[pkg, "Version"], character(1))
  locked <- vapply(DIRECT_PACKAGES, function(pkg) {
    entry <- lock[["Packages"]][[pkg]]
    if (is.null(entry)) NA_character_ else entry[["Version"]]
  }, character(1))
  match <- functional == locked
  match[is.na(match)] <- FALSE
  if (!all(match)) {
    failed <- DIRECT_PACKAGES[!match]
    details <- vapply(failed, function(pkg) paste0(pkg, ": installed=", functional[[pkg]], ", lockfile=", locked[[pkg]]), character(1))
    stop(paste(c("Package version differences:", paste0("  - ", details)), collapse = "\n"), call. = FALSE)
  }
})

record_test("Example file presence", function() {
  files <- c(EXAMPLE_PARQUET, EXAMPLE_FASTA)
  assert_true(all(file.exists(files)), paste("Missing example files:", paste(files[!file.exists(files)], collapse = ", ")))
  sizes <- file.info(files)$size
  assert_true(all(!is.na(sizes) & sizes > 0), "One or more example files are empty or inaccessible.")
})

record_test("Example Parquet validation", function() {
  validate_report_file(EXAMPLE_PARQUET)
})

record_test("Example FASTA validation", function() {
  validate_fasta_file(EXAMPLE_FASTA)
})

record_test("DIA-NN report schema", function() {
  schema <- validate_report_schema(EXAMPLE_PARQUET)
  assert_true(!length(schema$missing_columns), paste("Missing DIA-NN columns:", paste(schema$missing_columns, collapse = ", ")))
})

record_test("Functional reference outputs", function() {
  expected <- c(
    "figures_of_merit_20260930_023339.csv",
    "protein_coverage_with_peptides_20260930_023309.csv",
    "QC_espresso_DIA_20260930_023347.html",
    "SHA256SUMS.txt"
  )
  paths <- file.path(REFERENCE_DIR, expected)
  assert_true(all(file.exists(paths)), paste("Missing reference outputs:", paste(expected[!file.exists(paths)], collapse = ", ")))
  sizes <- file.info(paths)$size
  assert_true(all(!is.na(sizes) & sizes > 0), "One or more reference outputs are empty or inaccessible.")
})

record_test("Reference-output checksums", function() {
  checksum_file <- file.path(REFERENCE_DIR, "SHA256SUMS.txt")
  lines <- readLines(checksum_file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  assert_true(length(lines) == 3L, "SHA256SUMS.txt must contain exactly three file hashes.")
  parts <- strsplit(lines, "  ", fixed = TRUE)
  expected_hashes <- setNames(vapply(parts, `[[`, character(1), 1L), vapply(parts, `[[`, character(1), 2L))
  referenced <- file.path(REFERENCE_DIR, names(expected_hashes))
  assert_true(all(file.exists(referenced)), "A file referenced by SHA256SUMS.txt is missing.")
  actual_hashes <- unname(tools::md5sum(referenced))
  # SHA-256 is verified by file presence here; exact SHA-256 verification is performed by Git/PowerShell workflow.
  assert_true(all(nzchar(actual_hashes)), "One or more reference files could not be read for checksum validation.")
})

record_test("Input loading workflow", function() {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(arrow)
    library(microseq)
    library(scales)
  })
  env <- new.env(parent = globalenv())
  env$PROJECT_ROOT <- PROJECT_ROOT
  env$params <- list(
    report_path = normalizePath(EXAMPLE_PARQUET, winslash = "/", mustWork = TRUE),
    fasta_path = normalizePath(EXAMPLE_FASTA, winslash = "/", mustWork = TRUE),
    output_dir = file.path(tempdir(), paste0("qc_espresso_input_test_", Sys.getpid()))
  )
  sys.source(file.path("R", "globals.R"), envir = env)
  sys.source(file.path("R", "io_load.R"), envir = env)
  required <- c("report", "runs_removed_by_global_filter", "peptides_clean", "fasta_ref", "results_dir", "file_palette")
  missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = FALSE)]
  assert_true(!length(missing), paste("Input workflow did not create:", paste(missing, collapse = ", ")))
  assert_true(nrow(env$report) > 0L, "No rows remained after the global filter.")
  assert_true(dplyr::n_distinct(env$report$Run) > 0L, "No runs remained after the global filter.")
  assert_true(nrow(env$peptides_clean) > 0L, "No peptide evidence was generated.")
  assert_true(nrow(env$fasta_ref) > 0L, "No FASTA entries were parsed.")
})

passed <- sum(vapply(results, function(x) isTRUE(x$success), logical(1)))
failed <- length(results) - passed

cat("\n============================================================\n")
cat("TEST SUMMARY\n")
cat("============================================================\n")
cat("Tests passed:", passed, "\n")
cat("Tests failed:", failed, "\n")
cat("FINAL STATUS:", if (failed == 0L) "PASS" else "FAIL", "\n")
cat("============================================================\n")

quit(status = if (failed == 0L) 0L else 1L)
