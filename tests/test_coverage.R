###############################################################################
# File: tests/test_coverage.R
# Purpose: Validate protein sequence coverage and low-data plot behavior.
###############################################################################

results <- list()

record_test <- function(name, fun) {
  result <- tryCatch(
    {
      fun()
      list(success = TRUE, message = "")
    },
    error = function(e) list(success = FALSE, message = conditionMessage(e))
  )
  results[[length(results) + 1L]] <<- result
  cat(sprintf("%-58s %s\n", name, if (result$success) "PASS" else "FAIL"))
  if (!result$success) cat("  Reason: ", result$message, "\n", sep = "")
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

resolve_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(path), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- resolve_root()
setwd(ROOT)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tibble)
})

load_module <- function(fasta_ref, peptides_clean) {
  env <- new.env(parent = globalenv())
  env$fasta_ref <- fasta_ref
  env$peptides_clean <- peptides_clean
  env$results_dir <- tempfile("coverage_results_")
  dir.create(env$results_dir, recursive = TRUE)
  env$AR_DENSITY <- 1

  messages <- character()
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    unlink(pdf_file, force = TRUE)
    unlink(env$results_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  withCallingHandlers(
    sys.source(file.path("R", "qc_coverage.R"), envir = env),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  list(env = env, messages = messages)
}

helper_env <- new.env(parent = globalenv())
helper_env$coverage_available <- FALSE
sys.source(file.path("R", "qc_coverage.R"), envir = helper_env)
calculate_coverage <- helper_env$calculate_coverage

record_test("No peptide evidence gives zero coverage", function() {
  assert_true(calculate_coverage("ABCDEFGHIJ", NULL) == 0L, "Expected zero coverage")
})

record_test("Single peptide occurrence", function() {
  assert_true(calculate_coverage("ABCDEFGHIJ", "CDE") == 3L, "Expected three covered residues")
})

record_test("Repeated peptide occurrences are all covered", function() {
  assert_true(calculate_coverage("ABCABC", "ABC") == 6L, "Expected all six residues to be covered")
})

record_test("Overlapping peptides are not double counted", function() {
  observed <- calculate_coverage("ABCDEFGH", c("ABCD", "CDEF"))
  assert_true(observed == 6L, paste("Expected six covered residues; observed", observed))
})

record_test("Absent peptide does not add coverage", function() {
  assert_true(calculate_coverage("ABCDEFGHIJ", "XYZ") == 0L, "Absent peptide changed coverage")
})

record_test("Coverage stays within zero and one hundred percent", function() {
  fasta_ref <- tibble(
    Protein.Group = c("P1", "P2"),
    Protein.Name = c("Protein 1", "Protein 2"),
    Sequence = c("ABCDEFGHIJ", "KLMNOPQRST"),
    Total_Len = c(10L, 10L)
  )
  peptides_clean <- tibble(
    Protein.Group = c("P1", "P2"),
    Peptide_List = list(c("ABCDE", "DEFGH"), c("KLM", "QRST"))
  )
  result <- load_module(fasta_ref, peptides_clean)
  values <- result$env$coverage_results$Coverage_Pct
  assert_true(all(values >= 0 & values <= 100), "Coverage is outside 0-100%")
})

record_test("Shared evidence is mapped per accession without duplicate rows", function() {
  fasta_ref <- tibble(
    Protein.Group = c("A", "B"),
    Protein.Name = c("Protein A", "Protein B"),
    Sequence = c("PEPTIDEAAAA", "PEPTIDEBBBB"),
    Total_Len = c(11L, 11L)
  )
  peptides_clean <- tibble(
    Protein.Group = c("A", "B"),
    Peptide_List = list("PEPTIDE", "PEPTIDE")
  )
  result <- load_module(fasta_ref, peptides_clean)
  coverage <- result$env$coverage_results
  assert_true(nrow(coverage) == 2L, "Expected one row per FASTA accession")
  assert_true(!anyDuplicated(coverage$Protein.Group), "Coverage rows were duplicated")
  assert_true(all(coverage$Aminoacids_Covered == 7L), "Shared peptide was not mapped independently")
})

record_test("Zero detected proteins skip density plot", function() {
  fasta_ref <- tibble(
    Protein.Group = c("P1", "P2"),
    Protein.Name = c("Protein 1", "Protein 2"),
    Sequence = c("AAAA", "BBBB"),
    Total_Len = c(4L, 4L)
  )
  peptides_clean <- tibble(
    Protein.Group = character(),
    Peptide_List = list()
  )
  result <- load_module(fasta_ref, peptides_clean)
  assert_true(nrow(result$env$coverage_results_detected) == 0L, "Expected no detected proteins")
  assert_true(is.null(result$env$p_coverage), "Density plot was created without positive coverage")
  assert_true(any(grepl("fewer than two proteins", result$messages, fixed = TRUE)), "Skip message was not emitted")
})

record_test("One detected protein skips density plot", function() {
  fasta_ref <- tibble(
    Protein.Group = c("P1", "P2"),
    Protein.Name = c("Protein 1", "Protein 2"),
    Sequence = c("AAAA", "BBBB"),
    Total_Len = c(4L, 4L)
  )
  peptides_clean <- tibble(
    Protein.Group = "P1",
    Peptide_List = list("AA")
  )
  result <- load_module(fasta_ref, peptides_clean)
  assert_true(nrow(result$env$coverage_results_detected) == 1L, "Expected one detected protein")
  assert_true(is.null(result$env$p_coverage), "Density plot was created for one observation")
  assert_true(any(grepl("fewer than two proteins", result$messages, fixed = TRUE)), "Skip message was not emitted")
})

record_test("Two detected proteins generate density plot", function() {
  fasta_ref <- tibble(
    Protein.Group = c("P1", "P2"),
    Protein.Name = c("Protein 1", "Protein 2"),
    Sequence = c("AAAA", "BBBB"),
    Total_Len = c(4L, 4L)
  )
  peptides_clean <- tibble(
    Protein.Group = c("P1", "P2"),
    Peptide_List = list("AA", "BB")
  )
  result <- load_module(fasta_ref, peptides_clean)
  assert_true(inherits(result$env$p_coverage, "ggplot"), "Density plot was not created")
})

passed <- sum(vapply(results, function(x) isTRUE(x$success), logical(1)))
failed <- length(results) - passed

cat("\n============================================================\n")
cat("COVERAGE TEST SUMMARY\n")
cat("============================================================\n")
cat("Tests passed: ", passed, "\n", sep = "")
cat("Tests failed: ", failed, "\n", sep = "")
cat("FINAL STATUS: ", if (failed == 0L) "PASS" else "FAIL", "\n", sep = "")
cat("============================================================\n")

quit(status = if (failed == 0L) 0L else 1L)
