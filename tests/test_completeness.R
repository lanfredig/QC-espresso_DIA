###############################################################################
# File: tests/test_completeness.R
# Purpose:
#   Validate equivalence between the previous dense completeness calculation
#   and the new aggregated long-format implementation.
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
  cat(sprintf("%-56s %s\n", name, if (result$success) "PASS" else "FAIL"))
  if (!result$success) cat("  Reason: ", result$message, "\n", sep = "")
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

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

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

calculate_dense <- function(input_report) {
  matrix <- input_report %>%
    dplyr::select(Run, Precursor.Id, Precursor.Quantity) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = Precursor.Id,
      values_from = Precursor.Quantity
    )

  n_runs <- dplyr::n_distinct(matrix$Run)
  n_precursors <- ncol(matrix) - 1L
  total_possible <- n_runs * n_precursors
  observed <- sum(!is.na(matrix[-1]))

  per_run <- matrix %>%
    dplyr::mutate(
      Observed = rowSums(!is.na(dplyr::across(-Run))),
      Total = n_precursors,
      Completeness_Pct = 100 * Observed / Total
    ) %>%
    dplyr::select(Run, Observed, Total, Completeness_Pct) %>%
    dplyr::arrange(Run)

  list(
    n_runs = n_runs,
    n_precursors = n_precursors,
    total_possible = total_possible,
    observed = observed,
    completeness = 100 * observed / total_possible,
    per_run = per_run
  )
}

calculate_aggregated <- function(input_report) {
  env <- new.env(parent = globalenv())
  env$report <- input_report
  sys.source(file.path("R", "qc_completeness.R"), envir = env)

  values <- setNames(
    env$completeness_global_df$Value,
    env$completeness_global_df$Metric
  )

  list(
    n_runs = unname(values[["Number of runs"]]),
    n_precursors = unname(values[["Number of precursors"]]),
    total_possible = unname(values[["Total possible values (Run × Precursor)"]]),
    observed = unname(values[["Observed values"]]),
    completeness = unname(values[["Data completeness (%)"]]),
    per_run = env$completeness_per_run_df %>% dplyr::arrange(Run),
    data_complt = env$data_complt
  )
}

record_test("Controlled complete matrix", function() {
  report <- tidyr::crossing(
    Run = c("A", "B"),
    Precursor.Id = c("P1", "P2", "P3")
  ) %>%
    dplyr::mutate(Precursor.Quantity = seq_len(dplyr::n()))

  dense <- calculate_dense(report)
  aggregated <- calculate_aggregated(report)

  assert_true(aggregated$observed == 6, "Expected six observed values")
  assert_true(aggregated$completeness == 100, "Expected 100% completeness")
  assert_true(isTRUE(all.equal(dense$per_run, aggregated$per_run)), "Per-run results differ")
})

record_test("Controlled incomplete matrix", function() {
  report <- tibble::tibble(
    Run = c("A", "A", "A", "B", "B", "B"),
    Precursor.Id = c("P1", "P2", "P3", "P1", "P2", "P3"),
    Precursor.Quantity = c(10, NA, 30, NA, 20, NA)
  )

  dense <- calculate_dense(report)
  aggregated <- calculate_aggregated(report)

  assert_true(dense$observed == aggregated$observed, "Observed counts differ")
  assert_true(isTRUE(all.equal(dense$completeness, aggregated$completeness)), "Global completeness differs")
  assert_true(isTRUE(all.equal(dense$per_run, aggregated$per_run)), "Per-run results differ")
})

record_test("Duplicate quantitative rows are collapsed", function() {
  report <- tibble::tibble(
    Run = c("A", "A", "A", "B"),
    Precursor.Id = c("P1", "P1", "P2", "P1"),
    Precursor.Quantity = c(10, 10, NA, 20)
  )

  aggregated <- calculate_aggregated(report)
  duplicates <- aggregated$data_complt %>%
    dplyr::count(Run, Precursor.Id) %>%
    dplyr::filter(n != 1L)

  assert_true(nrow(duplicates) == 0L, "Duplicate Run + Precursor.Id rows remain")
  assert_true(aggregated$observed == 2, "Expected two observed pairs")
})

record_test("Missing and blank identifiers are excluded", function() {
  report <- tibble::tibble(
    Run = c("A", "A", "A", NA_character_),
    Precursor.Id = c("P1", NA_character_, "", "P2"),
    Precursor.Quantity = c(10, 20, 30, 40)
  )

  aggregated <- calculate_aggregated(report)
  assert_true(aggregated$n_runs == 1, "Expected one valid run")
  assert_true(aggregated$n_precursors == 1, "Expected one valid precursor")
  assert_true(aggregated$observed == 1, "Expected one observed value")
})

report <- arrow::open_dataset(
  file.path("example-data", "dataset_4_contaminans.parquet"),
  format = "parquet"
) %>%
  dplyr::filter(Decoy == 0, PG.Q.Value <= 0.01) %>%
  dplyr::collect()

record_test("Article example dense and aggregated equivalence", function() {
  dense <- calculate_dense(report)
  aggregated <- calculate_aggregated(report)

  assert_true(dense$n_runs == aggregated$n_runs, "Run counts differ")
  assert_true(dense$n_precursors == aggregated$n_precursors, "Precursor counts differ")
  assert_true(dense$total_possible == aggregated$total_possible, "Total possible counts differ")
  assert_true(dense$observed == aggregated$observed, "Observed counts differ")
  assert_true(isTRUE(all.equal(dense$completeness, aggregated$completeness, tolerance = 0)), "Global completeness differs")
  assert_true(isTRUE(all.equal(dense$per_run, aggregated$per_run, tolerance = 0)), "Per-run completeness differs")
})

record_test("Implementation contains no pivot_wider", function() {
  text <- paste(readLines(file.path("R", "qc_completeness.R"), warn = FALSE), collapse = "\n")
  assert_true(!grepl("pivot_wider", text, fixed = TRUE), "qc_completeness.R still creates a dense matrix")
})

passed <- sum(vapply(results, function(x) isTRUE(x$success), logical(1)))
failed <- length(results) - passed

cat("\n============================================================\n")
cat("COMPLETENESS TEST SUMMARY\n")
cat("============================================================\n")
cat("Tests passed: ", passed, "\n", sep = "")
cat("Tests failed: ", failed, "\n", sep = "")
cat("FINAL STATUS: ", if (failed == 0L) "PASS" else "FAIL", "\n", sep = "")
cat("============================================================\n")

quit(status = if (failed == 0L) 0L else 1L)
