###############################################################################
# File: tests/test_missed_cleavages.R
# Purpose: Validate the corrected missed-tryptic-cleavage rule consistently in
# qc_peptides.R and qc_figures_of_merit.R.
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
  library(stringr)
})

count_missed_tryptic_cleavages <- function(sequence) {
  cleavage_sites <- stringr::str_count(sequence, "(?<=[KR])(?!P)")
  terminal_tryptic_site <- as.integer(stringr::str_detect(sequence, "[KR]$"))
  pmax(0L, cleavage_sites - terminal_tryptic_site)
}

record_test("Controlled peptide sequences", function() {
  sequences <- c(
    "PEPTIDEK",
    "PEPTIDER",
    "PEPKTIDER",
    "PEPKTIDE",
    "PEPRTIDE",
    "PEPKTIDKR",
    "PEKPTIDE",
    "PERPTIDE"
  )
  expected <- c(0L, 0L, 1L, 1L, 1L, 2L, 0L, 0L)
  observed <- count_missed_tryptic_cleavages(sequences)
  assert_true(
    identical(observed, expected),
    paste0(
      "Expected ", paste(expected, collapse = ", "),
      "; observed ", paste(observed, collapse = ", ")
    )
  )
})

record_test("Both modules use the shared corrected function", function() {
  peptide_text <- paste(readLines("R/qc_peptides.R", warn = FALSE), collapse = "\n")
  fom_text <- paste(readLines("R/qc_figures_of_merit.R", warn = FALSE), collapse = "\n")
  assert_true(grepl("count_missed_tryptic_cleavages", peptide_text, fixed = TRUE), "qc_peptides.R does not use the corrected function")
  assert_true(grepl("count_missed_tryptic_cleavages", fom_text, fixed = TRUE), "qc_figures_of_merit.R does not use the corrected function")
  assert_true(!grepl("str_count(\n        Stripped.Sequence,\n        \"(?<=[KR])(?!P)\"\n      ) - 1", peptide_text, fixed = TRUE), "Legacy subtraction remains in qc_peptides.R")
  assert_true(!grepl("str_count(\n          Stripped.Sequence,\n          \"(?<=[KR])(?!P)\"\n        ) - 1", fom_text, fixed = TRUE), "Legacy subtraction remains in qc_figures_of_merit.R")
})

record_test("Missed-cleavage plot has no fixed 50-percent limit", function() {
  text <- paste(readLines("R/qc_peptides.R", warn = FALSE), collapse = "\n")
  assert_true(!grepl("limits = c(0, 0.5)", text, fixed = TRUE), "Fixed 0-50% limit remains in qc_peptides.R")
})

report <- arrow::open_dataset(
  file.path("example-data", "dataset_4_contaminans.parquet"),
  format = "parquet"
) |>
  dplyr::filter(Decoy == 0, PG.Q.Value <= 0.01) |>
  dplyr::collect()

comparison <- report |>
  dplyr::filter(!is.na(Stripped.Sequence), nzchar(Stripped.Sequence)) |>
  dplyr::distinct(Run, Stripped.Sequence) |>
  dplyr::mutate(
    cleavage_sites = stringr::str_count(Stripped.Sequence, "(?<=[KR])(?!P)"),
    current = pmax(0L, cleavage_sites - 1L),
    corrected = count_missed_tryptic_cleavages(Stripped.Sequence),
    changed = current != corrected
  )

record_test("Article example impact matches diagnostic baseline", function() {
  assert_true(nrow(comparison) == 25803L, paste("Unexpected run-peptide count:", nrow(comparison)))
  assert_true(sum(comparison$changed) == 79L, paste("Unexpected changed-observation count:", sum(comparison$changed)))
  assert_true(dplyr::n_distinct(comparison$Stripped.Sequence[comparison$changed]) == 48L, "Unexpected affected-sequence count")
  affected_runs <- comparison |>
    dplyr::group_by(Run) |>
    dplyr::summarise(changed = any(changed), .groups = "drop") |>
    dplyr::filter(changed) |>
    nrow()
  assert_true(affected_runs == 14L, paste("Unexpected affected-run count:", affected_runs))
})

per_run <- comparison |>
  dplyr::group_by(Run) |>
  dplyr::summarise(
    corrected_pct = 100 * mean(corrected > 0L),
    .groups = "drop"
  )

record_test("Corrected article percentages match expected range", function() {
  maximum <- max(per_run$corrected_pct)
  assert_true(abs(maximum - 28.17365) < 1e-4, paste("Unexpected maximum corrected percentage:", maximum))
  assert_true(all(per_run$corrected_pct <= 50), "The article example contains a corrected value above 50%")
})

passed <- sum(vapply(results, function(x) isTRUE(x$success), logical(1)))
failed <- length(results) - passed

cat("\n============================================================\n")
cat("MISSED-CLEAVAGE TEST SUMMARY\n")
cat("============================================================\n")
cat("Tests passed: ", passed, "\n", sep = "")
cat("Tests failed: ", failed, "\n", sep = "")
cat("FINAL STATUS: ", if (failed == 0L) "PASS" else "FAIL", "\n", sep = "")
cat("============================================================\n")

quit(status = if (failed == 0L) 0L else 1L)
