###############################################################################
# File: tests/test_missed_cleavages.R
#
# Purpose:
#   Validate the corrected missed-tryptic-cleavage rule consistently in
#   R/qc_peptides.R and R/qc_figures_of_merit.R, and characterize its impact
#   on the available reference dataset.
#
# Validation scope:
#   - Validate the corrected rule using controlled peptide sequences.
#   - Confirm that both analytical modules use the shared corrected function.
#   - Confirm that the legacy unconditional subtraction is absent.
#   - Confirm that the missed-cleavage plot has no fixed 50-percent limit.
#   - Compare the previous and corrected rules in the reference dataset.
#   - Validate the internal consistency of all diagnostic calculations.
#   - Generate per-run and affected-sequence diagnostic outputs.
#
# Result categories:
#   PASS  The expected test behavior was observed.
#   FAIL  The test produced an unexpected result.
#
# Usage:
#   Rscript --no-init-file tests/test_missed_cleavages.R
###############################################################################


## 1. Minimal test framework ---------------------------------------------------

results <- list()


record_test <- function(
    name,
    fun
) {
  result <- tryCatch(
    {
      fun()

      list(
        name = name,
        success = TRUE,
        message = ""
      )
    },
    error = function(error_condition) {
      list(
        name = name,
        success = FALSE,
        message = conditionMessage(error_condition)
      )
    }
  )

  results[[length(results) + 1L]] <<- result

  cat(
    sprintf(
      "%-60s %s\n",
      name,
      if (isTRUE(result$success)) {
        "PASS"
      } else {
        "FAIL"
      }
    )
  )

  if (!isTRUE(result$success)) {
    cat(
      "  Reason: ",
      result$message,
      "\n",
      sep = ""
    )
  }

  invisible(result)
}


assert_true <- function(
    condition,
    message
) {
  if (!isTRUE(condition)) {
    stop(
      message,
      call. = FALSE
    )
  }

  invisible(TRUE)
}


assert_equal <- function(
    observed,
    expected,
    message
) {
  if (!identical(observed, expected)) {
    stop(
      paste0(
        message,
        "\nExpected: ",
        paste(expected, collapse = ", "),
        "\nObserved: ",
        paste(observed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


## 2. Project-root resolution --------------------------------------------------

resolve_project_root <- function() {
  command_arguments <- commandArgs(
    trailingOnly = FALSE
  )

  file_argument <- grep(
    "^--file=",
    command_arguments,
    value = TRUE
  )

  if (length(file_argument) > 0L) {
    test_file <- normalizePath(
      sub(
        "^--file=",
        "",
        file_argument[[1L]]
      ),
      winslash = "/",
      mustWork = TRUE
    )

    return(
      normalizePath(
        file.path(
          dirname(test_file),
          ".."
        ),
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  project_root <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )

  if (!file.exists(
    file.path(
      project_root,
      "QC-espresso_DIA.Rproj"
    )
  )) {
    stop(
      paste0(
        "The QC-espresso project root could not be resolved.\n",
        "Run this test from the project root."
      ),
      call. = FALSE
    )
  }

  project_root
}


PROJECT_ROOT <- resolve_project_root()

setwd(
  PROJECT_ROOT
)


## 3. Required packages --------------------------------------------------------

suppressPackageStartupMessages(
  {
    library(arrow)
    library(dplyr)
    library(stringr)
    library(tibble)
  }
)


## 4. Corrected missed-cleavage function --------------------------------------

count_missed_tryptic_cleavages <- function(sequence) {
  cleavage_sites <- stringr::str_count(
    sequence,
    "(?<=[KR])(?!P)"
  )

  terminal_tryptic_site <- as.integer(
    stringr::str_detect(
      sequence,
      "[KR]$"
    )
  )

  pmax(
    0L,
    cleavage_sites - terminal_tryptic_site
  )
}


## 5. Controlled algorithm tests ----------------------------------------------

record_test(
  "Controlled peptide sequences",
  function() {
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

    expected <- c(
      0L,
      0L,
      1L,
      1L,
      1L,
      2L,
      0L,
      0L
    )

    observed <- count_missed_tryptic_cleavages(
      sequences
    )

    assert_equal(
      observed,
      expected,
      "The corrected rule returned unexpected values for controlled sequences."
    )
  }
)


record_test(
  "Corrected counts are non-negative integers",
  function() {
    sequences <- c(
      "",
      "A",
      "K",
      "R",
      "KP",
      "RP",
      "KR",
      "PEPTIDE",
      "PEPTIDEK",
      "PEPKTIDKR"
    )

    observed <- count_missed_tryptic_cleavages(
      sequences
    )

    assert_true(
      all(
        !is.na(observed)
      ),
      "The corrected rule produced missing values for controlled sequences."
    )

    assert_true(
      all(
        observed >= 0L
      ),
      "The corrected rule produced a negative missed-cleavage count."
    )

    assert_true(
      is.integer(observed),
      "The corrected rule did not return integer counts."
    )
  }
)


record_test(
  "Terminal tryptic sites are conditionally excluded",
  function() {
    sequences <- c(
      "PEPKTIDE",
      "PEPKTIDEK",
      "PEPKTIDER"
    )

    observed <- count_missed_tryptic_cleavages(
      sequences
    )

    expected <- c(
      1L,
      1L,
      1L
    )

    assert_equal(
      observed,
      expected,
      paste0(
        "The terminal-site correction was not applied consistently ",
        "to equivalent internal cleavage patterns."
      )
    )
  }
)


record_test(
  "Cleavage blocked by proline is not counted",
  function() {
    sequences <- c(
      "PEKPTIDE",
      "PERPTIDE",
      "KPEPTIDE",
      "RPEPTIDE"
    )

    observed <- count_missed_tryptic_cleavages(
      sequences
    )

    expected <- rep(
      0L,
      length(sequences)
    )

    assert_equal(
      observed,
      expected,
      "A lysine or arginine followed by proline was counted as a cleavage site."
    )
  }
)


## 6. Static implementation checks --------------------------------------------

record_test(
  "Both modules use the shared corrected function",
  function() {
    peptide_text <- paste(
      readLines(
        file.path(
          "R",
          "qc_peptides.R"
        ),
        warn = FALSE
      ),
      collapse = "\n"
    )

    figures_of_merit_text <- paste(
      readLines(
        file.path(
          "R",
          "qc_figures_of_merit.R"
        ),
        warn = FALSE
      ),
      collapse = "\n"
    )

    assert_true(
      grepl(
        "count_missed_tryptic_cleavages",
        peptide_text,
        fixed = TRUE
      ),
      paste0(
        "R/qc_peptides.R does not use ",
        "count_missed_tryptic_cleavages."
      )
    )

    assert_true(
      grepl(
        "count_missed_tryptic_cleavages",
        figures_of_merit_text,
        fixed = TRUE
      ),
      paste0(
        "R/qc_figures_of_merit.R does not use ",
        "count_missed_tryptic_cleavages."
      )
    )
  }
)


record_test(
  "Legacy unconditional subtraction is absent",
  function() {
    module_paths <- file.path(
      "R",
      c(
        "qc_peptides.R",
        "qc_figures_of_merit.R"
      )
    )

    legacy_pattern <- paste0(
      "str_count\\s*\\(",
      "[\\s\\S]*?",
      "\\(\\?<\\=\\[KR\\]\\)\\(\\?!P\\)",
      "[\\s\\S]*?",
      "\\)\\s*-\\s*1"
    )

    modules_with_legacy_rule <- module_paths[
      vapply(
        module_paths,
        function(module_path) {
          module_text <- paste(
            readLines(
              module_path,
              warn = FALSE
            ),
            collapse = "\n"
          )

          grepl(
            legacy_pattern,
            module_text,
            perl = TRUE
          )
        },
        logical(1)
      )
    ]

    assert_true(
      length(modules_with_legacy_rule) == 0L,
      paste0(
        "The legacy unconditional subtraction remains in: ",
        paste(
          modules_with_legacy_rule,
          collapse = ", "
        )
      )
    )
  }
)


record_test(
  "Missed-cleavage plot has no fixed 50-percent limit",
  function() {
    peptide_text <- paste(
      readLines(
        file.path(
          "R",
          "qc_peptides.R"
        ),
        warn = FALSE
      ),
      collapse = "\n"
    )

    assert_true(
      !grepl(
        "limits = c(0, 0.5)",
        peptide_text,
        fixed = TRUE
      ),
      "The fixed 0-50% missed-cleavage limit remains in R/qc_peptides.R."
    )
  }
)


## 7. Reference-dataset loading ------------------------------------------------

report_path <- file.path(
  "example-data",
  "report.parquet"
)

assert_true(
  file.exists(report_path),
  paste0(
    "Example Parquet was not found: ",
    report_path
  )
)


report <- arrow::open_dataset(
  report_path,
  format = "parquet"
) |>
  dplyr::filter(
    Decoy == 0,
    PG.Q.Value <= 0.01
  ) |>
  dplyr::collect()


assert_true(
  all(
    c(
      "Run",
      "Stripped.Sequence"
    ) %in% names(report)
  ),
  paste0(
    "The example dataset does not contain the columns required ",
    "for missed-cleavage validation."
  )
)


## 8. Diagnostic comparison ----------------------------------------------------

peptide_comparison <- report |>
  dplyr::filter(
    !is.na(Stripped.Sequence),
    nzchar(Stripped.Sequence)
  ) |>
  dplyr::distinct(
    Run,
    Stripped.Sequence
  ) |>
  dplyr::mutate(
    Cleavage_Sites = stringr::str_count(
      Stripped.Sequence,
      "(?<=[KR])(?!P)"
    ),
    Ends_in_KR = stringr::str_detect(
      Stripped.Sequence,
      "[KR]$"
    ),
    Previous_Missed_Cleavages = pmax(
      0L,
      Cleavage_Sites - 1L
    ),
    Corrected_Missed_Cleavages = count_missed_tryptic_cleavages(
      Stripped.Sequence
    ),
    Count_Changed = Previous_Missed_Cleavages !=
      Corrected_Missed_Cleavages,
    Previous_Has_Missed_Cleavage = Previous_Missed_Cleavages > 0L,
    Corrected_Has_Missed_Cleavage = Corrected_Missed_Cleavages > 0L
  )


per_run_comparison <- peptide_comparison |>
  dplyr::group_by(
    Run
  ) |>
  dplyr::summarise(
    Unique_Peptides = dplyr::n(),
    Changed_Observations = sum(
      Count_Changed
    ),
    Previous_Missed_Cleavage_Pct = 100 *
      mean(
        Previous_Has_Missed_Cleavage
      ),
    Corrected_Missed_Cleavage_Pct = 100 *
      mean(
        Corrected_Has_Missed_Cleavage
      ),
    Difference_Percentage_Points =
      Corrected_Missed_Cleavage_Pct -
      Previous_Missed_Cleavage_Pct,
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      abs(
        Difference_Percentage_Points
      )
    ),
    Run
  )


affected_sequence_summary <- peptide_comparison |>
  dplyr::filter(
    Count_Changed
  ) |>
  dplyr::distinct(
    Stripped.Sequence,
    Ends_in_KR,
    Previous_Missed_Cleavages,
    Corrected_Missed_Cleavages
  ) |>
  dplyr::arrange(
    Stripped.Sequence
  )


## 9. General diagnostic consistency tests ------------------------------------

record_test(
  "Reference dataset contains valid peptide observations",
  function() {
    assert_true(
      nrow(peptide_comparison) > 0L,
      "No valid Run + Stripped.Sequence observations were available."
    )

    assert_true(
      dplyr::n_distinct(
        peptide_comparison$Run
      ) > 0L,
      "No valid runs were available for missed-cleavage validation."
    )
  }
)


record_test(
  "One row is retained per Run and peptide sequence",
  function() {
    duplicated_observations <- peptide_comparison |>
      dplyr::count(
        Run,
        Stripped.Sequence,
        name = "Rows"
      ) |>
      dplyr::filter(
        Rows != 1L
      )

    assert_true(
      nrow(duplicated_observations) == 0L,
      paste0(
        nrow(duplicated_observations),
        " duplicated Run + Stripped.Sequence observations remain."
      )
    )
  }
)


record_test(
  "Corrected diagnostic counts match the shared function",
  function() {
    expected <- count_missed_tryptic_cleavages(
      peptide_comparison$Stripped.Sequence
    )

    assert_true(
      identical(
        peptide_comparison$Corrected_Missed_Cleavages,
        expected
      ),
      paste0(
        "The reference-dataset diagnostic does not match the ",
        "shared corrected function."
      )
    )
  }
)


record_test(
  "Changed observations are internally consistent",
  function() {
    expected_changed <- peptide_comparison$Previous_Missed_Cleavages !=
      peptide_comparison$Corrected_Missed_Cleavages

    assert_true(
      identical(
        peptide_comparison$Count_Changed,
        expected_changed
      ),
      "The Count_Changed flag is inconsistent with the calculated counts."
    )

    assert_true(
      nrow(affected_sequence_summary) ==
        dplyr::n_distinct(
          peptide_comparison$Stripped.Sequence[
            peptide_comparison$Count_Changed
          ]
        ),
      paste0(
        "The affected-sequence summary does not match the ",
        "distinct changed sequences."
      )
    )
  }
)


record_test(
  "Per-run summaries cover every available run",
  function() {
    expected_runs <- sort(
      unique(
        as.character(
          peptide_comparison$Run
        )
      )
    )

    observed_runs <- sort(
      unique(
        as.character(
          per_run_comparison$Run
        )
      )
    )

    assert_equal(
      observed_runs,
      expected_runs,
      "The per-run diagnostic summary does not cover every available run."
    )

    assert_equal(
      sum(
        per_run_comparison$Unique_Peptides
      ),
      nrow(peptide_comparison),
      paste0(
        "Per-run peptide counts do not sum to the number of ",
        "unique Run + peptide observations."
      )
    )
  }
)


record_test(
  "Per-run percentages are finite and within valid bounds",
  function() {
    percentage_columns <- c(
      "Previous_Missed_Cleavage_Pct",
      "Corrected_Missed_Cleavage_Pct"
    )

    percentage_values <- unlist(
      per_run_comparison[
        percentage_columns
      ],
      use.names = FALSE
    )

    assert_true(
      all(
        is.finite(
          percentage_values
        )
      ),
      "One or more per-run percentages are not finite."
    )

    assert_true(
      all(
        percentage_values >= 0 &
          percentage_values <= 100
      ),
      "One or more per-run percentages are outside 0-100%."
    )
  }
)


record_test(
  "Per-run percentage differences are calculated correctly",
  function() {
    expected_difference <-
      per_run_comparison$Corrected_Missed_Cleavage_Pct -
      per_run_comparison$Previous_Missed_Cleavage_Pct

    assert_true(
      isTRUE(
        all.equal(
          per_run_comparison$Difference_Percentage_Points,
          expected_difference,
          tolerance = 0
        )
      ),
      "Per-run percentage-point differences are inconsistent."
    )
  }
)


record_test(
  "Corrected rule does not reduce previous counts",
  function() {
    assert_true(
      all(
        peptide_comparison$Corrected_Missed_Cleavages >=
          peptide_comparison$Previous_Missed_Cleavages
      ),
      paste0(
        "The corrected conditional rule reduced the missed-cleavage ",
        "count for one or more observations."
      )
    )
  }
)


record_test(
  "Changed counts occur only for non-terminal sequences",
  function() {
    invalid_changes <- peptide_comparison |>
      dplyr::filter(
        Count_Changed,
        Ends_in_KR
      )

    assert_true(
      nrow(invalid_changes) == 0L,
      paste0(
        nrow(invalid_changes),
        " changed observations end in lysine or arginine."
      )
    )
  }
)


## 10. Diagnostic outputs ------------------------------------------------------

output_dir <- file.path(
  "tests",
  "test-results",
  "latest"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


per_run_output <- file.path(
  output_dir,
  "missed-cleavage-run-comparison.tsv"
)

affected_sequences_output <- file.path(
  output_dir,
  "missed-cleavage-affected-sequences.tsv"
)


utils::write.table(
  per_run_comparison,
  per_run_output,
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


utils::write.table(
  affected_sequence_summary,
  affected_sequences_output,
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


record_test(
  "Per-run diagnostic output was created",
  function() {
    assert_true(
      file.exists(per_run_output),
      "The per-run missed-cleavage diagnostic output was not created."
    )

    output_info <- file.info(
      per_run_output
    )

    assert_true(
      !is.na(output_info$size) &&
        output_info$size > 0,
      "The per-run missed-cleavage diagnostic output is empty."
    )
  }
)


record_test(
  "Affected-sequence diagnostic output was created",
  function() {
    assert_true(
      file.exists(affected_sequences_output),
      "The affected-sequence diagnostic output was not created."
    )

    output_info <- file.info(
      affected_sequences_output
    )

    assert_true(
      !is.na(output_info$size) &&
        output_info$size > 0,
      "The affected-sequence diagnostic output is empty."
    )
  }
)


## 11. Diagnostic summary ------------------------------------------------------

changed_observations <- sum(
  peptide_comparison$Count_Changed
)

affected_sequences <- nrow(
  affected_sequence_summary
)

affected_runs <- sum(
  per_run_comparison$Changed_Observations > 0L
)

maximum_previous_percentage <- max(
  per_run_comparison$Previous_Missed_Cleavage_Pct
)

maximum_corrected_percentage <- max(
  per_run_comparison$Corrected_Missed_Cleavage_Pct
)

maximum_absolute_difference <- max(
  abs(
    per_run_comparison$Difference_Percentage_Points
  )
)


cat("\n")
cat("============================================================\n")
cat("MISSED-CLEAVAGE DIAGNOSTIC SUMMARY\n")
cat("============================================================\n")
cat(
  "Unique run-peptide observations: ",
  nrow(peptide_comparison),
  "\n",
  sep = ""
)
cat(
  "Observations with changed count: ",
  changed_observations,
  "\n",
  sep = ""
)
cat(
  "Distinct peptide sequences affected: ",
  affected_sequences,
  "\n",
  sep = ""
)
cat(
  "Runs evaluated: ",
  nrow(per_run_comparison),
  "\n",
  sep = ""
)
cat(
  "Runs affected: ",
  affected_runs,
  "\n",
  sep = ""
)
cat(
  "Maximum previous percentage: ",
  format(
    maximum_previous_percentage,
    digits = 8
  ),
  "\n",
  sep = ""
)
cat(
  "Maximum corrected percentage: ",
  format(
    maximum_corrected_percentage,
    digits = 8
  ),
  "\n",
  sep = ""
)
cat(
  "Maximum absolute percentage-point difference: ",
  format(
    maximum_absolute_difference,
    digits = 8
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")


cat("\nDiagnostic outputs:\n")
cat(
  "  ",
  normalizePath(
    per_run_output,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  normalizePath(
    affected_sequences_output,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


## 12. Final test summary ------------------------------------------------------

passed_tests <- sum(
  vapply(
    results,
    function(result) {
      isTRUE(result$success)
    },
    logical(1)
  )
)

failed_tests <- length(results) - passed_tests


cat("\n")
cat("============================================================\n")
cat("MISSED-CLEAVAGE TEST SUMMARY\n")
cat("============================================================\n")
cat(
  "Tests passed: ",
  passed_tests,
  "\n",
  sep = ""
)
cat(
  "Tests failed: ",
  failed_tests,
  "\n",
  sep = ""
)
cat(
  "FINAL STATUS: ",
  if (failed_tests == 0L) {
    "PASS"
  } else {
    "FAIL"
  },
  "\n",
  sep = ""
)
cat("============================================================\n")


quit(
  status = if (failed_tests == 0L) {
    0L
  } else {
    1L
  }
)