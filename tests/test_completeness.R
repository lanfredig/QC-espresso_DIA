###############################################################################
# File: tests/test_completeness.R
#
# Purpose:
#   Validate the aggregated long-format completeness calculation implemented
#   in R/qc_completeness.R against controlled expectations and a dense
#   reference calculation.
#
# Validation scope:
#   - Validate complete and incomplete controlled datasets.
#   - Confirm that duplicate quantitative rows are collapsed.
#   - Exclude missing and blank run or precursor identifiers.
#   - Distinguish missing quantities from observed zero quantities.
#   - Validate global and per-run completeness dynamically.
#   - Compare the aggregated implementation with the dense reference
#     calculation using the available example DIA-NN report.
#   - Confirm that the production implementation does not use pivot_wider.
#
# Result categories:
#   PASS  The expected completeness behavior was observed.
#   FAIL  The completeness behavior was not observed.
#
# Usage:
#   Rscript --no-init-file tests/test_completeness.R
###############################################################################


## 1. Minimal test framework ---------------------------------------------------

results <- list()


record_test <- function(
    test_name,
    test_function
) {
  result <- tryCatch(
    {
      test_function()

      list(
        name = test_name,
        success = TRUE,
        message = ""
      )
    },
    error = function(error_condition) {
      list(
        name = test_name,
        success = FALSE,
        message = conditionMessage(
          error_condition
        )
      )
    }
  )

  results[[length(results) + 1L]] <<- result

  cat(
    sprintf(
      "%-62s %s\n",
      test_name,
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
        paste(
          expected,
          collapse = ", "
        ),
        "\nObserved: ",
        paste(
          observed,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


assert_equal_numeric <- function(
    observed,
    expected,
    message,
    tolerance = sqrt(.Machine$double.eps)
) {
  comparison <- isTRUE(
    all.equal(
      observed,
      expected,
      tolerance = tolerance,
      check.attributes = FALSE
    )
  )

  if (!comparison) {
    stop(
      paste0(
        message,
        "\nExpected: ",
        paste(
          expected,
          collapse = ", "
        ),
        "\nObserved: ",
        paste(
          observed,
          collapse = ", "
        )
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
    library(tidyr)
    library(tibble)
  }
)


## 4. Completeness calculations -----------------------------------------------

calculate_dense_reference <- function(
    input_report
) {
  required_columns <- c(
    "Run",
    "Precursor.Id",
    "Precursor.Quantity"
  )

  missing_columns <- setdiff(
    required_columns,
    names(input_report)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "The dense reference input is missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  valid_report <- input_report |>
    dplyr::filter(
      !is.na(Run),
      nzchar(
        as.character(Run)
      ),
      !is.na(Precursor.Id),
      nzchar(
        as.character(Precursor.Id)
      )
    ) |>
    dplyr::transmute(
      Run = as.character(Run),
      Precursor.Id = as.character(
        Precursor.Id
      ),
      Precursor.Quantity
    ) |>
    dplyr::group_by(
      Run,
      Precursor.Id
    ) |>
    dplyr::summarise(
      Precursor.Quantity = if (
        any(
          !is.na(Precursor.Quantity)
        )
      ) {
        dplyr::first(
          Precursor.Quantity[
            !is.na(Precursor.Quantity)
          ]
        )
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  if (nrow(valid_report) == 0L) {
    return(
      list(
        n_runs = 0L,
        n_precursors = 0L,
        total_possible = 0,
        observed = 0L,
        completeness = NA_real_,
        per_run = tibble::tibble(
          Run = character(),
          Observed = integer(),
          Total = integer(),
          Completeness_Pct = numeric()
        )
      )
    )
  }

  dense_matrix <- valid_report |>
    tidyr::pivot_wider(
      names_from = Precursor.Id,
      values_from = Precursor.Quantity
    )

  n_runs <- dplyr::n_distinct(
    dense_matrix$Run
  )

  n_precursors <- ncol(
    dense_matrix
  ) - 1L

  total_possible <- n_runs *
    n_precursors

  observed <- sum(
    !is.na(
      dense_matrix[
        setdiff(
          names(dense_matrix),
          "Run"
        )
      ]
    )
  )

  per_run <- dense_matrix |>
    dplyr::mutate(
      Observed = rowSums(
        !is.na(
          dplyr::across(
            -Run
          )
        )
      ),
      Total = n_precursors,
      Completeness_Pct = if (
        n_precursors > 0L
      ) {
        100 *
          Observed /
          Total
      } else {
        NA_real_
      }
    ) |>
    dplyr::select(
      Run,
      Observed,
      Total,
      Completeness_Pct
    ) |>
    dplyr::arrange(
      Run
    )

  list(
    n_runs = n_runs,
    n_precursors = n_precursors,
    total_possible = total_possible,
    observed = observed,
    completeness = if (
      total_possible > 0
    ) {
      100 *
        observed /
        total_possible
    } else {
      NA_real_
    },
    per_run = per_run
  )
}


calculate_aggregated <- function(
    input_report
) {
  execution_environment <- new.env(
    parent = globalenv()
  )

  execution_environment$report <- input_report

  sys.source(
    file.path(
      "R",
      "qc_completeness.R"
    ),
    envir = execution_environment
  )

  required_objects <- c(
    "completeness_global_df",
    "completeness_per_run_df",
    "data_complt"
  )

  missing_objects <- required_objects[
    !vapply(
      required_objects,
      exists,
      logical(1),
      envir = execution_environment,
      inherits = FALSE
    )
  ]

  if (length(missing_objects) > 0L) {
    stop(
      paste0(
        "qc_completeness.R did not create required objects: ",
        paste(
          missing_objects,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  global_table <- get(
    "completeness_global_df",
    envir = execution_environment,
    inherits = FALSE
  )

  per_run_table <- get(
    "completeness_per_run_df",
    envir = execution_environment,
    inherits = FALSE
  )

  completeness_data <- get(
    "data_complt",
    envir = execution_environment,
    inherits = FALSE
  )

  global_values <- setNames(
    global_table$Value,
    global_table$Metric
  )

  required_metrics <- c(
    "Number of runs",
    "Number of precursors",
    "Total possible values (Run × Precursor)",
    "Observed values",
    "Data completeness (%)"
  )

  missing_metrics <- setdiff(
    required_metrics,
    names(global_values)
  )

  if (length(missing_metrics) > 0L) {
    stop(
      paste0(
        "The global completeness table is missing metrics: ",
        paste(
          missing_metrics,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  list(
    n_runs = unname(
      global_values[["Number of runs"]]
    ),
    n_precursors = unname(
      global_values[["Number of precursors"]]
    ),
    total_possible = unname(
      global_values[["Total possible values (Run × Precursor)"]]
    ),
    observed = unname(
      global_values[["Observed values"]]
    ),
    completeness = unname(
      global_values[["Data completeness (%)"]]
    ),
    per_run = per_run_table |>
      dplyr::mutate(
        Run = as.character(Run)
      ) |>
      dplyr::arrange(
        Run
      ),
    data_complt = completeness_data
  )
}


compare_completeness_results <- function(
    dense_result,
    aggregated_result
) {
  assert_equal_numeric(
    aggregated_result$n_runs,
    dense_result$n_runs,
    "The number of runs differs between implementations.",
    tolerance = 0
  )

  assert_equal_numeric(
    aggregated_result$n_precursors,
    dense_result$n_precursors,
    "The number of precursors differs between implementations.",
    tolerance = 0
  )

  assert_equal_numeric(
    aggregated_result$total_possible,
    dense_result$total_possible,
    "The total possible count differs between implementations.",
    tolerance = 0
  )

  assert_equal_numeric(
    aggregated_result$observed,
    dense_result$observed,
    "The observed-value count differs between implementations.",
    tolerance = 0
  )

  assert_equal_numeric(
    aggregated_result$completeness,
    dense_result$completeness,
    "Global completeness differs between implementations."
  )

  assert_true(
    isTRUE(
      all.equal(
        dense_result$per_run,
        aggregated_result$per_run,
        tolerance = sqrt(
          .Machine$double.eps
        ),
        check.attributes = FALSE
      )
    ),
    "Per-run completeness differs between implementations."
  )

  invisible(TRUE)
}


## 5. Controlled complete dataset ---------------------------------------------

record_test(
  "Controlled complete matrix",
  function() {
    report <- tidyr::crossing(
      Run = c(
        "A",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P2",
        "P3"
      )
    ) |>
      dplyr::mutate(
        Precursor.Quantity = seq_len(
          dplyr::n()
        )
      )

    dense_result <- calculate_dense_reference(
      report
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    assert_equal_numeric(
      aggregated_result$observed,
      6,
      "Expected six observed values.",
      tolerance = 0
    )

    assert_equal_numeric(
      aggregated_result$completeness,
      100,
      "Expected 100% completeness."
    )

    compare_completeness_results(
      dense_result,
      aggregated_result
    )
  }
)


## 6. Controlled incomplete dataset -------------------------------------------

record_test(
  "Controlled incomplete matrix",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "A",
        "A",
        "B",
        "B",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P2",
        "P3",
        "P1",
        "P2",
        "P3"
      ),
      Precursor.Quantity = c(
        10,
        NA_real_,
        30,
        NA_real_,
        20,
        NA_real_
      )
    )

    dense_result <- calculate_dense_reference(
      report
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    assert_equal_numeric(
      aggregated_result$observed,
      3,
      "Expected three observed Run + precursor values.",
      tolerance = 0
    )

    assert_equal_numeric(
      aggregated_result$total_possible,
      6,
      "Expected six possible Run + precursor values.",
      tolerance = 0
    )

    assert_equal_numeric(
      aggregated_result$completeness,
      50,
      "Expected 50% global completeness."
    )

    compare_completeness_results(
      dense_result,
      aggregated_result
    )
  }
)


## 7. Duplicate-row behavior ---------------------------------------------------

record_test(
  "Duplicate quantitative rows are collapsed",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "A",
        "A",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P1",
        "P2",
        "P1"
      ),
      Precursor.Quantity = c(
        10,
        10,
        NA_real_,
        20
      )
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    duplicate_pairs <- aggregated_result$data_complt |>
      dplyr::count(
        Run,
        Precursor.Id,
        name = "Rows"
      ) |>
      dplyr::filter(
        Rows != 1L
      )

    assert_true(
      nrow(duplicate_pairs) == 0L,
      paste0(
        nrow(duplicate_pairs),
        " duplicated Run + Precursor.Id pairs remain."
      )
    )

    assert_equal_numeric(
      aggregated_result$observed,
      2,
      "Expected two observed Run + precursor pairs.",
      tolerance = 0
    )
  }
)


record_test(
  "Missing and observed duplicates resolve to observed",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "A",
        "B",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P1",
        "P1",
        "P1"
      ),
      Precursor.Quantity = c(
        NA_real_,
        10,
        NA_real_,
        NA_real_
      )
    )

    dense_result <- calculate_dense_reference(
      report
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    assert_equal_numeric(
      aggregated_result$observed,
      1,
      paste0(
        "A Run + precursor pair with one observed duplicate ",
        "should be counted once."
      ),
      tolerance = 0
    )

    compare_completeness_results(
      dense_result,
      aggregated_result
    )
  }
)


## 8. Identifier validation ----------------------------------------------------

record_test(
  "Missing and blank identifiers are excluded",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "A",
        "A",
        NA_character_,
        ""
      ),
      Precursor.Id = c(
        "P1",
        NA_character_,
        "",
        "P2",
        "P3"
      ),
      Precursor.Quantity = c(
        10,
        20,
        30,
        40,
        50
      )
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    assert_equal_numeric(
      aggregated_result$n_runs,
      1,
      "Expected one valid run.",
      tolerance = 0
    )

    assert_equal_numeric(
      aggregated_result$n_precursors,
      1,
      "Expected one valid precursor.",
      tolerance = 0
    )

    assert_equal_numeric(
      aggregated_result$observed,
      1,
      "Expected one observed value.",
      tolerance = 0
    )
  }
)


## 9. Quantity semantics -------------------------------------------------------

record_test(
  "Observed zero quantity contributes to completeness",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P1"
      ),
      Precursor.Quantity = c(
        0,
        NA_real_
      )
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    assert_equal_numeric(
      aggregated_result$observed,
      1,
      paste0(
        "A zero quantity is present data and should count as ",
        "one observed value."
      ),
      tolerance = 0
    )

    assert_equal_numeric(
      aggregated_result$completeness,
      50,
      "Expected 50% completeness when one of two values is observed."
    )
  }
)


record_test(
  "Negative quantity contributes to completeness",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P1"
      ),
      Precursor.Quantity = c(
        -5,
        NA_real_
      )
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    assert_equal_numeric(
      aggregated_result$observed,
      1,
      paste0(
        "Completeness measures presence rather than positivity; ",
        "a non-missing negative value should count as observed."
      ),
      tolerance = 0
    )
  }
)


## 10. Dynamic table consistency ----------------------------------------------

record_test(
  "Per-run table contains one row per valid run",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "A",
        "B",
        "B",
        "C",
        "C"
      ),
      Precursor.Id = c(
        "P1",
        "P2",
        "P1",
        "P2",
        "P1",
        "P2"
      ),
      Precursor.Quantity = c(
        10,
        20,
        30,
        NA_real_,
        NA_real_,
        NA_real_
      )
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    expected_runs <- sort(
      unique(
        report$Run
      )
    )

    observed_runs <- sort(
      unique(
        as.character(
          aggregated_result$per_run$Run
        )
      )
    )

    assert_equal(
      observed_runs,
      expected_runs,
      "The per-run completeness table does not cover every valid run."
    )

    assert_equal_numeric(
      nrow(
        aggregated_result$per_run
      ),
      length(
        expected_runs
      ),
      "Expected one completeness row per valid run.",
      tolerance = 0
    )
  }
)


record_test(
  "Per-run percentages remain within valid bounds",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "A",
        "B",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P2",
        "P1",
        "P2"
      ),
      Precursor.Quantity = c(
        10,
        NA_real_,
        20,
        30
      )
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    completeness_values <-
      aggregated_result$per_run$Completeness_Pct

    assert_true(
      all(
        is.finite(
          completeness_values
        )
      ),
      "One or more per-run completeness percentages are not finite."
    )

    assert_true(
      all(
        completeness_values >= 0 &
          completeness_values <= 100
      ),
      "One or more per-run completeness percentages are outside 0-100%."
    )
  }
)


record_test(
  "Global counts are internally consistent",
  function() {
    report <- tibble::tibble(
      Run = c(
        "A",
        "A",
        "B",
        "B"
      ),
      Precursor.Id = c(
        "P1",
        "P2",
        "P1",
        "P2"
      ),
      Precursor.Quantity = c(
        10,
        NA_real_,
        20,
        30
      )
    )

    aggregated_result <- calculate_aggregated(
      report
    )

    assert_equal_numeric(
      aggregated_result$total_possible,
      aggregated_result$n_runs *
        aggregated_result$n_precursors,
      paste0(
        "Total possible values do not equal runs multiplied ",
        "by precursors."
      ),
      tolerance = 0
    )

    expected_completeness <- 100 *
      aggregated_result$observed /
      aggregated_result$total_possible

    assert_equal_numeric(
      aggregated_result$completeness,
      expected_completeness,
      "Global completeness is inconsistent with global counts."
    )
  }
)


## 11. Example-dataset equivalence --------------------------------------------

example_report_path <- file.path(
  "example-data",
  "report.parquet"
)


record_test(
  "Example DIA-NN report is available",
  function() {
    assert_true(
      file.exists(
        example_report_path
      ),
      paste0(
        "Example DIA-NN report was not found: ",
        example_report_path
      )
    )

    report_info <- file.info(
      example_report_path
    )

    assert_true(
      !is.na(report_info$size) &&
        report_info$size > 0,
      "The example DIA-NN report is empty or inaccessible."
    )
  }
)


example_report <- arrow::open_dataset(
  example_report_path,
  format = "parquet"
) |>
  dplyr::filter(
    Decoy == 0,
    PG.Q.Value <= 0.01
  ) |>
  dplyr::collect()


record_test(
  "Example dense and aggregated calculations are equivalent",
  function() {
    dense_result <- calculate_dense_reference(
      example_report
    )

    aggregated_result <- calculate_aggregated(
      example_report
    )

    compare_completeness_results(
      dense_result,
      aggregated_result
    )
  }
)


record_test(
  "Example completeness results are internally valid",
  function() {
    aggregated_result <- calculate_aggregated(
      example_report
    )

    assert_true(
      aggregated_result$n_runs > 0,
      "The example produced no valid runs."
    )

    assert_true(
      aggregated_result$n_precursors > 0,
      "The example produced no valid precursors."
    )

    assert_true(
      aggregated_result$total_possible > 0,
      "The example produced no possible quantitative values."
    )

    assert_true(
      aggregated_result$observed >= 0 &&
        aggregated_result$observed <=
        aggregated_result$total_possible,
      "The example observed count is outside the possible range."
    )

    assert_true(
      is.finite(
        aggregated_result$completeness
      ) &&
        aggregated_result$completeness >= 0 &&
        aggregated_result$completeness <= 100,
      "The example global completeness is outside 0-100%."
    )
  }
)


## 12. Static implementation validation ---------------------------------------

record_test(
  "Production implementation contains no pivot_wider",
  function() {
    module_text <- paste(
      readLines(
        file.path(
          "R",
          "qc_completeness.R"
        ),
        warn = FALSE
      ),
      collapse = "\n"
    )

    assert_true(
      !grepl(
        "pivot_wider",
        module_text,
        fixed = TRUE
      ),
      paste0(
        "R/qc_completeness.R still creates a dense matrix ",
        "with pivot_wider."
      )
    )
  }
)


## 13. Dynamic example summary -------------------------------------------------

example_result <- calculate_aggregated(
  example_report
)


cat("\n")
cat("============================================================\n")
cat("COMPLETENESS EXAMPLE SUMMARY\n")
cat("============================================================\n")
cat(
  "Runs: ",
  example_result$n_runs,
  "\n",
  sep = ""
)
cat(
  "Precursors: ",
  example_result$n_precursors,
  "\n",
  sep = ""
)
cat(
  "Total possible values: ",
  example_result$total_possible,
  "\n",
  sep = ""
)
cat(
  "Observed values: ",
  example_result$observed,
  "\n",
  sep = ""
)
cat(
  "Global completeness (%): ",
  format(
    example_result$completeness,
    digits = 8
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")


## 14. Final summary -----------------------------------------------------------

passed_tests <- sum(
  vapply(
    results,
    function(result) {
      isTRUE(
        result$success
      )
    },
    logical(1)
  )
)

failed_tests <- length(results) -
  passed_tests


cat("\n")
cat("============================================================\n")
cat("COMPLETENESS TEST SUMMARY\n")
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