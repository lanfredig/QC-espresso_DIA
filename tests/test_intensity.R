###############################################################################
# File: tests/test_intensity.R
#
# Purpose:
#   Validate protein-level LFQ deduplication, precursor-intensity filtering,
#   consistency checks, and low-data behavior implemented in R/qc_intensity.R.
#
# Validation scope:
#   - Collapse repeated identical PG.MaxLFQ values per Run + Protein.Group.
#   - Preserve independent protein groups within each run.
#   - Reject multiple distinct positive PG.MaxLFQ values for the same pair.
#   - Exclude non-positive and non-finite quantitative values.
#   - Skip plots when insufficient positive values are available.
#   - Generate plots when sufficient valid values are available.
#   - Validate general PG.MaxLFQ consistency in the example DIA-NN report.
#
# Result categories:
#   PASS  The expected intensity behavior was observed.
#   FAIL  The intensity behavior was not observed.
#
# Usage:
#   Rscript --no-init-file tests/test_intensity.R
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
      "%-68s %s\n",
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
    library(dplyr)
    library(ggplot2)
    library(tibble)
  }
)


## 4. Shared test configuration ------------------------------------------------

file_palette <- c(
  A = "#0072B2",
  B = "#E69F00"
)


## 5. Isolated module execution ------------------------------------------------

run_intensity_module <- function(
    input_report
) {
  execution_environment <- new.env(
    parent = globalenv()
  )

  execution_environment$report <- input_report
  execution_environment$file_palette <- file_palette
  execution_environment$AR_DENSITY <- 1
  execution_environment$AR_BOXPLOT <- 1

  messages <- character()
  warnings <- character()
  error_condition <- NULL

  plot_file <- tempfile(
    fileext = ".pdf"
  )

  grDevices::pdf(
    plot_file
  )

  tryCatch(
    withCallingHandlers(
      sys.source(
        file.path(
          "R",
          "qc_intensity.R"
        ),
        envir = execution_environment
      ),
      message = function(message_condition) {
        messages <<- c(
          messages,
          conditionMessage(
            message_condition
          )
        )

        invokeRestart(
          "muffleMessage"
        )
      },
      warning = function(warning_condition) {
        warnings <<- c(
          warnings,
          conditionMessage(
            warning_condition
          )
        )

        invokeRestart(
          "muffleWarning"
        )
      }
    ),
    error = function(condition) {
      error_condition <<- condition
    }
  )

  try(
    grDevices::dev.off(),
    silent = TRUE
  )

  list(
    environment = execution_environment,
    messages = messages,
    warnings = warnings,
    error = error_condition,
    plot_file = plot_file
  )
}


cleanup_module_result <- function(
    module_result
) {
  if (
    !is.null(module_result$plot_file) &&
    file.exists(module_result$plot_file)
  ) {
    unlink(
      module_result$plot_file,
      force = TRUE
    )
  }

  invisible(NULL)
}


object_exists <- function(
    module_result,
    object_name
) {
  exists(
    object_name,
    envir = module_result$environment,
    inherits = FALSE
  )
}


get_object <- function(
    module_result,
    object_name
) {
  get(
    object_name,
    envir = module_result$environment,
    inherits = FALSE
  )
}


message_mentions <- function(
    module_result,
    expected_text
) {
  any(
    grepl(
      expected_text,
      module_result$messages,
      fixed = TRUE
    )
  )
}


## 6. PG.MaxLFQ deduplication --------------------------------------------------

record_test(
  "Repeated identical PG.MaxLFQ values collapse to one row",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A",
          "A"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P1",
        "P1"
      ),
      PG.MaxLFQ = c(
        100,
        100,
        100
      ),
      Precursor.Quantity = c(
        10,
        20,
        30
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    assert_true(
      object_exists(
        module_result,
        "intensity_boxplot"
      ),
      "intensity_boxplot was not created."
    )

    intensity_boxplot <- get_object(
      module_result,
      "intensity_boxplot"
    )

    assert_equal(
      nrow(intensity_boxplot),
      1L,
      "Expected one protein-level intensity row."
    )
  }
)


record_test(
  "Different proteins remain independent",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P2"
      ),
      PG.MaxLFQ = c(
        100,
        200
      ),
      Precursor.Quantity = c(
        10,
        20
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    intensity_boxplot <- get_object(
      module_result,
      "intensity_boxplot"
    )

    assert_equal(
      nrow(intensity_boxplot),
      2L,
      "Expected two independent protein-level intensity rows."
    )

    assert_true(
      setequal(
        intensity_boxplot$Protein.Group,
        c(
          "P1",
          "P2"
        )
      ),
      "The protein-level intensity inventory changed."
    )
  }
)


record_test(
  "Same protein in different runs remains independent",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "B"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P1"
      ),
      PG.MaxLFQ = c(
        100,
        200
      ),
      Precursor.Quantity = c(
        10,
        20
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    intensity_boxplot <- get_object(
      module_result,
      "intensity_boxplot"
    )

    assert_equal(
      nrow(intensity_boxplot),
      2L,
      "Expected one protein-level row for each run."
    )

    observed_runs <- sort(
      unique(
        as.character(
          intensity_boxplot$Run
        )
      )
    )

    assert_equal(
      observed_runs,
      c(
        "A",
        "B"
      ),
      "One or more runs were lost from the protein-level intensity table."
    )
  }
)


## 7. PG.MaxLFQ consistency checks --------------------------------------------

record_test(
  "Distinct positive PG.MaxLFQ values for one pair stop execution",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P1"
      ),
      PG.MaxLFQ = c(
        100,
        120
      ),
      Precursor.Quantity = c(
        10,
        20
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      !is.null(
        module_result$error
      ),
      "Expected an inconsistency error for distinct positive PG.MaxLFQ values."
    )

    assert_true(
      grepl(
        "Multiple distinct positive PG.MaxLFQ values",
        conditionMessage(
          module_result$error
        ),
        fixed = TRUE
      ),
      "The inconsistency error did not explain the PG.MaxLFQ conflict."
    )
  }
)


record_test(
  "Repeated positive value with missing duplicates is accepted",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A",
          "A"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P1",
        "P1"
      ),
      PG.MaxLFQ = c(
        100,
        NA_real_,
        100
      ),
      Precursor.Quantity = c(
        10,
        20,
        30
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    intensity_boxplot <- get_object(
      module_result,
      "intensity_boxplot"
    )

    assert_equal(
      nrow(intensity_boxplot),
      1L,
      "Missing duplicate values prevented valid PG.MaxLFQ deduplication."
    )
  }
)


record_test(
  "Non-positive duplicate values do not create inconsistency",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A",
          "A"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P1",
        "P1"
      ),
      PG.MaxLFQ = c(
        100,
        0,
        -10
      ),
      Precursor.Quantity = c(
        10,
        20,
        30
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    intensity_boxplot <- get_object(
      module_result,
      "intensity_boxplot"
    )

    assert_equal(
      nrow(intensity_boxplot),
      1L,
      "Non-positive PG.MaxLFQ values changed the valid protein-level row."
    )
  }
)


## 8. Quantitative filtering ---------------------------------------------------

record_test(
  "Non-positive values are excluded",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A",
          "B"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P2",
        "P3"
      ),
      PG.MaxLFQ = c(
        100,
        0,
        -10
      ),
      Precursor.Quantity = c(
        10,
        0,
        -5
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    intensity_boxplot <- get_object(
      module_result,
      "intensity_boxplot"
    )

    precursor_density <- get_object(
      module_result,
      "prec_density"
    )

    assert_equal(
      nrow(intensity_boxplot),
      1L,
      "Non-positive PG.MaxLFQ values were retained."
    )

    assert_equal(
      nrow(precursor_density),
      1L,
      "Non-positive precursor quantities were retained."
    )

    assert_true(
      all(
        intensity_boxplot$PG.MaxLFQ > 0
      ),
      "The protein-level intensity table contains non-positive values."
    )

    assert_true(
      all(
        precursor_density$Precursor.Quantity > 0
      ),
      "The precursor-density table contains non-positive values."
    )
  }
)


record_test(
  "Non-finite values are excluded",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A",
          "A",
          "B"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P2",
        "P3",
        "P4"
      ),
      PG.MaxLFQ = c(
        100,
        Inf,
        -Inf,
        NA_real_
      ),
      Precursor.Quantity = c(
        10,
        Inf,
        -Inf,
        NA_real_
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    intensity_boxplot <- get_object(
      module_result,
      "intensity_boxplot"
    )

    precursor_density <- get_object(
      module_result,
      "prec_density"
    )

    assert_true(
      all(
        is.finite(
          intensity_boxplot$PG.MaxLFQ
        )
      ),
      "The protein-level intensity table contains non-finite values."
    )

    assert_true(
      all(
        is.finite(
          precursor_density$Precursor.Quantity
        )
      ),
      "The precursor-density table contains non-finite values."
    )

    assert_equal(
      nrow(intensity_boxplot),
      1L,
      "Expected one finite positive PG.MaxLFQ observation."
    )

    assert_equal(
      nrow(precursor_density),
      1L,
      "Expected one finite positive precursor intensity observation."
    )
  }
)


## 9. Low-data plot behavior ---------------------------------------------------

record_test(
  "No positive precursor quantities skip density plot",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "B"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P2"
      ),
      PG.MaxLFQ = c(
        100,
        200
      ),
      Precursor.Quantity = c(
        0,
        NA_real_
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    assert_true(
      object_exists(
        module_result,
        "p_density"
      ),
      "p_density was not defined."
    )

    assert_true(
      is.null(
        get_object(
          module_result,
          "p_density"
        )
      ),
      "The raw-intensity density plot was created without positive quantities."
    )

    assert_true(
      message_mentions(
        module_result,
        "Skipping raw precursor intensity density"
      ),
      "The expected raw-intensity skip message was not emitted."
    )
  }
)


record_test(
  "No positive PG.MaxLFQ values skip protein boxplot",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "B"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P2"
      ),
      PG.MaxLFQ = c(
        0,
        NA_real_
      ),
      Precursor.Quantity = c(
        10,
        20
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    assert_true(
      object_exists(
        module_result,
        "p_box"
      ),
      "p_box was not defined."
    )

    assert_true(
      is.null(
        get_object(
          module_result,
          "p_box"
        )
      ),
      "The protein boxplot was created without positive PG.MaxLFQ values."
    )

    assert_true(
      message_mentions(
        module_result,
        "Skipping protein-level LFQ intensity boxplot"
      ),
      "The expected protein-intensity skip message was not emitted."
    )
  }
)


record_test(
  "Positive precursor quantities generate density plot",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A",
          "B",
          "B"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P2",
        "P3",
        "P4"
      ),
      PG.MaxLFQ = c(
        100,
        200,
        300,
        400
      ),
      Precursor.Quantity = c(
        10,
        20,
        30,
        40
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    density_plot <- get_object(
      module_result,
      "p_density"
    )

    assert_true(
      inherits(
        density_plot,
        "ggplot"
      ),
      "The raw precursor intensity density plot was not created."
    )
  }
)


record_test(
  "Positive PG.MaxLFQ values generate protein boxplot",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "A",
          "A",
          "B",
          "B"
        ),
        levels = c(
          "A",
          "B"
        )
      ),
      Protein.Group = c(
        "P1",
        "P2",
        "P3",
        "P4"
      ),
      PG.MaxLFQ = c(
        100,
        200,
        300,
        400
      ),
      Precursor.Quantity = c(
        10,
        20,
        30,
        40
      )
    )

    module_result <- run_intensity_module(
      report
    )

    on.exit(
      cleanup_module_result(
        module_result
      ),
      add = TRUE
    )

    assert_true(
      is.null(
        module_result$error
      ),
      if (is.null(module_result$error)) {
        ""
      } else {
        conditionMessage(
          module_result$error
        )
      }
    )

    boxplot_object <- get_object(
      module_result,
      "p_box"
    )

    assert_true(
      inherits(
        boxplot_object,
        "ggplot"
      ),
      "The protein-level LFQ boxplot was not created."
    )
  }
)


## 10. Example-dataset validation ---------------------------------------------

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


record_test(
  "Example has one positive PG.MaxLFQ value per run and group",
  function() {
    suppressPackageStartupMessages(
      library(arrow)
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

    required_columns <- c(
      "Run",
      "Protein.Group",
      "PG.MaxLFQ"
    )

    missing_columns <- setdiff(
      required_columns,
      names(
        example_report
      )
    )

    assert_true(
      length(missing_columns) == 0L,
      paste0(
        "The example report is missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )

    inconsistent_pairs <- example_report |>
      dplyr::filter(
        !is.na(PG.MaxLFQ),
        is.finite(PG.MaxLFQ),
        PG.MaxLFQ > 0
      ) |>
      dplyr::distinct(
        Run,
        Protein.Group,
        PG.MaxLFQ
      ) |>
      dplyr::count(
        Run,
        Protein.Group,
        name = "Distinct_Values"
      ) |>
      dplyr::filter(
        Distinct_Values > 1L
      )

    assert_true(
      nrow(inconsistent_pairs) == 0L,
      paste0(
        nrow(inconsistent_pairs),
        " Run + Protein.Group pairs contain multiple distinct ",
        "positive PG.MaxLFQ values."
      )
    )
  }
)


record_test(
  "Example contains valid positive intensity observations",
  function() {
    suppressPackageStartupMessages(
      library(arrow)
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

    positive_pg_values <- example_report |>
      dplyr::filter(
        !is.na(PG.MaxLFQ),
        is.finite(PG.MaxLFQ),
        PG.MaxLFQ > 0
      )

    positive_precursor_values <- example_report |>
      dplyr::filter(
        !is.na(Precursor.Quantity),
        is.finite(Precursor.Quantity),
        Precursor.Quantity > 0
      )

    assert_true(
      nrow(positive_pg_values) > 0L,
      "The example contains no positive finite PG.MaxLFQ observations."
    )

    assert_true(
      nrow(positive_precursor_values) > 0L,
      "The example contains no positive finite precursor quantities."
    )
  }
)


## 11. Final summary -----------------------------------------------------------

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
cat("INTENSITY TEST SUMMARY\n")
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