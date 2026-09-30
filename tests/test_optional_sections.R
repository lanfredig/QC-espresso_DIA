###############################################################################
# File: tests/test_optional_sections.R
#
# Purpose:
#   Characterize the current behavior of QC-espresso modules when optional
#   data are absent, incomplete, encoded as zero, or insufficient for fitting.
#
# Status categories:
#   PASS                     Expected behavior was observed.
#   OPTIONAL SECTION SKIPPED Optional data were unavailable and the module
#                            skipped the corresponding output without failing.
#   REQUIRED COLUMN ERROR    A mandatory column was absent and the module
#                            stopped, as expected for the current architecture.
#   BEHAVIOR TO CORRECT      The current behavior is not robust enough for the
#                            scientific execution mode.
#   FAIL                     An unexpected result occurred in the test itself.
#
# Usage from the project root:
#   Rscript --no-init-file tests/test_optional_sections.R
###############################################################################


## 1. Test framework -----------------------------------------------------------

results <- list()

record_result <- function(
  category,
  scenario,
  status,
  details = ""
) {
  allowed_status <- c(
    "PASS",
    "OPTIONAL SECTION SKIPPED",
    "REQUIRED COLUMN ERROR",
    "BEHAVIOR TO CORRECT",
    "FAIL"
  )

  if (!status %in% allowed_status) {
    stop(
      paste0(
        "Invalid test status: ",
        status
      ),
      call. = FALSE
    )
  }

  result <- data.frame(
    Category = category,
    Scenario = scenario,
    Status = status,
    Details = details,
    stringsAsFactors = FALSE
  )

  results[[length(results) + 1L]] <<- result

  cat(
    sprintf(
      "%-24s %-48s %s\n",
      category,
      scenario,
      status
    )
  )

  if (nzchar(details) && status != "PASS") {
    cat(
      "  ",
      details,
      "\n",
      sep = ""
    )
  }

  invisible(result)
}


run_scenario <- function(
  category,
  scenario,
  test_function
) {
  tryCatch(
    test_function(),
    error = function(error_condition) {
      record_result(
        category = category,
        scenario = scenario,
        status = "FAIL",
        details = conditionMessage(error_condition)
      )
    }
  )
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
      "Run this test from the QC-espresso project root.",
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
    library(tidyr)
    library(ggplot2)
    library(tibble)
    library(rlang)
    library(scales)
  }
)


## 4. Isolated module execution ------------------------------------------------

source_module <- function(
  module_name,
  objects
) {
  module_path <- file.path(
    PROJECT_ROOT,
    "R",
    module_name
  )

  if (!file.exists(module_path)) {
    stop(
      paste0(
        "Module not found: ",
        module_path
      ),
      call. = FALSE
    )
  }

  execution_environment <- list2env(
    objects,
    parent = globalenv()
  )

  plot_file <- tempfile(
    fileext = ".pdf"
  )

  grDevices::pdf(
    plot_file
  )

  on.exit(
    {
      try(
        grDevices::dev.off(),
        silent = TRUE
      )

      unlink(
        plot_file,
        force = TRUE
      )
    },
    add = TRUE
  )

  messages <- character(0)
  warnings <- character(0)

  error_condition <- NULL

  tryCatch(
    withCallingHandlers(
      sys.source(
        module_path,
        envir = execution_environment
      ),
      message = function(condition) {
        messages <<- c(
          messages,
          conditionMessage(condition)
        )

        invokeRestart(
          "muffleMessage"
        )
      },
      warning = function(condition) {
        warnings <<- c(
          warnings,
          conditionMessage(condition)
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

  list(
    environment = execution_environment,
    messages = messages,
    warnings = warnings,
    error = error_condition
  )
}


object_exists <- function(
  execution_result,
  object_name
) {
  exists(
    object_name,
    envir = execution_result$environment,
    inherits = FALSE
  )
}


error_mentions <- function(
  execution_result,
  expected_text
) {
  !is.null(execution_result$error) &&
    grepl(
      expected_text,
      conditionMessage(execution_result$error),
      fixed = TRUE
    )
}


message_mentions <- function(
  execution_result,
  expected_text
) {
  any(
    grepl(
      expected_text,
      execution_result$messages,
      fixed = TRUE
    )
  )
}


## 5. Shared synthetic inputs --------------------------------------------------

base_runs <- factor(
  c(
    "Run01",
    "Run02"
  ),
  levels = c(
    "Run01",
    "Run02"
  )
)

file_palette <- c(
  Run01 = "#0072B2",
  Run02 = "#E69F00"
)

legend_color <- function() {
  ggplot2::theme()
}

PANELS_PER_ROW <- 3L
AR_DENSITY <- 1


## 6. Internal trypsin standards ----------------------------------------------

run_scenario(
  "Internal standards",
  "No trypsin precursors",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Precursor.Id = c(
        "OTHER1",
        "OTHER2"
      ),
      RT = c(
        10,
        11
      ),
      Precursor.Quantity = c(
        100,
        110
      )
    )

    result <- source_module(
      "qc_chromatography_internalstd.R",
      list(
        report = report,
        legend_color = legend_color
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      message_mentions(
        result,
        "No trypsin autolysis peptides detected for RT analysis."
      ),
      "The expected RT skip message was not emitted."
    )

    assert_true(
      message_mentions(
        result,
        "No trypsin autolysis peptides detected for quantitative analysis."
      ),
      "The expected quantitative skip message was not emitted."
    )

    assert_true(
      !object_exists(
        result,
        "p_rt_std"
      ) &&
        !object_exists(
          result,
          "p_quant_std"
        ),
      "Trypsin plots were created without matching precursors."
    )

    record_result(
      "Internal standards",
      "No trypsin precursors",
      "OPTIONAL SECTION SKIPPED",
      "Both RT and quantitative panels were skipped with informative messages."
    )
  }
)


run_scenario(
  "Internal standards",
  "Partial trypsin precursor set",
  function() {
    report <- tibble::tibble(
      Run = factor(
        c(
          "Run01",
          "Run02"
        ),
        levels = levels(base_runs)
      ),
      Precursor.Id = c(
        "VATVSLPR2",
        "VATVSLPR2"
      ),
      RT = c(
        20,
        20.3
      ),
      Precursor.Quantity = c(
        1000,
        1100
      )
    )

    result <- source_module(
      "qc_chromatography_internalstd.R",
      list(
        report = report,
        legend_color = legend_color
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      object_exists(
        result,
        "p_rt_std"
      ) &&
        object_exists(
          result,
          "p_quant_std"
        ),
      "Expected plots were not created for the partial trypsin set."
    )

    record_result(
      "Internal standards",
      "Partial trypsin precursor set",
      "PASS"
    )
  }
)


run_scenario(
  "Internal standards",
  "All trypsin precursors",
  function() {
    report <- tidyr::crossing(
      Run = base_runs,
      Precursor.Id = c(
        "VATVSLPR2",
        "LSSPATLNSR2",
        "LGEHNIDVLEGNEQFINAAK3"
      )
    ) %>%
      dplyr::mutate(
        RT = rep(
          c(
            20,
            30,
            40
          ),
          times = 2
        ),
        Precursor.Quantity = seq_len(
          dplyr::n()
        ) *
          100
      )

    result <- source_module(
      "qc_chromatography_internalstd.R",
      list(
        report = report,
        legend_color = legend_color
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      object_exists(
        result,
        "p_rt_std"
      ) &&
        object_exists(
          result,
          "p_quant_std"
        ),
      "Expected plots were not created for all trypsin precursors."
    )

    record_result(
      "Internal standards",
      "All trypsin precursors",
      "PASS"
    )
  }
)


run_scenario(
  "Internal standards",
  "Mandatory column absent",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Precursor.Id = c(
        "VATVSLPR2",
        "VATVSLPR2"
      ),
      RT = c(
        20,
        21
      )
    )

    result <- source_module(
      "qc_chromatography_internalstd.R",
      list(
        report = report,
        legend_color = legend_color
      )
    )

    assert_true(
      !is.null(result$error),
      "The module accepted a report without Precursor.Quantity."
    )

    record_result(
      "Internal standards",
      "Mandatory column absent",
      "REQUIRED COLUMN ERROR",
      conditionMessage(result$error)
    )
  }
)


## 7. Ion mobility -------------------------------------------------------------

run_scenario(
  "Ion mobility",
  "No mobility column",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Precursor.Mz = c(
        400,
        500
      )
    )

    result <- source_module(
      "qc_ion_mobility.R",
      list(
        report = report,
        PANELS_PER_ROW = PANELS_PER_ROW,
        AR_DENSITY = AR_DENSITY
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      message_mentions(
        result,
        "[ION MOBILITY] Skipping: no ion mobility column found in report"
      ),
      "The expected missing-column skip message was not emitted."
    )

    record_result(
      "Ion mobility",
      "No mobility column",
      "OPTIONAL SECTION SKIPPED",
      "The complete ion-mobility section was skipped."
    )
  }
)


run_scenario(
  "Ion mobility",
  "Mobility encoded only as zero",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Precursor.Mz = c(
        400,
        500
      ),
      IM = c(
        0,
        0
      )
    )

    result <- source_module(
      "qc_ion_mobility.R",
      list(
        report = report,
        PANELS_PER_ROW = PANELS_PER_ROW,
        AR_DENSITY = AR_DENSITY
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      message_mentions(
        result,
        "no non-zero ion mobility values available for plotting"
      ),
      "The expected zero-mobility skip message was not emitted."
    )

    record_result(
      "Ion mobility",
      "Mobility encoded only as zero",
      "OPTIONAL SECTION SKIPPED",
      "The module recognized zero as unavailable ion mobility."
    )
  }
)


run_scenario(
  "Ion mobility",
  "Valid IM without predicted iIM",
  function() {
    report <- tibble::tibble(
      Run = factor(
        rep(
          levels(base_runs),
          each = 12
        ),
        levels = levels(base_runs)
      ),
      Precursor.Mz = rep(
        seq(
          400,
          950,
          length.out = 12
        ),
        times = 2
      ),
      IM = c(
        seq(
          0.8,
          1.4,
          length.out = 12
        ),
        seq(
          0.82,
          1.42,
          length.out = 12
        )
      )
    )

    result <- source_module(
      "qc_ion_mobility.R",
      list(
        report = report,
        PANELS_PER_ROW = PANELS_PER_ROW,
        AR_DENSITY = AR_DENSITY
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      object_exists(
        result,
        "p_im_density"
      ),
      "The ion-mobility density plot was not created."
    )

    assert_true(
      !object_exists(
        result,
        "p_trend"
      ),
      "A predicted trend plot was created without iIM values."
    )

    record_result(
      "Ion mobility",
      "Valid IM without predicted iIM",
      "PASS",
      "Density was generated and observed-versus-predicted trend was skipped."
    )
  }
)


run_scenario(
  "Ion mobility",
  "Valid IM with predicted iIM",
  function() {
    report <- tibble::tibble(
      Run = factor(
        rep(
          levels(base_runs),
          each = 12
        ),
        levels = levels(base_runs)
      ),
      Precursor.Mz = rep(
        seq(
          400,
          950,
          length.out = 12
        ),
        times = 2
      ),
      IM = c(
        seq(
          0.8,
          1.4,
          length.out = 12
        ),
        seq(
          0.82,
          1.42,
          length.out = 12
        )
      ),
      iIM = c(
        seq(
          0.81,
          1.39,
          length.out = 12
        ),
        seq(
          0.83,
          1.41,
          length.out = 12
        )
      )
    )

    result <- source_module(
      "qc_ion_mobility.R",
      list(
        report = report,
        PANELS_PER_ROW = PANELS_PER_ROW,
        AR_DENSITY = AR_DENSITY
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      object_exists(
        result,
        "p_im_density"
      ) &&
        object_exists(
          result,
          "p_trend"
        ),
      "Both ion-mobility plots were not created."
    )

    record_result(
      "Ion mobility",
      "Valid IM with predicted iIM",
      "PASS"
    )
  }
)


## 8. Chromatography -----------------------------------------------------------

run_scenario(
  "Chromatography",
  "Valid required columns",
  function() {
    report <- tibble::tibble(
      Run = factor(
        rep(
          levels(base_runs),
          each = 10
        ),
        levels = levels(base_runs)
      ),
      RT = rep(
        seq(
          1,
          10,
          length.out = 10
        ),
        times = 2
      ),
      FWHM = rep(
        seq(
          0.05,
          0.1,
          length.out = 10
        ),
        times = 2
      ),
      RT.Start = rep(
        seq(
          0.9,
          9.9,
          length.out = 10
        ),
        times = 2
      ),
      RT.Stop = rep(
        seq(
          1.1,
          10.1,
          length.out = 10
        ),
        times = 2
      )
    )

    result <- source_module(
      "qc_chromatography.R",
      list(
        report = report,
        file_palette = file_palette
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      object_exists(result, "p_rt") &&
        object_exists(result, "p_fwhm") &&
        object_exists(result, "p_base"),
      "One or more chromatography plots were not created."
    )

    record_result(
      "Chromatography",
      "Valid required columns",
      "PASS"
    )
  }
)


run_scenario(
  "Chromatography",
  "Required FWHM column absent",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      RT = c(1, 2),
      RT.Start = c(0.9, 1.9),
      RT.Stop = c(1.1, 2.1)
    )

    result <- source_module(
      "qc_chromatography.R",
      list(
        report = report,
        file_palette = file_palette
      )
    )

    assert_true(
      !is.null(result$error),
      "The module accepted a report without FWHM."
    )

    record_result(
      "Chromatography",
      "Required FWHM column absent",
      "REQUIRED COLUMN ERROR",
      conditionMessage(result$error)
    )
  }
)


run_scenario(
  "Chromatography",
  "Required RT.Start column absent",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      RT = c(1, 2),
      FWHM = c(0.1, 0.1),
      RT.Stop = c(1.1, 2.1)
    )

    result <- source_module(
      "qc_chromatography.R",
      list(
        report = report,
        file_palette = file_palette
      )
    )

    assert_true(
      !is.null(result$error),
      "The module accepted a report without RT.Start."
    )

    record_result(
      "Chromatography",
      "Required RT.Start column absent",
      "REQUIRED COLUMN ERROR",
      conditionMessage(result$error)
    )
  }
)


run_scenario(
  "Chromatography",
  "Required RT.Stop column absent",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      RT = c(1, 2),
      FWHM = c(0.1, 0.1),
      RT.Start = c(0.9, 1.9)
    )

    result <- source_module(
      "qc_chromatography.R",
      list(
        report = report,
        file_palette = file_palette
      )
    )

    assert_true(
      !is.null(result$error),
      "The module accepted a report without RT.Stop."
    )

    record_result(
      "Chromatography",
      "Required RT.Stop column absent",
      "REQUIRED COLUMN ERROR",
      conditionMessage(result$error)
    )
  }
)


run_scenario(
  "Chromatography",
  "All chromatography values missing",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      RT = c(NA_real_, NA_real_),
      FWHM = c(NA_real_, NA_real_),
      RT.Start = c(NA_real_, NA_real_),
      RT.Stop = c(NA_real_, NA_real_)
    )

    result <- source_module(
      "qc_chromatography.R",
      list(
        report = report,
        file_palette = file_palette
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    expected_messages <- c(
      "Skipping retention-time identification density",
      "Skipping FWHM distribution",
      "Skipping peak base width distribution"
    )

    missing_messages <- expected_messages[
      !vapply(
        expected_messages,
        message_mentions,
        logical(1),
        execution_result = result
      )
    ]

    assert_true(
      length(missing_messages) == 0L,
      paste0(
        "Missing controlled skip message(s): ",
        paste(missing_messages, collapse = ", ")
      )
    )

    assert_true(
      !object_exists(result, "p_rt") &&
        !object_exists(result, "p_fwhm") &&
        !object_exists(result, "p_base"),
      "One or more chromatography plots were created without valid values."
    )

    record_result(
      "Chromatography",
      "All chromatography values missing",
      "OPTIONAL SECTION SKIPPED",
      "All three chromatography subsections were skipped with informative messages."
    )
  }
)


## 9. Mass error ---------------------------------------------------------------

run_scenario(
  "Mass error",
  "Valid required columns",
  function() {
    report <- tibble::tibble(
      Run = factor(
        rep(
          levels(base_runs),
          each = 20
        ),
        levels = levels(base_runs)
      ),
      Precursor.Mz = rep(
        seq(
          400,
          1000,
          length.out = 20
        ),
        times = 2
      ),
      Precursor.Quantity = rep(
        1000,
        40
      ),
      Ms1.Apex.Mz.Delta = rep(
        seq(
          -0.001,
          0.001,
          length.out = 20
        ),
        times = 2
      )
    )

    result <- source_module(
      "qc_mass_error.R",
      list(
        report = report
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      object_exists(result, "p_mass") &&
        object_exists(result, "p_bin"),
      "One or more mass-error plots were not created."
    )

    record_result(
      "Mass error",
      "Valid required columns",
      "PASS"
    )
  }
)


run_scenario(
  "Mass error",
  "No valid positive quantities",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Precursor.Mz = c(500, 600),
      Precursor.Quantity = c(0, NA_real_),
      Ms1.Apex.Mz.Delta = c(0.001, 0.001)
    )

    result <- source_module(
      "qc_mass_error.R",
      list(
        report = report
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      message_mentions(
        result,
        "[MASS ERROR] Skipping: no valid positive-quantity mass-error observations available"
      ),
      "The expected mass-error skip message was not emitted."
    )

    assert_true(
      !object_exists(result, "p_mass") &&
        !object_exists(result, "p_bin"),
      "Mass-error plots were created without valid observations."
    )

    record_result(
      "Mass error",
      "No valid positive quantities",
      "OPTIONAL SECTION SKIPPED",
      "Both mass-error plots were skipped with an informative message."
    )
  }
)


run_scenario(
  "Mass error",
  "Required delta column absent",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Precursor.Mz = c(500, 600),
      Precursor.Quantity = c(100, 100)
    )

    result <- source_module(
      "qc_mass_error.R",
      list(
        report = report
      )
    )

    assert_true(
      !is.null(result$error),
      "The module accepted a report without Ms1.Apex.Mz.Delta."
    )

    record_result(
      "Mass error",
      "Required delta column absent",
      "REQUIRED COLUMN ERROR",
      conditionMessage(result$error)
    )
  }
)


## 10. DIA-NN quality metrics --------------------------------------------------

run_scenario(
  "DIA-NN quality",
  "Valid quality columns",
  function() {
    report <- tibble::tibble(
      Run = factor(
        rep(
          levels(base_runs),
          each = 12
        ),
        levels = levels(base_runs)
      ),
      Quantity.Quality = rep(
        seq(0.1, 0.95, length.out = 12),
        times = 2
      ),
      Empirical.Quality = rep(
        seq(0.15, 0.98, length.out = 12),
        times = 2
      ),
      Protein.Group = rep(
        paste0("P", seq_len(6)),
        times = 4
      ),
      Precursor.Id = paste0("PREC", seq_len(24)),
      PG.MaxLFQ.Quality = rep(
        seq(0.4, 0.95, length.out = 12),
        times = 2
      )
    )

    result <- source_module(
      "qc_diann_quality_metrics.R",
      list(
        report = report,
        PANELS_PER_ROW = PANELS_PER_ROW,
        file_palette = file_palette,
        AR_DENSITY = AR_DENSITY
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    assert_true(
      object_exists(result, "p_joint") &&
        object_exists(result, "p_pg"),
      "One or more DIA-NN quality plots were not created."
    )

    record_result(
      "DIA-NN quality",
      "Valid quality columns",
      "PASS"
    )
  }
)


run_scenario(
  "DIA-NN quality",
  "All quality values missing",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Quantity.Quality = c(NA_real_, NA_real_),
      Empirical.Quality = c(NA_real_, NA_real_),
      Protein.Group = c("P1", "P2"),
      Precursor.Id = c("PREC1", "PREC2"),
      PG.MaxLFQ.Quality = c(NA_real_, NA_real_)
    )

    result <- source_module(
      "qc_diann_quality_metrics.R",
      list(
        report = report,
        PANELS_PER_ROW = PANELS_PER_ROW,
        file_palette = file_palette,
        AR_DENSITY = AR_DENSITY
      )
    )

    assert_true(
      is.null(result$error),
      conditionMessage(result$error)
    )

    expected_messages <- c(
      "Skipping joint precursor reliability space",
      "Skipping protein quantification stability"
    )

    missing_messages <- expected_messages[
      !vapply(
        expected_messages,
        message_mentions,
        logical(1),
        execution_result = result
      )
    ]

    assert_true(
      length(missing_messages) == 0L,
      paste0(
        "Missing controlled skip message(s): ",
        paste(missing_messages, collapse = ", ")
      )
    )

    assert_true(
      !object_exists(result, "p_joint") &&
        !object_exists(result, "p_pg"),
      "DIA-NN quality plots were created without valid quality values."
    )

    record_result(
      "DIA-NN quality",
      "All quality values missing",
      "OPTIONAL SECTION SKIPPED",
      "Both DIA-NN quality subsections were skipped with informative messages."
    )
  }
)


run_scenario(
  "DIA-NN quality",
  "Required quality column absent",
  function() {
    report <- tibble::tibble(
      Run = base_runs,
      Quantity.Quality = c(0.5, 0.6),
      Protein.Group = c("P1", "P2"),
      Precursor.Id = c("PREC1", "PREC2"),
      PG.MaxLFQ.Quality = c(0.7, 0.8)
    )

    result <- source_module(
      "qc_diann_quality_metrics.R",
      list(
        report = report,
        PANELS_PER_ROW = PANELS_PER_ROW,
        file_palette = file_palette,
        AR_DENSITY = AR_DENSITY
      )
    )

    assert_true(
      !is.null(result$error),
      "The module accepted a report without Empirical.Quality."
    )

    record_result(
      "DIA-NN quality",
      "Required quality column absent",
      "REQUIRED COLUMN ERROR",
      conditionMessage(result$error)
    )
  }
)


## 11. Consolidated outputs ----------------------------------------------------

result_table <- do.call(
  rbind,
  results
)

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

utils::write.table(
  result_table,
  file.path(
    output_dir,
    "optional-section-results.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

status_levels <- c(
  "PASS",
  "OPTIONAL SECTION SKIPPED",
  "REQUIRED COLUMN ERROR",
  "BEHAVIOR TO CORRECT",
  "FAIL"
)

status_counts <- table(
  factor(
    result_table$Status,
    levels = status_levels
  )
)

unexpected_failures <- unname(
  status_counts[
    "FAIL"
  ]
)

behavior_gaps <- unname(
  status_counts[
    "BEHAVIOR TO CORRECT"
  ]
)

cat("\n")
cat("============================================================\n")
cat("OPTIONAL-SECTION TEST SUMMARY\n")
cat("============================================================\n")

for (status_name in status_levels) {
  cat(
    status_name,
    ": ",
    unname(
      status_counts[
        status_name
      ]
    ),
    "\n",
    sep = ""
  )
}

cat("\n")
cat(
  "FINAL STATUS: ",
  if (unexpected_failures == 0L) {
    if (behavior_gaps == 0L) {
      "PASS"
    } else {
      "PASS WITH DOCUMENTED GAPS"
    }
  } else {
    "FAIL"
  },
  "\n",
  sep = ""
)
cat("============================================================\n")
cat(
  "Detailed results: ",
  normalizePath(
    file.path(
      output_dir,
      "optional-section-results.tsv"
    ),
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)

quit(
  status = if (unexpected_failures == 0L) {
    0L
  } else {
    1L
  }
)
