###############################################################################
# File: tests/test_runner_report.R
#
# Purpose:
#   Validate the QC-espresso scientific runner and Quarto report through two
#   isolated end-to-end executions using the available example inputs.
#
# Validation scope:
#   - Confirm that the runner and scientific QMD are available and readable.
#   - Validate required QMD parameters and project-version integration.
#   - Confirm that obsolete parameters and global setwd calls are absent.
#   - Resolve the Quarto executable and read its version.
#   - Execute the scientific runner twice under a path containing spaces.
#   - Confirm that each execution uses an isolated output directory.
#   - Validate result-object structure and execution metadata.
#   - Validate the rendered HTML report and all declared generated files.
#   - Validate recognized scientific CSV outputs when they are generated.
#
# Result categories:
#   PASS  The expected runner or report behavior was observed.
#   FAIL  The runner or report behavior was not observed.
#
# Usage:
#   Rscript --no-init-file tests/test_runner_report.R
###############################################################################


## 1. Project-root resolution --------------------------------------------------

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


## 2. Minimal test framework ---------------------------------------------------

results <- list()


record_test <- function(
    test_name,
    test_function
) {
  started_at <- Sys.time()

  outcome <- tryCatch(
    {
      details <- test_function()

      list(
        success = TRUE,
        details = if (
          is.character(details) &&
          length(details) > 0L
        ) {
          paste(
            details,
            collapse = " "
          )
        } else {
          ""
        },
        elapsed_seconds = as.numeric(
          difftime(
            Sys.time(),
            started_at,
            units = "secs"
          )
        )
      )
    },
    error = function(error_condition) {
      list(
        success = FALSE,
        details = conditionMessage(
          error_condition
        ),
        elapsed_seconds = as.numeric(
          difftime(
            Sys.time(),
            started_at,
            units = "secs"
          )
        )
      )
    }
  )

  results[[length(results) + 1L]] <<- data.frame(
    Test = test_name,
    Success = outcome$success,
    Details = outcome$details,
    Elapsed_Seconds = outcome$elapsed_seconds,
    stringsAsFactors = FALSE
  )

  cat(
    sprintf(
      "%-60s %s\n",
      test_name,
      if (outcome$success) {
        "PASS"
      } else {
        "FAIL"
      }
    )
  )

  if (!outcome$success) {
    cat(
      "  Reason: ",
      outcome$details,
      "\n",
      sep = ""
    )
  }

  invisible(outcome)
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


## 3. Paths and shared objects -------------------------------------------------

runner_path <- file.path(
  PROJECT_ROOT,
  "run_qc.R"
)

qmd_path <- file.path(
  PROJECT_ROOT,
  "report",
  "qc_espresso_scientific.qmd"
)

report_path <- file.path(
  PROJECT_ROOT,
  "example-data",
  "report.parquet"
)

fasta_path <- file.path(
  PROJECT_ROOT,
  "example-data",
  "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta"
)


runner_environment <- new.env(
  parent = globalenv()
)

runner_result_1 <- NULL
runner_result_2 <- NULL


## 4. Static runner and QMD validation -----------------------------------------

record_test(
  "Scientific runner exists",
  function() {
    assert_true(
      file.exists(runner_path),
      paste0(
        "Scientific runner was not found: ",
        runner_path
      )
    )

    runner_info <- file.info(
      runner_path
    )

    assert_true(
      !is.na(runner_info$size) &&
        runner_info$size > 0,
      "The scientific runner is empty or inaccessible."
    )
  }
)


record_test(
  "Scientific QMD exists",
  function() {
    assert_true(
      file.exists(qmd_path),
      paste0(
        "Scientific QMD was not found: ",
        qmd_path
      )
    )

    qmd_info <- file.info(
      qmd_path
    )

    assert_true(
      !is.na(qmd_info$size) &&
        qmd_info$size > 0,
      "The scientific QMD is empty or inaccessible."
    )
  }
)


record_test(
  "Runner syntax and QMD content are readable",
  function() {
    parse(
      file = runner_path
    )

    qmd_lines <- readLines(
      qmd_path,
      warn = FALSE
    )

    assert_true(
      length(qmd_lines) > 0L &&
        any(
          nzchar(
            trimws(qmd_lines)
          )
        ),
      "The scientific QMD contains no readable content."
    )
  }
)


record_test(
  "Scientific QMD declares required parameters",
  function() {
    qmd_text <- paste(
      readLines(
        qmd_path,
        warn = FALSE
      ),
      collapse = "\n"
    )

    required_parameters <- c(
      "report_path:",
      "fasta_path:",
      "output_dir:"
    )

    missing_parameters <- required_parameters[
      !vapply(
        required_parameters,
        grepl,
        logical(1),
        x = qmd_text,
        fixed = TRUE
      )
    ]

    assert_true(
      length(missing_parameters) == 0L,
      paste0(
        "The scientific QMD is missing required parameters: ",
        paste(
          missing_parameters,
          collapse = ", "
        )
      )
    )
  }
)


record_test(
  "Scientific QMD excludes obsolete quarto_exec parameter",
  function() {
    qmd_text <- paste(
      readLines(
        qmd_path,
        warn = FALSE
      ),
      collapse = "\n"
    )

    assert_true(
      !grepl(
        "quarto_exec:",
        qmd_text,
        fixed = TRUE
      ),
      "The obsolete quarto_exec parameter remains in the scientific QMD."
    )
  }
)


record_test(
  "Scientific QMD reads project version",
  function() {
    qmd_text <- paste(
      readLines(
        qmd_path,
        warn = FALSE
      ),
      collapse = "\n"
    )

    assert_true(
      grepl(
        "VERSION",
        qmd_text,
        fixed = TRUE
      ),
      "The scientific QMD does not read the project VERSION file."
    )
  }
)


record_test(
  "Scientific QMD avoids global setwd",
  function() {
    qmd_text <- paste(
      readLines(
        qmd_path,
        warn = FALSE
      ),
      collapse = "\n"
    )

    assert_true(
      !grepl(
        "setwd(",
        qmd_text,
        fixed = TRUE
      ),
      "The scientific QMD changes the global working directory."
    )
  }
)


## 5. Example-input validation -------------------------------------------------

record_test(
  "Example DIA-NN report is available",
  function() {
    assert_true(
      file.exists(report_path),
      paste0(
        "Example DIA-NN report was not found: ",
        report_path
      )
    )

    report_info <- file.info(
      report_path
    )

    assert_true(
      !is.na(report_info$size) &&
        report_info$size > 0,
      "The example DIA-NN report is empty or inaccessible."
    )
  }
)


record_test(
  "Example FASTA reference is available",
  function() {
    assert_true(
      file.exists(fasta_path),
      paste0(
        "Example FASTA reference was not found: ",
        fasta_path
      )
    )

    fasta_info <- file.info(
      fasta_path
    )

    assert_true(
      !is.na(fasta_info$size) &&
        fasta_info$size > 0,
      "The example FASTA reference is empty or inaccessible."
    )
  }
)


## 6. Runner loading and Quarto discovery -------------------------------------

record_test(
  "Scientific runner loads required functions",
  function() {
    sys.source(
      runner_path,
      envir = runner_environment
    )

    required_functions <- c(
      "run_qc",
      "resolve_quarto_executable",
      "read_quarto_version"
    )

    missing_functions <- required_functions[
      !vapply(
        required_functions,
        exists,
        logical(1),
        envir = runner_environment,
        mode = "function",
        inherits = FALSE
      )
    ]

    assert_true(
      length(missing_functions) == 0L,
      paste0(
        "The runner did not define required functions: ",
        paste(
          missing_functions,
          collapse = ", "
        )
      )
    )
  }
)


record_test(
  "Quarto executable is available",
  function() {
    quarto_path <- runner_environment$resolve_quarto_executable(
      PROJECT_ROOT
    )

    assert_true(
      is.character(quarto_path) &&
        length(quarto_path) == 1L &&
        !is.na(quarto_path) &&
        nzchar(quarto_path),
      "The resolved Quarto path is missing or invalid."
    )

    assert_true(
      file.exists(quarto_path),
      paste0(
        "The resolved Quarto executable does not exist: ",
        quarto_path
      )
    )

    quarto_version <- runner_environment$read_quarto_version(
      quarto_path
    )

    assert_true(
      length(quarto_version) == 1L &&
        !is.na(quarto_version) &&
        nzchar(quarto_version),
      "The Quarto version could not be read."
    )

    quarto_version
  }
)


## 7. Isolated end-to-end output root -----------------------------------------

output_root <- file.path(
  tempdir(),
  paste0(
    "QC espresso runner test ",
    Sys.getpid()
  )
)

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

on.exit(
  unlink(
    output_root,
    recursive = TRUE,
    force = TRUE
  ),
  add = TRUE
)


record_test(
  "Output root contains spaces",
  function() {
    assert_true(
      grepl(
        " ",
        output_root,
        fixed = TRUE
      ),
      "The runner integration test does not use a path containing spaces."
    )

    assert_true(
      dir.exists(output_root),
      "The isolated output root was not created."
    )
  }
)


## 8. End-to-end executions ----------------------------------------------------

record_test(
  "First end-to-end scientific render",
  function() {
    runner_result_1 <<- runner_environment$run_qc(
      report_path = report_path,
      fasta_path = fasta_path,
      output_root = output_root,
      open_report = FALSE
    )

    assert_true(
      inherits(
        runner_result_1,
        "qc_espresso_result"
      ),
      "The runner did not return a qc_espresso_result object."
    )

    assert_true(
      isTRUE(runner_result_1$success),
      "The first runner execution did not report success."
    )

    assert_true(
      file.exists(
        runner_result_1$html_report
      ),
      paste0(
        "The first rendered HTML report was not found: ",
        runner_result_1$html_report
      )
    )
  }
)


record_test(
  "Second end-to-end scientific render",
  function() {
    runner_result_2 <<- runner_environment$run_qc(
      report_path = report_path,
      fasta_path = fasta_path,
      output_root = output_root,
      open_report = FALSE
    )

    assert_true(
      inherits(
        runner_result_2,
        "qc_espresso_result"
      ),
      "The second runner execution did not return qc_espresso_result."
    )

    assert_true(
      isTRUE(runner_result_2$success),
      "The second runner execution did not report success."
    )

    assert_true(
      file.exists(
        runner_result_2$html_report
      ),
      paste0(
        "The second rendered HTML report was not found: ",
        runner_result_2$html_report
      )
    )
  }
)


record_test(
  "Two renders use isolated output directories",
  function() {
    first_output_directory <- normalizePath(
      runner_result_1$output_directory,
      winslash = "/",
      mustWork = TRUE
    )

    second_output_directory <- normalizePath(
      runner_result_2$output_directory,
      winslash = "/",
      mustWork = TRUE
    )

    assert_true(
      !identical(
        first_output_directory,
        second_output_directory
      ),
      "The two runner executions used the same output directory."
    )

    assert_true(
      dir.exists(first_output_directory) &&
        dir.exists(second_output_directory),
      "One or both execution output directories do not exist."
    )
  }
)


## 9. Result-object validation -------------------------------------------------

required_result_fields <- c(
  "success",
  "qc_espresso_version",
  "r_version",
  "quarto_version",
  "report_path",
  "fasta_path",
  "output_directory",
  "html_report",
  "generated_files",
  "started_at",
  "completed_at",
  "elapsed_seconds"
)


validate_result_object <- function(
    result_object,
    execution_label
) {
  missing_fields <- setdiff(
    required_result_fields,
    names(result_object)
  )

  assert_true(
    length(missing_fields) == 0L,
    paste0(
      execution_label,
      " result is missing required fields: ",
      paste(
        missing_fields,
        collapse = ", "
      )
    )
  )

  assert_true(
    isTRUE(result_object$success),
    paste0(
      execution_label,
      " result does not report success."
    )
  )

  assert_true(
    is.character(result_object$output_directory) &&
      length(result_object$output_directory) == 1L &&
      !is.na(result_object$output_directory) &&
      nzchar(result_object$output_directory) &&
      dir.exists(result_object$output_directory),
    paste0(
      execution_label,
      " output_directory is missing or invalid."
    )
  )

  assert_true(
    is.character(result_object$html_report) &&
      length(result_object$html_report) == 1L &&
      !is.na(result_object$html_report) &&
      nzchar(result_object$html_report) &&
      file.exists(result_object$html_report),
    paste0(
      execution_label,
      " html_report is missing or invalid."
    )
  )

  assert_true(
    is.character(result_object$generated_files),
    paste0(
      execution_label,
      " generated_files is not a character vector."
    )
  )

  assert_true(
    is.numeric(result_object$elapsed_seconds) &&
      length(result_object$elapsed_seconds) == 1L &&
      is.finite(result_object$elapsed_seconds) &&
      result_object$elapsed_seconds >= 0,
    paste0(
      execution_label,
      " elapsed_seconds is missing or invalid."
    )
  )

  assert_true(
    !is.na(result_object$started_at) &&
      !is.na(result_object$completed_at),
    paste0(
      execution_label,
      " execution timestamps are missing."
    )
  )

  assert_true(
    result_object$completed_at >=
      result_object$started_at,
    paste0(
      execution_label,
      " completed_at precedes started_at."
    )
  )

  invisible(TRUE)
}


record_test(
  "First result contains execution metadata",
  function() {
    validate_result_object(
      runner_result_1,
      "First execution"
    )
  }
)


record_test(
  "Second result contains execution metadata",
  function() {
    validate_result_object(
      runner_result_2,
      "Second execution"
    )
  }
)


record_test(
  "Result input paths match requested inputs",
  function() {
    expected_report_path <- normalizePath(
      report_path,
      winslash = "/",
      mustWork = TRUE
    )

    expected_fasta_path <- normalizePath(
      fasta_path,
      winslash = "/",
      mustWork = TRUE
    )

    observed_report_path <- normalizePath(
      runner_result_1$report_path,
      winslash = "/",
      mustWork = TRUE
    )

    observed_fasta_path <- normalizePath(
      runner_result_1$fasta_path,
      winslash = "/",
      mustWork = TRUE
    )

    assert_true(
      identical(
        observed_report_path,
        expected_report_path
      ),
      "The result object records a different DIA-NN report path."
    )

    assert_true(
      identical(
        observed_fasta_path,
        expected_fasta_path
      ),
      "The result object records a different FASTA path."
    )
  }
)


record_test(
  "Result versions are non-empty",
  function() {
    version_fields <- c(
      "qc_espresso_version",
      "r_version",
      "quarto_version"
    )

    invalid_fields <- version_fields[
      !vapply(
        version_fields,
        function(field_name) {
          field_value <- runner_result_1[[field_name]]

          is.character(field_value) &&
            length(field_value) == 1L &&
            !is.na(field_value) &&
            nzchar(
              trimws(field_value)
            )
        },
        logical(1)
      )
    ]

    assert_true(
      length(invalid_fields) == 0L,
      paste0(
        "The result contains missing or invalid version fields: ",
        paste(
          invalid_fields,
          collapse = ", "
        )
      )
    )
  }
)


## 10. Rendered HTML validation ------------------------------------------------

validate_html_report <- function(
    html_path,
    execution_label
) {
  html_info <- file.info(
    html_path
  )

  assert_true(
    !is.na(html_info$size) &&
      html_info$size > 10000,
    paste0(
      execution_label,
      " HTML report is unexpectedly small."
    )
  )

  html_lines <- readLines(
    html_path,
    warn = FALSE
  )

  assert_true(
    length(html_lines) > 0L,
    paste0(
      execution_label,
      " HTML report is empty."
    )
  )

  html_text <- paste(
    html_lines,
    collapse = "\n"
  )

  assert_true(
    grepl(
      "QC",
      html_text,
      fixed = TRUE
    ),
    paste0(
      execution_label,
      " HTML report does not contain the expected QC title marker."
    )
  )

  forbidden_markers <- c(
    "Stack trace",
    "Browse[1]",
    "Execution halted"
  )

  present_forbidden_markers <- forbidden_markers[
    vapply(
      forbidden_markers,
      grepl,
      logical(1),
      x = html_text,
      fixed = TRUE
    )
  ]

  assert_true(
    length(present_forbidden_markers) == 0L,
    paste0(
      execution_label,
      " HTML report contains failure markers: ",
      paste(
        present_forbidden_markers,
        collapse = ", "
      )
    )
  )

  invisible(TRUE)
}


record_test(
  "First rendered HTML is valid",
  function() {
    validate_html_report(
      runner_result_1$html_report,
      "First rendered"
    )
  }
)


record_test(
  "Second rendered HTML is valid",
  function() {
    validate_html_report(
      runner_result_2$html_report,
      "Second rendered"
    )
  }
)


record_test(
  "Rendered HTML reports are distinct files",
  function() {
    first_html <- normalizePath(
      runner_result_1$html_report,
      winslash = "/",
      mustWork = TRUE
    )

    second_html <- normalizePath(
      runner_result_2$html_report,
      winslash = "/",
      mustWork = TRUE
    )

    assert_true(
      !identical(
        first_html,
        second_html
      ),
      "The two executions resolved to the same HTML report file."
    )
  }
)


## 11. Generated-file validation ----------------------------------------------

validate_generated_files <- function(
    result_object,
    execution_label
) {
  generated_files <- unique(
    result_object$generated_files
  )

  generated_files <- generated_files[
    !is.na(generated_files) &
      nzchar(generated_files)
  ]

  assert_true(
    length(generated_files) > 0L,
    paste0(
      execution_label,
      " did not declare any generated files."
    )
  )

  missing_files <- generated_files[
    !file.exists(generated_files)
  ]

  assert_true(
    length(missing_files) == 0L,
    paste0(
      execution_label,
      " declared generated files that do not exist: ",
      paste(
        basename(missing_files),
        collapse = ", "
      )
    )
  )

  generated_info <- file.info(
    generated_files
  )

  empty_files <- generated_files[
    is.na(generated_info$size) |
      generated_info$size <= 0
  ]

  assert_true(
    length(empty_files) == 0L,
    paste0(
      execution_label,
      " generated empty or inaccessible files: ",
      paste(
        basename(empty_files),
        collapse = ", "
      )
    )
  )

  normalized_output_directory <- normalizePath(
    result_object$output_directory,
    winslash = "/",
    mustWork = TRUE
  )

  normalized_generated_files <- normalizePath(
    generated_files,
    winslash = "/",
    mustWork = TRUE
  )

  output_prefix <- paste0(
    normalized_output_directory,
    "/"
  )

  files_outside_output <- normalized_generated_files[
    !startsWith(
      normalized_generated_files,
      output_prefix
    ) &
      normalized_generated_files !=
      normalized_output_directory
  ]

  assert_true(
    length(files_outside_output) == 0L,
    paste0(
      execution_label,
      " generated files outside its isolated output directory: ",
      paste(
        basename(files_outside_output),
        collapse = ", "
      )
    )
  )

  invisible(TRUE)
}


record_test(
  "First execution generated files are valid",
  function() {
    validate_generated_files(
      runner_result_1,
      "First execution"
    )
  }
)


record_test(
  "Second execution generated files are valid",
  function() {
    validate_generated_files(
      runner_result_2,
      "Second execution"
    )
  }
)


record_test(
  "Generated-file inventories are isolated",
  function() {
    first_generated_files <- normalizePath(
      runner_result_1$generated_files,
      winslash = "/",
      mustWork = TRUE
    )

    second_generated_files <- normalizePath(
      runner_result_2$generated_files,
      winslash = "/",
      mustWork = TRUE
    )

    shared_files <- intersect(
      first_generated_files,
      second_generated_files
    )

    assert_true(
      length(shared_files) == 0L,
      paste0(
        "The two executions share generated files: ",
        paste(
          basename(shared_files),
          collapse = ", "
        )
      )
    )
  }
)


## 12. Scientific CSV validation ----------------------------------------------

validate_csv_file <- function(
    csv_path
) {
  csv_info <- file.info(
    csv_path
  )

  assert_true(
    !is.na(csv_info$size) &&
      csv_info$size > 0,
    paste0(
      "Generated CSV is empty or inaccessible: ",
      basename(csv_path)
    )
  )

  csv_header <- names(
    utils::read.csv(
      csv_path,
      nrows = 1L,
      check.names = FALSE
    )
  )

  assert_true(
    length(csv_header) > 0L,
    paste0(
      "Generated CSV has no readable columns: ",
      basename(csv_path)
    )
  )

  invisible(TRUE)
}


record_test(
  "Generated scientific CSV files are readable",
  function() {
    generated_csv_files <- runner_result_1$generated_files[
      tolower(
        tools::file_ext(
          runner_result_1$generated_files
        )
      ) == "csv"
    ]

    if (length(generated_csv_files) == 0L) {
      return(
        "No scientific CSV was generated for this example dataset."
      )
    }

    for (csv_path in generated_csv_files) {
      validate_csv_file(
        csv_path
      )
    }

    invisible(TRUE)
  }
)


record_test(
  "Recognized scientific CSV schemas are valid",
  function() {
    generated_files <- runner_result_1$generated_files

    generated_names <- basename(
      generated_files
    )

    figures_of_merit_files <- generated_files[
      grepl(
        "^figures_of_merit_.*[.]csv$",
        generated_names
      )
    ]

    coverage_files <- generated_files[
      grepl(
        "^protein_coverage_with_peptides_.*[.]csv$",
        generated_names
      )
    ]

    if (length(figures_of_merit_files) > 0L) {
      figures_columns <- names(
        utils::read.csv(
          figures_of_merit_files[[1L]],
          nrows = 1L,
          check.names = FALSE
        )
      )

      assert_true(
        "Run" %in% figures_columns,
        "The Figures of Merit CSV does not contain the Run column."
      )
    }

    if (length(coverage_files) > 0L) {
      coverage_columns <- names(
        utils::read.csv(
          coverage_files[[1L]],
          nrows = 1L,
          check.names = FALSE
        )
      )

      required_coverage_columns <- c(
        "Protein.Group",
        "Coverage_Pct"
      )

      missing_coverage_columns <- setdiff(
        required_coverage_columns,
        coverage_columns
      )

      assert_true(
        length(missing_coverage_columns) == 0L,
        paste0(
          "The protein coverage CSV is missing required columns: ",
          paste(
            missing_coverage_columns,
            collapse = ", "
          )
        )
      )
    }

    invisible(TRUE)
  }
)


## 13. Dynamic execution summary ----------------------------------------------

first_file_names <- basename(
  runner_result_1$generated_files
)

second_file_names <- basename(
  runner_result_2$generated_files
)


cat("\n")
cat("============================================================\n")
cat("RUNNER AND REPORT EXECUTION SUMMARY\n")
cat("============================================================\n")
cat(
  "First execution directory: ",
  runner_result_1$output_directory,
  "\n",
  sep = ""
)
cat(
  "Second execution directory: ",
  runner_result_2$output_directory,
  "\n",
  sep = ""
)
cat(
  "First execution elapsed seconds: ",
  format(
    runner_result_1$elapsed_seconds,
    digits = 8
  ),
  "\n",
  sep = ""
)
cat(
  "Second execution elapsed seconds: ",
  format(
    runner_result_2$elapsed_seconds,
    digits = 8
  ),
  "\n",
  sep = ""
)
cat(
  "First execution generated files: ",
  length(first_file_names),
  "\n",
  sep = ""
)

for (file_name in first_file_names) {
  cat(
    "  ",
    file_name,
    "\n",
    sep = ""
  )
}

cat(
  "Second execution generated files: ",
  length(second_file_names),
  "\n",
  sep = ""
)

for (file_name in second_file_names) {
  cat(
    "  ",
    file_name,
    "\n",
    sep = ""
  )
}

cat("============================================================\n")


## 14. Consolidated results ----------------------------------------------------

result_table <- do.call(
  rbind,
  results
)

passed_tests <- sum(
  result_table$Success
)

failed_tests <- sum(
  !result_table$Success
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


results_output <- file.path(
  output_dir,
  "runner-report-results.tsv"
)


utils::write.table(
  result_table,
  results_output,
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


cat("\n")
cat("============================================================\n")
cat("RUNNER AND REPORT TEST SUMMARY\n")
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
cat(
  "Detailed results: ",
  normalizePath(
    results_output,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


quit(
  status = if (failed_tests == 0L) {
    0L
  } else {
    1L
  }
)