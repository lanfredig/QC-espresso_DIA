###############################################################################
# File: tests/run_audit.R
#
# Purpose:
#   Execute and consolidate the complete QC-espresso validation suite without
#   duplicating the scientific, structural, or integration checks implemented
#   by the individual test scripts.
#
# Validation scope:
#   - Confirm that every required validation suite is available and non-empty.
#   - Execute each validation suite once in an independent R process.
#   - Preserve the console output of every suite in an individual log file.
#   - Record exit status, execution time, and failure details for each suite.
#   - Consolidate all suite results into audit TSV and text-summary outputs.
#   - Return a non-zero process status when one or more suites fail.
#
# Result categories:
#   PASS  The validation suite completed with exit status zero.
#   FAIL  The validation suite was unavailable or returned a non-zero status.
#
# Usage:
#   Rscript --no-init-file tests/run_audit.R
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
    audit_file <- normalizePath(
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
          dirname(audit_file),
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
        "Run this audit from the project root."
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


## 2. Audit configuration ------------------------------------------------------

output_dir <- file.path(
  "tests",
  "test-results",
  "latest"
)

log_dir <- file.path(
  output_dir,
  "audit-logs"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


rscript_path <- file.path(
  R.home("bin"),
  "Rscript"
)

if (.Platform$OS.type == "windows") {
  rscript_path <- paste0(
    rscript_path,
    ".exe"
  )
}

rscript_path <- normalizePath(
  rscript_path,
  winslash = "/",
  mustWork = TRUE
)


suite_inventory <- data.frame(
  Order = seq_len(9L),
  Category = c(
    "Baseline",
    "Science",
    "Science",
    "Robustness",
    "Science",
    "Science",
    "Science",
    "Integration",
    "Scientific mode"
  ),
  Label = c(
    "Repository and input baseline",
    "Multi-protein abundance validation",
    "Missed-cleavage validation",
    "Optional-section validation",
    "Completeness validation",
    "Coverage validation",
    "Intensity validation",
    "Scientific core integration",
    "End-to-end scientific render"
  ),
  Script = c(
    "run_tests.R",
    "test_multiprotein_abundance.R",
    "test_missed_cleavages.R",
    "test_optional_sections.R",
    "test_completeness.R",
    "test_coverage.R",
    "test_intensity.R",
    "test_scientific_core.R",
    "test_runner_report.R"
  ),
  stringsAsFactors = FALSE
)


suite_inventory$Path <- file.path(
  "tests",
  suite_inventory$Script
)

suite_inventory$Log <- file.path(
  log_dir,
  paste0(
    tools::file_path_sans_ext(
      suite_inventory$Script
    ),
    ".log"
  )
)


## 3. Minimal audit framework --------------------------------------------------

audit_results <- list()


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


record_result <- function(
    category,
    test,
    script,
    status,
    exit_status,
    elapsed_seconds,
    log_file,
    details = ""
) {
  allowed_statuses <- c(
    "PASS",
    "FAIL"
  )

  if (!status %in% allowed_statuses) {
    stop(
      paste0(
        "Invalid audit status: ",
        status
      ),
      call. = FALSE
    )
  }

  result <- data.frame(
    Category = category,
    Test = test,
    Script = script,
    Status = status,
    Exit_Status = exit_status,
    Elapsed_Seconds = elapsed_seconds,
    Log_File = log_file,
    Details = details,
    stringsAsFactors = FALSE
  )

  audit_results[[length(audit_results) + 1L]] <<- result

  cat(
    sprintf(
      "%-18s %-44s %s\n",
      category,
      test,
      status
    )
  )

  if (
    identical(
      status,
      "FAIL"
    ) &&
    nzchar(details)
  ) {
    cat(
      "  ",
      details,
      "\n",
      sep = ""
    )
  }

  invisible(result)
}


extract_failure_details <- function(
    log_lines
) {
  if (length(log_lines) == 0L) {
    return(
      "The suite returned no console output."
    )
  }

  important_patterns <- c(
    "FINAL STATUS:",
    "Tests failed:",
    "Reason:",
    "FAIL",
    "Error",
    "Execution halted"
  )

  important_lines <- log_lines[
    vapply(
      log_lines,
      function(log_line) {
        any(
          vapply(
            important_patterns,
            grepl,
            logical(1),
            x = log_line,
            ignore.case = TRUE,
            fixed = TRUE
          )
        )
      },
      logical(1)
    )
  ]

  important_lines <- important_lines[
    nzchar(
      trimws(
        important_lines
      )
    )
  ]

  if (length(important_lines) == 0L) {
    important_lines <- tail(
      log_lines,
      20L
    )
  }

  important_lines <- unique(
    important_lines
  )

  paste(
    tail(
      important_lines,
      20L
    ),
    collapse = "\n"
  )
}


## 4. Audit header -------------------------------------------------------------

audit_started_at <- Sys.time()


cat("\n")
cat("============================================================\n")
cat("QC-ESPRESSO CONSOLIDATED VALIDATION AUDIT\n")
cat("============================================================\n")
cat(
  "Project: ",
  PROJECT_ROOT,
  "\n",
  sep = ""
)
cat(
  "R version: ",
  R.version.string,
  "\n",
  sep = ""
)
cat(
  "Rscript: ",
  rscript_path,
  "\n",
  sep = ""
)
cat(
  "Suites scheduled: ",
  nrow(
    suite_inventory
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")
cat("\n")


## 5. Validate suite inventory -------------------------------------------------

missing_suite_files <- suite_inventory$Path[
  !file.exists(
    suite_inventory$Path
  )
]

empty_suite_files <- character()

existing_suite_paths <- suite_inventory$Path[
  file.exists(
    suite_inventory$Path
  )
]

if (length(existing_suite_paths) > 0L) {
  existing_suite_sizes <- file.info(
    existing_suite_paths
  )$size

  empty_suite_files <- existing_suite_paths[
    is.na(existing_suite_sizes) |
      existing_suite_sizes <= 0
  ]
}


if (length(missing_suite_files) > 0L) {
  stop(
    paste0(
      "Required validation suites are missing: ",
      paste(
        missing_suite_files,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


if (length(empty_suite_files) > 0L) {
  stop(
    paste0(
      "Required validation suites are empty or inaccessible: ",
      paste(
        empty_suite_files,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


deprecated_script <- file.path(
  "tests",
  "diagnose_missed_cleavages.R"
)

assert_true(
  !file.exists(
    deprecated_script
  ),
  paste0(
    "The deprecated diagnostic script remains in the repository: ",
    deprecated_script
  )
)


## 6. Execute each suite once --------------------------------------------------

for (
  suite_index in
  seq_len(
    nrow(
      suite_inventory
    )
  )
) {
  suite_category <- suite_inventory$Category[
    suite_index
  ]

  suite_label <- suite_inventory$Label[
    suite_index
  ]

  suite_script <- suite_inventory$Script[
    suite_index
  ]

  suite_path <- suite_inventory$Path[
    suite_index
  ]

  suite_log <- suite_inventory$Log[
    suite_index
  ]

  if (file.exists(suite_log)) {
    unlink(
      suite_log,
      force = TRUE
    )
  }

  suite_started_at <- Sys.time()

  execution_error <- NULL

  exit_status <- tryCatch(
    system2(
      command = rscript_path,
      args = c(
        "--no-init-file",
        shQuote(
          suite_path
        )
      ),
      stdout = suite_log,
      stderr = suite_log,
      wait = TRUE
    ),
    error = function(error_condition) {
      execution_error <<- conditionMessage(
        error_condition
      )

      1L
    }
  )

  suite_elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      suite_started_at,
      units = "secs"
    )
  )

  if (
    is.null(exit_status) ||
    length(exit_status) == 0L ||
    is.na(exit_status)
  ) {
    exit_status <- 0L
  }

  log_lines <- if (file.exists(suite_log)) {
    readLines(
      suite_log,
      warn = FALSE
    )
  } else {
    character()
  }

  suite_status <- if (
    identical(
      as.integer(exit_status),
      0L
    )
  ) {
    "PASS"
  } else {
    "FAIL"
  }

  suite_details <- if (
    identical(
      suite_status,
      "PASS"
    )
  ) {
    ""
  } else if (
    !is.null(execution_error) &&
    nzchar(execution_error)
  ) {
    execution_error
  } else {
    extract_failure_details(
      log_lines
    )
  }

  record_result(
    category = suite_category,
    test = suite_label,
    script = suite_script,
    status = suite_status,
    exit_status = as.integer(
      exit_status
    ),
    elapsed_seconds = suite_elapsed_seconds,
    log_file = normalizePath(
      suite_log,
      winslash = "/",
      mustWork = file.exists(
        suite_log
      )
    ),
    details = suite_details
  )
}


## 7. Consolidated audit artifacts --------------------------------------------

result_table <- do.call(
  rbind,
  audit_results
)


passed_suites <- sum(
  result_table$Status ==
    "PASS"
)

failed_suites <- sum(
  result_table$Status ==
    "FAIL"
)


audit_completed_at <- Sys.time()

audit_elapsed_seconds <- as.numeric(
  difftime(
    audit_completed_at,
    audit_started_at,
    units = "secs"
  )
)


final_status <- if (
  failed_suites == 0L
) {
  "PASS"
} else {
  "FAIL"
}


audit_results_file <- file.path(
  output_dir,
  "audit-results.tsv"
)

audit_summary_file <- file.path(
  output_dir,
  "audit-summary.txt"
)


utils::write.table(
  result_table,
  audit_results_file,
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


summary_lines <- c(
  "QC-espresso consolidated validation audit",
  paste(
    "Generated:",
    format(
      audit_completed_at,
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste(
    "Project:",
    PROJECT_ROOT
  ),
  paste(
    "R:",
    R.version.string
  ),
  paste(
    "Suites executed:",
    nrow(
      result_table
    )
  ),
  paste(
    "Suites passed:",
    passed_suites
  ),
  paste(
    "Suites failed:",
    failed_suites
  ),
  paste(
    "Elapsed seconds:",
    format(
      audit_elapsed_seconds,
      digits = 10
    )
  ),
  paste(
    "FINAL STATUS:",
    final_status
  )
)


writeLines(
  summary_lines,
  audit_summary_file
)


## 8. Console summary ----------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("AUDIT SUMMARY\n")
cat("============================================================\n")
cat(
  "Suites executed: ",
  nrow(
    result_table
  ),
  "\n",
  sep = ""
)
cat(
  "Suites passed: ",
  passed_suites,
  "\n",
  sep = ""
)
cat(
  "Suites failed: ",
  failed_suites,
  "\n",
  sep = ""
)
cat(
  "Elapsed seconds: ",
  format(
    audit_elapsed_seconds,
    digits = 10
  ),
  "\n",
  sep = ""
)
cat(
  "FINAL STATUS: ",
  final_status,
  "\n",
  sep = ""
)
cat("============================================================\n")
cat(
  "Detailed results: ",
  normalizePath(
    audit_results_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "Summary: ",
  normalizePath(
    audit_summary_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "Suite logs: ",
  normalizePath(
    log_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


quit(
  status = if (
    identical(
      final_status,
      "PASS"
    )
  ) {
    0L
  } else {
    1L
  }
)