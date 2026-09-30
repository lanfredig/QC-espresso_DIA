###############################################################################
# File: tests/run_audit.R
# Purpose:
#   Consolidated post-refactoring audit for QC-espresso.
#
# Usage:
#   Rscript --no-init-file tests/run_audit.R
###############################################################################

resolve_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg)) {
    script <- normalizePath(
      sub("^--file=", "", file_arg[[1L]]),
      winslash = "/",
      mustWork = TRUE
    )

    return(
      normalizePath(
        file.path(dirname(script), ".."),
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- resolve_root()
setwd(ROOT)
results <- list()

record <- function(category, test, status, details = "") {
  allowed <- c("PASS", "FAIL", "KNOWN ISSUE", "NOT IMPLEMENTED")

  if (!status %in% allowed) {
    stop("Invalid audit status: ", status, call. = FALSE)
  }

  results[[length(results) + 1L]] <<- data.frame(
    Category = category,
    Test = test,
    Status = status,
    Details = details,
    stringsAsFactors = FALSE
  )

  cat(sprintf("%-18s %-50s %s\n", category, test, status))

  if (nzchar(details) && status != "PASS") {
    cat("  ", details, "\n", sep = "")
  }
}

run_check <- function(category, test, fun) {
  tryCatch(
    {
      value <- fun()
      details <- if (is.character(value) && length(value) == 1L) value else ""
      record(category, test, "PASS", details)
    },
    error = function(e) {
      record(category, test, "FAIL", conditionMessage(e))
    }
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

run_suite <- function(category, label, script_name) {
  script_path <- file.path("tests", script_name)

  if (!file.exists(script_path)) {
    record(category, label, "FAIL", paste("Test file not found:", script_path))
    return(invisible(NULL))
  }

  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = c("--no-init-file", shQuote(script_path)),
    stdout = log_file,
    stderr = log_file
  )

  log_lines <- readLines(log_file, warn = FALSE)

  if (identical(status, 0L)) {
    record(category, label, "PASS")
  } else {
    important_lines <- log_lines[
      grepl(
        "FAIL|Reason:|FINAL STATUS|Tests failed|Error|Execution halted",
        log_lines,
        ignore.case = TRUE
      )
    ]

    if (!length(important_lines)) {
      important_lines <- tail(log_lines, 15L)
    }

    record(
      category,
      label,
      "FAIL",
      paste(important_lines, collapse = "\n")
    )
  }

  invisible(NULL)
}

cat("\n============================================================\n")
cat("QC-espresso consolidated post-refactoring audit\n")
cat("============================================================\n")
cat("Project: ", ROOT, "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")
cat("============================================================\n\n")

## 1. Repository and environment ----------------------------------------------

run_check("Repository", "Required scientific files", function() {
  required <- c(
    "QC-espresso_DIA.Rproj",
    "VERSION",
    "run_qc.R",
    "renv.lock",
    ".Rprofile",
    ".gitignore",
    file.path("report", "qc_espresso_scientific.qmd")
  )

  missing <- required[!file.exists(required)]
  assert(!length(missing), paste("Missing:", paste(missing, collapse = ", ")))
})

run_check("Repository", "No tracked launcher artifacts", function() {
  tracked <- system2("git", "ls-files", stdout = TRUE, stderr = TRUE)
  forbidden <- tracked[grepl("[.](exe|bat|ps1)$", tracked, ignore.case = TRUE)]
  assert(!length(forbidden), paste("Tracked launcher artifacts:", paste(forbidden, collapse = ", ")))
})

run_check("Repository", "No portable runtime in repository", function() {
  assert(
    !dir.exists("runtime") && !dir.exists(file.path("R", "runtime")),
    "Portable runtime is inside QC-espresso_DIA."
  )
})

run_check("Repository", "No tracked R or Quarto residue", function() {
  tracked <- system2("git", "ls-files", stdout = TRUE, stderr = TRUE)
  forbidden <- tracked[
    basename(tracked) %in% c(".Rhistory", ".RData") |
      grepl("(^|/)[.]quarto(/|$)", tracked)
  ]
  assert(!length(forbidden), paste("Tracked residue:", paste(forbidden, collapse = ", ")))
})

run_check("Environment", "R 4.6.0", function() {
  observed <- paste(R.version$major, R.version$minor, sep = ".")
  assert(identical(observed, "4.6.0"), paste("Observed", observed))
})

run_check("Environment", "Lockfile records R 4.6.0", function() {
  assert(requireNamespace("jsonlite", quietly = TRUE), "jsonlite unavailable")
  lock <- jsonlite::read_json("renv.lock", simplifyVector = FALSE)
  assert(
    identical(lock[["R"]][["Version"]], "4.6.0"),
    "renv.lock does not record R 4.6.0"
  )
})

run_check("Environment", "Portable Quarto is discoverable", function() {
  runner_environment <- new.env(parent = globalenv())
  sys.source("run_qc.R", envir = runner_environment)
  quarto <- runner_environment$resolve_quarto_executable(ROOT)
  assert(file.exists(quarto), paste("Resolved Quarto was not found:", quarto))
  version <- runner_environment$read_quarto_version(quarto)
  assert(!is.na(version) && nzchar(version), "Quarto version could not be read")
})

## 2. Static implementation checks --------------------------------------------

run_check("Science", "A1 repeated PG.MaxLFQ formula absent", function() {
  text <- paste(readLines(file.path("R", "qc_abundance.R"), warn = FALSE), collapse = "\n")
  assert(
    !grepl("PG.MaxLFQ / sum(PG.MaxLFQ", text, fixed = TRUE),
    "Legacy repeated-PG.MaxLFQ proportion formula remains"
  )
})

run_check("Science", "A3 multiprotein annotation implemented", function() {
  text <- paste(readLines(file.path("R", "qc_abundance.R"), warn = FALSE), collapse = "\n")
  required <- c(
    "MULTIPLE_ANNOTATIONS",
    "PARTIAL_FASTA_MATCH",
    "relationship = \"many-to-one\"",
    "Signal_Difference"
  )
  missing <- required[!vapply(required, grepl, logical(1), x = text, fixed = TRUE)]
  assert(!length(missing), paste("Missing implementation markers:", paste(missing, collapse = ", ")))
})

run_check("Science", "Missed-cleavage conditional rule implemented", function() {
  for (file in c("qc_peptides.R", "qc_figures_of_merit.R")) {
    text <- paste(readLines(file.path("R", file), warn = FALSE), collapse = "\n")
    assert(
      grepl("count_missed_tryptic_cleavages", text, fixed = TRUE),
      paste("Corrected helper not used in", file)
    )
  }
})

run_check("Science", "Completeness avoids dense pivot", function() {
  text <- paste(readLines(file.path("R", "qc_completeness.R"), warn = FALSE), collapse = "\n")
  assert(!grepl("pivot_wider", text, fixed = TRUE), "Dense pivot remains in qc_completeness.R")
})

run_check("Visualization", "Missed-cleavage axis is not truncated", function() {
  text <- paste(readLines(file.path("R", "qc_peptides.R"), warn = FALSE), collapse = "\n")
  assert(!grepl("limits = c(0, 0.5)", text, fixed = TRUE), "Fixed 0-50% limit remains")
})

run_check("Visualization", "Mass-error visual range intentionally retained", function() {
  text <- paste(readLines(file.path("R", "qc_mass_error.R"), warn = FALSE), collapse = "\n")
  assert(
    grepl("xlim = c(-10, 10)", text, fixed = TRUE),
    "The approved -10 to 10 ppm visual range is absent"
  )
})

run_check("Scientific mode", "Scientific QMD implemented", function() {
  qmd <- file.path("report", "qc_espresso_scientific.qmd")
  assert(file.exists(qmd), "Scientific QMD is absent")
  text <- paste(readLines(qmd, warn = FALSE), collapse = "\n")
  assert(!grepl("quarto_exec:", text, fixed = TRUE), "Obsolete quarto_exec parameter remains")
  assert(!grepl("setwd(", text, fixed = TRUE), "Scientific QMD changes the global working directory")
})

## 3. Consolidated executable suites ------------------------------------------

suite_inventory <- list(
  c("Baseline", "General automated validation", "run_tests.R"),
  c("Science", "Multiprotein abundance validation", "test_multiprotein_abundance.R"),
  c("Science", "Missed-cleavage validation", "test_missed_cleavages.R"),
  c("Robustness", "Optional-section validation", "test_optional_sections.R"),
  c("Science", "Completeness validation", "test_completeness.R"),
  c("Science", "Coverage validation", "test_coverage.R"),
  c("Science", "Intensity validation", "test_intensity.R"),
  c("Integration", "Scientific core validation", "test_scientific_core.R"),
  c("Scientific mode", "End-to-end scientific render", "test_runner_report.R")
)

for (suite in suite_inventory) {
  run_suite(suite[[1]], suite[[2]], suite[[3]])
}

## 4. Consolidated audit artifacts --------------------------------------------

result_table <- do.call(rbind, results)
output_dir <- file.path("tests", "test-results", "latest")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  result_table,
  file.path(output_dir, "audit-results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

status_levels <- c("PASS", "FAIL", "KNOWN ISSUE", "NOT IMPLEMENTED")
counts <- table(factor(result_table$Status, levels = status_levels))

final_status <- if (
  counts[["FAIL"]] == 0L &&
    counts[["KNOWN ISSUE"]] == 0L &&
    counts[["NOT IMPLEMENTED"]] == 0L
) {
  "PASS"
} else if (counts[["FAIL"]] == 0L) {
  "PASS WITH DOCUMENTED GAPS"
} else {
  "FAIL"
}

summary_lines <- c(
  "QC-espresso consolidated post-refactoring audit",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Project:", ROOT),
  paste("R:", R.version.string),
  "",
  paste("PASS:", counts[["PASS"]]),
  paste("FAIL:", counts[["FAIL"]]),
  paste("KNOWN ISSUE:", counts[["KNOWN ISSUE"]]),
  paste("NOT IMPLEMENTED:", counts[["NOT IMPLEMENTED"]]),
  "",
  paste("FINAL STATUS:", final_status)
)

writeLines(
  summary_lines,
  file.path(output_dir, "audit-summary.txt")
)

cat("\n============================================================\n")
cat("AUDIT SUMMARY\n")
cat("============================================================\n")
cat(paste(summary_lines[6:length(summary_lines)], collapse = "\n"), "\n")
cat("============================================================\n")
cat(
  "Detailed results: ",
  normalizePath(
    file.path(output_dir, "audit-results.tsv"),
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)

quit(status = if (identical(final_status, "PASS")) 0L else 1L)
