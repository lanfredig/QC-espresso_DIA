###############################################################################
# File: tests/test_runner_report.R
# Purpose:
#   Validate the scientific runner and Quarto report using the article example.
###############################################################################

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
results <- list()

record <- function(name, fun) {
  outcome <- tryCatch(
    {
      details <- fun()
      list(success = TRUE, details = if (is.character(details)) paste(details, collapse = " ") else "")
    },
    error = function(e) list(success = FALSE, details = conditionMessage(e))
  )
  results[[length(results) + 1L]] <<- data.frame(
    Test = name,
    Success = outcome$success,
    Details = outcome$details,
    stringsAsFactors = FALSE
  )
  cat(sprintf("%-58s %s\n", name, if (outcome$success) "PASS" else "FAIL"))
  if (!outcome$success) cat("  Reason: ", outcome$details, "\n", sep = "")
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

qmd_path <- file.path("report", "qc_espresso_scientific.qmd")
runner_path <- "run_qc.R"
report_path <- file.path("example-data", "dataset_4_contaminans.parquet")
fasta_path <- file.path("example-data", "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta")

record("Scientific runner exists", function() {
  assert(file.exists(runner_path), "run_qc.R was not found")
})

record("Scientific QMD exists", function() {
  assert(file.exists(qmd_path), "report/qc_espresso_scientific.qmd was not found")
})

record("Runner and QMD syntax is readable", function() {
  parse(runner_path)
  text <- readLines(qmd_path, warn = FALSE)
  assert(length(text) > 0L, "Scientific QMD is empty")
})

record("Scientific QMD declares required parameters", function() {
  text <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  for (token in c("report_path:", "fasta_path:", "output_dir:")) {
    assert(grepl(token, text, fixed = TRUE), paste("Missing QMD parameter:", token))
  }
})

record("Scientific QMD does not retain obsolete quarto_exec parameter", function() {
  text <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  assert(!grepl("quarto_exec:", text, fixed = TRUE), "Obsolete quarto_exec parameter remains in the QMD")
})

record("Scientific QMD reads version from VERSION", function() {
  text <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  assert(grepl("VERSION", text, fixed = TRUE), "The QMD does not read the project VERSION file")
})

record("Scientific QMD avoids project-global setwd", function() {
  text <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  assert(!grepl("setwd(", text, fixed = TRUE), "The QMD still changes the global working directory")
})

record("Quarto executable is available", function() {
  source(runner_path, local = TRUE)
  quarto <- resolve_quarto_executable()
  assert(file.exists(quarto), paste("Resolved Quarto does not exist:", quarto))
  version <- read_quarto_version(quarto)
  assert(!is.na(version) && nzchar(version), "Quarto version could not be read")
  version
})

# End-to-end execution uses a root with spaces and runs twice to verify isolation.
output_root <- file.path(tempdir(), paste0("QC espresso runner test ", Sys.getpid()))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(output_root, recursive = TRUE, force = TRUE), add = TRUE)

runner_result_1 <- NULL
runner_result_2 <- NULL

record("First end-to-end scientific render", function() {
  source(runner_path, local = globalenv())
  runner_result_1 <<- run_qc(
    report_path = report_path,
    fasta_path = fasta_path,
    output_root = output_root,
    open_report = FALSE
  )
  assert(inherits(runner_result_1, "qc_espresso_result"), "Runner did not return qc_espresso_result")
  assert(isTRUE(runner_result_1$success), "Runner did not report success")
  assert(file.exists(runner_result_1$html_report), "Rendered HTML was not found")
})

record("Second render creates an isolated directory", function() {
  runner_result_2 <<- run_qc(
    report_path = report_path,
    fasta_path = fasta_path,
    output_root = output_root,
    open_report = FALSE
  )
  assert(file.exists(runner_result_2$html_report), "Second rendered HTML was not found")
  assert(
    !identical(runner_result_1$output_directory, runner_result_2$output_directory),
    "Two executions used the same output directory"
  )
})

record("Rendered HTML is non-empty and self-contained", function() {
  info <- file.info(runner_result_1$html_report)
  assert(!is.na(info$size) && info$size > 10000, "Rendered HTML is unexpectedly small")
  html <- paste(readLines(runner_result_1$html_report, warn = FALSE), collapse = "\n")
  assert(grepl("QC", html, fixed = TRUE), "Rendered HTML does not contain the report title")
  assert(!grepl("Stack trace", html, fixed = TRUE), "Rendered HTML contains a stack trace")
  assert(!grepl("Browse[1]", html, fixed = TRUE), "Rendered HTML contains an R browser prompt")
})

record("Expected scientific CSV exports exist", function() {
  files <- basename(runner_result_1$generated_files)
  assert(any(grepl("^figures_of_merit_.*[.]csv$", files)), "Figures of Merit CSV was not generated")
  assert(any(grepl("^protein_coverage_with_peptides_.*[.]csv$", files)), "Coverage CSV was not generated")
})

record("Result object contains execution metadata", function() {
  required <- c(
    "success", "qc_espresso_version", "r_version", "quarto_version",
    "report_path", "fasta_path", "output_directory", "html_report",
    "generated_files", "started_at", "completed_at", "elapsed_seconds"
  )
  assert(all(required %in% names(runner_result_1)), "Runner result is missing required fields")
})

result_table <- do.call(rbind, results)
passed <- sum(result_table$Success)
failed <- sum(!result_table$Success)

output_dir <- file.path("tests", "test-results", "latest")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  result_table,
  file.path(output_dir, "runner-report-results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

cat("\n============================================================\n")
cat("RUNNER AND REPORT TEST SUMMARY\n")
cat("============================================================\n")
cat("Tests passed: ", passed, "\n", sep = "")
cat("Tests failed: ", failed, "\n", sep = "")
cat("FINAL STATUS: ", if (failed == 0L) "PASS" else "FAIL", "\n", sep = "")
cat("============================================================\n")

quit(status = if (failed == 0L) 0L else 1L)
