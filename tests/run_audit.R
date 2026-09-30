###############################################################################
# File: tests/run_audit.R
# Purpose:
#   Consolidated pre-refactoring audit for QC-espresso.
#
# Status categories:
#   PASS            Requirement is satisfied.
#   FAIL            Unexpected failure that blocks progression.
#   KNOWN ISSUE     Reproduced issue already identified in the technical review.
#   NOT IMPLEMENTED Planned scientific-mode component not yet available.
#
# Usage:
#   Rscript --no-init-file tests/run_audit.R
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

record <- function(category, test, status, details = "") {
  allowed <- c("PASS", "FAIL", "KNOWN ISSUE", "NOT IMPLEMENTED")
  if (!status %in% allowed) stop("Invalid audit status: ", status, call. = FALSE)
  results[[length(results) + 1L]] <<- data.frame(
    Category = category,
    Test = test,
    Status = status,
    Details = details,
    stringsAsFactors = FALSE
  )
  cat(sprintf("%-16s %-46s %s\n", category, test, status))
  if (nzchar(details) && status != "PASS") cat("  ", details, "\n", sep = "")
}

run_check <- function(category, test, fun) {
  tryCatch(
    {
      value <- fun()
      details <- if (is.character(value) && length(value) == 1L) value else ""
      record(category, test, "PASS", details)
    },
    error = function(e) record(category, test, "FAIL", conditionMessage(e))
  )
}

assert <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

read_script <- function(name) paste(readLines(file.path("R", name), warn = FALSE), collapse = "\n")

cat("\n============================================================\n")
cat("QC-espresso consolidated pre-refactoring audit\n")
cat("============================================================\n")
cat("Project:", ROOT, "\n")
cat("R version:", R.version.string, "\n")
cat("============================================================\n\n")

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(microseq)
  library(stringr)
  library(tibble)
})

source(file.path("R", "validate_inputs.R"), local = globalenv())
source(file.path("R", "globals.R"), local = globalenv())

PARQUET <- file.path("example-data", "dataset_4_contaminans.parquet")
FASTA <- file.path("example-data", "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta")

# -----------------------------------------------------------------------------
# Repository and environment
# -----------------------------------------------------------------------------

run_check("Repository", "Required scientific files", function() {
  files <- c("QC-espresso_DIA.Rproj", "VERSION", "run_qc.R", "renv.lock", ".Rprofile", ".gitignore")
  missing <- files[!file.exists(files)]
  assert(!length(missing), paste("Missing:", paste(missing, collapse = ", ")))
})

run_check("Repository", "No tracked user-launcher artifacts", function() {
  tracked <- system2("git", c("ls-files"), stdout = TRUE, stderr = TRUE)
  forbidden <- tracked[grepl("[.](exe|bat|ps1)$", tracked, ignore.case = TRUE)]
  assert(!length(forbidden), paste("Tracked launcher artifacts:", paste(forbidden, collapse = ", ")))
})

run_check("Repository", "No portable runtime in repository", function() {
  assert(!dir.exists("runtime") && !dir.exists(file.path("R", "runtime")), "Portable runtime is inside the repository.")
})

run_check("Repository", "No tracked R or Quarto residue", function() {
  tracked <- system2("git", c("ls-files"), stdout = TRUE, stderr = TRUE)
  forbidden_names <- c(".Rhistory", ".RData")
  forbidden <- tracked[basename(tracked) %in% forbidden_names | grepl("(^|/)[.]quarto(/|$)", tracked)]
  assert(!length(forbidden), paste("Tracked residue:", paste(forbidden, collapse = ", ")))
})

run_check("Environment", "R 4.6.0", function() {
  observed <- paste(R.version$major, R.version$minor, sep = ".")
  assert(identical(observed, "4.6.0"), paste("Observed", observed))
})

run_check("Environment", "Lockfile R version", function() {
  assert(requireNamespace("jsonlite", quietly = TRUE), "jsonlite unavailable")
  lock <- jsonlite::read_json("renv.lock", simplifyVector = FALSE)
  assert(identical(lock[["R"]][["Version"]], "4.6.0"), "renv.lock does not record R 4.6.0")
})

# -----------------------------------------------------------------------------
# Input validation negative cases
# -----------------------------------------------------------------------------

run_check("Validation", "Missing Parquet rejected", function() {
  ok <- inherits(try(validate_report_file(file.path(tempdir(), "missing.parquet")), silent = TRUE), "try-error")
  assert(ok, "Missing Parquet was accepted")
})

run_check("Validation", "Wrong Parquet extension rejected", function() {
  path <- tempfile(fileext = ".txt"); writeLines("x", path)
  on.exit(unlink(path), add = TRUE)
  ok <- inherits(try(validate_report_file(path), silent = TRUE), "try-error")
  assert(ok, "Wrong Parquet extension was accepted")
})

run_check("Validation", "Empty Parquet rejected", function() {
  path <- tempfile(fileext = ".parquet"); file.create(path)
  on.exit(unlink(path), add = TRUE)
  ok <- inherits(try(validate_report_file(path), silent = TRUE), "try-error")
  assert(ok, "Empty Parquet was accepted")
})

run_check("Validation", "Invalid Parquet content rejected", function() {
  path <- tempfile(fileext = ".parquet"); writeLines("not parquet", path)
  on.exit(unlink(path), add = TRUE)
  path <- validate_report_file(path)
  ok <- inherits(try(validate_report_schema(path), silent = TRUE), "try-error")
  assert(ok, "Invalid Parquet content was accepted")
})

run_check("Validation", "Missing required column reported", function() {
  path <- tempfile(fileext = ".parquet")
  arrow::write_parquet(data.frame(Run = "A", PG.Q.Value = 0.001), path)
  on.exit(unlink(path), add = TRUE)
  message <- tryCatch({ validate_report_schema(path); "" }, error = function(e) conditionMessage(e))
  assert(grepl("Decoy", message, fixed = TRUE), paste("Unexpected message:", message))
})

run_check("Validation", "Empty FASTA rejected", function() {
  path <- tempfile(fileext = ".fasta"); file.create(path)
  on.exit(unlink(path), add = TRUE)
  ok <- inherits(try(validate_fasta_file(path), silent = TRUE), "try-error")
  assert(ok, "Empty FASTA was accepted")
})

run_check("Validation", "FASTA without header rejected", function() {
  path <- tempfile(fileext = ".fasta"); writeLines("PEPTIDE", path)
  on.exit(unlink(path), add = TRUE)
  ok <- inherits(try(validate_fasta_file(path), silent = TRUE), "try-error")
  assert(ok, "FASTA without header was accepted")
})

# -----------------------------------------------------------------------------
# I/O equivalence: current materialized method vs proposed Arrow pushdown
# -----------------------------------------------------------------------------

run_check("I/O", "Old and Arrow-filtered row counts", function() {
  old <- arrow::read_parquet(PARQUET) |> filter(Decoy == 0, PG.Q.Value <= 0.01)
  new <- arrow::open_dataset(PARQUET, format = "parquet") |> filter(Decoy == 0, PG.Q.Value <= 0.01) |> collect()
  assert(nrow(old) == nrow(new), paste("old=", nrow(old), " new=", nrow(new), sep = ""))
  assert(nrow(old) == 28748L, paste("Baseline changed:", nrow(old)))
})

run_check("I/O", "Old and Arrow-filtered columns", function() {
  old <- arrow::read_parquet(PARQUET) |> filter(Decoy == 0, PG.Q.Value <= 0.01)
  new <- arrow::open_dataset(PARQUET, format = "parquet") |> filter(Decoy == 0, PG.Q.Value <= 0.01) |> collect()
  assert(identical(names(old), names(new)), "Column names or order differ")
  old_classes <- vapply(old, function(x) paste(class(x), collapse = "/"), character(1))
  new_classes <- vapply(new, function(x) paste(class(x), collapse = "/"), character(1))
  assert(identical(old_classes, new_classes), "Column classes differ")
})

run_check("I/O", "Old and Arrow-filtered values", function() {
  old <- arrow::read_parquet(PARQUET) |> filter(Decoy == 0, PG.Q.Value <= 0.01)
  new <- arrow::open_dataset(PARQUET, format = "parquet") |> filter(Decoy == 0, PG.Q.Value <= 0.01) |> collect()
  assert(isTRUE(all.equal(old, new, check.attributes = TRUE, check.names = TRUE)), "Filtered tables are not equivalent")
})

run_check("I/O", "Global filter invariants", function() {
  filtered <- arrow::open_dataset(PARQUET, format = "parquet") |> filter(Decoy == 0, PG.Q.Value <= 0.01) |> collect()
  assert(all(filtered$Decoy == 0, na.rm = TRUE), "A decoy row was retained")
  assert(all(filtered$PG.Q.Value <= 0.01, na.rm = TRUE), "A row above PG.Q.Value 0.01 was retained")
  assert(n_distinct(filtered$Run) == 16L, "Retained-run count changed")
})

run_check("I/O", "Removed-run identity", function() {
  raw_runs <- arrow::open_dataset(PARQUET, format = "parquet") |> select(Run) |> distinct() |> collect() |> pull(Run)
  kept_runs <- arrow::open_dataset(PARQUET, format = "parquet") |> filter(Decoy == 0, PG.Q.Value <= 0.01) |> select(Run) |> distinct() |> collect() |> pull(Run)
  removed <- sort(setdiff(raw_runs, kept_runs))
  assert(length(removed) == 3L, paste("Observed", length(removed), "removed runs"))
})

run_check("I/O", "Run order equivalence", function() {
  old <- arrow::read_parquet(PARQUET) |> filter(Decoy == 0, PG.Q.Value <= 0.01)
  new <- arrow::open_dataset(PARQUET, format = "parquet") |> filter(Decoy == 0, PG.Q.Value <= 0.01) |> collect()
  assert(identical(unique(old$Run), unique(new$Run)), "Run order differs")
})

run_check("I/O", "Peptide-evidence equivalence", function() {
  build <- function(report) report |> select(Stripped.Sequence, Protein.Group) |> separate_rows(Protein.Group, sep = ";") |> mutate(Protein.Group = trimws(Protein.Group)) |> filter(!is.na(Protein.Group), nzchar(Protein.Group)) |> distinct() |> arrange(Protein.Group, Stripped.Sequence)
  old <- build(arrow::read_parquet(PARQUET) |> filter(Decoy == 0, PG.Q.Value <= 0.01))
  new <- build(arrow::open_dataset(PARQUET, format = "parquet") |> filter(Decoy == 0, PG.Q.Value <= 0.01) |> collect())
  assert(isTRUE(all.equal(old, new, check.attributes = TRUE)), "Peptide evidence differs")
  assert(n_distinct(old$Protein.Group) == 1872L, "Protein-group baseline changed")
})

run_check("I/O", "FASTA parser baseline", function() {
  fasta <- microseq::readFasta(FASTA)
  assert(nrow(fasta) == 20416L, paste("Observed", nrow(fasta), "FASTA rows"))
})

run_check("I/O", "Single FASTA read implementation", function() {
  scripts <- list.files("R", pattern = "[.]R$", full.names = TRUE)
  occurrences <- sum(vapply(scripts, function(path) sum(grepl("microseq::readFasta", readLines(path, warn = FALSE), fixed = TRUE)), integer(1)))
  assert(occurrences == 1L, paste("Found", occurrences, "readFasta calls"))
})

# -----------------------------------------------------------------------------
# Scientific diagnostics from the technical review
# -----------------------------------------------------------------------------

filtered_report <- arrow::read_parquet(PARQUET) |> filter(Decoy == 0, PG.Q.Value <= 0.01)

run_check("Science", "A1 repeated PG.MaxLFQ formula absent", function() {
  text <- read_script("qc_abundance.R")
  assert(!grepl("PG.MaxLFQ / sum(PG.MaxLFQ", text, fixed = TRUE), "Legacy repeated-PG.MaxLFQ proportion formula is present")
})

run_check("Science", "A2 one PG.MaxLFQ value per group/run", function() {
  multiplicity <- filtered_report |> filter(!is.na(PG.MaxLFQ), PG.MaxLFQ > 0) |> distinct(Run, Protein.Group, PG.MaxLFQ) |> count(Run, Protein.Group) |> filter(n > 1)
  assert(nrow(multiplicity) == 0L, paste(nrow(multiplicity), "Run + Protein.Group pairs contain multiple PG.MaxLFQ values"))
})

# A3 is diagnosed, not silently accepted.
{
  multi_groups <- filtered_report |> distinct(Protein.Group) |> filter(grepl(";", Protein.Group, fixed = TRUE))
  if (nrow(multi_groups) == 0L) {
    record("Science", "A3 multi-protein group annotation", "PASS", "No multi-protein groups in the example")
  } else {
    record("Science", "A3 multi-protein group annotation", "KNOWN ISSUE", paste(nrow(multi_groups), "multi-protein groups require explicit annotation policy"))
  }
}

# A5 is diagnosed by the shared static expression in both modules.
{
  affected <- c("qc_peptides.R", "qc_figures_of_merit.R")
  pattern <- "str_count("
  minus_one <- "- 1"
  present <- vapply(affected, function(file) {
    text <- read_script(file)
    grepl(pattern, text, fixed = TRUE) && grepl(minus_one, text, fixed = TRUE)
  }, logical(1))
  if (any(present)) {
    record("Science", "A5 missed-cleavage terminal rule", "KNOWN ISSUE", paste("Current subtraction rule detected in", paste(affected[present], collapse = ", ")))
  } else {
    record("Science", "A5 missed-cleavage terminal rule", "PASS")
  }
}

run_check("Science", "Completeness dense/aggregate equivalence", function() {
  dense <- filtered_report |> select(Run, Precursor.Id, Precursor.Quantity) |> distinct() |> pivot_wider(names_from = Precursor.Id, values_from = Precursor.Quantity)
  dense_observed <- sum(!is.na(dense[-1]))
  pairs <- filtered_report |> select(Run, Precursor.Id, Precursor.Quantity) |> distinct()
  aggregate_observed <- sum(!is.na(pairs$Precursor.Quantity))
  assert(dense_observed == aggregate_observed, paste("dense=", dense_observed, " aggregate=", aggregate_observed, sep = ""))
})

# Static review risks that are intentionally classified before refactoring.
{
  text <- read_script("qc_peptides.R")
  if (grepl("limits = c(0, 0.5)", text, fixed = TRUE)) {
    record("Visualization", "Missed-cleavage axis truncation", "KNOWN ISSUE", "Fixed 0-50% axis can remove pathological observations")
  } else {
    record("Visualization", "Missed-cleavage axis truncation", "PASS")
  }
}

{
  qmd <- file.path("report", "qc_espresso_scientific.qmd")
  if (file.exists(qmd)) record("Scientific mode", "Scientific QMD present", "PASS") else record("Scientific mode", "Scientific QMD present", "NOT IMPLEMENTED")
}

{
  qmd <- file.path("report", "qc_espresso_scientific.qmd")
  if (!file.exists(qmd)) {
    record("Scientific mode", "End-to-end scientific render", "NOT IMPLEMENTED")
  } else {
    record("Scientific mode", "End-to-end scientific render", "FAIL", "QMD exists but render test has not yet been added")
  }
}

# -----------------------------------------------------------------------------
# Output consolidated audit artifacts
# -----------------------------------------------------------------------------

result_table <- do.call(rbind, results)
output_dir <- file.path("tests", "test-results", "latest")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(result_table, file.path(output_dir, "audit-results.tsv"), sep = "\t", row.names = FALSE, quote = TRUE, na = "")

counts <- table(factor(result_table$Status, levels = c("PASS", "FAIL", "KNOWN ISSUE", "NOT IMPLEMENTED")))
summary_lines <- c(
  "QC-espresso consolidated pre-refactoring audit",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Project:", ROOT),
  paste("R:", R.version.string),
  "",
  paste("PASS:", counts[["PASS"]]),
  paste("FAIL:", counts[["FAIL"]]),
  paste("KNOWN ISSUE:", counts[["KNOWN ISSUE"]]),
  paste("NOT IMPLEMENTED:", counts[["NOT IMPLEMENTED"]]),
  "",
  paste("FINAL STATUS:", if (counts[["FAIL"]] == 0L) "PASS WITH DOCUMENTED GAPS" else "FAIL")
)
writeLines(summary_lines, file.path(output_dir, "audit-summary.txt"))

cat("\n============================================================\n")
cat("AUDIT SUMMARY\n")
cat("============================================================\n")
cat(paste(summary_lines[6:length(summary_lines)], collapse = "\n"), "\n")
cat("============================================================\n")
cat("Detailed results:", normalizePath(file.path(output_dir, "audit-results.tsv"), winslash = "/", mustWork = TRUE), "\n")

quit(status = if (counts[["FAIL"]] == 0L) 0L else 1L)
