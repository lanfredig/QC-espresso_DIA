###############################################################################
# File: tests/test_intensity.R
# Purpose: Validate deduplication, consistency checks, and low-data behavior in
# R/qc_intensity.R.
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

  cat(
    sprintf(
      "%-60s %s\n",
      name,
      if (result$success) "PASS" else "FAIL"
    )
  )

  if (!result$success) {
    cat(
      "  Reason: ",
      result$message,
      "\n",
      sep = ""
    )
  }
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

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

  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}

ROOT <- resolve_root()
setwd(ROOT)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

file_palette <- c(
  A = "#0072B2",
  B = "#E69F00"
)

run_module <- function(report) {
  env <- new.env(parent = globalenv())
  env$report <- report
  env$file_palette <- file_palette
  env$AR_DENSITY <- 1
  env$AR_BOXPLOT <- 1

  messages <- character()
  error_condition <- NULL

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)

  on.exit(
    {
      try(grDevices::dev.off(), silent = TRUE)
      unlink(plot_file, force = TRUE)
    },
    add = TRUE
  )

  tryCatch(
    withCallingHandlers(
      sys.source(
        file.path("R", "qc_intensity.R"),
        envir = env
      ),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_condition <<- e
    }
  )

  list(
    env = env,
    messages = messages,
    error = error_condition
  )
}

record_test("Repeated identical PG.MaxLFQ values collapse to one row", function() {
  report <- tibble(
    Run = factor(c("A", "A", "A"), levels = c("A", "B")),
    Protein.Group = c("P1", "P1", "P1"),
    PG.MaxLFQ = c(100, 100, 100),
    Precursor.Quantity = c(10, 20, 30)
  )

  result <- run_module(report)

  assert_true(is.null(result$error), conditionMessage(result$error))
  assert_true(nrow(result$env$intensity_boxplot) == 1L, "Expected one protein-level row")
})

record_test("Different proteins remain independent", function() {
  report <- tibble(
    Run = factor(c("A", "A"), levels = c("A", "B")),
    Protein.Group = c("P1", "P2"),
    PG.MaxLFQ = c(100, 200),
    Precursor.Quantity = c(10, 20)
  )

  result <- run_module(report)

  assert_true(is.null(result$error), conditionMessage(result$error))
  assert_true(nrow(result$env$intensity_boxplot) == 2L, "Expected two protein-level rows")
})

record_test("Distinct PG.MaxLFQ values for one pair stop execution", function() {
  report <- tibble(
    Run = factor(c("A", "A"), levels = c("A", "B")),
    Protein.Group = c("P1", "P1"),
    PG.MaxLFQ = c(100, 120),
    Precursor.Quantity = c(10, 20)
  )

  result <- run_module(report)

  assert_true(!is.null(result$error), "Expected an inconsistency error")
  assert_true(
    grepl(
      "Multiple distinct positive PG.MaxLFQ values",
      conditionMessage(result$error),
      fixed = TRUE
    ),
    "The error did not explain the PG.MaxLFQ inconsistency"
  )
})

record_test("Non-positive values are excluded", function() {
  report <- tibble(
    Run = factor(c("A", "A", "B"), levels = c("A", "B")),
    Protein.Group = c("P1", "P2", "P3"),
    PG.MaxLFQ = c(100, 0, -10),
    Precursor.Quantity = c(10, 0, -5)
  )

  result <- run_module(report)

  assert_true(is.null(result$error), conditionMessage(result$error))
  assert_true(nrow(result$env$intensity_boxplot) == 1L, "Non-positive PG.MaxLFQ values were retained")
  assert_true(nrow(result$env$prec_density) == 1L, "Non-positive precursor quantities were retained")
})

record_test("No positive precursor quantities skip density plot", function() {
  report <- tibble(
    Run = factor(c("A", "B"), levels = c("A", "B")),
    Protein.Group = c("P1", "P2"),
    PG.MaxLFQ = c(100, 200),
    Precursor.Quantity = c(0, NA_real_)
  )

  result <- run_module(report)

  assert_true(is.null(result$error), conditionMessage(result$error))
  assert_true(is.null(result$env$p_density), "Raw density plot was created without positive quantities")
  assert_true(
    any(grepl("Skipping raw precursor intensity density", result$messages, fixed = TRUE)),
    "The expected raw-intensity skip message was not emitted"
  )
})

record_test("No positive PG.MaxLFQ values skip protein boxplot", function() {
  report <- tibble(
    Run = factor(c("A", "B"), levels = c("A", "B")),
    Protein.Group = c("P1", "P2"),
    PG.MaxLFQ = c(0, NA_real_),
    Precursor.Quantity = c(10, 20)
  )

  result <- run_module(report)

  assert_true(is.null(result$error), conditionMessage(result$error))
  assert_true(is.null(result$env$p_box), "Protein boxplot was created without positive PG.MaxLFQ values")
  assert_true(
    any(grepl("Skipping protein-level LFQ intensity boxplot", result$messages, fixed = TRUE)),
    "The expected protein-intensity skip message was not emitted"
  )
})

record_test("Article example has one distinct value per Run and group", function() {
  suppressPackageStartupMessages(library(arrow))

  report <- arrow::open_dataset(
    file.path("example-data", "dataset_4_contaminans.parquet"),
    format = "parquet"
  ) %>%
    dplyr::filter(
      Decoy == 0,
      PG.Q.Value <= 0.01
    ) %>%
    dplyr::collect()

  inconsistent <- report %>%
    dplyr::filter(
      !is.na(PG.MaxLFQ),
      is.finite(PG.MaxLFQ),
      PG.MaxLFQ > 0
    ) %>%
    dplyr::distinct(
      Run,
      Protein.Group,
      PG.MaxLFQ
    ) %>%
    dplyr::count(
      Run,
      Protein.Group
    ) %>%
    dplyr::filter(n > 1L)

  assert_true(nrow(inconsistent) == 0L, "The article example contains inconsistent PG.MaxLFQ values")
})

passed <- sum(vapply(results, function(x) isTRUE(x$success), logical(1)))
failed <- length(results) - passed

cat("\n============================================================\n")
cat("INTENSITY TEST SUMMARY\n")
cat("============================================================\n")
cat("Tests passed: ", passed, "\n", sep = "")
cat("Tests failed: ", failed, "\n", sep = "")
cat("FINAL STATUS: ", if (failed == 0L) "PASS" else "FAIL", "\n", sep = "")
cat("============================================================\n")

quit(status = if (failed == 0L) 0L else 1L)
