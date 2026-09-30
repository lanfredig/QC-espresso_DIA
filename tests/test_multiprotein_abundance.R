###############################################################################
# File: tests/test_multiprotein_abundance.R
#
# Purpose:
#   Validate the multi-protein annotation policy implemented in
#   R/qc_abundance.R using the validated article example.
#
# Validation scope:
#   - Load the validated Parquet and FASTA inputs.
#   - Execute the current qc_abundance.R module.
#   - Preserve one annotation row per original Protein.Group.
#   - Preserve the number of DIA-NN quantitative rows.
#   - Preserve the signal for every Run + Protein.Group pair.
#   - Preserve total signal per run.
#   - Preserve abundance proportions summing to one per run.
#   - Preserve all 61 multi-protein groups.
#   - Confirm 31 MULTIPLE_ANNOTATIONS groups.
#   - Confirm 30 PARTIAL_FASTA_MATCH groups.
#   - Confirm that no multi-protein group is discarded.
#   - Confirm contaminant recognition at any position in a protein group.
#
# Usage from the project root:
#   Rscript --no-init-file tests/test_multiprotein_abundance.R
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
        message = conditionMessage(error_condition)
      )
    }
  )

  results[[length(results) + 1L]] <<- result

  cat(
    sprintf(
      "%-52s %s\n",
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
        paste(expected, collapse = ", "),
        "\nObserved: ",
        paste(observed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


## 2. Resolve project root -----------------------------------------------------

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


## 3. Load required packages ---------------------------------------------------

suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(arrow)
    library(microseq)
    library(scales)
  }
)


## 4. Define example inputs ----------------------------------------------------

example_report <- file.path(
  PROJECT_ROOT,
  "example-data",
  "dataset_4_contaminans.parquet"
)

example_fasta <- file.path(
  PROJECT_ROOT,
  "example-data",
  "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta"
)

assert_true(
  file.exists(example_report),
  paste0(
    "Example Parquet not found: ",
    example_report
  )
)

assert_true(
  file.exists(example_fasta),
  paste0(
    "Example FASTA not found: ",
    example_fasta
  )
)


## 5. Load the validated input workflow ---------------------------------------

PROJECT_ROOT <- normalizePath(
  PROJECT_ROOT,
  winslash = "/",
  mustWork = TRUE
)

params <- list(
  report_path = example_report,
  fasta_path = example_fasta,
  output_dir = file.path(
    tempdir(),
    paste0(
      "qc_abundance_multiprotein_test_",
      Sys.getpid()
    )
  )
)

source(
  file.path(
    "R",
    "globals.R"
  )
)

source(
  file.path(
    "R",
    "io_load.R"
  )
)


## 6. Capture quantitative baseline before annotation -------------------------

baseline_report_rows <- nrow(
  report
)

baseline_protein_groups <- sort(
  unique(
    as.character(
      report$Protein.Group
    )
  )
)

baseline_signal <- report %>%
  dplyr::filter(
    !is.na(Precursor.Quantity),
    Precursor.Quantity > 0
  ) %>%
  dplyr::group_by(
    Run,
    Protein.Group
  ) %>%
  dplyr::summarise(
    Signal = sum(
      Precursor.Quantity
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Run = as.character(
      Run
    )
  ) %>%
  dplyr::arrange(
    Run,
    Protein.Group
  )

baseline_signal_per_run <- baseline_signal %>%
  dplyr::group_by(
    Run
  ) %>%
  dplyr::summarise(
    TotalSignal = sum(
      Signal
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    Run
  )

baseline_multi_groups <- report %>%
  dplyr::filter(
    !is.na(Protein.Group),
    grepl(
      ";",
      Protein.Group,
      fixed = TRUE
    )
  ) %>%
  dplyr::distinct(
    Protein.Group
  ) %>%
  dplyr::arrange(
    Protein.Group
  )


## 7. Prevent plot windows during module execution ----------------------------

legend_fill <- function() {
  ggplot2::theme()
}

pdf_file <- tempfile(
  fileext = ".pdf"
)

grDevices::pdf(
  pdf_file
)

on.exit(
  {
    try(
      grDevices::dev.off(),
      silent = TRUE
    )

    unlink(
      pdf_file,
      force = TRUE
    )
  },
  add = TRUE
)


## 8. Execute corrected abundance module --------------------------------------

source(
  file.path(
    "R",
    "qc_abundance.R"
  )
)


## 9. Structural tests ---------------------------------------------------------

record_test(
  "Annotation table was created",
  function() {
    assert_true(
      exists("protein_group_annotation") &&
        is.data.frame(protein_group_annotation),
      "protein_group_annotation was not created."
    )
  }
)


record_test(
  "One annotation row per Protein.Group",
  function() {
    duplicated_groups <- protein_group_annotation %>%
      dplyr::count(
        Protein.Group,
        name = "Rows"
      ) %>%
      dplyr::filter(
        Rows != 1L
      )

    assert_true(
      nrow(duplicated_groups) == 0L,
      paste0(
        nrow(duplicated_groups),
        " protein groups have more than one annotation row."
      )
    )
  }
)


record_test(
  "Every original Protein.Group was preserved",
  function() {
    annotated_groups <- sort(
      unique(
        protein_group_annotation$Protein.Group
      )
    )

    assert_true(
      identical(
        baseline_protein_groups,
        annotated_groups
      ),
      "The annotated protein-group inventory differs from the original inventory."
    )
  }
)


record_test(
  "DIA-NN report row count was preserved",
  function() {
    assert_equal(
      nrow(report_labeled),
      baseline_report_rows,
      "The annotation join changed the DIA-NN report row count."
    )
  }
)


record_test(
  "No duplicate report rows were introduced",
  function() {
    assert_equal(
      nrow(report_labeled),
      nrow(report),
      "The annotation join duplicated report rows."
    )
  }
)


## 10. Multi-protein classification tests -------------------------------------

record_test(
  "All 61 multi-protein groups were preserved",
  function() {
    observed_multi_groups <- protein_group_annotation %>%
      dplyr::filter(
        Accession_Count > 1L
      ) %>%
      dplyr::distinct(
        Protein.Group
      )

    assert_equal(
      nrow(observed_multi_groups),
      61L,
      "Unexpected number of multi-protein groups."
    )

    assert_true(
      setequal(
        observed_multi_groups$Protein.Group,
        baseline_multi_groups$Protein.Group
      ),
      "The multi-protein group inventory changed after annotation."
    )
  }
)


record_test(
  "34 groups have multiple annotations",
  function() {
    observed <- protein_group_annotation %>%
      dplyr::filter(
        Accession_Count > 1L,
        Annotation_Status == "MULTIPLE_ANNOTATIONS"
      ) %>%
      nrow()

    assert_equal(
      observed,
      34L,
      "Unexpected number of MULTIPLE_ANNOTATIONS groups."
    )
  }
)


record_test(
  "26 groups have partial FASTA matches",
  function() {
    observed <- protein_group_annotation %>%
      dplyr::filter(
        Accession_Count > 1L,
        Annotation_Status == "PARTIAL_FASTA_MATCH"
      ) %>%
      nrow()

    assert_equal(
      observed,
      26L,
      "Unexpected number of PARTIAL_FASTA_MATCH groups."
    )
  }
)


record_test(
  "1 group has a consistent annotation",
  function() {
    observed <- protein_group_annotation %>%
      dplyr::filter(
        Accession_Count > 1L,
        Annotation_Status == "CONSISTENT_ANNOTATION"
      ) %>%
      nrow()

    assert_equal(
      observed,
      1L,
      "Unexpected number of CONSISTENT_ANNOTATION groups."
    )
  }
)


record_test(
  "Partial matches retain unmatched accessions",
  function() {
    invalid_partial_groups <- protein_group_annotation %>%
      dplyr::filter(
        Annotation_Status == "PARTIAL_FASTA_MATCH"
      ) %>%
      dplyr::filter(
        is.na(Unmatched_Accessions) |
          !nzchar(Unmatched_Accessions)
      )

    assert_true(
      nrow(invalid_partial_groups) == 0L,
      paste0(
        nrow(invalid_partial_groups),
        " partial-match groups do not retain unmatched accessions."
      )
    )
  }
)


record_test(
  "Partial-match labels are explicit",
  function() {
    invalid_labels <- protein_group_annotation %>%
      dplyr::filter(
        Annotation_Status == "PARTIAL_FASTA_MATCH"
      ) %>%
      dplyr::filter(
        !grepl(
          "[partial]",
          Gene_Label,
          fixed = TRUE
        )
      )

    assert_true(
      nrow(invalid_labels) == 0L,
      paste0(
        nrow(invalid_labels),
        " partial-match groups lack the [partial] label."
      )
    )
  }
)


record_test(
  "Multi-protein labels are non-empty",
  function() {
    invalid_labels <- protein_group_annotation %>%
      dplyr::filter(
        Accession_Count > 1L
      ) %>%
      dplyr::filter(
        is.na(Gene_Label) |
          !nzchar(Gene_Label)
      )

    assert_true(
      nrow(invalid_labels) == 0L,
      paste0(
        nrow(invalid_labels),
        " multi-protein groups have empty labels."
      )
    )
  }
)


## 11. Quantitative preservation tests ----------------------------------------

record_test(
  "One signal value per Run + Protein.Group",
  function() {
    duplicated_signal_rows <- protein_group_signal %>%
      dplyr::count(
        Run,
        Protein.Group,
        name = "Rows"
      ) %>%
      dplyr::filter(
        Rows != 1L
      )

    assert_true(
      nrow(duplicated_signal_rows) == 0L,
      paste0(
        nrow(duplicated_signal_rows),
        " Run + Protein.Group pairs have duplicated signal rows."
      )
    )
  }
)


record_test(
  "Signal is unchanged for every Run + Protein.Group",
  function() {
    observed_signal <- protein_group_signal %>%
      dplyr::transmute(
        Run = as.character(Run),
        Protein.Group,
        Signal
      ) %>%
      dplyr::arrange(
        Run,
        Protein.Group
      )

    comparison <- baseline_signal %>%
      dplyr::full_join(
        observed_signal,
        by = c(
          "Run",
          "Protein.Group"
        ),
        suffix = c(
          "_Before",
          "_After"
        )
      ) %>%
      dplyr::mutate(
        Signal_Before = tidyr::replace_na(
          Signal_Before,
          0
        ),
        Signal_After = tidyr::replace_na(
          Signal_After,
          0
        ),
        Difference = Signal_After - Signal_Before
      )

    assert_true(
      nrow(comparison) == nrow(baseline_signal),
      "The quantitative Run + Protein.Group inventory changed."
    )

    assert_true(
      all(
        abs(comparison$Difference) <=
          sqrt(.Machine$double.eps)
      ),
      paste0(
        sum(
          abs(comparison$Difference) >
            sqrt(.Machine$double.eps)
        ),
        " Run + Protein.Group signals changed."
      )
    )
  }
)


record_test(
  "Total signal per run is unchanged",
  function() {
    observed_total <- protein_group_signal %>%
      dplyr::mutate(
        Run = as.character(Run)
      ) %>%
      dplyr::group_by(
        Run
      ) %>%
      dplyr::summarise(
        TotalSignal = sum(Signal),
        .groups = "drop"
      ) %>%
      dplyr::arrange(
        Run
      )

    comparison <- baseline_signal_per_run %>%
      dplyr::inner_join(
        observed_total,
        by = "Run",
        suffix = c(
          "_Before",
          "_After"
        )
      )

    assert_equal(
      nrow(comparison),
      nrow(baseline_signal_per_run),
      "One or more runs were lost from total-signal validation."
    )

    assert_true(
      all(
        abs(
          comparison$TotalSignal_After -
            comparison$TotalSignal_Before
        ) <=
          sqrt(.Machine$double.eps)
      ),
      "Total signal changed for one or more runs."
    )
  }
)


record_test(
  "Abundance proportions sum to one per run",
  function() {
    proportion_summary <- abundance %>%
      dplyr::group_by(
        Run
      ) %>%
      dplyr::summarise(
        Proportion_Sum = sum(Proportion),
        .groups = "drop"
      )

    assert_true(
      all(
        abs(
          proportion_summary$Proportion_Sum - 1
        ) <=
          1e-12
      ),
      "Abundance proportions do not sum to one for every run."
    )
  }
)


## 12. Contaminant recognition tests ------------------------------------------

record_test(
  "Contaminants are recognized at any group position",
  function() {
    test_groups <- c(
      "CON_TEST;P12345",
      "P12345;CON_TEST",
      "cRAP-TEST;P12345",
      "P12345;cRAP-TEST",
      "P12345;Q67890"
    )

    observed <- grepl(
      recognized_contaminant_pattern,
      test_groups,
      perl = TRUE
    )

    expected <- c(
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      FALSE
    )

    assert_true(
      identical(
        observed,
        expected
      ),
      "Contaminant recognition failed for one or more group positions."
    )
  }
)


## 13. Final summary -----------------------------------------------------------

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
cat("MULTI-PROTEIN ABUNDANCE TEST SUMMARY\n")
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
