###############################################################################
# File: tests/test_multiprotein_abundance.R
#
# Purpose:
#   Validate the multi-protein annotation and abundance policy implemented in
#   R/qc_abundance.R using the available example DIA-NN report and FASTA.
#
# Validation scope:
#   - Preserve one annotation row per original Protein.Group.
#   - Preserve every original Protein.Group after FASTA annotation.
#   - Preserve the number of DIA-NN quantitative rows.
#   - Prevent row duplication during annotation joins.
#   - Preserve every multi-protein group present in the input.
#   - Validate multi-protein annotation statuses dynamically.
#   - Retain unmatched accessions in partial FASTA matches.
#   - Preserve signal for every Run + Protein.Group pair.
#   - Preserve total signal per run.
#   - Preserve abundance proportions summing to one per run.
#   - Recognize contaminants at any position in a protein group.
#
# Result categories:
#   PASS  The expected test behavior was observed.
#   FAIL  The test produced an unexpected result.
#
# Usage:
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
        paste(expected, collapse = ", "),
        "\nObserved: ",
        paste(observed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


object_exists <- function(
    object_name
) {
  exists(
    object_name,
    envir = globalenv(),
    inherits = FALSE
  )
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
    library(tidyverse)
    library(arrow)
    library(microseq)
    library(scales)
  }
)


## 4. Example inputs -----------------------------------------------------------

example_report <- file.path(
  PROJECT_ROOT,
  "example-data",
  "report.parquet"
)

example_fasta <- file.path(
  PROJECT_ROOT,
  "example-data",
  "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta"
)


record_test(
  "Example DIA-NN report is available",
  function() {
    assert_true(
      file.exists(example_report),
      paste0(
        "Example Parquet was not found: ",
        example_report
      )
    )

    report_info <- file.info(
      example_report
    )

    assert_true(
      !is.na(report_info$size) &&
        report_info$size > 0,
      "The example Parquet is empty or inaccessible."
    )
  }
)


record_test(
  "Example FASTA reference is available",
  function() {
    assert_true(
      file.exists(example_fasta),
      paste0(
        "Example FASTA was not found: ",
        example_fasta
      )
    )

    fasta_info <- file.info(
      example_fasta
    )

    assert_true(
      !is.na(fasta_info$size) &&
        fasta_info$size > 0,
      "The example FASTA is empty or inaccessible."
    )
  }
)


## 5. Validated input workflow -------------------------------------------------

test_output_dir <- file.path(
  tempdir(),
  paste0(
    "qc_abundance_multiprotein_test_",
    Sys.getpid()
  )
)

dir.create(
  test_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

on.exit(
  unlink(
    test_output_dir,
    recursive = TRUE,
    force = TRUE
  ),
  add = TRUE
)


params <- list(
  report_path = example_report,
  fasta_path = example_fasta,
  output_dir = test_output_dir
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


## 6. Input requirements -------------------------------------------------------

required_report_columns <- c(
  "Run",
  "Protein.Group",
  "Precursor.Quantity"
)

missing_report_columns <- setdiff(
  required_report_columns,
  names(report)
)

assert_true(
  length(missing_report_columns) == 0L,
  paste0(
    "The loaded report is missing required columns: ",
    paste(
      missing_report_columns,
      collapse = ", "
    )
  )
)


assert_true(
  nrow(report) > 0L,
  "The loaded report contains no retained rows."
)


assert_true(
  nrow(fasta_ref) > 0L,
  "The parsed FASTA reference contains no entries."
)


## 7. Quantitative baseline before annotation ---------------------------------

baseline_report_rows <- nrow(
  report
)


baseline_protein_groups <- sort(
  unique(
    as.character(
      report$Protein.Group[
        !is.na(report$Protein.Group) &
          nzchar(report$Protein.Group)
      ]
    )
  )
)


baseline_signal <- report |>
  dplyr::filter(
    !is.na(Run),
    !is.na(Protein.Group),
    nzchar(as.character(Run)),
    nzchar(Protein.Group),
    !is.na(Precursor.Quantity),
    is.finite(Precursor.Quantity),
    Precursor.Quantity > 0
  ) |>
  dplyr::group_by(
    Run,
    Protein.Group
  ) |>
  dplyr::summarise(
    Signal = sum(
      Precursor.Quantity
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Run = as.character(
      Run
    )
  ) |>
  dplyr::arrange(
    Run,
    Protein.Group
  )


baseline_signal_per_run <- baseline_signal |>
  dplyr::group_by(
    Run
  ) |>
  dplyr::summarise(
    TotalSignal = sum(
      Signal
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    Run
  )


baseline_multi_groups <- report |>
  dplyr::filter(
    !is.na(Protein.Group),
    nzchar(Protein.Group),
    grepl(
      ";",
      Protein.Group,
      fixed = TRUE
    )
  ) |>
  dplyr::distinct(
    Protein.Group
  ) |>
  dplyr::arrange(
    Protein.Group
  )


baseline_single_groups <- report |>
  dplyr::filter(
    !is.na(Protein.Group),
    nzchar(Protein.Group),
    !grepl(
      ";",
      Protein.Group,
      fixed = TRUE
    )
  ) |>
  dplyr::distinct(
    Protein.Group
  ) |>
  dplyr::arrange(
    Protein.Group
  )


## 8. Prevent interactive plot windows ----------------------------------------

legend_fill <- function() {
  ggplot2::theme()
}


pdf_file <- tempfile(
  fileext = ".pdf"
)

grDevices::pdf(
  pdf_file
)

graphics_device_open <- TRUE

on.exit(
  {
    if (isTRUE(graphics_device_open)) {
      try(
        grDevices::dev.off(),
        silent = TRUE
      )
    }

    unlink(
      pdf_file,
      force = TRUE
    )
  },
  add = TRUE
)


## 9. Execute abundance module -------------------------------------------------

source(
  file.path(
    "R",
    "qc_abundance.R"
  )
)


## 10. Required module outputs -------------------------------------------------

record_test(
  "Annotation table was created",
  function() {
    assert_true(
      object_exists(
        "protein_group_annotation"
      ),
      "protein_group_annotation was not created."
    )

    assert_true(
      is.data.frame(
        protein_group_annotation
      ),
      "protein_group_annotation is not a data frame."
    )

    assert_true(
      nrow(protein_group_annotation) > 0L,
      "protein_group_annotation is empty."
    )
  }
)


record_test(
  "Labeled report was created",
  function() {
    assert_true(
      object_exists(
        "report_labeled"
      ),
      "report_labeled was not created."
    )

    assert_true(
      is.data.frame(
        report_labeled
      ),
      "report_labeled is not a data frame."
    )
  }
)


record_test(
  "Protein-group signal table was created",
  function() {
    assert_true(
      object_exists(
        "protein_group_signal"
      ),
      "protein_group_signal was not created."
    )

    assert_true(
      is.data.frame(
        protein_group_signal
      ),
      "protein_group_signal is not a data frame."
    )
  }
)


record_test(
  "Abundance table was created",
  function() {
    assert_true(
      object_exists(
        "abundance"
      ),
      "abundance was not created."
    )

    assert_true(
      is.data.frame(
        abundance
      ),
      "abundance is not a data frame."
    )
  }
)


record_test(
  "Contaminant recognition pattern was created",
  function() {
    assert_true(
      object_exists(
        "recognized_contaminant_pattern"
      ),
      "recognized_contaminant_pattern was not created."
    )

    assert_true(
      is.character(
        recognized_contaminant_pattern
      ) &&
        length(
          recognized_contaminant_pattern
        ) == 1L &&
        nzchar(
          recognized_contaminant_pattern
        ),
      "recognized_contaminant_pattern is missing or invalid."
    )
  }
)


## 11. Annotation structure ----------------------------------------------------

record_test(
  "Annotation table contains required columns",
  function() {
    required_columns <- c(
      "Protein.Group",
      "Accession_Count",
      "Annotation_Status",
      "Gene_Label",
      "Unmatched_Accessions"
    )

    missing_columns <- setdiff(
      required_columns,
      names(protein_group_annotation)
    )

    assert_true(
      length(missing_columns) == 0L,
      paste0(
        "protein_group_annotation is missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
)


record_test(
  "One annotation row exists per Protein.Group",
  function() {
    duplicated_groups <- protein_group_annotation |>
      dplyr::count(
        Protein.Group,
        name = "Rows"
      ) |>
      dplyr::filter(
        Rows != 1L
      )

    assert_true(
      nrow(duplicated_groups) == 0L,
      paste0(
        nrow(duplicated_groups),
        " protein groups have an annotation-row count different from one."
      )
    )
  }
)


record_test(
  "Every original Protein.Group was preserved",
  function() {
    annotated_groups <- sort(
      unique(
        as.character(
          protein_group_annotation$Protein.Group
        )
      )
    )

    assert_true(
      identical(
        baseline_protein_groups,
        annotated_groups
      ),
      paste0(
        "The annotated Protein.Group inventory differs from ",
        "the original report inventory."
      )
    )
  }
)


record_test(
  "Accession counts match Protein.Group identifiers",
  function() {
    expected_accession_count <- stringr::str_count(
      protein_group_annotation$Protein.Group,
      stringr::fixed(";")
    ) + 1L

    invalid_rows <- which(
      is.na(
        protein_group_annotation$Accession_Count
      ) |
        protein_group_annotation$Accession_Count !=
        expected_accession_count
    )

    assert_true(
      length(invalid_rows) == 0L,
      paste0(
        length(invalid_rows),
        " protein groups have an Accession_Count inconsistent ",
        "with the semicolon-delimited identifier."
      )
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
      "The annotation join duplicated one or more report rows."
    )
  }
)


## 12. Protein-group inventory -------------------------------------------------

record_test(
  "Every input multi-protein group was preserved",
  function() {
    observed_multi_groups <- protein_group_annotation |>
      dplyr::filter(
        Accession_Count > 1L
      ) |>
      dplyr::distinct(
        Protein.Group
      ) |>
      dplyr::arrange(
        Protein.Group
      )

    assert_equal(
      nrow(observed_multi_groups),
      nrow(baseline_multi_groups),
      "The number of multi-protein groups changed after annotation."
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
  "Single-protein and multi-protein inventories are consistent",
  function() {
    observed_single_groups <- protein_group_annotation |>
      dplyr::filter(
        Accession_Count == 1L
      ) |>
      dplyr::distinct(
        Protein.Group
      )

    observed_multi_groups <- protein_group_annotation |>
      dplyr::filter(
        Accession_Count > 1L
      ) |>
      dplyr::distinct(
        Protein.Group
      )

    assert_equal(
      nrow(observed_single_groups),
      nrow(baseline_single_groups),
      "The number of single-protein groups changed after annotation."
    )

    assert_true(
      setequal(
        observed_single_groups$Protein.Group,
        baseline_single_groups$Protein.Group
      ),
      "The single-protein group inventory changed after annotation."
    )

    assert_equal(
      nrow(observed_single_groups) +
        nrow(observed_multi_groups),
      length(baseline_protein_groups),
      paste0(
        "The annotated single-protein and multi-protein inventories ",
        "do not partition the original Protein.Group inventory."
      )
    )
  }
)


## 13. Annotation-status validation -------------------------------------------

multi_protein_statuses <- c(
  "CONSISTENT_ANNOTATION",
  "MULTIPLE_ANNOTATIONS",
  "PARTIAL_FASTA_MATCH"
)


record_test(
  "Multi-protein statuses are non-missing and recognized",
  function() {
    invalid_multi_groups <- protein_group_annotation |>
      dplyr::filter(
        Accession_Count > 1L
      ) |>
      dplyr::filter(
        is.na(Annotation_Status) |
          !nzchar(
            trimws(
              Annotation_Status
            )
          ) |
          !Annotation_Status %in%
          multi_protein_statuses
      )

    assert_true(
      nrow(invalid_multi_groups) == 0L,
      paste0(
        nrow(invalid_multi_groups),
        " multi-protein groups have a missing or unsupported ",
        "Annotation_Status."
      )
    )
  }
)


record_test(
  "Single-protein statuses are internally consistent",
  function() {
    single_group_annotation <- protein_group_annotation |>
      dplyr::filter(
        Accession_Count == 1L
      )

    populated_single_statuses <- single_group_annotation |>
      dplyr::filter(
        !is.na(Annotation_Status),
        nzchar(
          trimws(
            Annotation_Status
          )
        )
      )

    invalid_single_statuses <- populated_single_statuses |>
      dplyr::filter(
        Annotation_Status %in%
          multi_protein_statuses
      )

    assert_true(
      nrow(invalid_single_statuses) == 0L,
      paste0(
        nrow(invalid_single_statuses),
        " single-protein groups use a status reserved for ",
        "multi-protein annotation."
      )
    )
  }
)


record_test(
  "Multi-protein classification accounts for every group",
  function() {
    multi_status_summary <- protein_group_annotation |>
      dplyr::filter(
        Accession_Count > 1L
      ) |>
      dplyr::count(
        Annotation_Status,
        name = "Groups"
      )

    classified_multi_groups <- sum(
      multi_status_summary$Groups
    )

    assert_equal(
      classified_multi_groups,
      nrow(baseline_multi_groups),
      paste0(
        "The annotation-status counts do not account for every ",
        "multi-protein group."
      )
    )
  }
)


record_test(
  "Partial FASTA matches retain unmatched accessions",
  function() {
    invalid_partial_groups <- protein_group_annotation |>
      dplyr::filter(
        Annotation_Status ==
          "PARTIAL_FASTA_MATCH"
      ) |>
      dplyr::filter(
        is.na(Unmatched_Accessions) |
          !nzchar(
            trimws(
              Unmatched_Accessions
            )
          )
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
    invalid_labels <- protein_group_annotation |>
      dplyr::filter(
        Annotation_Status ==
          "PARTIAL_FASTA_MATCH"
      ) |>
      dplyr::filter(
        is.na(Gene_Label) |
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
        " partial-match groups lack an explicit [partial] label."
      )
    )
  }
)


record_test(
  "Multi-protein labels are non-empty",
  function() {
    invalid_labels <- protein_group_annotation |>
      dplyr::filter(
        Accession_Count > 1L
      ) |>
      dplyr::filter(
        is.na(Gene_Label) |
          !nzchar(
            trimws(
              Gene_Label
            )
          )
      )

    assert_true(
      nrow(invalid_labels) == 0L,
      paste0(
        nrow(invalid_labels),
        " multi-protein groups have empty gene labels."
      )
    )
  }
)


## 14. Quantitative preservation ----------------------------------------------

record_test(
  "One signal value exists per Run + Protein.Group",
  function() {
    duplicated_signal_rows <- protein_group_signal |>
      dplyr::count(
        Run,
        Protein.Group,
        name = "Rows"
      ) |>
      dplyr::filter(
        Rows != 1L
      )

    assert_true(
      nrow(duplicated_signal_rows) == 0L,
      paste0(
        nrow(duplicated_signal_rows),
        " Run + Protein.Group pairs have a signal-row count ",
        "different from one."
      )
    )
  }
)


record_test(
  "Signal inventory matches the positive quantitative baseline",
  function() {
    observed_pairs <- protein_group_signal |>
      dplyr::transmute(
        Run = as.character(
          Run
        ),
        Protein.Group
      ) |>
      dplyr::arrange(
        Run,
        Protein.Group
      )

    expected_pairs <- baseline_signal |>
      dplyr::select(
        Run,
        Protein.Group
      ) |>
      dplyr::arrange(
        Run,
        Protein.Group
      )

    assert_true(
      identical(
        observed_pairs,
        expected_pairs
      ),
      paste0(
        "The Run + Protein.Group signal inventory differs from ",
        "the positive quantitative baseline."
      )
    )
  }
)


record_test(
  "Signal is unchanged for every Run + Protein.Group",
  function() {
    observed_signal <- protein_group_signal |>
      dplyr::transmute(
        Run = as.character(
          Run
        ),
        Protein.Group,
        Signal
      ) |>
      dplyr::arrange(
        Run,
        Protein.Group
      )

    comparison <- baseline_signal |>
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
      )

    assert_equal(
      nrow(comparison),
      nrow(baseline_signal),
      "The quantitative Run + Protein.Group inventory changed."
    )

    assert_true(
      all(
        !is.na(
          comparison$Signal_Before
        ) &
          !is.na(
            comparison$Signal_After
          )
      ),
      "One or more quantitative pairs are missing before or after annotation."
    )

    tolerance <- sqrt(
      .Machine$double.eps
    )

    relative_difference <- abs(
      comparison$Signal_After -
        comparison$Signal_Before
    ) /
      pmax(
        1,
        abs(
          comparison$Signal_Before
        ),
        abs(
          comparison$Signal_After
        )
      )

    assert_true(
      all(
        relative_difference <=
          tolerance
      ),
      "Signal changed for one or more Run + Protein.Group pairs."
    )
  }
)


record_test(
  "Total signal per run is unchanged",
  function() {
    observed_total <- protein_group_signal |>
      dplyr::mutate(
        Run = as.character(
          Run
        )
      ) |>
      dplyr::group_by(
        Run
      ) |>
      dplyr::summarise(
        TotalSignal = sum(
          Signal
        ),
        .groups = "drop"
      ) |>
      dplyr::arrange(
        Run
      )

    comparison <- baseline_signal_per_run |>
      dplyr::full_join(
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
      "One or more runs were lost or added during total-signal validation."
    )

    assert_true(
      all(
        !is.na(
          comparison$TotalSignal_Before
        ) &
          !is.na(
            comparison$TotalSignal_After
          )
      ),
      "One or more runs are absent before or after abundance annotation."
    )

    tolerance <- sqrt(
      .Machine$double.eps
    )

    relative_difference <- abs(
      comparison$TotalSignal_After -
        comparison$TotalSignal_Before
    ) /
      pmax(
        1,
        abs(
          comparison$TotalSignal_Before
        ),
        abs(
          comparison$TotalSignal_After
        )
      )

    assert_true(
      all(
        relative_difference <=
          tolerance
      ),
      "Total signal changed for one or more runs."
    )
  }
)


record_test(
  "Abundance table covers every quantified run",
  function() {
    abundance_runs <- sort(
      unique(
        as.character(
          abundance$Run
        )
      )
    )

    expected_runs <- sort(
      unique(
        baseline_signal$Run
      )
    )

    assert_true(
      identical(
        abundance_runs,
        expected_runs
      ),
      paste0(
        "The abundance-table run inventory does not match ",
        "the quantified baseline."
      )
    )
  }
)


record_test(
  "Abundance proportions are finite and bounded",
  function() {
    assert_true(
      "Proportion" %in%
        names(abundance),
      "The abundance table does not contain the Proportion column."
    )

    assert_true(
      all(
        is.finite(
          abundance$Proportion
        )
      ),
      "One or more abundance proportions are not finite."
    )

    tolerance <- sqrt(
      .Machine$double.eps
    )

    assert_true(
      all(
        abundance$Proportion >=
          -tolerance &
          abundance$Proportion <=
          1 + tolerance
      ),
      "One or more abundance proportions are outside 0-1."
    )
  }
)


record_test(
  "Abundance proportions sum to one per quantified run",
  function() {
    proportion_summary <- abundance |>
      dplyr::group_by(
        Run
      ) |>
      dplyr::summarise(
        Proportion_Sum = sum(
          Proportion
        ),
        .groups = "drop"
      )

    assert_true(
      nrow(proportion_summary) > 0L,
      "No per-run abundance proportions were available."
    )

    assert_true(
      all(
        abs(
          proportion_summary$Proportion_Sum -
            1
        ) <=
          1e-12
      ),
      "Abundance proportions do not sum to one for every quantified run."
    )
  }
)


## 15. Contaminant recognition -------------------------------------------------

record_test(
  "Contaminants are recognized at any group position",
  function() {
    test_groups <- c(
      "CON_TEST;P12345",
      "P12345;CON_TEST",
      "cRAP-TEST;P12345",
      "P12345;cRAP-TEST",
      "CON_A;P12345;Q67890",
      "P12345;CON_A;Q67890",
      "P12345;Q67890;CON_A",
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
      paste0(
        "Contaminant recognition failed for one or more ",
        "positions in a protein group."
      )
    )
  }
)


record_test(
  "Ordinary protein groups are not classified as contaminants",
  function() {
    ordinary_groups <- c(
      "P12345",
      "P12345;Q67890",
      "Q9TEST;A0A000",
      "PROTEIN_CONSERVED"
    )

    observed <- grepl(
      recognized_contaminant_pattern,
      ordinary_groups,
      perl = TRUE
    )

    assert_true(
      !any(observed),
      "One or more ordinary protein groups were classified as contaminants."
    )
  }
)


## 16. Dynamic annotation summary ---------------------------------------------

multi_annotation_summary <- protein_group_annotation |>
  dplyr::filter(
    Accession_Count > 1L
  ) |>
  dplyr::count(
    Annotation_Status,
    name = "Groups"
  ) |>
  dplyr::arrange(
    Annotation_Status
  )


single_status_summary <- protein_group_annotation |>
  dplyr::filter(
    Accession_Count == 1L
  ) |>
  dplyr::mutate(
    Annotation_Status = dplyr::if_else(
      is.na(Annotation_Status) |
        !nzchar(
          trimws(
            Annotation_Status
          )
        ),
      "NOT_APPLICABLE",
      Annotation_Status
    )
  ) |>
  dplyr::count(
    Annotation_Status,
    name = "Groups"
  ) |>
  dplyr::arrange(
    Annotation_Status
  )


cat("\n")
cat("============================================================\n")
cat("MULTI-PROTEIN ANNOTATION SUMMARY\n")
cat("============================================================\n")
cat(
  "Original Protein.Group values: ",
  length(
    baseline_protein_groups
  ),
  "\n",
  sep = ""
)
cat(
  "Single-protein groups: ",
  nrow(
    baseline_single_groups
  ),
  "\n",
  sep = ""
)
cat(
  "Multi-protein groups: ",
  nrow(
    baseline_multi_groups
  ),
  "\n",
  sep = ""
)

if (nrow(multi_annotation_summary) > 0L) {
  cat("Multi-protein statuses:\n")

  for (
    summary_index in
    seq_len(
      nrow(
        multi_annotation_summary
      )
    )
  ) {
    cat(
      "  ",
      multi_annotation_summary$Annotation_Status[
        summary_index
      ],
      ": ",
      multi_annotation_summary$Groups[
        summary_index
      ],
      "\n",
      sep = ""
    )
  }
}

if (nrow(single_status_summary) > 0L) {
  cat("Single-protein status values:\n")

  for (
    summary_index in
    seq_len(
      nrow(
        single_status_summary
      )
    )
  ) {
    cat(
      "  ",
      single_status_summary$Annotation_Status[
        summary_index
      ],
      ": ",
      single_status_summary$Groups[
        summary_index
      ],
      "\n",
      sep = ""
    )
  }
}

cat(
  "Quantified Run + Protein.Group pairs: ",
  nrow(
    baseline_signal
  ),
  "\n",
  sep = ""
)
cat(
  "Quantified runs: ",
  nrow(
    baseline_signal_per_run
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")


## 17. Graphics-device completion ---------------------------------------------

if (isTRUE(graphics_device_open)) {
  grDevices::dev.off()
  graphics_device_open <- FALSE
}


## 18. Final summary -----------------------------------------------------------

passed_tests <- sum(
  vapply(
    results,
    function(result) {
      isTRUE(result$success)
    },
    logical(1)
  )
)

failed_tests <- length(results) -
  passed_tests


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