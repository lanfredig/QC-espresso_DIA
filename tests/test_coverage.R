###############################################################################
# File: tests/test_coverage.R
#
# Purpose:
#   Validate protein sequence coverage calculations and low-data visualization
#   behavior implemented in R/qc_coverage.R.
#
# Validation scope:
#   - Validate peptide mapping with absent, single, repeated, and overlapping
#     peptide occurrences.
#   - Confirm that covered residues are not counted more than once.
#   - Confirm that coverage remains between zero and one hundred percent.
#   - Validate independent mapping of shared peptide evidence per accession.
#   - Confirm one coverage row per FASTA accession.
#   - Validate behavior with zero, one, and multiple detected proteins.
#   - Confirm that coverage outputs are generated dynamically when possible.
#
# Result categories:
#   PASS  The expected coverage behavior was observed.
#   FAIL  The coverage behavior was not observed.
#
# Usage:
#   Rscript --no-init-file tests/test_coverage.R
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
      "%-66s %s\n",
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
    library(purrr)
    library(stringr)
    library(tibble)
  }
)


## 4. Isolated module execution ------------------------------------------------

load_coverage_module <- function(
    fasta_ref,
    peptides_clean
) {
  execution_environment <- new.env(
    parent = globalenv()
  )

  execution_environment$fasta_ref <- fasta_ref
  execution_environment$peptides_clean <- peptides_clean
  execution_environment$results_dir <- tempfile(
    "coverage_results_"
  )
  execution_environment$AR_DENSITY <- 1

  dir.create(
    execution_environment$results_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

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
          "qc_coverage.R"
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
    plot_file = plot_file,
    results_dir = execution_environment$results_dir
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

  if (
    !is.null(module_result$results_dir) &&
    dir.exists(module_result$results_dir)
  ) {
    unlink(
      module_result$results_dir,
      recursive = TRUE,
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


## 5. Load the coverage helper -------------------------------------------------

helper_fasta_ref <- tibble::tibble(
  Protein.Group = "HELPER_PROTEIN",
  Protein.Name = "Helper protein",
  Sequence = "PEPTIDE",
  Total_Len = 7L
)

helper_peptides_clean <- tibble::tibble(
  Protein.Group = "HELPER_PROTEIN",
  Peptide_List = list(
    "PEPTIDE"
  )
)

helper_result <- load_coverage_module(
  fasta_ref = helper_fasta_ref,
  peptides_clean = helper_peptides_clean
)


assert_true(
  is.null(
    helper_result$error
  ),
  if (is.null(helper_result$error)) {
    ""
  } else {
    conditionMessage(
      helper_result$error
    )
  }
)


assert_true(
  object_exists(
    helper_result,
    "calculate_coverage"
  ),
  "qc_coverage.R did not define calculate_coverage."
)


calculate_coverage <- get_object(
  helper_result,
  "calculate_coverage"
)


cleanup_module_result(
  helper_result
)


## 6. Controlled coverage calculations ----------------------------------------

record_test(
  "No peptide evidence gives zero coverage",
  function() {
    observed <- calculate_coverage(
      "ABCDEFGHIJ",
      NULL
    )

    assert_equal(
      observed,
      0L,
      "Expected zero covered residues without peptide evidence."
    )
  }
)


record_test(
  "Empty peptide vector gives zero coverage",
  function() {
    observed <- calculate_coverage(
      "ABCDEFGHIJ",
      character()
    )

    assert_equal(
      observed,
      0L,
      "Expected zero covered residues for an empty peptide vector."
    )
  }
)


record_test(
  "Single peptide occurrence is counted",
  function() {
    observed <- calculate_coverage(
      "ABCDEFGHIJ",
      "CDE"
    )

    assert_equal(
      observed,
      3L,
      "Expected three covered residues."
    )
  }
)


record_test(
  "Repeated peptide occurrences are all covered",
  function() {
    observed <- calculate_coverage(
      "ABCABC",
      "ABC"
    )

    assert_equal(
      observed,
      6L,
      "Expected all repeated peptide occurrences to be covered."
    )
  }
)


record_test(
  "Overlapping peptides are not double counted",
  function() {
    observed <- calculate_coverage(
      "ABCDEFGH",
      c(
        "ABCD",
        "CDEF"
      )
    )

    assert_equal(
      observed,
      6L,
      "Expected six unique covered residues."
    )
  }
)


record_test(
  "Duplicate peptide evidence is not double counted",
  function() {
    observed <- calculate_coverage(
      "ABCDEFGHIJ",
      c(
        "CDE",
        "CDE"
      )
    )

    assert_equal(
      observed,
      3L,
      "Duplicate peptide evidence changed the number of covered residues."
    )
  }
)


record_test(
  "Absent peptide does not add coverage",
  function() {
    observed <- calculate_coverage(
      "ABCDEFGHIJ",
      "XYZ"
    )

    assert_equal(
      observed,
      0L,
      "An absent peptide changed the coverage."
    )
  }
)


record_test(
  "Mixed present and absent peptides are handled correctly",
  function() {
    observed <- calculate_coverage(
      "ABCDEFGHIJ",
      c(
        "ABC",
        "XYZ",
        "HIJ"
      )
    )

    assert_equal(
      observed,
      6L,
      "Present and absent peptides were not handled independently."
    )
  }
)


record_test(
  "Complete sequence peptide gives full residue coverage",
  function() {
    sequence <- "ABCDEFGHIJ"

    observed <- calculate_coverage(
      sequence,
      sequence
    )

    assert_equal(
      observed,
      nchar(sequence),
      "A complete-sequence peptide did not cover the complete sequence."
    )
  }
)


record_test(
  "Covered residues never exceed sequence length",
  function() {
    sequence <- "ABCABC"

    observed <- calculate_coverage(
      sequence,
      c(
        "ABC",
        "BCA",
        "CAB",
        "ABC"
      )
    )

    assert_true(
      observed >= 0L &&
        observed <= nchar(sequence),
      paste0(
        "Covered residues are outside the valid sequence range: ",
        observed,
        "."
      )
    )
  }
)


## 7. Coverage table structure -------------------------------------------------

record_test(
  "Coverage stays within zero and one hundred percent",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Protein.Name = c(
        "Protein 1",
        "Protein 2"
      ),
      Sequence = c(
        "ABCDEFGHIJ",
        "KLMNOPQRST"
      ),
      Total_Len = c(
        10L,
        10L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Peptide_List = list(
        c(
          "ABCDE",
          "DEFGH"
        ),
        c(
          "KLM",
          "QRST"
        )
      )
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_results <- get_object(
      module_result,
      "coverage_results"
    )

    assert_true(
      all(
        is.finite(
          coverage_results$Coverage_Pct
        )
      ),
      "One or more coverage percentages are not finite."
    )

    assert_true(
      all(
        coverage_results$Coverage_Pct >= 0 &
          coverage_results$Coverage_Pct <= 100
      ),
      "One or more coverage percentages are outside 0-100%."
    )
  }
)


record_test(
  "Covered residues do not exceed protein length",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Protein.Name = c(
        "Protein 1",
        "Protein 2"
      ),
      Sequence = c(
        "ABCABC",
        "DEFGHI"
      ),
      Total_Len = c(
        6L,
        6L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Peptide_List = list(
        c(
          "ABC",
          "BCA",
          "ABC"
        ),
        c(
          "DEF",
          "FGH"
        )
      )
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_results <- get_object(
      module_result,
      "coverage_results"
    )

    assert_true(
      all(
        coverage_results$Aminoacids_Covered >= 0L
      ),
      "One or more proteins have negative covered-residue counts."
    )

    assert_true(
      all(
        coverage_results$Aminoacids_Covered <=
          coverage_results$Total_Len
      ),
      "One or more covered-residue counts exceed the protein length."
    )
  }
)


record_test(
  "Shared evidence is mapped independently per accession",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "A",
        "B"
      ),
      Protein.Name = c(
        "Protein A",
        "Protein B"
      ),
      Sequence = c(
        "PEPTIDEAAAA",
        "PEPTIDEBBBB"
      ),
      Total_Len = c(
        11L,
        11L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = c(
        "A",
        "B"
      ),
      Peptide_List = list(
        "PEPTIDE",
        "PEPTIDE"
      )
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_results <- get_object(
      module_result,
      "coverage_results"
    )

    assert_equal(
      nrow(coverage_results),
      2L,
      "Expected one coverage row per FASTA accession."
    )

    assert_true(
      !anyDuplicated(
        coverage_results$Protein.Group
      ),
      "Coverage rows were duplicated."
    )

    assert_true(
      all(
        coverage_results$Aminoacids_Covered == 7L
      ),
      "Shared peptide evidence was not mapped independently."
    )
  }
)


record_test(
  "FASTA accessions without evidence remain in coverage table",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2",
        "P3"
      ),
      Protein.Name = c(
        "Protein 1",
        "Protein 2",
        "Protein 3"
      ),
      Sequence = c(
        "AAAA",
        "BBBB",
        "CCCC"
      ),
      Total_Len = c(
        4L,
        4L,
        4L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = "P1",
      Peptide_List = list(
        "AA"
      )
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_results <- get_object(
      module_result,
      "coverage_results"
    )

    assert_equal(
      nrow(coverage_results),
      nrow(fasta_ref),
      "The coverage table does not contain one row per FASTA accession."
    )

    assert_true(
      setequal(
        coverage_results$Protein.Group,
        fasta_ref$Protein.Group
      ),
      "The coverage-table accession inventory differs from the FASTA."
    )

    zero_coverage_groups <- coverage_results |>
      dplyr::filter(
        Protein.Group %in%
          c(
            "P2",
            "P3"
          )
      )

    assert_true(
      all(
        zero_coverage_groups$Aminoacids_Covered == 0L
      ),
      "FASTA accessions without peptide evidence have non-zero coverage."
    )
  }
)


## 8. Low-data visualization behavior -----------------------------------------

record_test(
  "Zero detected proteins skip density plot",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Protein.Name = c(
        "Protein 1",
        "Protein 2"
      ),
      Sequence = c(
        "AAAA",
        "BBBB"
      ),
      Total_Len = c(
        4L,
        4L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = character(),
      Peptide_List = list()
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_detected <- get_object(
      module_result,
      "coverage_results_detected"
    )

    assert_equal(
      nrow(coverage_detected),
      0L,
      "Expected no detected proteins."
    )

    assert_true(
      object_exists(
        module_result,
        "p_coverage"
      ),
      "p_coverage was not defined."
    )

    assert_true(
      is.null(
        get_object(
          module_result,
          "p_coverage"
        )
      ),
      "A density plot was created without detected proteins."
    )

    assert_true(
      message_mentions(
        module_result,
        "fewer than two proteins"
      ),
      "The expected low-data skip message was not emitted."
    )
  }
)


record_test(
  "One detected protein skips density plot",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Protein.Name = c(
        "Protein 1",
        "Protein 2"
      ),
      Sequence = c(
        "AAAA",
        "BBBB"
      ),
      Total_Len = c(
        4L,
        4L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = "P1",
      Peptide_List = list(
        "AA"
      )
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_detected <- get_object(
      module_result,
      "coverage_results_detected"
    )

    assert_equal(
      nrow(coverage_detected),
      1L,
      "Expected one detected protein."
    )

    assert_true(
      is.null(
        get_object(
          module_result,
          "p_coverage"
        )
      ),
      "A density plot was created for one observation."
    )

    assert_true(
      message_mentions(
        module_result,
        "fewer than two proteins"
      ),
      "The expected low-data skip message was not emitted."
    )
  }
)


record_test(
  "Two detected proteins generate density plot",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Protein.Name = c(
        "Protein 1",
        "Protein 2"
      ),
      Sequence = c(
        "AAAA",
        "BBBB"
      ),
      Total_Len = c(
        4L,
        4L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Peptide_List = list(
        "AA",
        "BB"
      )
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_plot <- get_object(
      module_result,
      "p_coverage"
    )

    assert_true(
      inherits(
        coverage_plot,
        "ggplot"
      ),
      "The density plot was not created for two detected proteins."
    )
  }
)


## 9. Output-file behavior -----------------------------------------------------

record_test(
  "Coverage CSV is generated and readable",
  function() {
    fasta_ref <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Protein.Name = c(
        "Protein 1",
        "Protein 2"
      ),
      Sequence = c(
        "AAAA",
        "BBBB"
      ),
      Total_Len = c(
        4L,
        4L
      )
    )

    peptides_clean <- tibble::tibble(
      Protein.Group = c(
        "P1",
        "P2"
      ),
      Peptide_List = list(
        "AA",
        "BB"
      )
    )

    module_result <- load_coverage_module(
      fasta_ref,
      peptides_clean
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

    coverage_files <- list.files(
      module_result$results_dir,
      pattern = "^protein_coverage_with_peptides_.*[.]csv$",
      full.names = TRUE
    )

    assert_true(
      length(coverage_files) > 0L,
      "The coverage module did not generate a coverage CSV."
    )

    coverage_info <- file.info(
      coverage_files
    )

    valid_files <- coverage_files[
      !is.na(
        coverage_info$size
      ) &
        coverage_info$size > 0
    ]

    assert_true(
      length(valid_files) > 0L,
      "The generated coverage CSV is empty or inaccessible."
    )

    coverage_columns <- names(
      utils::read.csv(
        valid_files[[1L]],
        nrows = 1L,
        check.names = FALSE
      )
    )

    required_columns <- c(
      "Protein.Group",
      "Aminoacids_Covered",
      "Total_Len",
      "Coverage_Pct"
    )

    missing_columns <- setdiff(
      required_columns,
      coverage_columns
    )

    assert_true(
      length(missing_columns) == 0L,
      paste0(
        "The generated coverage CSV is missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
)


## 10. Final summary -----------------------------------------------------------

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
cat("COVERAGE TEST SUMMARY\n")
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