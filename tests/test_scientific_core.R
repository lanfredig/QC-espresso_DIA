###############################################################################
# File: tests/test_scientific_core.R
#
# Purpose:
#   Execute the QC-espresso scientific modules in their intended sequence and
#   validate the integration between the loaded inputs, analytical modules,
#   derived objects, and generated scientific outputs.
#
# Validation scope:
#   - Load a valid example DIA-NN report and FASTA reference.
#   - Execute each scientific module once in the intended sequence.
#   - Confirm that shared input objects remain available after execution.
#   - Validate dynamic relationships between runs and per-run outputs.
#   - Confirm that completeness, coverage, and Figures of Merit outputs exist.
#   - Confirm that generated scientific output files are readable and non-empty.
#   - Avoid executing focused test suites already managed by run_audit.R.
#
# Result categories:
#   PASS  The module or integration invariant completed successfully.
#   FAIL  The module failed or an integration invariant was not satisfied.
#
# Usage:
#   Rscript --no-init-file tests/test_scientific_core.R
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


record_result <- function(
    name,
    success,
    details = ""
) {
  result <- data.frame(
    Test = name,
    Success = isTRUE(success),
    Details = details,
    stringsAsFactors = FALSE
  )

  results[[length(results) + 1L]] <<- result

  cat(
    sprintf(
      "%-58s %s\n",
      name,
      if (isTRUE(success)) {
        "PASS"
      } else {
        "FAIL"
      }
    )
  )

  if (!isTRUE(success) && nzchar(details)) {
    cat(
      "  Reason: ",
      details,
      "\n",
      sep = ""
    )
  }

  invisible(result)
}


run_check <- function(
    name,
    test_function
) {
  outcome <- tryCatch(
    {
      test_function()

      list(
        success = TRUE,
        details = ""
      )
    },
    error = function(error_condition) {
      list(
        success = FALSE,
        details = conditionMessage(error_condition)
      )
    }
  )

  record_result(
    name = name,
    success = outcome$success,
    details = outcome$details
  )

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


## 3. Required packages --------------------------------------------------------

suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(arrow)
    library(microseq)
    library(scales)
    library(ggrepel)
    library(RColorBrewer)
  }
)


## 4. Example inputs and isolated output --------------------------------------

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


run_check(
  "Example DIA-NN report is available",
  function() {
    assert_true(
      file.exists(example_report),
      paste0(
        "Example DIA-NN report was not found: ",
        example_report
      )
    )

    report_info <- file.info(
      example_report
    )

    assert_true(
      !is.na(report_info$size) &&
        report_info$size > 0,
      "The example DIA-NN report is empty or inaccessible."
    )
  }
)


run_check(
  "Example FASTA reference is available",
  function() {
    assert_true(
      file.exists(example_fasta),
      paste0(
        "Example FASTA reference was not found: ",
        example_fasta
      )
    )

    fasta_info <- file.info(
      example_fasta
    )

    assert_true(
      !is.na(fasta_info$size) &&
        fasta_info$size > 0,
      "The example FASTA reference is empty or inaccessible."
    )
  }
)


integration_output_dir <- file.path(
  tempdir(),
  paste0(
    "qc_espresso_integrated_",
    Sys.getpid()
  )
)

dir.create(
  integration_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

on.exit(
  unlink(
    integration_output_dir,
    recursive = TRUE,
    force = TRUE
  ),
  add = TRUE
)


params <- list(
  report_path = example_report,
  fasta_path = example_fasta,
  output_dir = integration_output_dir
)


## 5. Graphics device ----------------------------------------------------------

plot_file <- file.path(
  integration_output_dir,
  "integrated_plots.pdf"
)

grDevices::pdf(
  plot_file,
  width = 12,
  height = 8
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
  },
  add = TRUE
)


## 6. Shared input loading -----------------------------------------------------

input_loading_success <- TRUE


run_check(
  "Shared input workflow",
  function() {
    tryCatch(
      {
        sys.source(
          file.path(
            "R",
            "globals.R"
          ),
          envir = globalenv()
        )

        sys.source(
          file.path(
            "R",
            "io_load.R"
          ),
          envir = globalenv()
        )
      },
      error = function(error_condition) {
        input_loading_success <<- FALSE
        stop(
          conditionMessage(error_condition),
          call. = FALSE
        )
      }
    )

    required_objects <- c(
      "report",
      "runs_removed_by_global_filter",
      "peptides_clean",
      "fasta_ref",
      "results_dir",
      "file_palette"
    )

    missing_objects <- required_objects[
      !vapply(
        required_objects,
        exists,
        logical(1),
        envir = globalenv(),
        inherits = FALSE
      )
    ]

    assert_true(
      length(missing_objects) == 0L,
      paste0(
        "The input workflow did not create: ",
        paste(
          missing_objects,
          collapse = ", "
        )
      )
    )
  }
)


if (!isTRUE(input_loading_success)) {
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
      "scientific-core-results.tsv"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = ""
  )

  quit(
    status = 1L
  )
}


## 7. Dynamic input baseline ---------------------------------------------------

baseline_report_rows <- nrow(
  report
)

baseline_runs <- sort(
  unique(
    as.character(
      report$Run
    )
  )
)

baseline_run_count <- length(
  baseline_runs
)

baseline_peptide_group_count <- nrow(
  peptides_clean
)

baseline_fasta_entry_count <- nrow(
  fasta_ref
)


run_check(
  "Filtered report contains observations",
  function() {
    assert_true(
      is.data.frame(report),
      "The filtered report is not a data frame."
    )

    assert_true(
      baseline_report_rows > 0L,
      "The filtered report contains no observations."
    )
  }
)


run_check(
  "Filtered report contains valid runs",
  function() {
    assert_true(
      "Run" %in% names(report),
      "The filtered report does not contain the Run column."
    )

    assert_true(
      baseline_run_count > 0L,
      "The filtered report contains no valid runs."
    )

    assert_true(
      all(
        !is.na(baseline_runs) &
          nzchar(baseline_runs)
      ),
      "The filtered report contains missing or blank run identifiers."
    )
  }
)


run_check(
  "Peptide evidence is available",
  function() {
    assert_true(
      is.data.frame(peptides_clean),
      "peptides_clean is not a data frame."
    )

    assert_true(
      baseline_peptide_group_count > 0L,
      "No protein groups with peptide evidence are available."
    )
  }
)


run_check(
  "FASTA reference is available",
  function() {
    assert_true(
      is.data.frame(fasta_ref),
      "fasta_ref is not a data frame."
    )

    assert_true(
      baseline_fasta_entry_count > 0L,
      "The parsed FASTA reference contains no entries."
    )
  }
)


## 8. Defensive plotting defaults ---------------------------------------------

if (!exists(
  "PANELS_PER_ROW",
  inherits = FALSE
)) {
  PANELS_PER_ROW <- 3L
}

if (!exists(
  "AR_DENSITY",
  inherits = FALSE
)) {
  AR_DENSITY <- 1
}

if (!exists(
  "AR_BOXPLOT",
  inherits = FALSE
)) {
  AR_BOXPLOT <- 1
}

if (!exists(
  "legend_fill",
  mode = "function",
  inherits = FALSE
)) {
  legend_fill <- function() {
    ggplot2::theme()
  }
}

if (!exists(
  "legend_color",
  mode = "function",
  inherits = FALSE
)) {
  legend_color <- function() {
    ggplot2::theme()
  }
}


## 9. Integrated scientific module execution ----------------------------------

modules <- c(
  "qc_identification.R",
  "qc_abundance.R",
  "qc_chromatography.R",
  "qc_chromatography_internalstd.R",
  "qc_completeness.R",
  "qc_coverage.R",
  "qc_diann_quality_metrics.R",
  "qc_intensity.R",
  "qc_ion_mobility.R",
  "qc_mass_error.R",
  "qc_peptides.R",
  "qc_figures_of_merit.R"
)


cat("\n")
cat("============================================================\n")
cat("INTEGRATED SCIENTIFIC MODULE EXECUTION\n")
cat("============================================================\n")


module_execution_status <- setNames(
  rep(
    FALSE,
    length(modules)
  ),
  modules
)


for (module in modules) {
  module_path <- file.path(
    "R",
    module
  )

  run_check(
    module,
    function() {
      assert_true(
        file.exists(module_path),
        paste0(
          "Scientific module was not found: ",
          module_path
        )
      )

      sys.source(
        module_path,
        envir = globalenv()
      )

      module_execution_status[[module]] <<- TRUE

      invisible(TRUE)
    }
  )
}


## 10. Shared-object preservation ---------------------------------------------

run_check(
  "Filtered report remains available",
  function() {
    assert_true(
      exists(
        "report",
        envir = globalenv(),
        inherits = FALSE
      ),
      "The filtered report was removed during module execution."
    )

    assert_true(
      is.data.frame(report),
      "The filtered report is no longer a data frame."
    )

    assert_true(
      nrow(report) == baseline_report_rows,
      paste0(
        "The number of report rows changed during module execution. ",
        "Before: ",
        baseline_report_rows,
        "; after: ",
        nrow(report),
        "."
      )
    )
  }
)


run_check(
  "Run inventory remains unchanged",
  function() {
    observed_runs <- sort(
      unique(
        as.character(
          report$Run
        )
      )
    )

    assert_true(
      identical(
        observed_runs,
        baseline_runs
      ),
      "The run inventory changed during scientific module execution."
    )
  }
)


run_check(
  "Peptide evidence remains available",
  function() {
    assert_true(
      exists(
        "peptides_clean",
        envir = globalenv(),
        inherits = FALSE
      ),
      "peptides_clean was removed during module execution."
    )

    assert_true(
      is.data.frame(peptides_clean),
      "peptides_clean is no longer a data frame."
    )

    assert_true(
      nrow(peptides_clean) == baseline_peptide_group_count,
      paste0(
        "The peptide-evidence row count changed during module execution. ",
        "Before: ",
        baseline_peptide_group_count,
        "; after: ",
        nrow(peptides_clean),
        "."
      )
    )
  }
)


run_check(
  "FASTA reference remains available",
  function() {
    assert_true(
      exists(
        "fasta_ref",
        envir = globalenv(),
        inherits = FALSE
      ),
      "fasta_ref was removed during module execution."
    )

    assert_true(
      is.data.frame(fasta_ref),
      "fasta_ref is no longer a data frame."
    )

    assert_true(
      nrow(fasta_ref) == baseline_fasta_entry_count,
      paste0(
        "The FASTA-reference row count changed during module execution. ",
        "Before: ",
        baseline_fasta_entry_count,
        "; after: ",
        nrow(fasta_ref),
        "."
      )
    )
  }
)


## 11. Completeness integration invariants ------------------------------------

run_check(
  "Completeness global table was created",
  function() {
    assert_true(
      exists(
        "completeness_global_df",
        envir = globalenv(),
        inherits = FALSE
      ),
      "completeness_global_df was not created."
    )

    assert_true(
      is.data.frame(completeness_global_df),
      "completeness_global_df is not a data frame."
    )

    assert_true(
      nrow(completeness_global_df) > 0L,
      "completeness_global_df is empty."
    )

    required_columns <- c(
      "Metric",
      "Value"
    )

    assert_true(
      all(
        required_columns %in%
          names(completeness_global_df)
      ),
      paste0(
        "completeness_global_df is missing required columns: ",
        paste(
          setdiff(
            required_columns,
            names(completeness_global_df)
          ),
          collapse = ", "
        )
      )
    )
  }
)


run_check(
  "Completeness per-run table matches run inventory",
  function() {
    assert_true(
      exists(
        "completeness_per_run_df",
        envir = globalenv(),
        inherits = FALSE
      ),
      "completeness_per_run_df was not created."
    )

    assert_true(
      is.data.frame(completeness_per_run_df),
      "completeness_per_run_df is not a data frame."
    )

    assert_true(
      "Run" %in% names(completeness_per_run_df),
      "completeness_per_run_df does not contain the Run column."
    )

    observed_runs <- sort(
      unique(
        as.character(
          completeness_per_run_df$Run
        )
      )
    )

    assert_true(
      identical(
        observed_runs,
        baseline_runs
      ),
      paste0(
        "The completeness table run inventory does not match ",
        "the filtered report."
      )
    )

    assert_true(
      nrow(completeness_per_run_df) ==
        baseline_run_count,
      paste0(
        "Expected one completeness row per run. ",
        "Expected: ",
        baseline_run_count,
        "; observed: ",
        nrow(completeness_per_run_df),
        "."
      )
    )
  }
)


run_check(
  "Completeness percentages are within valid bounds",
  function() {
    assert_true(
      "Completeness_Pct" %in%
        names(completeness_per_run_df),
      paste0(
        "completeness_per_run_df does not contain ",
        "Completeness_Pct."
      )
    )

    completeness_values <-
      completeness_per_run_df$Completeness_Pct

    assert_true(
      all(
        is.finite(
          completeness_values
        )
      ),
      "One or more per-run completeness values are not finite."
    )

    assert_true(
      all(
        completeness_values >= 0 &
          completeness_values <= 100
      ),
      "One or more per-run completeness values are outside 0-100%."
    )
  }
)


## 12. Coverage integration invariants ----------------------------------------

run_check(
  "Coverage table was created",
  function() {
    assert_true(
      exists(
        "coverage_results",
        envir = globalenv(),
        inherits = FALSE
      ),
      "coverage_results was not created."
    )

    assert_true(
      is.data.frame(coverage_results),
      "coverage_results is not a data frame."
    )

    assert_true(
      nrow(coverage_results) > 0L,
      "coverage_results is empty."
    )

    required_columns <- c(
      "Protein.Group",
      "Coverage_Pct"
    )

    assert_true(
      all(
        required_columns %in%
          names(coverage_results)
      ),
      paste0(
        "coverage_results is missing required columns: ",
        paste(
          setdiff(
            required_columns,
            names(coverage_results)
          ),
          collapse = ", "
        )
      )
    )
  }
)


run_check(
  "Coverage values are within valid bounds",
  function() {
    coverage_values <- coverage_results$Coverage_Pct

    valid_values <- coverage_values[
      !is.na(coverage_values)
    ]

    assert_true(
      length(valid_values) > 0L,
      "No non-missing protein coverage values were produced."
    )

    assert_true(
      all(
        is.finite(
          valid_values
        )
      ),
      "One or more protein coverage values are not finite."
    )

    assert_true(
      all(
        valid_values >= 0 &
          valid_values <= 100
      ),
      "One or more protein coverage values are outside 0-100%."
    )
  }
)


run_check(
  "Coverage table contains no duplicated protein rows",
  function() {
    assert_true(
      !anyDuplicated(
        coverage_results$Protein.Group
      ),
      "coverage_results contains duplicated Protein.Group rows."
    )
  }
)


run_check(
  "Coverage CSV was created",
  function() {
    coverage_files <- list.files(
      results_dir,
      pattern = "^protein_coverage_with_peptides_.*[.]csv$",
      recursive = TRUE,
      full.names = TRUE
    )

    assert_true(
      length(coverage_files) > 0L,
      paste0(
        "No protein coverage CSV matching ",
        "protein_coverage_with_peptides_*.csv was found in: ",
        results_dir
      )
    )

    coverage_files <- coverage_files[
      file.exists(coverage_files)
    ]

    assert_true(
      length(coverage_files) > 0L,
      "The detected protein coverage CSV does not exist."
    )

    coverage_file_info <- file.info(
      coverage_files
    )

    valid_coverage_files <- coverage_files[
      !is.na(coverage_file_info$size) &
        coverage_file_info$size > 0
    ]

    assert_true(
      length(valid_coverage_files) > 0L,
      "No readable and non-empty protein coverage CSV was generated."
    )

    coverage_csv_columns <- names(
      utils::read.csv(
        valid_coverage_files[[1L]],
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
      coverage_csv_columns
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
)


## 13. Figures of Merit integration invariants --------------------------------

run_check(
  "Figures of Merit table was created",
  function() {
    assert_true(
      exists(
        "figures_of_merit_df",
        envir = globalenv(),
        inherits = FALSE
      ),
      "figures_of_merit_df was not created."
    )

    assert_true(
      is.data.frame(figures_of_merit_df),
      "figures_of_merit_df is not a data frame."
    )

    assert_true(
      nrow(figures_of_merit_df) > 0L,
      "figures_of_merit_df is empty."
    )

    assert_true(
      "Run" %in% names(figures_of_merit_df),
      "figures_of_merit_df does not contain the Run column."
    )
  }
)


run_check(
  "Figures of Merit table matches run inventory",
  function() {
    figures_of_merit_runs <- sort(
      unique(
        as.character(
          figures_of_merit_df$Run
        )
      )
    )

    assert_true(
      identical(
        figures_of_merit_runs,
        baseline_runs
      ),
      paste0(
        "The Figures of Merit run inventory does not match ",
        "the filtered report."
      )
    )

    assert_true(
      nrow(figures_of_merit_df) ==
        baseline_run_count,
      paste0(
        "Expected one Figures of Merit row per run. ",
        "Expected: ",
        baseline_run_count,
        "; observed: ",
        nrow(figures_of_merit_df),
        "."
      )
    )
  }
)


run_check(
  "Figures of Merit CSV was created",
  function() {
    assert_true(
      exists(
        "figures_of_merit_csv",
        envir = globalenv(),
        inherits = FALSE
      ),
      "figures_of_merit_csv was not created."
    )

    assert_true(
      length(figures_of_merit_csv) == 1L &&
        !is.na(figures_of_merit_csv) &&
        nzchar(figures_of_merit_csv),
      "The Figures of Merit CSV path is missing or invalid."
    )

    assert_true(
      file.exists(figures_of_merit_csv),
      paste0(
        "The Figures of Merit CSV was not found: ",
        figures_of_merit_csv
      )
    )

    figures_of_merit_file_info <- file.info(
      figures_of_merit_csv
    )

    assert_true(
      !is.na(figures_of_merit_file_info$size) &&
        figures_of_merit_file_info$size > 0,
      "The Figures of Merit CSV is empty or inaccessible."
    )
  }
)


## 14. Generated-output integration invariants --------------------------------

run_check(
  "Scientific results directory exists",
  function() {
    assert_true(
      dir.exists(results_dir),
      paste0(
        "The scientific results directory does not exist: ",
        results_dir
      )
    )
  }
)


run_check(
  "Scientific modules generated output files",
  function() {
    generated_files <- list.files(
      results_dir,
      recursive = TRUE,
      full.names = TRUE
    )

    generated_files <- generated_files[
      file.exists(generated_files) &
        !dir.exists(generated_files)
    ]

    assert_true(
      length(generated_files) > 0L,
      "The scientific modules did not generate any output files."
    )

    generated_file_sizes <- file.info(
      generated_files
    )$size

    assert_true(
      all(
        !is.na(generated_file_sizes) &
          generated_file_sizes > 0
      ),
      "One or more generated scientific output files are empty."
    )
  }
)


## 15. Graphics-device completion ---------------------------------------------

if (isTRUE(graphics_device_open)) {
  grDevices::dev.off()
  graphics_device_open <- FALSE
}


run_check(
  "Integrated plot file was created",
  function() {
    assert_true(
      file.exists(plot_file),
      "The integrated plot file was not created."
    )

    plot_file_info <- file.info(
      plot_file
    )

    assert_true(
      !is.na(plot_file_info$size) &&
        plot_file_info$size > 0,
      "The integrated plot file is empty or inaccessible."
    )
  }
)


## 16. Consolidated results ----------------------------------------------------

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
  "scientific-core-results.tsv"
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
cat("SCIENTIFIC CORE TEST SUMMARY\n")
cat("============================================================\n")
cat(
  "Modules evaluated: ",
  length(modules),
  "\n",
  sep = ""
)
cat(
  "Runs evaluated: ",
  baseline_run_count,
  "\n",
  sep = ""
)
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