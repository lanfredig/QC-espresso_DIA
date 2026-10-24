###############################################################################
# File: tests/run_tests.R
#
# Purpose:
#   Validate the QC-espresso repository structure, computational environment,
#   example inputs, reference-output inventory, and input-loading baseline.
#
# Validation scope:
#   - Confirm that required repository files and R modules are available.
#   - Parse every repository R script to detect syntax errors.
#   - Validate the QC-espresso and R environment versions.
#   - Confirm availability and lockfile registration of direct dependencies.
#   - Validate the example DIA-NN report and FASTA reference.
#   - Validate the DIA-NN report schema.
#   - Validate the reference-output manifest and referenced files.
#   - Execute the input-loading workflow in an isolated environment.
#   - Validate dynamic relationships among loaded report, run, peptide, and
#     FASTA objects without requiring dataset-specific absolute counts.
#
# Result categories:
#   PASS  The expected repository or baseline condition was satisfied.
#   FAIL  The repository or baseline condition was not satisfied.
#
# Usage:
#   Rscript --no-init-file tests/run_tests.R
###############################################################################


## 1. Minimal test framework ---------------------------------------------------

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
        name = test_name,
        success = TRUE,
        message = if (
          is.character(details) &&
          length(details) == 1L
        ) {
          details
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
        name = test_name,
        success = FALSE,
        message = conditionMessage(
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

  results[[length(results) + 1L]] <<- outcome

  cat(
    sprintf(
      "%-46s %s\n",
      test_name,
      if (isTRUE(outcome$success)) {
        "PASS"
      } else {
        "FAIL"
      }
    )
  )

  if (!isTRUE(outcome$success)) {
    cat(
      "  Reason: ",
      outcome$message,
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


## 3. Repository expectations --------------------------------------------------

EXPECTED_R_SCRIPTS <- sort(
  c(
    "globals.R",
    "io_load.R",
    "qc_abundance.R",
    "qc_chromatography.R",
    "qc_chromatography_internalstd.R",
    "qc_completeness.R",
    "qc_coverage.R",
    "qc_diann_quality_metrics.R",
    "qc_figures_of_merit.R",
    "qc_identification.R",
    "qc_intensity.R",
    "qc_ion_mobility.R",
    "qc_mass_error.R",
    "qc_peptides.R",
    "qc_session.R",
    "validate_inputs.R"
  )
)


REQUIRED_ROOT_FILES <- c(
  "QC-espresso_DIA.Rproj",
  "VERSION",
  "run_qc.R",
  "renv.lock",
  ".Rprofile",
  ".gitignore"
)


REQUIRED_SCIENTIFIC_FILES <- c(
  file.path(
    "report",
    "qc_espresso_scientific.qmd"
  )
)


EXAMPLE_PARQUET <- file.path(
  "example-data",
  "report.parquet"
)


EXAMPLE_FASTA <- file.path(
  "example-data",
  "uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta"
)


REFERENCE_DIR <- file.path(
  "tests",
  "reference-output",
  "v1.3.1"
)


REFERENCE_MANIFEST <- file.path(
  REFERENCE_DIR,
  "SHA256SUMS.txt"
)


DIRECT_PACKAGES <- c(
  "arrow",
  "dplyr",
  "forcats",
  "ggplot2",
  "ggrepel",
  "gt",
  "gtExtras",
  "knitr",
  "microseq",
  "purrr",
  "RColorBrewer",
  "rlang",
  "scales",
  "stringr",
  "tibble",
  "tidyr",
  "tidyverse"
)


## 4. Input-validation functions ----------------------------------------------

source(
  file.path(
    "R",
    "validate_inputs.R"
  ),
  local = globalenv()
)


## 5. Execution header ---------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("QC-ESPRESSO REPOSITORY AND BASELINE VALIDATION\n")
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
cat("============================================================\n")
cat("\n")


## 6. Repository structure -----------------------------------------------------

record_test(
  "Repository root files",
  function() {
    missing_files <- REQUIRED_ROOT_FILES[
      !file.exists(
        REQUIRED_ROOT_FILES
      )
    ]

    assert_true(
      length(missing_files) == 0L,
      paste0(
        "Missing required root files: ",
        paste(
          missing_files,
          collapse = ", "
        )
      )
    )

    file_sizes <- file.info(
      REQUIRED_ROOT_FILES
    )$size

    assert_true(
      all(
        !is.na(file_sizes) &
          file_sizes > 0
      ),
      "One or more required root files are empty or inaccessible."
    )
  }
)


record_test(
  "Scientific report files",
  function() {
    missing_files <- REQUIRED_SCIENTIFIC_FILES[
      !file.exists(
        REQUIRED_SCIENTIFIC_FILES
      )
    ]

    assert_true(
      length(missing_files) == 0L,
      paste0(
        "Missing scientific report files: ",
        paste(
          missing_files,
          collapse = ", "
        )
      )
    )

    file_sizes <- file.info(
      REQUIRED_SCIENTIFIC_FILES
    )$size

    assert_true(
      all(
        !is.na(file_sizes) &
          file_sizes > 0
      ),
      "One or more scientific report files are empty or inaccessible."
    )
  }
)


record_test(
  "R module inventory",
  function() {
    observed_scripts <- sort(
      list.files(
        "R",
        pattern = "[.]R$",
        full.names = FALSE
      )
    )

    assert_equal(
      observed_scripts,
      EXPECTED_R_SCRIPTS,
      paste0(
        "The R module inventory differs from the expected ",
        "repository inventory."
      )
    )

    script_sizes <- file.info(
      file.path(
        "R",
        observed_scripts
      )
    )$size

    assert_true(
      all(
        !is.na(script_sizes) &
          script_sizes > 0
      ),
      "One or more R modules are empty or inaccessible."
    )
  }
)


record_test(
  "No deprecated diagnostic script",
  function() {
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
  }
)


## 7. R syntax validation ------------------------------------------------------

record_test(
  "Repository R script syntax",
  function() {
    repository_scripts <- unique(
      c(
        list.files(
          "R",
          pattern = "[.]R$",
          full.names = TRUE
        ),
        list.files(
          "tests",
          pattern = "[.]R$",
          full.names = TRUE,
          recursive = TRUE
        ),
        "run_qc.R"
      )
    )

    repository_scripts <- repository_scripts[
      file.exists(
        repository_scripts
      )
    ]

    assert_true(
      length(repository_scripts) > 0L,
      "No repository R scripts were found for syntax validation."
    )

    syntax_errors <- character()

    for (script_path in repository_scripts) {
      syntax_message <- tryCatch(
        {
          parse(
            file = script_path
          )

          ""
        },
        error = function(error_condition) {
          conditionMessage(
            error_condition
          )
        }
      )

      if (nzchar(syntax_message)) {
        syntax_errors[
          script_path
        ] <- syntax_message
      }
    }

    assert_true(
      length(syntax_errors) == 0L,
      paste(
        c(
          "Syntax errors were detected:",
          paste0(
            "  ",
            names(syntax_errors),
            ": ",
            syntax_errors
          )
        ),
        collapse = "\n"
      )
    )
  }
)


## 8. Project and environment versions ----------------------------------------

record_test(
  "QC-espresso version",
  function() {
    version_lines <- readLines(
      "VERSION",
      warn = FALSE
    )

    assert_true(
      length(version_lines) > 0L,
      "The VERSION file is empty."
    )

    project_version <- trimws(
      version_lines[[1L]]
    )

    assert_true(
      nzchar(
        project_version
      ),
      "The QC-espresso version is blank."
    )

    assert_true(
      grepl(
        "^[0-9]+[.][0-9]+[.][0-9]+([.-][A-Za-z0-9]+)?$",
        project_version
      ),
      paste0(
        "The QC-espresso version is not a recognized semantic ",
        "version: ",
        project_version
      )
    )

    project_version
  }
)


record_test(
  "Reference R version",
  function() {
    observed_version <- paste(
      R.version$major,
      R.version$minor,
      sep = "."
    )

    expected_version <- "4.6.0"

    assert_equal(
      observed_version,
      expected_version,
      "Tests must run with the validated R environment."
    )
  }
)


record_test(
  "Required package availability",
  function() {
    availability <- vapply(
      DIRECT_PACKAGES,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )

    missing_packages <- DIRECT_PACKAGES[
      !availability
    ]

    assert_true(
      length(missing_packages) == 0L,
      paste0(
        "Missing required packages: ",
        paste(
          missing_packages,
          collapse = ", "
        )
      )
    )
  }
)


record_test(
  "Package-lockfile equivalence",
  function() {
    assert_true(
      requireNamespace(
        "jsonlite",
        quietly = TRUE
      ),
      "jsonlite is required to inspect renv.lock."
    )

    lockfile <- jsonlite::read_json(
      "renv.lock",
      simplifyVector = FALSE
    )

    locked_r_version <- lockfile[["R"]][["Version"]]

    assert_true(
      !is.null(
        locked_r_version
      ) &&
        length(
          locked_r_version
        ) == 1L &&
        nzchar(
          locked_r_version
        ),
      "renv.lock does not record a valid R version."
    )

    observed_r_version <- paste(
      R.version$major,
      R.version$minor,
      sep = "."
    )

    assert_equal(
      locked_r_version,
      observed_r_version,
      paste0(
        "The running R version differs from the version ",
        "recorded in renv.lock."
      )
    )

    installed_packages <- installed.packages()

    unavailable_installed_packages <- setdiff(
      DIRECT_PACKAGES,
      rownames(
        installed_packages
      )
    )

    assert_true(
      length(
        unavailable_installed_packages
      ) == 0L,
      paste0(
        "Required packages are absent from installed.packages(): ",
        paste(
          unavailable_installed_packages,
          collapse = ", "
        )
      )
    )

    functional_versions <- vapply(
      DIRECT_PACKAGES,
      function(package_name) {
        installed_packages[
          package_name,
          "Version"
        ]
      },
      character(1)
    )

    locked_versions <- vapply(
      DIRECT_PACKAGES,
      function(package_name) {
        package_entry <- lockfile[["Packages"]][[package_name]]

        if (is.null(package_entry)) {
          return(
            NA_character_
          )
        }

        package_version <- package_entry[["Version"]]

        if (
          is.null(package_version) ||
          length(package_version) != 1L
        ) {
          return(
            NA_character_
          )
        }

        package_version
      },
      character(1)
    )

    matching_versions <-
      functional_versions ==
      locked_versions

    matching_versions[
      is.na(
        matching_versions
      )
    ] <- FALSE

    if (!all(matching_versions)) {
      differing_packages <- DIRECT_PACKAGES[
        !matching_versions
      ]

      difference_details <- vapply(
        differing_packages,
        function(package_name) {
          paste0(
            package_name,
            ": installed=",
            functional_versions[[package_name]],
            ", lockfile=",
            locked_versions[[package_name]]
          )
        },
        character(1)
      )

      stop(
        paste(
          c(
            "Package version differences were detected:",
            paste0(
              "  ",
              difference_details
            )
          ),
          collapse = "\n"
        ),
        call. = FALSE
      )
    }

    invisible(TRUE)
  }
)


## 9. Example-input validation -------------------------------------------------

record_test(
  "Example file presence",
  function() {
    example_files <- c(
      EXAMPLE_PARQUET,
      EXAMPLE_FASTA
    )

    missing_files <- example_files[
      !file.exists(
        example_files
      )
    ]

    assert_true(
      length(missing_files) == 0L,
      paste0(
        "Missing example files: ",
        paste(
          missing_files,
          collapse = ", "
        )
      )
    )

    file_sizes <- file.info(
      example_files
    )$size

    assert_true(
      all(
        !is.na(file_sizes) &
          file_sizes > 0
      ),
      "One or more example files are empty or inaccessible."
    )
  }
)


record_test(
  "Example Parquet validation",
  function() {
    validate_report_file(
      EXAMPLE_PARQUET
    )
  }
)


record_test(
  "Example FASTA validation",
  function() {
    validate_fasta_file(
      EXAMPLE_FASTA
    )
  }
)


record_test(
  "DIA-NN report schema",
  function() {
    schema_result <- validate_report_schema(
      EXAMPLE_PARQUET
    )

    assert_true(
      is.list(
        schema_result
      ),
      "validate_report_schema did not return the expected result."
    )

    assert_true(
      "missing_columns" %in%
        names(
          schema_result
        ),
      paste0(
        "The report-schema validation result does not contain ",
        "missing_columns."
      )
    )

    assert_true(
      length(
        schema_result$missing_columns
      ) == 0L,
      paste0(
        "Missing required DIA-NN columns: ",
        paste(
          schema_result$missing_columns,
          collapse = ", "
        )
      )
    )
  }
)


## 10. Reference-output inventory ---------------------------------------------

record_test(
  "Reference-output directory",
  function() {
    assert_true(
      dir.exists(
        REFERENCE_DIR
      ),
      paste0(
        "Reference-output directory was not found: ",
        REFERENCE_DIR
      )
    )

    assert_true(
      file.exists(
        REFERENCE_MANIFEST
      ),
      paste0(
        "Reference-output manifest was not found: ",
        REFERENCE_MANIFEST
      )
    )

    manifest_info <- file.info(
      REFERENCE_MANIFEST
    )

    assert_true(
      !is.na(
        manifest_info$size
      ) &&
        manifest_info$size > 0,
      "The reference-output manifest is empty or inaccessible."
    )
  }
)


reference_manifest_entries <- NULL


record_test(
  "Reference-output manifest format",
  function() {
    manifest_lines <- readLines(
      REFERENCE_MANIFEST,
      warn = FALSE
    )

    manifest_lines <- trimws(
      manifest_lines
    )

    manifest_lines <- manifest_lines[
      nzchar(
        manifest_lines
      )
    ]

    assert_true(
      length(manifest_lines) > 0L,
      "SHA256SUMS.txt contains no reference-output entries."
    )

    manifest_pattern <-
      "^([A-Fa-f0-9]{64})[[:space:]]{2}(.+)$"

    valid_lines <- grepl(
      manifest_pattern,
      manifest_lines,
      perl = TRUE
    )

    assert_true(
      all(
        valid_lines
      ),
      paste0(
        "One or more SHA256SUMS.txt entries do not follow ",
        "the expected SHA-256 and filename format."
      )
    )

    reference_hashes <- sub(
      manifest_pattern,
      "\\1",
      manifest_lines,
      perl = TRUE
    )

    reference_files <- sub(
      manifest_pattern,
      "\\2",
      manifest_lines,
      perl = TRUE
    )

    assert_true(
      !anyDuplicated(
        reference_files
      ),
      "SHA256SUMS.txt contains duplicated filenames."
    )

    reference_manifest_entries <<- data.frame(
      SHA256 = tolower(
        reference_hashes
      ),
      File = reference_files,
      stringsAsFactors = FALSE
    )

    invisible(TRUE)
  }
)


record_test(
  "Referenced output files",
  function() {
    assert_true(
      !is.null(
        reference_manifest_entries
      ),
      paste0(
        "Reference-output manifest entries were not available ",
        "for file validation."
      )
    )

    referenced_paths <- file.path(
      REFERENCE_DIR,
      reference_manifest_entries$File
    )

    missing_files <- reference_manifest_entries$File[
      !file.exists(
        referenced_paths
      )
    ]

    assert_true(
      length(missing_files) == 0L,
      paste0(
        "Files referenced by SHA256SUMS.txt are missing: ",
        paste(
          missing_files,
          collapse = ", "
        )
      )
    )

    referenced_sizes <- file.info(
      referenced_paths
    )$size

    assert_true(
      all(
        !is.na(referenced_sizes) &
          referenced_sizes > 0
      ),
      paste0(
        "One or more files referenced by SHA256SUMS.txt ",
        "are empty or inaccessible."
      )
    )

    referenced_extensions <- tolower(
      tools::file_ext(
        referenced_paths
      )
    )

    supported_extensions <- c(
      "csv",
      "html",
      "htm",
      "tsv",
      "txt"
    )

    unsupported_files <-
      reference_manifest_entries$File[
        !referenced_extensions %in%
          supported_extensions
      ]

    assert_true(
      length(unsupported_files) == 0L,
      paste0(
        "Unsupported reference-output file extensions: ",
        paste(
          unsupported_files,
          collapse = ", "
        )
      )
    )
  }
)


record_test(
  "Reference-output readability",
  function() {
    assert_true(
      !is.null(
        reference_manifest_entries
      ),
      "Reference-output manifest entries are unavailable."
    )

    referenced_paths <- file.path(
      REFERENCE_DIR,
      reference_manifest_entries$File
    )

    readable_files <- vapply(
      referenced_paths,
      function(reference_path) {
        tryCatch(
          {
            first_byte <- readBin(
              reference_path,
              what = "raw",
              n = 1L
            )

            length(first_byte) == 1L
          },
          error = function(error_condition) {
            FALSE
          }
        )
      },
      logical(1)
    )

    unreadable_files <- basename(
      referenced_paths[
        !readable_files
      ]
    )

    assert_true(
      length(unreadable_files) == 0L,
      paste0(
        "Unreadable reference-output files: ",
        paste(
          unreadable_files,
          collapse = ", "
        )
      )
    )
  }
)


## 11. Isolated input-loading workflow ----------------------------------------

input_environment <- new.env(
  parent = globalenv()
)


input_output_dir <- file.path(
  tempdir(),
  paste0(
    "qc_espresso_input_test_",
    Sys.getpid()
  )
)

dir.create(
  input_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

on.exit(
  unlink(
    input_output_dir,
    recursive = TRUE,
    force = TRUE
  ),
  add = TRUE
)


input_environment$PROJECT_ROOT <- PROJECT_ROOT

input_environment$params <- list(
  report_path = normalizePath(
    EXAMPLE_PARQUET,
    winslash = "/",
    mustWork = TRUE
  ),
  fasta_path = normalizePath(
    EXAMPLE_FASTA,
    winslash = "/",
    mustWork = TRUE
  ),
  output_dir = input_output_dir
)


record_test(
  "Input loading workflow",
  function() {
    suppressPackageStartupMessages(
      {
        library(tidyverse)
        library(arrow)
        library(microseq)
        library(scales)
      }
    )

    sys.source(
      file.path(
        "R",
        "globals.R"
      ),
      envir = input_environment
    )

    sys.source(
      file.path(
        "R",
        "io_load.R"
      ),
      envir = input_environment
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
        envir = input_environment,
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


## 12. Dynamic input baseline --------------------------------------------------

record_test(
  "Loaded report baseline",
  function() {
    loaded_report <- get(
      "report",
      envir = input_environment,
      inherits = FALSE
    )

    assert_true(
      is.data.frame(
        loaded_report
      ),
      "The loaded report is not a data frame."
    )

    assert_true(
      nrow(
        loaded_report
      ) > 0L,
      "The loaded report contains no retained observations."
    )

    required_columns <- c(
      "Run",
      "Protein.Group",
      "Precursor.Id",
      "Precursor.Quantity"
    )

    missing_columns <- setdiff(
      required_columns,
      names(
        loaded_report
      )
    )

    assert_true(
      length(missing_columns) == 0L,
      paste0(
        "The loaded report is missing baseline columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
)


record_test(
  "Loaded run inventory",
  function() {
    loaded_report <- get(
      "report",
      envir = input_environment,
      inherits = FALSE
    )

    run_values <- unique(
      as.character(
        loaded_report$Run
      )
    )

    assert_true(
      length(run_values) > 0L,
      "The loaded report contains no runs."
    )

    assert_true(
      all(
        !is.na(run_values) &
          nzchar(run_values)
      ),
      "The loaded report contains missing or blank run identifiers."
    )

    file_palette <- get(
      "file_palette",
      envir = input_environment,
      inherits = FALSE
    )

    assert_true(
      is.character(
        file_palette
      ),
      "file_palette is not a character vector."
    )

    assert_true(
      length(file_palette) > 0L,
      "file_palette is empty."
    )

    assert_true(
      all(
        !is.na(file_palette) &
          nzchar(file_palette)
      ),
      "file_palette contains missing or blank color values."
    )

    palette_names <- names(
      file_palette
    )

    palette_is_named <-
      !is.null(palette_names) &&
      length(palette_names) ==
      length(file_palette) &&
      all(
        !is.na(palette_names) &
          nzchar(palette_names)
      )

    if (palette_is_named) {
      assert_true(
        all(
          run_values %in%
            palette_names
        ),
        paste0(
          "The named file palette does not contain a color ",
          "for every retained run."
        )
      )
    } else {
      assert_true(
        length(file_palette) >=
          length(run_values),
        paste0(
          "The positional file palette does not contain enough ",
          "colors for all retained runs. Runs: ",
          length(run_values),
          "; colors: ",
          length(file_palette),
          "."
        )
      )
    }
  }
)


record_test(
  "Removed-run inventory",
  function() {
    removed_runs <- get(
      "runs_removed_by_global_filter",
      envir = input_environment,
      inherits = FALSE
    )

    assert_true(
      is.data.frame(
        removed_runs
      ),
      "runs_removed_by_global_filter is not a data frame."
    )

    if (nrow(removed_runs) > 0L) {
      assert_true(
        "Run" %in%
          names(
            removed_runs
          ),
        paste0(
          "runs_removed_by_global_filter contains observations ",
          "but no Run column."
        )
      )

      retained_report <- get(
        "report",
        envir = input_environment,
        inherits = FALSE
      )

      retained_runs <- unique(
        as.character(
          retained_report$Run
        )
      )

      removed_run_values <- unique(
        as.character(
          removed_runs$Run
        )
      )

      assert_true(
        length(
          intersect(
            retained_runs,
            removed_run_values
          )
        ) == 0L,
        paste0(
          "One or more runs are simultaneously classified as ",
          "retained and removed."
        )
      )
    }
  }
)


record_test(
  "Peptide evidence baseline",
  function() {
    peptide_evidence <- get(
      "peptides_clean",
      envir = input_environment,
      inherits = FALSE
    )

    assert_true(
      is.data.frame(
        peptide_evidence
      ),
      "peptides_clean is not a data frame."
    )

    assert_true(
      nrow(
        peptide_evidence
      ) > 0L,
      "No protein groups with peptide evidence were generated."
    )

    assert_true(
      "Protein.Group" %in%
        names(
          peptide_evidence
        ),
      "peptides_clean does not contain Protein.Group."
    )

    assert_true(
      !anyDuplicated(
        peptide_evidence$Protein.Group
      ),
      "peptides_clean contains duplicated Protein.Group rows."
    )
  }
)


record_test(
  "FASTA parsing baseline",
  function() {
    parsed_fasta <- get(
      "fasta_ref",
      envir = input_environment,
      inherits = FALSE
    )

    assert_true(
      is.data.frame(
        parsed_fasta
      ),
      "fasta_ref is not a data frame."
    )

    assert_true(
      nrow(
        parsed_fasta
      ) > 0L,
      "The parsed FASTA reference contains no entries."
    )

    required_columns <- c(
      "Protein.Group",
      "Sequence"
    )

    missing_columns <- setdiff(
      required_columns,
      names(
        parsed_fasta
      )
    )

    assert_true(
      length(missing_columns) == 0L,
      paste0(
        "fasta_ref is missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )

    assert_true(
      all(
        !is.na(
          parsed_fasta$Protein.Group
        ) &
          nzchar(
            parsed_fasta$Protein.Group
          )
      ),
      "The parsed FASTA contains missing or blank accessions."
    )

    assert_true(
      all(
        !is.na(
          parsed_fasta$Sequence
        ) &
          nzchar(
            parsed_fasta$Sequence
          )
      ),
      "The parsed FASTA contains missing or blank sequences."
    )
  }
)


record_test(
  "Results directory baseline",
  function() {
    loaded_results_dir <- get(
      "results_dir",
      envir = input_environment,
      inherits = FALSE
    )

    assert_true(
      is.character(
        loaded_results_dir
      ) &&
        length(
          loaded_results_dir
        ) == 1L &&
        !is.na(
          loaded_results_dir
        ) &&
        nzchar(
          loaded_results_dir
        ),
      "results_dir is missing or invalid."
    )

    assert_true(
      dir.exists(
        loaded_results_dir
      ),
      paste0(
        "The input workflow did not create results_dir: ",
        loaded_results_dir
      )
    )
  }
)


## 13. Dynamic baseline summary ------------------------------------------------

loaded_report <- get(
  "report",
  envir = input_environment,
  inherits = FALSE
)

removed_runs <- get(
  "runs_removed_by_global_filter",
  envir = input_environment,
  inherits = FALSE
)

peptide_evidence <- get(
  "peptides_clean",
  envir = input_environment,
  inherits = FALSE
)

parsed_fasta <- get(
  "fasta_ref",
  envir = input_environment,
  inherits = FALSE
)


cat("\n")
cat("============================================================\n")
cat("DYNAMIC INPUT BASELINE SUMMARY\n")
cat("============================================================\n")
cat(
  "Retained DIA-NN rows: ",
  nrow(
    loaded_report
  ),
  "\n",
  sep = ""
)
cat(
  "Retained runs: ",
  dplyr::n_distinct(
    loaded_report$Run
  ),
  "\n",
  sep = ""
)
cat(
  "Runs removed by global filter: ",
  nrow(
    removed_runs
  ),
  "\n",
  sep = ""
)
cat(
  "Protein groups with peptide evidence: ",
  nrow(
    peptide_evidence
  ),
  "\n",
  sep = ""
)
cat(
  "Parsed FASTA entries: ",
  nrow(
    parsed_fasta
  ),
  "\n",
  sep = ""
)
cat(
  "Reference-output files: ",
  nrow(
    reference_manifest_entries
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")


## 14. Consolidated outputs ----------------------------------------------------

result_table <- data.frame(
  Test = vapply(
    results,
    function(result) {
      result$name
    },
    character(1)
  ),
  Success = vapply(
    results,
    function(result) {
      isTRUE(
        result$success
      )
    },
    logical(1)
  ),
  Details = vapply(
    results,
    function(result) {
      result$message
    },
    character(1)
  ),
  Elapsed_Seconds = vapply(
    results,
    function(result) {
      result$elapsed_seconds
    },
    numeric(1)
  ),
  stringsAsFactors = FALSE
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
  "repository-baseline-results.tsv"
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
cat("REPOSITORY AND BASELINE TEST SUMMARY\n")
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