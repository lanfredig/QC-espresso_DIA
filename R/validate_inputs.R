###############################################################################
# File: validate_inputs.R
#
# Purpose:
#   Validate the scientific QC-espresso execution environment, input files,
#   DIA-NN report schema, FASTA structure, and output directory before the
#   analytical workflow is started.
#
# Notes:
#   - This file defines validation functions only.
#   - No validation is performed automatically when this file is sourced.
#   - No analytical result is calculated or modified by these functions.
###############################################################################


## 1. Required package definitions --------------------------------------------

qc_required_packages <- function() {
  c(
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
}


## 2. Required DIA-NN report columns ------------------------------------------

qc_required_report_columns <- function() {
  c(
    "Run",
    "Decoy",
    "PG.Q.Value",
    "Protein.Group",
    "Precursor.Id",
    "Modified.Sequence",
    "Stripped.Sequence",
    "Precursor.Charge",
    "Precursor.Quantity",
    "PG.MaxLFQ",
    "Quantity.Quality",
    "Empirical.Quality",
    "PG.MaxLFQ.Quality",
    "RT",
    "RT.Start",
    "RT.Stop",
    "FWHM",
    "Precursor.Mz",
    "Ms1.Apex.Mz.Delta"
  )
}


## 3. Optional DIA-NN report columns ------------------------------------------

qc_optional_report_columns <- function() {
  c(
    "IM",
    "Ion.Mobility",
    "Precursor.Ion.Mobility",
    "1/K0",
    "iIM"
  )
}


## 4. Scalar path validation ---------------------------------------------------

validate_scalar_path <- function(
    path,
    argument_name
) {
  if (
    is.null(path) ||
    length(path) != 1L ||
    is.na(path) ||
    !is.character(path) ||
    !nzchar(trimws(path))
  ) {
    stop(
      paste0(
        "`",
        argument_name,
        "` must contain exactly one valid file path."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


## 5. DIA-NN Parquet file validation ------------------------------------------

validate_report_file <- function(report_path) {
  validate_scalar_path(
    path = report_path,
    argument_name = "report_path"
  )
  
  normalized_path <- normalizePath(
    report_path,
    winslash = "/",
    mustWork = FALSE
  )
  
  if (!file.exists(normalized_path)) {
    stop(
      paste0(
        "The DIA-NN report was not found:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  if (dir.exists(normalized_path)) {
    stop(
      paste0(
        "The supplied DIA-NN report path refers to a directory, not a file:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  extension <- tolower(
    tools::file_ext(normalized_path)
  )
  
  if (!identical(extension, "parquet")) {
    stop(
      paste0(
        "The DIA-NN report must have the .parquet extension:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  file_size <- file.info(normalized_path)$size
  
  if (
    is.na(file_size) ||
    file_size <= 0
  ) {
    stop(
      paste0(
        "The DIA-NN Parquet file is empty or inaccessible:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  normalized_path
}


## 6. FASTA file validation ----------------------------------------------------

validate_fasta_file <- function(fasta_path) {
  validate_scalar_path(
    path = fasta_path,
    argument_name = "fasta_path"
  )
  
  normalized_path <- normalizePath(
    fasta_path,
    winslash = "/",
    mustWork = FALSE
  )
  
  if (!file.exists(normalized_path)) {
    stop(
      paste0(
        "The FASTA reference file was not found:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  if (dir.exists(normalized_path)) {
    stop(
      paste0(
        "The supplied FASTA path refers to a directory, not a file:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  extension <- tolower(
    tools::file_ext(normalized_path)
  )
  
  allowed_extensions <- c(
    "fasta",
    "fa",
    "faa"
  )
  
  if (!extension %in% allowed_extensions) {
    stop(
      paste0(
        "The reference database must use one of these extensions: ",
        paste(allowed_extensions, collapse = ", "),
        ".\nReceived:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  file_size <- file.info(normalized_path)$size
  
  if (
    is.na(file_size) ||
    file_size <= 0
  ) {
    stop(
      paste0(
        "The FASTA file is empty or inaccessible:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  preview_lines <- readLines(
    normalized_path,
    n = 100L,
    warn = FALSE
  )
  
  preview_lines <- trimws(
    preview_lines
  )
  
  preview_lines <- preview_lines[
    nzchar(preview_lines)
  ]
  
  if (length(preview_lines) == 0L) {
    stop(
      paste0(
        "No readable content was found in the FASTA file:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  if (!any(startsWith(preview_lines, ">"))) {
    stop(
      paste0(
        "The file does not contain a recognizable FASTA header ",
        "within the first 100 lines:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  normalized_path
}


## 7. Required package validation ---------------------------------------------

validate_required_packages <- function(
    packages = qc_required_packages()
) {
  package_status <- vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
  
  missing_packages <- packages[
    !package_status
  ]
  
  if (length(missing_packages) > 0L) {
    stop(
      paste0(
        "The QC-espresso R environment is incomplete.\n",
        "Missing package(s):\n",
        paste(
          paste0("  - ", missing_packages),
          collapse = "\n"
        ),
        "\n\nRestore the project environment with `renv::restore()` ",
        "before running the analysis."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


## 8. DIA-NN report schema validation -----------------------------------------

validate_report_schema <- function(
    report_path,
    required_columns = qc_required_report_columns()
) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(
      paste0(
        "The `arrow` package is required to inspect the DIA-NN ",
        "Parquet report."
      ),
      call. = FALSE
    )
  }
  
  report_dataset <- tryCatch(
    arrow::open_dataset(
      report_path,
      format = "parquet"
    ),
    error = function(error_condition) {
      stop(
        paste0(
          "The DIA-NN Parquet report could not be opened.\n",
          "File: ",
          report_path,
          "\nReason: ",
          conditionMessage(error_condition)
        ),
        call. = FALSE
      )
    }
  )
  
  available_columns <- names(
    report_dataset
  )
  
  missing_columns <- setdiff(
    required_columns,
    available_columns
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "The DIA-NN report is missing required column(s):\n",
        paste(
          paste0("  - ", missing_columns),
          collapse = "\n"
        ),
        "\n\nAvailable columns were read successfully, but the report ",
        "does not contain all fields required by the current QC workflow."
      ),
      call. = FALSE
    )
  }
  
  optional_columns <- qc_optional_report_columns()
  
  available_optional_columns <- intersect(
    optional_columns,
    available_columns
  )
  
  list(
    available_columns = available_columns,
    required_columns = required_columns,
    missing_columns = missing_columns,
    available_optional_columns = available_optional_columns
  )
}


## 9. Output-directory validation ---------------------------------------------

validate_output_directory <- function(output_dir) {
  validate_scalar_path(
    path = output_dir,
    argument_name = "output_dir"
  )
  
  normalized_path <- normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )
  
  if (!dir.exists(normalized_path)) {
    directory_created <- dir.create(
      normalized_path,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (
      !isTRUE(directory_created) &&
      !dir.exists(normalized_path)
    ) {
      stop(
        paste0(
          "The output directory could not be created:\n",
          normalized_path
        ),
        call. = FALSE
      )
    }
  }
  
  test_file <- file.path(
    normalized_path,
    paste0(
      ".qc_espresso_write_test_",
      Sys.getpid(),
      ".tmp"
    )
  )
  
  write_success <- tryCatch(
    {
      writeLines(
        "QC-espresso output-directory write test",
        con = test_file,
        useBytes = TRUE
      )
      
      TRUE
    },
    error = function(error_condition) {
      FALSE
    }
  )
  
  if (!isTRUE(write_success)) {
    stop(
      paste0(
        "The output directory is not writable:\n",
        normalized_path
      ),
      call. = FALSE
    )
  }
  
  if (file.exists(test_file)) {
    unlink(
      test_file,
      force = TRUE
    )
  }
  
  normalized_path
}


## 10. Combined scientific input validation -----------------------------------

validate_qc_inputs <- function(
    report_path,
    fasta_path,
    output_dir,
    validate_packages = TRUE,
    validate_schema = TRUE
) {
  if (isTRUE(validate_packages)) {
    validate_required_packages()
  }
  
  validated_report_path <- validate_report_file(
    report_path
  )
  
  validated_fasta_path <- validate_fasta_file(
    fasta_path
  )
  
  validated_output_dir <- validate_output_directory(
    output_dir
  )
  
  schema_result <- NULL
  
  if (isTRUE(validate_schema)) {
    schema_result <- validate_report_schema(
      validated_report_path
    )
  }
  
  result <- list(
    report_path = validated_report_path,
    fasta_path = validated_fasta_path,
    output_dir = validated_output_dir,
    report_schema = schema_result
  )
  
  class(result) <- c(
    "qc_espresso_validated_inputs",
    class(result)
  )
  
  result
}