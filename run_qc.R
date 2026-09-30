###############################################################################
# File: run_qc.R
#
# Purpose:
#   Provide the scientific, cross-platform entry point for QC-espresso.
#   The runner validates user inputs, creates an isolated execution directory,
#   invokes the scientific Quarto report, and preserves all previous results.
#
# Usage:
#   source("run_qc.R")
#
#   result <- run_qc(
#     report_path = "example-data/report.parquet",
#     fasta_path = "example-data/reference.fasta"
#   )
#
# Notes:
#   - The RStudio project must be opened before this file is sourced.
#   - The Quarto CLI must be available through RStudio or the system PATH.
#   - Project dependencies must be restored with renv::restore().
###############################################################################


## 1. Project-root resolution --------------------------------------------------

resolve_qc_project_root <- function() {
  current_directory <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
  
  required_root_files <- c(
    "QC-espresso_DIA.Rproj",
    "VERSION"
  )
  
  missing_root_files <- required_root_files[
    !file.exists(
      file.path(
        current_directory,
        required_root_files
      )
    )
  ]
  
  if (length(missing_root_files) > 0L) {
    stop(
      paste0(
        "QC-espresso must be executed from the project root.\n",
        "Open `QC-espresso_DIA.Rproj` in RStudio before running the analysis.\n",
        "Current directory:\n",
        current_directory,
        "\n\nMissing project file(s):\n",
        paste(
          paste0("  - ", missing_root_files),
          collapse = "\n"
        )
      ),
      call. = FALSE
    )
  }
  
  current_directory
}


## 2. Software-version resolution ---------------------------------------------

read_qc_version <- function(project_root) {
  version_file <- file.path(
    project_root,
    "VERSION"
  )
  
  version_value <- readLines(
    version_file,
    n = 1L,
    warn = FALSE
  )
  
  version_value <- trimws(
    version_value
  )
  
  if (
    length(version_value) != 1L ||
    !nzchar(version_value)
  ) {
    stop(
      paste0(
        "The QC-espresso VERSION file is empty or invalid:\n",
        version_file
      ),
      call. = FALSE
    )
  }
  
  version_value
}


## 3. Quarto executable resolution --------------------------------------------

resolve_quarto_executable <- function() {
  quarto_path <- Sys.which(
    "quarto"
  )
  
  if (nzchar(quarto_path)) {
    return(
      normalizePath(
        quarto_path,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }
  
  rstudio_quarto <- Sys.getenv(
    "RSTUDIO_QUARTO",
    unset = ""
  )
  
  if (nzchar(rstudio_quarto)) {
    candidates <- c(
      rstudio_quarto,
      file.path(
        rstudio_quarto,
        "bin",
        "quarto"
      ),
      file.path(
        rstudio_quarto,
        "bin",
        "quarto.exe"
      ),
      file.path(
        rstudio_quarto,
        "quarto"
      ),
      file.path(
        rstudio_quarto,
        "quarto.exe"
      )
    )
    
    candidates <- unique(
      candidates[
        file.exists(candidates)
      ]
    )
    
    if (length(candidates) > 0L) {
      return(
        normalizePath(
          candidates[[1L]],
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }
  
  stop(
    paste0(
      "Quarto could not be located.\n",
      "Install Quarto or use an RStudio installation that includes Quarto.\n",
      "The `quarto --version` command must be available before running ",
      "QC-espresso."
    ),
    call. = FALSE
  )
}


## 4. Quarto-version resolution -----------------------------------------------

read_quarto_version <- function(quarto_executable) {
  version_output <- tryCatch(
    system2(
      command = quarto_executable,
      args = "--version",
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(error_condition) {
      character(0)
    }
  )
  
  if (length(version_output) == 0L) {
    return(NA_character_)
  }
  
  trimws(
    version_output[[1L]]
  )
}


## 5. Execution-directory creation --------------------------------------------

create_qc_execution_directory <- function(
    project_root,
    output_root = "output"
) {
  if (!grepl("^([A-Za-z]:|/)", output_root)) {
    output_root <- file.path(
      project_root,
      output_root
    )
  }
  
  output_root <- normalizePath(
    output_root,
    winslash = "/",
    mustWork = FALSE
  )
  
  if (!dir.exists(output_root)) {
    directory_created <- dir.create(
      output_root,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (
      !isTRUE(directory_created) &&
      !dir.exists(output_root)
    ) {
      stop(
        paste0(
          "The QC-espresso output root could not be created:\n",
          output_root
        ),
        call. = FALSE
      )
    }
  }
  
  execution_timestamp <- format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )
  
  execution_basename <- paste0(
    "QC_espresso_DIA_",
    execution_timestamp
  )
  
  execution_directory <- file.path(
    output_root,
    execution_basename
  )
  
  numerical_suffix <- 1L
  
  while (dir.exists(execution_directory)) {
    execution_directory <- file.path(
      output_root,
      paste0(
        execution_basename,
        "_",
        numerical_suffix
      )
    )
    
    numerical_suffix <- numerical_suffix + 1L
  }
  
  creation_success <- dir.create(
    execution_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (
    !isTRUE(creation_success) &&
    !dir.exists(execution_directory)
  ) {
    stop(
      paste0(
        "The execution directory could not be created:\n",
        execution_directory
      ),
      call. = FALSE
    )
  }
  
  normalizePath(
    execution_directory,
    winslash = "/",
    mustWork = TRUE
  )
}


## 6. Quarto argument construction --------------------------------------------

build_quarto_arguments <- function(
    qmd_path,
    execution_directory,
    report_path,
    fasta_path
) {
  c(
    "render",
    shQuote(qmd_path),
    "--to",
    "html",
    "--output-dir",
    shQuote(execution_directory),
    "-P",
    shQuote(
      paste0(
        "report_path=",
        report_path
      )
    ),
    "-P",
    shQuote(
      paste0(
        "fasta_path=",
        fasta_path
      )
    ),
    "-P",
    shQuote(
      paste0(
        "output_dir=",
        execution_directory
      )
    )
  )
}


## 7. Main scientific runner ---------------------------------------------------

run_qc <- function(
    report_path,
    fasta_path,
    output_root = "output",
    open_report = interactive()
) {
  execution_started_at <- Sys.time()
  
  project_root <- resolve_qc_project_root()
  
  validation_file <- file.path(
    project_root,
    "R",
    "validate_inputs.R"
  )
  
  if (!file.exists(validation_file)) {
    stop(
      paste0(
        "The QC-espresso validation module was not found:\n",
        validation_file
      ),
      call. = FALSE
    )
  }
  
  source(
    validation_file,
    local = environment()
  )
  
  qc_version <- read_qc_version(
    project_root
  )
  
  quarto_executable <- resolve_quarto_executable()
  
  quarto_version <- read_quarto_version(
    quarto_executable
  )
  
  execution_directory <- create_qc_execution_directory(
    project_root = project_root,
    output_root = output_root
  )
  
  validated_inputs <- tryCatch(
    validate_qc_inputs(
      report_path = report_path,
      fasta_path = fasta_path,
      output_dir = execution_directory,
      validate_packages = TRUE,
      validate_schema = TRUE
    ),
    error = function(error_condition) {
      unlink(
        execution_directory,
        recursive = TRUE,
        force = TRUE
      )
      
      stop(
        conditionMessage(error_condition),
        call. = FALSE
      )
    }
  )
  
  scientific_qmd <- file.path(
    project_root,
    "report",
    "qc_espresso_scientific.qmd"
  )
  
  if (!file.exists(scientific_qmd)) {
    unlink(
      execution_directory,
      recursive = TRUE,
      force = TRUE
    )
    
    stop(
      paste0(
        "The scientific QC-espresso report was not found:\n",
        scientific_qmd
      ),
      call. = FALSE
    )
  }
  
  quarto_arguments <- build_quarto_arguments(
    qmd_path = scientific_qmd,
    execution_directory = validated_inputs$output_dir,
    report_path = validated_inputs$report_path,
    fasta_path = validated_inputs$fasta_path
  )
  
  cat("\n")
  cat("============================================================\n")
  cat("QC-espresso scientific execution\n")
  cat("============================================================\n")
  cat("QC-espresso version: ", qc_version, "\n", sep = "")
  cat("R version: ", R.version.string, "\n", sep = "")
  cat(
    "Quarto version: ",
    ifelse(
      is.na(quarto_version),
      "Not detected",
      quarto_version
    ),
    "\n",
    sep = ""
  )
  cat("DIA-NN report: ", validated_inputs$report_path, "\n", sep = "")
  cat("FASTA reference: ", validated_inputs$fasta_path, "\n", sep = "")
  cat("Output directory: ", validated_inputs$output_dir, "\n", sep = "")
  cat("============================================================\n")
  cat("\nRendering QC report...\n\n")
  
  render_output <- system2(
    command = quarto_executable,
    args = quarto_arguments,
    stdout = TRUE,
    stderr = TRUE
  )
  
  render_status <- attr(
    render_output,
    "status"
  )
  
  if (is.null(render_status)) {
    render_status <- 0L
  }
  
  if (length(render_output) > 0L) {
    cat(
      paste(
        render_output,
        collapse = "\n"
      ),
      "\n"
    )
  }
  
  if (!identical(as.integer(render_status), 0L)) {
    stop(
      paste0(
        "Quarto rendering failed with exit code ",
        render_status,
        ".\n",
        "The execution directory was preserved for inspection:\n",
        validated_inputs$output_dir
      ),
      call. = FALSE
    )
  }
  
  rendered_html <- file.path(
    validated_inputs$output_dir,
    "qc_espresso_scientific.html"
  )
  
  if (!file.exists(rendered_html)) {
    stop(
      paste0(
        "Quarto completed without a reported error, but the expected ",
        "HTML report was not found:\n",
        rendered_html
      ),
      call. = FALSE
    )
  }
  
  execution_completed_at <- Sys.time()
  
  elapsed_seconds <- as.numeric(
    difftime(
      execution_completed_at,
      execution_started_at,
      units = "secs"
    )
  )
  
  generated_files <- list.files(
    validated_inputs$output_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )
  
  cat("\n")
  cat("============================================================\n")
  cat("QC-espresso completed successfully\n")
  cat("============================================================\n")
  cat("HTML report: ", rendered_html, "\n", sep = "")
  cat("Generated files: ", length(generated_files), "\n", sep = "")
  cat(
    "Elapsed time: ",
    round(elapsed_seconds, 1),
    " seconds\n",
    sep = ""
  )
  cat("============================================================\n")
  
  if (
    isTRUE(open_report) &&
    interactive()
  ) {
    utils::browseURL(
      rendered_html
    )
  }
  
  result <- list(
    success = TRUE,
    qc_espresso_version = qc_version,
    r_version = R.version.string,
    quarto_version = quarto_version,
    report_path = validated_inputs$report_path,
    fasta_path = validated_inputs$fasta_path,
    output_directory = validated_inputs$output_dir,
    html_report = rendered_html,
    generated_files = generated_files,
    started_at = execution_started_at,
    completed_at = execution_completed_at,
    elapsed_seconds = elapsed_seconds
  )
  
  class(result) <- c(
    "qc_espresso_result",
    class(result)
  )
  
  invisible(result)
}


## 8. Runner load message ------------------------------------------------------

cat(
  paste0(
    "\nQC-espresso scientific runner loaded.\n\n",
    "Example execution:\n\n",
    "result <- run_qc(\n",
    "  report_path = \"example-data/report.parquet\",\n",
    "  fasta_path = \"example-data/reference.fasta\"\n",
    ")\n"
  )
)