###############################################################################
# File: run_qc.R
#
# Purpose:
#   Provide the scientific, cross-platform entry point for QC-espresso.
#   The runner validates user inputs, creates an isolated execution directory,
#   invokes the scientific Quarto report, and preserves all previous results.
###############################################################################

resolve_qc_project_root <- function() {
  current_directory <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  required_root_files <- c("QC-espresso_DIA.Rproj", "VERSION")
  missing_root_files <- required_root_files[
    !file.exists(file.path(current_directory, required_root_files))
  ]

  if (length(missing_root_files) > 0L) {
    stop(
      paste0(
        "QC-espresso must be executed from the project root.\n",
        "Open `QC-espresso_DIA.Rproj` in RStudio before running the analysis.\n",
        "Current directory:\n", current_directory,
        "\n\nMissing project file(s):\n",
        paste(paste0("  - ", missing_root_files), collapse = "\n")
      ),
      call. = FALSE
    )
  }

  current_directory
}

read_qc_version <- function(project_root) {
  version_file <- file.path(project_root, "VERSION")
  version_value <- trimws(readLines(version_file, n = 1L, warn = FALSE))

  if (length(version_value) != 1L || !nzchar(version_value)) {
    stop(
      paste0("The QC-espresso VERSION file is empty or invalid:\n", version_file),
      call. = FALSE
    )
  }

  version_value
}

is_quarto_executable <- function(path) {
  length(path) == 1L &&
    !is.na(path) &&
    nzchar(path) &&
    file.exists(path) &&
    !dir.exists(path)
}

expand_quarto_candidate <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(character(0))
  }

  path <- path.expand(path)

  if (file.exists(path) && !dir.exists(path)) {
    return(path)
  }

  c(
    file.path(path, "quarto"),
    file.path(path, "quarto.exe"),
    file.path(path, "bin", "quarto"),
    file.path(path, "bin", "quarto.exe")
  )
}

resolve_quarto_executable <- function(project_root = resolve_qc_project_root()) {
  path_quarto <- unname(Sys.which("quarto"))

  ## Prefer the portable Quarto distributed beside the project runtime.
  ## For the repository layout, project_root is QC-espresso_DIA and the
  ## shared portable runtime is located at ../runtime.
  portable_runtime_roots <- unique(c(
    file.path(project_root, "runtime"),
    file.path(project_root, "..", "runtime"),
    file.path(project_root, "R", "runtime"),
    file.path(project_root, "..", "QC-espresso_dev_JPR", "R", "runtime")
  ))

  portable_candidates <- unlist(
    lapply(
      portable_runtime_roots,
      function(runtime_root) {
        c(
          file.path(runtime_root, "quarto", "bin", "quarto.exe"),
          file.path(runtime_root, "quarto", "bin", "quarto"),
          file.path(runtime_root, "quarto", "quarto.exe"),
          file.path(runtime_root, "quarto", "quarto")
        )
      }
    ),
    use.names = FALSE
  )

  environment_candidates <- c(
    Sys.getenv("QUARTO_PATH", unset = ""),
    Sys.getenv("RSTUDIO_QUARTO", unset = "")
  )

  windows_roots <- unique(c(
    Sys.getenv("ProgramFiles", unset = ""),
    Sys.getenv("ProgramW6432", unset = ""),
    Sys.getenv("ProgramFiles(x86)", unset = ""),
    Sys.getenv("LOCALAPPDATA", unset = "")
  ))
  windows_roots <- windows_roots[nzchar(windows_roots)]

  windows_candidates <- character(0)

  if (.Platform$OS.type == "windows") {
    independent_candidates <- unlist(
      lapply(
        windows_roots,
        function(root) {
          c(
            file.path(root, "Quarto", "bin", "quarto.exe"),
            file.path(root, "Programs", "Quarto", "bin", "quarto.exe")
          )
        }
      ),
      use.names = FALSE
    )

    rstudio_installations <- unlist(
      lapply(
        windows_roots,
        function(root) {
          c(
            file.path(root, "RStudio"),
            file.path(root, "Programs", "RStudio")
          )
        }
      ),
      use.names = FALSE
    )

    rstudio_candidates <- unlist(
      lapply(
        rstudio_installations,
        function(root) {
          c(
            file.path(root, "resources", "app", "bin", "quarto", "bin", "quarto.exe"),
            file.path(root, "resources", "app", "quarto", "bin", "quarto.exe"),
            file.path(root, "bin", "quarto", "bin", "quarto.exe")
          )
        }
      ),
      use.names = FALSE
    )

    windows_candidates <- c(independent_candidates, rstudio_candidates)
  }

  unix_candidates <- c(
    "/usr/local/bin/quarto",
    "/usr/bin/quarto",
    "/opt/quarto/bin/quarto",
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto",
    "/Applications/RStudio.app/Contents/Resources/app/bin/quarto/bin/quarto"
  )

  candidates <- unique(c(
    portable_candidates,
    unlist(lapply(environment_candidates, expand_quarto_candidate), use.names = FALSE),
    path_quarto,
    windows_candidates,
    unix_candidates
  ))
  candidates <- candidates[nzchar(candidates)]

  valid_candidates <- candidates[
    vapply(candidates, is_quarto_executable, logical(1))
  ]

  if (length(valid_candidates) > 0L) {
    for (candidate in valid_candidates) {
      normalized_candidate <- normalizePath(
        candidate,
        winslash = "/",
        mustWork = TRUE
      )

      version_status <- suppressWarnings(
        tryCatch(
          system2(
            command = normalized_candidate,
            args = "--version",
            stdout = FALSE,
            stderr = FALSE
          ),
          error = function(error_condition) 1L
        )
      )

      if (identical(as.integer(version_status), 0L)) {
        return(normalized_candidate)
      }
    }
  }

  checked_locations <- paste(
    paste0("  - ", candidates),
    collapse = "\n"
  )

  stop(
    paste0(
      "Quarto could not be located.\n",
      "Install Quarto or use an RStudio installation that includes Quarto.\n",
      "QC-espresso checked the system PATH, QUARTO_PATH, RSTUDIO_QUARTO, ",
      "standard standalone installations, and bundled RStudio locations.\n\n",
      "Checked locations:\n",
      checked_locations
    ),
    call. = FALSE
  )
}

read_quarto_version <- function(quarto_executable) {
  version_output <- tryCatch(
    system2(
      command = quarto_executable,
      args = "--version",
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(error_condition) character(0)
  )

  if (length(version_output) == 0L) {
    return(NA_character_)
  }

  trimws(version_output[[1L]])
}

create_qc_execution_directory <- function(project_root, output_root = "output") {
  if (!grepl("^([A-Za-z]:|/)", output_root)) {
    output_root <- file.path(project_root, output_root)
  }

  output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)

  if (!dir.exists(output_root)) {
    directory_created <- dir.create(
      output_root,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (!isTRUE(directory_created) && !dir.exists(output_root)) {
      stop(
        paste0("The QC-espresso output root could not be created:\n", output_root),
        call. = FALSE
      )
    }
  }

  execution_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  execution_basename <- paste0("QC_espresso_DIA_", execution_timestamp)
  execution_directory <- file.path(output_root, execution_basename)
  numerical_suffix <- 1L

  while (dir.exists(execution_directory)) {
    execution_directory <- file.path(
      output_root,
      paste0(execution_basename, "_", numerical_suffix)
    )
    numerical_suffix <- numerical_suffix + 1L
  }

  creation_success <- dir.create(
    execution_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!isTRUE(creation_success) && !dir.exists(execution_directory)) {
    stop(
      paste0("The execution directory could not be created:\n", execution_directory),
      call. = FALSE
    )
  }

  normalizePath(execution_directory, winslash = "/", mustWork = TRUE)
}

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
    shQuote(paste0("report_path=", report_path)),
    "-P",
    shQuote(paste0("fasta_path=", fasta_path)),
    "-P",
    shQuote(paste0("output_dir=", execution_directory))
  )
}

run_qc <- function(
    report_path,
    fasta_path,
    output_root = "output",
    open_report = interactive()
) {
  execution_started_at <- Sys.time()
  project_root <- resolve_qc_project_root()
  validation_file <- file.path(project_root, "R", "validate_inputs.R")

  if (!file.exists(validation_file)) {
    stop(
      paste0("The QC-espresso validation module was not found:\n", validation_file),
      call. = FALSE
    )
  }

  source(validation_file, local = environment())

  qc_version <- read_qc_version(project_root)
  quarto_executable <- resolve_quarto_executable(project_root)
  quarto_version <- read_quarto_version(quarto_executable)
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
      unlink(execution_directory, recursive = TRUE, force = TRUE)
      stop(conditionMessage(error_condition), call. = FALSE)
    }
  )

  scientific_qmd <- file.path(project_root, "report", "qc_espresso_scientific.qmd")

  if (!file.exists(scientific_qmd)) {
    unlink(execution_directory, recursive = TRUE, force = TRUE)
    stop(
      paste0("The scientific QC-espresso report was not found:\n", scientific_qmd),
      call. = FALSE
    )
  }

  quarto_arguments <- build_quarto_arguments(
    qmd_path = scientific_qmd,
    execution_directory = validated_inputs$output_dir,
    report_path = validated_inputs$report_path,
    fasta_path = validated_inputs$fasta_path
  )

  cat("\n============================================================\n")
  cat("QC-espresso scientific execution\n")
  cat("============================================================\n")
  cat("QC-espresso version: ", qc_version, "\n", sep = "")
  cat("R version: ", R.version.string, "\n", sep = "")
  cat("Quarto executable: ", quarto_executable, "\n", sep = "")
  cat("Quarto version: ", ifelse(is.na(quarto_version), "Not detected", quarto_version), "\n", sep = "")
  cat("DIA-NN report: ", validated_inputs$report_path, "\n", sep = "")
  cat("FASTA reference: ", validated_inputs$fasta_path, "\n", sep = "")
  cat("Output directory: ", validated_inputs$output_dir, "\n", sep = "")
  cat("============================================================\n")
  cat("\nRendering QC report...\n\n")

  ## Ensure that Quarto uses the same R installation and package library as
  ## the current QC-espresso process. QUARTO_R is the authoritative Quarto
  ## setting for the Rscript executable used by the knitr engine.
  current_r_bin <- normalizePath(
    R.home("bin"),
    winslash = "/",
    mustWork = TRUE
  )

  current_path <- Sys.getenv("PATH", unset = "")
  current_r_libs_user <- Sys.getenv("R_LIBS_USER", unset = "")

  render_path <- paste(
    current_r_bin,
    current_path,
    sep = .Platform$path.sep
  )

  current_rscript <- file.path(
    current_r_bin,
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )

  if (!file.exists(current_rscript)) {
    stop(
      paste0(
        "The Rscript executable for the active R session was not found:\n",
        current_rscript
      ),
      call. = FALSE
    )
  }

  previous_environment <- Sys.getenv(
    c(
      "PATH",
      "R_HOME",
      "R_LIBS_USER",
      "QUARTO_R",
      "QUARTO_KNITR_RSCRIPT_ARGS",
      "RENV_CONFIG_AUTOLOADER_ENABLED"
    ),
    unset = NA_character_
  )

  restore_render_environment <- function() {
    for (variable_name in names(previous_environment)) {
      previous_value <- previous_environment[[variable_name]]

      if (is.na(previous_value)) {
        Sys.unsetenv(variable_name)
      } else {
        do.call(
          Sys.setenv,
          stats::setNames(
            list(previous_value),
            variable_name
          )
        )
      }
    }
  }

  on.exit(
    restore_render_environment(),
    add = TRUE
  )

  Sys.setenv(
    PATH = render_path,
    R_HOME = R.home(),
    QUARTO_R = current_rscript,
    QUARTO_KNITR_RSCRIPT_ARGS = "--no-init-file",
    RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE"
  )

  if (nzchar(current_r_libs_user)) {
    Sys.setenv(
      R_LIBS_USER = current_r_libs_user
    )
  }

  cat("Render R executable: ", current_rscript, "\n", sep = "")
  cat(
    "Render R library: ",
    ifelse(
      nzchar(current_r_libs_user),
      current_r_libs_user,
      "R default library paths"
    ),
    "\n",
    sep = ""
  )

  render_output <- system2(
    command = quarto_executable,
    args = quarto_arguments,
    stdout = TRUE,
    stderr = TRUE
  )

  restore_render_environment()

  render_status <- attr(render_output, "status")
  if (is.null(render_status)) render_status <- 0L

  if (length(render_output) > 0L) {
    cat(paste(render_output, collapse = "\n"), "\n")
  }

  if (!identical(as.integer(render_status), 0L)) {
    stop(
      paste0(
        "Quarto rendering failed with exit code ", render_status, ".\n",
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
        "Quarto completed without a reported error, but the expected HTML ",
        "report was not found:\n",
        rendered_html
      ),
      call. = FALSE
    )
  }

  execution_completed_at <- Sys.time()
  elapsed_seconds <- as.numeric(
    difftime(execution_completed_at, execution_started_at, units = "secs")
  )

  generated_files <- list.files(
    validated_inputs$output_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )

  cat("\n============================================================\n")
  cat("QC-espresso completed successfully\n")
  cat("============================================================\n")
  cat("HTML report: ", rendered_html, "\n", sep = "")
  cat("Generated files: ", length(generated_files), "\n", sep = "")
  cat("Elapsed time: ", round(elapsed_seconds, 1), " seconds\n", sep = "")
  cat("============================================================\n")

  if (isTRUE(open_report) && interactive()) {
    utils::browseURL(rendered_html)
  }

  result <- list(
    success = TRUE,
    qc_espresso_version = qc_version,
    r_version = R.version.string,
    quarto_executable = quarto_executable,
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

  class(result) <- c("qc_espresso_result", class(result))
  invisible(result)
}

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
