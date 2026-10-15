###############################################################################
# File: qc_session.R
#
# Purpose:
#   Record and report execution-environment details used to generate the QC
#   report, including the execution timestamp, R and Quarto versions, and
#   loaded R package versions.
#
# Inputs (expected to exist before sourcing):
#   - params : Quarto parameters list, optionally containing quarto_exec
#
# Outputs (objects defined):
#   - session_date   : Date and time at which session metadata were captured
#   - r_version      : R version used to generate the report
#   - quarto_exec    : Quarto executable path or command obtained from params
#   - quarto_valid   : Availability flag for the specified Quarto executable
#   - quarto_version : Detected Quarto version or NA when unavailable
#   - pkg_info       : Loaded R package information returned by sessionInfo()
#   - packages       : Loaded R package names and versions
#   - session_summary: Execution-environment summary table
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tibble
#   - gt
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - The params object is optional; quarto_exec remains NULL when params or
#     params$quarto_exec is unavailable.
#   - The Quarto executable is validated as either an existing file path or a
#     command available through the system PATH.
#   - Failure to detect or query Quarto does not interrupt report generation.
#   - Only packages recorded in sessionInfo()$otherPkgs are included in the
#     loaded-package table.
#   - The

## 1. Capture execution metadata -----------------------------------------------

session_date <- Sys.time()
r_version    <- R.version.string

## Safe extraction of quarto_exec from params

quarto_exec <- NULL

if (exists("params") && !is.null(params$quarto_exec)) {
  quarto_exec <- params$quarto_exec
}

## Validate quarto_exec

quarto_valid <- FALSE

if (!is.null(quarto_exec) &&
    is.character(quarto_exec) &&
    length(quarto_exec) == 1 &&
    nzchar(quarto_exec)) {
  
  # check if executable exists (path or in PATH)
  if (file.exists(quarto_exec) || nzchar(Sys.which(quarto_exec))) {
    quarto_valid <- TRUE
  }
}

## Retrieve Quarto version safely

quarto_version <- NA_character_

if (quarto_valid) {
  quarto_version <- tryCatch(
    {
      system2(quarto_exec, "--version", stdout = TRUE, stderr = TRUE)[1]
    },
    warning = function(w) NA_character_,
    error   = function(e) NA_character_
  )
}


## 2. Loaded R packages ---------------------------------------------------------

pkg_info <- sessionInfo()$otherPkgs

if (is.null(pkg_info)) {
  packages <- tibble::tibble(
    Package = character(0),
    Version = character(0)
  )
} else {
  packages <- tibble::tibble(
    Package = names(pkg_info),
    Version = vapply(
      pkg_info,
      function(x) x$Version,
      character(1)
    )
  )
}


## 3. Execution environment summary --------------------------------------------

session_summary <- tibble::tibble(
  Item = c(
    "Report generated on",
    "R version",
    "Quarto version",
    "Number of loaded packages"
  ),
  Value = c(
    format(session_date, "%Y-%m-%d %H:%M:%S"),
    r_version,
    ifelse(is.na(quarto_version), "Not detected", quarto_version),
    nrow(packages)
  )
)

print(
  session_summary %>% 
    gt::gt() %>% 
    gt::tab_header(
      title = "Execution environment summary"
    ) %>% 
    gt::cols_label(
      Item  = "Item",
      Value = "Value"
    ) %>% 
    gt::tab_options(
      table.width = gt::pct(100)
    )
)


## 4. Loaded packages and versions ---------------------------------------------

print(
  packages %>% 
    gt::gt() %>% 
    gt::tab_header(
      title = "Loaded R packages and versions"
    ) %>% 
    gt::tab_options(
      table.width          = gt::pct(100),
      container.height     = gt::px(450),
      container.overflow.y = "scroll"
    )
)