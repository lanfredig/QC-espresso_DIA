###############################################################################
# File: qc_completeness.R
#
# Purpose:
#   Assess precursor-level quantitative data completeness without materializing
#   a dense run-by-precursor matrix.
#
# Inputs (expected to exist before sourcing):
#   - report : Filtered DIA-NN report data frame
#
# Outputs (objects defined):
#   - completeness_input      : Valid run, precursor, and quantitative entries
#   - run_reference           : Distinct runs included in the assessment
#   - precursor_reference     : Distinct precursor identifiers included in the assessment
#   - n_runs                  : Number of retained runs
#   - n_precursors            : Number of retained precursor identifiers
#   - total_possible          : Number of possible run–precursor combinations
#   - data_complt             : Long-format observation status per run and precursor
#   - observed                : Number of observed run–precursor combinations
#   - completeness_global_pct : Global quantitative data completeness percentage
#   - completeness_global_df  : Global data completeness summary table
#   - completeness_per_run_df : Per-run data completeness summary table
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tibble
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - Completeness is based on the presence of non-missing
#     Precursor.Quantity values.
#   - The denominator is the Cartesian space defined by all retained runs and
#     all retained precursor identifiers.
#   - Duplicate Run and Precursor.Id combinations are collapsed before
#     counting.
#   - A run–precursor combination is considered observed when at least one
#     corresponding Precursor.Quantity value is non-missing.
#   - No dense run-by-precursor matrix is created.
###############################################################################
## 1. Validate identifiers -----------------------------------------------------

completeness_input <- report %>%
  dplyr::select(
    Run,
    Precursor.Id,
    Precursor.Quantity
  ) %>%
  dplyr::filter(
    !is.na(Run),
    nzchar(trimws(as.character(Run))),
    !is.na(Precursor.Id),
    nzchar(trimws(as.character(Precursor.Id)))
  )

run_reference <- completeness_input %>%
  dplyr::distinct(Run)

precursor_reference <- completeness_input %>%
  dplyr::distinct(Precursor.Id)

n_runs <- nrow(run_reference)
n_precursors <- nrow(precursor_reference)
total_possible <- n_runs * n_precursors

## 2. Collapse to one observation per Run + Precursor.Id -----------------------

## This long-format object replaces the previous dense matrix. A pair is
## considered observed when at least one corresponding quantitative value is
## non-missing.
data_complt <- completeness_input %>%
  dplyr::group_by(
    Run,
    Precursor.Id
  ) %>%
  dplyr::summarise(
    Observed = any(!is.na(Precursor.Quantity)),
    .groups = "drop"
  )

observed <- sum(data_complt$Observed)

## 3. Global data completeness -------------------------------------------------

completeness_global_pct <- if (total_possible > 0) {
  100 * observed / total_possible
} else {
  NA_real_
}

completeness_global_df <- tibble::tibble(
  Metric = c(
    "Number of runs",
    "Number of precursors",
    "Total possible values (Run × Precursor)",
    "Observed values",
    "Data completeness (%)"
  ),
  Value = c(
    n_runs,
    n_precursors,
    total_possible,
    observed,
    completeness_global_pct
  )
)

## 4. Run-level data completeness ----------------------------------------------

completeness_per_run_df <- run_reference %>%
  dplyr::left_join(
    data_complt %>%
      dplyr::group_by(Run) %>%
      dplyr::summarise(
        Observed = sum(Observed),
        .groups = "drop"
      ),
    by = "Run"
  ) %>%
  dplyr::mutate(
    Observed = dplyr::coalesce(Observed, 0L),
    Total = n_precursors,
    Completeness_Pct = dplyr::if_else(
      Total > 0,
      100 * Observed / Total,
      NA_real_
    )
  ) %>%
  dplyr::select(
    Run,
    Observed,
    Total,
    Completeness_Pct
  )

## 5. Clean exit ---------------------------------------------------------------

invisible(NULL)
