###############################################################################
# File: qc_completeness.R
#
# Purpose:
#   Assess precursor-level data completeness by constructing a run × precursor
#   matrix and summarizing both global and per-run completeness metrics.
#
# Inputs (expected to exist before sourcing):
#   - report : Filtered DIA-NN report data frame
#
# Outputs (objects defined):
#   - data_complt             : Run × precursor quantitative matrix
#   - completeness_global_df  : Global data completeness summary table
#   - completeness_per_run_df : Per-run data completeness summary table
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tidyr
#   - tibble
#
# Notes:
#   - Completeness is computed using presence/absence of quantitative values
#     only and is independent of downstream normalization or filtering steps.
#   - This script performs no file I/O and does not generate plots directly.
###############################################################################


## 1. Construct run × precursor matrix -----------------------------------------

data_complt <- report %>% 
  dplyr::select(
    Run,
    Precursor.Id,
    Precursor.Quantity
  ) %>% 
  dplyr::distinct() %>% 
  tidyr::pivot_wider(
    names_from  = Precursor.Id,
    values_from = Precursor.Quantity
  )


## 2. Global data completeness --------------------------------------------------

n_runs       <- dplyr::n_distinct(data_complt$Run)
n_precursors <- ncol(data_complt) - 1

total_possible <- n_runs * n_precursors
observed       <- sum(!is.na(data_complt[-1]))

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
    (observed / total_possible) * 100
  )
)


## 3. Run-level data completeness -----------------------------------------------

completeness_per_run_df <- data_complt %>% 
  dplyr::mutate(
    Observed         = rowSums(!is.na(dplyr::across(-Run))),
    Total            = n_precursors,
    Completeness_Pct = (Observed / Total) * 100
  ) %>% 
  dplyr::select(
    Run,
    Observed,
    Total,
    Completeness_Pct
  )


## 4. Clean exit ---------------------------------------------------------------

invisible(NULL)