###############################################################################
# File: qc_figures_of_merit.R
#
# Purpose:
#   Generate a consolidated run-level QC summary table (Figures of Merit),
#   export the table as CSV, and expose it for rendering in the QMD report.
#
# Inputs (expected to exist before sourcing):
#   - report      : Filtered DIA-NN report data frame
#   - results_dir : Canonical output directory
#
# Outputs:
#   - figures_of_merit_df
#   - figures_of_merit_csv
#
# Exported metrics:
#   - Precursors
#   - Peptides
#   - Protein_Groups
#   - Median_FWHM_s
#   - Median_PeakBase_s
#   - Median_RT_Deviation_min
#   - Median_Mass_Error_ppm
#   - Median_Peptide_Length
#   - Missed_Cleavage_Pct
#   - Charge2_Pct
#   - Charge3_Pct
#   - Contaminant_Pct
#
# Notes:
#   - This script does not compute data completeness.
#   - This script does not include intensity-derived metrics.
#   - This script does not include DIA-NN quality metrics.
#   - This script does not include internal standard signal metrics.
#   - Retention-time reproducibility is computed globally using shared
#     precursors observed in at least 50% of runs.
###############################################################################


## 1. Input validation ---------------------------------------------------------

fom_available <- TRUE

figures_of_merit_df  <- NULL
figures_of_merit_csv <- NULL

if (!exists("report")) {
  cat("[FOM] Skipping: report not found\n")
  fom_available <- FALSE
}

if (!exists("results_dir")) {
  cat("[FOM] Skipping: results_dir not found\n")
  fom_available <- FALSE
}


## 2. Compute Figures of Merit -------------------------------------------------

if (isTRUE(fom_available)) {
  
  ## --------------------------------------------------------------------------
  ## Run reference
  ## --------------------------------------------------------------------------
  
  run_reference <- report %>%
    dplyr::distinct(Run) %>%
    dplyr::arrange(Run)
  
  
  ## --------------------------------------------------------------------------
  ## Identification depth
  ## --------------------------------------------------------------------------
  
  identification_df <- report %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Precursors     = dplyr::n_distinct(Precursor.Id),
      Peptides       = dplyr::n_distinct(Modified.Sequence),
      Protein_Groups = dplyr::n_distinct(Protein.Group),
      .groups = "drop"
    )
  
  
  ## --------------------------------------------------------------------------
  ## Chromatographic metrics
  ## --------------------------------------------------------------------------
  
  chromatography_df <- report %>%
    dplyr::transmute(
      Run,
      FWHM_s = FWHM * 60,
      PeakBase_s = (RT.Stop - RT.Start) * 60
    ) %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Median_FWHM_s = median(FWHM_s, na.rm = TRUE),
      Median_PeakBase_s = median(PeakBase_s, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  ## --------------------------------------------------------------------------
  ## Global retention-time reproducibility
  ## --------------------------------------------------------------------------
  ## Strategy:
  ##   1. Collapse RT to one value per Run + Precursor.Id.
  ##   2. Keep precursors observed in at least 50% of runs.
  ##   3. Compute the global median RT per precursor.
  ##   4. Compute absolute RT deviation from the global median.
  ##   5. Summarise the median absolute RT deviation per run.
  ## --------------------------------------------------------------------------
  
  rt_reference <- report %>%
    dplyr::select(
      Run,
      Precursor.Id,
      RT
    ) %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Precursor.Id),
      !is.na(RT)
    ) %>%
    dplyr::group_by(
      Run,
      Precursor.Id
    ) %>%
    dplyr::summarise(
      RT_Run = median(RT, na.rm = TRUE),
      .groups = "drop"
    )
  
  n_runs_rt <- dplyr::n_distinct(rt_reference$Run)
  
  min_runs_required <- max(
    2,
    ceiling(0.5 * n_runs_rt)
  )
  
  shared_precursors <- rt_reference %>%
    dplyr::group_by(Precursor.Id) %>%
    dplyr::summarise(
      n_runs = dplyr::n_distinct(Run),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      n_runs >= min_runs_required
    ) %>%
    dplyr::pull(Precursor.Id)
  
  if (length(shared_precursors) > 0) {
    
    rt_df <- rt_reference %>%
      dplyr::filter(
        Precursor.Id %in% shared_precursors
      ) %>%
      dplyr::group_by(Precursor.Id) %>%
      dplyr::mutate(
        RT_Global_Median = median(RT_Run, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        RT_Deviation = abs(RT_Run - RT_Global_Median)
      ) %>%
      dplyr::group_by(Run) %>%
      dplyr::summarise(
        Median_RT_Deviation_min = median(
          RT_Deviation,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
  } else {
    
    rt_df <- run_reference %>%
      dplyr::mutate(
        Median_RT_Deviation_min = NA_real_
      )
    
    message(
      "[FOM] Median_RT_Deviation_min could not be computed: ",
      "no shared precursors passed the minimum run-frequency threshold"
    )
  }
  
  
  ## --------------------------------------------------------------------------
  ## Mass accuracy
  ## --------------------------------------------------------------------------
  
  mass_df <- report %>%
    dplyr::filter(
      !is.na(Precursor.Quantity),
      Precursor.Quantity > 0,
      !is.na(Precursor.Mz),
      Precursor.Mz > 0,
      !is.na(Ms1.Apex.Mz.Delta)
    ) %>%
    dplyr::mutate(
      Mass_Error_PPM =
        (Ms1.Apex.Mz.Delta / Precursor.Mz) * 1e6
    ) %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Median_Mass_Error_ppm = median(
        Mass_Error_PPM,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  ## --------------------------------------------------------------------------
  ## Peptide metrics
  ## --------------------------------------------------------------------------
  
  peptide_df <- report %>%
    dplyr::select(
      Run,
      Stripped.Sequence
    ) %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Stripped.Sequence),
      nzchar(Stripped.Sequence)
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Peptide_Length = nchar(Stripped.Sequence),
      MissedCleavages = pmax(
        0,
        stringr::str_count(
          Stripped.Sequence,
          "(?<=[KR])(?!P)"
        ) - 1
      )
    ) %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Median_Peptide_Length = median(
        Peptide_Length,
        na.rm = TRUE
      ),
      Missed_Cleavage_Pct =
        100 * mean(MissedCleavages > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  ## --------------------------------------------------------------------------
  ## Charge-state distribution
  ## --------------------------------------------------------------------------
  
  charge_df <- report %>%
    dplyr::select(
      Run,
      Precursor.Id,
      Precursor.Charge
    ) %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Precursor.Id),
      !is.na(Precursor.Charge)
    ) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Charge2_Pct =
        100 * mean(Precursor.Charge == 2, na.rm = TRUE),
      Charge3_Pct =
        100 * mean(Precursor.Charge == 3, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  ## --------------------------------------------------------------------------
  ## Contaminant contribution
  ## --------------------------------------------------------------------------
  
  contaminant_df <- report %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Protein.Group),
      !is.na(Precursor.Quantity),
      Precursor.Quantity > 0
    ) %>%
    dplyr::group_by(
      Run,
      Protein.Group
    ) %>%
    dplyr::summarise(
      Signal = sum(Precursor.Quantity),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Is_Contaminant = grepl(
        "^(cRAP-|CON_)",
        Protein.Group
      )
    ) %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Total_Signal = sum(Signal),
      Contaminant_Signal = sum(Signal[Is_Contaminant]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Contaminant_Pct = dplyr::if_else(
        Total_Signal > 0,
        100 * Contaminant_Signal / Total_Signal,
        0
      )
    ) %>%
    dplyr::select(
      Run,
      Contaminant_Pct
    )
  
  
  ## --------------------------------------------------------------------------
  ## Merge all metrics
  ## --------------------------------------------------------------------------
  
  figures_of_merit_df <- run_reference %>%
    dplyr::left_join(
      identification_df,
      by = "Run"
    ) %>%
    dplyr::left_join(
      chromatography_df,
      by = "Run"
    ) %>%
    dplyr::left_join(
      rt_df,
      by = "Run"
    ) %>%
    dplyr::left_join(
      mass_df,
      by = "Run"
    ) %>%
    dplyr::left_join(
      peptide_df,
      by = "Run"
    ) %>%
    dplyr::left_join(
      charge_df,
      by = "Run"
    ) %>%
    dplyr::left_join(
      contaminant_df,
      by = "Run"
    ) %>%
    dplyr::arrange(Run)
  
  
  ## --------------------------------------------------------------------------
  ## CSV export
  ## --------------------------------------------------------------------------
  
  figures_of_merit_csv <- file.path(
    results_dir,
    paste0(
      "figures_of_merit_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".csv"
    )
  )
  
  write.csv(
    figures_of_merit_df,
    figures_of_merit_csv,
    row.names = FALSE
  )
  
  message(
    "Figures of Merit CSV saved to: ",
    figures_of_merit_csv
  )
  
} else {
  
  message(
    "[FOM] Section skipped due to missing inputs"
  )
}


## 3. Clean exit ---------------------------------------------------------------

invisible(NULL)