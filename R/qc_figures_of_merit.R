###############################################################################
# File: qc_figures_of_merit.R
#
# Purpose:
#   Generate a consolidated run-level QC summary table containing
#   identification, chromatography, retention-time deviation, mass-error,
#   peptide, precursor-charge, and contaminant-signal metrics, export the
#   table as CSV, and expose it for rendering in the QMD report.
#
# Inputs (expected to exist before sourcing):
#   - report      : Filtered DIA-NN report data frame
#   - results_dir : Canonical output directory created by io_load.R
#
# Outputs (objects defined):
#   - count_missed_tryptic_cleavages : Helper function for counting missed
#                                      tryptic cleavages, if not already defined
#   - fom_available                  : Availability flag for the required inputs
#   - run_reference                 : Distinct runs included in the summary
#   - identification_df             : Identification counts per run
#   - chromatography_df             : Median FWHM and peak base width per run
#   - rt_reference                  : Median retention time per run and precursor
#   - n_runs_rt                     : Number of runs represented in rt_reference
#   - min_runs_required             : Minimum run-frequency threshold for
#                                     retention-time deviation assessment
#   - shared_precursors             : Precursors meeting the run-frequency
#                                     threshold
#   - rt_df                         : Median absolute retention-time deviation
#                                     per run
#   - mass_df                       : Median MS1 apex mass error per run
#   - peptide_df                    : Median peptide length and percentage of
#                                     peptides with missed cleavages per run
#   - charge_df                     : Percentages of charge 2 and charge 3
#                                     precursors per run
#   - contaminant_df                : Contaminant signal percentage per run
#   - figures_of_merit_df           : Consolidated run-level QC summary table
#   - figures_of_merit_csv          : Path to the exported Figures of Merit CSV
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - stringr
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - The missed-cleavage helper is defined only when a function with the same
#     name is not already available.
#   - Identification counts are based on distinct Precursor.Id,
#     Modified.Sequence, and Protein.Group values per run.
#   - Chromatographic widths are converted from minutes to seconds.
#   - Retention-time deviation is calculated for precursors observed in at
#     least half of the represented runs, with a minimum of two runs.
#   - Mass error is calculated from Ms1.Apex.Mz.Delta relative to
#     Precursor.Mz and expressed in parts per million.
#   - Missed-cleavage and peptide-length metrics are calculated from distinct
#     Stripped.Sequence values per run.
#   - Charge-state percentages are calculated from distinct combinations of
#     Run, Precursor.Id, and Precursor.Charge.
#   - Contaminants are recognized when CON_ or cRAP- occurs at the beginning
#     of any accession within a protein group.
###############################################################################

## 0. Shared missed-cleavage calculation --------------------------------------

if (!exists("count_missed_tryptic_cleavages", mode = "function")) {
  count_missed_tryptic_cleavages <- function(sequence) {
    cleavage_sites <- stringr::str_count(
      sequence,
      "(?<=[KR])(?!P)"
    )

    terminal_tryptic_site <- as.integer(
      stringr::str_detect(
        sequence,
        "[KR]$"
      )
    )

    pmax(
      0L,
      cleavage_sites - terminal_tryptic_site
    )
  }
}

## 1. Input validation ---------------------------------------------------------

fom_available <- TRUE
figures_of_merit_df <- NULL
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
  run_reference <- report %>%
    dplyr::distinct(Run) %>%
    dplyr::arrange(Run)

  identification_df <- report %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Precursors = dplyr::n_distinct(Precursor.Id),
      Peptides = dplyr::n_distinct(Modified.Sequence),
      Protein_Groups = dplyr::n_distinct(Protein.Group),
      .groups = "drop"
    )

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

  rt_reference <- report %>%
    dplyr::select(Run, Precursor.Id, RT) %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Precursor.Id),
      !is.na(RT)
    ) %>%
    dplyr::group_by(Run, Precursor.Id) %>%
    dplyr::summarise(
      RT_Run = median(RT, na.rm = TRUE),
      .groups = "drop"
    )

  n_runs_rt <- dplyr::n_distinct(rt_reference$Run)
  min_runs_required <- max(2, ceiling(0.5 * n_runs_rt))

  shared_precursors <- rt_reference %>%
    dplyr::group_by(Precursor.Id) %>%
    dplyr::summarise(
      n_runs = dplyr::n_distinct(Run),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_runs >= min_runs_required) %>%
    dplyr::pull(Precursor.Id)

  if (length(shared_precursors) > 0) {
    rt_df <- rt_reference %>%
      dplyr::filter(Precursor.Id %in% shared_precursors) %>%
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

  mass_df <- report %>%
    dplyr::filter(
      !is.na(Precursor.Quantity),
      Precursor.Quantity > 0,
      !is.na(Precursor.Mz),
      Precursor.Mz > 0,
      !is.na(Ms1.Apex.Mz.Delta)
    ) %>%
    dplyr::mutate(
      Mass_Error_PPM = (Ms1.Apex.Mz.Delta / Precursor.Mz) * 1e6
    ) %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Median_Mass_Error_ppm = median(
        Mass_Error_PPM,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  peptide_df <- report %>%
    dplyr::select(Run, Stripped.Sequence) %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Stripped.Sequence),
      nzchar(Stripped.Sequence)
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Peptide_Length = nchar(Stripped.Sequence),
      MissedCleavages = count_missed_tryptic_cleavages(
        Stripped.Sequence
      )
    ) %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Median_Peptide_Length = median(
        Peptide_Length,
        na.rm = TRUE
      ),
      Missed_Cleavage_Pct = 100 * mean(
        MissedCleavages > 0,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  charge_df <- report %>%
    dplyr::select(Run, Precursor.Id, Precursor.Charge) %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Precursor.Id),
      !is.na(Precursor.Charge)
    ) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      Charge2_Pct = 100 * mean(Precursor.Charge == 2, na.rm = TRUE),
      Charge3_Pct = 100 * mean(Precursor.Charge == 3, na.rm = TRUE),
      .groups = "drop"
    )

  contaminant_df <- report %>%
    dplyr::filter(
      !is.na(Run),
      !is.na(Protein.Group),
      !is.na(Precursor.Quantity),
      Precursor.Quantity > 0
    ) %>%
    dplyr::group_by(Run, Protein.Group) %>%
    dplyr::summarise(
      Signal = sum(Precursor.Quantity),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Is_Contaminant = grepl(
        "(^|;\\s*)(cRAP-|CON_)",
        Protein.Group,
        perl = TRUE
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
    dplyr::select(Run, Contaminant_Pct)

  figures_of_merit_df <- run_reference %>%
    dplyr::left_join(identification_df, by = "Run") %>%
    dplyr::left_join(chromatography_df, by = "Run") %>%
    dplyr::left_join(rt_df, by = "Run") %>%
    dplyr::left_join(mass_df, by = "Run") %>%
    dplyr::left_join(peptide_df, by = "Run") %>%
    dplyr::left_join(charge_df, by = "Run") %>%
    dplyr::left_join(contaminant_df, by = "Run") %>%
    dplyr::arrange(Run)

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
