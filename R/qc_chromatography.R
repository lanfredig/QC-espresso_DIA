###############################################################################
# File: qc_chromatography.R
#
# Purpose:
#   Generate chromatography-focused QC plots across runs, covering retention-
#   time identification density, precursor peak width at half maximum (FWHM),
#   and precursor peak base width distributions.
#
# Inputs:
#   - report
#   - file_palette
#
# Outputs:
#   - Printed ggplot objects when valid values are available
#   - Informative messages when a subsection has no valid values
###############################################################################

## 1. Identification density over retention time ------------------------------

id_vs_rt <- report %>%
  dplyr::select(Run, RT) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(RT),
    is.finite(RT)
  ) %>%
  dplyr::mutate(RT = round(RT, 1))

if (nrow(id_vs_rt) > 0L) {
  runs <- unique(id_vs_rt$Run)
  run_groups <- split(runs, ceiling(seq_along(runs) / 3))

  for (grp in run_groups) {
    df_subset <- id_vs_rt %>%
      dplyr::filter(Run %in% grp)

    p_rt <- ggplot2::ggplot(
      df_subset,
      ggplot2::aes(RT, fill = Run)
    ) +
      ggplot2::geom_histogram(
        binwidth = 1,
        color = "black"
      ) +
      ggplot2::facet_wrap(
        ~ Run,
        scales = "free_y",
        ncol = 3
      ) +
      ggplot2::scale_fill_manual(values = file_palette) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::labs(
        title = "Identifications over Retention Time per Run",
        subtitle = "Distribution of peptide identifications across the LC gradient",
        x = "Retention Time (min)",
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

    print(p_rt)
  }
} else {
  message(
    "[CHROMATOGRAPHY] Skipping retention-time identification density: ",
    "no valid RT values available"
  )
}

## 2. Peak width at half maximum ----------------------------------------------

FWHM <- report %>%
  dplyr::select(Run, FWHM) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(FWHM),
    is.finite(FWHM)
  ) %>%
  dplyr::mutate(FWHM = 60 * round(FWHM, 3))

if (nrow(FWHM) > 0L) {
  FWHM_stat <- FWHM %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      median_val = median(FWHM, na.rm = TRUE),
      .groups = "drop"
    )

  runs <- unique(FWHM$Run)
  run_groups <- split(runs, ceiling(seq_along(runs) / 3))

  for (grp in run_groups) {
    df_subset <- FWHM %>%
      dplyr::filter(Run %in% grp)

    stat_subset <- FWHM_stat %>%
      dplyr::filter(Run %in% grp)

    p_fwhm <- ggplot2::ggplot(
      df_subset,
      ggplot2::aes(FWHM)
    ) +
      ggplot2::geom_histogram(
        binwidth = 0.5,
        fill = "grey70",
        color = "black"
      ) +
      ggplot2::facet_wrap(
        ~ Run,
        scales = "free_y",
        ncol = 3
      ) +
      ggplot2::geom_vline(
        data = stat_subset,
        ggplot2::aes(xintercept = median_val),
        color = "red",
        linetype = "dashed"
      ) +
      ggplot2::geom_text(
        data = stat_subset,
        ggplot2::aes(
          x = median_val,
          y = Inf,
          label = paste0("Median = ", round(median_val, 2))
        ),
        nudge_x = 0.5,
        vjust = 1.5,
        hjust = 0,
        color = "red",
        inherit.aes = FALSE
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::labs(
        title = "FWHM Distribution per Run",
        subtitle = "Peak width at half maximum",
        x = "FWHM (s)",
        y = "Count"
      ) +
      ggplot2::theme_minimal()

    print(p_fwhm)
  }
} else {
  message(
    "[CHROMATOGRAPHY] Skipping FWHM distribution: ",
    "no valid FWHM values available"
  )
}

## 3. Peak base width ----------------------------------------------------------

peak_base <- report %>%
  dplyr::select(Run, RT.Start, RT.Stop) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(RT.Start),
    !is.na(RT.Stop),
    is.finite(RT.Start),
    is.finite(RT.Stop)
  ) %>%
  dplyr::mutate(
    Peakbase = (RT.Stop - RT.Start) * 60
  ) %>%
  dplyr::filter(
    is.finite(Peakbase)
  )

if (nrow(peak_base) > 0L) {
  peak_base_stat <- peak_base %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      median_val = median(Peakbase, na.rm = TRUE),
      .groups = "drop"
    )

  runs <- unique(peak_base$Run)
  run_groups <- split(runs, ceiling(seq_along(runs) / 3))

  for (grp in run_groups) {
    df_subset <- peak_base %>%
      dplyr::filter(Run %in% grp)

    stat_subset <- peak_base_stat %>%
      dplyr::filter(Run %in% grp)

    p_base <- ggplot2::ggplot(
      df_subset,
      ggplot2::aes(Peakbase)
    ) +
      ggplot2::geom_histogram(
        binwidth = 0.1,
        fill = "grey70",
        color = "black"
      ) +
      ggplot2::facet_wrap(
        ~ Run,
        scales = "free_y",
        ncol = 3
      ) +
      ggplot2::geom_vline(
        data = stat_subset,
        ggplot2::aes(xintercept = median_val),
        color = "red",
        linetype = "dashed"
      ) +
      ggplot2::geom_text(
        data = stat_subset,
        ggplot2::aes(
          x = median_val,
          y = Inf,
          label = paste0("Median = ", round(median_val, 2))
        ),
        nudge_x = 0.5,
        vjust = 1.5,
        hjust = 0,
        color = "red",
        inherit.aes = FALSE
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::labs(
        title = "Peak Base Width per Run",
        subtitle = "Distribution of peak base width",
        x = "Base Width (s)",
        y = "Count"
      ) +
      ggplot2::theme_minimal()

    print(p_base)
  }
} else {
  message(
    "[CHROMATOGRAPHY] Skipping peak base width distribution: ",
    "no valid RT.Start and RT.Stop pairs available"
  )
}

invisible(NULL)
