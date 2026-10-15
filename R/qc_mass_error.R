###############################################################################
# File: qc_mass_error.R
#
# Purpose:
#   Calculate precursor mass error in parts per million from the DIA-NN
#   Ms1.Apex.Mz.Delta and Precursor.Mz fields, summarize its median per run,
#   and visualize its distribution and variation across precursor m/z.
#
# Inputs (expected to exist before sourcing):
#   - report : Filtered DIA-NN report data frame
#
# Outputs (objects defined):
#   - mass_error             : Valid positive-quantity observations with
#                              calculated precursor mass error in ppm
#   - mass_error_stat        : Median precursor mass error per run
#   - p_mass_runs            : Runs included in the mass-error distributions
#   - p_mass_run_groups      : Run groups used for facetted distribution plots
#   - p_mass                 : Mass-error distribution plot for the current
#                              run group
#   - binned                 : Median mass error across precursor m/z bins
#   - p_bin                  : Binned mass-error plot for the current run group
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - ggplot2
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - Mass error is calculated as Ms1.Apex.Mz.Delta divided by Precursor.Mz
#     and converted to parts per million.
#   - Only observations with positive Precursor.Quantity and positive
#     Precursor.Mz values are included.
#   - Mass error is expressed consistently in parts per million across all
#     plots.
#   - The mass-error distribution is displayed within a fixed visual range
#     from -10 to 10 ppm.
#   - The relationship between mass error and precursor m/z is summarized
#     using 50 precursor m/z bins and the median mass error within each bin.
#   - Runs are displayed in groups of up to three facets.
#   - Plot generation is skipped when no valid mass-error observations are
#     available.
###############################################################################

## 1. Compute precursor mass error --------------------------------------------

mass_error <- report %>%
  dplyr::select(
    Run,
    Precursor.Mz,
    Precursor.Quantity,
    Ms1.Apex.Mz.Delta
  ) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(Precursor.Quantity),
    Precursor.Quantity > 0,
    !is.na(Precursor.Mz),
    Precursor.Mz > 0,
    !is.na(Ms1.Apex.Mz.Delta),
    is.finite(Precursor.Mz),
    is.finite(Ms1.Apex.Mz.Delta)
  ) %>%
  dplyr::mutate(
    Mass_Error_PPM = (Ms1.Apex.Mz.Delta / Precursor.Mz) * 1e6
  ) %>%
  dplyr::filter(
    !is.na(Mass_Error_PPM),
    is.finite(Mass_Error_PPM)
  )

if (nrow(mass_error) == 0L) {
  message(
    "[MASS ERROR] Skipping: no valid positive-quantity mass-error ",
    "observations available"
  )
} else {
  ## 2. Summary statistics per run --------------------------------------------

  mass_error_stat <- mass_error %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
      median_val = median(Mass_Error_PPM, na.rm = TRUE),
      .groups = "drop"
    )

  ## 3. Mass-error distributions ----------------------------------------------

  p_mass_runs <- unique(mass_error$Run)
  p_mass_run_groups <- split(
    p_mass_runs,
    ceiling(seq_along(p_mass_runs) / 3)
  )

  for (grp in p_mass_run_groups) {
    mass_error_subset <- mass_error %>%
      dplyr::filter(Run %in% grp)

    mass_error_stat_subset <- mass_error_stat %>%
      dplyr::filter(Run %in% grp)

    p_mass <- ggplot2::ggplot(
      mass_error_subset,
      ggplot2::aes(x = Mass_Error_PPM)
    ) +
      ggplot2::geom_histogram(
        binwidth = 0.2,
        fill = "steelblue",
        color = "white"
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "red",
        linewidth = 1
      ) +
      ggplot2::geom_vline(
        data = mass_error_stat_subset,
        ggplot2::aes(xintercept = median_val),
        linetype = "dashed",
        color = "darkred",
        linewidth = 1
      ) +
      ggplot2::geom_text(
        data = mass_error_stat_subset,
        ggplot2::aes(
          x = median_val,
          y = Inf,
          label = paste0(
            "Median = ",
            round(median_val, 2),
            " ppm"
          )
        ),
        vjust = 2,
        color = "darkred",
        fontface = "bold",
        inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(
        ~ Run,
        scales = "free_y",
        ncol = 3
      ) +
      ggplot2::coord_cartesian(
        xlim = c(-10, 10)
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq(-10, 10, 2)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Mass Error Distribution per Run",
        subtitle = "Deviation from theoretical precursor m/z (target = 0 ppm)",
        x = "Mass error (ppm)",
        y = "Count"
      )

    print(p_mass)
  }

## 4. Binned median mass error across precursor m/z -------------------------
  binned <- mass_error %>%
    dplyr::mutate(
      bin = cut(Precursor.Mz, breaks = 50)
    ) %>%
    dplyr::group_by(Run, bin) %>%
    dplyr::summarise(
      mz = mean(Precursor.Mz, na.rm = TRUE),
      median_ppm = median(Mass_Error_PPM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      !is.na(mz),
      is.finite(mz),
      !is.na(median_ppm),
      is.finite(median_ppm)
    )

  if (nrow(binned) > 0L) {
    runs <- unique(binned$Run)
    run_groups <- split(
      runs,
      ceiling(seq_along(runs) / 3)
    )

    for (grp in run_groups) {
      df_subset <- binned %>%
        dplyr::filter(Run %in% grp)

      p_bin <- ggplot2::ggplot(
        df_subset,
        ggplot2::aes(mz, median_ppm)
      ) +
        ggplot2::geom_line(color = "blue", linewidth = 1) +
        ggplot2::geom_point(size = 1) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = "red"
        ) +
        ggplot2::facet_wrap(
          ~ Run,
          scales = "fixed",
          ncol = 3
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Median Mass Error Across Precursor m/z Bins",
          x = "Precursor m/z",
          y = "Median mass error (ppm)"
        )

      print(p_bin)
    }
  } else {
    message(
      "[MASS ERROR] Skipping binned m/z trend: ",
      "no valid binned mass-error observations available"
    )
  }
}

invisible(NULL)
