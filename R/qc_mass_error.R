###############################################################################
# File: qc_mass_error.R
#
# Purpose:
#   Compute and visualize precursor mass measurement error (ppm) across runs
#   in order to assess mass accuracy, centering, and stability of the DIA
#   acquisition.
#
# Notes:
#   - Mass error is consistently expressed in ppm across all plots.
###############################################################################

## 1. Compute precursor mass error (ppm) ---------------------------------------

mass_error <- report %>% 
  dplyr::select(
    Run,
    Precursor.Mz,
    Precursor.Quantity,
    Ms1.Apex.Mz.Delta
  ) %>% 
  dplyr::filter(
    !is.na(Precursor.Quantity),
    Precursor.Quantity > 0,
    !is.na(Precursor.Mz),
    Precursor.Mz > 0
  ) %>% 
  dplyr::mutate(
    Mass_Error_PPM = (Ms1.Apex.Mz.Delta / Precursor.Mz) * 1e6
  ) %>% 
  dplyr::filter(
    !is.na(Mass_Error_PPM)
  )


## 2. Summary statistics per run -----------------------------------------------

mass_error_stat <- mass_error %>% 
  dplyr::group_by(Run) %>% 
  dplyr::summarise(
    median_val = median(Mass_Error_PPM, na.rm = TRUE),
    .groups = "drop"
  )


## 3. Visualization: mass error distributions ----------------------------------

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
      fill     = "steelblue",
      color    = "white"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      color      = "red",
      linewidth  = 1
    ) +
    ggplot2::geom_vline(
      data = mass_error_stat_subset,
      ggplot2::aes(xintercept = median_val),
      linetype   = "dashed",
      color      = "darkred",
      linewidth  = 1
    ) +
    ggplot2::geom_text(
      data = mass_error_stat_subset,
      ggplot2::aes(
        x     = median_val,
        y     = Inf,
        label = paste0(
          "Median = ",
          round(median_val, 2),
          " ppm"
        )
      ),
      vjust       = 2,
      color       = "darkred",
      fontface    = "bold",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(
      ~ Run,
      scales = "free_y",
      ncol   = 3
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-10, 10)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(-10, 10, 2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = "Mass Error Distribution per Run",
      subtitle = "Deviation from theoretical precursor m/z (target = 0 ppm)",
      x        = "Mass error (ppm)",
      y        = "Count"
    )
  
  print(p_mass)
}



## 4. Visualization: mass error vs precursor m/z (ppm-consistent) -------------

binned <- mass_error %>%
  dplyr::mutate(
    Mass_Error_PPM = (Ms1.Apex.Mz.Delta / Precursor.Mz) * 1e6,
    bin = cut(Precursor.Mz, breaks = 50)
  ) %>%
  dplyr::group_by(Run, bin) %>%
  dplyr::summarise(
    mz = mean(Precursor.Mz),
    median_ppm = median(Mass_Error_PPM),
    .groups = "drop"
  )

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
      linetype   = "dashed",
      color      = "red"
    ) +
    ggplot2::facet_wrap(
      ~ Run,
      scales = "fixed",
      ncol   = 3
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Mass Error Bias Across m/z (Binned Median)",
      x = "Precursor m/z",
      y = "Median mass error (ppm)"
    )
  
  print(p_bin)
}