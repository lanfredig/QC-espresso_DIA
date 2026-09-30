###############################################################################
# File: qc_chromatography.R
#
# Purpose:
#   Generate chromatography-focused QC plots across runs, covering retention-
#   time identification density, precursor peak width at half maximum (FWHM),
#   and precursor peak base width distributions.
#
# Inputs (expected to exist before sourcing):
#   - report:
#       Data frame containing at least the following columns:
#       Run, RT, FWHM, RT.Start, RT.Stop
#   - file_palette:
#       Named color vector used to keep run colors consistent where applicable
#
# Outputs:
#   - Printed ggplot objects for:
#       1. Identification density over retention time per run
#       2. FWHM distribution per run, with median annotation
#       3. Peak base width distribution per run, with median annotation
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - ggplot2
#   - magrittr
#
# Notes:
#   - Runs are split into groups of three so that faceted plots remain readable
#     within the fixed figure size used by the QMD report
#   - FWHM and peak base width are converted to seconds before plotting
#   - Median values are computed per run and displayed as dashed red reference
#     lines with text annotations
#   - This script assumes all required objects are already available in the
#     environment before sourcing
###############################################################################


## 1- Identification density over retention time

id_vs_rt <- report %>% 
  dplyr::select(Run, RT) %>% 
  dplyr::mutate(RT = round(RT, 1))

runs <- unique(id_vs_rt$Run)

# split runs into groups of 3
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
      ncol   = 3
    ) +
    ggplot2::scale_fill_manual(values = file_palette) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      title    = "Identifications over Retention Time per Run",
      subtitle = "Distribution of peptide identifications across the LC gradient",
      x        = "Retention Time (min)",
      y        = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none"
    )
  
  print(p_rt)
}


## 2. Peak width at half maximum (FWHM)

FWHM <- report %>% 
  dplyr::mutate(FWHM = 60 * round(FWHM, 3))

FWHM_stat <- FWHM %>% 
  dplyr::group_by(Run) %>% 
  dplyr::summarise(
    median_val = median(FWHM, na.rm = TRUE),
    .groups = "drop"
  )

runs <- unique(FWHM$Run)

# split runs into groups of 3
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
      fill  = "grey70",
      color = "black"
    ) +
    ggplot2::facet_wrap(
      ~ Run,
      scales = "free_y",
      ncol   = 3
    ) +
    ggplot2::geom_vline(
      data = stat_subset,
      ggplot2::aes(xintercept = median_val),
      color    = "red",
      linetype = "dashed"
    ) +
    ggplot2::geom_text(
      data = stat_subset,
      ggplot2::aes(
        x     = median_val,
        y     = Inf,
        label = paste0("Median = ", round(median_val, 2))
      ),
      nudge_x     = 0.5,
      vjust       = 1.5,
      hjust       = 0,
      color       = "red",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      title    = "FWHM Distribution per Run",
      subtitle = "Peak width at half maximum",
      x        = "FWHM (s)",
      y        = "Count"
    ) +
    ggplot2::theme_minimal()
  
  print(p_fwhm)
}

## 3. Peak base width

peak_base <- report %>% 
  dplyr::select(Run, RT.Start, RT.Stop) %>% 
  dplyr::mutate(
    Peakbase = (RT.Stop - RT.Start) * 60
  )

peak_base_stat <- peak_base %>% 
  dplyr::group_by(Run) %>% 
  dplyr::summarise(
    median_val = median(Peakbase, na.rm = TRUE),
    .groups = "drop"
  )

runs <- unique(peak_base$Run)

# split runs into groups of 3
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
      fill     = "grey70",
      color    = "black"
    ) +
    ggplot2::facet_wrap(
      ~ Run,
      scales = "free_y",
      ncol   = 3
    ) +
    ggplot2::geom_vline(
      data = stat_subset,
      ggplot2::aes(xintercept = median_val),
      color    = "red",
      linetype = "dashed"
    ) +
    ggplot2::geom_text(
      data = stat_subset,
      ggplot2::aes(
        x     = median_val,
        y     = Inf,
        label = paste0("Median = ", round(median_val, 2))
      ),
      nudge_x     = 0.5,
      vjust       = 1.5,
      hjust       = 0,
      color       = "red",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      title    = "Peak Base Width per Run",
      subtitle = "Distribution of peak base width",
      x        = "Base Width (s)",
      y        = "Count"
    ) +
    ggplot2::theme_minimal()
  
  print(p_base)
}
