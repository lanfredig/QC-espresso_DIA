
###############################################################################
# File: qc_diann_quality_metrics.R
#
# Purpose:
#   Visualization-only QC for DIA-NN quantification quality
#
# Plots included:
#   1. Joint precursor reliability space
#   2. Protein-group MaxLFQ quality vs precursor support
#
# Inputs (expected to exist before sourcing):
#   - report          : Filtered DIA-NN report data frame
#   - PANELS_PER_ROW  : Number of panels per row
#   - file_palette    : Named color vector for Run
#   - AR_DENSITY      : Aspect ratio for density-like plots
#
# Required packages (loaded by QMD):
#   - dplyr
#   - ggplot2
###############################################################################



## 1. Joint precursor reliability space ---------------------------------------


qc_precursor <- report %>% 
  dplyr::select(
    Run,
    Quantity.Quality,
    Empirical.Quality
  ) %>% 
  dplyr::filter(
    !is.na(Quantity.Quality),
    !is.na(Empirical.Quality)
  )

runs <- unique(qc_precursor$Run)

# split runs into groups of 3
run_groups <- split(runs, ceiling(seq_along(runs) / 3))

for (grp in run_groups) {
  
  df_subset <- qc_precursor %>% 
    dplyr::filter(Run %in% grp)
  
  p_joint <- ggplot2::ggplot(
    df_subset,
    ggplot2::aes(
      x = Quantity.Quality,
      y = Empirical.Quality
    )
  ) +
    ggplot2::geom_bin2d(bins = 80) +
    ggplot2::scale_fill_gradientn(
      colours = c("navy", "turquoise4", "yellow3"),
      name    = "Precursor count",
      trans   = "sqrt"
    ) +
    ggplot2::geom_vline(
      xintercept = 0.5,
      linetype   = "dashed",
      color      = "red"
    ) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      linetype   = "dashed",
      color      = "red"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    ggplot2::facet_wrap(
      ~ Run,
      ncol = 3
    ) +
    ggplot2::labs(
      title    = "Joint Precursor Reliability Space",
      subtitle = "Quantity.Quality vs Empirical.Quality",
      x        = "Quantity.Quality",
      y        = "Empirical.Quality"
    ) +
    ggplot2::theme_minimal()
  
  print(p_joint)
}



## 2. Protein-group MaxLFQ quality vs precursor support -----------------------


qc_protein <- report %>% 
  dplyr::select(
    Run,
    Protein.Group,
    Precursor.Id,
    PG.MaxLFQ.Quality
  ) %>% 
  dplyr::filter(!is.na(PG.MaxLFQ.Quality)) %>% 
  dplyr::group_by(
    Run,
    Protein.Group,
    PG.MaxLFQ.Quality
  ) %>% 
  dplyr::summarise(
    n_precursors = dplyr::n_distinct(Precursor.Id),
    .groups = "drop"
  )

runs <- unique(qc_protein$Run)

# split runs into groups of 3
run_groups <- split(runs, ceiling(seq_along(runs) / 3))

for (grp in run_groups) {
  
  df_subset <- qc_protein %>% 
    dplyr::filter(Run %in% grp)
  
  p_pg <- ggplot2::ggplot(
    df_subset,
    ggplot2::aes(
      x     = n_precursors,
      y     = PG.MaxLFQ.Quality,
      color = Run
    )
  ) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::geom_hline(
      yintercept = 0.7,
      linetype   = "dashed",
      color      = "blue"
    ) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      linetype   = "dashed",
      color      = "red"
    ) +
    ggplot2::facet_wrap(
      ~ Run,
      ncol = 3
    ) +
    ggplot2::scale_color_manual(
      values = file_palette,
      drop   = FALSE
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title    = "Protein Quantification Stability",
      subtitle = "MaxLFQ consistency versus precursor support",
      x        = "Number of precursors (log scale)",
      y        = "PG.MaxLFQ.Quality"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  print(p_pg)
}
