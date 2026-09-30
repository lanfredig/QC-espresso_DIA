###############################################################################
# File: qc_diann_quality_metrics.R
#
# Purpose:
#   Visualization-only QC for DIA-NN quantification quality.
#
# Plots:
#   1. Joint precursor reliability space
#   2. Protein-group MaxLFQ quality versus precursor support
#
# Notes:
#   - The two subsections are evaluated independently.
#   - A subsection is skipped with an informative message when its columns are
#     present but contain no valid observations.
###############################################################################

## 1. Joint precursor reliability space ---------------------------------------

qc_precursor <- report %>%
  dplyr::select(
    Run,
    Quantity.Quality,
    Empirical.Quality
  ) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(Quantity.Quality),
    !is.na(Empirical.Quality),
    is.finite(Quantity.Quality),
    is.finite(Empirical.Quality)
  )

if (nrow(qc_precursor) > 0L) {
  runs <- unique(qc_precursor$Run)
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
        name = "Precursor count",
        trans = "sqrt"
      ) +
      ggplot2::geom_vline(
        xintercept = 0.5,
        linetype = "dashed",
        color = "red"
      ) +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        color = "red"
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
        title = "Joint Precursor Reliability Space",
        subtitle = "Quantity.Quality vs Empirical.Quality",
        x = "Quantity.Quality",
        y = "Empirical.Quality"
      ) +
      ggplot2::theme_minimal()

    print(p_joint)
  }
} else {
  message(
    "[DIA-NN QUALITY] Skipping joint precursor reliability space: ",
    "no valid Quantity.Quality and Empirical.Quality pairs available"
  )
}

## 2. Protein-group MaxLFQ quality versus precursor support --------------------

qc_protein <- report %>%
  dplyr::select(
    Run,
    Protein.Group,
    Precursor.Id,
    PG.MaxLFQ.Quality
  ) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(Protein.Group),
    !is.na(Precursor.Id),
    !is.na(PG.MaxLFQ.Quality),
    is.finite(PG.MaxLFQ.Quality)
  ) %>%
  dplyr::group_by(
    Run,
    Protein.Group,
    PG.MaxLFQ.Quality
  ) %>%
  dplyr::summarise(
    n_precursors = dplyr::n_distinct(Precursor.Id),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_precursors > 0
  )

if (nrow(qc_protein) > 0L) {
  runs <- unique(qc_protein$Run)
  run_groups <- split(runs, ceiling(seq_along(runs) / 3))

  for (grp in run_groups) {
    df_subset <- qc_protein %>%
      dplyr::filter(Run %in% grp)

    p_pg <- ggplot2::ggplot(
      df_subset,
      ggplot2::aes(
        x = n_precursors,
        y = PG.MaxLFQ.Quality,
        color = Run
      )
    ) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::geom_hline(
        yintercept = 0.7,
        linetype = "dashed",
        color = "blue"
      ) +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        color = "red"
      ) +
      ggplot2::facet_wrap(
        ~ Run,
        ncol = 3
      ) +
      ggplot2::scale_color_manual(
        values = file_palette,
        drop = FALSE
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = "Protein Quantification Stability",
        subtitle = "MaxLFQ consistency versus precursor support",
        x = "Number of precursors (log scale)",
        y = "PG.MaxLFQ.Quality"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    print(p_pg)
  }
} else {
  message(
    "[DIA-NN QUALITY] Skipping protein quantification stability: ",
    "no valid PG.MaxLFQ.Quality observations available"
  )
}

invisible(NULL)
