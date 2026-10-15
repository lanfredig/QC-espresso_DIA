###############################################################################
# File: qc_diann_quality_metrics.R
#
# Purpose:
#   Visualize DIA-NN-derived quantification quality fields at the precursor
#   and protein-group levels, including the joint distribution of
#   Quantity.Quality and Empirical.Quality and the relationship between
#   PG.MaxLFQ.Quality and the number of distinct supporting precursors.
#
# Inputs (expected to exist before sourcing):
#   - report       : Filtered DIA-NN report data frame
#   - file_palette : Run-level color palette created by io_load.R
#
# Outputs (objects defined):
#   - qc_precursor : Valid Quantity.Quality and Empirical.Quality observations
#   - p_joint      : Joint precursor quality-metric plot for the current run group
#   - qc_protein   : Protein-group quality values and distinct precursor counts
#   - p_pg         : Protein-group quality versus precursor-count plot for the
#                    current run group
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - ggplot2
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - Quantity.Quality, Empirical.Quality, and PG.MaxLFQ.Quality are
#     DIA-NN-derived QuantUMS-related quality fields.
#   - The precursor-level plot shows the joint distribution of
#     Quantity.Quality and Empirical.Quality using two-dimensional bins.
#   - The protein-group-level plot relates PG.MaxLFQ.Quality to the number
#     of distinct Precursor.Id values per run and protein group.
#   - The two subsections are evaluated independently.
#   - Runs are displayed in groups of up to three facets.
#   - A subsection is skipped when its required fields contain no valid
#     observations.
###############################################################################

## 1. Joint precursor quality-metric space ------------------------------------
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
        title = "Joint Precursor Quality-Metric Space",
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

## 2. Protein-group MaxLFQ quality versus precursor count ----------------------
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
        title = "Protein-Group MaxLFQ Quality versus Precursor Count",
        subtitle = "PG.MaxLFQ.Quality versus  number of distinct precursors per protein group",
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
