###############################################################################
# File: qc_intensity.R
#
# Purpose:
#   Evaluate signal intensity distributions at multiple levels, including raw
#   precursor intensities and protein-level LFQ-normalized intensities, to
#   assess quantitative consistency and normalization performance across runs.
#
# Inputs:
#   - report
#   - file_palette
#   - AR_DENSITY
#   - AR_BOXPLOT
#
# Outputs:
#   - prec_density
#   - intensity_boxplot
#   - p_density
#   - p_box
#
# Notes:
#   - Raw precursor distributions use positive Precursor.Quantity values only.
#   - Protein-level distributions use one positive PG.MaxLFQ value per
#     Run + Protein.Group.
#   - Repeated identical PG.MaxLFQ values are collapsed.
#   - Multiple distinct positive PG.MaxLFQ values for the same Run +
#     Protein.Group are treated as an input inconsistency and stop execution.
#   - A plot is skipped with an informative message when no valid positive
#     values are available for that subsection.
###############################################################################


## 1. Initialize output objects ------------------------------------------------

prec_density <- NULL
intensity_boxplot <- NULL
p_density <- NULL
p_box <- NULL


## 2. Raw precursor intensity distribution ------------------------------------

prec_density <- report %>%
  dplyr::select(
    Run,
    Precursor.Quantity
  ) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(Precursor.Quantity),
    is.finite(Precursor.Quantity),
    Precursor.Quantity > 0
  )


if (nrow(prec_density) > 0L) {

  legend_columns <- max(
    1L,
    ceiling(
      dplyr::n_distinct(
        prec_density$Run
      ) /
        20
    )
  )


  p_density <- ggplot2::ggplot(
    prec_density,
    ggplot2::aes(
      x = Precursor.Quantity,
      fill = Run
    )
  ) +
    ggplot2::geom_density(
      alpha = 0.25,
      color = NA
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_manual(
      values = file_palette
    ) +
    ggplot2::labs(
      title = "Precursor Intensity Density (Raw Signal)",
      x = "Precursor quantity (log10)",
      y = "Density"
    ) +
    ggplot2::theme_classic() +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        ncol = legend_columns,
        byrow = FALSE
      )
    ) +
    ggplot2::theme(
      aspect.ratio = AR_DENSITY,
      legend.position = "right",
      legend.text = ggplot2::element_text(
        size = 9
      ),
      legend.key.size = grid::unit(
        0.5,
        "cm"
      ),
      legend.title = ggplot2::element_blank()
    )


  print(
    p_density
  )

} else {

  message(
    paste0(
      "[INTENSITY] Skipping raw precursor intensity density: ",
      "no valid positive Precursor.Quantity values available"
    )
  )
}


## 3. Validate protein-level LFQ values ----------------------------------------

pg_maxlfq_values <- report %>%
  dplyr::select(
    Run,
    Protein.Group,
    PG.MaxLFQ
  ) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(Protein.Group),
    nzchar(
      trimws(
        Protein.Group
      )
    ),
    !is.na(PG.MaxLFQ),
    is.finite(PG.MaxLFQ),
    PG.MaxLFQ > 0
  ) %>%
  dplyr::distinct(
    Run,
    Protein.Group,
    PG.MaxLFQ
  )


pg_maxlfq_consistency <- pg_maxlfq_values %>%
  dplyr::group_by(
    Run,
    Protein.Group
  ) %>%
  dplyr::summarise(
    Distinct_PG_MaxLFQ_Values = dplyr::n_distinct(
      PG.MaxLFQ
    ),
    PG_MaxLFQ_Values = paste(
      sort(
        unique(
          PG.MaxLFQ
        )
      ),
      collapse = ";"
    ),
    .groups = "drop"
  )


inconsistent_pg_maxlfq <- pg_maxlfq_consistency %>%
  dplyr::filter(
    Distinct_PG_MaxLFQ_Values > 1L
  )


if (nrow(inconsistent_pg_maxlfq) > 0L) {

  example_inconsistencies <- inconsistent_pg_maxlfq %>%
    dplyr::slice_head(
      n = 5L
    ) %>%
    dplyr::transmute(
      Description = paste0(
        "Run=",
        as.character(Run),
        ", Protein.Group=",
        Protein.Group,
        ", PG.MaxLFQ=",
        PG_MaxLFQ_Values
      )
    ) %>%
    dplyr::pull(
      Description
    )


  stop(
    paste0(
      "[INTENSITY] Multiple distinct positive PG.MaxLFQ values were found ",
      "for ",
      nrow(inconsistent_pg_maxlfq),
      " Run + Protein.Group pair(s). A protein-level intensity cannot be ",
      "selected without an explicit aggregation rule. Examples: ",
      paste(
        example_inconsistencies,
        collapse = " | "
      )
    ),
    call. = FALSE
  )
}


## 4. Protein-level LFQ intensity distribution ---------------------------------

intensity_boxplot <- pg_maxlfq_values %>%
  dplyr::select(
    Run,
    Protein.Group,
    PG.MaxLFQ
  )


if (nrow(intensity_boxplot) > 0L) {

  duplicated_protein_values <- intensity_boxplot %>%
    dplyr::count(
      Run,
      Protein.Group,
      name = "Rows"
    ) %>%
    dplyr::filter(
      Rows != 1L
    )


  if (nrow(duplicated_protein_values) > 0L) {
    stop(
      paste0(
        "[INTENSITY] Internal consistency failure: ",
        nrow(duplicated_protein_values),
        " Run + Protein.Group pair(s) remain duplicated after PG.MaxLFQ ",
        "deduplication."
      ),
      call. = FALSE
    )
  }


  p_box <- ggplot2::ggplot(
    intensity_boxplot,
    ggplot2::aes(
      x = Run,
      y = log10(PG.MaxLFQ),
      fill = Run
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.alpha = 0.25,
      width = 0.7,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = file_palette
    ) +
    ggplot2::labs(
      title = "Protein Intensity Distribution (LFQ-normalized)",
      x = "Run",
      y = "PG.MaxLFQ (log10)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      aspect.ratio = AR_BOXPLOT,
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none"
    )


  print(
    p_box
  )

} else {

  message(
    paste0(
      "[INTENSITY] Skipping protein-level LFQ intensity boxplot: ",
      "no valid positive PG.MaxLFQ values available"
    )
  )
}


## 5. Clean exit ---------------------------------------------------------------

invisible(NULL)
