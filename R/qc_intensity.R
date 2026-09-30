###############################################################################
# File: qc_intensity.R
#
# Purpose:
#   Evaluate signal intensity distributions at multiple levels, including
#   raw precursor intensities and protein-level LFQ-normalized intensities,
#   to assess quantitative consistency and normalization performance across
#   runs.
#
# Inputs (expected to exist before sourcing):
#   - report       : Filtered DIA-NN report data frame
#   - file_palette : Run-level color palette
#   - AR_DENSITY   : Aspect ratio constant for density plots
#   - AR_BOXPLOT   : Aspect ratio constant for boxplots
#
# Outputs:
#   - Density plots of raw precursor intensities
#   - Boxplots of protein-level LFQ intensity distributions
###############################################################################


## 2. Raw precursor intensity distribution ------------------------------------

prec_density <- report %>%
  dplyr::select(Run, Precursor.Quantity) %>%
  dplyr::filter(
    !is.na(Precursor.Quantity),
    Precursor.Quantity > 0
  )

legend_columns <- max(
  1,
  ceiling(dplyr::n_distinct(prec_density$Run) / 20)
)

p_density <- ggplot2::ggplot(
  prec_density,
  ggplot2::aes(
    x    = Precursor.Quantity,
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
    x     = "Precursor quantity (log10)",
    y     = "Density"
  ) +
  ggplot2::theme_classic() +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      ncol  = legend_columns,
      byrow = FALSE
    )
  ) +
  ggplot2::theme(
    aspect.ratio     = AR_DENSITY,
    legend.position  = "right",
    legend.text      = ggplot2::element_text(size = 9),
    legend.key.size  = grid::unit(0.5, "cm"),
    legend.title     = ggplot2::element_blank()
  )

print(p_density)


## 3. Protein-level LFQ intensity distribution --------------------------------
## IMPORTANT: ensure one value per Run + Protein.Group

intensity_boxplot <- report %>% 
  dplyr::select(Run, Protein.Group, PG.MaxLFQ) %>% 
  dplyr::filter(
    !is.na(PG.MaxLFQ),
    PG.MaxLFQ > 0
  ) %>% 
  dplyr::distinct(Run, Protein.Group, PG.MaxLFQ)

## OPTIONAL CONSISTENCY CHECK (can be commented if too verbose)
dup_check <- report %>% 
  dplyr::filter(!is.na(PG.MaxLFQ), PG.MaxLFQ > 0) %>% 
  dplyr::count(Run, Protein.Group) %>% 
  dplyr::filter(n > 1)

if (nrow(dup_check) > 0) {
  message("[INTENSITY] Collapsing duplicated PG.MaxLFQ per Run + Protein.Group")
}

p_box <- ggplot2::ggplot(
  intensity_boxplot,
  ggplot2::aes(
    x    = Run,
    y    = log10(PG.MaxLFQ),
    fill = Run
  )
) +
  ggplot2::geom_boxplot(
    outlier.alpha = 0.25,
    width         = 0.7,
    color         = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = file_palette
  ) +
  ggplot2::labs(
    title = "Protein Intensity Distribution (LFQ-normalized)",
    x     = "Run",
    y     = "PG.MaxLFQ (log10)"
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

print(p_box)
