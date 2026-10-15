###############################################################################
# File: globals.R
#
# Purpose:
#   Define global visualization settings used by the downstream QC modules,
#   including the shared color palette, panel layout configuration,
#   aspect-ratio presets, and standardized ggplot2 legend helpers.
#
# Inputs (expected to exist before sourcing):
#   - None
#
# Outputs (objects defined):
#   - cb_palette    : Global color palette for QC visualizations
#   - PANELS_PER_ROW: Number of plot panels displayed per row
#   - AR_DENSITY    : Aspect-ratio preset for density plots
#   - AR_BOXPLOT    : Aspect-ratio preset for boxplots
#   - AR_HIST       : Aspect-ratio preset for histograms
#   - legend_fill   : Helper function for standardized fill legends
#   - legend_color  : Helper function for standardized color legends
#
# Required R packages (loaded by QMD):
#   - ggplot2
#   - grid
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - This script has no data-dependent inputs.
#   - The defined objects are shared across downstream QC modules.
###############################################################################

## 1. Global color palette -----------------------------------------------------

cb_palette <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)


## 2. Panel layout configuration -----------------------------------------------

PANELS_PER_ROW <- 3


## 3. Aspect-ratio presets -----------------------------------------------------

AR_DENSITY <- 0.8
AR_BOXPLOT <- 0.35
AR_HIST    <- 0.7


## 4. Standardized ggplot2 legend helpers --------------------------------------

legend_fill <- function() {
  
  list(
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        ncol  = 1,
        byrow = TRUE
      )
    ),
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 9),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_blank()
    )
  )
}


legend_color <- function() {
  
  list(
    ggplot2::guides(
      color = ggplot2::guide_legend(
        ncol  = 1,
        byrow = TRUE
      )
    ),
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 9),
      legend.key.size = grid::unit(0.5, "cm")
    )
  )
}