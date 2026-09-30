###############################################################################
# File: globals.R
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