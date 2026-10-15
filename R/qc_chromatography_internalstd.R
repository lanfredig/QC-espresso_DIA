###############################################################################
# File: qc_chromatography_internalstd.R
#
# Purpose:
#   Generate QC plots for predefined trypsin autolysis precursors, covering
#   retention-time reproducibility and quantitative signal stability across
#   runs.
#
# Inputs (expected to exist before sourcing):
#   - report       : Filtered DIA-NN report data frame containing Run,
#                    Precursor.Id, RT, and Precursor.Quantity
#   - legend_color : Standardized color-legend helper defined in globals.R
#
# Outputs (objects defined):
#   - TRYPSIN_PRECURSORS : Predefined trypsin autolysis precursor identifiers
#   - TRYPSIN_COLORS     : Colors assigned to the monitored precursors
#   - trypsin_rt         : Retention-time values for monitored precursors
#   - p_rt_std           : Retention-time reproducibility plot
#   - trypsin_quant      : Quantitative values for monitored precursors
#   - p_quant_std        : Quantitative signal-stability plot
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - ggplot2
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - The script monitors the predefined trypsin autolysis precursors
#     VATVSLPR2, LSSPATLNSR2, and LGEHNIDVLEGNEQFINAAK3.
#   - Colors are explicitly assigned to each monitored precursor for
#     consistency across plots.
#   - Plot objects are created only when matching precursor entries are
#     available in the filtered DIA-NN report.
###############################################################################

## 1. Internal standard RT and signal stability (trypsin) ---------------------


TRYPSIN_PRECURSORS <- c(
  "VATVSLPR2",
  "LSSPATLNSR2",
  "LGEHNIDVLEGNEQFINAAK3"
)

TRYPSIN_COLORS <- c(
  "VATVSLPR2"             = "#0072B2",
  "LSSPATLNSR2"           = "#E69F00",
  "LGEHNIDVLEGNEQFINAAK3" = "#009E73"
)

## 2. Retention‑time reproducibility ------------------------------------------
trypsin_rt <- report %>% 
  dplyr::filter(Precursor.Id %in% TRYPSIN_PRECURSORS) %>% 
  dplyr::select(Run, Precursor.Id, RT)

if (nrow(trypsin_rt) > 0) {
  
  rt_max    <- ceiling(max(report$RT, na.rm = TRUE) * 1.05)
  rt_breaks <- seq(0, rt_max, by = 5)
  
  p_rt_std <- ggplot2::ggplot(
    trypsin_rt,
    ggplot2::aes(
      x     = Run,
      y     = RT,
      group = Precursor.Id,
      color = Precursor.Id
    )
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3, alpha = 0.5) +
    ggplot2::scale_y_continuous(
      limits = c(0, rt_max),
      breaks = rt_breaks
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(add = 1)
    ) +
    ggplot2::scale_color_manual(values = TRYPSIN_COLORS) +
    ggplot2::labs(
      title = "RT Reproducibility — Trypsin Autolysis Precursorss",
      x     = "Run",
      y     = "Retention Time (min)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      plot.margin      = ggplot2::margin(t = 10, r = 10, b = 40, l = 20)
    ) +
    legend_color() +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    ) +
    ggplot2::theme(
      legend.position  = "top",
      legend.direction = "horizontal"
    )
  
  print(p_rt_std)
  
} else {
  message("No trypsin autolysis peptides detected for RT analysis.")
}



## 3. Quantitative signal stability -------------------------------------------

trypsin_quant <- report %>% 
  dplyr::filter(Precursor.Id %in% TRYPSIN_PRECURSORS) %>% 
  dplyr::select(Run, Precursor.Id, Precursor.Quantity)

if (nrow(trypsin_quant) > 0) {
  
  p_quant_std <- ggplot2::ggplot(
    trypsin_quant,
    ggplot2::aes(
      x    = Run,
      y    = Precursor.Quantity,
      fill = Precursor.Id
    )
  ) +
    ggplot2::geom_col(color = "black", width = 0.7) +
    ggplot2::facet_wrap(
      ~ Precursor.Id,
      scales = "free_y",
      ncol   = 1
    ) +
    ggplot2::scale_fill_manual(values = TRYPSIN_COLORS) +
    ggplot2::labs(
      title = "Quantitative Stability — Trypsin Autolysis Precursors",
      x     = "Run",
      y     = "Precursor Quantity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  print(p_quant_std)
  
} else {
  message("No trypsin autolysis peptides detected for quantitative analysis.")
}

