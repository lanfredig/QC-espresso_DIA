###############################################################################
# File: qc_ion_mobility.R
#
# Purpose:
#   Evaluate ion mobility performance across runs using:
#   1. Precursor m/z vs ion mobility density structure
#   2. Observed vs predicted ion mobility trends with per-run QC metrics
#
# Inputs (expected to exist before sourcing):
#   - report                  : Filtered DIA-NN report data frame
#   - PANELS_PER_ROW          : Number of panels per row in facetted plots
#   - AR_DENSITY              : Aspect ratio constant for density-like plots
#   - HAS_VALID_ION_MOBILITY  : Logical flag computed upstream in QMD
#
# Outputs:
#   - Ion mobility density plot rendered to the active graphics device
#   - Observed vs library predicted ion mobility trend plot rendered when applicable
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tidyr
#   - ggplot2
#   - tibble
#   - rlang
#
# Notes:
#   - This script assumes ion mobility validity has already been determined
#     upstream in the QMD.
#   - The script is source()-safe and does not use return() at top level.
#   - Plot generation is skipped only if no ion mobility column is present,
#     as a defensive safeguard against accidental standalone execution.
#   - IM values equal to zero are excluded from plotting because datasets from
#     instruments without ion mobility may encode absent IM as 0.
###############################################################################


## 1. Detect ion mobility column -----------------------------------------------

im_candidates <- c(
  "IM",
  "Ion.Mobility",
  "Precursor.Ion.Mobility",
  "1/K0"
)

im_col <- im_candidates[im_candidates %in% colnames(report)][1]


## 2. Execute only when an ion mobility column is present ----------------------

if (!is.na(im_col)) {
  
  ## 2.1 Prepare data ----------------------------------------------------------
  
  im_df <- report %>% 
    dplyr::select(
      Run,
      Precursor.Mz,
      !!rlang::sym(im_col),
      dplyr::any_of("iIM")
    ) %>% 
    dplyr::rename(IM = !!rlang::sym(im_col)) %>% 
    dplyr::filter(
      !is.na(IM),
      IM != 0,
      !is.na(Precursor.Mz)
    )
  
  if (!"iIM" %in% colnames(im_df)) {
    im_df <- im_df %>% 
      dplyr::mutate(iIM = NA_real_)
  }
  
  
  ## 2.2 Helper: determine whether a run supports trend fitting ----------------
  
  run_has_valid_im_trend <- function(df) {
    
    df_valid <- df %>% 
      dplyr::filter(
        !is.na(Precursor.Mz),
        !is.na(IM),
        IM != 0,
        !is.na(iIM)
      )
    
    if (nrow(df_valid) < 5) {
      return(FALSE)
    }
    
    if (dplyr::n_distinct(df_valid$Precursor.Mz) < 2) {
      return(FALSE)
    }
    
    if (dplyr::n_distinct(df_valid$IM) < 2) {
      return(FALSE)
    }
    
    if (dplyr::n_distinct(df_valid$iIM) < 2) {
      return(FALSE)
    }
    
    TRUE
  }
  
  
  ## 2.3 Plot 1: ion mobility vs precursor m/z density -------------------------
  
  if (nrow(im_df) > 0) {
    
    p_im_density <- ggplot2::ggplot(
      im_df,
      ggplot2::aes(
        x = Precursor.Mz,
        y = IM
      )
    ) +
      ggplot2::geom_bin2d(bins = 120) +
      ggplot2::scale_fill_gradientn(
        colours = c(
          "#56B4E9",
          "white",
          "yellow",
          "red"
        ),
        name  = "Precursor count",
        trans = "sqrt"
      ) +
      ggplot2::facet_wrap(
        ~ Run,
        ncol = PANELS_PER_ROW
      ) +
      ggplot2::labs(
        title    = "Ion Mobility vs Precursor m/z",
        subtitle = "Precursor density in ion mobility separation space",
        x        = "Precursor m/z",
        y        = "Ion Mobility"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        aspect.ratio = AR_DENSITY
      )
    
    print(p_im_density)
    
    
    ## 2.4 Identify runs valid for observed vs predicted trend comparison ------
    
    valid_trend_runs_df <- im_df %>% 
      dplyr::group_by(Run) %>% 
      dplyr::group_modify(~{
        tibble::tibble(
          valid_for_trend = run_has_valid_im_trend(.x)
        )
      }) %>% 
      dplyr::ungroup()
    
    trend_runs_vec <- valid_trend_runs_df %>% 
      dplyr::filter(valid_for_trend) %>% 
      dplyr::pull(Run)
    
    im_trend_df <- im_df %>% 
      dplyr::filter(Run %in% trend_runs_vec)
    
    
    ## 2.5 Plot 2: observed vs predicted ion mobility trends -------------------
    
    if (nrow(im_trend_df) > 0) {
      
      df_long <- im_trend_df %>% 
        dplyr::select(Run, Precursor.Mz, IM, iIM) %>% 
        tidyr::pivot_longer(
          cols      = c("IM", "iIM"),
          names_to  = "Type",
          values_to = "Value"
        )
      
      metrics_df <- im_trend_df %>% 
        dplyr::group_by(Run) %>% 
        dplyr::group_modify(~{
          
          fit_obs  <- stats::lm(IM ~ Precursor.Mz, data = .x)
          fit_pred <- stats::lm(iIM ~ Precursor.Mz, data = .x)
          
          slope_obs      <- unname(stats::coef(fit_obs)[2])
          intercept_obs  <- unname(stats::coef(fit_obs)[1])
          slope_pred     <- unname(stats::coef(fit_pred)[2])
          intercept_pred <- unname(stats::coef(fit_pred)[1])
          
          r2 <- suppressWarnings(
            stats::cor(.x$IM, .x$iIM, use = "complete.obs")^2
          )
          
          label <- paste0(
            "Observed: y = ",
            signif(slope_obs, 3), "x + ",
            signif(intercept_obs, 3), "\n",
            "Predicted: y = ",
            signif(slope_pred, 3), "x + ",
            signif(intercept_pred, 3), "\n",
            "R² = ",
            signif(r2, 3)
          )
          
          x_min <- min(.x$Precursor.Mz, na.rm = TRUE)
          x_max <- max(.x$Precursor.Mz, na.rm = TRUE)
          
          y_all <- c(.x$IM, .x$iIM)
          y_min <- min(y_all, na.rm = TRUE)
          y_max <- max(y_all, na.rm = TRUE)
          
          tibble::tibble(
            label = label,
            x = x_min + 0.02 * (x_max - x_min),
            y = y_max - 0.02 * (y_max - y_min)
          )
        }) %>% 
        dplyr::ungroup()
      
      p_trend <- ggplot2::ggplot(
        df_long,
        ggplot2::aes(
          x = Precursor.Mz,
          y = Value,
          colour = Type
        )
      ) +
        ggplot2::geom_smooth(
          method = "lm",
          se = FALSE,
          linewidth = 1
        ) +
        ggplot2::scale_colour_manual(
          values = c(
            "IM" = "black",
            "iIM" = "red"
          ),
          labels = c(
            "IM" = "Observed",
            "iIM" = "Predicted"
          ),
          name = NULL
        ) +
        ggplot2::geom_text(
          data = metrics_df,
          ggplot2::aes(
            x = x,
            y = y,
            label = label
          ),
          inherit.aes = FALSE,
          hjust = 0,
          vjust = 1,
          size = 3
        ) +
        ggplot2::facet_wrap(
          ~ Run,
          ncol = PANELS_PER_ROW
        ) +
        ggplot2::labs(
          title    = "Observed vs Predicted IM Trends",
          subtitle = "Linear trend comparison",
          x        = "Precursor m/z",
          y        = "Ion Mobility"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          aspect.ratio = AR_DENSITY,
          legend.position = "top"
        )
      
      print(p_trend)
    }
    
  } else {
    message("[ION MOBILITY] Skipping: no non-zero ion mobility values available for plotting")
  }
  
} else {
  message("[ION MOBILITY] Skipping: no ion mobility column found in report")
}


## 3. Clean exit ---------------------------------------------------------------

invisible(NULL)
