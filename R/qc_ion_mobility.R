###############################################################################
# File: qc_ion_mobility.R
#
# Purpose:
#   Visualize ion mobility patterns across runs using the joint distribution
#   of precursor m/z and the available ion mobility field, and compare the
#   linear trends of observed ion mobility and DIA-NN iIM values when iIM is
#   available.
#
# Inputs (expected to exist before sourcing):
#   - report         : Filtered DIA-NN report data frame
#   - PANELS_PER_ROW : Number of panels per row defined in globals.R
#   - AR_DENSITY     : Aspect-ratio preset for density-like plots defined in
#                      globals.R
#
# Outputs (objects defined):
#   - im_candidates           : Candidate names for the observed ion mobility field
#   - im_col                  : Ion mobility column selected from the report
#   - im_df                   : Valid precursor m/z, ion mobility, and iIM values
#   - run_has_valid_im_trend  : Helper function for assessing whether a run
#                               supports linear trend fitting
#   - p_im_density            : Precursor m/z versus ion mobility density plot
#   - valid_trend_runs_df     : Trend-fitting availability by run
#   - trend_runs_vec          : Runs meeting the trend-fitting requirements
#   - im_trend_df             : Ion mobility values from runs valid for trend fitting
#   - df_long                 : Long-format observed ion mobility and iIM values
#   - metrics_df              : Linear-model coefficients and squared correlation
#                               between observed ion mobility and iIM by run
#   - p_trend                 : Observed ion mobility versus iIM trend plot
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tidyr
#   - ggplot2
#   - tibble
#   - rlang
#   - stats
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - The observed ion mobility field is selected from IM, Ion.Mobility,
#     Precursor.Ion.Mobility, or 1/K0, using the first available field in that
#     order.
#   - DIA-NN iIM values are analyzed only when the iIM column is available.
#   - The script does not use Predicted.IM.
#   - Ion mobility values equal to zero are excluded because datasets from
#     instruments without ion mobility may encode unavailable values as zero.
#   - Trend fitting requires at least five complete observations and at least
#     two distinct values for Precursor.Mz, observed ion mobility, and iIM.
#   - The reported R² is the squared Pearson correlation between observed ion
#     mobility and iIM, rather than the coefficient of determination from a
#     regression of one field against the other.
#   - Plot generation is skipped when no supported ion mobility column or no
#     valid non-zero ion mobility values are available.
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
    
    
    ## 2.4 Identify runs valid for observed ion mobility and iIM comparison ----
    
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
    
    
    ## 2.5 Plot 2: observed ion mobility and iIM trends ------------------------
    
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
          title    = "Observed Ion Mobility and DIA-NN iIM Trends",
          subtitle = "Separate linear trends across precursor m/z",
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
