###############################################################################
# File: qc_identification.R
#
# Purpose:
#   Quantify identification depth across runs for distinct modified sequences,
#   precursors, and protein groups, and visualize per-run counts with global
#   mean, maximum, and total reference values.
#
# Inputs (expected to exist before sourcing):
#   - report : Filtered DIA-NN report data frame
#
# Outputs (objects defined):
#   - id_depth_long    : Long-format identification counts per run and level
#   - avgmax           : Global mean and maximum per-run counts by level
#   - total_ids        : Total distinct identifiers across all runs by level
#   - collapsed_labels : Reference statistics and annotation positions by level
#   - mid_runs         : Runs used to position reference-statistic labels
#   - p_ident          : Identification-depth bar plot
#   - id_depth_summary : Per-run identification summary table
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tidyr
#   - ggplot2
#   - scales
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - Identification depth is calculated independently for distinct
#     Modified.Sequence, Precursor.Id, and Protein.Group values.
#   - Per-run mean and maximum counts are calculated across all retained runs.
#   - Total counts represent distinct identifiers across the complete filtered
#     DIA-NN report.
#   - The three identification levels are reported independently and do not
#     imply hierarchy or aggregation between identifier types.
###############################################################################

## 1. Identification counts (long format) --------------------------------------

id_depth_long <- report  %>% 
  dplyr::select(
    Run,
    Precursor.Id,
    Modified.Sequence,
    Protein.Group
  ) %>%
  tidyr::pivot_longer(
    cols       = -Run,
    names_to   = "Level",
    values_to  = "ID"
  ) %>%
  dplyr::group_by(Run, Level) %>%
  dplyr::distinct(ID) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

id_depth_long$Level <- factor(
  id_depth_long$Level,
  levels = c(
    "Modified.Sequence",
    "Precursor.Id",
    "Protein.Group"
  )
)


## 2. Reference statistics across runs ------------------------------------------

avgmax <- id_depth_long %>%
  dplyr::group_by(Level) %>%
  dplyr::summarise(
    mean_n = mean(n),
    max_n  = max(n),
    .groups = "drop"
  )


## 2b. Global TOTAL identification counts ---------------------------------------

total_ids <- report %>%
  dplyr::summarise(
    `Modified.Sequence` = dplyr::n_distinct(Modified.Sequence),
    `Precursor.Id`      = dplyr::n_distinct(Precursor.Id),
    `Protein.Group`     = dplyr::n_distinct(Protein.Group)
  ) %>%
  tidyr::pivot_longer(
    cols      = everything(),
    names_to  = "Level",
    values_to = "total_n"
  )

total_ids$Level <- factor(
  total_ids$Level,
  levels = levels(id_depth_long$Level)
)


collapsed_labels <- avgmax %>%
  dplyr::left_join(total_ids, by = "Level") %>%
  dplyr::mutate(
    y_base = pmax(mean_n, max_n, total_n)
  )


## 3. Visualization -------------------------------------------------------------
mid_runs <- id_depth_long %>% 
  dplyr::group_by(Level) %>% 
  dplyr::summarise(
    middle_run = dplyr::first(Run[ceiling(dplyr::n() / 2)]),
    .groups = "drop"
  )

collapsed_labels <- collapsed_labels %>% 
  dplyr::left_join(mid_runs, by = "Level")

p_ident <- ggplot2::ggplot(
  id_depth_long,
  ggplot2::aes(x = Run, y = n)
) +
  ggplot2::geom_col(
    fill  = "grey85",
    color = "black"
  ) +
  ggplot2::facet_wrap(
    ~ Level,
    scales = "free_y",
    ncol   = 1
  ) +
  
  ggplot2::geom_hline(
    data = avgmax,
    ggplot2::aes(yintercept = mean_n),
    color     = "blue",
    linetype  = "dotted",
    linewidth = 1
  ) +
  ggplot2::geom_hline(
    data = avgmax,
    ggplot2::aes(yintercept = max_n),
    color     = "red",
    linetype  = "dotted",
    linewidth = 0.8
  ) +
  ggplot2::geom_hline(
    data = total_ids,
    ggplot2::aes(yintercept = total_n),
    color     = "purple",
    linewidth = 1
  ) +
  
  ggplot2::geom_text(
    data = collapsed_labels,
    ggplot2::aes(
      x = middle_run,
      y = y_base * 1.1,
      label = paste0("Avg: ", scales::comma(round(mean_n)))
    ),
    hjust = 0.5,
    color = "blue",
    size  = 3,
    inherit.aes = FALSE
  ) +
  
  ggplot2::geom_text(
    data = collapsed_labels,
    ggplot2::aes(
      x = middle_run,
      y = y_base * 1.25,
      label = paste0("Max: ", scales::comma(max_n))
    ),
    hjust = 0.5,
    color = "red",
    size  = 3,
    inherit.aes = FALSE
  ) +
  
  ggplot2::geom_text(
    data = collapsed_labels,
    ggplot2::aes(
      x = middle_run,
      y = y_base * 1.4,
      label = paste0("Total: ", scales::comma(total_n))
    ),
    hjust = 0.5,
    color = "purple",
    size  = 3,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.22))
  ) +
  ggplot2::labs(
    title = "Identification Depth per Run",
    x     = "Run",
    y     = "Number of Identifications"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    plot.margin = ggplot2::margin(
      t = 5,
      r = 5,
      b = 20,
      l = 25
    )
  )

print(p_ident)

## 4. Identification depth summary table ---------------------------------------

id_depth_summary <- report %>%
  dplyr::select(
    Run,
    Precursor.Id,
    Modified.Sequence,
    Protein.Group
  ) %>%
  dplyr::distinct() %>%
  dplyr::group_by(Run) %>%
  dplyr::summarise(
    Precursors     = dplyr::n_distinct(Precursor.Id),
    Peptides       = dplyr::n_distinct(Modified.Sequence),
    Protein.Groups = dplyr::n_distinct(Protein.Group),
    .groups = "drop"
  )


## 5. Clean exit ---------------------------------------------------------------

invisible(NULL)