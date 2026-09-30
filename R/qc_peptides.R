###############################################################################
# File: qc_peptides.R
#
# Purpose:
#   Summarize peptide-related quality characteristics derived from identified
#   precursors and peptide sequences, including precursor charge-state
#   distribution, peptide length distribution, and digestion efficiency as
#   assessed by missed tryptic cleavages.
###############################################################################

## 0. Shared missed-cleavage calculation --------------------------------------

if (!exists("count_missed_tryptic_cleavages", mode = "function")) {
  count_missed_tryptic_cleavages <- function(sequence) {
    cleavage_sites <- stringr::str_count(
      sequence,
      "(?<=[KR])(?!P)"
    )

    terminal_tryptic_site <- as.integer(
      stringr::str_detect(
        sequence,
        "[KR]$"
      )
    )

    pmax(
      0L,
      cleavage_sites - terminal_tryptic_site
    )
  }
}

## 1. Centralized peptide-related analysis tables ------------------------------

precursor_level <- report %>%
  dplyr::select(
    Run,
    Precursor.Id,
    Precursor.Charge
  ) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(Precursor.Id),
    !is.na(Precursor.Charge)
  ) %>%
  dplyr::distinct()

peptide_level <- report %>%
  dplyr::select(
    Run,
    Stripped.Sequence
  ) %>%
  dplyr::filter(
    !is.na(Run),
    !is.na(Stripped.Sequence),
    nzchar(Stripped.Sequence)
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    seq_length = nchar(Stripped.Sequence),
    MissedCleavages = count_missed_tryptic_cleavages(
      Stripped.Sequence
    )
  )

## 2. Precursor charge-state distribution -------------------------------------

chrg <- precursor_level %>%
  dplyr::group_by(Run, Precursor.Charge) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop_last"
  ) %>%
  dplyr::mutate(
    pct = n / sum(n)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Precursor.Charge = factor(Precursor.Charge)
  )

charge_colors <- setNames(
  grDevices::colorRampPalette(cb_palette)(
    length(levels(chrg$Precursor.Charge))
  ),
  levels(chrg$Precursor.Charge)
)

runs <- unique(chrg$Run)
run_groups <- split(runs, ceiling(seq_along(runs) / 3))

for (grp in run_groups) {
  df_subset <- chrg %>%
    dplyr::filter(Run %in% grp)

  p_charge <- ggplot2::ggplot(
    df_subset,
    ggplot2::aes(
      Precursor.Charge,
      pct,
      fill = Precursor.Charge
    )
  ) +
    ggplot2::geom_col(color = "black") +
    ggplot2::facet_wrap(
      ~ Run,
      scales = "free_y",
      ncol = 3
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent
    ) +
    ggplot2::scale_fill_manual(
      values = charge_colors
    ) +
    ggplot2::labs(
      title = "Precursor Charge-State Distribution per Run",
      subtitle = "Relative frequency computed from unique precursors per run",
      x = "Precursor charge",
      y = "Percentage (%)",
      fill = "Charge state"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none"
    )

  print(p_charge)
}

## 3. Peptide length distribution ---------------------------------------------

pep_len <- peptide_level %>%
  dplyr::select(Run, Stripped.Sequence, seq_length)

pep_len_stat <- pep_len %>%
  dplyr::group_by(Run) %>%
  dplyr::summarise(
    median_val = median(seq_length, na.rm = TRUE),
    .groups = "drop"
  )

pep_len_ymax <- pep_len %>%
  dplyr::count(Run, seq_length) %>%
  dplyr::group_by(Run) %>%
  dplyr::summarise(
    y_top = max(n),
    .groups = "drop"
  )

pep_len_stat <- pep_len_stat %>%
  dplyr::left_join(pep_len_ymax, by = "Run")

runs <- unique(pep_len$Run)
run_groups <- split(runs, ceiling(seq_along(runs) / 3))

for (grp in run_groups) {
  df_subset <- pep_len %>%
    dplyr::filter(Run %in% grp)

  stat_subset <- pep_len_stat %>%
    dplyr::filter(Run %in% grp)

  p_len <- ggplot2::ggplot(
    df_subset,
    ggplot2::aes(seq_length, fill = Run)
  ) +
    ggplot2::geom_histogram(
      binwidth = 1,
      color = "black"
    ) +
    ggplot2::geom_vline(
      data = stat_subset,
      ggplot2::aes(xintercept = median_val),
      color = "red",
      linetype = "dashed"
    ) +
    ggplot2::geom_text(
      data = stat_subset,
      ggplot2::aes(
        x = median_val,
        y = y_top,
        label = paste0("Median = ", round(median_val, 0))
      ),
      hjust = 0,
      nudge_x = 0.3,
      vjust = -0.3,
      color = "red",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(
      ~ Run,
      scales = "free_y",
      ncol = 3
    ) +
    ggplot2::scale_fill_manual(
      values = file_palette,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      title = "Peptide Length Distribution per Run",
      subtitle = "Computed from unique stripped peptide sequences per run",
      x = "Peptide length (aa)",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  print(p_len)
}

## 4. Digestion efficiency: missed cleavages ----------------------------------

missed_cleavages <- peptide_level %>%
  dplyr::group_by(Run) %>%
  dplyr::summarise(
    Total_sequences = dplyr::n(),
    With_missed = sum(MissedCleavages > 0),
    pct_Missed_Tryptic_Cleavages = With_missed / Total_sequences,
    .groups = "drop"
  )

p_miss <- ggplot2::ggplot(
  missed_cleavages,
  ggplot2::aes(
    x = Run,
    y = pct_Missed_Tryptic_Cleavages
  )
) +
  ggplot2::geom_col(
    fill = "steelblue",
    color = "black",
    width = 0.7
  ) +
  ggplot2::geom_hline(
    yintercept = median(
      missed_cleavages$pct_Missed_Tryptic_Cleavages
    ),
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey40"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  ggplot2::labs(
    title = "Digestion Efficiency per Run",
    subtitle = "Fraction of unique peptide sequences with at least one missed tryptic cleavage",
    x = "Run",
    y = "Percentage (%)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_miss)
