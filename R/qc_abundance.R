###############################################################################
# File: qc_abundance.R
#
# Purpose:
#   Quantify and visualize protein-level signal composition across runs and
#   generate contaminant-specific profiles when recognized annotations exist.
#
# Inputs:
#   - report
#   - fasta_ref
#   - cb_palette
#   - legend_fill
#
# Outputs:
#   - p_abundance
#   - p_contaminants
#   - CONTAMINANT_PREFIXES_AVAILABLE
#   - CONTAMINANT_SIGNAL_AVAILABLE
#   - contaminant_fraction
#
# Notes:
#   - Contaminants are recognized by Protein.Group prefixes CON_ or cRAP-.
#   - If recognized contaminant entries are unavailable, the contaminant plot
#     is replaced by a Quarto warning callout.
###############################################################################

## 1. Freeze run order ---------------------------------------------------------

run_levels <- unique(report$Run)

## 2. Gene annotation from centralized FASTA object ----------------------------

gene_map_available <- (
  exists("fasta_ref") &&
    is.data.frame(fasta_ref) &&
    nrow(fasta_ref) > 0 &&
    all(c("Protein.Group", "Header") %in% names(fasta_ref))
)

if (isTRUE(gene_map_available)) {
  gene_map <- fasta_ref %>%
    dplyr::transmute(
      Protein.Group = trimws(Protein.Group),
      Gene_from_fasta = dplyr::if_else(
        grepl("GN=", Header),
        sub(".*GN=([^ ]+).*", "\\1", Header),
        NA_character_
      )
    ) %>%
    dplyr::mutate(
      Gene_from_fasta = trimws(Gene_from_fasta),
      Gene_from_fasta = dplyr::if_else(
        !is.na(Gene_from_fasta) & nzchar(Gene_from_fasta),
        Gene_from_fasta,
        NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Protein.Group), nzchar(Protein.Group)) %>%
    dplyr::distinct(Protein.Group, .keep_all = TRUE)

  message("Gene annotation derived from centralized fasta_ref object")
} else {
  gene_map <- data.frame(
    Protein.Group = character(0),
    Gene_from_fasta = character(0),
    stringsAsFactors = FALSE
  )

  message(
    "[ABUNDANCE] fasta_ref not available or incomplete; using Protein.Group as fallback label"
  )
}

report_labeled <- report %>%
  dplyr::left_join(gene_map, by = "Protein.Group") %>%
  dplyr::mutate(
    Gene_Label = dplyr::if_else(
      !is.na(Gene_from_fasta),
      Gene_from_fasta,
      Protein.Group
    )
  ) %>%
  dplyr::select(-Gene_from_fasta)

## 3. Protein-group signal aggregation -----------------------------------------

protein_group_signal <- report_labeled %>%
  dplyr::filter(!is.na(Precursor.Quantity), Precursor.Quantity > 0) %>%
  dplyr::group_by(Run, Protein.Group, Gene_Label) %>%
  dplyr::summarise(Signal = sum(Precursor.Quantity), .groups = "drop") %>%
  dplyr::mutate(Run = factor(Run, levels = run_levels))

## 4. Global top contributors --------------------------------------------------

top_gene_labels <- protein_group_signal %>%
  dplyr::group_by(Gene_Label) %>%
  dplyr::summarise(TotalSignal = sum(Signal), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(TotalSignal)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(Gene_Label)

## 5. Relative abundance composition ------------------------------------------

abundance <- protein_group_signal %>%
  dplyr::mutate(
    Plot_Label = dplyr::if_else(
      Gene_Label %in% top_gene_labels,
      Gene_Label,
      "Others"
    )
  ) %>%
  dplyr::group_by(Run, Plot_Label) %>%
  dplyr::summarise(Signal = sum(Signal), .groups = "drop") %>%
  tidyr::complete(Run, Plot_Label, fill = list(Signal = 0)) %>%
  dplyr::group_by(Run) %>%
  dplyr::mutate(
    TotalSignal_Run = sum(Signal),
    Proportion = dplyr::if_else(
      TotalSignal_Run > 0,
      Signal / TotalSignal_Run,
      0
    ),
    Stack_Order = dplyr::if_else(
      Plot_Label == "Others",
      Inf,
      dplyr::dense_rank(dplyr::desc(Signal))
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-TotalSignal_Run) %>%
  dplyr::mutate(
    Run = factor(Run, levels = run_levels),
    Plot_Label = factor(Plot_Label),
    Plot_Label = forcats::fct_relevel(Plot_Label, "Others", after = Inf)
  )

gene_palette <- c(
  setNames(
    grDevices::colorRampPalette(cb_palette)(length(top_gene_labels)),
    top_gene_labels
  ),
  Others = "grey85"
)

## 6. Plot: relative signal composition ----------------------------------------

p_abundance <- ggplot2::ggplot(
  abundance,
  ggplot2::aes(
    x = Run,
    y = Proportion,
    fill = Plot_Label,
    group = interaction(Run, Stack_Order)
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_stack(reverse = TRUE)
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(values = gene_palette) +
  ggplot2::labs(
    title = "Protein and Gene Signal Composition per Run",
    subtitle = "Largest contributions stacked at the bottom; aggregated remainder on top",
    x = "Run",
    y = "Fraction of total signal (%)",
    fill = "Protein / Gene"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) +
  legend_fill()

print(p_abundance)

## 7. Detect recognized contaminant annotations --------------------------------

recognized_contaminant_pattern <- "^(cRAP-|CON_)"

CONTAMINANT_PREFIXES_AVAILABLE <- report %>%
  dplyr::filter(
    !is.na(Protein.Group),
    grepl(recognized_contaminant_pattern, Protein.Group)
  ) %>%
  nrow() > 0

contaminant_pg_signal <- protein_group_signal %>%
  dplyr::filter(
    !is.na(Protein.Group),
    grepl(recognized_contaminant_pattern, Protein.Group)
  )

CONTAMINANT_SIGNAL_AVAILABLE <- (
  isTRUE(CONTAMINANT_PREFIXES_AVAILABLE) &&
    nrow(contaminant_pg_signal) > 0 &&
    any(
      !is.na(contaminant_pg_signal$Signal) &
        contaminant_pg_signal$Signal > 0
    )
)

p_contaminants <- NULL
contaminant_abundance <- data.frame()
contaminant_fraction <- data.frame()
top_contaminants <- character(0)

## 8. Contaminant composition or conditional warning ---------------------------

if (!isTRUE(CONTAMINANT_PREFIXES_AVAILABLE)) {

  cat(
    paste0(
      '\n<div class="alert alert-warning" role="alert" ',
      'style="margin-top: 1rem; margin-bottom: 1rem;">\n',
      '<h4 class="alert-heading">Contaminant assessment unavailable</h4>\n',
      '<p>',
      'No protein groups with recognized <code>CON_</code> or ',
      '<code>cRAP-</code> contaminant prefixes were found in the filtered ',
      'DIA-NN report. Contaminants may not have been included in the search ',
      'database, may use different accession prefixes, or may not have passed ',
      'the initial confidence filter.',
      '</p>\n',
      '<p class="mb-0">',
      'The contaminant signal-composition plot was therefore not generated.',
      '</p>\n',
      '</div>\n\n'
    )
  )

} else if (!isTRUE(CONTAMINANT_SIGNAL_AVAILABLE)) {

  cat(
    paste0(
      '\n<div class="alert alert-warning" role="alert" ',
      'style="margin-top: 1rem; margin-bottom: 1rem;">\n',
      '<h4 class="alert-heading">',
      'Contaminant signal composition unavailable',
      '</h4>\n',
      '<p>',
      'Protein groups with recognized <code>CON_</code> or ',
      '<code>cRAP-</code> prefixes were found in the filtered DIA-NN report, ',
      'but no positive <code>Precursor.Quantity</code> values were available ',
      'for contaminant signal aggregation.',
      '</p>\n',
      '<p class="mb-0">',
      'The contaminant signal-composition plot was therefore not generated.',
      '</p>\n',
      '</div>\n\n'
    )
  )

} else {
  ## 8.1 Global top contaminant contributors ----------------------------------

  top_contaminants <- contaminant_pg_signal %>%
    dplyr::group_by(Gene_Label) %>%
    dplyr::summarise(TotalSignal = sum(Signal), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(TotalSignal)) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::pull(Gene_Label)

  contaminant_levels <- c(top_contaminants, "Others")

  ## 8.2 Relative contaminant composition -------------------------------------

  contaminant_abundance <- contaminant_pg_signal %>%
    dplyr::mutate(
      Plot_Label = dplyr::if_else(
        Gene_Label %in% top_contaminants,
        Gene_Label,
        "Others"
      )
    ) %>%
    dplyr::group_by(Run, Plot_Label) %>%
    dplyr::summarise(Signal = sum(Signal), .groups = "drop") %>%
    tidyr::complete(
      Run = factor(run_levels, levels = run_levels),
      Plot_Label = contaminant_levels,
      fill = list(Signal = 0)
    ) %>%
    dplyr::group_by(Run) %>%
    dplyr::mutate(
      TotalContaminantSignal = sum(Signal),
      Proportion = dplyr::if_else(
        TotalContaminantSignal > 0,
        Signal / TotalContaminantSignal,
        0
      ),
      Stack_Order = dplyr::if_else(
        Plot_Label == "Others",
        Inf,
        dplyr::dense_rank(dplyr::desc(Signal))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-TotalContaminantSignal) %>%
    dplyr::mutate(
      Run = factor(Run, levels = run_levels),
      Plot_Label = factor(Plot_Label),
      Plot_Label = forcats::fct_relevel(Plot_Label, "Others", after = Inf)
    )

  ## 8.3 Contaminant fraction per run -----------------------------------------

  total_signal <- protein_group_signal %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(TotalSignal = sum(Signal), .groups = "drop")

  contaminant_fraction <- total_signal %>%
    dplyr::left_join(
      contaminant_pg_signal %>%
        dplyr::group_by(Run) %>%
        dplyr::summarise(
          ContaminantSignal = sum(Signal),
          .groups = "drop"
        ),
      by = "Run"
    ) %>%
    dplyr::mutate(
      ContaminantSignal = tidyr::replace_na(ContaminantSignal, 0),
      ContaminantFraction = dplyr::if_else(
        TotalSignal > 0,
        ContaminantSignal / TotalSignal,
        0
      ),
      Label = dplyr::if_else(
        ContaminantSignal > 0,
        scales::percent(ContaminantFraction, accuracy = 0.1),
        "No recognized contaminant signal"
      ),
      Run = factor(Run, levels = run_levels)
    )

  ## 8.4 Contaminant color palette --------------------------------------------

  gene_palette_cont <- c(
    setNames(
      grDevices::colorRampPalette(cb_palette)(length(top_contaminants)),
      top_contaminants
    ),
    Others = "grey80"
  )

  ## 8.5 Plot: contaminant composition ----------------------------------------

  p_contaminants <- ggplot2::ggplot(
    contaminant_abundance,
    ggplot2::aes(
      y = Run,
      x = Proportion,
      fill = Plot_Label,
      group = interaction(Run, Stack_Order)
    )
  ) +
    ggplot2::geom_col(
      width = 0.6,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::geom_text(
      data = contaminant_fraction,
      ggplot2::aes(y = Run, x = 1.15, label = Label),
      inherit.aes = FALSE,
      hjust = 1,
      size = 3
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.3, 0.6, 0.9),
      labels = scales::percent
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1.2)) +
    ggplot2::scale_fill_manual(
      values = gene_palette_cont,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Contaminant Signal Composition per Run",
      subtitle = "Largest contributors stacked at the bottom; aggregated remainder on top",
      x = "Relative abundance within contaminants (%)",
      y = "Run",
      fill = "Contaminant"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = grid::unit(0.3, "cm")
    ) +
    legend_fill()

  print(p_contaminants)
}

## 9. Clean exit ---------------------------------------------------------------

invisible(NULL)
