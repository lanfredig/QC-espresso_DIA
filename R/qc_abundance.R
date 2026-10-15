###############################################################################
# File: qc_abundance.R
#
# Purpose:
#   Quantify and visualize protein-group signal composition across runs,
#   annotate single- and multi-protein groups using the centralized FASTA
#   reference, and generate contaminant-specific signal profiles when
#   recognized contaminant annotations are available.
#
# Inputs (expected to exist before sourcing):
#   - report      : Filtered DIA-NN report data frame
#   - fasta_ref   : Parsed FASTA reference table created by io_load.R
#   - cb_palette  : Global color palette defined in globals.R
#   - legend_fill : Standardized fill-legend helper defined in globals.R
#
# Outputs (objects defined):
#   - protein_group_annotation           : Protein-group annotation summary
#   - report_labeled                     : DIA-NN report with protein-group annotations
#   - multi_protein_annotation_summary   : Annotation summary for multi-protein groups
#   - protein_group_signal               : Aggregated protein-group signal per run
#   - top_gene_labels                    : Top annotation labels by total signal
#   - abundance                          : Relative signal composition per run
#   - p_abundance                        : Protein-group signal-composition plot
#   - CONTAMINANT_PREFIXES_AVAILABLE     : Availability flag for recognized contaminant prefixes
#   - CONTAMINANT_SIGNAL_AVAILABLE       : Availability flag for positive contaminant signal
#   - contaminant_abundance              : Relative contaminant signal composition per run
#   - contaminant_fraction               : Contaminant fraction of total signal per run
#   - top_contaminants                   : Top contaminant annotation labels by total signal
#   - p_contaminants                     : Contaminant signal-composition plot or NULL
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tidyr
#   - forcats
#   - ggplot2
#   - scales
#   - grDevices
#   - grid
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - Quantitative values remain assigned to the original DIA-NN
#     Protein.Group.
#   - Multi-protein groups are expanded only for FASTA annotation and are
#     aggregated back to one annotation row before joining to quantitative
#     data.
#   - Quantitative signal is never duplicated across individual accessions.
#   - Contaminants are recognized when CON_ or cRAP- occurs at the beginning
#     of any accession within a protein group.
#   - The contaminant plot is generated only when recognized contaminant
#     prefixes and positive contaminant Precursor.Quantity values are present.
###############################################################################

## 1. Freeze run order ---------------------------------------------------------
run_levels <- unique(report$Run)

## 2. Build accession-level FASTA annotation map -------------------------------
gene_map_available <- (
  exists("fasta_ref") &&
    is.data.frame(fasta_ref) &&
    nrow(fasta_ref) > 0 &&
    all(c("Protein.Group", "Header") %in% names(fasta_ref))
)

if (isTRUE(gene_map_available)) {
  gene_map <- fasta_ref %>%
    dplyr::transmute(
      Accession = trimws(Protein.Group),
      FASTA_Matched = TRUE,
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
    dplyr::filter(!is.na(Accession), nzchar(Accession)) %>%
    dplyr::distinct(Accession, .keep_all = TRUE)

  message("Gene annotation derived from centralized fasta_ref object")
} else {
  gene_map <- data.frame(
    Accession = character(0),
    FASTA_Matched = logical(0),
    Gene_from_fasta = character(0),
    stringsAsFactors = FALSE
  )

  message(
    paste0(
      "[ABUNDANCE] fasta_ref not available or incomplete; ",
      "using Protein.Group as fallback label"
    )
  )
}

## 3. Annotate single- and multi-protein groups --------------------------------
protein_group_annotation <- report %>%
  dplyr::filter(
    !is.na(Protein.Group),
    nzchar(trimws(Protein.Group))
  ) %>%
  dplyr::distinct(Protein.Group) %>%
  dplyr::mutate(Accession = Protein.Group) %>%
  tidyr::separate_rows(Accession, sep = ";") %>%
  dplyr::mutate(Accession = trimws(Accession)) %>%
  dplyr::filter(!is.na(Accession), nzchar(Accession)) %>%
  dplyr::left_join(gene_map, by = "Accession") %>%
  dplyr::mutate(
    FASTA_Matched = tidyr::replace_na(FASTA_Matched, FALSE)
  ) %>%
  dplyr::group_by(Protein.Group) %>%
  dplyr::summarise(
    Accession_Count = dplyr::n_distinct(Accession),
    Accessions = paste(sort(unique(Accession)), collapse = ";"),
    Matched_Accession_Count = dplyr::n_distinct(
      Accession[FASTA_Matched]
    ),
    Matched_Accessions = paste(
      sort(unique(Accession[FASTA_Matched])),
      collapse = ";"
    ),
    Unmatched_Accessions = paste(
      sort(unique(Accession[!FASTA_Matched])),
      collapse = ";"
    ),
    Gene_Annotations = paste(
      sort(unique(Gene_from_fasta[
        !is.na(Gene_from_fasta) & nzchar(Gene_from_fasta)
      ])),
      collapse = " | "
    ),
    Distinct_Gene_Count = dplyr::n_distinct(
      Gene_from_fasta[
        !is.na(Gene_from_fasta) & nzchar(Gene_from_fasta)
      ]
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Annotation_Status = dplyr::case_when(
      Matched_Accession_Count == 0L ~ "NO_FASTA_MATCH",
      Matched_Accession_Count < Accession_Count ~ "PARTIAL_FASTA_MATCH",
      Accession_Count == 1L ~ "SINGLE_PROTEIN",
      Distinct_Gene_Count <= 1L ~ "CONSISTENT_ANNOTATION",
      TRUE ~ "MULTIPLE_ANNOTATIONS"
    ),
    Annotation_Label = dplyr::case_when(
      Annotation_Status == "NO_FASTA_MATCH" ~ Protein.Group,
      Annotation_Status == "PARTIAL_FASTA_MATCH" & nzchar(Gene_Annotations) ~
        paste0(Gene_Annotations, " [partial]"),
      Annotation_Status == "PARTIAL_FASTA_MATCH" ~
        paste0(Matched_Accessions, " [partial]"),
      nzchar(Gene_Annotations) ~ Gene_Annotations,
      nzchar(Matched_Accessions) ~ Matched_Accessions,
      TRUE ~ Protein.Group
    ),
    Gene_Label = dplyr::if_else(
      !is.na(Annotation_Label) & nzchar(Annotation_Label),
      Annotation_Label,
      Protein.Group
    )
  )

## 4. Join annotation without duplicating quantitative rows --------------------
report_labeled <- report %>%
  dplyr::left_join(
    protein_group_annotation,
    by = "Protein.Group",
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    Accession_Count = tidyr::replace_na(Accession_Count, 1L),
    Matched_Accession_Count = tidyr::replace_na(
      Matched_Accession_Count,
      0L
    ),
    Annotation_Status = dplyr::if_else(
      is.na(Annotation_Status),
      "NO_FASTA_MATCH",
      Annotation_Status
    ),
    Gene_Label = dplyr::if_else(
      !is.na(Gene_Label) & nzchar(Gene_Label),
      Gene_Label,
      Protein.Group
    )
  )

## 5. Verify annotation join invariants ----------------------------------------
if (nrow(report_labeled) != nrow(report)) {
  stop(
    paste0(
      "[ABUNDANCE] Protein-group annotation changed the number of DIA-NN ",
      "report rows. Quantitative values may have been duplicated."
    )
  )
}

if (
  dplyr::n_distinct(report_labeled$Protein.Group) !=
    dplyr::n_distinct(report$Protein.Group)
) {
  stop(
    paste0(
      "[ABUNDANCE] One or more protein groups were lost during FASTA ",
      "annotation."
    )
  )
}

## 6. Report multi-protein annotation summary ----------------------------------
multi_protein_annotation_summary <- protein_group_annotation %>%
  dplyr::filter(Accession_Count > 1L) %>%
  dplyr::count(Annotation_Status, name = "Protein_Group_Count") %>%
  dplyr::arrange(Annotation_Status)

if (nrow(multi_protein_annotation_summary) > 0L) {
  message(
    "Multi-protein groups annotated: ",
    sum(multi_protein_annotation_summary$Protein_Group_Count)
  )

  for (summary_index in seq_len(nrow(multi_protein_annotation_summary))) {
    message(
      "  ",
      multi_protein_annotation_summary$Annotation_Status[summary_index],
      ": ",
      multi_protein_annotation_summary$Protein_Group_Count[summary_index]
    )
  }
}

## 7. Protein-group signal aggregation -----------------------------------------
protein_group_signal <- report_labeled %>%
  dplyr::filter(!is.na(Precursor.Quantity), Precursor.Quantity > 0) %>%
  dplyr::group_by(
    Run,
    Protein.Group,
    Gene_Label,
    Annotation_Status
  ) %>%
  dplyr::summarise(
    Signal = sum(Precursor.Quantity),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Run = factor(Run, levels = run_levels))

## 8. Verify signal preservation -----------------------------------------------
signal_before_annotation <- report %>%
  dplyr::filter(!is.na(Precursor.Quantity), Precursor.Quantity > 0) %>%
  dplyr::group_by(Run, Protein.Group) %>%
  dplyr::summarise(
    Signal_Before = sum(Precursor.Quantity),
    .groups = "drop"
  )

signal_after_annotation <- protein_group_signal %>%
  dplyr::group_by(Run, Protein.Group) %>%
  dplyr::summarise(
    Signal_After = sum(Signal),
    .groups = "drop"
  )

signal_validation <- signal_before_annotation %>%
  dplyr::full_join(
    signal_after_annotation,
    by = c("Run", "Protein.Group")
  ) %>%
  dplyr::mutate(
    Signal_Before = tidyr::replace_na(Signal_Before, 0),
    Signal_After = tidyr::replace_na(Signal_After, 0),
    Signal_Difference = Signal_After - Signal_Before
  )

if (
  nrow(signal_validation) == 0L ||
    any(abs(signal_validation$Signal_Difference) > sqrt(.Machine$double.eps))
) {
  stop(
    paste0(
      "[ABUNDANCE] Quantitative signal changed during protein-group ",
      "annotation."
    )
  )
}

rm(
  signal_before_annotation,
  signal_after_annotation,
  signal_validation
)

## 9. Global top contributors --------------------------------------------------
top_gene_labels <- protein_group_signal %>%
  dplyr::group_by(Gene_Label) %>%
  dplyr::summarise(
    TotalSignal = sum(Signal),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(TotalSignal), Gene_Label) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(Gene_Label)

## 10. Relative abundance composition ------------------------------------------
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
    Plot_Label = forcats::fct_relevel(
      Plot_Label,
      "Others",
      after = Inf
    )
  )

gene_palette <- c(
  setNames(
    grDevices::colorRampPalette(cb_palette)(length(top_gene_labels)),
    top_gene_labels
  ),
  Others = "grey85"
)

## 11. Plot: relative signal composition ---------------------------------------
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
    subtitle = paste0(
      "Largest contributions stacked at the bottom; ",
      "aggregated remainder on top"
    ),
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

## 12. Detect recognized contaminant annotations -------------------------------
recognized_contaminant_pattern <- "(^|;\\s*)(cRAP-|CON_)"

CONTAMINANT_PREFIXES_AVAILABLE <- report %>%
  dplyr::filter(
    !is.na(Protein.Group),
    grepl(
      recognized_contaminant_pattern,
      Protein.Group,
      perl = TRUE
    )
  ) %>%
  nrow() >
  0

contaminant_pg_signal <- protein_group_signal %>%
  dplyr::filter(
    !is.na(Protein.Group),
    grepl(
      recognized_contaminant_pattern,
      Protein.Group,
      perl = TRUE
    )
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

## 13. Contaminant composition or conditional warning --------------------------
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
  ## 13.1 Global top contaminant contributors ---------------------------------
  top_contaminants <- contaminant_pg_signal %>%
    dplyr::group_by(Gene_Label) %>%
    dplyr::summarise(TotalSignal = sum(Signal), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(TotalSignal), Gene_Label) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::pull(Gene_Label)

  contaminant_levels <- c(top_contaminants, "Others")

  ## 13.2 Relative contaminant composition ------------------------------------
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
      Plot_Label = factor(Plot_Label, levels = contaminant_levels),
      Plot_Label = forcats::fct_relevel(
        Plot_Label,
        "Others",
        after = Inf
      )
    )

  ## 13.3 Contaminant fraction per run ----------------------------------------
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

  ## 13.4 Contaminant color palette -------------------------------------------
  gene_palette_cont <- c(
    setNames(
      grDevices::colorRampPalette(cb_palette)(length(top_contaminants)),
      top_contaminants
    ),
    Others = "grey80"
  )

  ## 13.5 Plot: contaminant composition ---------------------------------------
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
      subtitle = paste0(
        "Largest contributors stacked at the bottom; ",
        "aggregated remainder on top"
      ),
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

## 14. Clean exit --------------------------------------------------------------
invisible(NULL)
