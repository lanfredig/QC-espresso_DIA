###############################################################################
# File: qc_coverage.R
#
# Purpose:
#   Compute protein sequence coverage from peptide evidence, export coverage
#   results for proteins with mapped peptide evidence, prepare a Top 100
#   coverage summary, and visualize the distribution of positive sequence
#   coverage values when sufficient data are available.
#
# Inputs (expected to exist before sourcing):
#   - fasta_ref      : Parsed FASTA reference table created by io_load.R
#   - peptides_clean : Peptide lists grouped by protein group, created by io_load.R
#   - results_dir    : Canonical output directory created by io_load.R
#
# Outputs (objects defined):
#   - calculate_coverage         : Helper function for calculating covered residues
#   - coverage_available         : Availability flag for the required input objects
#   - coverage_results           : Sequence coverage results for all FASTA proteins
#   - coverage_results_detected  : Coverage results for proteins with mapped residues
#   - coverage_results_top100    : Top 100 proteins ranked by sequence coverage
#   - coverage_results_graph     : Positive sequence coverage values used for plotting
#   - output_csv                 : Path to the exported protein coverage CSV file
#   - p_coverage                 : Protein sequence coverage density plot or NULL
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - purrr
#   - stringr
#   - ggplot2
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - Peptides are mapped to protein sequences by exact fixed-string matching.
#   - All occurrences of each unique peptide are included when calculating
#     covered residues.
#   - Overlapping peptide mappings do not duplicate covered residue counts.
#   - The complete in-memory coverage table includes all proteins in fasta_ref.
#   - The exported CSV includes only proteins with at least one covered residue.
#   - The Top 100 summary is ranked by descending sequence coverage percentage.
#   - The density plot includes only proteins with positive sequence coverage
#     and is generated when at least two such proteins are available.
###############################################################################
`
## 1. Protein sequence coverage helper ----------------------------------------

calculate_coverage <- function(seq, peps) {
  if (is.null(seq) || length(seq) != 1L || is.na(seq) || !nzchar(seq)) {
    return(0L)
  }

  if (is.null(peps) || length(peps) == 0L) {
    return(0L)
  }

  mask <- integer(nchar(seq))

  for (p in unique(peps)) {
    if (is.na(p) || !nzchar(p)) {
      next
    }

    locations <- stringr::str_locate_all(
      seq,
      stringr::fixed(p)
    )[[1]]

    if (!is.null(locations) && nrow(locations) > 0L) {
      for (location_index in seq_len(nrow(locations))) {
        mask[
          locations[location_index, 1]:locations[location_index, 2]
        ] <- 1L
      }
    }
  }

  sum(mask)
}

## 2. Input validation ---------------------------------------------------------

coverage_available <- TRUE
coverage_results <- NULL
coverage_results_detected <- NULL
coverage_results_top100 <- NULL
output_csv <- NULL
p_coverage <- NULL

if (!exists("fasta_ref")) {
  cat("[COVERAGE] Skipping: fasta_ref not found\n")
  coverage_available <- FALSE
}

if (!exists("peptides_clean")) {
  cat("[COVERAGE] Skipping: peptides_clean not found\n")
  coverage_available <- FALSE
}

if (!exists("results_dir")) {
  cat("[COVERAGE] Skipping: results_dir not found\n")
  coverage_available <- FALSE
}

## 3. Compute, export, and summarize protein sequence coverage ------------------------------------

if (isTRUE(coverage_available)) {
  coverage_results <- fasta_ref %>%
    dplyr::left_join(
      peptides_clean,
      by = "Protein.Group",
      relationship = "one-to-one"
    ) %>%
    dplyr::mutate(
      Aminoacids_Covered = purrr::map2_int(
        Sequence,
        Peptide_List,
        calculate_coverage
      ),
      Coverage_Pct = dplyr::if_else(
        Total_Len > 0,
        pmin(
          100,
          pmax(
            0,
            100 * Aminoacids_Covered / Total_Len
          )
        ),
        0
      ),
      Peptides = purrr::map_chr(
        Peptide_List,
        ~ if (is.null(.x) || length(.x) == 0L) {
          NA_character_
        } else {
          paste(sort(unique(.x)), collapse = ";")
        }
      )
    ) %>%
    dplyr::select(
      Protein.Group,
      Protein.Name,
      Total_Len,
      Aminoacids_Covered,
      Coverage_Pct,
      Peptides
    )

  if (nrow(coverage_results) == 0L) {
    stop(
      "Protein coverage table is empty; computation failed",
      call. = FALSE
    )
  }

  coverage_results_detected <- coverage_results %>%
    dplyr::filter(Aminoacids_Covered > 0L)

  output_csv <- file.path(
    results_dir,
    paste0(
      "protein_coverage_with_peptides_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".csv"
    )
  )

  write.csv(
    coverage_results_detected,
    output_csv,
    row.names = FALSE
  )

  message("Protein coverage CSV saved to: ", output_csv)

  coverage_results_top100 <- coverage_results %>%
    dplyr::select(-Peptides) %>%
    dplyr::arrange(
      dplyr::desc(Coverage_Pct),
      Protein.Group
    ) %>%
    dplyr::slice_head(n = 100)

  coverage_results_graph <- coverage_results %>%
    dplyr::filter(
      !is.na(Coverage_Pct),
      Coverage_Pct > 0
    )

  if (nrow(coverage_results_graph) >= 2L) {
    p_coverage <- ggplot2::ggplot(
      coverage_results_graph,
      ggplot2::aes(x = Coverage_Pct)
    ) +
      ggplot2::geom_density(
        fill = "#56B4E9",
        color = "black",
        alpha = 0.8
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq(0, 100, 10),
        limits = c(-1, 101)
      ) +
      ggplot2::labs(
        title = "Proteome Sequence Coverage Distribution",
        x = "Sequence coverage (%)",
        y = "Density"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA),
        plot.title = ggplot2::element_text(face = "bold"),
        aspect.ratio = 0.8
      )

    print(p_coverage)
  } else {
    message(
      "[COVERAGE] Skipping coverage density plot: fewer than two proteins ",
      "have positive sequence coverage"
    )
  }
} else {
  message("[COVERAGE] Section skipped due to missing inputs")
}

invisible(NULL)
