###############################################################################
# File: qc_coverage.R
#
# Purpose:
#   Compute protein sequence coverage using peptide evidence, export a complete
#   coverage table including peptide mappings, prepare a Top 100 summary table
#   object for downstream rendering in QMD, and visualize the global
#   distribution of sequence coverage.
#
# Inputs (expected to exist before sourcing):
#   - fasta_ref      : Parsed FASTA reference table
#   - peptides_clean : Peptide lists grouped by protein group
#   - results_dir    : Canonical output directory
#   - AR_DENSITY     : Aspect ratio constant for density plots
#
# Outputs:
#   - coverage_results              : Full protein coverage table
#   - coverage_results_detected     : Coverage table filtered to detected proteins
#   - coverage_results_top100       : Top 100 coverage summary table for QMD rendering
#   - output_csv                    : CSV file with peptide-resolved coverage information
#   - p_coverage                    : Density plot of protein sequence coverage
###############################################################################


## 1. Protein sequence coverage helper -----------------------------------------

calculate_coverage <- function(seq, peps) {
  
  if (is.null(seq) || is.na(seq)) {
    return(0L)
  }
  
  if (is.null(peps) || length(peps) == 0) {
    return(0L)
  }
  
  mask <- integer(nchar(seq))
  
  for (p in peps) {
    
    if (is.na(p) || !nzchar(p)) {
      next
    }
    
    res <- stringr::str_locate_all(
      seq,
      stringr::fixed(p)
    )[[1]]
    
    if (!is.null(res) && nrow(res) > 0) {
      for (i in seq_len(nrow(res))) {
        mask[res[i, 1]:res[i, 2]] <- 1L
      }
    }
  }
  
  sum(mask)
}


## 2. Input validation (SAFE MODE FOR SOURCE) ----------------------------------

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


## 3–5. Execute only if inputs are available -----------------------------------

if (isTRUE(coverage_available)) {
  
  ## 3. Protein sequence coverage computation ----------------------------------
  
  coverage_results <- fasta_ref %>% 
    dplyr::left_join(peptides_clean, by = "Protein.Group") %>% 
    dplyr::mutate(
      Aminoacids_Covered = purrr::map2_int(
        Sequence,
        Peptide_List,
        calculate_coverage
      ),
      Coverage_Pct = dplyr::if_else(
        Total_Len > 0,
        (Aminoacids_Covered / Total_Len) * 100,
        0
      ),
      Peptides = purrr::map_chr(
        Peptide_List,
        ~ if (is.null(.x) || length(.x) == 0) {
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
  
  if (nrow(coverage_results) == 0) {
    stop("Protein coverage table is empty — computation failed")
  }
  
  
  ## 4. CSV export with peptide evidence ---------------------------------------
  
  coverage_results_detected <- coverage_results %>% 
    dplyr::filter(Aminoacids_Covered > 0)
  
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
  
  
  ## 5. Top 100 summary object for QMD rendering -------------------------------
  
  coverage_results_top100 <- coverage_results %>% 
    dplyr::select(-Peptides) %>% 
    dplyr::arrange(dplyr::desc(Coverage_Pct)) %>% 
    dplyr::slice_head(n = 100)
  
  
  ## 6. Global protein sequence coverage distribution --------------------------
  
  coverage_results_graph <- coverage_results %>% 
    dplyr::filter(Coverage_Pct > 0)
  
  p_coverage <- ggplot2::ggplot(
    coverage_results_graph,
    ggplot2::aes(x = Coverage_Pct)
  ) +
    ggplot2::geom_density(
      fill  = "#56B4E9",
      color = "black",
      alpha = 0.8
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 100, 10),
      limits = c(-1, 101)
    ) +
    ggplot2::labs(
      title = "Proteome Sequence Coverage Distribution",
      x     = "Sequence coverage (%)",
      y     = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA),
      plot.title       = ggplot2::element_text(face = "bold"),
      aspect.ratio     = 0.8
    )
  
  print(p_coverage)
  
} else {
  message("[COVERAGE] Section skipped due to missing inputs")
}