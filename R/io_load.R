###############################################################################
# File: io_load.R
#
# Purpose:
#   Resolve all input and output paths deterministically, load and filter the
#   DIA-NN report, identify runs completely removed by the global confidence
#   filter, construct peptide evidence per protein group, locate and load the
#   FASTA reference, and perform data-dependent global setup required by
#   downstream QC modules.
#
# Inputs (expected to exist before sourcing):
#   - PROJECT_ROOT : Absolute project root directory
#   - params       : Quarto parameters list
#   - cb_palette   : Global color palette defined in globals.R
#
# Outputs (objects defined):
#   - results_dir                   : Canonical output directory for results
#   - report_path                   : Resolved DIA-NN report file path
#   - report                        : Filtered DIA-NN report data frame
#   - runs_removed_by_global_filter : Runs completely removed by initial filters
#   - peptides_clean                : Peptide lists grouped by protein group
#   - fasta_path                    : Selected FASTA reference file path
#   - fasta_ref                     : Parsed FASTA reference table
#   - file_palette                  : Run-level color palette derived from cb_palette
#
# Required R packages (loaded by QMD):
#   - dplyr
#   - tidyr
#   - arrow
#   - microseq
#   - scales
#   - grDevices
#
# Notes:
#   - No library() or require() calls are made in this script.
#   - No setwd() calls are made in this script.
#   - Designed to be fully portable across local, CI, and containerized
#     execution environments.
#   - The global analytical filter remains unchanged: Decoy == 0 and
#     PG.Q.Value <= 0.01.
#   - Runs completely removed by this filter are recorded so that the QMD can
#     display a conditional notice in the final report.
###############################################################################


## 0. FASTA header parsing helpers ---------------------------------------------

parse_fasta_header <- function(header) {

  header_clean <- trimws(sub("^>", "", header))

  if (is.na(header_clean) || !nzchar(header_clean)) {
    return(
      data.frame(
        Header_Original = NA_character_,
        Protein.Group   = NA_character_,
        Protein.Name    = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }

  header_parts <- strsplit(header_clean, "\\s+", perl = TRUE)[[1]]
  first_token  <- header_parts[1]
  description  <- if (length(header_parts) > 1) {
    paste(header_parts[-1], collapse = " ")
  } else {
    NA_character_
  }

  description_clean <- description

  if (!is.na(description_clean) && nzchar(description_clean)) {
    description_clean <- sub(
      "\\s+(OS|OX|GN|PE|SV)=.*$",
      "",
      description_clean,
      perl = TRUE
    )
    description_clean <- trimws(description_clean)
    if (!nzchar(description_clean)) {
      description_clean <- NA_character_
    }
  }

  accession <- first_token
  entry_name <- NA_character_

  if (grepl("\\|", first_token, perl = TRUE)) {

    pipe_parts <- strsplit(first_token, "\\|", perl = TRUE)[[1]]
    pipe_parts <- pipe_parts[nzchar(pipe_parts)]

    if (length(pipe_parts) >= 2) {

      db_tag <- pipe_parts[1]

      if (db_tag %in% c("sp", "tr", "ref", "gb", "emb", "dbj", "pir", "prf")) {
        accession <- pipe_parts[2]
        if (length(pipe_parts) >= 3) {
          entry_name <- pipe_parts[3]
        }
      } else if (db_tag == "gi" && length(pipe_parts) >= 4) {
        accession <- pipe_parts[4]
        if (length(pipe_parts) >= 5) {
          entry_name <- pipe_parts[5]
        }
      } else {
        accession <- pipe_parts[length(pipe_parts)]
        if (length(pipe_parts) >= 2) {
          entry_name <- pipe_parts[length(pipe_parts) - 1]
        }
      }
    }
  }

  accession <- trimws(accession)
  if (!nzchar(accession)) {
    accession <- NA_character_
  }

  protein_name_candidates <- c(description_clean, entry_name, accession)
  protein_name_candidates <- protein_name_candidates[
    !is.na(protein_name_candidates) & nzchar(protein_name_candidates)
  ]

  protein_name <- if (length(protein_name_candidates) > 0) {
    protein_name_candidates[1]
  } else {
    NA_character_
  }

  data.frame(
    Header_Original = header_clean,
    Protein.Group   = accession,
    Protein.Name    = protein_name,
    stringsAsFactors = FALSE
  )
}


## 1. Canonical output directory -----------------------------------------------

if (!is.null(params$output_dir)) {
  results_dir <- normalizePath(
    params$output_dir,
    winslash = "/",
    mustWork = FALSE
  )
} else {
  results_dir <- file.path(PROJECT_ROOT, "output", "Results")
}

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
message("Results directory: ", results_dir)


## 2. DIA-NN report path resolution --------------------------------------------

if (!is.null(params$report_path)) {

  report_path <- normalizePath(
    params$report_path,
    winslash = "/",
    mustWork = TRUE
  )

} else {

  report_path <- normalizePath(
    file.path(PROJECT_ROOT, "output", "DIANN-out", "report.parquet"),
    winslash = "/",
    mustWork = TRUE
  )
}

if (!file.exists(report_path)) {
  stop("report.parquet not found at: ", report_path)
}

message("Using DIA-NN report: ", report_path)


## 3. Load and filter DIA-NN report --------------------------------------------

report_raw <- arrow::read_parquet(report_path)

required_report_columns <- c(
  "Run",
  "Decoy",
  "PG.Q.Value"
)

missing_report_columns <- setdiff(
  required_report_columns,
  colnames(report_raw)
)

if (length(missing_report_columns) > 0) {
  stop(
    "Required DIA-NN report column(s) not found: ",
    paste(missing_report_columns, collapse = ", ")
  )
}

raw_run_names <- report_raw %>%
  dplyr::filter(
    !is.na(Run),
    nzchar(trimws(as.character(Run)))
  ) %>%
  dplyr::transmute(
    Run = as.character(Run)
  ) %>%
  dplyr::distinct()

report <- report_raw %>%
  dplyr::filter(
    Decoy == 0,
    PG.Q.Value <= 0.01
  )

retained_run_names <- report %>%
  dplyr::filter(
    !is.na(Run),
    nzchar(trimws(as.character(Run)))
  ) %>%
  dplyr::transmute(
    Run = as.character(Run)
  ) %>%
  dplyr::distinct()

runs_removed_by_global_filter <- raw_run_names %>%
  dplyr::anti_join(
    retained_run_names,
    by = "Run"
  ) %>%
  dplyr::arrange(Run)

if (nrow(report) == 0) {
  stop(
    paste0(
      "No entries remained after removal of decoys and application of ",
      "PG.Q.Value <= 0.01."
    )
  )
}

message(
  "DIA-NN report loaded: ",
  nrow(report),
  " retained rows across ",
  dplyr::n_distinct(report$Run),
  " runs"
)

if (nrow(runs_removed_by_global_filter) > 0) {
  message(
    nrow(runs_removed_by_global_filter),
    " run(s) were completely removed by the global filter"
  )
}

rm(
  report_raw,
  raw_run_names,
  retained_run_names,
  required_report_columns,
  missing_report_columns
)


## 4. Peptide evidence per protein group ---------------------------------------

peptides_clean <- report %>%
  dplyr::select(Stripped.Sequence, Protein.Group) %>%
  tidyr::separate_rows(Protein.Group, sep = ";") %>%
  dplyr::mutate(
    Protein.Group = trimws(Protein.Group)
  ) %>%
  dplyr::filter(!is.na(Protein.Group), nzchar(Protein.Group)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(Protein.Group) %>%
  dplyr::summarise(
    Peptide_List = list(Stripped.Sequence),
    .groups = "drop"
  )

message(
  "Peptide list generated for ",
  nrow(peptides_clean),
  " protein groups"
)


## 5. FASTA reference discovery ------------------------------------------------

fasta_candidates <- character(0)

if (!is.null(params$fasta_path)) {

  fasta_param <- normalizePath(
    params$fasta_path,
    winslash = "/",
    mustWork = TRUE
  )

  if (dir.exists(fasta_param)) {
    fasta_candidates <- list.files(
      fasta_param,
      pattern = "\\.(fa|fasta|faa)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  } else {
    fasta_candidates <- fasta_param
  }

  message("FASTA path provided via parameter")

} else {

  fasta_dirs <- c(
    file.path(PROJECT_ROOT, "data", "FASTA"),
    file.path(PROJECT_ROOT, "FASTA"),
    file.path(dirname(report_path), "FASTA")
  )

  for (d in fasta_dirs) {

    d <- normalizePath(d, winslash = "/", mustWork = FALSE)

    if (dir.exists(d)) {
      fasta_candidates <- c(
        fasta_candidates,
        list.files(
          d,
          pattern = "\\.(fa|fasta|faa)$",
          full.names = TRUE,
          ignore.case = TRUE
        )
      )
    }
  }
}

if (length(fasta_candidates) == 0) {
  stop("No FASTA file found")
}


## 6. FASTA selection strategy -------------------------------------------------

non_contaminant <- !grepl(
  "crap|contaminant",
  basename(fasta_candidates),
  ignore.case = TRUE
)

if (any(non_contaminant)) {
  fasta_candidates <- fasta_candidates[non_contaminant]
}

fasta_path <- fasta_candidates[
  which.max(file.info(fasta_candidates)$mtime)
]

message("Using FASTA reference: ", fasta_path)


## 7. Load and parse FASTA reference -------------------------------------------

fasta_raw <- microseq::readFasta(fasta_path)

fasta_header_info <- dplyr::bind_rows(
  lapply(fasta_raw$Header, parse_fasta_header)
)

fasta_ref <- fasta_raw %>%
  dplyr::mutate(
    Header = trimws(sub("^>", "", Header))
  ) %>%
  dplyr::bind_cols(fasta_header_info) %>%
  dplyr::select(
    Header = Header_Original,
    Protein.Group,
    Protein.Name,
    Sequence
  ) %>%
  dplyr::mutate(
    Protein.Group = trimws(Protein.Group),
    Protein.Name  = trimws(Protein.Name),
    Total_Len     = nchar(Sequence)
  ) %>%
  dplyr::filter(
    !is.na(Protein.Group),
    nzchar(Protein.Group),
    !is.na(Sequence),
    nzchar(Sequence)
  )

dup_fasta <- fasta_ref %>%
  dplyr::count(Protein.Group, name = "n_entries") %>%
  dplyr::filter(n_entries > 1)

if (nrow(dup_fasta) > 0) {
  message(
    "Duplicated FASTA accessions detected: ",
    nrow(dup_fasta),
    " accessions. Keeping the first occurrence per Protein.Group."
  )

  fasta_ref <- fasta_ref %>%
    dplyr::distinct(Protein.Group, .keep_all = TRUE)
}

if (nrow(fasta_ref) == 0) {
  stop("FASTA parsing failed: no valid protein entries were recovered")
}

message("FASTA loaded: ", basename(fasta_path))
message("FASTA entries parsed: ", nrow(fasta_ref))


## 8. Data-dependent global setup ----------------------------------------------

file_palette <- grDevices::colorRampPalette(cb_palette)(
  dplyr::n_distinct(report$Run)
)

report <- report %>%
  dplyr::mutate(
    Run = factor(Run, levels = unique(Run))
  )

message("Input loading completed successfully")
