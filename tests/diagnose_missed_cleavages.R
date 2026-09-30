###############################################################################
# File: tests/diagnose_missed_cleavages.R
# Purpose: Compare the current and conditional missed-cleavage rules without
# modifying analytical scripts or generated reference outputs.
###############################################################################

resolve_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

root <- resolve_project_root()
setwd(root)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(stringr)
})

report_path <- file.path("example-data", "dataset_4_contaminans.parquet")
if (!file.exists(report_path)) stop("Example Parquet was not found.", call. = FALSE)

report <- arrow::open_dataset(report_path, format = "parquet") |>
  dplyr::filter(Decoy == 0, PG.Q.Value <= 0.01) |>
  dplyr::collect()

peptides <- report |>
  dplyr::filter(!is.na(Stripped.Sequence), nzchar(Stripped.Sequence)) |>
  dplyr::distinct(Run, Stripped.Sequence) |>
  dplyr::mutate(
    Cleavage_Sites = stringr::str_count(Stripped.Sequence, "(?<=[KR])(?!P)"),
    Ends_in_KR = stringr::str_detect(Stripped.Sequence, "[KR]$"),
    Current_Missed_Cleavages = pmax(0L, Cleavage_Sites - 1L),
    Conditional_Missed_Cleavages = pmax(
      0L,
      Cleavage_Sites - as.integer(Ends_in_KR)
    ),
    Count_Changed = Current_Missed_Cleavages != Conditional_Missed_Cleavages,
    Current_Has_Missed_Cleavage = Current_Missed_Cleavages > 0L,
    Conditional_Has_Missed_Cleavage = Conditional_Missed_Cleavages > 0L
  )

per_run <- peptides |>
  dplyr::group_by(Run) |>
  dplyr::summarise(
    Unique_Peptides = dplyr::n(),
    Changed_Observations = sum(Count_Changed),
    Current_Missed_Cleavage_Pct = 100 * mean(Current_Has_Missed_Cleavage),
    Conditional_Missed_Cleavage_Pct = 100 * mean(Conditional_Has_Missed_Cleavage),
    Difference_Percentage_Points = Conditional_Missed_Cleavage_Pct - Current_Missed_Cleavage_Pct,
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(abs(Difference_Percentage_Points)))

sequence_summary <- peptides |>
  dplyr::filter(Count_Changed) |>
  dplyr::distinct(Stripped.Sequence, Ends_in_KR, Current_Missed_Cleavages, Conditional_Missed_Cleavages) |>
  dplyr::arrange(Stripped.Sequence)

output_dir <- file.path("tests", "test-results", "latest")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  per_run,
  file.path(output_dir, "missed-cleavage-run-comparison.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

utils::write.table(
  sequence_summary,
  file.path(output_dir, "missed-cleavage-affected-sequences.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

cat("\n============================================================\n")
cat("MISSED-CLEAVAGE DIAGNOSTIC\n")
cat("============================================================\n")
cat("Unique run-peptide observations:", nrow(peptides), "\n")
cat("Observations with changed count:", sum(peptides$Count_Changed), "\n")
cat("Distinct peptide sequences affected:", nrow(sequence_summary), "\n")
cat("Runs affected:", sum(per_run$Difference_Percentage_Points != 0), "\n")
cat("Maximum current percentage:", max(per_run$Current_Missed_Cleavage_Pct), "\n")
cat("Maximum conditional percentage:", max(per_run$Conditional_Missed_Cleavage_Pct), "\n")
cat("Maximum absolute percentage-point difference:", max(abs(per_run$Difference_Percentage_Points)), "\n")
cat("Current values above 50%:", sum(per_run$Current_Missed_Cleavage_Pct > 50), "\n")
cat("Conditional values above 50%:", sum(per_run$Conditional_Missed_Cleavage_Pct > 50), "\n")
cat("============================================================\n\n")

print(per_run, n = Inf)

cat("\nOutputs:\n")
cat("  tests/test-results/latest/missed-cleavage-run-comparison.tsv\n")
cat("  tests/test-results/latest/missed-cleavage-affected-sequences.tsv\n")
