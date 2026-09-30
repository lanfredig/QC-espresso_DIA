# <img src="QC-espresso.ico" width="40" height="40" alt="QC-espresso logo"> QC-espresso

QC-espresso generates an HTML quality-control report and exportable summary tables for data-independent acquisition (DIA) proteomics experiments processed with DIA-NN.

The workflow is written in R, runs through RStudio or a regular R session, and uses the DIA-NN `report.parquet` file together with the same FASTA database used for data processing.

## Requirements

- **R 4.6.0**
- **RStudio** recommended for interactive use
- **Quarto** for HTML report generation
- R package dependencies recorded in `renv.lock`

QC-espresso automatically searches for Quarto in the system `PATH`, common standalone installation directories, recent RStudio installations, and a portable runtime located next to the repository.

### Input files

Two input files are required:

1. **DIA-NN report** in Parquet format, usually named `report.parquet`.
2. **FASTA database** used for the corresponding DIA-NN analysis.

## Installation

### 1. Clone or download the repository

```bash
git clone https://github.com/lanfredig/QC-espresso_DIA.git
cd QC-espresso_DIA
```

To use the validated technical branch:

```bash
git checkout jpr-technical-note
```

### 2. Open the RStudio project

Open:

```text
QC-espresso_DIA.Rproj
```

### 3. Restore the R environment

Run once in the R console:

```r
install.packages("renv")
renv::restore()
```

The lockfile records the validated R package environment. Quarto must also be available to generate the HTML report.

## Usage

### Recommended use in RStudio

Open the project, then run the following code in the R console or in a new R script:

```r
source("run_qc.R")

result <- run_qc(
  report_path = "path/to/report.parquet",
  fasta_path = "path/to/reference.fasta",
  output_root = "output",
  open_report = TRUE
)
```

Replace the two input paths with the locations of the DIA-NN report and FASTA database.

The analysis will:

1. validate the input files and required DIA-NN columns;
2. create a new timestamped output directory;
3. compute the quality-control metrics;
4. render the scientific HTML report;
5. export the coverage and Figures of Merit tables;
6. open the report in the default browser when `open_report = TRUE` and the session is interactive.

### Run the included example

```r
source("run_qc.R")

result <- run_qc(
  report_path = "example-data/dataset_4_contaminans.parquet",
  fasta_path = "example-data/uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta",
  output_root = "output",
  open_report = TRUE
)
```

## Output files

Each analysis is written to a separate directory, so previous results are not overwritten:

```text
output/
└── QC_espresso_DIA_YYYYMMDD_HHMMSS/
    ├── qc_espresso_scientific.html
    ├── figures_of_merit_YYYYMMDD_HHMMSS.csv
    └── protein_coverage_with_peptides_YYYYMMDD_HHMMSS.csv
```

### HTML report

The report includes:

- identification depth at precursor, peptide, and protein-group levels;
- relative abundance and contaminant contribution;
- protein sequence coverage;
- peptide charge, length, and missed-cleavage metrics;
- chromatographic performance and internal trypsin standards;
- precursor mass accuracy;
- DIA-NN quantification-quality metrics;
- raw precursor and protein-level LFQ intensity distributions;
- precursor-level data completeness;
- run-level Figures of Merit;
- ion-mobility metrics when valid ion-mobility data are available;
- reproducibility and session information.

### CSV exports

- `figures_of_merit_*.csv` contains the consolidated run-level QC metrics.
- `protein_coverage_with_peptides_*.csv` contains protein sequence coverage and peptide evidence.

## Important analytical behavior

- The original DIA-NN `Protein.Group` remains the quantitative unit. Accessions in multiprotein groups are expanded only for annotation, without duplicating quantitative signal.
- Missed tryptic cleavages are counted using the corrected terminal K/R rule.
- Protein-level intensity plots use one positive `PG.MaxLFQ` value per `Run + Protein.Group`. Multiple distinct values for the same pair are reported as an input inconsistency.
- Protein coverage is calculated independently for each compatible FASTA sequence. Shared peptide evidence may contribute to more than one sequence, but intensity is not duplicated.
- Optional sections are skipped with informative messages when the required measurements are unavailable.
- The mass-error distribution retains an intentional visual range of -10 to +10 ppm; this display range does not remove observations from the calculations.

## Reproducibility and validation

The repository includes automated tests for input validation, multiprotein annotation, missed cleavages, optional sections, completeness, coverage, protein-level intensity, integrated scientific execution, and end-to-end Quarto rendering.

To run the consolidated audit from the project root:

```r
system2(
  file.path(R.home("bin"), "Rscript"),
  c("--no-init-file", "tests/run_audit.R")
)
```

The validated audit result for the technical branch is:

```text
PASS: 23
FAIL: 0
KNOWN ISSUE: 0
NOT IMPLEMENTED: 0
FINAL STATUS: PASS
```

The audit is intended for development, review, and release validation. Routine users only need to call `run_qc()`.

## Troubleshooting

### Quarto is not found

Install Quarto or use a recent RStudio installation that includes it. The runner also detects a portable Quarto installation placed in a sibling `runtime/quarto` directory.

Check a standalone installation with:

```bash
quarto --version
```

### A required DIA-NN column is missing

QC-espresso validates the report schema before rendering. Confirm that the selected file is the DIA-NN report in Parquet format and that it contains the required identification, quantification, chromatography, and quality columns.

### The report is not generated

Read the error printed by `run_qc()`. If Quarto starts but rendering fails, the timestamped execution directory is preserved for inspection.

### A section does not appear

Some sections are conditional. For example, ion-mobility plots are omitted when the input does not contain valid nonzero ion-mobility measurements.

## Repository structure

```text
QC-espresso_DIA/
├── R/                              # Scientific and input-processing modules
├── report/                         # Quarto scientific report
├── tests/                          # Automated validation and audit suites
├── example-data/                   # Validated example inputs
├── run_qc.R                        # User-facing scientific runner
├── renv.lock                       # Reproducible R package environment
├── VERSION                         # QC-espresso version
└── QC-espresso_DIA.Rproj           # RStudio project
```

## Development status

The `jpr-technical-note` branch contains the validated scientific refactoring and end-to-end rendering workflow. Further work should focus on validation from a clean clone, formal regression of exported files, reviewer-oriented documentation, and release preparation.

## Feedback

Issues and suggestions can be submitted through the repository issue tracker, particularly for:

- reproducibility on additional operating systems;
- DIA-NN schema variants;
- clearer diagnostics for invalid inputs;
- additional QC metrics or export formats.
