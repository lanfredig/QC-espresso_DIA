# <img src="QC-espresso.ico" width="40" height="40" alt="QC-espresso logo"> QC-espresso

QC-espresso generates an HTML quality-control report and exportable summary tables for data-independent acquisition (DIA) proteomics experiments processed with DIA-NN.

The scientific workflow is written in R and can be executed through RStudio or a regular R session.

A portable Windows executable is available for [download here](https://github.com/lanfredig/QC-espresso_DIA/releases/download/v1.3.1/QC-espresso_v1.3.1_Windows_portable.zip) as a secondary option for routine use. Both execution modes use the DIA-NN `report.parquet` file together with the same FASTA database used for data processing and call the same scientific workflow.

## Requirements

### 🔬 Scientific and reproducible execution

- **R 4.6.0**
- **RStudio** recommended for interactive use
- **Quarto** for HTML report generation
- R package dependencies recorded in `renv.lock`
- **Rtools45 on Windows**, using the toolchain compatible with R 4.6.x

> [!IMPORTANT]
> Windows users must install Rtools45 before running `renv::restore()`.
> The lockfile may require packages to be installed from source. Packages
> containing compiled C, C++, or Fortran code, including dependencies such
> as `arrow` and `data.table`, require a working compiler toolchain and build
> utilities such as `make` and `sh`. These tools are not included in a
> standard Windows installation and are provided by Rtools45.
>
> Rtools44 is intended for R 4.4.x and must not be used with the validated
> R 4.6.0 environment.

QC-espresso automatically searches for Quarto in the system `PATH`, common standalone installation directories, recent RStudio installations, and a portable runtime located next to the repository.

### ▶ Portable Windows execution

The portable Windows Release includes:

- R 4.6.0;
- the validated R package library;
- Quarto;
- the QC-espresso scientific scripts and report;
- the `QC-espresso.exe` graphical launcher.

The portable package does not require separate installations of R, RStudio, Quarto, Rtools, or R packages because the validated package library is distributed with the application.

### Input files

Both execution modes require:

1. a **DIA-NN report** in Parquet format, usually named `report.parquet`;
2. the **same FASTA database** used for the corresponding DIA-NN analysis.

## Installation

### 1. Scientific and reproducible installation

Clone the repository:

```bash
git clone https://github.com/lanfredig/QC-espresso_DIA.git
cd QC-espresso_DIA
```

#### Windows: install Rtools45

Download and install Rtools45 from the official CRAN page:

https://cran.r-project.org/bin/windows/Rtools/

Use the default installation location:

```text
C:\rtools45
```

The default installation is recommended because R 4.6.x can detect the toolchain without additional manual configuration.

After installation, restart R or RStudio and verify that the required build tools are available:

```r
Sys.which(c("make", "sh", "gcc", "g++"))
```

Each command should return a valid path. Empty values indicate that the build toolchain is not available to the current R session.

Optionally, if the `pkgbuild` package is available, run:

```r
pkgbuild::has_build_tools(debug = TRUE)
```

The result should be `TRUE` before attempting to restore the environment.

#### Open the RStudio project

Open:

```text
QC-espresso_DIA.Rproj
```

#### Restore the recorded R environment

Run once in the R console:

```r
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::restore(prompt = FALSE)
```

The restoration process may compile packages from source. On Windows, do not proceed if Rtools45 is unavailable or if the build-tool checks above fail.

After restoration, confirm that the compiled packages required by the workflow are available and that `arrow` can read the validated Parquet example:

```r
required_packages <- c("arrow", "data.table")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "The following required packages are unavailable: ",
    paste(missing_packages, collapse = ", ")
  )
}

arrow::read_parquet(
  "example-data/report.parquet",
  as_data_frame = FALSE
)
```

Quarto must also be available to generate the HTML report. Verify the installation from a terminal with:

```bash
quarto --version
```

### 2. Portable Windows installation

1. Download `QC-espresso_v1.3.1_Windows_portable.zip` from the GitHub Releases page.
2. Extract the complete archive to a local directory.
3. Keep the `runtime` directory next to the `QC-espresso_DIA` directory.
4. Do not move individual runtime folders or application files after extraction.

The portable package has the following structure:

```text
QC-espresso_v1.3.1_Windows_portable/
├── START_HERE.txt
├── runtime/
│   ├── R-4.6.0/
│   ├── library/
│   └── quarto/
└── QC-espresso_DIA/
    ├── example-data/
    ├── launcher/
    │   └── QC-espresso.exe
    ├── output/
    ├── R/
    ├── report/
    │   ├── QC-espresso.png
    │   └── qc_espresso_scientific.qmd
    ├── run_qc.R
    └── VERSION
```

Rtools45 is not required for the portable Windows Release because the validated R package library is already included. Rtools45 is required only when reproducing the scientific environment from `renv.lock` and source compilation is necessary.

## Usage

### 1. Scientific and reproducible execution

Open the RStudio project and run the following code in the R console or in an R script:

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

#### Run the included example

The repository includes a validated example dataset:

```r
source("run_qc.R")

result <- run_qc(
  report_path = "example-data/report.parquet",
  fasta_path = "example-data/uniprotkb_proteome_UP000005640_AND_revi_2026_05_17.fasta",
  output_root = "output",
  open_report = TRUE
)
```

The example data are included both in the scientific repository and in the portable Windows Release, allowing users to test the complete workflow before analyzing their own data.

### 2. Portable Windows executable

The executable provides a graphical interface for routine Windows use while calling the same validated scientific runner and Quarto report.

1. Open `QC-espresso_DIA\launcher`.
2. Double-click `QC-espresso.exe`.
3. Select the DIA-NN `.parquet` report when prompted.
4. Select the FASTA database used in the DIA-NN analysis.
5. Follow the progress displayed in the console window.
6. Wait for the HTML report to open automatically in the default browser.

The console displays the major analysis stages:

```text
[Loading libraries]
[Loading inputs]
[Runs removed by global filtering]
[Identification numbers]
[Relative abundance]
[Sequence coverage]
[Peptide metrics]
[Chromatography quality]
[Mass accuracy]
[Quantification quality]
[Intensity and LFQ distribution]
[Data completeness]
[Figures of Merit]
[Ion mobility]
[Session information]
```

The execution is complete only after the launcher displays:

```text
[SUCCESS] Report generated:
```

followed by:

```text
PROCESS COMPLETED SUCCESSFULLY
```

## Output files

Both execution modes call the same scientific workflow and generate the same output structure. Each analysis is written to a separate timestamped directory, so previous results are not overwritten.

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

- `figures_of_merit_*.csv` contains consolidated run-level QC metrics.
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

The validated audit result for version 1.3.1 is:

```text
PASS: 23
FAIL: 0
KNOWN ISSUE: 0
NOT IMPLEMENTED: 0
FINAL STATUS: PASS
```

The audit is intended for development, review, and release validation. Scientific and reproducible execution should use `run_qc()` directly. The portable executable is a secondary interface for routine Windows use.

## Troubleshooting

### Rtools45 or Windows build tools are not detected

Confirm that Rtools45 is installed in `C:\rtools45`, restart R or RStudio, and run:

```r
Sys.which(c("make", "sh", "gcc", "g++"))
```

If any returned value is empty, reinstall Rtools45 from the official CRAN page and keep the default installation path. Do not substitute Rtools44 for the validated R 4.6.0 environment.

### `arrow` or another compiled package fails during `renv::restore()`

Confirm that Rtools45 is available before retrying the restore. Remove only the failed partial package installation if instructed by `renv`, then run:

```r
renv::restore(prompt = FALSE)
```

Review the complete installation log if the failure persists. Compiled dependencies may report the underlying missing compiler, build utility, header, or library explicitly.

### Quarto is not found during scientific execution

Install Quarto or use a recent RStudio installation that includes it. The runner also detects a portable Quarto installation placed in a sibling `runtime/quarto` directory.

Check a standalone installation with:

```bash
quarto --version
```

### The portable executable does not start

Confirm that:

- the complete portable archive was extracted before execution;
- `QC-espresso.exe` remains inside `QC-espresso_DIA\launcher`;
- the `runtime` directory remains next to `QC-espresso_DIA`;
- security software has not quarantined or blocked the executable.

### The portable runtime is not found

Confirm that the extracted package still has this structure:

```text
QC-espresso_v1.3.1_Windows_portable/
├── runtime/
└── QC-espresso_DIA/
```

The `runtime` directory must contain:

```text
runtime/
├── R-4.6.0/
├── library/
└── quarto/
```

### A required DIA-NN column is missing

QC-espresso validates the report before rendering. Confirm that the selected file is the DIA-NN report in Parquet format and contains the required identification, quantification, chromatography, and quality columns.

### The report is not generated

Keep the R or executable console open and read the displayed error. If rendering begins but fails, the timestamped execution directory is preserved for inspection.

### A section does not appear

Some sections are conditional. For example, ion-mobility plots are omitted when the input does not contain valid nonzero ion-mobility measurements.

## Repository structure

```text
QC-espresso_DIA/
├── R/                              # Scientific and input-processing modules
├── report/                         # Quarto scientific report and visual assets
├── tests/                          # Automated validation and audit suites
├── example-data/                   # Validated example inputs
├── run_qc.R                        # Cross-platform scientific runner
├── renv.lock                       # Reproducible R package environment
├── VERSION                         # QC-espresso version
└── QC-espresso_DIA.Rproj           # RStudio project
```

## Distribution and source availability

The R scripts, Quarto report, automated tests, example data, and reproducible environment specification are maintained in the GitHub repository.

The portable Windows Release is intended for routine end-user execution and includes the compiled executable, the required scientific scripts, the Quarto report, and the portable runtime. The executable is an optional Windows interface that calls the same scientific workflow implemented by `run_qc.R`.

## Feedback

Issues and suggestions can be submitted through the repository issue tracker, particularly for:

- reproducibility on additional operating systems;
- DIA-NN schema variants;
- clearer diagnostics for invalid inputs;
- additional QC metrics or export formats.
