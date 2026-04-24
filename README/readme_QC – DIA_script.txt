Proteomics Data QC – DIA
Script Logic and Workflow Description
Version: v0.2.7
Author: Guilherme Lanfredi

============================================================
OVERVIEW
============================================================

This script generates a comprehensive Quality Control (QC) report for
Data-Independent Acquisition (DIA) proteomics datasets processed with
DIA-NN. The report is implemented as a Quarto (.qmd) document that
combines data processing, visualization, and summary tables into a
single, fully self-contained HTML output.

The script is designed to support portable, reproducible execution
across heterogeneous computing environments while maintaining a clear,
deterministic QC workflow.

This README describes the internal logic and workflow of the QC script,
including data flow, intermediate objects, and the rationale underlying
each QC metric and visualization.

Details related to execution orchestration, environment preparation, and
portable runtime handling are documented separately in the launcher
documentation.


============================================================
DEVELOPMENT STATUS AND FEEDBACK
============================================================

This script is under active development.

The current version is fully functional and suitable for routine DIA QC
assessment. Version v0.2.7 represents a major step toward execution
portability and robustness, alongside refinements to aggregation logic
and visualization determinism.

Feedback, suggestions, and adaptations are welcome, particularly with
regard to:
- Additional QC metrics or diagnostic visualizations
- Alternative representations of existing metrics
- Performance optimizations for large DIA-NN reports
- Improved robustness across different DIA-NN configurations
- Improved clarity of metric definitions and aggregation logic

Users are encouraged to report issues, suggest enhancements, or adapt
the script to their specific workflows as needed.


============================================================
AUTHORSHIP AND DEVELOPMENT NOTE
============================================================

The overall concept, QC strategy, analytical logic, and structural
design of this script are original and were conceived and implemented by
the author.

AI-assisted tools were used during development as programming aids to
support code drafting, refactoring, and documentation, in a role similar
to interactive development environments or code suggestion systems.
All methodological decisions, QC definitions, and final code structure
were defined, reviewed, and validated by the author.

Responsibility for the correctness, interpretation, and intended use of
the script rests fully with the author.


============================================================
WHAT IS NEW IN VERSION v0.2.7
============================================================

Version v0.2.7 introduces a major architectural update focused on
execution portability and context-aware robustness, while preserving the
overall scope and analytical intent of the QC report.

1) Portable and Context-Aware Execution (PRIMARY CHANGE)

The QC script is now explicitly designed to support fully portable
execution across heterogeneous environments, including launcher-based,
command-line, and interactive usage.

Key changes include:

- Robust project root discovery based on structural markers rather than
  hard-coded paths.
- Context-aware path resolution that behaves consistently when executed
  via Quarto, directly from RStudio, or through an external launcher.
- Parameterized input handling, allowing DIA-NN report and FASTA paths to
  be injected externally while maintaining sensible defaults.
- Deterministic use of project-local directories for inputs and outputs,
  independent of the active working directory.
- Explicit recording of execution environment metadata to ensure
  traceability and reproducibility.

These changes enable reliable report generation on different machines
without requiring user-specific configuration, fixed installation paths,
or manual path adjustments.

2) Relative Abundance of Top Genes (UPDATED)

- Relative abundance is computed using explicit gene-level aggregation
  prior to ranking.
- Protein intensities are normalized per sample to total signal.
- Normalized intensities are summed at the gene level.
- Genes are ranked by cumulative normalized abundance within each run.
- The ten most abundant genes per run are retained individually.
- All remaining genes are collapsed into a single "Others" category.
- Relative contributions always sum to 100% per sample.
- Visualization uses stacked bar plots with explicit stacking order,
  placing "Others" on top for improved interpretability.
- The y-axis is fixed from 0 to 100% with 10% increments.

This update ensures that ranking reflects true gene-level dominance and
eliminates ambiguities caused by row-level ranking.

3) Protein Sequence Coverage Assessment (introduced in v0.2.6)

- Identified peptides are mapped to full-length protein sequences from a
  FASTA reference.
- Coverage is computed at the amino-acid level as the percentage of the
  sequence covered by observed peptides.
- For visualization purposes, only the most highly covered proteins are
  displayed in the HTML report.
- The complete coverage dataset, including peptide evidence, is exported
  automatically as a timestamped CSV file during report generation.

4) Trypsin Autolysis Peptide Monitoring (introduced in v0.2.6)

- Known trypsin autolysis peptides are monitored as internal QC
  landmarks.
- Retention time reproducibility and quantitative signal stability are
  evaluated when these peptides are detected.
- These peptides serve as endogenous internal standards for LC stability,
  RT reproducibility, injection consistency, and signal stability.

All analyses involving autolysis peptides are conditionally executed and
skipped gracefully when the peptides are absent.


============================================================
SCRIPT EXECUTION FLOW
============================================================

The script follows a linear, modular QC workflow:

1. Load required R libraries and define global plotting parameters
2. Resolve the project root depending on execution context
3. Load DIA-NN results from a Parquet report
4. Apply global QC filters (decoy removal and protein group FDR)
5. Resolve and load a FASTA reference
6. Construct intermediate data objects reused across QC sections
7. Compute QC metrics by category
8. Generate diagnostic plots and summary tables
9. Export selected results (e.g., protein coverage tables)
10. Record execution environment and session metadata


============================================================
CORE DATA OBJECTS
============================================================

report
- Central data frame loaded from DIA-NN report.parquet
- Filtered to remove decoys and enforce protein group confidence
- Serves as the single source of truth for all QC metrics

fasta_ref
- Reference protein sequences loaded from FASTA
- FASTA file is selected automatically, prioritizing recent files
- Provides sequences and lengths for coverage calculation

peptides_clean
- Maps identified peptide sequences to protein groups
- Handles multi-protein assignments explicitly
- Used for protein sequence coverage estimation


============================================================
QC METRIC SECTIONS (SUMMARY)
============================================================

Identification Metrics
- Counts unique precursors and protein groups per run

Relative Abundance
- Computes per-run relative gene abundance with explicit gene-level
  aggregation

Protein Coverage
- Computes amino-acid-level protein sequence coverage
- Exports full coverage results as CSV

Charge State Distribution
- Evaluates precursor charge-state patterns per run

Peptide Length
- Assesses digestion suitability and peptide size distribution

Identifications vs Retention Time
- Evaluates peptide identification density across the LC gradient

Chromatographic Performance
- Reports FWHM and peak base width statistics

Trypsin Autolysis Peptides
- Monitors internal standard RT and signal stability

Mass Accuracy
- Computes and visualizes mass error distributions

Digestion Efficiency
- Estimates missed tryptic cleavages

Intensity Distribution
- Compares raw and normalized intensity distributions

Data Completeness
- Quantifies run × precursor matrix coverage


============================================================
DESIGN PRINCIPLES
============================================================

- Single source of truth for all QC metrics
- Deterministic, explicit aggregation and ordering
- Clear separation between computation and visualization
- Defensive execution for optional inputs
- Emphasis on portability, reproducibility, and traceability


============================================================
THIRD-PARTY SOFTWARE
============================================================

This script relies on third-party software and R packages, including:

- R (GNU General Public License)
- Quarto (MIT License)

These components are used without modification and are governed by their
respective licenses. The QC methodology and analytical logic implemented
in this script are independent of these tools.


============================================================
INTENDED USE
============================================================

This script is intended for technical QC assessment of DIA proteomics
datasets, focusing on instrument performance, acquisition quality, and
data suitability for downstream analysis.

It is not intended for biological interpretation, hypothesis testing, or
statistical inference.

============================================================
END OF FILE
============================================================