Proteomics Data QC – DIA
Launcher Script (Windows Batch)
Version: v0.2.7
Author: Guilherme Lanfredi

============================================================
OVERVIEW
============================================================

This Windows batch script provides a portable launcher for executing the
DIA proteomics Quality Control (QC) report implemented as a Quarto (.qmd)
document.

The launcher is responsible for preparing the execution environment,
validating required inputs, rendering the QC report via Quarto, and
managing the generated HTML output. It is designed for routine use by
analysts and users who may not be familiar with command-line or Quarto
workflows.


============================================================
DEVELOPMENT STATUS AND FEEDBACK
============================================================

This launcher is under active development.

The current version is fully functional and intended for routine QC
report execution. Recent updates focused on improving robustness and
determinism across execution contexts, and on supporting a portable
runtime model for R and Quarto, without modifying the analytical
behavior of the QC report.

Suggestions and feedback are welcome, particularly regarding:
- Improved robustness across different Windows environments
- Additional pre-flight validation checks
- Clearer diagnostics for failure scenarios
- Support for alternative execution modes


============================================================
DESIGN PHILOSOPHY
============================================================

The launcher follows these design principles:

- Portability: no hard-coded installation paths
- Portable-first execution: bundled R and Quarto are used when available
- Safety: required inputs are validated before execution
- Transparency: resolved paths for Quarto and R are printed to the console
- Determinism: the QC script is always executed from a known directory
- Separation of concerns: execution logic is kept independent of QC logic

The launcher does not modify or control the analytical behavior of the QC
script itself.


============================================================
DIRECTORY ASSUMPTIONS
============================================================

The launcher assumes the following directory structure relative to the
location of the .bat file:

<ProjectRoot>/
│
├── Rscripts/
│   └── DIA_QC.qmd
│
├── DIANN-out/
│   └── report.parquet
│
├── FASTA/
│   └── reference FASTA files (.fa, .fasta, .faa)
│
├── Results/
│   └── generated reports (created automatically if missing)
│
├── runtime/
│   ├── R/             (portable R distribution)
│   ├── quarto/        (portable Quarto distribution)
│   └── .Renviron      (project-local R environment configuration)

The directory containing the launcher script is always treated as the
project root.


============================================================
INPUT VALIDATION
============================================================

Before rendering, the launcher performs explicit validation of all
required inputs:

- Presence of the Quarto QC script (DIA_QC.qmd)
- Presence of the DIA-NN report file (report.parquet)
- Presence of a FASTA directory
- Presence of at least one FASTA file with a supported extension

If any required input is missing, execution stops immediately with a
clear error message. Rendering is not attempted unless all inputs are
successfully validated.


============================================================
EXECUTION FLOW
============================================================

The launcher follows a fixed, sequential execution flow:

1) Environment initialization
   - Disable command echo
   - Enable delayed variable expansion
   - Display a versioned startup banner

2) Path resolution
   - Determine project root from the launcher location
   - Define paths for:
     - Rscripts directory
     - Results directory
     - QC script
     - Input data folders

3) Input validation
   - Verify presence of required input files and directories
   - Abort execution if validation fails

4) Tool discovery (portable-first)
   - Prefer bundled Quarto and R from the runtime/ directory
   - Fallback to PATH, registry, RStudio, and standard install locations
   - Abort execution if required tools are not found

5) Environment configuration
   - Force Quarto to use the resolved Rscript executable
   - Activate a project-local .Renviron file
   - Use a portable R library location defined via R_HOME

6) Report rendering
   - Change working directory to Rscripts/
   - Invoke `quarto render DIA_QC.qmd --to html`
   - Capture and evaluate render exit status

7) Output management
   - Generate a timestamp
   - Rename and move the rendered HTML file to Results/
   - Prevent overwriting of previous reports
   - Automatically open the final report

8) Termination
   - Display completion banner
   - Pause to allow inspection of messages


============================================================
WORKING DIRECTORY AND PATH RESOLUTION MODEL
============================================================

The launcher enforces execution from the Rscripts/ directory before
rendering the Quarto document.

Within the QC script, project paths are resolved relative to the Quarto
document location during rendering. When the script is executed
interactively or outside a Quarto context, paths are resolved relative to
the active working directory.

This dual-resolution model ensures:

- Consistent Quarto rendering regardless of invocation method
- Safe interactive execution without a Quarto context
- No reliance on hard-coded absolute paths
- Compatibility with RStudio, command-line, and automated environments


============================================================
OUTPUT HANDLING
============================================================

During rendering, Quarto produces an HTML file in the Rscripts/ directory.
The launcher subsequently:

- Appends a timestamp to the filename
- Moves the file to the Results/ directory
- Ensures previous reports are preserved
- Opens the final report automatically for inspection

Responsibility for HTML output management is intentionally centralized
in the launcher.


============================================================
ERROR HANDLING
============================================================

The launcher implements explicit failure checks at all critical stages:

- Missing input files or directories
- Quarto not found
- Rscript not found
- Quarto rendering failure

In each failure scenario, execution stops immediately and a clear error
message is printed to the console, followed by a pause for user review.


============================================================
LIMITATIONS
============================================================

- The launcher is designed for Windows environments only.
- R and Quarto are required for execution.
- When a portable runtime is provided in the runtime/ directory, bundled
  versions of R and Quarto are used. Otherwise, the launcher falls back to
  system-installed versions.
- The launcher validates the presence, but not the contents, of input
  files.
- Analytical parameters and QC thresholds are controlled exclusively by
  the QC script.

These limitations are intentional to keep the launcher focused on
execution control and portability.


============================================================
THIRD-PARTY SOFTWARE
============================================================

This package may include unmodified binary distributions of third-party
software for execution convenience, including:

- R (GNU General Public License)
- Quarto (MIT License)

These components are distributed in accordance with their respective
licenses and remain the property of their original authors.


============================================================
INTENDED USE
============================================================

This launcher is intended for:

- Routine execution of DIA QC reports
- Use by analysts without command-line experience
- Ensuring consistent execution context across machines
- Archival of timestamped QC reports

It is not intended to replace direct Quarto or RStudio usage during
development or debugging.

============================================================
END OF FILE
============================================================