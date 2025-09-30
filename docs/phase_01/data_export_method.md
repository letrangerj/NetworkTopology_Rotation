# Data Export Method: HDF5 → RDS Workflow (updated)

Overview
- This document explains the actual, current workflow used in this project to convert MATLAB `.mat` inputs into R-friendly serialized objects (`.rds`) for fast downstream analysis.
- The pipeline is implemented with a mix of MATLAB (for initial `.mat` → HDF5 conversion when needed) and R scripts (for HDF5 → R list + `.rds` export). The R scripts are the canonical, repeatable step used in Phase 1.

Why this workflow
- HDF5 provides a cross-language, inspectable container for numeric matrices and byte-encoded strings produced by MATLAB.
- The R-side workflow decodes byte-encoded string fields (uint8 char arrays), reconstructs nested lists for complex fields (sometimes using `jsonlite`), and then writes R-native objects to `.rds` for fast reload with `readRDS()`.

Where files live
- Raw MATLAB inputs (immutable): `data/raw/` (project convention; do not modify).
- HDF5 intermediate files (if produced): `data/interim/` (e.g., `syn.h5`, `rec.h5`, `Syn_PWA_inGroup.h5`, `Rec_MSA_inGroup.h5`).
- RDS outputs (fast-loading intermediates used by downstream R analysis): `data/interim/rds/`.

RDS files currently produced in this repository (examples created during Phase 1)
- `data/interim/rds/syn.rds`           — synthetase (BGC) data exported as an R list
- `data/interim/rds/rec.rds`           — receptor data exported as an R list
- `data/interim/rds/syn_cluster.rds`   — synthetase clustering (`Syn_PWA_inGroup`) data
- `data/interim/rds/rec_cluster.rds`   — receptor clustering (`Rec_MSA_inGroup`) data
- `data/interim/rds/pairing_result.rds`— pairing / final_pairs data (MAT → R conversion produced an R object saved as .rds)

Canonical R scripts that implement the HDF5→RDS conversion
- `scripts/phase_01_explore_validate/save_syn.R`
  - Reads `syn.h5` with `rhdf5::h5read()`
  - Converts uint8 char matrices into R character vectors
  - Decodes JSON-encoded fields (where present) via `jsonlite::fromJSON()`
  - Assembles a named list `syn` and `saveRDS(syn, ".../syn.rds")`
- `scripts/phase_01_explore_validate/save_rec.R`
  - Reads `rec.h5` and converts string and numeric fields similarly
  - Assembles `rec` list and `saveRDS(rec, ".../rec.rds")`
- `scripts/phase_01_explore_validate/save_mats.R` (a small orchestrator / converter)
  - Handles conversion of MATLAB v7 `pairing_result_v7.mat` (via `R.matlab::readMat()`) when that file is present
  - Loads `Syn_PWA_inGroup.h5` / `Rec_MSA_inGroup.h5` cluster files and saves them as `.rds`
  - Writes cluster R objects to `data/interim/rds/`

Key implementation notes and conventions
- String fields stored in HDF5 are typically uint8 char matrices. The R conversion routines:
  - Convert each numeric byte to raw via `as.raw()` and then `rawToChar()`
  - Reconstruct per-row strings with `apply(..., paste, collapse = "")` and `trimws()` to remove padding spaces
- Complex cell arrays that are not trivially representable in HDF5 are either:
  - Converted in MATLAB to JSON blobs (and then parsed in R with `jsonlite::fromJSON()`), or
  - Saved as MAT blobs and loaded with `R.matlab::readMat()` on the R side if available
- After successful HDF5 parsing, each object is written using `saveRDS(..., file = "<path>")`. Use `readRDS()` to load them back quickly for analysis.

Quick usage (run from project root)
- If you already have HDF5 files in `data/interim/`:
  - Run the R scripts (example): use an R session or `Rscript` to execute the specific converter scripts for the file you need.
  - Example loading in R after conversion:
    - `syn <- readRDS("data/interim/rds/syn.rds")`
    - `rec <- readRDS("data/interim/rds/rec.rds")`
- If you have a MATLAB-only `.mat` (legacy v7) pairing file:
  - The helper R script will attempt `R.matlab::readMat()` and then write an R object to `data/interim/rds/`.

Validation and inspection recommendations
- After loading any `*.rds` object, inspect the top-level structure:
  - `str(my_object, max.level = 2)` to confirm fields and types
  - `head(my_object$strainName)` (for `syn`) or `head(my_object$recname)` (for `rec`) to validate string decoding
- For clustering objects, confirm:
  - `dim(cluster$distm)` matches expectation and `length(cluster$dataname)` aligns with `nrow(cluster$distm)`.

Error handling and reproducibility
- The R scripts use try/catch-style patterns to isolate read/write failures and continue processing other files where possible.
- All conversion scripts explicitly write their output paths; confirm the `data/interim/rds/` directory after running conversion.
- Record the R session info (`sessionInfo()`) in `docs/` if you need to later reproduce conversions under the same package versions.

Notes for future work
- If additional complex nested MATLAB cell arrays are encountered, prefer converting them in MATLAB to JSON before HDF5 export to make the R-side parsing deterministic.
- Consider a lightweight R wrapper script that performs all conversions sequentially and emits a short manifest (`data/interim/rds/manifest.json`) listing created files and their checksums.

If you want, I can:
- Add a small wrapper R script that runs the three converter scripts in order and writes a manifest,
- Or update/standardize the naming of the `.rds` outputs (for example: use consistent suffixes like `*_cluster.rds` / `pairing_result.rds`) and adjust callers downstream to match.