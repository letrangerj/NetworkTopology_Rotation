# Debug diagnostics for Phase 3 Step 1

This directory contains a focused set of diagnostic scripts that inspect and validate the agent-set outputs produced in Phase 2 and Step 1 of Phase 3. The diagnostics are intentionally compact and output only the essential CSV/text artifacts needed to decide how to proceed with pyoverdine pair mapping and adjacency construction.

Purpose
- Verify functional-group (FG) collapse integrity and size distribution.
- Summarize validated pyoverdine component coverage vs. raw production groups.
- Map component pyoverdine IDs to pairing rows (lock-key pairs) and report pair-level activity.
- Provide a concise, actionable summary and a small set of CSVs for downstream review.

Scripts
- `diagnostics_step1_basic.R`
  - Quick checks: top-10 FG sizes, counts of raw vs. validated pyoverdine components, and a short plain-text summary.
  - Produces a minimal set of files for rapid inspection.

- `diagnostics_step1_detailed.R`
  - Builds `component -> pair_id` mapping from `pairing_result`.
  - Produces pair-level activity/coverage summaries and a trimmed deep inspection of the largest FG (`FG_067`) (member list and component/pair frequencies).
  - Produces only the essential CSVs needed to evaluate pair activation and the influence of large FGs.

How to run
From the project root (where `data/` and `scripts/` live) run:
- `Rscript scripts/phase_03_network/debug/diagnostics_step1_basic.R`
- `Rscript scripts/phase_03_network/debug/diagnostics_step1_detailed.R`

Notes:
- `diagnostics_step1_detailed.R` requires the `pairing_result` RDS produced by Phase 1/2 to be present in `data/interim/rds/` (or an expected location). If that file is not present the script will exit with a clear message.
- Both scripts expect the Phase 2/Step1 outputs in `data/interim/`:
  - `nodes_functional_groups_conservative.rds`
  - `nodes_strains_conservative.rds`
  - `strain_repertoires.rds`
  - `receptor_to_pyoverdine.rds`

Essential outputs (written to `scripts/phase_03_network/debug/`)
- `fg_top10.csv` — concise table of the top 10 functional groups by size.
- `pyov_summary.csv` — counts: raw production groups, validated components, pair rows, and active components.
- `diagnostics_step1_basic_summary.txt` — short human-readable summary for quick review.
- `pair_activity_summary.csv` — per-pair component activity (n components / n active components).
- `pair_level_coverage_summary.csv` — per-pair FG & strain coverage summary (how many FGs/strains produce or utilize each pair).
- `active_pair_ids.csv` / `inactive_pair_ids.csv` — lists of pair rows that are active/inactive.
- If present, trimmed FG diagnostics for `FG_067`:
  - `FG_067_member_strains_full.csv`
  - `FG_067_receptor_group_freq.csv`
  - `FG_067_production_pairid_freq.csv`
  - `FG_067_util_pairid_freq.csv`

Quick interpretation guide
- If a single FG (e.g., `FG_067`) contains a very large fraction of strains (>20–30%), expect large influence on modularity / nestedness; plan sensitivity tests (remove or downsample that FG).
- If many validated `pair_id` rows are inactive, re-check pairing mappings and strain repertoires; if all pair rows are active, proceed to build pair-aggregated adjacency matrices.
- Use the `pair_level_coverage_summary.csv` to decide whether to aggregate to pair-level nodes (recommended) and to mark `U_active` vs. `U_all`.

Recommended next steps after diagnostics
1. Review the two summary files (`diagnostics_step1_basic_summary.txt`, `diagnostics_step1_detailed_summary.txt`) and the top-10 FG list.
2. If you confirm the pair-level activity is sufficient, run Step 2 (pyoverdine index mapping) to create the `nodes_pyoverdines_conservative` dictionary and the contiguous `matrix_index` for adjacency building.
3. Plan sensitivity analyses for any very large FG before running large-scale topology computations.

Contact / provenance
- These scripts are deterministic (seeded where sampling is used) and intentionally produce a minimal set of outputs to simplify review and reduce disk noise.
- If you need additional fields or different aggregation levels, edit the scripts in this folder to include the desired exports.
