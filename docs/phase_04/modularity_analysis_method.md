# Modularity Analysis Plan for Bipartite Siderophore Networks (Phase 04 — Revised)

## Purpose & scope
This document defines a minimal, biologically-aware pipeline for bipartite modularity analysis of siderophore (pyoverdine) networks in Phase 04. It prioritizes layer-specific modularity (production and utilization) as the primary analyses, treats union-incidence as an optional exploratory analysis, and prescribes layer-appropriate null models and stability checks. The design emphasizes correctness, reproducibility, and compact outputs.

Key goals
- Compute bipartite modularity separately for production and utilization layers (primary).
- Assess replicate stability and statistical significance using degree-preserving nulls per layer.
- Save only the minimal programmatic artifacts required for validation and downstream use.
- Keep runtime and storage modest (pilot replicates/nulls; escalate only if necessary).

---

## Network model and inputs

Network layers (directional)
- Production layer: FG → PYO (agent produces pyoverdine)
- Utilization layer: PYO → FG (agent can utilize a pyoverdine)

Primary data (FG-level)
- `data/interim/adj_production_FG_agentsxpyov_conservative.rds`  (FG rows × PYO cols)
- `data/interim/adj_utilization_FG_pyovxagents_conservative.rds` (PYO rows × FG cols; transpose before use)
- Node reference tables (for labels and provenance):
  - `data/interim/nodes_functional_groups_conservative.rds`
  - `data/interim/nodes_pyoverdines_conservative.rds`

Matrix conventions
- For each layer, use a binary incidence matrix:
  - Production-incidence P = 1{adj_production_FG > 0} (FG × PYO)
  - Utilization-incidence U = 1{t(adj_utilization_FG) > 0} (FG × PYO)
- Preprocessing:
  - Remove zero-degree rows and zero-degree columns (record counts removed and reasons).
  - Keep orientation consistent: rows = FGs, cols = PYOs for both layers.
  - Do not conflate production and utilization during primary analyses.

STR-level
- Same approach applies for strain-level networks (STR × PYO), but run STR analyses only after FG-level pilot due to higher runtime/memory.

---

## Primary modularity definition and optimizer

Bipartite modularity
- Use Barber-style bipartite modularity (Qb) as implemented in bipartite-aware algorithms.
- Optimizer: `bipartite::computeModules()` (or equivalent bipartite-optimized routine). The optimizer operates on incidence matrices and evaluates Qb-like objectives appropriate to bipartite data.

Why layer-specific modularity
- Production vs. utilization encode different biological strategies (synthesis vs. uptake); they must be analyzed separately to preserve ecological interpretation.
- Degree-preserving nulls and interpretation must be defined with respect to the analyzed layer.

---

## Replicates and representative partition selection

Replicate policy
- For each layer run R replicates with different seeds to capture optimization stochasticity.
  - Default pilot: R = 50.
  - If stability is inadequate, increase to R = 100.
- Use random initial conditions between replicates (seeded) to avoid order-dependent bias.

Stability metrics (computed per layer)
- Q statistics: mean(Q), sd(Q), min, max, coefficient of variation CV = sd(Q) / mean(Q).
  - Warning threshold: CV > 0.10 indicates optimizer instability; increase R or inspect partitions.
- Partition similarity: pairwise ARI (or Variation of Information, VI) across replicates.
  - Report mean ARI and IQR.
  - Warning threshold: mean ARI < 0.60 indicates low partition consistency.

Representative partition selection
- Medoid partition: the replicate partition with the smallest mean VI (or largest mean ARI) to all other replicates. Use medoid as the representative for downstream reporting and saving.
- Also record the best-Q replicate (for transparency), but prefer medoid for stability.

Determinism
- Use a master seed from `docs/phase_04/parameters_manifest.json` (default 2025). Derive per-replicate seeds deterministically (e.g., master_seed + replicate_index).
- Do not force deterministic row/column ordering before modularity optimization. Controlled seeded randomness across replicates mitigates ordering bias while allowing reproducibility.

---

## Null model design (layer-specific)

General principle
- Nulls must preserve the biologically relevant degree sequence for the layer under test.

Production nulls
- Preserve production-row-degrees (number of pyovs produced per FG) and production-column-degrees (number of FGs producing each pyov).
- Generate null incidence matrices via a degree-preserving bipartite swap method (e.g., curveball / quasiswap / swap.web). Document method choice and swap counts.

Utilization nulls
- Preserve utilization-row-degrees (number of pyovs utilized per FG) and utilization-column-degrees (number of FGs utilizing each pyov).
- Use the same degree-preserving machinery as for production but applied to utilization-incidence.

Union-incidence nulls
- Not recommended by default: the union degree sequence mixes production and utilization degrees and lacks a clear biological constraint.
- Compute union nulls only upon explicit request and interpret results with caution (document biological assumption being preserved).

Null generation parameters (compact defaults)
- Pilot null count N = 100. Increase only if p-values are borderline or estimates unstable.
- Swaps per null (if using swap-based): heuristic 10 × E (E = number of ones in the matrix). Document the number used.
- For `oecosimu`-style wrappers that generate independent nulls, ensure unique seeds per null to avoid serial correlation.

Null testing procedure
- For each null matrix, compute modularity Q (single run of computeModules() is sufficient for the null).
- Aggregate null Qs to produce null mean, sd, min, max.
- Compute empirical p-value with +1 correction:
  - p = (count(null_Q >= obs_Q) + 1) / (N + 1)  for one-tailed positive test (or two-tailed logic if desired).
- Compute SES = (obs_Q − mean(null_Q)) / sd(null_Q).
- Save only summary statistics (no full null matrices) to keep storage minimal.

---

## Required outputs (minimal but sufficient)

Per-layer outputs (production and utilization)
- Programmatic (small):
  - `results/phase_04/modularity/module_assignments_<layer>.rds`  
    - Representative partition: a table with node_id, node_type (FG or PYOV), module_id, and counts.
  - `results/phase_04/modularity/replicate_Q_<layer>.csv`  
    - Columns: replicate_index, seed, Q_value
- Visual (compact):
  - `figures/network_topology/modularity_<layer>_composition.pdf`  
    - Stacked bars showing module composition (FG vs PYOV) and module sizes.
  - Optional heatmap ordered by module if size permits: `figures/network_topology/modularity_<layer>_incidence_heatmap.pdf`
- Logs:
  - `docs/phase_04/logs/step03_modularity_summary.txt` (consolidated summary for both layers; includes stability metrics)
  - `docs/phase_04/logs/step03_session_info.txt` (R session info and package versions)

Union (optional)
- Text-only summary under `docs/phase_04/logs/03b_union_modularity.txt` (dimensions, fill, Q, #modules, and cautions). No RDS by default.

Null diagnostics (debug)
- Text-only summary under `scripts/phase_04_topology/debug/logs/03c_null_modularity_summary.txt` containing:
  - Observed Q, null mean ± sd, null range, p-value, SES, N_valid_nulls
  - Short observed stability diagnostics (e.g., Q CV, mean ARI over R = 20)
- No storage of null matrices or full replicate partitions.

Rationale for minimal artifacts
- Representative partition RDS and replicate Q vector are the minimal artifacts needed to:
  - Reproduce module assignments
  - Allow downstream validation (taxonomic, metadata joins)
  - Re-run focused follow-up analyses without storing large intermediate datasets
- Avoid saving thousands of null matrices or full replicate partitions to control size.

---

## Scripts and responsibilities

1) Main script (production + utilization; primary)
- Path:
  - `scripts/phase_04_topology/03_modularity_analysis.R` (replaces existing)
- Responsibilities:
  - Load adjacency matrices · validate shapes and names
  - Build P and U incidence matrices (binary)
  - For each layer:
    - Remove zero-degree rows/cols (log counts removed)
    - Run R replicates (default R = 50) of `computeModules()` with seeded randomness
    - Compute stability metrics (Q mean/sd/CV, pairwise ARI distribution)
    - Select medoid representative partition and best-Q replicate
    - Save `module_assignments_<layer>.rds` and `replicate_Q_<layer>.csv`
    - Save the compact module composition figure
  - Write consolidated summary and session info to `docs/phase_04/logs/step03_modularity_summary.txt`

2) Additional (optional exploratory union)
- Path:
  - `scripts/phase_04_topology/additional/03b_union_modularity.R`
- Responsibilities:
  - Build union incidence I = 1{P > 0 OR U > 0}
  - Compute a small number of replicates (R = 10) or a single run for a quick exploratory Q
  - Save a short, cautious text summary at `docs/phase_04/logs/03b_union_modularity.txt`
  - Do not run nulls by default

3) Debug (null tests; layer-specific; compact)
- Path:
  - `scripts/phase_04_topology/debug/03c_null_model_modularity.R`
- Responsibilities:
  - Accept a `layer` argument (production or utilization) or run both serially
  - Generate N degree-preserving nulls (default N = 100) for the given layer
  - For each null, compute Q via `computeModules()` (single run) and collect Qs
  - Report null summary statistics, p-value (+1 correction), SES, and a small observed stability check (R = 20)
  - Save only a text summary at `scripts/phase_04_topology/debug/logs/03c_null_modularity_summary.txt`

---

## Parameters, seeds, and reproducibility

Parameter source
- Primary configuration read from:
  - `docs/phase_04/parameters_manifest.json`
- Required fields (with defaults):
  - `master_seed`: default 2025
  - `modularity_replicates.pilot`: default 50
  - `modularity_replicates.final`: default 100
  - `null_model_N.pilot`: default 100
  - `swaps_per_null`: default expression "10 * E" (document actual numeric choice)

Seed management
- Derive per-task seeds from `master_seed` deterministically (e.g., master_seed + task_index + replicate_index).
- Document all seeds used in replicate_Q files and in null summaries for traceability.

Session info
- Each script must write `docs/phase_04/logs/step03_session_info.txt` (or script-specific session info) containing R version and package versions used.

---

## Quality control checks (per run)

Mandatory checks (fail early and log)
- Input alignment: row/col names match between adjacency matrices and node tables.
- Incidence sanity: sum of ones > 0; no NA entries.
- Degree check: confirm row/col degree distributions before and after zero-degree pruning.

Stability/reporting checks
- CV(Q) reported and flagged if > 0.10.
- mean ARI reported and flagged if < 0.60.
- Module composition: report count of modules that contain both FG and PYOV nodes (bipartite-mixed).
- If instability flagged:
  - Increase replicates to final R (e.g., 100) and re-evaluate.
  - Consider alternative optimization parameters or examine problematic high-degree nodes.

Null-model diagnostics
- For swap-based nulls:
  - Document swaps per null and verify degree sequences preserved (spot-check).
  - Report number of valid nulls (non-NA Qs) and success rate.
- Independence check:
  - For N >= 20, compute difference in mean Q between first and second halves of nulls; warn if large deviation.

---

## Limitations and caveats
- `computeModules()` is stochastic; representative partition selection and stability reporting mitigate but do not eliminate optimizer uncertainty.
- Union-incidence aggregation mixes biologically distinct interactions and should be interpreted cautiously; it is retained here only as an optional exploratory complement.
- Null models preserve layer-specific degrees but cannot account for latent biological covariates (phylogeny, habitat) — such covariate analyses belong in later phases (Phase 06).
- STR-level analyses are computationally heavier; pilot on FG-level first.

---

## Minimal recommended workflow (quick)

1. Run main script:
   - `Rscript scripts/phase_04_topology/03_modularity_analysis.R`
   - Inspect `docs/phase_04/logs/step03_modularity_summary.txt` and composition figures.
2. If partitions are unstable (CV > 0.1 or mean ARI < 0.6), rerun with increased replicates (R = 100).
3. If you need statistical testing, run debug null script for the relevant layer:
   - `Rscript scripts/phase_04_topology/debug/03c_null_model_modularity.R --layer production`
4. Optionally run union exploratory script for contrast:
   - `Rscript scripts/phase_04_topology/additional/03b_union_modularity.R`

---

## Next steps for implementation
- Update `scripts/phase_04_topology/03_modularity_analysis.R` to implement the per-layer replicate and medoid selection scheme and to save the minimal artifacts listed above.
- Create the additional and debug scripts with compact behavior and careful documentation of seeds and null settings.
- Ensure `docs/phase_04/parameters_manifest.json` contains the necessary fields and defaults.

---

## Contact and provenance
- Every produced summary must include: timestamp, master_seed, software version (session info), input file paths and their modification timestamps (or file hashes), and parameter values used. This is required for reproducibility and later inspection.