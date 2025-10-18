# Modularity Analysis Method (Phase 04) — Updated

This document describes the modularity analysis pipeline and validation procedures used for bipartite siderophore networks in Phase 04. It includes: a high-level summary of the main analysis script (`03a_modularity_analysis_sk.py`), descriptions of two debug null-test scripts (`03g_edge_rewiring_null_test.py` and `03g_edge_rewiring_null_test_production.py`), orientation conventions for each layer, the null-generation methods used, BRIM optimizer settings and rationale, expected outputs, and a validation checklist with recommended parameters.

Table of contents
- Purpose & scope
- Network conventions & orientation
- Summary: main analysis script (`03a_modularity_analysis_sk.py`)
- Summary: debug null-test scripts
  - `03g_edge_rewiring_null_test.py` (utilization)
  - `03g_edge_rewiring_null_test_production.py` (production)
- Null model methods (swap-based, configuration/stub-matching)
- BRIM optimization: settings and recommendations
- Outputs produced by the pipeline
- Validation guidance and diagnostics
- Reproducibility, seeds, and logging
- Recommended default parameters
- Quick workflow and decision rules
- Contact / provenance

---

Purpose & scope
- Quantify bipartite modular structure separately in the `production` and `utilization` layers of FG–PYO networks.
- Provide robust statistical testing using degree-preserving null models.
- Save compact artifacts sufficient for downstream analyses and reproducibility checks (representative partitions, replicate Q values, null summaries).

Network conventions & orientation
- Canonical orientation used throughout Phase 04:
  - Rows = Functional groups (FG_*) (producers/agents)
  - Columns = Pyoverdines (PYO_*)
- Layer-specific orientation in practice:
  - `production` layer: use original edges as `source -> target` → rows = `source` (FG), cols = `target` (PYO).
  - `utilization` layer: depending on raw data representation, scripts may transpose (flip) the input so rows = FG and cols = PYO. The debug/util scripts explicitly handle this inversion where necessary.
- Always remove zero-degree rows/columns before analysis and log the counts removed.

---

Summary: main analysis script (`03a_modularity_analysis_sk.py`)
Purpose
- The central script that runs the complete modularity analysis for both layers, selects representative partitions, and produces compact outputs and logs.

Key capabilities
- Builds binary incidence matrices from input edge lists (deduplicated).
- Implements `simple_brim()` and `run_brim()` to perform BRIM-like optimization and to compute `Barber` bipartite modularity (`Q`).
- Runs `R` replicates per layer (seeded) and:
  - Records per-replicate `Q`s and seeds.
  - Computes replicate stability metrics (mean, sd, CV).
  - Computes pairwise partition similarity (e.g., ARI or VI), and selects a medoid partition as representative.
  - Records best-Q replicate for transparency.
- Generates null models via `_bipartite_config_model()` (stub-matching implementation) for statistical testing and compares observed `Q` to null `Q` distribution.
- Produces small output artifacts:
  - `results/phase_04/modularity/module_assignments_<layer>_sk.csv`
  - `results/phase_04/modularity/replicate_Q_<layer>_sk.csv`
  - Null replicate file summary and aggregated null summary CSVs
- Contains safety/diagnostic checks (degree sums, zero-degree pruning, input file checks).

Notes on implementation
- The script exposes deterministic seed derivation from a `master_seed` to obtain reproducible per-replicate seeds.
- The `_bipartite_config_model()` routine is a fast stub-matching generator used to cross-check swap-based nulls in debug scripts.

---

Summary: debug null-test scripts

Purpose
- Provide lightweight, script-level cross-validation of the null testing implemented in the main script.
- Compare multiple null generation strategies and test whether `Q` distribution behavior is layer-specific.

1) `03g_edge_rewiring_null_test.py` (utilization debug)
- Designed for the `utilization` layer; it flips input orientation so rows = FG, cols = PYO (matching main-script convention).
- Implements a swap-rewiring null: Maslov–Sneppen pairwise edge swaps that preserve row/column degree sequences.
- Improvements over naive swap implementations:
  - `n_swaps` is interpreted as the number of *successful* swaps (not multiplied by |E| by default).
  - `swaps_per_edge` option available (toggles target swaps relative to edge count).
  - `edge_rewiring_model(..., return_stats=True)` returns actual realized swaps for diagnostics.
  - BRIM (`simple_brim`) includes:
    - Multi-restart support (`n_restarts`).
    - Minimum module enforcement (`min_modules`) to avoid degeneracy (single-module collapse).
    - Convergence checks and deterministic seeding per replicate.
  - Optional per-replicate logging of swaps and `Q`.
- Use: run for varying `n_swaps` or `swaps_per_edge` to examine mixing/stabilization of null `Q` distribution.

2) `03g_edge_rewiring_null_test_production.py` (production debug)
- Same logic as utilization debug but keeps production orientation (rows = `source` FG, cols = `target` PYO); do *not* flip orientation for production.
- Provides identical configuration options (absolute `n_swaps`, `swaps_per_edge`, `n_restarts`, `min_modules`, `log_replicate_stats`).
- Purpose: cross-validate main-script nulls and to show layer-dependent differences (e.g., production null `mean` can be > 0 and stabilize, utilization null mean may be ≈ 0 if the layer is dense and degree-driven).

Why both scripts exist
- Orientation matters and must be handled carefully — the utilization script flips to match the canonical FG rows × PYO cols arrangement, production script does not.
- Having both scripts allows side-by-side comparisons and checks that null behavior is not an artifact of orientation or implementation details.

---

Null model methods (summary & when to use)
Two complementary degree-preserving null approaches are used and cross-compared:

1) Swap-based rewiring (Maslov–Sneppen)
- Iteratively pick two edges (r1,c1) and (r2,c2) and swap to (r1,c2),(r2,c1) if the new edges do not already exist.
- Preserves exact row/column degree sequences.
- Key parameters:
  - `n_swaps` (target successful swaps) or `swaps_per_edge` (target = swaps_per_edge × E).
  - `max_attempts` per replicate to prevent infinite loops (set heuristics like `max(10 * n_swaps, 1000)`).
- Diagnostics:
  - Track realized number of successful swaps.
  - Verify degree preservation for a sample of nulls.
  - Check Jaccard distance or edge-set diversity across null replicates.

2) Configuration / stub-matching (exact degree-preserving generator)
- Build stubs per row/col according to degrees and randomly match row-stubs to column-stubs (permuting column stubs).
- Fast, direct and often used to validate swap-based nulls.
- Drawback: naive stub matching can generate duplicate edges; scripts usually rebuild unique edges and eliminate duplicates—verify degree preservation afterwards.

Guidance
- Use both swap-based and configuration-based nulls for cross-validation. If they agree (same null mean & sd for `Q`), it increases confidence that null inference is robust.
- For swap-based nulls, choose `swaps_per_edge` large enough for mixing (starting point: 10–20; increase if null stats change with more swaps).
- For configuration stub-matching, ensure degree sums match and check for duplicate removal side-effects.

---

BRIM optimization: settings and recommendations
Implementation notes
- `simple_brim()` used across scripts is a BRIM-like alternating label propagation tailored for binary bipartite incidence matrices.
- It computes Barber modularity `Qb` after labels converge.

Important parameters
- `n_iter`: maximum iterations per BRIM run (typical default: 20k in production scripts when maximizing thoroughly).
- `n_restarts`: number of randomized initializations for BRIM; choose 20–100 depending on runtime. Multi-starts mitigate local optima.
- `min_modules`: minimum number of valid (non-empty, both-row-and-col containing) modules required for a restart to be considered valid. Using `min_modules=2` prevents trivial single-cluster solutions from contaminating null variance.
- Tie-breaking and stochasticity:
  - Include randomness in tie-breaking during label assignment to avoid deterministic convergence that may bias `sd_null` downwards.
  - Seeding: pass explicit seeds for reproducibility but decouple seeds for rewiring and BRIM.

Recommendations
- For observed data: run `n_restarts >= 20` (increase to 50–100 if partitions are unstable).
- For nulls: run `n_restarts >= 10` or more to get a robust maximized Q per null matrix.
- If BRIM collapses to one module frequently on randomized matrices, enforce `min_modules=2` during null runs to get meaningful variance estimates for `Q`.

Interpreting BRIM results
- `Q_obs` alone is not a universal signal of structure; use null comparisons (z-score and empirical p-value) to quantify significance.
- Under a well-mixed degree-preserving null, the expected `mean(Q_null)` often approaches 0; however, when optimizing over partitions, the maximum modularity found on random graphs can still be positive (positive null means are common). Use the observed `Q_obs` relative to `sd_null` and empirical p-values for inference.

---

Outputs produced by the pipeline
Per layer outputs (compact & reproducible)
- `results/phase_04/modularity/module_assignments_<layer>_sk.csv`  
  - Columns: `node_id`, `node_type` (`FG` or `PYOV`), `module_id`
- `results/phase_04/modularity/replicate_Q_<layer>_sk.csv`  
  - Columns: `replicate`, `Q`, `random_state`, `brim_iterations`
- Null outputs:
  - `scripts/phase_04_topology/debug/logs/null_replicates_<layer>_sk.csv` (debug scripts) with per-null `Q` and seeds
  - `scripts/phase_04_topology/debug/logs/null_summary_<layer>_sk.csv` with aggregated null summary (mean, sd, z-score, p-value, n_valid_nulls)
- Log files:
  - `docs/phase_04/logs/step03_modularity_summary.txt` — consolidated summary of both layers and diagnostic flags
  - `docs/phase_04/logs/step03_session_info.txt` — environment/session info and parameter provenance

Minimal visualization (optional)
- `figures/network_topology/modularity_<layer>_composition.pdf` — module composition and sizes
- `figures/network_topology/modularity_<layer>_incidence_heatmap.pdf` — ordered incidence heatmap (large matrices may skip heatmap)

---

Validation guidance and diagnostics (practical checklist)
1. Degree preservation check
   - For the first few nulls (both swap-based and stub-matching), confirm that the row/col degree sequences match the original. If not, investigate duplicate-edge removal or bug in stub matching.

2. Null mixing / stabilization
   - Run nulls at increasing `swaps_per_edge` (e.g., 5, 10, 20, 50) and plot `mean(Q_null)` and `sd(Q_null)`. When statistics plateau, mixing is adequate.
   - For stub-matching, compare with swap-based nulls. Close agreement indicates both methods are sampling similar ensembles.

3. BRIM optimizer robustness
   - Increase `n_restarts` and confirm `Q_obs` stabilizes and representative partition (medoid) is consistent.
   - Compute replicate CV = `sd(Q)/mean(Q)`. If `CV > 0.10`, increase `n_restarts` or examine high-degree nodes.

4. Partition stability
   - Compute pairwise ARI (or VI) across replicates. Flag if mean ARI < 0.60.
   - Use medoid partition (minimum mean VI) as representative; save best-Q for transparency.

5. Null inference and interpretation
   - Compute empirical p-value using +1 correction:
     - `p = (count(null_Q >= Q_obs) + 1) / (N + 1)`
   - Compute SES (z-score) if `sd_null` is well-estimated:
     - `z = (Q_obs - mean(Q_null)) / sd(Q_null)`
   - Interpret results in context:
     - `mean(Q_null)` ≈ 0 is expected for a well-mixed degree-preserving null, but `sd(Q_null)` drives significance.
     - Non-zero `mean(Q_null)` (stabilized) can occur when optimizer finds positive structure in random graphs (depends on size/density/heterogeneity). This is expected behavior and must be used when computing z/p.

6. Layer comparison sanity
   - Compare `density`, degree variances, and extreme degrees across layers — dense or skewed degree distributions can explain differences in null behavior (e.g., utilization near-zero vs production positive null mean).

7. Additional diagnostics (optional)
   - Jaccard distances between null edge sets to assess diversity of null ensemble.
   - Report number of valid nulls (non-`NA` `Q` values).
   - Report realized successful swaps per null replicate when `log_replicate_stats=True`.

---

Reproducibility, seeds, and logging
- Master seed is recorded in `docs/phase_04/parameters_manifest.json`.
- Derive per-layer and per-replicate seeds deterministically (e.g., `master_seed + layer_idx * 10^7 + rep_idx`).
- Log seeds in `replicate_Q_<layer>.csv` and in null replicate files.
- Capture software environment (R/Python/NumPy/Scipy/SciPy versions) in `docs/phase_04/logs/step03_session_info.txt`.

---

Recommended default parameters (pilot)
- BRIM:
  - `brim_iterations = 20000` (max per run)
  - `n_restarts = 20` (pilot); increase to `50–100` if instability observed
  - `min_modules = 2` for null runs to avoid degenerate single-cluster outcomes
- Nulls:
  - `n_null_replicates = 100` (pilot)
  - `swaps_per_edge = 10` (start) or `n_swaps = 100` absolute; increase if stabilization not reached
  - `max_attempts` heuristic: `max(10 * n_swaps, 1000)`
- Diagnostics:
  - `R` replicates per layer: pilot `R = 50`; final `R = 100` if needed for stability checks
  - Use +1 correction for empirical p-values

---

Quick recommended workflow
1. Run main modularity analysis (`03a_modularity_analysis_sk.py`) for both layers with pilot replicates.
2. Inspect `docs/phase_04/logs/step03_modularity_summary.txt` for stability flags (CV, ARI).
3. If instability is detected or if null results are surprising, run debug null-test scripts:
   - `03g_edge_rewiring_null_test.py` for utilization (note: flips orientation internally).
   - `03g_edge_rewiring_null_test_production.py` for production.
4. Validate null statistics across both swap-based and configuration-based nulls. If they agree, accept null summaries. If they disagree, increase swap counts and re-evaluate.
5. Save medoid representative partition and best-Q replicate; report both in logs and export minimal artifacts.

---

Contact & provenance
- Every summary must include: timestamp, `master_seed` used, script version/path, input file paths (and modification timestamps/hashes), and parameter values used for the run.
- Store this provenance information in `docs/phase_04/logs/step03_session_info.txt` and include links to `results/phase_04/modularity` artifacts.

---

Appendix: interpretation notes (short)
- A near-zero `mean(Q_null)` is not an indication of an invalid null — it often indicates that, under degree-preserving constraints, there is no systematic modular structure beyond degree effects. The importance is the relative difference between `Q_obs` and null distribution spread (`sd_null`) and the empirical p-value.
- Conversely, a consistently positive `mean(Q_null)` after optimization is common and expected; this reflects the optimizer's ability to find positive deviations even in random graphs — use empirical testing with many nulls to assess significance.

--- 

End of document.