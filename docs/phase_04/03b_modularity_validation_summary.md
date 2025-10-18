## Summary: Modularity Analysis and Null Testing Validation

### What the Main Script Does (03a_modularity_analysis_sk.py)
- **Modularity method**: Barber bipartite modularity (Qb) for a bipartite incidence matrix (rows = FG, cols = PYO). Implemented in `_barber_modularity(biadj, labels_row, labels_col)`; computes intra-module edge weight minus degree-based expectation and normalizes by total weight.
- **Optimizer**: A BRIM-like alternating-label majority vote implemented as `simple_brim` that alternately updates row labels by majority neighbor column labels, then column labels by majority neighbor row labels, up to `brim_iterations` or convergence. In the committed `03a_modularity_analysis_sk.py`, `simple_brim` is single-start and returns `(labels_row, labels_col, quality)`; the pipeline achieves multi-start behavior by running repeated replicates (`n_replicates`) and selecting the best-Q replicate (via `run_brim` and the outer replicate loop). The debug scripts provide an alternative `simple_brim` with internal `n_restarts` and `min_modules`.
- **Stability analysis**: Runs multiple replicates (configured `n_replicates = 100` in the committed script) with seeded randomness. Reports per-replicate Q and basic statistics (mean Q, sd, min, max) in `replicate_Q_<layer>_sk.csv`. **Does not compute ARI or medoid selection**; selects the best-Q replicate partition as the representative.
- **Null model**: Degree-preserving bipartite configuration model (`_bipartite_config_model`): exact row/column degree sequence is preserved via stub-matching (row and column stub lists are created and column stubs shuffled and matched). The script reports degree-preservation checks for the first null replicate.
  - Implementation nuance: stub-matching can pair multiple stubs to the same (row,col) pair, producing integer (>1) entries in the CSR at that position (i.e., parallel counts). Degree sums are preserved; if strictly binary adjacency (0/1) is required, binarization would change degree sums and thus constitute a different null model.
- **Output**: Membership tables (`module_assignments_<layer>_sk.csv`), per-replicate Q CSV (`replicate_Q_<layer>_sk.csv`), null replicate CSV (`null_replicates_<layer>_sk.csv`), null summary CSV (`null_summary_<layer>_sk.csv`), and summary text logs. Uses deterministic per-task seeds derived from a master seed for reproducibility.
- **Orientation handling**:
  - Production: uses original direction (rows = source (FG), cols = target (PYO)).
  - Utilization: flips the direction so rows = FG, cols = PYO before building the bipartite adjacency.

### What Seemed Odd
When running the main script, the null distribution for the utilization layer appeared to trend toward a mean near 0 as swaps/mixing increased, while the production null stabilized at a positive mean. This prompted verification to confirm whether this was an implementation artifact or a real layer difference.

### How We Used the Two Debug Scripts to Verify

#### 1) Debug Setup
- Both debug scripts (`03g_edge_rewiring_null_test.py` for utilization and `03g_edge_rewiring_null_test_production.py` for production) were written to:
  - Re-implement a different null generation method (Maslov–Sneppen edge swaps) with explicit control over the number of successful swaps (`n_swaps` or `swaps_per_edge`).
  - Use a stronger BRIM optimization in the debug code with internal multiple restarts and a minimum-modules constraint (`n_restarts = 20`, `min_modules = 2`) to avoid degenerate single-cluster partitions.
  - Decouple RNG seeds per replicate step (rewiring vs. BRIM).
  - Output compact summaries to `logs/03g_edge_rewiring_null_test_results.txt` (utilization) and `logs/03g_edge_rewiring_null_test_production_results.txt` (production).

#### 2) Key Verifications
- Two null generators were cross-compared:
  - Swap-based edge rewiring (degree-preserving).
  - Configuration/stub-matching generator used in the main script.
- Layer-specific orientation was handled correctly (utilization scripts invert orientation where needed).
- Mixing intensity was varied (`n_swaps` from low to high or `swaps_per_edge`) and stabilized null statistics were observed.

#### 3) Findings
- **Utilization**:
  - Swap-based rewiring and stub-matching nulls agree: null mean trends to ~0 as mixing increases, and sd stabilizes.
  - This is expected for dense or degree-heterogeneous bipartite graphs where optimized Q under the configuration model is small.
  - After increasing BRIM restarts and enforcing `min_modules=2` in debug runs, null variance is well-behaved and empirical p-values are sensible.
  - Parity note: the utilization debug script computed the observed Q (`q_obs`) with a single-start BRIM call while the debug nulls used multi-restarts; for strict parity compute `q_obs` with the same multi-restart settings used for null optimization.
- **Production**:
  - Both null generators agree on a stabilized positive null mean (not zero), which is expected: the maximum achievable Barber Q on random bipartite graphs under degree constraints can be positive, so the null ensemble mean can be > 0.
  - This layer difference appears structural rather than an implementation artifact.

### Conclusion
- The observed differences between utilization and production null distributions are consistent across null generators and are interpretable in terms of network density and degree heterogeneity.
- No implementation bug was found in the configuration-model null or BRIM implementations; the debug scripts validated the main script's behavior when accounting for the noted differences (function signatures, multi-start treatment, stub-matching nuance).
- Recommendation: ensure observed-Q optimization uses the same optimization regime (single-start vs multi-start, and `min_modules` setting) as null replicates for a strict comparison; choose swap-based rewiring if you require strictly binary nulls that preserve degree sequences without stub-duplication.

### Artifacts
- Main analysis: `03a_modularity_analysis_sk.py` (replicate Q CSVs, membership CSVs, summary logs).
- Debug validation:
  - `scripts/phase_04_topology/debug/03g_edge_rewiring_null_test.py` → `logs/03g_edge_rewiring_null_test_results.txt`
  - `scripts/phase_04_topology/debug/03g_edge_rewiring_null_test_production.py` → `logs/03g_edge_rewiring_null_test_production_results.txt`

### Critical Stats Comparison (from debug runs)

| Layer | Q_obs (best) | Null Mean (stub/swaps) | Null SD | Z-score | Empirical p | Network notes |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| Utilization | 0.323641 | 0.001433 | 0.011840 | 27.212581 | 0.000000 | size: 173 x 93; edges: 2238; density: 0.139101; null range: 0.000000 - 0.118386; n_null_replicates: 100 |
| Production  | 0.857273 | 0.832903 | 0.010498 | 2.321522 | 0.010000 | size: 134 x 93; edges: 175; density: 0.014043; null range: 0.802906 - 0.859918; n_null_replicates: 100 |

---