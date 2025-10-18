#!/usr/bin/env python3
"""
Production layer edge-rewiring null test for bipartite modularity (Barber Q).

This script mirrors the utilization-layer script but keeps the original
row/column orientation for the production layer:
- Rows: source (FG_*)
- Columns: target (PYO_*)

Features:
- Degree-preserving edge rewiring (Maslov–Sneppen swaps)
- Either absolute number of successful swaps (n_swaps) or swaps_per_edge
- BRIM-style label propagation with multiple restarts; enforces minimum modules
- Decoupled RNG seeds per replicate step (rewiring vs. BRIM)
- Optional per-replicate logging of realized swaps and Q

Outputs summary statistics:
- Observed Q
- Null mean/sd, z-score, empirical p-value
- Null Q range
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from numpy.random import default_rng
from typing import Optional
from tqdm import tqdm


def build_biadjacency(edges: pd.DataFrame, src_col: str, tgt_col: str):
    """
    Build a bipartite adjacency (rows=src nodes, cols=tgt nodes) from an edge list.
    Returns (csr_matrix, row_labels, col_labels).
    """
    edges = edges[[src_col, tgt_col]].drop_duplicates().copy()
    rows = sorted(edges[src_col].unique())
    cols = sorted(edges[tgt_col].unique())
    r_index = {v: i for i, v in enumerate(rows)}
    c_index = {v: i for i, v in enumerate(cols)}
    r_idx = edges[src_col].map(r_index).to_numpy()
    c_idx = edges[tgt_col].map(c_index).to_numpy()
    data = np.ones(len(edges), dtype=np.float64)
    biadj = csr_matrix((data, (r_idx, c_idx)), shape=(len(rows), len(cols)))
    return biadj, rows, cols


def edge_rewiring_model(
    biadj: csr_matrix,
    n_swaps: int = 10,
    random_state: Optional[int] = None,
    return_stats: bool = False,
):
    """
    Maslov–Sneppen edge rewiring for bipartite graphs:
    - Preserves row/column degree sequences exactly
    - n_swaps is the target number of successful swaps (absolute count)
    - Prevents creation of duplicate edges; bipartite "self-loop" checks are unnecessary

    If return_stats is True, returns (randomized_matrix, successful_swaps),
    otherwise returns randomized_matrix.
    """
    n_rows, n_cols = biadj.shape
    if biadj.nnz == 0 or n_rows == 0 or n_cols == 0 or n_swaps <= 0:
        return (biadj.copy(), 0) if return_stats else biadj.copy()

    rng = default_rng(random_state)

    # Extract current edge list
    rows, cols = biadj.nonzero()
    edges = list(zip(rows, cols))
    current_edges = set(edges)

    successful_swaps = 0
    # Allow enough attempts to reach the target in sparse or constrained graphs
    max_attempts = max(10 * n_swaps, 1000)

    for _ in range(max_attempts):
        if successful_swaps >= n_swaps or len(edges) < 2:
            break

        # Choose two distinct edges at random
        idx1, idx2 = rng.choice(len(edges), size=2, replace=False)
        r1, c1 = edges[idx1]
        r2, c2 = edges[idx2]

        # Proposed swap: (r1,c1),(r2,c2) -> (r1,c2),(r2,c1)
        new_edge1 = (r1, c2)
        new_edge2 = (r2, c1)

        # Disallow if new edges already exist (keeps matrix simple)
        if new_edge1 in current_edges or new_edge2 in current_edges:
            continue

        # Perform swap
        old_edge1 = edges[idx1]
        old_edge2 = edges[idx2]
        edges[idx1] = new_edge1
        edges[idx2] = new_edge2

        current_edges.remove(old_edge1)
        current_edges.remove(old_edge2)
        current_edges.add(new_edge1)
        current_edges.add(new_edge2)

        successful_swaps += 1

    if len(edges) == 0:
        result = csr_matrix((n_rows, n_cols), dtype=np.float64)
        return (result, successful_swaps) if return_stats else result

    rows, cols = zip(*edges)
    data = np.ones(len(rows), dtype=np.float64)
    result = csr_matrix((data, (rows, cols)), shape=(n_rows, n_cols))
    return (result, successful_swaps) if return_stats else result


def simple_brim(
    biadj: csr_matrix,
    n_iter: int = 100,
    random_state: Optional[int] = None,
    n_restarts: int = 1,
    min_modules: int = 1,
) -> float:
    """
    Simplified BRIM-like co-clustering for bipartite modularity (Barber Q):
    - Alternating majority-label updates between row and column partitions
    - Multiple restarts; returns best Q across restarts
    - Enforces a minimum number of non-empty modules (with both rows and cols)

    Returns: best modularity score (float)
    """
    rng = default_rng(random_state)
    n_rows, n_cols = biadj.shape
    if n_rows == 0 or n_cols == 0:
        return 0.0

    biadj = biadj.tocsr()
    biadj_t = biadj.transpose().tocsr()

    best_q = -np.inf

    for _restart in range(max(1, int(n_restarts))):
        # Initialize number of modules and labels
        min_side = max(1, min(n_rows, n_cols))
        if min_side < 3:
            n_init_modules = min_side
        else:
            n_init_modules = min(15, max(2, min_side // 3))
        n_init_modules = max(1, min(n_init_modules, min_side))

        labels_row = rng.integers(0, n_init_modules, n_rows)
        labels_col = rng.integers(0, n_init_modules, n_cols)

        # Iterate until convergence or max iterations
        for _ in range(max(1, n_iter)):
            old_row = labels_row.copy()
            old_col = labels_col.copy()

            # Update row labels
            for i in range(n_rows):
                start, end = biadj.indptr[i], biadj.indptr[i + 1]
                neighbours = biadj.indices[start:end]
                if neighbours.size == 0:
                    continue
                label_counts = np.bincount(labels_col[neighbours])
                labels_row[i] = int(np.argmax(label_counts))

            # Update column labels
            for j in range(n_cols):
                start, end = biadj_t.indptr[j], biadj_t.indptr[j + 1]
                neighbours = biadj_t.indices[start:end]
                if neighbours.size == 0:
                    continue
                label_counts = np.bincount(labels_row[neighbours])
                labels_col[j] = int(np.argmax(label_counts))

            # Convergence
            if np.array_equal(labels_row, old_row) and np.array_equal(
                labels_col, old_col
            ):
                break

        # Compute Barber Q
        total_weight = float(biadj.sum())
        if total_weight == 0:
            q_val = 0.0
        else:
            degrees_row = np.asarray(biadj.sum(axis=1)).ravel()
            degrees_col = np.asarray(biadj.sum(axis=0)).ravel()

            communities = {}
            for idx, label in enumerate(labels_row):
                if label not in communities:
                    communities[label] = {"rows": [], "cols": []}
                communities[label]["rows"].append(idx)
            for idx, label in enumerate(labels_col):
                if label not in communities:
                    communities[label] = {"rows": [], "cols": []}
                communities[label]["cols"].append(idx)

            valid_modules = 0
            quality = 0.0
            for group in communities.values():
                row_idx = np.asarray(group["rows"], dtype=int)
                col_idx = np.asarray(group["cols"], dtype=int)
                if row_idx.size == 0 or col_idx.size == 0:
                    continue
                valid_modules += 1
                sub = biadj[row_idx][:, col_idx]
                intra = float(sub.sum())
                deg_r = float(degrees_row[row_idx].sum())
                deg_c = float(degrees_col[col_idx].sum())
                quality += intra - (deg_r * deg_c) / total_weight

            if valid_modules < max(1, int(min_modules)):
                q_val = -np.inf
            else:
                q_val = quality / total_weight

        if q_val > best_q:
            best_q = q_val

    return 0.0 if not np.isfinite(best_q) else float(best_q)


def main():
    # Configuration
    edges_file = "data/interim/edges_functional_groups_conservative.csv"
    random_state = 42

    # Rewiring configuration
    n_swaps = 1000  # absolute successful swaps (used if swaps_per_edge is None)
    swaps_per_edge = None  # e.g., 10.0 => target_swaps ≈ 10 * nnz; None to use n_swaps

    # BRIM configuration
    brim_iterations = 20000
    n_restarts = 20  # multi-start BRIM per replicate
    min_modules = 2  # enforce at least 2 non-empty modules

    # Null sampling
    n_null_replicates = 100
    log_replicate_stats = False

    print("=== Production网络边重排 Null Test ===")

    # Load data
    edges = pd.read_csv(edges_file)
    prod_edges = edges.query("edge_type == 'production'")[["source", "target"]]
    if prod_edges.empty:
        raise RuntimeError("No production edges available")

    # For production layer, keep original orientation: rows=source(FG), cols=target(PYO)
    biadj, rows, cols = build_biadjacency(prod_edges, "source", "target")

    # Observed modularity (with multi-start)
    master_rng = default_rng(random_state)
    rs_obs = int(master_rng.integers(0, 2**31 - 1))
    q_obs = simple_brim(
        biadj,
        n_iter=brim_iterations,
        random_state=rs_obs,
        n_restarts=n_restarts,
        min_modules=min_modules,
    )

    # Basic stats
    n_edges = int(biadj.nnz)
    max_edges = biadj.shape[0] * biadj.shape[1]
    density = (n_edges / max_edges) if max_edges > 0 else 0.0

    print(f"网络规模: {biadj.shape[0]} x {biadj.shape[1]}")
    print(f"边数: {n_edges}")
    print(f"网络密度: {density:.4f}")
    print(f"原始模块性 (Q_obs): {q_obs:.6f}")
    print(f"BRIM重启次数: {n_restarts}, 最小模块数: {min_modules}")

    # Null model runs
    null_q_values = []
    target_swaps = (
        int(round(swaps_per_edge * n_edges))
        if swaps_per_edge is not None
        else int(n_swaps)
    )
    print(
        f"\n运行 {n_null_replicates} 个边重排null模型，目标交换数 = {target_swaps}"
        + (f"（约等于每边 {swaps_per_edge} 次）" if swaps_per_edge is not None else "")
        + " ..."
    )

    for rep in tqdm(range(n_null_replicates), desc="Null模型进度"):
        rs_rewire = int(master_rng.integers(0, 2**31 - 1))
        rs_brim = int(master_rng.integers(0, 2**31 - 1))

        if log_replicate_stats:
            randomized, swaps_done = edge_rewiring_model(
                biadj, n_swaps=target_swaps, random_state=rs_rewire, return_stats=True
            )
        else:
            randomized = edge_rewiring_model(
                biadj, n_swaps=target_swaps, random_state=rs_rewire, return_stats=False
            )
            swaps_done = None

        q_null = simple_brim(
            randomized,
            n_iter=brim_iterations,
            random_state=rs_brim,
            n_restarts=n_restarts,
            min_modules=min_modules,
        )
        null_q_values.append(q_null)
        if log_replicate_stats:
            print(f"[replicate {rep + 1}] swaps={swaps_done} Q_null={q_null:.6f}")

    # Statistics
    null_q_values = np.array(null_q_values, dtype=float)
    mean_q = float(np.mean(null_q_values)) if null_q_values.size > 0 else float("nan")
    sd_q = (
        float(np.std(null_q_values, ddof=1)) if null_q_values.size > 1 else float("nan")
    )
    if np.isfinite(sd_q) and sd_q > 0.0:
        z_score = float((q_obs - mean_q) / sd_q)
    else:
        z_score = float("nan")
    p_value = (
        float(np.mean(null_q_values >= q_obs))
        if null_q_values.size > 0
        else float("nan")
    )

    print(f"\n=== 结果 ===")
    print(f"Null模型Q值均值: {mean_q:.6f}")
    print(f"Null模型Q值标准差: {sd_q:.6f}")
    print(f"Z-score: {z_score:.6f}")
    print(f"P值: {p_value:.6f}")
    if null_q_values.size > 0:
        print(
            f"Null模型Q值范围: {float(np.min(null_q_values)):.6f} - {float(np.max(null_q_values)):.6f}"
        )

    # --- write the same printed summary to a compact log file ---
    try:
        log_dir = os.path.join(os.path.dirname(__file__), "logs")
    except NameError:
        log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(
        log_dir, "03g_edge_rewiring_null_test_production_results.txt"
    )

    now = datetime.now().isoformat()
    lines = []
    lines.append(f"Timestamp: {now}")
    lines.append("=== Production 网络边重排 Null Test Results ===")
    lines.append(f"Network size: {biadj.shape[0]} x {biadj.shape[1]}")
    lines.append(f"Edges: {n_edges}")
    lines.append(f"Density: {density:.6f}")
    lines.append(f"Observed Q (Q_obs): {q_obs:.6f}")
    lines.append("")
    lines.append("Null model summary:")
    lines.append(f"  Null Q mean: {mean_q:.6f}")
    lines.append(f"  Null Q sd: {sd_q:.6f}")
    lines.append(f"  Z-score: {z_score:.6f}")
    lines.append(f"  P-value: {p_value:.6f}")
    lines.append(f"  BRIM restarts: {n_restarts}, min_modules: {min_modules}")
    if null_q_values.size > 0:
        lines.append(
            f"  Null Q range: {float(np.min(null_q_values)):.6f} - {float(np.max(null_q_values)):.6f}"
        )
    lines.append("")
    lines.append("Parameters:")
    lines.append(f"  n_null_replicates: {n_null_replicates}")
    lines.append(f"  target_swaps: {target_swaps}")
    lines.append(f"  swaps_per_edge: {swaps_per_edge}")
    lines.append(f"  log_replicate_stats: {log_replicate_stats}")

    # write file
    try:
        with open(log_file, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
    except Exception as exc:
        print(f"[WARNING] Failed to write log file {log_file}: {exc}")


if __name__ == "__main__":
    main()
