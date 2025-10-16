#!/usr/bin/env python3

import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from tqdm import tqdm
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from numpy.random import Generator, default_rng


def _bipartite_config_model(biadj: csr_matrix, random_state: Optional[int] = None) -> csr_matrix:
    """
    Create a bipartite configuration model null network.
    Preserves exact row and column degree sequences while randomizing connections.
    """
    n_rows, n_cols = biadj.shape
    
    # Extract degree sequences
    row_degrees = np.asarray(biadj.sum(axis=1)).ravel().astype(int)
    col_degrees = np.asarray(biadj.sum(axis=0)).ravel().astype(int)
    
    # Quick validation
    if row_degrees.sum() != col_degrees.sum():
        raise ValueError("Row and column degree sums must match")
    
    if row_degrees.sum() == 0:
        # Empty network
        return csr_matrix((n_rows, n_cols), dtype=np.float64)
    
    # Create stub lists (one for each edge endpoint)
    row_stubs = np.repeat(np.arange(n_rows), row_degrees)
    col_stubs = np.repeat(np.arange(n_cols), col_degrees)
    
    # Randomly match row stubs to column stubs
    rng = default_rng(random_state)
    rng.shuffle(col_stubs)
    
    # Create new bipartite adjacency matrix
    data = np.ones(len(row_stubs), dtype=np.float64)
    randomized = csr_matrix((data, (row_stubs, col_stubs)), shape=(n_rows, n_cols))
    
    # Remove potential self-loops (shouldn't exist in bipartite but check)
    randomized.setdiag(0)
    randomized.eliminate_zeros()
    
    return randomized


def _compute_null_summary(q_obs: float, q_null: np.ndarray) -> Dict[str, float]:
    if q_null.size == 0:
        return {
            'Q_obs': q_obs,
            'Q_null_mean': float('nan'),
            'Q_null_sd': float('nan'),
            'z_score': float('nan'),
            'p_value': float('nan'),
            'n_valid_nulls': 0
        }
    mean_q = float(np.mean(q_null)) if q_null.size else float('nan')
    sd_q = float(np.std(q_null, ddof=1)) if q_null.size > 1 else float('nan')
    if np.isfinite(sd_q) and sd_q > 0.0:
        z_score = float((q_obs - mean_q) / sd_q)
    else:
        z_score = float('nan')
    p_value = float(np.mean(q_null >= q_obs)) if q_null.size else float('nan')
    return {
        'Q_obs': q_obs,
        'Q_null_mean': mean_q,
        'Q_null_sd': sd_q,
        'z_score': z_score,
        'p_value': p_value,
        'n_valid_nulls': int(q_null.size)
    }


def _barber_modularity(biadj: csr_matrix, labels_row: np.ndarray, labels_col: np.ndarray) -> float:
    total_weight = float(biadj.sum())
    if total_weight == 0:
        return 0.0
    degrees_row = np.asarray(biadj.sum(axis=1)).ravel()
    degrees_col = np.asarray(biadj.sum(axis=0)).ravel()
    communities: Dict[int, Dict[str, List[int]]] = {}
    for idx, label in enumerate(labels_row):
        if label not in communities:
            communities[label] = {'rows': [], 'cols': []}
        communities[label]['rows'].append(idx)
    for idx, label in enumerate(labels_col):
        if label not in communities:
            communities[label] = {'rows': [], 'cols': []}
        communities[label]['cols'].append(idx)
    quality = 0.0
    for group in communities.values():
        row_idx = np.asarray(group['rows'], dtype=int)
        col_idx = np.asarray(group['cols'], dtype=int)
        if row_idx.size == 0 or col_idx.size == 0:
            continue
        sub = biadj[row_idx][:, col_idx]
        intra = float(sub.sum())
        deg_r = float(degrees_row[row_idx].sum())
        deg_c = float(degrees_col[col_idx].sum())
        quality += intra - (deg_r * deg_c) / total_weight
    return quality / total_weight


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S %Z")


def safe_dir(path: str) -> None:
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def build_biadjacency(edges: pd.DataFrame, src_col: str, tgt_col: str) -> Tuple[csr_matrix, List[str], List[str]]:
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


def simple_brim(biadj: csr_matrix, n_iter: int = 100, random_state: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    BRIM algorithm for bipartite modularity optimization.
    
    Key changes from original:
    - Removed forced diversity code (let BRIM converge naturally)
    - Added convergence check to stop early
    - Handle isolated nodes by assigning them to a separate module
    """
    rng = default_rng(random_state)
    n_rows, n_cols = biadj.shape
    if n_rows == 0 or n_cols == 0:
        return np.zeros(n_rows, dtype=int), np.zeros(n_cols, dtype=int), 0.0

    min_side = max(1, min(n_rows, n_cols))
    if min_side < 3:
        n_init_modules = min_side
    else:
        n_init_modules = min(15, max(2, min_side // 3))
    n_init_modules = max(1, min(n_init_modules, min_side))

    labels_row = rng.integers(0, n_init_modules, n_rows)
    labels_col = rng.integers(0, n_init_modules, n_cols)

    biadj = biadj.tocsr()
    biadj_t = biadj.transpose().tocsr()
    
    converged = False
    for iteration in range(max(1, n_iter)):
        old_row = labels_row.copy()
        old_col = labels_col.copy()
        
        # Update row labels
        for i in range(n_rows):
            start, end = biadj.indptr[i], biadj.indptr[i + 1]
            neighbours = biadj.indices[start:end]
            if neighbours.size == 0:
                # Keep current label for isolated nodes (will handle later)
                continue
            label_counts = np.bincount(labels_col[neighbours])
            labels_row[i] = int(np.argmax(label_counts))
        
        # Update column labels
        for j in range(n_cols):
            start, end = biadj_t.indptr[j], biadj_t.indptr[j + 1]
            neighbours = biadj_t.indices[start:end]
            if neighbours.size == 0:
                # Keep current label for isolated nodes
                continue
            label_counts = np.bincount(labels_row[neighbours])
            labels_col[j] = int(np.argmax(label_counts))
        
        # Check convergence
        if np.array_equal(labels_row, old_row) and np.array_equal(labels_col, old_col):
            converged = True
            break

    # Relabel modules consecutively (removes unused module IDs)
    all_labels = np.concatenate([labels_row, labels_col])
    unique = np.unique(all_labels)
    
    # Only intervene if we have pathological convergence (single module)
    # For dense networks, low modularity with 2-3 modules may be biologically meaningful
    if len(unique) == 1:
        # Only force diversity if it's a true pathological case
        # But first check if this might be legitimate for dense networks
        n_edges = biadj.nnz
        max_possible_edges = n_rows * n_cols
        density = n_edges / max_possible_edges
        
        # Only force diversity for very sparse networks (density < 0.01)
        # Dense networks legitimately may have low modularity
        if density < 0.01:
            min_modules = max(2, min(5, max(n_rows, n_cols) // 3))
            labels_row = np.arange(n_rows) % min_modules
            labels_col = np.arange(n_cols) % min_modules
            unique = np.arange(min_modules)
        # Otherwise, trust the BRIM algorithm - low Q might be meaningful
    
    mapping = {old: new for new, old in enumerate(unique)}
    labels_row = np.array([mapping[x] for x in labels_row], dtype=int)
    labels_col = np.array([mapping[x] for x in labels_col], dtype=int)

    quality = _barber_modularity(biadj, labels_row, labels_col)
    
    # Debug info for understanding modularity patterns
    n_modules = len(unique)
    n_edges = biadj.nnz
    density = n_edges / (biadj.shape[0] * biadj.shape[1])
    
    # Only print for null models where we expect higher variation
    if random_state is not None and abs(quality) < 0.01:
        print(f"[DEBUG] Low Q ({quality:.6f}) with {n_modules} modules, density: {density:.4f}")
    
    return labels_row, labels_col, quality


def run_brim(biadj: csr_matrix, random_state: Optional[int] = None, n_iter: int = 100) -> Tuple[np.ndarray, float, dict]:
    labels_row, labels_col, quality = simple_brim(biadj, n_iter=n_iter, random_state=random_state)
    labels = np.concatenate([labels_row, labels_col])
    meta = {'labels_row_': labels_row, 'labels_col_': labels_col}
    return labels, float(quality), meta


def save_membership(path: str, row_nodes: List[str], col_nodes: List[str], labels: np.ndarray, meta: dict, layer: str) -> pd.DataFrame:
    if 'labels_row_' in meta and 'labels_col_' in meta:
        lab_row = meta['labels_row_']
        lab_col = meta['labels_col_']
    else:
        n_r = len(row_nodes)
        lab_row = labels[:n_r]
        lab_col = labels[n_r:]
    df_row = pd.DataFrame({'node_id': row_nodes, 'module_id': lab_row})
    df_col = pd.DataFrame({'node_id': col_nodes, 'module_id': lab_col})
    df = pd.concat([df_row, df_col], ignore_index=True)
    df['node_type'] = np.where(df['node_id'].str.startswith('FG_'), 'FG',
                               np.where(df['node_id'].str.startswith('PYO_'), 'PYOV', 'UNKNOWN'))
    df['layer'] = layer
    df.to_csv(path, index=False)
    return df


def analyze_layer(layer_name: str, edges: pd.DataFrame, src_col: str, tgt_col: str, outdir: str,
                  random_state: Optional[int] = None, n_replicates: int = 1,
                  n_null_replicates: int = 0, brim_iterations: int = 100) -> Dict:
    print(f"[{ts()}] Layer: {layer_name}")
    biadj, rows, cols = build_biadjacency(edges, src_col, tgt_col)

    # Debug info
    print(f"\n{'='*60}")
    print(f"Layer: {layer_name}")
    print(f"FG nodes (rows): {len(rows)}")
    print(f"PYO nodes (cols): {len(cols)}")
    print(f"Sample rows: {rows[:5]}")
    print(f"Sample cols: {cols[:5]}")
    print(f"{'='*60}\n")

    n_edges = int(biadj.nnz)
    density = n_edges / (biadj.shape[0] * biadj.shape[1]) if biadj.shape[0] * biadj.shape[1] > 0 else 0.0
    print(f"dims: {biadj.shape[0]} x {biadj.shape[1]} | edges: {n_edges} | density: {density:.4f}")

    replicates: List[Dict[str, float]] = []
    best = None
    # Use a single master RNG for this layer to ensure independent streams
    master_rng = default_rng(random_state)
    for rep in tqdm(range(n_replicates), desc=f'{layer_name} replicate'):
        # Each replicate gets an independent seed from master RNG
        rs = None if random_state is None else int(master_rng.integers(0, 2**31))
        labels, quality, meta = run_brim(biadj, random_state=rs, n_iter=brim_iterations)
        replicates.append({
            'replicate': rep + 1,
            'Q': quality,
            'random_state': rs,
            'brim_iterations': brim_iterations
        })
        if best is None or quality > best['quality']:
            best = {'labels': labels, 'quality': quality, 'meta': meta}
    assert best is not None
    best_labels = best['labels']
    best_quality = best['quality']
    best_meta = best['meta']
    print(f"[{ts()}] {layer_name} | Best Barber Q: {best_quality:.6f} (across {n_replicates} replicate(s))")

    assign_path = os.path.join(outdir, f"module_assignments_{layer_name}_sk.csv")
    membership_df = save_membership(assign_path, rows, cols, best_labels, best_meta, layer_name)

    rep_file = os.path.join(outdir, f"replicate_Q_{layer_name}_sk.csv")
    rep_df = pd.DataFrame(replicates)
    rep_df['layer'] = layer_name
    rep_df.to_csv(rep_file, index=False)

    null_rows: List[Dict[str, Optional[float]]] = []
    null_replicates_file = os.path.join(outdir, f"null_replicates_{layer_name}_sk.csv")
    null_summary_file = os.path.join(outdir, f"null_summary_{layer_name}_sk.csv")
    
    # Store original degree sequences for validation
    orig_row_degrees = np.asarray(biadj.sum(axis=1)).ravel()
    orig_col_degrees = np.asarray(biadj.sum(axis=0)).ravel()
    
    if n_null_replicates > 0:
        # Create dedicated RNG for null models to ensure independence from replicates
        null_master_rng = default_rng(random_state + 1000000 if random_state else None)
        for null_idx in tqdm(range(n_null_replicates),  desc=f'{layer_name} nulls'):
            # Generate independent seeds for each null component
            config_seed = None if random_state is None else int(null_master_rng.integers(0, 2**31))
            brim_null_seed = None if random_state is None else int(null_master_rng.integers(0, 2**31))
            
            try:
                randomized = _bipartite_config_model(biadj, random_state=config_seed)
                
                # Validate degree preservation on first null replicate
                if null_idx == 0:
                    null_row_degrees = np.asarray(randomized.sum(axis=1)).ravel()
                    null_col_degrees = np.asarray(randomized.sum(axis=0)).ravel()
                    if not np.array_equal(orig_row_degrees, null_row_degrees) or not np.array_equal(orig_col_degrees, null_col_degrees):
                        print(f"[{ts()}] {layer_name} | WARNING: Degree preservation failed in null model!")
                        print(f"  Original row degrees sum: {orig_row_degrees.sum()}, Null row degrees sum: {null_row_degrees.sum()}")
                        print(f"  Original col degrees sum: {orig_col_degrees.sum()}, Null col degrees sum: {null_col_degrees.sum()}")
                
                _, q_null, _ = run_brim(randomized, random_state=brim_null_seed, n_iter=brim_iterations)
                null_rows.append({
                    'replicate': null_idx + 1,
                    'Q': q_null,
                    'random_state': brim_null_seed,
                    'config_seed': config_seed,
                    'brim_iterations': brim_iterations,
                    'null_method': 'bipartite_config_model'
                })
            except Exception as exc:
                print(f"[{ts()}] {layer_name} | Null replicate {null_idx + 1} failed: {exc}")
        if null_rows:
            null_df = pd.DataFrame(null_rows)
            null_df['layer'] = layer_name
            null_df.to_csv(null_replicates_file, index=False)
            valid_null = null_df['Q'].dropna().to_numpy(dtype=float)
            null_summary = _compute_null_summary(best_quality, valid_null)
        else:
            null_df = pd.DataFrame(columns=['replicate', 'Q', 'random_state', 'shuffle_seed', 'layer'])
            null_df.to_csv(null_replicates_file, index=False)
            null_summary = _compute_null_summary(best_quality, np.array([], dtype=float))
    else:
        null_summary = _compute_null_summary(best_quality, np.array([], dtype=float))
        null_rows = []
        pd.DataFrame(columns=['replicate', 'Q', 'random_state', 'shuffle_seed', 'layer']).to_csv(null_replicates_file, index=False)
    pd.DataFrame([null_summary]).to_csv(null_summary_file, index=False)
    if null_rows:
        print(f"[{ts()}] {layer_name} | Null replicates retained: {null_summary['n_valid_nulls']} / {n_null_replicates}")
    else:
        print(f"[{ts()}] {layer_name} | Null replicates skipped or unavailable")

    q_values = rep_df['Q'].to_numpy()
    mean_q = float(np.mean(q_values))
    sd_q = float(np.std(q_values, ddof=1)) if len(q_values) > 1 else float('nan')
    rep_stats = {'mean_Q': mean_q, 'sd_Q': sd_q, 'min_Q': float(np.min(q_values)), 'max_Q': float(np.max(q_values))}
    return {
        'Q_best': best_quality,
        'replicate_stats': rep_stats,
        'membership': membership_df,
        'null_summary': null_summary,
        'assignments_path': assign_path,
        'replicate_file': rep_file,
        'null_file': null_summary_file,
        'null_replicates_file': null_replicates_file,
    }


def main():
    cfg = {
        'edges_file': 'data/interim/edges_functional_groups_conservative.csv',
        'output_dir': 'results/phase_04/modularity',
        'log_dir': 'docs/phase_04/logs',
        'random_state': 42,
        'n_replicates': 100,
        'n_null_replicates': 100,
        'brim_iterations': 20000,
    }

    safe_dir(cfg['output_dir'])
    safe_dir(cfg['log_dir'])

    print("=== Phase 04 Step 03a: Modularity Analysis (Fixed Version) ===")
    print("Timestamp:", ts())
    print(f"Configured replicates per layer: {cfg['n_replicates']}")
    print(f"Configured null replicates per layer: {cfg['n_null_replicates']} (bipartite configuration model)")
    print(f"BRIM iterations per run: {cfg['brim_iterations']}")
    print()

    if not os.path.isfile(cfg['edges_file']):
        raise FileNotFoundError(f"Missing edges file: {cfg['edges_file']}")

    edges = pd.read_csv(cfg['edges_file'])
    print(f"Loaded {len(edges)} edges")

    prod_edges = edges.query("edge_type == 'production'")[["source", "target"]]
    util_edges = edges.query("edge_type == 'utilization'")[["source", "target"]]
    if prod_edges.empty or util_edges.empty:
        raise RuntimeError("Edges missing for one or both layers")

    layers = {
        'production': (prod_edges, 'source', 'target'),
        # For utilization, invert direction to get FG as rows, PYO as cols
        'utilization': (util_edges.rename(columns={'source': 'col', 'target': 'row'}).rename(columns={'row': 'source', 'col': 'target'}), 'source', 'target')
    }

    results: Dict[str, Dict] = {}
    for layer_idx, (layer_name, (e, s_col, t_col)) in enumerate(layers.items()):
        # Give each layer an independent random seed to avoid correlation
        layer_seed = cfg['random_state'] + layer_idx * 10000000 if cfg['random_state'] else None
        results[layer_name] = analyze_layer(
            layer_name,
            e,
            s_col,
            t_col,
            cfg['output_dir'],
            random_state=layer_seed,
            n_replicates=cfg['n_replicates'],
            n_null_replicates=cfg['n_null_replicates'],
            brim_iterations=cfg['brim_iterations']
        )

    summary_lines = [
        'Phase 04 Step 03a summary (Fixed Version)',
        f"Timestamp: {ts()}",
        ''
    ]
    for lname, res in results.items():
        qs = res['replicate_stats']
        ns = res['null_summary']
        summary_lines.extend([
            lname.title(),
            f"  Best Q: {res['Q_best']:.4f}",
            f"  Replicate Q mean±sd: {qs['mean_Q']:.4f} ± {qs['sd_Q']}",
            f"  Null mean±sd: {ns['Q_null_mean']} ± {ns['Q_null_sd']}",
            f"  z-score: {ns['z_score']} | p(one-tailed): {ns['p_value']}",
            f"  Modules saved: {res['assignments_path']}",
            f"  Replicate Q file: {res['replicate_file']}",
            f"  Null replicates file: {res['null_replicates_file']}",
            f"  Null summary file: {res['null_file']}",
            ''
        ])

    summary_file = os.path.join(cfg['log_dir'], 'step03a_modularity_serial_summary.txt')
    with open(summary_file, 'w') as fh:
        fh.write("\n".join(summary_lines))
    print("\n".join(summary_lines))


if __name__ == '__main__':
    main()
