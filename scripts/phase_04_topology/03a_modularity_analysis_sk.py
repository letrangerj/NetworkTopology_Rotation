#!/usr/bin/env python3

import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from sknetwork.clustering import Louvain


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


def run_bipartite_louvain(biadj: csr_matrix, random_state: Optional[int] = None) -> Tuple[np.ndarray, float, dict]:
    algo = Louvain(random_state=random_state)
    labels = algo.fit_predict(biadj)
    if hasattr(algo, 'labels_row_') and hasattr(algo, 'labels_col_'):
        lab_row = np.array(algo.labels_row_, dtype=int)
        lab_col = np.array(algo.labels_col_, dtype=int)
    else:
        n_rows = biadj.shape[0]
        n_cols = biadj.shape[1]
        if len(labels) != n_rows + n_cols:
            raise ValueError("Unexpected label length from Louvain output for bipartite graph")
        lab_row = np.array(labels[:n_rows], dtype=int)
        lab_col = np.array(labels[n_rows:], dtype=int)
    quality = getattr(algo, 'modularity_', None)
    if quality is None:
        quality = getattr(algo, 'score_', None)
    if quality is None and hasattr(algo, 'score'):
        try:
            quality = algo.score(biadj, labels)
        except TypeError:
            quality = algo.score(biadj)
    if quality is None:
        try:
            from sknetwork.utils.metrics import modularity as modularity_fn
        except ImportError:
            try:
                from sknetwork.utils import modularity as modularity_fn  # type: ignore
            except ImportError:
                modularity_fn = None
        if modularity_fn is not None:
            try:
                quality = modularity_fn(biadj, labels)
            except TypeError:
                quality = modularity_fn(biadj, lab_row, lab_col)
    if quality is None:
        quality = _barber_modularity(biadj, lab_row, lab_col)
    quality = float(quality)
    meta = {}
    meta['labels_row_'] = lab_row
    meta['labels_col_'] = lab_col
    return np.array(labels, dtype=int), quality, meta


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
                  random_state: Optional[int] = None, n_replicates: int = 1) -> Dict:
    print(f"[{ts()}] Layer: {layer_name}")
    biadj, rows, cols = build_biadjacency(edges, src_col, tgt_col)
    n_edges = int(biadj.nnz)
    density = n_edges / (biadj.shape[0] * biadj.shape[1]) if biadj.shape[0] * biadj.shape[1] > 0 else 0.0
    print(f"dims: {biadj.shape[0]} x {biadj.shape[1]} | edges: {n_edges} | density: {density:.4f}")

    replicates: List[Dict[str, float]] = []
    best = None
    for rep in range(n_replicates):
        rs = None if random_state is None else random_state + rep
        labels, quality, meta = run_bipartite_louvain(biadj, random_state=rs)
        replicates.append({'replicate': rep + 1, 'Q': quality, 'random_state': rs})
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

    null_file = os.path.join(outdir, f"null_summary_{layer_name}_sk.csv")
    pd.DataFrame([{'Q_obs': best_quality, 'Q_null_mean': np.nan, 'Q_null_sd': np.nan, 'z_score': np.nan, 'p_value': np.nan, 'n_valid_nulls': 0}]).to_csv(null_file, index=False)

    q_values = rep_df['Q'].to_numpy()
    mean_q = float(np.mean(q_values))
    sd_q = float(np.std(q_values, ddof=1)) if len(q_values) > 1 else float('nan')
    rep_stats = {'mean_Q': mean_q, 'sd_Q': sd_q, 'min_Q': float(np.min(q_values)), 'max_Q': float(np.max(q_values))}
    return {
        'Q_best': best_quality,
        'replicate_stats': rep_stats,
        'membership': membership_df,
        'null_summary': {'Q_obs': best_quality, 'Q_null_mean': np.nan, 'Q_null_sd': np.nan, 'z_score': np.nan, 'p_value': np.nan, 'n_valid_nulls': 0},
        'assignments_path': assign_path,
        'replicate_file': rep_file,
        'null_file': null_file,
    }


def main():
    cfg = {
        'edges_file': 'data/interim/edges_functional_groups_conservative.csv',
        'output_dir': 'results/phase_04/modularity',
        'log_dir': 'docs/phase_04/logs',
        'random_state': 42,
        'n_replicates': 100,
    }

    safe_dir(cfg['output_dir'])
    safe_dir(cfg['log_dir'])

    print("=== Phase 04 Step 03a: Modularity Analysis (scikit-network) ===")
    print("Timestamp:", ts())
    print(f"Configured replicates per layer: {cfg['n_replicates']}")
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
    for layer_name, (e, s_col, t_col) in layers.items():
        results[layer_name] = analyze_layer(
            layer_name,
            e,
            s_col,
            t_col,
            cfg['output_dir'],
            random_state=cfg['random_state'],
            n_replicates=cfg['n_replicates']
        )

    summary_lines = [
        'Phase 04 Step 03a summary (scikit-network)',
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
            f"  Null summary file: {res['null_file']}",
            ''
        ])

    summary_file = os.path.join(cfg['log_dir'], 'step03a_modularity_serial_summary.txt')
    with open(summary_file, 'w') as fh:
        fh.write("\n".join(summary_lines))
    print("\n".join(summary_lines))


if __name__ == '__main__':
    main()
