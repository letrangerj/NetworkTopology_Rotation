#!/usr/bin/env python3

import os
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from numpy.random import default_rng
from typing import Optional
from tqdm import tqdm


def edge_rewiring_model(
    biadj: csr_matrix,
    n_swaps: int = 10,
    random_state: Optional[int] = None,
    return_stats: bool = False,
):
    """
    二分网络的边重排（Maslov-Sneppen式交换）：
    - 维持行、列度序不变
    - n_swaps 表示期望的“成功交换次数”（不再按边数缩放）
    - 不进行自环检查（对二分图无意义）；仅禁止产生重复边



    参数
    - biadj: 原始二分邻接矩阵（CSR）
    - n_swaps: 成功交换次数的目标值（绝对数量）
    - random_state: 随机种子
    - return_stats: 若为 True，则返回 (随机矩阵, 成功交换次数)


    返回

    - 随机重排后的二分邻接矩阵（CSR），或当 return_stats=True 时返回 (矩阵, 成功交换次数)

    """

    n_rows, n_cols = biadj.shape
    if biadj.nnz == 0 or n_rows == 0 or n_cols == 0 or n_swaps <= 0:
        return biadj.copy()

    rng = default_rng(random_state)

    # 提取边列表
    rows, cols = biadj.nonzero()
    edges = list(zip(rows, cols))
    current_edges = set(edges)

    successful_swaps = 0
    # 给每次交换留出足够尝试次数（如果网络较稀或结构限制多，成功率会偏低）
    max_attempts = max(10 * n_swaps, 1000)

    for _ in range(max_attempts):
        if successful_swaps >= n_swaps or len(edges) < 2:
            break

        # 随机选择两条不同的边
        idx1, idx2 = rng.choice(len(edges), size=2, replace=False)
        r1, c1 = edges[idx1]
        r2, c2 = edges[idx2]

        # 拟议交换： (r1,c1),(r2,c2) -> (r1,c2),(r2,c1)
        new_edge1 = (r1, c2)
        new_edge2 = (r2, c1)

        # 禁止产生已存在的边或重复集合（保持简单一致性）
        if new_edge1 in current_edges or new_edge2 in current_edges:
            continue

        # 执行交换并更新集合
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


def build_biadjacency(edges: pd.DataFrame, src_col: str, tgt_col: str):
    """构建二分邻接矩阵"""
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


def simple_brim(
    biadj: csr_matrix,
    n_iter: int = 100,
    random_state: Optional[int] = None,
    n_restarts: int = 1,
    min_modules: int = 1,
) -> float:
    """

    简化的BRIM算法，只返回模块性值（Barber bipartite modularity）

    - 支持多次随机初始化（n_restarts），返回最佳Q
    - 可要求最少模块数（min_modules），不足则视为无效尝试
    """

    rng = default_rng(random_state)
    n_rows, n_cols = biadj.shape
    if n_rows == 0 or n_cols == 0:
        return 0.0

    biadj = biadj.tocsr()
    biadj_t = biadj.transpose().tocsr()
    best_q = -np.inf

    for _restart in range(max(1, int(n_restarts))):
        # 随机初始化模块数量和标签
        min_side = max(1, min(n_rows, n_cols))

        if min_side < 3:
            n_init_modules = min_side
        else:
            n_init_modules = min(15, max(2, min_side // 3))

        n_init_modules = max(1, min(n_init_modules, min_side))
        labels_row = rng.integers(0, n_init_modules, n_rows)
        labels_col = rng.integers(0, n_init_modules, n_cols)

        # 交替更新直到收敛或达到迭代次数
        for _ in range(max(1, n_iter)):
            old_row = labels_row.copy()
            old_col = labels_col.copy()

            # 更新行标签

            for i in range(n_rows):
                start, end = biadj.indptr[i], biadj.indptr[i + 1]
                neighbours = biadj.indices[start:end]
                if neighbours.size == 0:
                    continue
                label_counts = np.bincount(labels_col[neighbours])
                labels_row[i] = int(np.argmax(label_counts))

            # 更新列标签

            for j in range(n_cols):
                start, end = biadj_t.indptr[j], biadj_t.indptr[j + 1]
                neighbours = biadj_t.indices[start:end]
                if neighbours.size == 0:
                    continue
                label_counts = np.bincount(labels_row[neighbours])
                labels_col[j] = int(np.argmax(label_counts))
            # 收敛检查

            if np.array_equal(labels_row, old_row) and np.array_equal(
                labels_col, old_col
            ):
                break

        # 计算Barber模块性 Q

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

            # 仅计数同时包含行和列的模块
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
                q_val = -np.inf  # 不足最小模块数时视为无效
            else:
                q_val = quality / total_weight

        if q_val > best_q:
            best_q = q_val

    return 0.0 if not np.isfinite(best_q) else float(best_q)


def main():
    # 配置参数

    edges_file = "data/interim/edges_functional_groups_conservative.csv"
    random_state = 42
    n_null_replicates = 100
    brim_iterations = 20000

    # 现在 n_swaps 表示“成功交换次数”的绝对值

    n_swaps = 1000
    # 可选：使用每条边的交换次数来设定混合强度。设为 None 时使用 n_swaps 的绝对次数
    swaps_per_edge = None  # 例如 10.0 表示每条边交换约10次；None 则使用 n_swaps
    # BRIM 多起点与最小模块数约束
    n_restarts = 20  # 每个replicate的随机初始化次数，取最佳Q
    min_modules = 2  # 最少有效模块数（同时包含行和列的模块）
    # 是否记录每个replicate的交换次数与Q
    log_replicate_stats = False

    print("=== Utilization网络边重排Null Test ===")

    # 加载数据
    edges = pd.read_csv(edges_file)
    util_edges = edges.query("edge_type == 'utilization'")[["source", "target"]]

    # 反转方向，使FG为行，PYO为列
    util_edges = util_edges.rename(columns={"source": "col", "target": "row"}).rename(
        columns={"row": "source", "col": "target"}
    )

    # 构建二分邻接矩阵
    biadj, rows, cols = build_biadjacency(util_edges, "source", "target")

    # 计算原始网络的模块性
    q_obs = simple_brim(biadj, n_iter=brim_iterations, random_state=random_state)

    # 计算网络密度
    n_edges = int(biadj.nnz)
    max_edges = biadj.shape[0] * biadj.shape[1]
    density = n_edges / max_edges if max_edges > 0 else 0.0

    print(f"网络规模: {biadj.shape[0]} x {biadj.shape[1]}")
    print(f"边数: {n_edges}")
    print(f"网络密度: {density:.4f}")
    print(f"原始模块性 (Q_obs): {q_obs:.6f}")

    # 运行null模型
    null_q_values = []
    master_rng = default_rng(random_state)

    # 计算每次重排的目标交换次数
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

    for idx in tqdm(range(n_null_replicates), desc="Null模型进度"):
        # 为每个 replicate 的重排和BRIM使用不同的种子，解耦两者
        rs_rewire = int(master_rng.integers(0, 2**31 - 1))
        rs_brim = int(master_rng.integers(0, 2**31 - 1))

        # 生成边重排网络（可返回成功交换次数）

        if log_replicate_stats:
            randomized, swaps_done = edge_rewiring_model(
                biadj, n_swaps=target_swaps, random_state=rs_rewire, return_stats=True
            )
        else:
            randomized = edge_rewiring_model(
                biadj, n_swaps=target_swaps, random_state=rs_rewire, return_stats=False
            )

        # 计算模块性（多起点，最小模块数约束）

        q_null = simple_brim(
            randomized,
            n_iter=brim_iterations,
            random_state=rs_brim,
            n_restarts=n_restarts,
            min_modules=min_modules,
        )
        null_q_values.append(q_null)

        if log_replicate_stats:
            print(f"[replicate {idx + 1}] swaps={swaps_done} Q_null={q_null:.6f}")

    # 统计
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
    print(f"BRIM重启次数: {n_restarts}, 最小模块数: {min_modules}")
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
    log_file = os.path.join(log_dir, "03g_edge_rewiring_null_test_results.txt")

    now = datetime.now().isoformat()
    lines = []
    lines.append(f"Timestamp: {now}")
    lines.append("=== Utilization 网络边重排 Null Test Results ===")
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
