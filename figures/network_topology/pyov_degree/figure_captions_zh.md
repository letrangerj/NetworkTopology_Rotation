# 图注 — Pyoverdine 度分布（步骤02b）
# Path: figures/network_topology/pyov_degree/figure_captions_zh.md
# Generated: 2025-10-13
# Author: automated-agent
#
# 说明:
# - 每个图注对应 `figures/network_topology/pyov_degree` 下的 PNG 文件。
# - 源脚本: [`scripts/phase_04_topology/02b_pyoverdine_degree_distributions.R`](scripts/phase_04_topology/02b_pyoverdine_degree_distributions.R:1)
# - 数据来源: `results/phase_04/degree/pyov_degree_summary.csv` 与 `data/interim/nodes_pyoverdines_conservative.csv`
#
# 总览（基于实际结果）
# - 节点数: 93 经验证 pyoverdine 节点（Step 2 将配对展开为 synthetase 级别）
# - 关键统计（示例）:
#   - 最大 FG 利用出度 (out_fg) = 61 : `pyov_53` (内部 `PYO_32`)
#   - 次高 FG 利用出度 = 58 : `pyov_27` (内部 `PYO_14`)
#   - 最大 FG 生产入度 (in_fg) = 12 : `pyov_23` (内部 `PYO_13`)
#   - 最大 STR 生产入度 (in_str) = 768 : `pyov_22` (内部 `PYO_12`)
#   - 最大 STR 利用出度 (out_str) = 869 : `pyov_22`/`pyov_23`（见表）
#
# 方法说明
# - CCDF 计算：脚本按降序排列非零度并计算 P(X ≥ x) = rank / n，图中在对数–对数坐标上绘制以检验重尾性
# - 解释要点：在 log–log CCDF 上线性尾部指示幂律样重尾；直线斜率与幂律指数 γ 相关（斜率 = 1-γ）
#
# 单图注
[`figures/network_topology/pyov_degree/pyov_degree_in_fg.png`](figures/network_topology/pyov_degree/pyov_degree_in_fg.png:1)
Caption:
功能组（FG）级别 — Pyoverdine 生产入度直方图。
- 内容: 展示每个 pyov 被多少功能组产生（`in_fg` 列）。
- 实际结果: 大多数 pyov 的 `in_fg` = 1，少数达到更高值；最大观测值为 12（`pyov_23` / `PYO_13`），提示某些合成簇在多个功能组中独立存在或被重复保留。
- 说明: 生产入度较为集中，表明合成能力在功能组层面受限，支持"生产较专一"的结论。

[`figures/network_topology/pyov_degree/pyov_degree_out_fg.png`](figures/network_topology/pyov_degree/pyov_degree_out_fg.png:1)
Caption:
功能组（FG）级别 — Pyoverdine 利用出度直方图。
- 内容: 展示每个 pyov 被多少功能组利用（`out_fg` 列）。
- 实际结果: 出度分布右偏，少数 pyov 被大量功能组利用；例子包括 `pyov_53`/`PYO_32`（out_fg=61）和 `pyov_27`/`PYO_14`（out_fg=58）。
- 说明: 该直方图与 CCDF 配合表明受体泛化在功能组层面存在显著不对称。

[`figures/network_topology/pyov_degree/pyov_degree_in_fg_ccdf.png`](figures/network_topology/pyov_degree/pyov_degree_in_fg_ccdf.png:1)
Caption:
功能组（FG）级别 — 生产入度 CCDF（对数–对数）。
- 内容: 绘制 P(X ≥ x) 对 degree 的 CCDF（log–log）。
- 实际结果: 入度尾部短且集中，尾部不显著线性，暗示生产入度不显著重尾。

[`figures/network_topology/pyov_degree/pyov_degree_out_fg_ccdf.png`](figures/network_topology/pyov_degree/pyov_degree_out_fg_ccdf.png:1)
Caption:
功能组（FG）級别 — 利用出度 CCDF（对数–对数）。
- 内容: 绘制 P(X ≥ x) 对 degree 的 CCDF（log–log）。
- 实际结果: CCDF 在尾部呈缓慢下降，接近线性段，支持重尾/幂律样行为；例如 P(X ≥ 61) 对应 `pyov_53`，尾部缓慢下降表明有少数高影响 pyov。
- 建议: 若需定量检验幂律，使用极大似然拟合与 KS 检验（Clauset et al., 2009）。

[`figures/network_topology/pyov_degree/pyov_degree_in_str.png`](figures/network_topology/pyov_degree/pyov_degree_in_str.png:1)
Caption:
菌株（STR）级别 — Pyoverdine 生产入度直方图。
- 内容: 展示每个 pyov 在菌株层面的生产计数（`in_str` 列）。
- 实际结果: 菌株层面变异更大；例如 `pyov_22`/`PYO_12` 的 `in_str` = 768，显示部分 pyov 在菌株范围内被广泛合成。
- 说明: 菌株尺度揭示了在不同菌株间合成能力的高度不均一。

[`figures/network_topology/pyov_degree/pyov_degree_out_str.png`](figures/network_topology/pyov_degree/pyov_degree_out_str.png:1)
Caption:
菌株（STR）级别 — Pyoverdine 利用出度直方图。
- 内容: 展示每个 pyov 在菌株层面的使用者计数（`out_str` 列）。
- 实际结果: 出度极端值显著；最大 `out_str` 可达 869（例如 `pyov_22` / `PYO_12` 或 `pyov_23` / `PYO_13`），表明在菌株尺度上存在广泛可利用的 pyov。
- 说明: 该图支持菌株级别的重尾行为觀察。

[`figures/network_topology/pyov_degree/pyov_degree_in_str_ccdf.png`](figures/network_topology/pyov_degree/pyov_degree_in_str_ccdf.png:1)
Caption:
菌株（STR）级别 — 生产入度 CCDF（对数–对数）。
- 内容: CCDF 显示菌株尺度的累积尾概率。
- 实际结果: 对数–对数图尾部近似线性，支持菌株级别强重尾假设，提示极少数 pyov 在菌株层面被大量合成。

[`figures/network_topology/pyov_degree/pyov_degree_out_str_ccdf.png`](figures/network_topology/pyov_degree/pyov_degree_out_str_ccdf.png:1)
Caption:
菌株（STR）級别 — 利用出度 CCDF（对数–对数）。
- 内容: CCDF 显示菌株尺度上 P(X ≥ x)。
- 实际结果: 尾部缓慢下降并接近线性段，匹配 `pyov_22`/`PYO_12`、`pyov_23`/`PYO_13` 等高连接节点，进一步支持重尾行为。

[`figures/network_topology/pyov_degree/pyov_degree_ccdf_fg_combined.png`](figures/network_topology/pyov_degree/pyov_degree_ccdf_fg_combined.png:1)
Caption:
FG 级别 CCDF 合成图 — 并列展示生产/利用 CCDF 以比较尾部强度。基于实际数据，利用（out_fg）显示明显更强的重尾性。

[`figures/network_topology/pyov_degree/pyov_degree_ccdf_str_combined.png`](figures/network_topology/pyov_degree/pyov_degree_ccdf_str_combined.png:1)
Caption:
STR 级别 CCDF 合成图 — 并列展示菌株级生产/利用 CCDF，强调菌株尺度的极端异质性和个别高连接 pyov 的主导作用。

# 参考与方法
# - 数据: `results/phase_04/degree/pyov_degree_summary.csv` 与 `data/interim/nodes_pyoverdines_conservative.csv`
# - CCDF 方法: 在脚本中以降序排列非零度并计算 P(X ≥ x) = rank / n（见 `02b_pyoverdine_degree_distributions.R`）。
# - 幂律拟合建议: 使用极大似然估计与 KS 检验（Clauset et al., 2009）；如需我可添加拟合代码与拟合统计。
#
# End of file.