# 图注 — Pyoverdine 介导的网络（Panel 3A）

**图3A：Pyoverdine介导的铁竞争二分网络。** 该网络揭示了铁载体相互作用的社会结构，由两类不同的节点组成：**功能组（FGs，圆形）**，代表了1809个细菌菌株根据其独特的生产和利用谱系压缩而成的177个群体；以及**Pyoverdine组（PYOs，菱形）**，代表了从基因数据中验证的93种不同的铁载体类型。有向边表示公共物品的流动：从FG指向PYO的深色箭头表示**生产**，而从PYO指向FG的浅色箭头表示**利用**。该网络包含175条生产边和2238条利用边。每个FG节点的大小与其所代表的菌株数量的平方根成正比，突显了主导的生态策略；而PYO节点的大小反映其总度数（生产者和消费者的总和），表明其作为公共物品的普遍性。FG节点的颜色标示了其基于生产/利用配置文件的生态策略：单受体生产者（绿色）、多受体生产者（蓝色）、非生产者（红色）以及其他专业类型。空间布局由Fruchterman-Reingold力导向算法确定，该算法将紧密互连的节点聚集在一起，以直观地揭示模块化结构并识别相互作用的核心枢纽。

---

### 详细分解与技术说明

描述：
该二分网络可视化展示功能组（FG，圆点）和 pyoverdine（PYO，六角形样标记）节点，及表示产生（FG → PYO）和利用（PYO → FG）的有向边。

数据与溯源：
- 脚本：[`scripts/phase_03_network/04_visualize_network.R`](scripts/phase_03_network/04_visualize_network.R:1)
- 输入文件：[`data/interim/edges_functional_groups_conservative.csv`](data/interim/edges_functional_groups_conservative.csv:1)、[`data/interim/nodes_pyoverdines_conservative.csv`](data/interim/nodes_pyoverdines_conservative.csv:1)、[`data/interim/nodes_functional_groups_conservative.csv`](data/interim/nodes_functional_groups_conservative.csv:1)
- 输出图像：[`figures/network/fg_pyov_network_panel3A.pdf`](figures/network/fg_pyov_network_panel3A.pdf:1)、[`figures/network/fg_pyov_network_panel3A.png`](figures/network/fg_pyov_network_panel3A.png:1)
- 随机种子：全局 set.seed(2025)，布局种子 layout seed = 42（见脚本）

可视编码：
- FG 节点：圆形，填充色表示策略类别：
  * 绿色 = Single-receptor producer（单受体生产者）
  * 橙色 = Multi-producer（多产者）
  * 蓝色 = Multi-receptor producer（多受体生产者）
  * 紫色 = Producer-only（仅生产，未利用）
  * 红色 = Nonproducer（非生产者）
- PYO 节点：六角形样标记（shape 23），填充橙色
- 节点大小：FG 节点大小 ∝ sqrt(strain_count)；PYO 节点大小 ∝ 使用/计数代理值
- 边：有向箭头段，生产边以 gray30 绘制，利用边以 gray50 绘制；段首尾缩短以避免箭头覆盖目标节点
- 布局：力导向 Fruchterman–Reingold（igraph::layout_with_fr，niter = 500），坐标通过种子固定以保证可复现

关键统计（FG 级网络，来自可视化度量）：
- 顶点数：270（177 个功能组 + 93 个 pyoverdine）
- 唯一边数：2,413（生产边 FG → PYO = 175；利用边 PYO → FG = 2,238）
- 弱连通分量：2（大型主分量 + 小分量）— 参见 [`docs/phase_03/phase_03_summary.md`](docs/phase_03/phase_03_summary.md:1)

如何解读该图：
- 枢纽（hubs）：大量外向利用边的 PYO（大六角）表示被广泛利用的 siderophore，可能作为公共物质。
- 生产枢纽：具有多个外向生产边的 FG（大圆）表示能产生多种 pyoverdine 的功能组。
- 策略区分：节点颜色展示生态策略类别，可用于识别模块是以生产者、非生产者或混合策略为主。
- 空间聚类：力导向布局的邻近性与边密度提示模块与枢纽－辐射结构；应与第4阶段的模块性与共分配稳定性结果联合解读。

可复现性与技术说明：
- 使用 ggplot2 绘图，边先绘制、节点覆盖于上以减少遮挡；对线段进行缩短处理以避让箭头和节点标记。
- 可视化度量保存于 [`docs/phase_03/visualization_metrics.csv`](docs/phase_03/visualization_metrics.csv:1)，文本摘要见 [`docs/phase_03/visualization_summary.txt`](docs/phase_03/visualization_summary.txt:1)。

生物学解读与建议用途：
- 该面板可作为保守 FG 级网络结构的描述性摘要：识别候选公共物质 pyoverdines、多产者枢纽以及候选通用消费者。
- 将此图与拓扑指标（模块性、嵌套性、度分布）及共分配稳定性分析结合，以支持关于模块稳健性與生态角色的论断。
- 稿件建议：在主图中提供颜色图例与简要编码说明，并在补充材料中给出共分配热图以证明模块稳定性。

交叉引用文件：
- 实现与溯源：[`scripts/phase_03_network/04_visualize_network.R`](scripts/phase_03_network/04_visualize_network.R:1)
- Phase3 汇总：[`docs/phase_03/phase_03_summary.md`](docs/phase_03/phase_03_summary.md:1)
