## Li & Hu 2014 实验参数参考

来源: *Global migration strategy with moving colony for hierarchical distributed evolutionary algorithms*, Soft Computing 18:2161–2176 (2014)

### 算法列表

| 缩写 | 全称 |
|---|---|
| Ring DEA | Ring topology distributed EA |
| Ring individual HDEA | Ring local + Ring individual global migration HDEA |
| Random individual HDEA | Ring local + Random individual global migration HDEA |
| Ring colony HDEA | Ring local + Ring colony global migration (本文提出) |
| Random colony HDEA | Ring local + Random colony global migration (本文提出) |

### 基础 EA（Inver-over GA for TSP）

- 算子: inver-over (变异/交叉自适应) + mapping 算子 (交叉)
- 临界速度 (critical velocity): 5000 — 低于此阈值时以给定概率执行 mapping
- 变异率 P: 初始 0.02，按 `P = 0.02 × (1 − g_n/g × 0.01)` 线性衰减
- 交叉率: `1 − P`（起始 ~0.98）
- 映射率: 0.05（固定）

### 并行/迁移参数

| 参数 | 值 |
|---|---|
| 子种群大小 | 100 |
| 子种群数量 | 64（8×8 组，每组 8 子种群；对比实验也用 16 和 36） |
| 局部迁移拓扑 | Ring |
| 局部迁移策略 | Random-worst（随机选迁出个体，替换目标子种群最差个体） |
| 局部迁移规模 | 1 个体 |
| 局部:全局迁移比 (HDEA) | 20:1（每 20 次局部迁移 → 1 次全局迁移） |
| Ring colony 全局迁移 | 子种群级，1 子种群/次，轮转 Ring 拓扑，无通信开销 |
| Random colony 全局迁移 | 子种群级，平均 2 子种群/次，随机对换 |
| Ring individual 全局迁移 | 个体级，1 最优个体/次，轮转选迁出 → 替换同子种群最差 |
| Random individual 全局迁移 | 个体级，平均 1 个体/次，随机 |

### 终止条件

- **DEA**: 与 HDEA 运行相同总代数
- **HDEA**: 完成 100 次全局迁移后终止

### 总代数计算

总代数 = 迁移间隔 × 20 × 100 = 迁移间隔 × 2000

### 各实例实测迁移间隔

| 实例 | 最优解 | 迁移间隔（5 档等差） | 对应总代数 |
|---|---|---|---|
| lin318 | 42,029 | 30k / 40k / 50k / 60k / 70k | 60M–140M |
| pcb442 | 50,778 | 100k / 150k / 200k / 250k / 300k | 200M–600M |
| u574 | 36,905 | 60k / 70k / 80k / 90k / 100k | 120M–200M |
| p654 | 34,643 | 50k / 60k / 70k / 80k / 90k | 100M–180M |
| d657 | 48,912 | 150k / 200k / 250k / 300k / 350k | 300M–700M |
| u724 | 41,910 | 150k / 200k / 250k / 300k / 350k | 300M–700M |
| rat783 | 8,806 | 150k / 200k / 250k / 300k / 350k | 300M–700M |
| dsj1000 | 18,659,688 | 600k / 800k / 1M / 1.2M / 1.4M | 1.2B–2.8B |
| pr1002 | 259,045 | 500k / 600k / 700k / 800k / 900k | 1B–1.8B |

### dEA 实验中的特殊设置

- dEA（非 HDEA）也使用相同总代数，迁移间隔为 `i=2` 代（每 2 代迁移一次）
- 每实例各算法独立运行 **30 次**，取平均值和标准差
- 显著性检验: **t-test with 95% confidence**
