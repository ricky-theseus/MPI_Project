# 并行计算课程设计 — 实现计划

## 目标

基于 TSP0.C（Inver-over EA），用 MPI 实现三种分布式进化算法求解 TSP，使结果逐级提升。

## 文件结构

```
course-project/
├── TSP_base.h                # 基类 EA 头文件（三算法共用）
├── TSP_dEA.cpp               # Alg1: 主从式 dEA
├── TSP_HdEA.cpp              # Alg2: 岛屿 + 层次精英
├── TSP_HdEA_MP.cpp           # Alg3: Moving Colony
├── run_experiments.ps1       # 自动化运行
├── analyze_results.py        # 统计 + t 检验
├── plot_results.py           # 三张图
├── pcb442.tsp
│
├── results/
│   ├── serial/               # 原版 TSP0.C ×10 次
│   ├── dEA/                  # Alg1 ×10 次
│   ├── HdEA/                 # Alg2 ×10 次
│   └── HdEA_MP/              # Alg3 ×10 次
│
report/figures/
├── boxplot.png
├── barplot.png
└── convergence.png
```

## Phase 0 — TSP_base.h

从 TSP0.C 提取核心逻辑，增加两个缺失功能：

### 1. mapping 算子（PMX 风格）

```
选两个不同个体 A(worse), B(better)
A 中随机选段，B 中找到对应段
用 B 段替换 A 段，PMX 调整剩余位置
```

### 2. critical velocity 机制（阈值 = 5000）

```
每 2000 代计算改进速度
velocity < 5000 且 rand() < 0.05 时触发 mapping
```

### 3. 封装接口

```c
void tsp_init(const char *tspFile);                              // 初始化
void tsp_evolve_range(int colony[][CITY], double dis[],          // inver-over
    int popSize, int start, int end);
void tsp_select(int colony[][CITY], double dis[], int popSize);  // 选择
int  tsp_best_idx(double dis[], int popSize);                    // 最优个体下标
int  tsp_worst_idx(double dis[], int popSize);                   // 最差个体下标
double tsp_best_distance(double dis[], int popSize);             // 最优距离
void tsp_log_convergence(FILE *fp, int gen, double best);        // 收敛日志
```

## Phase 1 — TSP_dEA.cpp（主从式）

### 架构

```
rank 0 (Master): 种群 100 个体
rank 1-3 (Worker): 每进程处理 ~33 个体

每代:
  Bcast(种群) → Worker inver_over → Gatherv(子代) → Master select1
  每 2000 代记收敛曲线
```

### 参数

| 参数 | 值 |
|------|-----|
| xColony | 100 |
| maxGen | 200000 |
| probab1 | 0.02 |
| 进程数 | 4 |

### 预期

解质量 ≈ 串行 (~51600)，速度提升 2-3x

## Phase 2 — TSP_HdEA.cpp（岛屿 + 层次精英）

### 架构

```
4 个 rank，每个维护 1 个岛屿 (50 个体)

每代: 各岛独立 inver_over + select
每 20 代:
  Worker → Send(best) → Master → 选出全局 best → Bcast → 替换最差
```

### 参数

| 参数 | 值 |
|------|-----|
| 每岛 xColony | 50 |
| 总种群 | 200 |
| maxGen | 200000 |
| MIG_INTERVAL | 20 |

### 与现有代码差异

- 去掉 Master 端 MASTER_GENS 额外进化
- 增加 mapping + critical velocity

## Phase 3 — TSP_HdEA_MP.cpp（Moving Colony）

### 架构

在 HdEA 基础上，每进程 50 个体划分为：

```
[0..24] = 核心种群 (25，始终保留)
[25..49] = 移民种群 (25，参与 moving colony)

每 20 代: 层次精英交换（同 HdEA）
每 400 代: 额外执行 Moving Colony

  MPI_Sendrecv(&colony[25], sendTo, recvBuf, recvFrom)
  memcpy(&colony[25], recvBuf, ...)  // 替换本地移民
```

### 参数

| 参数 | 值 |
|------|-----|
| 每岛 xColony | 50 |
| 核心/移民 | 25/25 |
| MIG_INTERVAL | 20 |
| MOVING_INTERVAL | 400 |
| maxGen | 200000 |

## Phase 4 — 实验运行

### run_experiments.ps1

```
流程:
  1. 串行基线: TSP0.exe %Runs% → results/serial/
  2. mpiexec -n 4 TSP_dEA.exe %Runs% → results/dEA/
  3. mpiexec -n 4 TSP_HdEA.exe %Runs% → results/HdEA/
  4. mpiexec -n 4 TSP_HdEA_MP.exe %Runs% → results/HdEA_MP/
```

### 每算法输出

```
results/{algo}/results.txt:   run \t best \t time
results/{algo}/convergence/run_XX.txt:   gen \t time \t best
```

## Phase 5 — 分析与可视化

### analyze_results.py

```
输出:
  mean, std, min, max for each algorithm
  Paired t-test between all pairs
```

### plot_results.py — 三张图

| 图 | 文件名 |
|----|--------|
| 箱线图 | boxplot.png |
| 柱状图 | barplot.png |
| 收敛曲线 | convergence.png |

## 预期结果梯度

```
Serial (~51600) — 基线
  │ 种群相同，只是并行评估
  ▼
dEA (~51600) — 速度提升 2-3x，质量持平
  │ 总种群翻倍 100→200，岛屿隔离保持多样性
  ▼
HdEA (~51450) — 质量优于 dEA，t 检验显著
  │ Moving Colony 半种群批量换血
  ▼
HdEA-MP (~51200) — 质量优于 HdEA，t 检验显著
```
