# 课程设计：TSP 并行进化算法（三种 MPI 算法实现与对比）

## 目录结构

| 目录 | 内容 |
|------|------|
| `01_分布式进化算法/` | **算法 1**：主从式分布式进化算法 (dEA) |
| `02_分层分布式进化算法/` | **算法 2**：分层分布式进化算法 (HdEA) |
| `03_移动种群分层分布式进化算法/` | **算法 3**：带迁移的分层分布式进化算法 (HdEA-MP) |
| `TSP0.C` | 串行参考实现 |
| `TSP_MPI.cpp` | MPI 并行版（算法 1 原型） |
| `pcb442.tsp` | TSPLIB 标准测试集，442 个城市 |
| `hierarchical.pdf` | 参考文献 |

每个算法目录下包含：
- 独立可编译的 `.cpp` 源码
- 已编译的 `.exe` 可执行文件
- 该算法的说明 README

## 问题描述

旅行商问题（TSP）：给定 n 个城市及两两间的距离，求经过每个城市恰好一次的最短回路。使用 `pcb442.tsp`（442 个城市，整数欧氏距离）。

## 实验环境

- **CPU**: AMD Ryzen 7 5800X (8C/16T)
- **RAM**: 32 GB
- **OS**: Windows 10
- **MPI**: Microsoft MPI v10.1
- **Compiler**: Visual Studio 2022 (MSVC 19.44) x64, Release 模式

## 三种算法架构

### 算法 1：分布式进化算法 (dEA)

**主从式（Master-Worker）细粒度并行。**

- rank 0（master）维护整个全局种群
- 每代广播完整种群到所有 worker
- 每个 worker 对其分配的路径子集执行反转变异
- 收集子代结果回 master
- master 执行锦标赛选择

**特点**：通信密集（每代 2 次全广播），全局同步，逻辑与串行版最接近。

### 算法 2：分层分布式进化算法 (HdEA)

**岛屿模型 + 双层层次结构。**

- **下层**：每个 rank 独立维护一个子种群（岛屿），独立进化
- **上层**：rank 0 维护精英种群，每 MIGRATE_INTERVAL 代收集各岛屿的最优个体，额外进化后广播全局最优回所有岛屿

**特点**：通信稀疏（每 MIGRATE_INTERVAL=50 代通信一次），适合大规模集群。

### 算法 3：移动种群分层分布式进化算法 (HdEA-MP)

在算法 2 的基础上增加**环状迁移**：

- 继承算法 2 的层次结构
- 每 MIGRATE_INTERVAL 代，除层次通信外，各岛屿将其 M=3 个最优个体发送给下一个 rank（环状拓扑）
- 接收到的迁移个体替换接收方的最差个体

**特点**：层次结构保证全局收敛压力，环状迁移增强种群多样性。

## 编译与运行

### 编译（使用 VS 2022 开发者命令提示符）

```cmd
cd course-project

:: 算法 1
cl /EHsc /I"%MSMPI_INC%" 01_分布式进化算法\TSP_dEA.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib

:: 算法 2
cl /EHsc /I"%MSMPI_INC%" 02_分层分布式进化算法\TSP_HdEA.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib

:: 算法 3
cl /EHsc /I"%MSMPI_INC%" 03_移动种群分层分布式进化算法\TSP_HdEA_MP.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
```

### 运行

```cmd
cd course-project
mpiexec -n 4 01_分布式进化算法\TSP_dEA.exe [runs] [maxGen]
mpiexec -n 4 02_分层分布式进化算法\TSP_HdEA.exe [runs] [maxGen]
mpiexec -n 4 03_移动种群分层分布式进化算法\TSP_HdEA_MP.exe [runs] [maxGen]
```

## 实验结果

> TODO: 在 10 次运行后填入

| 算法 | mean | std | min | max |
|------|------|-----|-----|-----|
| dEA | - | - | - | - |
| HdEA | - | - | - | - |
| HdEA-MP | - | - | - | - |

### T 检验（p < 0.05）

| 对比组 | t | p | 结论 |
|--------|---|---|------|
| dEA vs HdEA | - | - | - |
| dEA vs HdEA-MP | - | - | - |
| HdEA vs HdEA-MP | - | - | - |

### 分析

> 实验结果分析待补充。

## 参数配置

| 参数 | 默认值 | 说明 |
|------|--------|------|
| 种群大小 | 50/岛, 50/精英 | 算法 1: 100 全局; 算法 2/3: 50/岛 × 4 = 200 总 |
| 变异概率 | 0.02 | 随机变异 vs 交叉继承 |
| 迁移间隔 | 50 代 | 层次/环状通信频率 |
| 精英进化代数 | 20 | master 每轮进化精英种群的代数 |
| 迁移个体数 | 3 | 算法 3 的环状迁移中每次发送的个体数 |
