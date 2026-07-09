# 算法 2：分层分布式进化算法 (HdEA)

## 并行架构

**双层层次结构 (hierarchical, coarse-grained)**

```
Level 2 (Master):   rank 0 精英种群 (ELITE=50)
Level 1 (Islands):  rank 0..N  独立子种群 (POP=50)
```

每代：各岛屿独立进化（0 MPI 通信）。
每 MIGRATE_INTERVAL=50 代：

```
1. Send/Recv:  rank i -> rank 0   (各岛最优个体)
2. Evolution:  rank 0              (精英种群进化 MASTER_GENS 代)
3. Broadcast:  rank 0 -> ALL       (全局最优精英回灌)
```

## 与算法 1 的差异

| 方面 | dEA | HdEA |
|------|-----|------|
| 通信频率 | 每代 | 每 50 代 |
| 数据粒度 | 全种群广播 | 仅最优个体 |
| 进程角色 | master/worker | 每个 rank 都是独立进化岛屿 |
| 群体独立性 | 共享全局种群 | 各自独立进化 |

## 编译与运行

```cmd
cl /EHsc /I"%MSMPI_INC%" TSP_HdEA.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
mpiexec -n 4 TSP_HdEA.exe [runs] [maxGen]
```

## 关键参数

- `POP=50`：每个岛屿的种群大小
- `ELITE=50`：master 精英种群大小
- `MIGRATE=50`：层次通信间隔（代）
- `MASTER_GENS=20`：master 每轮精英进化的代数
