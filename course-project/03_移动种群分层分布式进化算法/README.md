# 算法 3：移动种群分层分布式进化算法 (HdEA-MP)

## 并行架构

**层次结构 + 环状迁移**

在算法 2 的层次结构上增加**同级岛屿间的环状迁移**：

```
每 MIGRATE_INTERVAL=50 代:
  Step 1: 层次通信 (同 HdEA)
    rank i -> rank 0:  send best
    rank 0 -> ALL:     broadcast global best (hierarchy)
    
  Step 2: 环状迁移 (新增)
    rank i -> rank (i+1)%N:  send M=3 best individuals
    rank i <- rank (i-1)%N:  recv M=3 individuals
    replace worst in island with migrants
```

## 与算法 2 的差异

| 方面 | HdEA | HdEA-MP |
|------|------|----------|
| 通信拓扑 | 星型（全部到 master） | 星型 + 环型 |
| 个体流动 | 仅上下层流动 | 上下层 + 同级流动 |
| 多样性维持 | 仅靠种群自身变异 | 迁移注入外部基因 |
| `sort_by_fitness` | — | √ |

## 关键参数

- `MIG_M=3`：每次迁移发送/接收的个体数
- 其他参数同 HdEA

## 编译与运行

```cmd
cl /EHsc /I"%MSMPI_INC%" TSP_HdEA_MP.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
mpiexec -n 4 TSP_HdEA_MP.exe [runs] [maxGen]
```
