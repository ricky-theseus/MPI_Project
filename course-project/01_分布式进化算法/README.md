# 算法 1：分布式进化算法 (dEA)

## 并行架构

**Master-Worker 模型 (fine-grained)**

主进程 (rank 0) 维护全局种群 `colony[POP][CITY]`。每代：

```
broadcast:  rank 0 -> ALL   (全局种群 + 适应度)
mutation:   each rank       (mutate its chunk)
gather:     ALL -> rank 0   (子代路径 + 适应度)
selection:  rank 0          (锦标赛选择)
```

## 与串行算法的差异

- 与串行逻辑**一致**：每代遍历所有个体，执行反转变异，再从父代和子代中择优
- `srand(time(NULL) + rank)` 使各进程获得不同随机序列
- 使用 `MPI_Gatherv` 灵活收集不等长子代数据

## 编译与运行

```cmd
cl /EHsc /I"%MSMPI_INC%" TSP_dEA.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
mpiexec -n 4 TSP_dEA.exe [runs] [maxGen]
```

## 关键代码分析

- `probab1 = 0.02`：以 2% 概率随机选变异起点，98% 概率从另一路径继承城市
- `invert()`：翻转城市序列中的一段，实现 2-opt 等效的反转变异
- `position()`：在排列中查找给定城市的索引
