# 并行计算课程设计 — 实现状态

## 目录结构 (course-project/)

```
04_dEA_island/
  TSP_dEA_island.cpp     # dEA 孤岛模型：Ring 拓扑 + Random-worst 迁移
  编译: cl /EHsc /I"%INC%" TSP_dEA_island.cpp /link msmpi.lib

05_HdEA/
  TSP_HdEA.cpp           # HdEA 4 变体：variant=0(RingIndiv) 1(RandIndiv) 2(RingColony) 3(RandColony)
  编译: cl /EHsc /I"%INC%" TSP_HdEA.cpp /link msmpi.lib

experiments/
  run_all.py             # 一键跑全部实验 (tiny/quick/paper 模式)
  analyze.py             # 生成 4 图 (收敛/箱线/t-test/耗时) + 报告
  output/                # 结果输出: <instance>/interval_<N>/<algorithm>/

tsp_instances/           # 9 个 TSP 实例数据
generate_tsp.py          # TSPLIB 下载/合成工具
```

## 使用

```bash
# 编译
cmd /c 04_dEA_island\compile.bat
cmd /c 05_HdEA\compile.bat

# 跑实验 (pcb442, 10轮, ~10分钟)
python experiments\run_all.py --instances pcb442 --runs 10 --intervals 100

# 分析出图
python experiments\analyze.py
```

## 参数设置 (Li & Hu 2014)

| 参数 | 值 |
|---|---|
| 子种群大小 | 100 |
| 变异率 (P) | 0.02, 线性衰减: P = 0.02×(1−g_n/g×0.01) |
| 交叉率 | 1−P |
| 映射率 | 0.05 |
| 临界速度 | 5000 |
| 局部迁移拓扑 | Ring |
| 局部迁移策略 | Random-worst |
| 局部:全局迁移比 | 20:1 |
| 终止条件 | 100 全局迁移轮次 |
| 总代数 | migInterval × 2000 |
