# 课程设计：TSP 并行计算解决方案

## 文件说明

| 文件 | 说明 |
|---|---|
| `TSP0.C` | 串行 TSP 求解器（群体搜索 + 反转变异） |
| `TSP_MPI.cpp` | MPI 并行版 TSP 求解器 |
| `pcb442.tsp` | TSPLIB 标准测试集，442 个 PCB 钻孔城市坐标 |
| `hierarchical.pdf` | 参考文献 |
| `TSP0.exe` | 串行可执行文件 |

## 问题描述

旅行商问题（TSP）：给定 n 个城市及两两间的距离，求经过每个城市恰好一次的最短回路。

测试数据 `pcb442.tsp` 来自 TSPLIB，包含 442 个城市，坐标为 PCB 钻孔位置，距离取四舍五入后的整数欧氏距离。

## 串行版（TSP0.C）

### 算法

维护 `xColony` 条随机路径（排列），每轮对每条路径执行**反转变异**（翻转两城市间的一段），变异来源有概率随机或从另一条路径继承。每 `xColony` 轮执行一次父子锦标赛选择。

### 编译与运行

打开 **Developer Command Prompt for VS 2022**（或执行 `VsDevCmd.bat`），然后：

```cmd
cd course-project
cl /TC /O2 TSP0.C
TSP0.exe
```

`/TC` 强制以 C 语言编译。输入文件 `pcb442.tsp` 放在当前目录。标准输出打印每代最佳距离：

```
init success!!!
2:720740.000000
...
2000:340177.000000
```

详细日志追加到 `tsp0.txt`，每行格式：`<代数>  <耗时秒>  <距离>`

### 参数调优

编辑 `TSP0.C` 顶部变量：

| 变量 | 默认值 | 含义 |
|---|---|---|
| `xColony` | 100 | 群体大小（路径数量） |
| `probab1` | 0.02 | 随机变异概率 |
| `maxGen` | 200000 | 最大迭代代数 |

## MPI 并行版（TSP_MPI.cpp）

### 并行策略

将 `xColony` 条路径均匀分配到所有 MPI 进程，每代流程：

1. **root 广播**当前群体到所有进程（保证跨进程的变异参考正确性）
2. **各进程独立**对其分块的路径执行反转变异 + 距离计算
3. **gather** 子代结果回 root
4. **root 执行**锦标赛选择，输出最佳距离
5. 循环至 `maxGen` 代

### 编译（MSBuild，与其他 homework 一致方式）

#### CMD

```cmd
cd course-project
msbuild TSP_MPI.vcxproj /p:Configuration=Release /p:Platform=x64
```

#### PowerShell

```powershell
Set-Location course-project
msbuild TSP_MPI.vcxproj /p:Configuration=Release /p:Platform=x64
```

或直接双击 `TSP_MPI.vcxproj` 在 Visual Studio 中打开编译。

### 运行

#### CMD

```cmd
cd course-project\x64\Release
mpiexec -n 4 TSP_MPI.exe
```

#### PowerShell

```powershell
Set-Location course-project\x64\Release
mpiexec -n 4 TSP_MPI.exe
```

进程数建议整除 `xColony`（100），例如 2、4、5、10 等，使负载均衡。

### 与串行版对比

| 方面 | TSP0.C | TSP_MPI.cpp |
|---|---|---|
| 语言 | C | C++（含 MPI） |
| 随机种子 | `srand(time(NULL))` | `srand(time(NULL) + rank)` 各进程不同 |
| 计时 | `clock()` | `MPI_Wtime()` |
| 输出 | stdout + tsp0.txt | 同串行 |
| 加速方式 | 单线程 | 多进程并行评估 |
