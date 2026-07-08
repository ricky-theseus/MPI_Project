# MPI 矩阵乘法 — homework3

## 题目

利用 MPI 实现 3×3 矩阵乘法 `C = A × B`，每个进程计算结果矩阵中的一个元素。

## 解题思路

### 算法设计

采用 9 个进程（3×3 网格），进程 rank `i * N + j` 负责计算 `C[i][j]`：

1. **初始化**：rank 0 硬编码矩阵 A 和 B
2. **广播**：`MPI_Bcast` 将 A、B 发送给所有进程
3. **局部计算**：每个进程计算 `C[i][j] = Σ A[i][k] * B[k][j]`，其中 `i = rank / N`，`j = rank % N`
4. **收集结果**：`MPI_Gather` 将各进程的结果汇集到 rank 0 并输出

### 矩阵

```
    A            B            C
1  2  3      9  8  7     30  24  18
4  5  6  ×   6  5  4  =  84  69  54
7  8  9      3  2  1    138 114  90
```

## 编译

```bash
MSBuild homework3.vcxproj /p:Configuration=Debug /p:Platform=x64
```

## 运行

```bash
mpiexec -n 9 homework3\x64\Debug\homework3.exe
```

## 示例输出

```
Matrix A:
   1   2   3
   4   5   6
   7   8   9

Matrix B:
   9   8   7
   6   5   4
   3   2   1

C = A x B:
  30  24  18
  84  69  54
 138 114  90
```
