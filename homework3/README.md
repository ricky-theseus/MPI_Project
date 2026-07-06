# MPI 矩阵乘法 — homework3

## 解题思路

**问题**：计算两个 3×3 矩阵 A × B 的乘积，使用 MPI 并行加速。

**并行策略**：

1. **任务分解**：结果矩阵 C 有 3×3 = 9 个元素，每个元素 C[i][j] 的计算互相独立，因此启动 9 个 MPI 进程，每个进程负责计算一个元素。

2. **数据分发**：
   - Rank 0 持有完整的 A、B 矩阵
   - 通过两次 `MPI_Bcast` 将 A、B 广播给所有 9 个进程
   - 每个进程收到完整矩阵后，根据自身 rank 确定负责的元素位置：
     - 行号 `i = rank / 3`
     - 列号 `j = rank % 3`

3. **本地计算**：每个进程独立计算 `C[i][j] = Σ(k=0..2) A[i][k] × B[k][j]`

4. **结果归集**：通过 `MPI_Gather` 将所有 9 个计算结果收集回 rank 0，按顺序拼成完整的 C 矩阵

5. **输出**：Rank 0 打印 A、B、C 三个矩阵

**图示**：

```
Rank 0: C[0][0]  Rank 1: C[0][1]  Rank 2: C[0][2]
Rank 3: C[1][0]  Rank 4: C[1][1]  Rank 5: C[1][2]
Rank 6: C[2][0]  Rank 7: C[2][1]  Rank 8: C[2][2]

广播: rank 0 ──MPI_Bcast──→ rank 0~8
计算: rank r → C[r/3][r%3]
收集: rank 0~8 ──MPI_Gather──→ rank 0
```

## 编译

用 MSBuild（或 Visual Studio 打开 `.vcxproj` 后生成）：

```bash
MSBuild homework3.vcxproj /p:Configuration=Debug /p:Platform=x64
```

## 运行

从项目根目录运行：

```bash
mpiexec -n 9 homework3\x64\Debug\homework3.exe
```

或进入输出目录：

```bash
cd homework3\x64\Debug
mpiexec -n 9 homework3.exe
```

**注意**：`-n 9` 不可省略，必须恰好 9 个进程，否则程序会报错退出。

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
