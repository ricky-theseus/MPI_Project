# 并行归并排序 — homework2

## 题目

用 MPI 实现并行排序。

- 长度为 `m` 的一个大数组被二分 `k` 次，得到 `2^k` 个子数组
- 每个子数组内部排序后两两归并排序
- 进程数 `n` 由 `m` 与 `k` 共同决定：`n = min(2^k, m)`
- `m` 和 `k` 在运行时输入

## 解题思路

### 进程数计算

二分 k 次 → 最多需要 `2^k` 个进程。当 `2^k > m` 时，只用 `m` 个进程即可（每进程至少 1 个元素）。

### 算法流程

1. **输入**：rank 0 读取 m、k 和数组内容，`MPI_Bcast` 广播
2. **进程分割**：`MPI_Comm_split` 将冗余进程 idle，实际参与进程数 `active_size = min(n, world_size)`
3. **数据分发**：`MPI_Scatterv` 不均等分配
4. **本地排序**：每个进程 `std::sort`
5. **树形归并**：`for step = 1; step < active_size; step <<= 1`
   - 偶数 rank 接收并归并
   - 奇数 rank 发送后退出
6. **验证**：结果与 `std::sort` 完整排序对比

### 示例（m=8, k=3, 4 进程）

```
输入: [3 1 4 1 5 9 2 6], k=3 需要 8 进程，实际只开 4 进程
  ↓ Scatterv
P0:[3 1]  P1:[4 1]  P2:[5 9]  P3:[2 6]
  ↓ sort
P0:[1 3]  P1:[1 4]  P2:[5 9]  P3:[2 6]
  ↓ step=1: P1→P0, P3→P2
P0:[1 1 3 4]  P2:[2 5 6 9]
  ↓ step=2: P2→P0
P0:[1 1 2 3 4 5 6 9]  ✓
```

## 编译

```bash
MSBuild homework2.vcxproj /p:Configuration=Debug /p:Platform=x64
```

## 运行

```bash
mpiexec -n 4 homework2\x64\Debug\homework2.exe
```

## 示例运行

```bash
mpiexec -n 4 homework2\x64\Debug\homework2.exe
```

运行时依次输入：

```
Enter m (array length): 8
Enter k (binary split count): 2
Enter 8 integers: 3 1 4 1 5 9 2 6
```

另一种运行方式（用管道一次传入）：

```bash
"8 2 3 1 4 1 5 9 2 6" | mpiexec -n 4 homework2\x64\Debug\homework2.exe
```

### 示例输出（k=2, 4 进程）

```
m=8, k=2, needed processes=4, actual=4
Rank 0 got chunk: 3 1
Rank 0 sorted: 1 3
Rank 0 merged from rank 1, now has: 1 1 3 4
Rank 0 merged from rank 2, now has: 1 1 2 3 4 5 6 9

Final sorted array: 1 1 2 3 4 5 6 9
PASS: matches std::sort
```

### 进程数对应关系

| m | k | 2^k | 需要进程数 | 实际启动 | 每进程分几块 |
|---|---|---|---|---|---|
| 8 | 1 | 2 | 2 | 2 | 4 个 |
| 8 | 2 | 4 | 4 | 4 | 2 个 |
| 8 | 3 | 8 | 8 | 建议 `-n 8` | 1 个 |
| 8 | 4 | 16 | 8（上限 m） | 8 | 1 个 |
