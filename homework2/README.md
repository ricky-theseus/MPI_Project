# 并行归并排序 — homework2

## 题目

利用 MPI 实现并行归并排序：将数组分割到多个进程，每个进程本地排序后逐层合并。

## 解题思路

### 算法流程

1. **硬编码数据**：rank 0 持有原始数组，通过 `MPI_Scatterv` 不均等分发给各进程
2. **本地排序**：每个进程对分到的子数组执行 `std::sort`
3. **树形归并**：`for step = 1; step < active_size; step <<= 1`
   - 偶数 rank 从 `rank + step` 接收数据并归并
   - 奇数 rank 将自身数据发送给 `rank - step` 后退出
4. **验证**：最终结果与 `std::sort` 对比，输出 PASS/FAIL

### 示例（4 进程，12 个元素）

```
原始: [3 1 4 1 5 9 2 6 5 3 8 0]
  ↓ Scatterv 不均等分发（每个进程 3 个）
P0:[3 1 4]  P1:[1 5 9]  P2:[2 6 5]  P3:[3 8 0]
  ↓ std::sort
P0:[1 3 4]  P1:[1 5 9]  P2:[2 5 6]  P3:[0 3 8]
  ↓ step=1: P1→P0, P3→P2
P0:[1 1 3 4 5 9]         P2:[0 2 3 5 6 8]
  ↓ step=2: P2→P0
P0:[0 1 1 2 3 3 4 5 5 6 8 9]  ✓
```

## 优化说明

- 改用硬编码数据（原版需手动输入）
- 增加 `std::sort` 全数组对比验证
- 每个进程打印分块、本地排序、归并过程
- 输出全英文，避免编码问题

## 编译

```bash
MSBuild homework2.vcxproj /p:Configuration=Debug /p:Platform=x64
```

## 运行

```bash
mpiexec -n 4 homework2\x64\Debug\homework2.exe
```

## 示例输出

```
Rank 0 got chunk: 3 1 4
Rank 0 sorted: 1 3 4
Rank 0 merged from rank 1, now has: 1 1 3 4 5 9
Rank 0 merged from rank 2, now has: 0 1 1 2 3 3 4 5 5 6 8 9

Final sorted array: 0 1 1 2 3 3 4 5 5 6 8 9
PASS: matches std::sort
```
