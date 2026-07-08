# 并行归并排序 — homework2

## 题目

利用 MPI 实现并行归并排序：将数组分割到多个进程，每个进程本地排序后逐层合并。

## 解题思路

### 算法流程

1. **进程分割**：主进程读取数组 `m`、二分次数 `k` 和数组内容，通过 `MPI_Bcast` 广播
2. **数据分发**：用 `MPI_Scatterv` 将数组不均等分发给各进程
3. **本地排序**：每个进程对分到的子数组执行 `std::sort`
4. **逐层合并**：按 `for step = 1; step < active_size; step <<= 1` 构建归并树：
   - 偶数 rank 的进程从 `rank + step` 接收排序后的子数组并归并
   - 奇数 rank 的进程将自身子数组发送给 `rank - step` 后退出
5. 最终结果归集到 rank 0 输出

### 进程数优化

函数 `calculate_required_processes(m, k)` 根据数组长度和二分次数计算所需进程数，避免浪费。

## 编译

```bash
MSBuild homework2.vcxproj /p:Configuration=Debug /p:Platform=x64
```

## 运行

```bash
mpiexec -n <进程数> homework2\x64\Debug\homework2.exe
```

## 示例

```
请输入 m: 8
请输入 k: 3
请输入 8 个整数作为数组内容: 3 1 4 1 5 9 2 6
原始数组: 3 1 4 1 5 9 2 6
排序结果: 1 1 2 3 4 5 6 9
```
