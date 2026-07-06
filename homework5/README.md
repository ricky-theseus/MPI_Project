# MPI Allgather 实现 — homework5

## 题目

基于 `MPI_Send()` 和 `MPI_Recv()` 实现 `MPI_Allgather()`，遵循 HMPI 规范——进程向自身通信时直接内存拷贝。

## 解题思路

### MPI_Allgather 语义

每个进程贡献 sendcount 个元素，所有进程最终都得到完整的 size × recvcount 数据。等同于先 Gather 再 Bcast。

### 环形算法（Ring Allgather）

采用标准环形通信模式，共 size-1 步：

1. **初始化**：每个进程用 `memcpy` 将自己的数据拷贝到 recvbuf 的对应位置（HMPI）
2. **循环**：for step = 0 to size-2：
   - 发送当前持有的某个 chunk 给右邻 `(rank+1) % size`
   - 从左邻 `(rank-1+size) % size` 接收一个新的 chunk

每次循环，每个进程获得一块新的数据，size-1 步后 recvbuf 完整。

### 示例（4 进程）

```
步 0: P0→P1:chunk0, P0←P3:chunk3   P1→P2:chunk1, P1←P0:chunk0  ...
步 1: P0→P1:chunk3, P0←P3:chunk2   ...
步 2: P0→P1:chunk2, P0←P3:chunk1   ... → 全部完成
```

## 编译

```bash
MSBuild homework5.vcxproj /p:Configuration=Debug /p:Platform=x64
```

## 运行

```bash
mpiexec -n 4 homework5\x64\Debug\homework5.exe
```

## 示例输出

```
Rank 0: my    = 1 2 3 4
Rank 0: ref   = 1 2 3 4
Rank 0: PASS
...
```
