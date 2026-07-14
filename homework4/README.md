# MPI Gather 实现 — homework4

## 题目

基于 `MPI_Send()` 和 `MPI_Recv()` 实现 `MPI_Gather()`，遵循 HMPI 规范——进程向自身通信时直接内存拷贝。

## 分布式打印演示

将字符串 `"姓名: 任天翌  学号: 20231003512"` 切分为 4 段，分发给 4 个进程：

| 进程 | 持有的段 |
|---|---|
| Rank 0 | "姓名: 任天" |
| Rank 1 | "翌  学号" |
| Rank 2 | ": 20231003" |
| Rank 3 | "512" |

调用 `my_MPI_Gather` 将所有段汇聚至 Rank 0，Rank 0 拼接后输出完整字符串。这展示了 Gather 的核心语义：将分散在各进程的数据按 rank 顺序汇聚到一处。

## 解题思路

### MPI_Gather 语义

- Root 进程收集所有进程的数据
- 每个进程（含 root）发送 sendcount 个元素
- Root 按进程 rank 顺序存入 recvbuf

### 实现策略

1. **Root 进程**：将自己的数据用 `memcpy` 直接拷贝到 recvbuf 的对应位置（HMPI 规范），然后循环接收其他进程的数据
2. **非 Root 进程**：直接 `MPI_Send` 将数据发送给 root

## 编译

```bash
MSBuild homework4.vcxproj /p:Configuration=Debug /p:Platform=x64
```

## 运行

```bash
mpiexec -n 4 homework4\x64\Debug\homework4.exe
```

## 示例输出

```
Rank 0 holds: [姓名: 任天]
Rank 1 holds: [翌  学号]
Rank 2 holds: [: 20231003]
Rank 3 holds: [512]

--- my_MPI_Gather: distributed concatenation ---
姓名: 任天翌  学号: 20231003512

--- standard MPI_Gather validation ---
  my_MPI_Gather: 1 2 3 4
  MPI_Gather:    1 2 3 4
  => PASS
```
