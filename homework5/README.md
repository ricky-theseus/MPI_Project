# MPI Allgather 实现 — homework5

## 题目

基于 `MPI_Send()` 和 `MPI_Recv()` 实现 `MPI_Allgather()`，遵循 HMPI 规范——进程向自身通信时直接内存拷贝。

## 分布式打印演示

将字符串 `"姓名: 任天翌  学号: 20231003512"` 切分为 4 段，分发给 4 个进程：

| 进程 | 持有的段 |
|---|---|
| Rank 0 | "姓名: 任天" |
| Rank 1 | "翌  学号" |
| Rank 2 | ": 20231003" |
| Rank 3 | "512" |

调用 `my_MPI_Allgather` 将所有段分发到每一个进程，每个进程各自拼接并打印完整字符串。这展示了 Allgather 与 Gather 的区别：**所有进程**都获得了完整数据，而非仅 root 进程。

## 解题思路

### MPI_Allgather 语义

每个进程贡献 sendcount 个元素，所有进程最终都得到完整的 size × recvcount 数据。等同于先 Gather 再 Bcast。

### 环形算法（Ring Allgather）

采用标准环形通信模式，共 size-1 步：

1. **初始化**：每个进程用 `memcpy` 将自己的数据拷贝到 recvbuf 的对应位置（HMPI）
2. **循环**：for step = 0 to size-2：发送当前持有的某个 chunk 给右邻，从左邻接收一个新的 chunk

每次循环，每个进程获得一块新的数据，size-1 步后 recvbuf 完整。

## 编译

```bash
cl /EHsc /I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" homework5.cpp ^
   /Fe:homework5.exe ^
   /link /LIBPATH:"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" msmpi.lib
```

## 运行

```bash
mpiexec -n 4 homework5.exe
```

## 示例输出

```
Rank 0 holds: [姓名: 任天]
Rank 1 holds: [翌  学号]
Rank 2 holds: [: 20231003]
Rank 3 holds: [512]
Rank 0 sees: 姓名: 任天翌  学号: 20231003512
Rank 1 sees: 姓名: 任天翌  学号: 20231003512
Rank 2 sees: 姓名: 任天翌  学号: 20231003512
Rank 3 sees: 姓名: 任天翌  学号: 20231003512
Rank 0 has all: 1 2 3 4  => PASS
Rank 1 has all: 1 2 3 4  => PASS
Rank 2 has all: 1 2 3 4  => PASS
Rank 3 has all: 1 2 3 4  => PASS
```
