# MPI Gather 实现 — homework4

## 题目

MPI 的复杂通信函数实际上就是用 6 个基础函数实现的。请基于 `MPI_Send()` 和 `MPI_Recv()` 实现 `MPI_Gather()`。

要求：遵循 HMPI（华为 MPI）规范——如果进程向自身收发消息（root 自己发给自己），不通过 `MPI_Send/Recv`，而是直接内存拷贝。

## 解题思路

### MPI_Gather 语义

- Root 进程收集所有进程的数据
- 每个进程（含 root）发送 sendcount 个元素
- Root 按进程 rank 顺序存入 recvbuf：进程 i 的数据放在 `recvbuf + i * recvcount * sizeof(type)`

### 实现策略

1. **Root 进程**：
   - 将自己的数据用 `memcpy` 直接拷贝到 recvbuf 的对应位置（HMPI 规范）
   - 循环接收其他所有进程的数据

2. **非 Root 进程**：
   - 直接 `MPI_Send` 将数据发送给 root

### 伪代码

```
if rank == root:
    memcpy(recvbuf, sendbuf, sendcount * typesize)
    for i in 0..size-1:
        if i == root: continue
        MPI_Recv(recvbuf + i * recvcount * typesize, ..., from i)
else:
    MPI_Send(sendbuf, sendcount, sendtype, to root)
```

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
my_MPI_Gather result:
 1 2 3 4
PASS: matches MPI_Gather
```
