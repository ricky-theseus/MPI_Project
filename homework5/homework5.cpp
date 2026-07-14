#include <mpi.h>

#include <cstring>
#include <iostream>
#include <vector>

void my_MPI_Allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                      void* recvbuf, int recvcount, MPI_Datatype recvtype,
                      MPI_Comm comm) {
    int rank = 0, size = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int typesize = 0;
    MPI_Type_size(sendtype, &typesize);
    int blocksize = recvcount * typesize;

    std::memcpy(static_cast<char*>(recvbuf) + rank * blocksize,
                sendbuf, sendcount * typesize);

    for (int step = 0; step < size - 1; ++step) {
        int target = (rank + 1) % size;
        int source = (rank - 1 + size) % size;
        int send_chunk = (rank - step + size) % size;
        int recv_chunk = (rank - step - 1 + size) % size;

        MPI_Send(static_cast<char*>(recvbuf) + send_chunk * blocksize,
                 recvcount, recvtype, target, 0, comm);
        MPI_Recv(static_cast<char*>(recvbuf) + recv_chunk * blocksize,
                 recvcount, recvtype, source, 0, comm, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Each process holds a segment of the name+ID string
    constexpr int MAX_SEG = 20;
    const char* segments[] = {
        "姓名: 任天",
        "\xE7\xBF\x8E  学号",    // "翌  学号"
        ": 20231003",
        "512"
    };
    char myseg[MAX_SEG] = {};
    strncpy(myseg, segments[rank], MAX_SEG - 1);

    std::cout << "Rank " << rank << " holds: [" << myseg << "]\n";

    // Use my_MPI_Allgather to distribute all segments to every process
    std::vector<char> full(MAX_SEG * size);
    my_MPI_Allgather(myseg, MAX_SEG, MPI_CHAR,
                     full.data(), MAX_SEG, MPI_CHAR,
                     MPI_COMM_WORLD);

    // Every process can now print the full string
    std::cout << "Rank " << rank << " sees: ";
    for (int i = 0; i < size; ++i)
        std::cout << (full.data() + i * MAX_SEG);
    std::cout << "\n";

    // Validation with standard MPI_Allgather (int version)
    int sendval = rank + 1;
    std::vector<int> mybuf(size), expected(size);
    my_MPI_Allgather(&sendval, 1, MPI_INT,
                     mybuf.data(), 1, MPI_INT,
                     MPI_COMM_WORLD);
    MPI_Allgather(&sendval, 1, MPI_INT,
                  expected.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    bool ok = true;
    for (int i = 0; i < size; ++i) {
        if (mybuf[i] != expected[i]) { ok = false; break; }
    }
    std::cout << "Rank " << rank << " has all:";
    for (int v : mybuf) std::cout << ' ' << v;
    std::cout << "  => " << (ok ? "PASS\n" : "FAIL\n");

    MPI_Finalize();
    return 0;
}
