#include <mpi.h>

#include <cstring>
#include <iostream>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

void my_MPI_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                   void* recvbuf, int recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm) {
    int rank = 0, size = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int type_size = 0;
    MPI_Type_size(sendtype, &type_size);

    if (rank == root) {
        std::memcpy(recvbuf, sendbuf, sendcount * type_size);

        for (int i = 0; i < size; ++i) {
            if (i == root) continue;
            MPI_Recv(static_cast<char*>(recvbuf) + i * recvcount * type_size,
                     recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
        }
    } else {
        MPI_Send(sendbuf, sendcount, sendtype, root, 0, comm);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    // Each process holds a segment of the name+ID string
    constexpr int MAX_SEG = 20;
    const char* segments[] = {
        "姓名: 任天",
        "\xE7\xBF\x8C  学号",    // "翌  学号"
        ": 20231003",
        "512"
    };
    char myseg[MAX_SEG] = {};
    strncpy(myseg, segments[rank], MAX_SEG - 1);

    std::cout << "Rank " << rank << " holds: [" << myseg << "]\n";

    // Use my_MPI_Gather to collect all segments at root
    std::vector<char> full(MAX_SEG * size);
    my_MPI_Gather(myseg, MAX_SEG, MPI_CHAR,
                  full.data(), MAX_SEG, MPI_CHAR,
                  0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n--- my_MPI_Gather: distributed concatenation ---\n";
        for (int i = 0; i < size; ++i)
            std::cout << (full.data() + i * MAX_SEG);
        std::cout << "\n";
    }

    // Validation with standard MPI_Gather (int version)
    int sendbuf = rank + 1;
    std::vector<int> recvbuf;
    if (rank == 0) recvbuf.resize(size);

    my_MPI_Gather(&sendbuf, 1, MPI_INT,
                  recvbuf.data(), 1, MPI_INT,
                  0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n--- standard MPI_Gather validation ---\n";
        std::vector<int> expected(size);
        MPI_Gather(&sendbuf, 1, MPI_INT,
                   expected.data(), 1, MPI_INT,
                   0, MPI_COMM_WORLD);
        bool ok = true;
        for (int i = 0; i < size; ++i) {
            if (recvbuf[i] != expected[i]) { ok = false; break; }
        }
        std::cout << "  my_MPI_Gather: ";
        for (int v : recvbuf) std::cout << v << " ";
        std::cout << "\n  MPI_Gather:    ";
        for (int v : expected) std::cout << v << " ";
        std::cout << "\n  => " << (ok ? "PASS\n" : "FAIL\n");
    } else {
        MPI_Gather(&sendbuf, 1, MPI_INT,
                   nullptr, 1, MPI_INT,
                   0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
