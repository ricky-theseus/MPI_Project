#include <mpi.h>

#include <cstring>
#include <iostream>
#include <vector>

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

    constexpr int N = 4;
    int sendbuf = rank + 1;

    std::cout << "Rank " << rank << " sendbuf = " << sendbuf << "\n";

    std::vector<int> recvbuf;
    if (rank == 0) {
        recvbuf.resize(size);
    }

    my_MPI_Gather(&sendbuf, 1, MPI_INT,
                  recvbuf.data(), 1, MPI_INT,
                  0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n--- my_MPI_Gather result ---\n";
        for (int i = 0; i < size; ++i) {
            std::cout << "  from rank " << i << " -> recvbuf[" << i << "] = " << recvbuf[i] << "\n";
        }

        std::vector<int> expected;
        expected.resize(size);
        MPI_Gather(&sendbuf, 1, MPI_INT,
                   expected.data(), 1, MPI_INT,
                   0, MPI_COMM_WORLD);

        bool ok = true;
        for (int i = 0; i < size; ++i) {
            if (recvbuf[i] != expected[i]) {
                ok = false;
                break;
            }
        }
        std::cout << "\n--- Comparison ---\n";
        std::cout << "  my_MPI_Gather: ";
        for (int v : recvbuf) std::cout << v << " ";
        std::cout << "\n  MPI_Gather:    ";
        for (int v : expected) std::cout << v << " ";
        std::cout << "\n  => " << (ok ? "PASS: results match\n" : "FAIL: mismatch\n");
    } else {
        MPI_Gather(&sendbuf, 1, MPI_INT,
                   nullptr, 1, MPI_INT,
                   0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
