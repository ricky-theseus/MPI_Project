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

    int sendval = rank + 1;

    std::vector<int> mybuf;
    mybuf.resize(size);

    my_MPI_Allgather(&sendval, 1, MPI_INT,
                     mybuf.data(), 1, MPI_INT,
                     MPI_COMM_WORLD);

    std::vector<int> expected;
    expected.resize(size);
    MPI_Allgather(&sendval, 1, MPI_INT,
                  expected.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    bool ok = true;
    for (int i = 0; i < size; ++i) {
        if (mybuf[i] != expected[i]) {
            ok = false;
            break;
        }
    }

    std::cout << "Rank " << rank << ": my    =";
    for (int v : mybuf) std::cout << ' ' << v;
    std::cout << "\nRank " << rank << ": ref   =";
    for (int v : expected) std::cout << ' ' << v;
    std::cout << "\nRank " << rank << ": " << (ok ? "PASS" : "FAIL") << "\n\n";

    MPI_Finalize();
    return 0;
}
