#include <mpi.h>

#include <iomanip>
#include <iostream>
#include <vector>

constexpr int N = 3;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank = 0, world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size != N * N) {
        if (world_rank == 0) {
            std::cerr << "错误：需要 " << N * N << " 个进程，当前 " << world_size << " 个\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::vector<int> A, B;
    if (world_rank == 0) {
        A = {
            1, 2, 3,
            4, 5, 6,
            7, 8, 9
        };
        B = {
            9, 8, 7,
            6, 5, 4,
            3, 2, 1
        };
    }

    A.resize(N * N);
    B.resize(N * N);
    MPI_Bcast(A.data(), N * N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B.data(), N * N, MPI_INT, 0, MPI_COMM_WORLD);

    const int i = world_rank / N;
    const int j = world_rank % N;

    int local_result = 0;
    for (int k = 0; k < N; ++k) {
        local_result += A[i * N + k] * B[k * N + j];
    }

    std::vector<int> results;
    if (world_rank == 0) {
        results.resize(N * N);
    }
    MPI_Gather(&local_result, 1, MPI_INT, results.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        std::cout << "矩阵 A:\n";
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << std::setw(4) << A[i * N + j];
            }
            std::cout << '\n';
        }

        std::cout << "\n矩阵 B:\n";
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << std::setw(4) << B[i * N + j];
            }
            std::cout << '\n';
        }

        std::cout << "\nC = A x B:\n";
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << std::setw(4) << results[i * N + j];
            }
            std::cout << '\n';
        }
    }

    MPI_Finalize();
    return 0;
}
