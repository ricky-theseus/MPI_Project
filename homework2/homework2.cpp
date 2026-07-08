#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <vector>

namespace {

std::vector<int> merge_sorted(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> merged;
    merged.reserve(a.size() + b.size());
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size())
        merged.push_back(a[i] <= b[j] ? a[i++] : b[j++]);
    while (i < a.size()) merged.push_back(a[i++]);
    while (j < b.size()) merged.push_back(b[j++]);
    return merged;
}

int calc_processes(int m, int k) {
    int n = 1;
    for (int i = 0; i < k; ++i) {
        if (n >= m) break;
        n *= 2;
    }
    return std::min(n, m);
}

}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m = 0, k = 0;
    std::vector<int> data;

    if (rank == 0) {
        std::cout << "Enter m (array length): ";
        std::cin >> m;
        std::cout << "Enter k (binary split count): ";
        std::cin >> k;
        data.resize(m);
        std::cout << "Enter " << m << " integers: ";
        for (int i = 0; i < m; ++i) std::cin >> data[i];
    }

    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) data.resize(m);
    if (m > 0) MPI_Bcast(data.data(), m, MPI_INT, 0, MPI_COMM_WORLD);

    int n = calc_processes(m, k);
    int active_size = std::min(n, size);

    if (rank == 0) {
        std::cout << "m=" << m << ", k=" << k << ", needed processes=" << n
                  << ", actual=" << active_size << "\n";
    }

    if (active_size > m) {
        if (rank == 0) std::cerr << "Error: too many processes\n";
        MPI_Finalize();
        return 1;
    }

    MPI_Comm active_comm = MPI_COMM_NULL;
    int color = (rank < active_size) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &active_comm);

    if (color == MPI_UNDEFINED) {
        MPI_Finalize();
        return 0;
    }

    int active_rank = 0;
    MPI_Comm_rank(active_comm, &active_rank);
    MPI_Comm_size(active_comm, &active_size);

    std::vector<int> sendcounts(active_size);
    std::vector<int> displs(active_size);
    int base = m / active_size;
    int rem = m % active_size;
    int offset = 0;
    for (int i = 0; i < active_size; ++i) {
        sendcounts[i] = base + (i < rem ? 1 : 0);
        displs[i] = offset;
        offset += sendcounts[i];
    }

    int local_count = sendcounts[active_rank];
    std::vector<int> local(local_count);
    MPI_Scatterv(data.data(), sendcounts.data(), displs.data(), MPI_INT,
                 local.data(), local_count, MPI_INT, 0, active_comm);

    std::cout << "Rank " << active_rank << " got chunk:";
    for (int v : local) std::cout << ' ' << v;
    std::cout << "\n";

    std::sort(local.begin(), local.end());

    std::cout << "Rank " << active_rank << " sorted:";
    for (int v : local) std::cout << ' ' << v;
    std::cout << "\n";

    for (int step = 1; step < active_size; step <<= 1) {
        if (active_rank % (step * 2) == 0) {
            int partner = active_rank + step;
            if (partner < active_size) {
                int recv_count = 0;
                MPI_Recv(&recv_count, 1, MPI_INT, partner, 0, active_comm, MPI_STATUS_IGNORE);
                std::vector<int> buf(recv_count);
                if (recv_count > 0)
                    MPI_Recv(buf.data(), recv_count, MPI_INT, partner, 1, active_comm, MPI_STATUS_IGNORE);
                local = merge_sorted(local, buf);
                std::cout << "Rank " << active_rank << " merged from rank " << partner
                          << ", now has:"; for (int v : local) std::cout << ' ' << v;
                std::cout << "\n";
            }
        } else if (active_rank % step == 0) {
            int partner = active_rank - step;
            int send_count = static_cast<int>(local.size());
            MPI_Send(&send_count, 1, MPI_INT, partner, 0, active_comm);
            if (send_count > 0)
                MPI_Send(local.data(), send_count, MPI_INT, partner, 1, active_comm);
            std::cout << "Rank " << active_rank << " sent data to rank " << partner << ", exiting\n";
            break;
        }
    }

    if (active_rank == 0) {
        std::cout << "\nFinal sorted array:";
        for (int v : local) std::cout << ' ' << v;
        std::cout << "\n";

        std::vector<int> expected = data;
        std::sort(expected.begin(), expected.end());
        bool ok = (local == expected);
        std::cout << (ok ? "PASS: matches std::sort\n" : "FAIL: mismatch\n");
    }

    MPI_Comm_free(&active_comm);
    MPI_Finalize();
    return 0;
}
