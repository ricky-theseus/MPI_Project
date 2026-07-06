#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {

	std::vector<int> merge_sorted_vectors(const std::vector<int>& left, const std::vector<int>& right) {
		std::vector<int> merged;
		merged.reserve(left.size() + right.size());

		auto left_it = left.begin();
		auto right_it = right.begin();
		while (left_it != left.end() && right_it != right.end()) {
			if (*left_it <= *right_it) {
				merged.push_back(*left_it++);
			}
			else {
				merged.push_back(*right_it++);
			}
		}

		merged.insert(merged.end(), left_it, left.end());
		merged.insert(merged.end(), right_it, right.end());
		return merged;
	}

	void print_vector(const std::vector<int>& values, const char* title) {
		std::cout << title;
		for (int value : values) {
			std::cout << ' ' << value;
		}
		std::cout << '\n';
	}

	size_t calculate_required_processes(size_t element_count, int split_times) {
		if (element_count == 0) {
			return 1;
		}

		size_t process_count = 1;
		for (int i = 0; i < split_times; ++i) {
			if (process_count >= element_count) {
				break;
			}
			if (process_count > element_count / 2) {
				process_count = element_count;
				break;
			}
			process_count *= 2;
		}

		return std::max<size_t>(1, std::min(process_count, element_count));
	}

}  // namespace

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size = 0;
	int world_rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	size_t m = 0;
	int k = 0;
	std::vector<int> input_values;

	if (world_rank == 0) {
		std::cout << "请输入 m: ";
		std::cin >> m;

		std::cout << "请输入 k: ";
		std::cin >> k;

		input_values.resize(m);
		std::cout << "请输入 " << m << " 个整数作为数组内容: ";
		for (size_t i = 0; i < m; ++i) {
			std::cin >> input_values[i];
		}
	}

	MPI_Bcast(&m, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (world_rank != 0) {
		input_values.resize(m);
	}
	if (m > 0) {
		MPI_Bcast(input_values.data(), static_cast<int>(m), MPI_INT, 0, MPI_COMM_WORLD);
	}

	if (m == 0) {
		if (world_rank == 0) {
			std::cout << "输入数组长度为 0，无需排序。\n";
		}
		MPI_Finalize();
		return 0;
	}

	const size_t required_processes = calculate_required_processes(m, k);
	const size_t active_processes = std::min<size_t>(required_processes, static_cast<size_t>(world_size));

	if (world_rank == 0) {
		std::cout << "数组长度 m = " << m << ", 二分次数 k = " << k << '\n';
		std::cout << "按算法需要的进程数 n = " << required_processes << '\n';
		std::cout << "当前实际启动进程数 = " << world_size << '\n';
		if (static_cast<size_t>(world_size) != required_processes) {
			std::cout << "提示：为完全按 k 次二分并行排序，建议使用 mpiexec -n " << required_processes << " 启动。\n";
		}
	}

	MPI_Comm active_comm = MPI_COMM_NULL;
	const int color = (static_cast<size_t>(world_rank) < active_processes) ? 0 : MPI_UNDEFINED;
	MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &active_comm);

	if (color == MPI_UNDEFINED) {
		MPI_Finalize();
		return 0;
	}

	int active_rank = 0;
	int active_size = 0;
	MPI_Comm_rank(active_comm, &active_rank);
	MPI_Comm_size(active_comm, &active_size);

	std::vector<int> global_values;
	std::vector<int> counts(active_size, 0);
	std::vector<int> displs(active_size, 0);

	const size_t base_count = m / static_cast<size_t>(active_size);
	const size_t remainder = m % static_cast<size_t>(active_size);
	for (int i = 0; i < active_size; ++i) {
		counts[i] = static_cast<int>(base_count + (static_cast<size_t>(i) < remainder ? 1 : 0));
		displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
	}

	if (active_rank == 0) {
		global_values = input_values;
		print_vector(global_values, "原始数组:");
	}

	std::vector<int> local_values(static_cast<size_t>(counts[active_rank]));
	MPI_Scatterv(active_rank == 0 ? global_values.data() : nullptr,
		counts.data(),
		displs.data(),
		MPI_INT,
		local_values.data(),
		counts[active_rank],
		MPI_INT,
		0,
		active_comm);

	std::sort(local_values.begin(), local_values.end());

	for (int step = 1; step < active_size; step <<= 1) {
		if (active_rank % (step * 2) == 0) {
			const int partner = active_rank + step;
			if (partner < active_size) {
				int received_count = 0;
				MPI_Recv(&received_count, 1, MPI_INT, partner, 0, active_comm, MPI_STATUS_IGNORE);
				std::vector<int> received_values(static_cast<size_t>(received_count));
				if (received_count > 0) {
					MPI_Recv(received_values.data(), received_count, MPI_INT, partner, 1, active_comm, MPI_STATUS_IGNORE);
				}
				local_values = merge_sorted_vectors(local_values, received_values);
			}
		}
		else if (active_rank % step == 0) {
			const int partner = active_rank - step;
			const int send_count = static_cast<int>(local_values.size());
			MPI_Send(&send_count, 1, MPI_INT, partner, 0, active_comm);
			if (send_count > 0) {
				MPI_Send(local_values.data(), send_count, MPI_INT, partner, 1, active_comm);
			}
			break;
		}
	}

	if (active_rank == 0) {
		print_vector(local_values, "排序结果:");
	}

	MPI_Comm_free(&active_comm);
	MPI_Finalize();
	return 0;
}
