#include <fstream>
#include <iostream>

#define ASSEMBLE_REDUCTION_MATRIX
#define INDICATE_PROGRESS
#define SKIP_APPARENT_PAIRS
#define USE_ARRAY_HASHMAP
#define USE_CELLS_WITHOUT_DIMENSION
#define SORT_COLUMNS_BY_PIVOT
#define COUNT_ONLY
// #define WITH_HDF5
// #define USE_COEFFICIENTS
// #define MANY_VERTICES

#include "../include/argparser.h"
#include "../include/parameters.h"
#include "../include/persistence.h"

#ifdef WITH_HDF5
#include "../include/output/hdf5.h"
#endif

//
// Compute directed flag complex homology
//

#include "../include/usage/flagser-count.h"

std::vector<size_t> count_cells(filtered_directed_graph_t& graph, const flagser_parameters& params) {
	// Aggregated counts
	std::vector<size_t> total_cell_count;

	int64_t total_euler_characteristic = 0;
	size_t total_max_dim = 0;

#ifdef WITH_HDF5
	hdf5_output_t* output = nullptr;
	output = new hdf5_output_t(params);
#endif

	std::vector<filtered_directed_graph_t> subgraphs{graph};
	if (params.split_into_connected_components) { subgraphs = graph.get_connected_subgraphs(2); }

	bool is_first_line = true;
	for (auto subgraph : subgraphs) {
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "computing the flag complex" << std::flush << "\r";
#endif

		directed_flag_complex_t complex(subgraph);
		struct cell_counter_t {
#ifdef WITH_HDF5
		private:
			hdf5_output_t* output;

		public:
			cell_counter_t(hdf5_output_t* _output = nullptr) : output(_output) {}
#endif
			void done() {}
			void operator()(vertex_index_t*
#ifdef WITH_HDF5
			                    first_vertex
#endif
			                ,
			                int size) {
				// Add (-1)^size to the Euler characteristic
				if (size & 1)
					ec++;
				else
					ec--;

				if (cell_counts.size() < size_t(size)) { cell_counts.resize(size, 0); }
				cell_counts[size - 1]++;

#ifdef WITH_HDF5
				if (output != nullptr) output->write_cell(first_vertex, size);
#endif
			}

			int64_t euler_characteristic() const { return ec; }
			std::vector<size_t> cell_count() const { return cell_counts; }

		private:
			int64_t ec = 0;
			std::vector<size_t> cell_counts;
		};

		std::vector<cell_counter_t> cell_counter(params.nb_threads);
		for (size_t i = 0; i < params.nb_threads; i++)
			cell_counter[i] = cell_counter_t(
#ifdef WITH_HDF5
			    output
#endif
			);

#ifdef WITH_HDF5
		if (output != nullptr)
			complex.for_each_cell(*cell_counter.data(), 0, 10000);
		else {
#endif
			complex.for_each_cell(cell_counter, 0, 10000);
#ifdef WITH_HDF5
		}
#endif
		int64_t euler_characteristic = 0;
		for (size_t i = 0; i < params.nb_threads; i++) euler_characteristic += cell_counter[i].euler_characteristic();

#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif

		if (is_first_line) std::cout << "# [euler_characteristic cell_count_dim_0 cell_count_dim_1 ...]" << std::endl;
		std::cout << euler_characteristic;

		std::vector<std::vector<size_t>> cell_counts(params.nb_threads);
		size_t max_dim = 0;
		for (size_t i = 0; i < params.nb_threads; i++) {
			cell_counts[i] = cell_counter[i].cell_count();
			size_t dim = cell_counts[i].size();
			max_dim = max_dim < dim ? dim : max_dim;
		}

		total_cell_count.resize(max_dim, 0);
		for (size_t dim = 0; dim < max_dim; dim++) {
			size_t size = 0;
			for (size_t i = 0; i < params.nb_threads; i++)
				size += cell_counts[i].size() > dim ? cell_counts[i][dim] : 0;
			std::cout << " " << size;
			total_cell_count[dim] += size;
		}
		std::cout << std::endl;
		total_euler_characteristic += euler_characteristic;
		total_max_dim = std::max(total_max_dim, max_dim);

		is_first_line = false;
	}

#ifdef WITH_HDF5
	if (output != nullptr) delete output;
#endif

	if (params.split_into_connected_components) {
		std::cout << std::endl << "# Total" << std::endl;
		std::cout << total_euler_characteristic;
		for (size_t dim = 0; dim < total_max_dim; dim++) std::cout << " " << total_cell_count[dim];
		std::cout << std::endl;
	}

	return total_cell_count;
}

int main(int argc, char** argv) {
	try {
		auto arguments = parse_arguments(argc, argv);

		auto positional_arguments = get_positional_arguments(arguments);
		auto named_arguments = get_named_arguments(arguments);
		auto params = flagser_parameters(named_arguments);
		named_arguments_t::const_iterator it;
		if (named_arguments.find("help") != named_arguments.end()) { print_usage_and_exit(-1); }

		if (positional_arguments.size() == 0) { print_usage_and_exit(-1); }
		const char* input_filename = positional_arguments[0];

		filtered_directed_graph_t graph = read_filtered_directed_graph(input_filename, params);

		auto cell_count = count_cells(graph, params);
	} catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}
