#include <fstream>
#include <iostream>

#define ASSEMBLE_REDUCTION_MATRIX
#define INDICATE_PROGRESS
#define SKIP_APPARENT_PAIRS
#define USE_ARRAY_HASHMAP
#define USE_CELLS_WITHOUT_DIMENSION
#define SORT_COLUMNS_BY_PIVOT
// #define WITH_HDF5
// #define KEEP_FLAG_COMPLEX_IN_MEMORY
// #define USE_COEFFICIENTS
// #define MANY_VERTICES

#include "../include/parameters.h"
#include "../include/persistence.h"

//
// Compute directed flag complex homology
//

#ifdef KEEP_FLAG_COMPLEX_IN_MEMORY
#include "../include/complex/directed_flag_complex_in_memory_computer.h"
#else
#include "../include/complex/directed_flag_complex_computer.h"
#endif

#include "../include/usage/flagser.h"

#ifdef KEEP_FLAG_COMPLEX_IN_MEMORY
typedef directed_flag_complex_in_memory_computer_t directed_flag_complex_compute_t;
#else
typedef directed_flag_complex_computer_t directed_flag_complex_compute_t;
#endif

#ifdef RETRIEVE_PERSISTENCE
std::vector<persistence_computer_t<directed_flag_complex_compute_t>>
#else
void
#endif
compute_homology(filtered_directed_graph_t& graph, const flagser_parameters& params) {

	std::vector<filtered_directed_graph_t> subgraphs{graph};
	if (params.split_into_connected_components) { subgraphs = graph.get_connected_subgraphs(2); }

#ifdef RETRIEVE_PERSISTENCE
	std::vector<persistence_computer_t<directed_flag_complex_compute_t>> complex_subgraphs;
#endif

	auto output = get_output<directed_flag_complex_compute_t>(params);
	size_t component_number = 1;
	for (auto subgraph : subgraphs) {
		directed_flag_complex_compute_t complex(subgraph, params);

		output->set_complex(&complex);
		if (params.split_into_connected_components) {
			if (component_number > 1) output->print("\n");
			output->print("## Path component number ");
			output->print(std::to_string(component_number));
			output->print("\n");

#ifdef INDICATE_PROGRESS
			std::cout << "\033[K";
#endif
			if (component_number > 1) std::cout << "\n";
			std::cout << "# Path component number " << component_number << std::endl;

			component_number++;
		}

#ifdef RETRIEVE_PERSISTENCE
		complex_subgraphs.push_back(
		    persistence_computer_t<decltype(complex)>(complex, output.get(), params.max_entries, params.modulus));
		complex_subgraphs.back().compute_persistence(params.min_dimension, params.max_dimension);
#else
		persistence_computer_t<decltype(complex)> persistence_computer(complex, output.get(), params.max_entries,
		                                                               params.modulus);
		persistence_computer.compute_persistence(params.min_dimension, params.max_dimension);
#endif
	}

	if (params.split_into_connected_components) { output->print("\n## Total\n"); }

	output->print_aggregated_results();

#ifdef RETRIEVE_PERSISTENCE
	return complex_subgraphs;
#endif
}

int main(int argc, char** argv) {
	auto arguments = parse_arguments(argc, argv);

	auto positional_arguments = get_positional_arguments(arguments);
	auto named_arguments = get_named_arguments(arguments);
	auto params = flagser_parameters(named_arguments);

	if (named_arguments.find("help") != named_arguments.end()) { print_usage_and_exit(-1); }

	if (positional_arguments.size() == 0) { print_usage_and_exit(-1); }
	const char* input_filename = positional_arguments[0];

	filtered_directed_graph_t graph = read_filtered_directed_graph(input_filename, named_arguments);

	compute_homology(graph, params);
}
