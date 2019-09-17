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

#include "../include/argparser.h"
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
typedef directed_flag_complex_in_memory_computer_t
    directed_flag_complex_compute_t;
#else
typedef directed_flag_complex_computer_t 
    directed_flag_complex_compute_t;
#endif

void compute_homology(filtered_directed_graph_t& graph, const named_arguments_t& named_arguments, size_t max_entries,
                      coefficient_t modulus) {

	unsigned short max_dimension = std::numeric_limits<unsigned short>::max();
	unsigned short min_dimension = 0;
	bool split_into_connected_components = named_arguments.find("components") != named_arguments.end();

	named_arguments_t::const_iterator it;
	if ((it = named_arguments.find("max-dim")) != named_arguments.end()) { max_dimension = atoi(it->second); }
	if ((it = named_arguments.find("min-dim")) != named_arguments.end()) { min_dimension = atoi(it->second); }

	std::vector<filtered_directed_graph_t> subgraphs{graph};
	if (split_into_connected_components) { subgraphs = graph.get_connected_subgraphs(2); }

	auto output = get_output<directed_flag_complex_compute_t>(named_arguments);

	size_t component_number = 1;
	for (auto subgraph : subgraphs) {
		directed_flag_complex_compute_t complex(subgraph, named_arguments);

		output->set_complex(&complex);
		if (split_into_connected_components) {
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

		persistence_computer_t<decltype(complex)> persistence_computer(complex, output, max_entries, modulus);
		persistence_computer.compute_persistence(min_dimension, max_dimension);
	}

	if (split_into_connected_components) { output->print("\n## Total\n"); }

	output->print_aggregated_results();
}

int main(int argc, char** argv) {
	auto arguments = parse_arguments(argc, argv);

	auto positional_arguments = get_positional_arguments(arguments);
	auto named_arguments = get_named_arguments(arguments);
	if (named_arguments.find("help") != named_arguments.end()) { print_usage_and_exit(-1); }

	if (positional_arguments.size() == 0) { print_usage_and_exit(-1); }
	const char* input_filename = positional_arguments[0];

	filtered_directed_graph_t graph = read_filtered_directed_graph(input_filename, named_arguments);

	size_t max_entries = std::numeric_limits<size_t>::max();
	coefficient_t modulus = 2;
	named_arguments_t::const_iterator it;
	if ((it = named_arguments.find("approximate")) != named_arguments.end()) { max_entries = atoi(it->second); }
#ifdef USE_COEFFICIENTS
	if ((it = named_arguments.find("modulus")) != named_arguments.end()) { modulus = atoi(it->second); }
#endif

	compute_homology(graph, named_arguments, max_entries, modulus);
}
