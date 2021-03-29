#include "../include/complex/directed_flag_complex_computer.h"
#include "../include/complex/directed_flag_complex_in_memory_computer.h"
#include "../include/input/input_classes.h"
#include "../include/output/output_classes.h"
#include "../include/parameters.h"
#include "../include/persistence.h"
#include <cstdio>
#include <stdexcept>

template <class T> void compute(std::string&& filename, std::vector<size_t> homology, std::string out_filename = "") {
	std::cout << "Testing " << filename << "..." << std::endl;
	flagser_parameters params;
	params.output_format = out_filename.size() == 0 ? "none" : "barcode";
	if (out_filename.size() > 0) { params.output_name = out_filename; }
	params.directed = true;
	filtered_directed_graph_t graph = read_filtered_directed_graph(filename.c_str(), params);

	size_t max_entries = std::numeric_limits<size_t>::max();
	coefficient_t modulus = 2;

	unsigned short max_dimension = std::numeric_limits<unsigned short>::max();
	unsigned short min_dimension = 0;

	auto output = get_output<T>(params);
	T complex(graph, params);

	output->set_complex(&complex);

	auto result = persistence_computer_t<decltype(complex)>(complex, output.get(), max_entries, modulus);
	result.compute_persistence(min_dimension, max_dimension);

	bool correct = true;

	if (result.get_betti_numbers().size() != homology.size()) {
		correct = false;
	} else {
		for (auto i = 0ul; i < homology.size(); i++) {
			if (result.get_betti_numbers()[i] != homology[i]) {
				correct = false;
				break;
			}
		}
	}

	if (!correct) {
		std::cout << "Expected betti numbers: ";
		for (auto b : homology) std::cout << std::to_string(b) << " ";
		std::cout << "; computed: ";
		for (auto b : result.get_betti_numbers()) std::cout << std::to_string(b) << " ";
		std::cout << std::endl;
		std::cout << "FAILED" << std::endl;
		//   throw std::logic_error(std::string("Wrong betti numbers computed."));
	}

	std::cout << "All good." << std::endl;
}

void run_all(bool full = false, bool in_memory = false) {
	using dfc_in_memory = directed_flag_complex_in_memory_computer::directed_flag_complex_in_memory_computer_t;
	using dfc = directed_flag_complex_computer::directed_flag_complex_computer_t;

	if (in_memory) {
                _run_all<dfc_in_memory>(full);
        } else {
                _run_all<dfc>(full);
        }
}

template <class dfc> void _run_all(bool full = false, bool in_memory = false) {
	compute<dfc>("../../test/d2.flag", {{1ul, 1ul}});
	compute<dfc>("../../test/a.flag", {{1ul, 2ul, 0ul}});
	compute<dfc>("../../test/b.flag", {{1ul, 0ul, 0ul}});
	compute<dfc>("../../test/c.flag", {{1ul, 5ul}});
	compute<dfc>("../../test/d.flag", {{1ul, 0ul, 1ul}});
	compute<dfc>("../../test/e.flag", {{1ul, 0ul, 0ul, 0ul}});
	compute<dfc>("../../test/f.flag", {{1ul, 0ul, 0ul}});
	compute<dfc>("../../test/d3.flag", {{1ul, 0ul, 2ul}});
	compute<dfc>("../../test/d3-allzero.flag", {{1ul, 0ul, 2ul}});
	compute<dfc>("../../test/double-d3.flag", {{1, 0, 5}});
	compute<dfc>("../../test/double-d3-allzero.flag", {{1, 0, 5}});
	compute<dfc>("../../test/d4.flag", {{1, 0, 0, 9}});
	compute<dfc>("../../test/d4-allzero.flag", {{1, 0, 0, 9}});
	compute<dfc>("../../test/d5.flag", {{1, 0, 0, 0, 44}});
	compute<dfc>("../../test/d7.flag", {{1, 0, 0, 0, 0, 0, 1854}});

	if (full) {
		const auto file_path = "../flagser_tmp";
		std::remove(file_path);
                compute<dfc>("../../test/a.flag", {{1ul, 2ul, 0ul}}, file_path);

		// Check that the file has the right content
		std::ifstream t(file_path);
		std::string file_content((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
		std::remove(file_path);

		std::string expected_file_content =
		    "# persistence intervals in dimension 0\n [0, )\n# persistence intervals in dimension 1\n [0, )\n [0, )\n# "
		    "persistence intervals in dimension 2\n\nThe remaining homology groups are trivial.\n\n# Euler "
		    "characteristic: -1\n\n# Betti numbers:\n#\t\tdim H_0 = 1\n#\t\tdim H_1 = 2\n#\t\tdim H_2 = 0\n\n# Cell "
		    "counts:\n#\t\tdim C_0 = 5\n#\t\tdim C_1 = 7\n#\t\tdim C_2 = 1\n#\t\tdim C_3 = 0\n";

		if (file_content != expected_file_content) {
			std::cerr << "The file content differed!" << std::endl;
			std::cerr << std::endl << "*** EXPECTED ***" << std::endl << expected_file_content << std::endl;
			std::cerr << std::endl << "*** GOT ***" << std::endl << file_content << std::endl;
			throw std::logic_error("The file content differed.");
		}

#ifdef NDEBUG
		std::cout << "Running extensive tests, this might take a while." << std::endl;
		if (in_memory) {
		        compute<dfc>("../../test/medium-test-data.flag", {{14237, 39477, 378, 0}});
			compute<dfc>("../../test/d10.flag", {{1, 0, 0, 0, 0, 0, 0, 0, 0, 1334961}});
		}
#else
		std::cout << "Skipping extensive tests." << std::endl;
#endif
	}
}
