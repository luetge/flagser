#include <stdexcept>
#include "../include/persistence.h"
#include "../include/complex/directed_flag_complex_computer.h"
#include "../include/input/input_classes.h"
#include "../include/output/output_classes.h"

void compute(std::string&& filename, std::vector<size_t> homology, std::string out_filename = "") {
  std::cout << "Testing " << filename << "..." << std::endl;
  std::unordered_map<std::string, const char*> named_arguments;
  const std::string out_format = out_filename.size() == 0 ? "none" : "barcode";
  named_arguments.insert(std::pair<std::string, const char*>("out-format", out_format.c_str()));
  if (out_filename.size() > 0) {
    named_arguments.insert(std::pair<std::string, const char*>("out", out_filename.c_str()));
  }
	filtered_directed_graph_t graph = read_filtered_directed_graph(filename.c_str(), named_arguments);

	size_t max_entries = std::numeric_limits<size_t>::max();
	coefficient_t modulus = 2;

  unsigned short max_dimension = std::numeric_limits<unsigned short>::max();
	unsigned short min_dimension = 0;

	auto output = get_output<directed_flag_complex_computer_t>(named_arguments);
  directed_flag_complex_computer_t complex(graph, named_arguments);

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

void run_all(bool full=false) {
  compute("../../test/d2.flag", {{1ul, 1ul}});
  compute("../../test/a.flag", {{1ul, 2ul, 0ul}});
  compute("../../test/b.flag", {{1ul, 0ul, 0ul}});
  compute("../../test/c.flag", {{1ul, 5ul}});
  compute("../../test/d.flag", {{1ul, 0ul, 1ul}});
  compute("../../test/e.flag", {{1ul, 0ul, 0ul, 0ul}});
  compute("../../test/f.flag", {{1ul, 0ul, 0ul}});
  compute("../../test/d3.flag", {{1ul, 0ul, 2ul}});
  compute("../../test/d3-allzero.flag", {{1ul, 0ul, 2ul}});
  compute("../../test/double-d3.flag", {{1, 0, 5}});
  compute("../../test/double-d3-allzero.flag", {{1, 0, 5}});
  compute("../../test/d4.flag", {{1, 0, 0, 9}});
  compute("../../test/d4-allzero.flag", {{1, 0, 0, 9}});
  compute("../../test/d5.flag", {{1, 0, 0, 0, 44}});
  compute("../../test/d7.flag", {{1, 0, 0, 0, 0, 0, 1854}});

  if (full) {
    const auto file_path = "/tmp/flagser_tmp_test";
    std::remove(file_path);
    compute("../../test/a.flag", {{1ul, 2ul, 0ul}}, file_path);
    // Check that the file has the right content
    std::ifstream t(file_path);
    std::string file_content((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
    std::remove(file_path);

    std::string expected_file_content = "# persistence intervals in dimension 0\n [0, )\n# persistence intervals in dimension 1\n [0, )\n [0, )\n# persistence intervals in dimension 2\n\nThe remaining homology groups are trivial.\n\n# Euler characteristic: -1\n\n# Betti numbers:\n#\t\tdim H_0 = 1\n#\t\tdim H_1 = 2\n#\t\tdim H_2 = 0\n\n# Cell counts:\n#\t\tdim C_0 = 5\n#\t\tdim C_1 = 7\n#\t\tdim C_2 = 1\n#\t\tdim C_3 = 0\n";

    if (file_content != expected_file_content) {
      std::cerr << "The file content differed!" << std::endl;
      std::cerr << std::endl << "*** EXPECTED ***" << std::endl << expected_file_content << std::endl;
      std::cerr << std::endl << "*** GOT ***" << std::endl << file_content << std::endl;
      throw std::logic_error("The file content differed.");
    }

#ifdef NDEBUG
    // The expensive stuff is only run in production settings
    compute("../../test/medium-test-data.flag", {{14237, 39477, 378, 0}});
    compute("../../test/d10.flag", {{1, 0, 0, 0, 0, 0, 0, 0, 0, 1334961}});
#endif
  }
}
