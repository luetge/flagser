#define ASSEMBLE_REDUCTION_MATRIX
#define INDICATE_PROGRESS
#define SKIP_APPARENT_PAIRS
#define USE_ARRAY_HASHMAP
#define USE_CELLS_WITHOUT_DIMENSION
#define SORT_COLUMNS_BY_PIVOT
#define RETRIEVE_PERSISTENCE
// #define WITH_HDF5
// #define KEEP_FLAG_COMPLEX_IN_MEMORY
// #define USE_COEFFICIENTS
// #define MANY_VERTICES


#include "../include/persistence.h"
#include "../include/complex/directed_flag_complex_computer.h"
#include "../include/input/input_classes.h"
#include "../include/output/output_classes.h"


void compute(std::string&& filename, std::vector<size_t> homology) {
  std::cout << "Testing " << filename << "..." << std::endl;
  std::unordered_map<std::string, const char*> named_arguments;
  const std::string out_format = "none";
  named_arguments.insert(std::pair<std::string, const char*>("out-format", out_format.c_str()));
	filtered_directed_graph_t graph = read_filtered_directed_graph(filename.c_str(), named_arguments);

	size_t max_entries = std::numeric_limits<size_t>::max();
	coefficient_t modulus = 2;

  unsigned short max_dimension = std::numeric_limits<unsigned short>::max();
	unsigned short min_dimension = 0;

	auto output = get_output<directed_flag_complex_computer_t>(named_arguments);
	size_t component_number = 1;
  directed_flag_complex_computer_t complex(graph, named_arguments);

  output->set_complex(&complex);

  auto result = persistence_computer_t<decltype(complex)>(complex, output, max_entries, modulus);
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
    std::cerr << "Expected betti numbers: ";
    for (auto b : homology) std::cerr << std::to_string(b) << " ";
    std::cerr << "; computed: ";
    for (auto b : result.get_betti_numbers()) std::cerr << std::to_string(b) << " ";
    std::cerr << std::endl;
    throw std::logic_error(std::string("Wrong betti numbers computed."));
  }

  std::cout << "All good." << std::endl;
}

int main(int argc, char** argv) {
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
  compute("../../test/medium-test-data.flag", {{14237, 39477, 378, 0}});
  compute("../../test/d10.flag", {{1, 0, 0, 0, 0, 0, 0, 0, 0, 1334961}});

  return 0;
}
