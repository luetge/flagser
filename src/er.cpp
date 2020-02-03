/**
 * Builds an Erdösz-Renyi graph.
 */

#include <ctime>
#include <fstream>
#include <iostream>

void print_usage_and_exit(int exit_code) {
	std::cerr << "Usage: "
	          << "er "
	          << "[options] n filename" << std::endl
	          << std::endl
	          << "Options:" << std::endl
	          << std::endl
	          << "  --density                    the average density" << std::endl
	          << "  --random-edge-filtration     generates random filtration "
	             "values in [0,1] for all edges"
	          << std::endl
	          << "  --undirected                 create a graph such that the directed and"
	          << " undirected flag complexes coincide" << std::endl
	          << "  --help                       print this screen" << std::endl
	          << std::endl;

	exit(exit_code);
}

int main(int argc, char** argv) {

	const char* filename = nullptr;
	int n = -1;
	float density = 1.0;
	bool edge_filtration = false;
	bool undirected = false;

	for (int i = 0; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--density") {
			density = float(atof(argv[++i]));
		} else if (arg == "--random-edge-filtration") {
			edge_filtration = true;
		} else if (arg == "--undirected") {
			undirected = true;
		} else {
			if (n > -1) {
				filename = argv[i];
			} else
				n = atoi(argv[i]);
		}
	}

	if (n == -1) { print_usage_and_exit(-1); }

	std::ofstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	srand((unsigned int)time(NULL));
	std::ostream& output_stream = filename ? file_stream : std::cout;
	output_stream << "dim 0" << std::endl;
	for (int i = 0; i < n; i++) output_stream << 0 << " ";
	output_stream << std::endl;

	output_stream << "dim 1" << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (density == 1 || (rand() % 1000) / 1000.0 < density) {
				output_stream << i << " " << j;
				if (edge_filtration) output_stream << " " << (1 + (rand() % 100) / 100.0);
				output_stream << std::endl;
			}

			if (undirected) continue;

			if (density == 1 || (rand() % 1000) / 1000.0 < density) {
				output_stream << j << " " << i;
				if (edge_filtration) output_stream << " " << (1 + (rand() % 100) / 100.0);
				output_stream << std::endl;
			}
		}
	}

	if (filename) {
		file_stream.close();
		std::cout << "Wrote directed graph to " << filename << "." << std::endl;
	}
}
