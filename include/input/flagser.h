#pragma once

#include <functional>

#ifdef _MSC_VER

// Windows compiler does not support `not`, `and`, etc. keywords out of the box
#include <iso646.h>

#endif

#include "../parameters.h"
#include "base.h"

// String to number conversion
unsigned int string_to_uint(std::string s) { return atoi(s.c_str()); }
float string_to_float(std::string s) { return float(atof(s.c_str())); }
template <typename t>
std::vector<t> split(const std::string& s, char delim, const std::function<t(std::string)>& transform) {
	std::vector<t> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) elems.push_back(transform(item));
	return elems;
}

enum HAS_EDGE_FILTRATION { TOO_EARLY_TO_DECIDE, MAYBE, YES, NO };
filtered_directed_graph_t read_graph_flagser(const std::string filename, const flagser_parameters& params) {
	std::string line;
	filtered_directed_graph_t graph{};
	int current_dimension = 0;
	std::vector<value_t> vertex_filtration;
	HAS_EDGE_FILTRATION has_edge_filtration = HAS_EDGE_FILTRATION::TOO_EARLY_TO_DECIDE;

	std::ifstream input_stream;
	open_file(filename, input_stream);

	while (not input_stream.eof()) {
		getline(input_stream, line);
		line = trim(line);
		if (line.length() == 0) continue;
		if (line[0] == 'd' && line[1] == 'i' && line[2] == 'm') {
			if (line[4] == '1') {
				graph = filtered_directed_graph_t(vertex_filtration, params.directed);
				current_dimension = 1;
				has_edge_filtration = HAS_EDGE_FILTRATION::MAYBE;
			}
			continue;
		}

		if (current_dimension == 0) {
			if (vertex_filtration.empty()) vertex_filtration = split<value_t>(line, ' ', string_to_float);
		} else {
			if (has_edge_filtration == HAS_EDGE_FILTRATION::MAYBE) {
				// If the edge has three components, then there are also filtration values, which we assume to come last
				int number_of_spaces = 0;
				for (auto& iter : line)
					if (iter == ' ') { number_of_spaces++; }
				has_edge_filtration = number_of_spaces == 1 ? HAS_EDGE_FILTRATION::NO : HAS_EDGE_FILTRATION::YES;
			}

			if (has_edge_filtration == NO) {
				std::vector<vertex_index_t> vertices = split<vertex_index_t>(line, ' ', string_to_uint);
				graph.add_edge(vertices[0], vertices[1]);
			} else {
				std::vector<value_t> vertices = split<value_t>(line, ' ', string_to_float);
				if (value_t(vertices[2]) <
				    std::max(vertex_filtration[size_t(vertices[0])], vertex_filtration[size_t(vertices[1])])) {
					std::string err_msg =
					    "The flagser file contained an edge filtration that contradicts the vertex "
					    "filtration, the edge (" +
					    std::to_string(vertices[0]) + ", " + std::to_string(vertices[1]) + ") has filtration value " +
					    std::to_string(vertices[2]) + ", which is lower than min(" +
					    std::to_string(vertex_filtration[size_t(vertices[0])]) + ", " +
					    std::to_string(vertex_filtration[size_t(vertices[1])]) + "), the filtrations of its vertices.";
					throw std::runtime_error(err_msg);
				}
				graph.add_filtered_edge((vertex_index_t)vertices[0], (vertex_index_t)vertices[1], vertices[2]);
			}
		}
	}

	return graph;
}
