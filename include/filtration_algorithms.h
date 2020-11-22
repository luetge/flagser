#pragma once

#include <algorithm>
#include <exception>
#include <math.h>

#include "complex/directed_flag_complex.h"
#include "directed_graph.h"

struct filtration_algorithm_t {
	filtration_algorithm_t() = default;
	virtual ~filtration_algorithm_t() = default;
	virtual inline value_t compute_filtration(unsigned short /* dimension */,
	                                          const directed_flag_complex_cell_t& /* cell */,
	                                          const filtered_directed_graph_t& /* graph */,
	                                          const value_t* /* boundary_filtration */) const {
		return 0.0f;
	}
	virtual inline bool overwrite_vertex_filtration() const { return false; }
	virtual inline bool overwrite_edge_filtration() const { return false; }
	virtual inline bool needs_face_filtration() const { return true; }
};

template <typename T> inline value_t max(T values, int from, int to) {
	value_t max = 0;
	for (int index = from; index <= to; index++) {
		value_t next_value = values[index];
		max = index == 0 || next_value > max ? next_value : max;
	}

	return max;
}

template <typename T> inline value_t min(T values, int from, int to) {
	value_t min = 0;
	for (int index = from; index <= to; index++) {
		value_t next_value = values[index];
		min = index == 0 || next_value < min ? next_value : min;
	}

	return min;
}

template <typename T> inline value_t sum(T values, int from, int to) {
	value_t sum = 0;
	for (int index = from; index <= to; index++) { sum += values[index]; }

	return sum;
}

template <typename T> value_t product(T values, int from, int to) {
	value_t product = 1;
	for (int index = from; index <= to; index++) { product *= values[index]; }

	return product;
}

#include "./algorithms.h"

filtration_algorithm_t* get_filtration_computer(std::string algorithm) {
	std::transform(algorithm.begin(), algorithm.end(), algorithm.begin(), ::tolower);

	// For the trivial filtration, skip the computation of filtration values altogether
	if (algorithm == "zero" || algorithm == "") return nullptr;

	filtration_algorithm_t* custom = get_custom_filtration_computer(algorithm);

	if (custom != nullptr) return custom;

	std::string err_msg = "The filtration algorithm \"" + algorithm + "\" could not be found.";
	throw std::invalid_argument(err_msg);
}
