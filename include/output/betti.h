#pragma once

#include <iostream>

#include "../definitions.h"
#include "../parameters.h"
#include "base.h"

template <typename Complex> class betti_output_t : public file_output_t<Complex> {
	Complex* complex;
	unsigned short current_dimension = 0;
	unsigned short min_dimension;
	unsigned short max_dimension;
	coefficient_t modulus;
	std::vector<size_t> betti;
	std::vector<size_t> skipped;
	index_t euler_characteristic = 0;
	bool approximate_computation;
	bool aggregate_results = false;

	static std::vector<size_t> total_betti;
	static std::vector<size_t> total_skipped;
	static std::vector<size_t> total_cell_count;

	static size_t total_top_dimension;

public:
	betti_output_t(const flagser_parameters& params)
	    : file_output_t<Complex>(params.output_name), min_dimension(params.min_dimension),
	      max_dimension(params.max_dimension), modulus(params.modulus),
	      approximate_computation(params.approximate_computation),
	      aggregate_results(!params.split_into_connected_components) {}

	virtual void finished(bool with_cell_counts = true) override;
	virtual void set_complex(Complex* _complex) override { complex = _complex; }
	virtual void computing_barcodes_in_dimension(unsigned short dimension) override { current_dimension = dimension; }
	virtual void betti_number(size_t _betti, size_t _skipped) override {
		betti.resize(current_dimension + 1, 0);
		skipped.resize(current_dimension + 1, 0);
		betti[current_dimension] = _betti;
		skipped[current_dimension] = _skipped;
		euler_characteristic += index_t((current_dimension & 1 ? -1 : 1) * _betti);
	}

	virtual void print_aggregated_results() override;
};

template <typename Complex> std::vector<size_t> betti_output_t<Complex>::total_betti;
template <typename Complex> std::vector<size_t> betti_output_t<Complex>::total_skipped;
template <typename Complex> std::vector<size_t> betti_output_t<Complex>::total_cell_count;

template <typename Complex> size_t betti_output_t<Complex>::total_top_dimension = 0;

void print_ordinary(std::ostream& outstream, int min_dimension, int max_dimension, int top_dimension,
                    std::vector<size_t> number_of_cells, std::vector<size_t> betti, std::vector<size_t> skipped,
                    bool approximate_computation, bool with_cell_counts) {
	if (with_cell_counts) {
		outstream << "# Cell counts (of dimensions between " << std::max(0, min_dimension - 1) << " and "
		          << top_dimension << "):" << std::endl;
		for (int i = std::max(0, min_dimension - 1); i <= top_dimension; i++)
			outstream << number_of_cells[i] << (i < top_dimension ? " " : "");
		outstream << std::endl;

		index_t cell_euler_characteristic = 0;
		for (int i = 0; i <= top_dimension; i++)
			cell_euler_characteristic += index_t((i % 2 == 1 ? -1 : 1) * number_of_cells[i]);

		bool computed_full_homology = min_dimension == 0 && max_dimension == std::numeric_limits<unsigned short>::max();
		if (computed_full_homology) {
			outstream << "# Euler characteristic:" << std::endl;
			outstream << cell_euler_characteristic << std::endl;
		}
	}

	if (approximate_computation) {
		outstream << "# Skipped columns of the coboundary matrix (in dimensions between "
		          << std::max(0, min_dimension - 1) << " and " << skipped.size() - 1 << "):" << std::endl;
		for (size_t idx = std::max(0, min_dimension - 1); idx < skipped.size(); idx++) {
			outstream << skipped[idx] << (idx < betti.size() ? " " : "");
		}
		outstream << std::endl;
	}

	outstream << "# Betti numbers (between " << min_dimension << " and " << (betti.size() - 1) << "):" << std::endl;
	for (size_t idx = min_dimension; idx < betti.size(); idx++) {
		outstream << betti[idx] << (idx < betti.size() ? " " : "");
	}

	outstream << std::endl;
}

template <typename Complex> void betti_output_t<Complex>::finished(bool with_cell_counts) {
	std::vector<size_t> number_of_cells;
	total_top_dimension = std::max(total_top_dimension, complex->top_dimension());

	if (with_cell_counts) {
		for (auto i = 0ul; i <= complex->top_dimension(); i++) number_of_cells.push_back(complex->number_of_cells(i));

		total_cell_count.resize(complex->top_dimension() + 1, 0);
		for (auto i = 0ul; i <= complex->top_dimension(); i++) total_cell_count[i] += complex->number_of_cells(i);
	} else {
		total_cell_count.resize(0);
	}

	print_ordinary(file_output_t<Complex>::outstream, int(min_dimension), int(max_dimension),
	               int(complex->top_dimension()), number_of_cells, betti, skipped, approximate_computation,
	               with_cell_counts);

	total_betti.resize(betti.size(), 0);
	for (size_t idx = 0; idx < betti.size(); idx++) total_betti[idx] += betti[idx];

	total_skipped.resize(skipped.size(), 0);
	for (size_t idx = 0; idx < skipped.size(); idx++) total_skipped[idx] += skipped[idx];
}

template <typename Complex> void betti_output_t<Complex>::print_aggregated_results() {
	if (!aggregate_results) return;

	print_ordinary(file_output_t<Complex>::outstream, int(min_dimension), int(max_dimension), int(total_top_dimension),
	               total_cell_count, total_betti, total_skipped, approximate_computation, total_cell_count.size() > 0);
}
