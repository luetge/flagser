#pragma once

#include <iostream>

#include "../definitions.h"
#include "base.h"

template <typename Complex> class trivial_output_t : public output_t<Complex> {
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
	trivial_output_t(const named_arguments_t& named_arguments)
	    : modulus(atoi(get_argument_or_default(named_arguments, "modulus", "2"))),
	      min_dimension(atoi(get_argument_or_default(named_arguments, "min-dim", "0"))),
	      max_dimension(atoi(get_argument_or_default(named_arguments, "max-dim", "65535"))),
	      approximate_computation(argument_was_passed(named_arguments, "approximate")),
	      aggregate_results(argument_was_passed(named_arguments, "components")) {}

	virtual void finished(bool with_cell_counts = true) override;
	virtual void set_complex(Complex* _complex) override { complex = _complex; }
	virtual void computing_barcodes_in_dimension(unsigned short dimension) override { current_dimension = dimension; }
	virtual void betti_number(size_t _betti, size_t _skipped) override {
		betti.resize(current_dimension + 1, 0);
		skipped.resize(current_dimension + 1, 0);
		betti[current_dimension] = _betti;
		skipped[current_dimension] = _skipped;
		euler_characteristic += (current_dimension & 1 ? -1 : 1) * _betti;
	}

	virtual void print_aggregated_results() override;
};

template <typename Complex> std::vector<size_t> trivial_output_t<Complex>::total_betti;
template <typename Complex> std::vector<size_t> trivial_output_t<Complex>::total_skipped;
template <typename Complex> std::vector<size_t> trivial_output_t<Complex>::total_cell_count;

template <typename Complex> size_t trivial_output_t<Complex>::total_top_dimension = 0;

template <typename Complex> void trivial_output_t<Complex>::finished(bool with_cell_counts) {}

template <typename Complex> void trivial_output_t<Complex>::print_aggregated_results() {}