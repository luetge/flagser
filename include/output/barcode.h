#pragma once

#include <iostream>

#include "../argparser.h"
#include "../definitions.h"
#include "../parameters.h"
#include "base.h"

template <typename Complex> class barcode_output_t : public file_output_t<Complex> {
	Complex* complex;
	unsigned short current_dimension = 0;
	unsigned short min_dimension;
	unsigned short max_dimension;
	coefficient_t modulus;
	std::vector<size_t> betti;
	std::vector<size_t> skipped;
	index_t euler_characteristic = 0;

public:
	barcode_output_t(const flagser_parameters& params)
	    : file_output_t<Complex>(params.output_name), min_dimension(params.min_dimension),
	      max_dimension(params.max_dimension), modulus(params.modulus) {}

	void set_complex(Complex* _complex) { complex = _complex; }

	void finished(bool with_cell_counts = true) {

		index_t cell_euler_characteristic = 0;
		if (with_cell_counts) {
			for (size_t i = 0; i <= complex->top_dimension(); i++)
				cell_euler_characteristic += (i % 2 == 1 ? -1 : 1) * index_t(complex->number_of_cells(index_t(i)));

			bool computed_full_homology =
			    min_dimension == 0 && max_dimension == std::numeric_limits<unsigned short>::max();
			if (computed_full_homology) {
				file_output_t<Complex>::outstream << std::endl;
				file_output_t<Complex>::outstream << "# Euler characteristic: " << cell_euler_characteristic
				                                  << std::endl;
			}
		}

		file_output_t<Complex>::outstream << std::endl;
		file_output_t<Complex>::outstream << "# Betti numbers:" << std::endl;
		int dim = 0;
		for (size_t idx = 0; idx < betti.size(); idx++) {
			file_output_t<Complex>::outstream << "#\t\tdim H_" << (min_dimension + dim++) << " = " << betti[idx];
			if (skipped[idx] > 0) file_output_t<Complex>::outstream << " (skipped " << skipped[idx] << ")";
			file_output_t<Complex>::outstream << std::endl;
		}

		if (with_cell_counts) {
			file_output_t<Complex>::outstream << std::endl;
			file_output_t<Complex>::outstream << "# Cell counts:" << std::endl;
			for (size_t i = size_t(std::max(0, min_dimension - 1)); i <= complex->top_dimension(); i++)
				file_output_t<Complex>::outstream << "#\t\tdim C_" << i << " = " << complex->number_of_cells(index_t(i))
				                                  << std::endl;
		}
	}
	void computing_barcodes_in_dimension(unsigned short dimension) {
		if (dimension >= min_dimension && dimension <= max_dimension)
			file_output_t<Complex>::outstream << "# persistence intervals in dimension " << dimension << std::endl;
		current_dimension = dimension;
	}
	inline void new_barcode(value_t birth, value_t death) {
		file_output_t<Complex>::outstream << " [" << birth << ", " << death << ")" << std::endl;
	}
	inline void new_infinite_barcode(value_t birth) {
		file_output_t<Complex>::outstream << " [" << birth << ", )" << std::endl;
	}
	inline void skipped_column(value_t birth) {
		file_output_t<Complex>::outstream << " [?" << birth << ", ?)" << std::endl;
	}
	void betti_number(size_t _betti, size_t _skipped) {
		const auto shifted_dimension = current_dimension - min_dimension;
		betti.resize(shifted_dimension + 1, 0);
		skipped.resize(shifted_dimension + 1, 0);
		betti[shifted_dimension] = _betti;
		skipped[shifted_dimension] = _skipped;
		euler_characteristic += (shifted_dimension & 1 ? -1 : 1) * index_t(_betti);
	}
	void remaining_homology_is_trivial() {
		file_output_t<Complex>::outstream << std::endl << "The remaining homology groups are trivial." << std::endl;
	}
};
