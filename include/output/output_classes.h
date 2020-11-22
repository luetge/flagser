#pragma once

#include <cstring>
#include <iostream>
#include <memory>

#include "../definitions.h"

#include "../argparser.h"
#include "barcode.h"
#include "base.h"
#include "betti.h"
#include "trivial.h"

#ifdef WITH_HDF5
#include "barcode_hdf5.h"
#endif

/* issues with compiler version and std::make_unique were discovered
 * see: https://github.com/luetge/flagser/pull/20
 */
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

// aparently clang defines also gnuc as macros, need to check first for clang
#if !defined(__clang__) && (defined(__GNUC__) || defined(__GNUG__)) && GCC_VERSION < 40900
template <typename T, typename... Args> std::unique_ptr<T> make_unique(Args&&... args) {
	return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#endif

template <typename Complex> std::unique_ptr<output_t<Complex>> get_output(const flagser_parameters& params) {
	// using std namespace because of std::make_unique workaround
	// This restrict the scope to inside the function
	using namespace std;

	if (params.output_format == "betti") return make_unique<betti_output_t<Complex>>(params);
	if (params.output_format == "barcode") return make_unique<barcode_output_t<Complex>>(params);

	if (params.output_format == "none") return make_unique<trivial_output_t<Complex>>(params);

#ifdef WITH_HDF5
	if (params.output_format == "barcode:hdf5") return make_unique<barcode_hdf5_output_t<Complex>>(params);
#endif

	std::string err_msg = "The output format \"" + params.output_format + "\" could not be found.";
	throw std::invalid_argument(err_msg);
}

std::vector<std::string> available_output_formats = {"barcode", "betti"
#ifdef WITH_HDF5
                                                     ,
                                                     "barcode:hdf5"
#endif
};
