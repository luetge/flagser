#pragma once

#include <iostream>

#include "../input/input_classes.h"

void print_input_usage() {
	std::cerr << "  --in-format       the format of the input file. Defaults to \"flagser\". Available options are:"
	          << std::endl;

	for (auto f : available_input_formats) std::cerr << "                         " << f << std::endl;

	std::cerr << "  --h5-type type     the type of data in the h5-file. The type can either be \"matrix\"" << std::endl
	          << "                     if at the given path in the HDF5-file there is the connectivity matrix or"
	          << std::endl
	          << "                     \"grouped\" if the connectivity matrices are grouped. To only" << std::endl
	          << "                     consider a subset of the groups you can list them after" << std::endl
	          << "                     \"grouped\", e.g. \"grouped:L1_DAC,L2*\". The star is a placeholder" << std::endl
	          << "                     for arbitrary characters. The type defaults to \"matrix\" and is" << std::endl
	          << "                     only relevant for the input type h5." << std::endl
#ifndef WITH_HDF5
	          << " [HDF5 library not found]" << std::endl
#endif
	          << "  --undirected       compute the *undirected* flag complex" << std::endl
	          << "  --components       compute the directed flag complex for each individual connected" << std::endl
	          << "                     component of the input graph. Warning: this currently only works" << std::endl
	          << "                     for the trivial filtration. Additionally, isolated vertices are" << std::endl
	          << "                     ignored." << std::endl
	          << std::endl;
}