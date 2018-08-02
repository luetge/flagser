#pragma once

#include <iostream>

#include "../output/output_classes.h"

void print_cell_output_usage() {
	std::cerr << "  --out filename     write lists of all cells into an HDF5 file" << std::endl;
}