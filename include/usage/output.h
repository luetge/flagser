#pragma once

#include <iostream>

#include "../output/output_classes.h"

void print_output_usage() {
	std::cerr << "  --out filename     write the barcodes to the given file" << std::endl
	          << "  --out-format       the format of the output file. Defaults to \"barcode\". Available options are:"
	          << std::endl;

	for (auto f : available_output_formats) std::cerr << "                         " << f << std::endl;
}