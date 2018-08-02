#pragma once
#include <iostream>

#include "cell-output.h"
#include "help.h"
#include "input.h"

void print_usage_and_exit(int exit_code) {
	std::cerr << "Computes the Euler characteristic and cell counts of the directed flag complex of a graph."
	          << std::endl
	          << "Usage: flagser-count [options] filename" << std::endl
	          << "Options:" << std::endl
	          << std::endl;

	print_input_usage();
	print_cell_output_usage();
	print_help_usage();

	exit(exit_code);
}