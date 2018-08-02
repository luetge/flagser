#pragma once

#include <iostream>

#include "../filtration_algorithms.h"
#include "../output/output_classes.h"
#include "help.h"
#include "homology.h"
#include "input.h"
#include "output.h"

void print_usage_and_exit(int exit_code) {
	std::cerr << "Computes the persistent homology of the directed flag complex of the given directed graph."
	          << std::endl
	          << "Usage: flagser [options] filename" << std::endl
	          << std::endl
	          << "Options:" << std::endl
	          << std::endl;

	print_output_usage();
	print_input_usage();
	print_homology_usage();
	print_help_usage();

	std::cerr << std::endl
	          << std::endl
	          << "If you updated the filtration algorithms in the file \"algorithms.math\" and they do not" << std::endl
	          << "show up here, try running \"make\" before running flagser." << std::endl
	          << std::endl;

	exit(exit_code);
}