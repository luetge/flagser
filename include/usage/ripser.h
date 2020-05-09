#pragma once

#include <iostream>

#include "../filtration_algorithms.h"
#include "../output/output_classes.h"
#include "help.h"
#include "homology.h"
#include "input.h"
#include "output.h"

void print_usage_and_exit(int exit_code) {
	std::cerr << "Computes the persistent homology of the Vietoris-Rips complex of a point cloud." << std::endl
	          << "Usage: ripser [options] filename" << std::endl
	          << std::endl
	          << "Options:" << std::endl
	          << std::endl;

	print_output_usage();
	std::cerr << "  --format         use the specified file format for the input. "
	             "Available options are:"
	          << std::endl
	          << "                         lower-distance (lower triangular distance "
	             "matrix; default)"
	          << std::endl
	          << "                         upper-distance (upper triangular distance "
	             "matrix)"
	          << std::endl
	          << "                         distance       (full distance matrix)" << std::endl
	          << "                         point-cloud    (point cloud in Euclidean space)" << std::endl
	          << "                         dipha          (distance matrix in DIPHA file "
	             "format)"
	          << std::endl
	          << "  --threshold <t>  compute Rips complexes up to filtration <t>" << std::endl;
	print_homology_usage(true);
	print_help_usage();

	exit(exit_code);
}