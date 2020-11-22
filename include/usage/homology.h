#pragma once

#include <iostream>

#include "../filtration_algorithms.h"

void print_homology_usage(bool no_filtration = false) {
	std::cerr << "  --cache folder     cache the coboundary matrix in the given folder" << std::endl;
	if (!no_filtration) {
		std::cerr << "  --filtration       use the specified algorithm to compute the filtration. Options are:"
		          << std::endl;

		for (auto f : custom_filtration_computer) std::cerr << "                         " << f << std::endl;
	}

	std::cerr << "  --max-dim          the maximal homology dimension to be computed" << std::endl
	          << "  --min-dim          the minimal homology dimension to be computed" << std::endl
#ifdef USE_COEFFICIENTS
	          << "  --modulus          compute homology with coefficients in the prime field Z/<p>Z" << std::endl
#endif
	          << "  --approximate n    skip all columns creating columns in the reduction matrix with" << std::endl
	          << "                     n non-trivial entries. Use this for hard problems, a good value" << std::endl
	          << "                     is often 100000. Increase for higher precision, decrease for faster computation." << std::endl
	          << "  --threads          number of threads to use for the computation (default: number of threads that your machine can execute simultaneously)" << std::endl
	          << std::endl;
}
