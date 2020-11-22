#pragma once

#include <cstring>
#include <iostream>
#include <string>

#include "../argparser.h"
#include "../definitions.h"
#include "../directed_graph.h"
#include "../parameters.h"

#include "flagser.h"
#include "h5.h"

filtered_directed_graph_t read_filtered_directed_graph(std::string input_filename, const flagser_parameters& params) {
#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	          << "reading in the graph" << std::flush << "\r";
#endif

	if (params.input_format == "h5" || params.hdf5_type.size() > 0 ||
	    (input_filename.rfind(".h5") != std::string::npos)) {
		return read_graph_h5(input_filename, params);
	}
	if (params.input_format == "flagser") return read_graph_flagser(input_filename, params);

#ifdef INDICATE_PROGRESS
	std::cout << "\033[K";
#endif

	std::string err_msg = "The input format \"" + params.input_format + "\" could not be found.";
	throw std::invalid_argument(err_msg);
}

std::vector<std::string> available_input_formats = {"flagser", "h5"};
