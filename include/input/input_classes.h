#pragma once

#include <iostream>
#include <string>
#include <cstring>

#include "../argparser.h"
#include "../definitions.h"
#include "../directed_graph.h"

#include "flagser.h"
#include "h5.h"

filtered_directed_graph_t read_filtered_directed_graph(std::string input_filename,
                                                       const named_arguments_t& named_arguments) {
	std::string input_name = "flagser";
	auto it = named_arguments.find("in-format");
	if (it != named_arguments.end()) { input_name = it->second; }

#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	          << "reading in the graph" << std::flush << "\r";
#endif

	if (input_name == "h5" || strlen(get_argument_or_default(named_arguments, "h5-type", "")) > 0 ||
	    (!argument_was_passed(named_arguments, "in-format") &&
	     input_filename.rfind(".h5") != std::string::npos)) {
		return read_graph_h5(input_filename, named_arguments);
	}
	if (input_name == "flagser") return read_graph_flagser(input_filename, named_arguments);

#ifdef INDICATE_PROGRESS
	std::cout << "\033[K";
#endif

	std::cerr << "The input format \"" << input_name << "\" could not be found." << std::endl;
	exit(1);
}

std::vector<std::string> available_input_formats = {"flagser", "h5"};