#pragma once

#include <cctype>
#include <fstream>
#include <string>

#include "../argparser.h"
#include "../definitions.h"
#include "../directed_graph.h"

void open_file(const std::string filename, std::ifstream& file_stream) {
	file_stream.open(filename);
	if (file_stream.fail()) {
		std::string err_msg = "couldn't open file " + filename;
		throw std::runtime_error(err_msg);
	}
}

inline std::string trim(const std::string& s) {
	auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) { return std::isspace(c); });
	auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) { return std::isspace(c); }).base();
	return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}
