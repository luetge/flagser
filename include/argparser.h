#pragma once

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

typedef std::vector<const char*> positional_arguments_t;
typedef std::unordered_map<std::string, const char*> named_arguments_t;
typedef std::pair<positional_arguments_t, named_arguments_t> cli_arguments_t;

const positional_arguments_t get_positional_arguments(const cli_arguments_t& args) { return args.first; }
const named_arguments_t get_named_arguments(const cli_arguments_t& args) { return args.second; }

const char* get_argument_or_default(const named_arguments_t& args, const std::string name, const char* def) {
	auto it = args.find(name);
	if (it == args.end()) return def;
	return it->second;
}

const char* get_argument_or_fail(const named_arguments_t& args, const std::string name, const std::string msg) {
	auto it = args.find(name);
	if (it != args.end()) return it->second;

	throw std::invalid_argument(msg.c_str());
}

bool argument_was_passed(const named_arguments_t& args, const std::string name) {
	return args.find(name) != args.end();
}

// We assume that every argument *except for the help argument* has extra data attached
const cli_arguments_t parse_arguments(int argc, char** argv) {
	positional_arguments_t positional_arguments;
	named_arguments_t named_arguments;

	for (int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);

		if (arg[0] != '-' || arg[1] != '-') {
			positional_arguments.push_back(argv[i]);
			continue;
		}

		// We have a named argument, check if we should have extra data attached to it
		arg = arg.substr(2);

		// For the help command, there is no additional data, and we can stop parsing
		if (arg == "help") {
			named_arguments.insert(std::make_pair(arg, ""));
			break;
		}

		// Also components and undirected has no data attached to it.
		// TODO: How to make this nicer without a lot of overhead
		if (arg == "components" || arg == "undirected") {
			named_arguments.insert(std::make_pair(arg, "true"));
			continue;
		}

		// We have extra data, so parse it
		named_arguments.insert(std::make_pair(arg, argv[++i]));
	}

	return std::make_pair(positional_arguments, named_arguments);
}
