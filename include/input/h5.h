#pragma once

#include "base.h"

#ifndef WITH_HDF5

filtered_directed_graph_t read_graph_h5(const std::string, const flagser_parameters&) {
	std::string err_msg = "Error: flagser was compiled without support for .h5-files. Please install the HDF5-library "
	                      "(https://support.hdfgroup.org/HDF5/) and rebuild flagser by running \"make\" again.";
	throw std::runtime_error(err_msg);
}

#else

#include <hdf5.h>
#include <regex>
#include <string>

#include "../parameters.h"

std::string replace_string(std::string subject, const std::string& search, const std::string& replace) {
	size_t pos = 0;
	while ((pos = subject.find(search, pos)) != std::string::npos) {
		subject.replace(pos, search.length(), replace);
		pos += replace.length();
	}
	return subject;
}

struct group_extractor_t {
	group_extractor_t(std::string groups) {
		std::vector<std::string> gs = split<std::string>(groups, ',', [](std::string s) { return s; });
		for (auto g : gs) {
			if (g == "_____all") {
				all_groups = true;
				break;
			}

			regexes.push_back(std::regex(replace_string(trim(g), "*", ".*")));
		}
	}

	herr_t operator()(const char* name) {
		if (!all_groups) {
			bool matches = false;
			for (auto reg : regexes) {
				if (std::regex_match(name, reg)) {
					matches = true;
					break;
				}
			}

			if (!matches) return 0;
		}

		group_names.push_back(name);
		return 0;
	}

	std::vector<std::string> group_names;

private:
	bool all_groups = false;
	std::vector<std::regex> regexes;
};

herr_t extract_groups(hid_t, const char* name, const H5L_info_t*, void* opdata) {
	group_extractor_t* extractor = (group_extractor_t*)opdata;
	return (*extractor)(name);
}

const filtered_directed_graph_t read_graph_h5(const std::string filename, const flagser_parameters& params) {
	filtered_directed_graph_t* graph;
	std::vector<value_t> vertex_filtration;
	std::string type = params.hdf5_type;
	bool directed = params.directed;
	bool with_filtration = params.filtration_algorithm.get() != nullptr;

	size_t h5_pos = filename.rfind(".h5");
	std::string fname = filename;
	std::string h5_path = "/";
	if (h5_pos != std::string::npos && h5_pos < filename.size() - 3) {
		fname = filename.substr(0, h5_pos + 3);
		h5_path = filename.substr(h5_pos + 3);
	}

	if (!H5Fis_hdf5(fname.c_str())) {
		std::string err_msg = "The file \"" + fname + "\" is not an HDF5-file. Please choose another input format.";
		throw std::invalid_argument(err_msg);
	}

	auto file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	if (type == "matrix") {
		auto dataset_id = H5Dopen(file_id, h5_path.c_str(), H5P_DEFAULT);

		if (dataset_id == -1) {
			std::string err_msg =
			    "\n\nThe connectivity data could not be found at \"" + h5_path + "\", please check the path again.";
			throw std::runtime_error(err_msg);
		}

		// Read the number of vertices
		size_t vertex_number;
		{
			auto space = H5Dget_space(dataset_id);
			hsize_t dim[2];
			hsize_t maxdim[2];
			H5Sget_simple_extent_dims(space, dim, maxdim);
			vertex_number = dim[0];
			H5Sclose(space);
		}

		// We now know the number of vertices
		graph = new filtered_directed_graph_t(std::vector<value_t>(vertex_number), directed);

		if (with_filtration) {
			value_t* entries = new value_t[vertex_number * vertex_number];
			H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, entries);

			// Add the edges
			for (auto x = 0ul; x < vertex_number; x++) {
				graph->vertex_filtration[x] = entries[x * vertex_number + x];

				for (auto y = 0ul; y < vertex_number; y++) {
					if (entries[x * vertex_number + y] != 0)
						graph->add_filtered_edge(x, y, entries[x * vertex_number + y]);
				}
			}
		} else {
			hid_t enumtype = H5Tcreate(H5T_ENUM, sizeof(index_t));
			index_t val = 0;
			H5Tenum_insert(enumtype, "FALSE", &val);
			val = 1;
			H5Tenum_insert(enumtype, "TRUE", &val);

			index_t* entries = new index_t[vertex_number * vertex_number];
			auto type = H5Dget_type(dataset_id);
			H5Dread(dataset_id, H5Tget_class(type) == H5T_ENUM ? enumtype : H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
			        H5P_DEFAULT, entries);
			H5Tclose(type);

			// Add the edges
			for (auto x = 0ul; x < vertex_number; x++) {
				for (auto y = 0ul; y < vertex_number; y++) {
					if (entries[x * vertex_number + y] != 0) graph->add_edge(x, y);
				}
			}

			delete[] entries;
			H5Tclose(enumtype);
		}

		H5Dclose(dataset_id);
	} else {
		// type == "grouped" or type == "grouped:L1*,..."
		std::string groups = type.length() >= 8 ? type.substr(8) : "";
		if (groups == "") groups = "_____all";

		auto connectivity_id = H5Gopen(file_id, h5_path.c_str(), H5P_DEFAULT);
		if (connectivity_id == -1) {
			std::string err_msg =
			    "\n\nThe connectivity data could not be found at \"" + h5_path + "\", please check the path again.";
			throw std::runtime_error(err_msg);
		}

		group_extractor_t group_extractor(groups);
		H5Literate(connectivity_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, &extract_groups, &group_extractor);

		std::vector<unsigned int> group_offsets;
		std::vector<unsigned int> group_sizes;
		group_offsets.reserve(group_extractor.group_names.size());
		unsigned int offset = 0;
		for (auto i = 0ul; i < group_extractor.group_names.size(); i++) {
			std::string name = group_extractor.group_names[i];
			name += "/";
			name += group_extractor.group_names[i];
			name += "/cMat";
			auto dataset_id = H5Dopen(connectivity_id, name.c_str(), H5P_DEFAULT);

			if (dataset_id == -1) {
				std::string err_msg =
				    "\n\nThe connectivity data could not be found at \"" + h5_path + "\", please check the path again.";
				throw std::runtime_error(err_msg);
			}

			auto space = H5Dget_space(dataset_id);
			hsize_t dim[2];
			hsize_t maxdim[2];
			H5Sget_simple_extent_dims(space, dim, maxdim);
			group_offsets.push_back(offset);
			group_sizes.push_back(dim[0]);
			offset += dim[0];
			H5Sclose(space);
			H5Dclose(dataset_id);
		}

		// We now know the number of vertices
		graph = new filtered_directed_graph_t(std::vector<value_t>(offset), directed);

		hid_t enumtype = H5Tcreate(H5T_ENUM, sizeof(index_t));
		index_t val = 0;
		H5Tenum_insert(enumtype, "FALSE", &val);
		val = 1;
		H5Tenum_insert(enumtype, "TRUE", &val);

		for (auto i = 0ul; i < group_extractor.group_names.size(); i++) {
			for (auto j = 0ul; j < group_extractor.group_names.size(); j++) {
				std::string name = group_extractor.group_names[i];
				name += "/";
				name += group_extractor.group_names[j];
				name += "/cMat";
				auto dataset_id = H5Dopen(connectivity_id, name.c_str(), H5P_DEFAULT);

				if (with_filtration) {
					value_t* entries = new value_t[group_sizes[i] * group_sizes[j]];
					H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, entries);

					// Add the edges
					for (auto x = 0ul; x < group_sizes[i]; x++) {
						graph->vertex_filtration[group_offsets[i] + x] = entries[x * group_sizes[j] + x];

						for (auto y = 0ul; y < group_sizes[j]; y++) {
							if (entries[x * group_sizes[j] + y] > 0) {
								graph->add_filtered_edge(group_offsets[i] + x, group_offsets[j] + y,
								                         entries[x * group_sizes[j] + y]);
							}
						}
					}
					delete[] entries;
				} else {
					index_t* entries = new index_t[group_sizes[i] * group_sizes[j]];
					auto type = H5Dget_type(dataset_id);
					H5Dread(dataset_id, H5Tget_class(type) == H5T_ENUM ? enumtype : H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
					        H5P_DEFAULT, entries);

					// Add the edges
					for (auto x = 0ul; x < group_sizes[i]; x++) {
						for (auto y = 0ul; y < group_sizes[j]; y++) {
							if (entries[x * group_sizes[j] + y] == 1) {
								graph->add_edge(group_offsets[i] + x, group_offsets[j] + y);
							}
						}
					}

					H5Tclose(type);
					delete[] entries;
				}
				H5Dclose(dataset_id);
			}
		}

		H5Tclose(enumtype);

		H5Gclose(connectivity_id);
	}

	H5Fclose(file_id);

	return *graph;
}
#endif
