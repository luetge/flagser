#pragma once

#include "../definitions.h"
#include "../input/flagser.h"

#include <hdf5.h>

const int BUFFER_SIZE = 10000;
#ifndef MANY_VERTICES
const auto VERTEX_INDEX_TYPE = H5T_NATIVE_USHORT;
#else
const auto VERTEX_INDEX_TYPE = H5T_NATIVE_UINT32;
#endif

std::pair<hid_t, hid_t> open_or_create_group(const std::string& filename) {
	size_t h5_pos = filename.rfind(".h5");
	std::string fname = filename;
	std::string h5_path = "/";
	if (h5_pos != std::string::npos && h5_pos < filename.size() - 3) {
		fname = filename.substr(0, h5_pos + 3);
		h5_path = filename.substr(h5_pos + 3);
	}

	hid_t file_id;
	std::ifstream f(fname);
	if (f.good()) {
		f.close();
		file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	} else {
		f.close();
		file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	}

	if (file_id == -1) {
		std::string err_msg = "The output file " + fname + " could not be opened.";
		throw std::invalid_argument(err_msg);
	}

	// Suppress errors
	H5E_auto_t old_func;
	void* old_client_data;
	H5Eget_auto2(H5P_DEFAULT, &old_func, &old_client_data);
	H5Eset_auto2(H5P_DEFAULT, NULL, NULL);

	hid_t group_id = file_id;
	std::vector<std::string> groups = split<std::string>(h5_path, '/', [](std::string s) { return s; });
	for (auto group : groups) {
		if (group.size() == 0) continue;
		hid_t new_group_id = H5Gopen2(group_id, group.c_str(), H5P_DEFAULT);
		if (new_group_id < 0) new_group_id = H5Gcreate2(group_id, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (new_group_id < 0) {
			std::string err_msg = "Could not create group " + group + ".";
			throw std::invalid_argument(err_msg);
		}
		group_id = new_group_id;
	}

	if (group_id >= 0)
		while (H5Ldelete_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, 0, H5P_DEFAULT) >= 0) {}

	// Show errors again
	H5Eset_auto2(H5P_DEFAULT, old_func, old_client_data);

	if (group_id < 0) {
		std::string err_msg = "The data could not be written to " + h5_path + ".";
		throw std::invalid_argument(err_msg);
	}

	return std::make_pair(file_id, group_id);
}
