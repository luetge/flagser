#pragma once

#ifdef WITH_HDF5
#include <iostream>

#include "../argparser.h"
#include "../definitions.h"
#include "../input/flagser.h"
#include "../parameters.h"
#include "./hdf5_helper.h"

#include <hdf5.h>
#include <regex>

class hdf5_output_t {
	hid_t file_id;
	hid_t group_id;

	std::vector<hid_t> datasets;
	std::vector<size_t> sizes;
	std::vector<std::vector<vertex_index_t>> buffer;

public:
	hdf5_output_t(const flagser_parameters& params) {
		auto ids = open_or_create_group(params.output_name);
		file_id = ids.first;
		group_id = ids.second;
	}
	~hdf5_output_t() {
		for (auto i = 0ul; i < datasets.size(); i++) flush(i + 2);
		for (auto id : datasets) H5Dclose(id);
		if (group_id != file_id) H5Gclose(group_id);
		H5Fclose(file_id);
	}

	void write_cell(vertex_index_t* first_vertex, int size) {
		if (size < 2) return;

		if (size - 2 >= long(datasets.size())) prepare(size);

		vertex_index_t* v = first_vertex;
		for (int i = 0; i < size; ++i, ++v) buffer[size - 2].push_back(*v);

		if (buffer[size - 2].size() == size_t(BUFFER_SIZE * size)) flush(size);
	}

private:
	void flush(int size) {
		int idx = size - 2;
		// Write data
		size_t new_size = buffer[idx].size() / size;
		hsize_t dims[2]{sizes[idx] + new_size, (hsize_t)size};
		H5Dset_extent(datasets[idx], dims);
		auto filespace = H5Dget_space(datasets[idx]);
		hsize_t offset[2]{sizes[idx], 0};
		hsize_t dimsext[2]{new_size, (hsize_t)size};
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dimsext, NULL);
		auto memspace = H5Screate_simple(2, dimsext, NULL);
		H5Dwrite(datasets[idx], VERTEX_INDEX_TYPE, memspace, filespace, H5P_DEFAULT, &buffer[idx][0]);
		H5Sclose(memspace);
		H5Sclose(filespace);

		sizes[idx] += buffer[idx].size() / size;
		buffer[idx].clear();
	}

	inline void prepare(int size) {
		while (size - 2 >= long(datasets.size())) {
			hsize_t dim = datasets.size() + 1;
			// Create an unlimited dataspace
			hsize_t dims[2]{0, 0};
			hsize_t max_dims[2]{H5S_UNLIMITED, H5S_UNLIMITED};
			hsize_t chunk_dims[2]{BUFFER_SIZE, dim + 1};
			auto dataspace = H5Screate_simple(2, dims, max_dims);
			auto prop = H5Pcreate(H5P_DATASET_CREATE);
			H5Pset_chunk(prop, 2, chunk_dims);
			std::string s = "Cells_";
			s += std::to_string(dim);
			datasets.push_back(
			    H5Dcreate2(group_id, s.c_str(), VERTEX_INDEX_TYPE, dataspace, H5P_DEFAULT, prop, H5P_DEFAULT));
			H5Pclose(prop);
			H5Sclose(dataspace);

			std::vector<vertex_index_t> b;
			b.reserve(BUFFER_SIZE + 50);
			buffer.push_back(b);
		}

		sizes.resize(size - 1, 0);
	}
};

#endif
