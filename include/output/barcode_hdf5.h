#pragma once

#ifdef WITH_HDF5
#include <iostream>

#include "../argparser.h"
#include "../definitions.h"
#include "../input/flagser.h"
#include "../parameters.h"
#include "hdf5_helper.h"

#include <hdf5.h>
#include <regex>

template <typename Complex> class barcode_hdf5_output_t : public output_t<Complex> {
	hid_t file_id;
	hid_t group_id;

	std::vector<hid_t> barcode_datasets;
	std::vector<hid_t> skipped_datasets;
	std::vector<value_t> barcode_buffer;
	std::vector<value_t> skipped_buffer;

	Complex* complex;
	unsigned short current_dimension = 0;
	unsigned short min_dimension;
	unsigned short max_dimension;
	coefficient_t modulus;
	std::vector<size_t> betti;
	std::vector<size_t> skipped;
	index_t euler_characteristic = 0;
	bool approximate_computation;
	bool aggregate_results = false;
	bool output_bars = true;

	static std::vector<size_t> total_betti;
	static std::vector<size_t> total_skipped;
	static std::vector<size_t> total_cell_count;
	static size_t total_top_dimension;

public:
	barcode_hdf5_output_t(const flagser_parameters& params)
	    : min_dimension(params.min_dimension), max_dimension(params.max_dimension), modulus(params.modulus),
	      approximate_computation(params.approximate_computation),
	      aggregate_results(!params.split_into_connected_components),
	      output_bars(params.filtration_algorithm.get() != nullptr) {
		const auto ids = open_or_create_group(params.output_name);
		file_id = ids.first;
		group_id = ids.second;
		barcode_buffer.reserve(BUFFER_SIZE + 50);
		skipped_buffer.reserve(BUFFER_SIZE + 50);
	}
	~barcode_hdf5_output_t() {
		for (auto i = 0ul; i < barcode_datasets.size(); i++) flush_if_necessary(true);
		for (auto id : barcode_datasets) H5Dclose(id);
		for (auto id : skipped_datasets) H5Dclose(id);
		if (group_id != file_id) H5Gclose(group_id);
		H5Fclose(file_id);
	}

	void set_complex(Complex* _complex) override { complex = _complex; }

	void finished(bool with_cell_counts = true) override {
		flush_if_necessary(true);
		total_top_dimension = std::max(total_top_dimension, complex->top_dimension());

		if (with_cell_counts) {
			std::vector<size_t> number_of_cells;
			for (auto i = 0ul; i <= complex->top_dimension(); i++)
				number_of_cells.push_back(complex->number_of_cells(i));
			// print_ordinary(file_output_t<Complex>::outstream, min_dimension, max_dimension, complex->top_dimension(),
			//  number_of_cells, betti, skipped, approximate_computation);

			total_cell_count.resize(complex->top_dimension() + 1, 0);
			for (auto i = 0ul; i <= complex->top_dimension(); i++) total_cell_count[i] += complex->number_of_cells(i);
		} else {
			total_cell_count.resize(0);
		}

		total_betti.resize(betti.size(), 0);
		for (size_t idx = 0; idx < betti.size(); idx++) total_betti[idx] += betti[idx];

		total_skipped.resize(skipped.size(), 0);
		for (size_t idx = 0; idx < skipped.size(); idx++) total_skipped[idx] += skipped[idx];
	}

	virtual void print_aggregated_results() override {
		bool computed_full_homology = total_cell_count.size() > 0 && min_dimension == 0 &&
		                              max_dimension == std::numeric_limits<unsigned short>::max();

		int first_dimension = std::max(0, min_dimension - 1);
		int last_dimension = total_top_dimension;

		const auto info_group_id = H5Gcreate2(group_id, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		{
			// Write dimension information
			hsize_t dims[1]{2};
			hsize_t max_dims[1]{2};
			const auto dataspace = H5Screate_simple(1, dims, max_dims);
			const auto dataset = H5Dcreate2(info_group_id, "dimension_range", H5T_NATIVE_INT, dataspace, H5P_DEFAULT,
			                                H5P_DEFAULT, H5P_DEFAULT);
			int data[2]{first_dimension, last_dimension};
			H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, data);
			H5Sclose(dataspace);
			H5Dclose(dataset);
		}

		if (total_cell_count.size() > 0) {
			// Write cell counts
			hsize_t dims[1]{static_cast<hsize_t>(last_dimension - first_dimension + 1)};
			hsize_t max_dims[1]{static_cast<hsize_t>(last_dimension - first_dimension + 1)};
			const auto dataspace = H5Screate_simple(1, dims, max_dims);
			const auto dataset = H5Dcreate2(info_group_id, "cell_counts", H5T_NATIVE_INT, dataspace, H5P_DEFAULT,
			                                H5P_DEFAULT, H5P_DEFAULT);
			std::vector<int> data(last_dimension - first_dimension + 1);
			for (int i = first_dimension; i <= last_dimension; i++) data[i] = total_cell_count[i];
			H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, data.data());
			H5Sclose(dataspace);
			H5Dclose(dataset);
		}

		{
			// Write betti numbers
			const hsize_t size = total_betti.size() - first_dimension;
			hsize_t dims[1]{size};
			hsize_t max_dims[1]{size};
			const auto dataspace = H5Screate_simple(1, dims, max_dims);
			const auto dataset = H5Dcreate2(info_group_id, "betti_numbers", H5T_NATIVE_INT, dataspace, H5P_DEFAULT,
			                                H5P_DEFAULT, H5P_DEFAULT);
			std::vector<int> data(last_dimension - first_dimension + 1);
			for (auto i = size_t(first_dimension); i < total_betti.size(); i++) data[i] = total_betti[i];
			H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, data.data());
			H5Sclose(dataspace);
			H5Dclose(dataset);
		}

		if (computed_full_homology) {
			index_t cell_euler_characteristic = 0;
			for (size_t i = 0; i <= total_top_dimension; i++)
				cell_euler_characteristic += (i % 2 == 1 ? -1 : 1) * total_cell_count[i];

			// Write euler characteristic
			hsize_t dims[1]{1};
			hsize_t max_dims[1]{1};
			const auto dataspace = H5Screate_simple(1, dims, max_dims);
			const auto dataset = H5Dcreate2(info_group_id, "euler_characteristic", H5T_NATIVE_INT, dataspace,
			                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			int data[1]{cell_euler_characteristic};
			H5Dwrite(dataset, H5T_NATIVE_INT, dataspace, dataspace, H5P_DEFAULT, data);
			H5Sclose(dataspace);
			H5Dclose(dataset);
		}
	}

	void computing_barcodes_in_dimension(unsigned short dimension) override {
		if (output_bars) {

			if (dimension > min_dimension) flush_if_necessary(true);

			const size_t index = dimension - std::max(0, min_dimension - 1);
			while (index >= barcode_datasets.size()) {
				std::string s = "dimension_";
				s += std::to_string(dimension);
				const auto current_dim_group_id =
				    H5Gcreate2(group_id, s.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

				{
					// Create an unlimited dataspace for the barcodes
					hsize_t dims[2]{0, 2};
					hsize_t max_dims[2]{H5S_UNLIMITED, 2};
					hsize_t chunk_dims[2]{BUFFER_SIZE, 2};
					const auto dataspace = H5Screate_simple(2, dims, max_dims);
					auto prop = H5Pcreate(H5P_DATASET_CREATE);
					H5Pset_chunk(prop, 2, chunk_dims);
					barcode_datasets.push_back(H5Dcreate2(current_dim_group_id, "barcodes", H5T_NATIVE_FLOAT, dataspace,
					                                      H5P_DEFAULT, prop, H5P_DEFAULT));
					H5Pclose(prop);
					H5Sclose(dataspace);
				}

				{
					// Create an unlimited dataspace for the skipped columns
					hsize_t dims[2]{0, 1};
					hsize_t max_dims[2]{H5S_UNLIMITED, 1};
					hsize_t chunk_dims[2]{BUFFER_SIZE, 1};
					const auto dataspace = H5Screate_simple(2, dims, max_dims);
					auto prop = H5Pcreate(H5P_DATASET_CREATE);
					H5Pset_chunk(prop, 2, chunk_dims);
					skipped_datasets.push_back(H5Dcreate2(current_dim_group_id, "skipped", H5T_NATIVE_FLOAT, dataspace,
					                                      H5P_DEFAULT, prop, H5P_DEFAULT));
					H5Pclose(prop);
					H5Sclose(dataspace);
				}
			}
		}

		current_dimension = dimension;
	}
	inline void new_barcode(value_t birth, value_t death) override {
		if (!output_bars) return;

		barcode_buffer.push_back(birth);
		barcode_buffer.push_back(death);

		flush_if_necessary();
	}
	inline void new_infinite_barcode(value_t birth) override {
		if (!output_bars) return;

		barcode_buffer.push_back(birth);
		barcode_buffer.push_back(std::numeric_limits<value_t>::infinity());

		flush_if_necessary();
	}
	inline void skipped_column(value_t birth) override {
		if (!output_bars) return;

		skipped_buffer.push_back(birth);

		flush_if_necessary();
	}
	void betti_number(size_t _betti, size_t _skipped) override {
		const auto shifted_dimension = current_dimension - min_dimension;
		betti.resize(shifted_dimension + 1, 0);
		skipped.resize(shifted_dimension + 1, 0);
		betti[shifted_dimension] = _betti;
		skipped[shifted_dimension] = _skipped;
		euler_characteristic += (shifted_dimension & 1 ? -1 : 1) * _betti;
	}
	void remaining_homology_is_trivial() override {}

private:
	inline size_t current_dataset_index() { return current_dimension - std::max(0, min_dimension - 1); }

	void flush_if_necessary(bool force = false) {
		if (!output_bars) return;

		if (force || barcode_buffer.size() >= 2 * BUFFER_SIZE) {
			// Write data
			size_t new_size = barcode_buffer.size() / 2;
			const auto dataset = barcode_datasets[current_dataset_index()];
			hsize_t dims[2];
			dims[0] = new_size;

			H5Dset_extent(dataset, dims);
			auto filespace = H5Dget_space(dataset);
			hsize_t offset[2]{dims[0] - new_size, 0};
			hsize_t dimsext[2]{new_size, 2};
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dimsext, NULL);
			auto memspace = H5Screate_simple(2, dimsext, NULL);
			H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, &barcode_buffer[0]);
			H5Sclose(memspace);
			H5Sclose(filespace);

			barcode_buffer.clear();
		}

		if (force || skipped_buffer.size() >= BUFFER_SIZE) {
			// Write data
			size_t new_size = skipped_buffer.size();
			const auto dataset = skipped_datasets[current_dataset_index()];

			hsize_t dims[2];
			dims[0] += new_size;
			H5Dset_extent(dataset, dims);
			auto filespace = H5Dget_space(dataset);
			hsize_t offset[2]{dims[0] - new_size, 0};
			hsize_t dimsext[2]{new_size, 1};
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dimsext, NULL);
			auto memspace = H5Screate_simple(2, dimsext, NULL);
			H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, &skipped_buffer[0]);
			H5Sclose(memspace);
			H5Sclose(filespace);

			skipped_buffer.clear();
		}
	}
};

template <typename Complex> std::vector<size_t> barcode_hdf5_output_t<Complex>::total_betti;
template <typename Complex> std::vector<size_t> barcode_hdf5_output_t<Complex>::total_skipped;
template <typename Complex> std::vector<size_t> barcode_hdf5_output_t<Complex>::total_cell_count;
template <typename Complex> size_t barcode_hdf5_output_t<Complex>::total_top_dimension = 0;

#endif
