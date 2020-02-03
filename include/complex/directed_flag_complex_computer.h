#pragma once

#define USE_CELLS_WITHOUT_DIMENSION

#include <algorithm>
#include <array>
#include <memory>

#include "../argparser.h"
#include "../directed_graph.h"
#include "../filtration_algorithms.h"
#include "../persistence.h"
#include "directed_flag_complex.h"

typedef hash_map<directed_flag_complex_cell_t, size_t, cell_hasher_t, cell_comparer_t> cell_hash_map_t;

template <typename Complex> class coboundary_iterator_t {
	const Complex* complex;
	short dimension;
	const compressed_sparse_matrix<entry_t>& matrix;
	const index_t index;
	const coefficient_t coefficient;
	const coefficient_t modulus;
	size_t current_offset = 0;

public:
	coboundary_iterator_t(const Complex* _complex, short _dimension, const compressed_sparse_matrix<entry_t>& _matrix,
	                      index_t _index, coefficient_t _coefficient, coefficient_t _modulus)
	    : complex(_complex), dimension(_dimension), matrix(_matrix), index(_index), coefficient(_coefficient),
	      modulus(_modulus) {}

	bool has_next() { return index != -1 && matrix.cbegin(index) + current_offset != matrix.cend(index); }

	filtration_entry_t next() {
		entry_t entry = *(matrix.cbegin(index) + current_offset++);
		coefficient_t coface_coefficient = get_coefficient(entry) * coefficient % modulus;
		return filtration_entry_t(complex->filtration(dimension + 1, get_index(entry)), get_index(entry),
		                          coface_coefficient);
	}
};

struct reorder_edges_t {
	reorder_edges_t(std::vector<vertex_index_t>& _new_edges, bool _reorder_filtration,
	                hash_map<directed_flag_complex_cell_t, size_t, cell_hasher_t, cell_comparer_t>& _cell_hash,
	                std::vector<value_t>& _old_filtration, std::vector<value_t>& _new_filtration)
	    : cell_hash(_cell_hash), old_filtration(_old_filtration), new_filtration(_new_filtration),
	      new_edges(_new_edges), reorder_filtration(_reorder_filtration) {}
	void done() {}
	void operator()(vertex_index_t* first_vertex, int) {
		directed_flag_complex_cell_t cell(first_vertex);
		if (reorder_filtration) new_filtration.push_back(old_filtration[cell_hash.find(cell)->second]);
		new_edges.push_back(cell.vertex(0));
		new_edges.push_back(cell.vertex(1));
	}

private:
	hash_map<directed_flag_complex_cell_t, size_t, cell_hasher_t, cell_comparer_t>& cell_hash;
	std::vector<value_t>& old_filtration;
	std::vector<value_t>& new_filtration;
	std::vector<vertex_index_t>& new_edges;
	bool reorder_filtration;
};

struct compute_filtration_t {
	compute_filtration_t(filtration_algorithm_t* _filtration_algorithm, const filtered_directed_graph_t& _graph,
	                     const std::vector<value_t>& _current_filtration, std::vector<cell_hash_map_t>& _cell_hash,
	                     size_t* _cell_hash_offsets = nullptr)
	    : graph(_graph), filtration_algorithm(_filtration_algorithm), current_filtration(_current_filtration),
	      cell_hash(_cell_hash), cell_hash_offsets(_cell_hash_offsets) {}
	void done() {}

	void operator()(vertex_index_t* first_vertex, int size) {
		// The index starts at -1, so the first cell gets index 0
		current_index++;
		directed_flag_complex_cell_t cell(first_vertex);
		if (cell_hasher_t::dimension != size - 2) cell_hasher_t::set_current_cell_dimension(size - 2);

		if (filtration_algorithm == nullptr) return;

		if (filtration_algorithm->needs_face_filtration()) {
			if (boundary_filtration == nullptr) boundary_filtration = new value_t[size];

			for (int i = 0; i < size; i++) {
				auto bdry = cell.boundary(i);

				if (!cell_hash.size()) {
					boundary_filtration[i] = current_filtration[bdry.vertex(i)];
				} else {
					// The threads are split by the first vertex
					short thread_index = bdry.vertex(0) % PARALLEL_THREADS;
					auto pair = cell_hash[thread_index].find(bdry);

					if (pair == cell_hash[thread_index].end()) {
						std::string err = "Could not find boundary ";
						err += cell.boundary(i).to_string(size - 2);
						err += " of ";
						err += cell.to_string(size - 1);
						err += ".\n";
						std::cerr << err;
						exit(-1);
					}
					boundary_filtration[i] = current_filtration[pair->second + cell_hash_offsets[thread_index]];
				}
			}

			next_filtration.push_back(
			    filtration_algorithm->compute_filtration(size - 1, cell, graph, boundary_filtration));
		} else {
			next_filtration.push_back(filtration_algorithm->compute_filtration(size - 1, cell, graph, nullptr));
		}
	}

	std::vector<value_t> filtration() const { return next_filtration; }

	index_t number_of_cells() const { return current_index + 1; }

private:
	index_t current_index = -1;
	const filtered_directed_graph_t& graph;

	value_t* boundary_filtration = nullptr;
	filtration_algorithm_t* filtration_algorithm;
	const std::vector<value_t>& current_filtration;
	std::vector<value_t> next_filtration;
	std::vector<cell_hash_map_t>& cell_hash;
	size_t* cell_hash_offsets;
};

template <typename Complex>
void prepare_graph_filtration(Complex& complex, filtered_directed_graph_t& graph,
                              filtration_algorithm_t* filtration_algorithm) {

	if (filtration_algorithm == nullptr) {
		// All vertices get the trivial filtration value
		graph.vertex_filtration = std::vector<value_t>(graph.vertex_filtration.size(), 0);
	} else {
		if (filtration_algorithm->overwrite_vertex_filtration()) {
			directed_flag_complex_cell_t cell;
			for (vertex_index_t v = 0; v < graph.number_of_vertices; v++) {
				cell.set_vertices(&v);
				graph.vertex_filtration[v] = filtration_algorithm->compute_filtration(0, cell, graph, nullptr);
			}
		}
	}

	bool computed_edge_filtration = false;
	if (filtration_algorithm != nullptr &&
	    (graph.edge_filtration.size() == 0 || filtration_algorithm->overwrite_edge_filtration())) {
		computed_edge_filtration = true;
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "computing the filtration of all edges" << std::flush << "\r";
#endif

        // Dummy parameter because we pass by reference the argument and
        // It cannot take a default value
        std::vector<cell_hash_map_t> tmp(0);
		std::vector<compute_filtration_t> compute_filtration(
		    PARALLEL_THREADS, compute_filtration_t(filtration_algorithm, graph, graph.vertex_filtration, tmp));
		complex.for_each_cell(compute_filtration, 1);
		graph.edge_filtration.clear();
		graph.edge_filtration.reserve(graph.edge_number());
		for (int i = 0; i < PARALLEL_THREADS; i++) {
			for (auto f : compute_filtration[i].filtration()) graph.edge_filtration.push_back(f);
		}
	}

	// Now, we reorder the edges and their filtration value such that it matches
	// the iteration order of "for_each_cell"
	std::vector<vertex_index_t> new_edges[PARALLEL_THREADS];
	std::vector<value_t> new_filtrations[PARALLEL_THREADS];

	hash_map<directed_flag_complex_cell_t, size_t, cell_hasher_t, cell_comparer_t> hash;
	if (filtration_algorithm != nullptr) {
		cell_hasher_t::set_current_cell_dimension(1);
		vertex_index_t* e = &(graph.edges[0]);
		directed_flag_complex_cell_t c;
		for (size_t i = 0; i < graph.edges.size() / 2; ++i, e += 2) {
			c.set_vertices(e);
			hash.insert(std::make_pair(c, i));
		}
	}

	std::vector<reorder_edges_t> reorder_filtration;
	for (int i = 0; i < PARALLEL_THREADS; i++) {
		reorder_filtration.push_back(reorder_edges_t(new_edges[i],
		                                             !computed_edge_filtration && filtration_algorithm != nullptr, hash,
		                                             graph.edge_filtration, new_filtrations[i]));
	}
	complex.for_each_cell(reorder_filtration, 1);

	hash.clear();
	hash_map<directed_flag_complex_cell_t, size_t, cell_hasher_t, cell_comparer_t> empty;
	std::swap(hash, empty);

	size_t current_filtration_index = 0;
	size_t current_edge_index = 0;
	for (int i = 0; i < PARALLEL_THREADS; i++) {
		if (!computed_edge_filtration)
			for (auto f : new_filtrations[i]) graph.edge_filtration[current_filtration_index++] = f;
		for (auto e : new_edges[i]) graph.edges[current_edge_index++] = e;
	}
}

class directed_flag_complex_computer_t {
	filtered_directed_graph_t& graph;
	std::unique_ptr<filtration_algorithm_t> filtration_algorithm;
	unsigned short max_dimension;
	unsigned short min_dimension;
	int current_dimension = 0;
	bool _is_top_dimension = false;
	std::vector<size_t> cell_count;
	std::string cache;

	directed_flag_complex_t flag_complex;

	// Filtration
	std::vector<value_t> next_filtration;

	// Coboundaries
	compressed_sparse_matrix<entry_t> coboundary_matrix[PARALLEL_THREADS];
	size_t coboundary_matrix_offsets[PARALLEL_THREADS];
	coefficient_t modulus;

public:
	directed_flag_complex_computer_t(filtered_directed_graph_t& _graph, const named_arguments_t& named_arguments)
	    : graph(_graph),
	      filtration_algorithm(get_filtration_computer(get_argument_or_default(named_arguments, "filtration", "zero"))),
	      max_dimension(atoi(get_argument_or_default(named_arguments, "max-dim", "65535"))),
	      min_dimension(atoi(get_argument_or_default(named_arguments, "min-dim", "0"))),
	      cache(get_argument_or_default(named_arguments, "cache", "")),
        flag_complex(graph),
	      modulus(atoi(get_argument_or_default(named_arguments, "modulus", "2"))) {
		cell_count.push_back(_graph.vertex_number());
		cell_count.push_back(_graph.edge_number());

		// Order the edges and compute the correct filtration
		if (min_dimension <= 1 || (filtration_algorithm != nullptr && filtration_algorithm->needs_face_filtration()))
			prepare_graph_filtration(flag_complex, graph, filtration_algorithm.get());
	}

	size_t number_of_cells(int dimension) const {
		assert(dimension < cell_count.size());
		return cell_count[dimension];
	}
	size_t top_dimension() const { return cell_count.size() - 1; }

	const std::vector<value_t> vertex_filtration() const { return graph.vertex_filtration; }

	inline value_t filtration(int dimension, index_t index) const {
		// Only return something if we have a non-trivial filtration algorithm
		if (filtration_algorithm == nullptr) return 0.0f;

		if (dimension == 0) return graph.vertex_filtration[index];

		assert(dimension == current_dimension + 1);
		if (dimension == 1) return graph.edge_filtration[index];

		return next_filtration[index];
	}

	std::pair<vertex_index_t, vertex_index_t> vertices_of_edge(index_t edge) const {
		return std::make_pair(graph.edges[2 * edge], graph.edges[2 * edge + 1]);
	}

	// Note: Gets called with consecutive dimensions, starting with zero
	void prepare_next_dimension(int dimension);

	inline coboundary_iterator_t<directed_flag_complex_computer_t> coboundary(filtration_entry_t cell) {
		if (current_dimension > max_dimension || is_top_dimension()) {
			return coboundary_iterator_t<directed_flag_complex_computer_t>(this, current_dimension,
			                                                               coboundary_matrix[0], -1, 1, modulus);
		}

		int i = 0;
		while (i < PARALLEL_THREADS - 1 && index_t(coboundary_matrix_offsets[i + 1]) <= get_index(cell)) { i++; }
		return coboundary_iterator_t<directed_flag_complex_computer_t>(this, current_dimension, coboundary_matrix[i],
		                                                               get_index(cell) - coboundary_matrix_offsets[i],
		                                                               get_coefficient(cell), modulus);
	}

	bool is_top_dimension() { return _is_top_dimension; }

	void computation_result(int, size_t, size_t) {}
	void finished() {}
};

struct add_cell_index_to_cache_t {
	add_cell_index_to_cache_t(cell_hash_map_t* _cache) : cache(_cache) {}
	~add_cell_index_to_cache_t() {
		// Cleanup memory
		for (auto p : *cache) delete[] p.first.get_vertex_array();
	}
	void done() {}
	void operator()(vertex_index_t* first_vertex, int size) {
		// The index starts at -1, so the first cell gets index 0
		current_index++;

		// Save the cell vertices on the heap
		vertex_index_t* c = new vertex_index_t[size];
		std::copy(first_vertex, first_vertex + size, c);

		directed_flag_complex_cell_t cell(c);
		cell_hasher_t::set_current_cell_dimension(size - 1);

		cache->insert(std::make_pair(cell, current_index));
	}

	index_t number_of_cells() const { return current_index + 1; }

private:
	cell_hash_map_t* cache;
	index_t current_index = -1;
};

struct store_coboundaries_in_cache_t {
	store_coboundaries_in_cache_t(compressed_sparse_matrix<entry_t>& _coboundary_matrix, int _current_dimension,
	                              const filtered_directed_graph_t& _graph, std::vector<cell_hash_map_t>& _cell_hash,
	                              size_t* _cell_hash_offsets, size_t _total_cell_number, bool _is_first,
	                              coefficient_t _modulus = 2)
	    : is_first(_is_first), current_dimension(_current_dimension), coboundary_matrix(_coboundary_matrix),
	      graph(_graph), cell_hash(_cell_hash), cell_hash_offsets(_cell_hash_offsets),
	      total_cell_number(_total_cell_number), modulus(_modulus) {}
	void done() {
#ifdef INDICATE_PROGRESS
		if (is_first)
			std::cout << "\033[K"
			          << "dimension " << current_dimension << ": computed almost all of the coboundaries" << std::flush
			          << "\r";
#endif
	}
	void operator()(vertex_index_t* first_vertex, int size) {
#ifdef INDICATE_PROGRESS
		if (is_first && (current_index + 1) % 10000 == 0) {
			std::cout << "\033[K"
			          << "dimension " << current_dimension << ": computed ca. "
			          << PARALLEL_THREADS * (current_index + 1);
			if (total_cell_number > 0) std::cout << "/" << total_cell_number;
			std::cout << " coboundaries" << std::flush << "\r";
		}
#endif
		directed_flag_complex_cell_t cell(first_vertex);
		coboundary_matrix.append_column();

		std::vector<size_t> vertex_offsets;
		for (int j = 0; j < size; j++) vertex_offsets.push_back(first_vertex[j] >> 6);

		for (int i = 0; i <= size; i++) {
			// Check intersections in chunks of 64
			for (size_t offset = 0; offset < graph.incidence_row_length; offset++) {
				size_t bits = -1; // All bits set

				for (int j = 0; bits > 0 && j < size; j++) {
					// Remove the vertices already making up the cellk
					if (vertex_offsets[j] == offset) bits &= ~(ONE_ << (first_vertex[j] - (vertex_offsets[j] << 6)));

					// Intersect with the outgoing/incoming edges of the current vertex
					bits &= j < i ? graph.get_outgoing_chunk(first_vertex[j], offset)
					              : graph.get_incoming_chunk(first_vertex[j], offset);
				}

				size_t vertex_offset = offset << 6;
				while (bits > 0) {
					// Get the least significant non-zero bit
					int b = __builtin_ctzl(bits);

					// Unset this bit
					bits &= ~(ONE_ << b);

					// Now insert the appropriate vertex at this position
					const auto& cb = cell.insert_vertex(i, vertex_offset + b);
					short thread_index = cb.vertex(0) % PARALLEL_THREADS;
					auto pair = cell_hash[thread_index].find(cb);
					if (pair == cell_hash[thread_index].end()) {
						std::string err = "Could not find coboundary ";
						err += cb.to_string(current_dimension + 1);
						err += ".\n";
						std::cerr << err;
						exit(-1);
					}
					coboundary_matrix.push_back(
					    make_entry(pair->second + cell_hash_offsets[thread_index], i & 1 ? -1 + modulus : 1));
				}
			}
		}

		current_index++;
	}

private:
	index_t current_index = 0;
	bool is_first;
	int current_dimension;
	compressed_sparse_matrix<entry_t>& coboundary_matrix;
	const filtered_directed_graph_t& graph;
	std::vector<cell_hash_map_t>& cell_hash;
	size_t* cell_hash_offsets;
	size_t total_cell_number;
	coefficient_t modulus;
};

void directed_flag_complex_computer_t::prepare_next_dimension(int dimension) {
	if (dimension == 0) return;

#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	          << "preparing dimension " << dimension << ": indexing " << (dimension) << "-dimensional cells"
	          << std::flush << "\r";
#endif

	assert(dimension == current_dimension + 1);
	current_dimension = dimension;
	cell_hasher_t::set_current_cell_dimension(dimension + 1);
	cell_count.resize(dimension + 2);

	// Clean up
	for (auto p : coboundary_matrix) { p.clear(); }

	if (dimension > max_dimension || _is_top_dimension) return;

	if (cache != "") {
#ifdef USE_COEFFICIENTS
		// TODO: Make this work
		std::cerr << "Sorry, caching does not work with coefficients yet." << std::endl;
		exit(1);
#endif
		bool loaded_from_file = false;
		for (int i = 0; i < PARALLEL_THREADS; i++) {
			coboundary_matrix[i].clear();
			std::string fname = cache;
			fname += "/matrix_";
			fname += std::to_string(current_dimension);
			fname += "_";
			fname += std::to_string(i);
			std::ifstream f(fname.c_str(), std::ios::binary);
			if (i == 0 && !f.good()) break;
			loaded_from_file = true;

			f.read((char*)&(cell_count[current_dimension + 1]), sizeof(size_t));
			f.read((char*)&(coboundary_matrix_offsets[i]), sizeof(size_t));

			size_t this_size;
			f.read((char*)&(this_size), sizeof(size_t));

			index_t next_value;
			while (this_size > 0 && f.read((char*)&next_value, sizeof(index_t))) {
				if (next_value == -1)
					coboundary_matrix[i].append_column();
				else
					coboundary_matrix[i].push_back(next_value);
			}
			f.close();
		}

		std::string fname = cache;
		fname += "/filtration_";
		fname += std::to_string(current_dimension);
		std::ifstream f(fname.c_str(), std::ios::binary);
		next_filtration.clear();
		value_t next_value;
		while (f.read((char*)&next_value, sizeof(value_t))) next_filtration.push_back(next_value);
		f.close();

		if (loaded_from_file) {
			_is_top_dimension = cell_count[current_dimension + 1] == 0;
			return;
		}
	}

	{
		if (filtration_algorithm != nullptr) {
			// Add the current and next cells into the hash and compute the new
			// filtration
			std::vector<cell_hash_map_t> _cache_current_cells(PARALLEL_THREADS);
			size_t _cache_current_cells_offsets[PARALLEL_THREADS];
			std::vector<add_cell_index_to_cache_t> init_current_cell_cache;
			for (int i = 0; i < PARALLEL_THREADS; i++) {
				init_current_cell_cache.push_back(add_cell_index_to_cache_t(&_cache_current_cells[i]));
			}

			cell_hasher_t::set_current_cell_dimension(dimension);

			// Add the current cells into the cache if the filtration algorithm needs
			// them
			if (filtration_algorithm->needs_face_filtration())
				flag_complex.for_each_cell(init_current_cell_cache, dimension);

			// If we will actually compute coboundaries, then compute the filtration.
			// Also if we need the face filtrations.
			if (dimension + 1 >= min_dimension || filtration_algorithm->needs_face_filtration()) {
#ifdef INDICATE_PROGRESS
				std::cout << "\033[K"
				          << "preparing dimension " << dimension << ": computing the filtration of all "
				          << (dimension + 1) << "-dimensional cells" << std::flush << "\r";
#endif
				size_t offset = 0;
				std::vector<compute_filtration_t> compute_filtration;
				for (int i = 0; i < PARALLEL_THREADS; i++) {
					_cache_current_cells_offsets[i] = offset;
					if (filtration_algorithm->needs_face_filtration()) 
                    {
                        offset += _cache_current_cells[i].size();
                    }
					compute_filtration.push_back(compute_filtration_t(
					    filtration_algorithm.get(), graph, dimension == 1 ? graph.edge_filtration : next_filtration,
					    _cache_current_cells, _cache_current_cells_offsets));
				}
				flag_complex.for_each_cell(compute_filtration, dimension + 1);

				size_t _cell_count = 0;
				for (int i = 0; i < PARALLEL_THREADS; i++) _cell_count += compute_filtration[i].number_of_cells();
				cell_count[dimension + 1] = _cell_count;
				if (_cell_count == 0) _is_top_dimension = true;

				// Combine the filtration
				next_filtration.clear();
				next_filtration.reserve(_cell_count);
				for (int i = 0; i < PARALLEL_THREADS; i++) {
					for (auto f : compute_filtration[i].filtration()) next_filtration.push_back(f);
				}
			}
		}

		// If we do not want the homology in the next degree, then we can stop
		// here
		if (dimension + 1 < min_dimension) return;

		cell_hasher_t::set_current_cell_dimension(dimension + 1);

		std::vector<cell_hash_map_t> _cache_next_cells(PARALLEL_THREADS);
		size_t _cache_next_cells_offsets[PARALLEL_THREADS];
		std::vector<add_cell_index_to_cache_t> init_next_cell_cache;
		for (int i = 0; i < PARALLEL_THREADS; i++) {
			init_next_cell_cache.push_back(add_cell_index_to_cache_t(&_cache_next_cells[i]));
		}
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "preparing dimension " << dimension << ": indexing " << (dimension + 1) << "-dimensional cells"
		          << std::flush << "\r";
#endif

		size_t _cell_count = 0;

		if (!is_top_dimension()) {
			flag_complex.for_each_cell(init_next_cell_cache, dimension + 1);

			for (int i = 0; i < PARALLEL_THREADS; i++) {
				_cache_next_cells_offsets[i] = _cell_count;
				_cell_count += init_next_cell_cache[i].number_of_cells();
			}
			cell_count[dimension + 1] = _cell_count;
			if (_cell_count == 0) _is_top_dimension = true;
		}

		// Now compute the coboundaries
		std::vector<store_coboundaries_in_cache_t> store_coboundaries;
		for (int i = 0; i < PARALLEL_THREADS; i++) coboundary_matrix[i] = compressed_sparse_matrix<entry_t>();
		for (int i = 0; i < PARALLEL_THREADS; i++) {
			store_coboundaries.push_back(store_coboundaries_in_cache_t(coboundary_matrix[i], dimension, graph,
			                                                           _cache_next_cells, _cache_next_cells_offsets,
			                                                           cell_count[dimension], i == 0, modulus));
		}
		flag_complex.for_each_cell(store_coboundaries, dimension);

		_cell_count = 0;
		for (int i = 0; i < PARALLEL_THREADS; i++) {
			coboundary_matrix_offsets[i] = _cell_count;
			_cell_count += coboundary_matrix[i].size();
		}
		cell_count[dimension] = _cell_count;

#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "preparing dimension " << dimension
		          << ": done computing coboundaries, now reducing memory consumption" << std::flush << "\r";
#endif
	}

	if (cache != "") {
		for (int i = 0; i < PARALLEL_THREADS; i++) {
			std::string fname = cache;
			fname += "/matrix_";
			fname += std::to_string(current_dimension);
			fname += "_";
			fname += std::to_string(i);
			std::ofstream o(fname.c_str(), std::ios::binary);
			index_t separator = -1;
			o.write((char*)&(cell_count[current_dimension + 1]), sizeof(size_t));
			o.write((char*)&(coboundary_matrix_offsets[i]), sizeof(size_t));
			size_t this_size = coboundary_matrix[i].size();
			o.write((char*)&this_size, sizeof(size_t));
			for (auto j = 0ul; j < this_size; j++) {
				o.write((char*)&(separator), sizeof(index_t));
				for (auto it = coboundary_matrix[i].cbegin(j); it != coboundary_matrix[i].cend(j); ++it) {
					o.write((char*)&(*it), sizeof(index_t));
				}
			}
			o.close();
		}

		std::string fname = cache;
		fname += "/filtration_";
		fname += std::to_string(current_dimension);
		std::ofstream o(fname.c_str(), std::ios::binary);
		for (auto f : next_filtration) o.write((char*)&f, sizeof(value_t));
		o.close();
	}
}
