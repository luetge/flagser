#pragma once

#define USE_CELLS_WITHOUT_DIMENSION

#include <algorithm>
#include <array>

#include "../argparser.h"
#include "../directed_graph.h"
#include "../filtration_algorithms.h"
#include "../persistence.h"

#include "directed_flag_complex_in_memory.h"

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
	reorder_edges_t(std::vector<vertex_index_t>& _new_edges, std::vector<value_t>& _new_filtration,
	                bool _reorder_filtration)
	    : new_edges(_new_edges), new_filtration(_new_filtration), reorder_filtration(_reorder_filtration) {}
	void done() {}
	void operator()(vertex_index_t* first_vertex, int, std::pair<index_t, value_t>& data) {
		new_edges.push_back(*first_vertex++);
		new_edges.push_back(*first_vertex);
		if (reorder_filtration) new_filtration.push_back(data.second);
	}

private:
	std::vector<vertex_index_t>& new_edges;
	std::vector<value_t>& new_filtration;
	bool reorder_filtration;
};

struct compute_filtration_t {
	compute_filtration_t(filtration_algorithm_t* _filtration_algorithm, filtered_directed_graph_t& _graph,
	                     directed_flag_complex_in_memory_t<std::pair<index_t, value_t>>& _complex)
	    : filtration_algorithm(_filtration_algorithm), graph(_graph), complex(_complex) {}
	void done() {}
	void operator()(vertex_index_t* first_vertex, int size, std::pair<index_t, value_t>&) {
		// The index starts at -1, so the first cell gets index 0
		current_index++;

		value_t nf = 0.0f;

		if (filtration_algorithm != nullptr) {
			directed_flag_complex_cell_t cell(first_vertex);
			if (cell_hasher_t::dimension != size - 2) cell_hasher_t::set_current_cell_dimension(size - 2);

			if (filtration_algorithm->needs_face_filtration()) {
				if (boundary_filtration == nullptr) boundary_filtration = new value_t[size];

				for (int i = 0; i < size; i++) {
					boundary_filtration[i] = complex.get_data(size - 2, cell.boundary(i)).second;
				}

				nf = filtration_algorithm->compute_filtration(size - 1, cell, graph, boundary_filtration);
			} else {
				nf = filtration_algorithm->compute_filtration(size - 1, cell, graph, nullptr);
			}
			next_filtration.push_back(nf);
		}
		complex.set_data(size - 1, first_vertex, std::make_pair(current_index, nf));
	}

	std::vector<value_t> filtration() const { return next_filtration; }

	index_t number_of_cells() const { return current_index + 1; }

private:
	index_t current_index = -1;
	filtration_algorithm_t* filtration_algorithm;
	filtered_directed_graph_t& graph;
	directed_flag_complex_in_memory_t<std::pair<index_t, value_t>>& complex;

	value_t* boundary_filtration = nullptr;
	std::vector<value_t> next_filtration;
};

template <typename Complex>
void prepare_graph_filtration(Complex& complex, filtered_directed_graph_t& graph,
                              filtration_algorithm_t* filtration_algorithm, const size_t nb_threads) {
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

		compute_filtration_t compute_filtration(filtration_algorithm, graph, complex);
		complex.for_each_cell(compute_filtration, 1);
		graph.edge_filtration.clear();
		graph.edge_filtration = compute_filtration.filtration();
	}

	// Reorder the edge filtration
	if (!computed_edge_filtration && filtration_algorithm != nullptr) {
		vertex_index_t* edge = &graph.edges[0];
		value_t* filtration = &graph.edge_filtration[0];
		for (size_t index = 0; index < graph.edge_number(); ++index, edge += 2, ++filtration) {
			complex.set_data(1, edge, std::make_pair(-1, *filtration));
		}
	}

	// Now reorder the edges
	std::vector<std::vector<vertex_index_t>> new_edges(nb_threads);
	std::vector<std::vector<value_t>> new_filtrations(nb_threads);
	std::vector<reorder_edges_t> reorder_filtration;
	for (size_t i = 0; i < nb_threads; i++)
		reorder_filtration.push_back(reorder_edges_t(new_edges[i], new_filtrations[i],
		                                             !computed_edge_filtration && filtration_algorithm != nullptr));
	complex.for_each_cell(reorder_filtration, 1);

	size_t current_filtration_index = 0;
	size_t current_edge_index = 0;
	for (size_t i = 0; i < nb_threads; i++) {
		for (auto e : new_edges[i]) graph.edges[current_edge_index++] = e;
		if (!computed_edge_filtration && filtration_algorithm != nullptr)
			for (auto f : new_filtrations[i]) graph.edge_filtration[current_filtration_index++] = f;
	}
}

class directed_flag_complex_in_memory_computer_t {
	filtered_directed_graph_t& graph;
	directed_flag_complex_in_memory_t<std::pair<index_t, value_t>> flag_complex;
	filtration_algorithm_t* filtration_algorithm;
	unsigned short min_dimension;
	unsigned short max_dimension;
	const char* cache;
	int current_dimension = 0;
	bool _is_top_dimension = false;
	std::vector<size_t> cell_count;
	const size_t nb_threads;

	// Filtration
	std::vector<value_t> next_filtration;

	// Coboundaries
	std::vector<compressed_sparse_matrix<entry_t>> coboundary_matrix;
	std::vector<size_t> coboundary_matrix_offsets;
	coefficient_t modulus;

public:
	directed_flag_complex_in_memory_computer_t(filtered_directed_graph_t& _graph, const flagser_parameters& params)
	    : graph(_graph), flag_complex(graph, params.nb_threads),
	      filtration_algorithm(params.filtration_algorithm.get()), min_dimension(params.min_dimension),
	      max_dimension(params.max_dimension), cache(params.cache.c_str()), nb_threads(params.nb_threads),
	      coboundary_matrix(nb_threads), coboundary_matrix_offsets(nb_threads, 0), modulus(params.modulus) {
		cell_count.push_back(_graph.vertex_number());
		cell_count.push_back(_graph.edge_number());

		// Order the edges and compute the correct filtration
		if (min_dimension <= 1 || (filtration_algorithm != nullptr && filtration_algorithm->needs_face_filtration()))
			prepare_graph_filtration(flag_complex, graph, filtration_algorithm, nb_threads);
	}

	size_t number_of_cells(int dimension) const {
		assert(size_t(dimension) < cell_count.size());
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

	inline coboundary_iterator_t<directed_flag_complex_in_memory_computer_t> coboundary(filtration_entry_t cell) {
		if (current_dimension > max_dimension || is_top_dimension()) {
			return coboundary_iterator_t<directed_flag_complex_in_memory_computer_t>(
			    this, current_dimension, coboundary_matrix[0], -1, 1, modulus);
		}

		size_t i = 0;
		while (i < nb_threads - 1 && index_t(coboundary_matrix_offsets[i + 1]) <= get_index(cell)) { i++; }
		return coboundary_iterator_t<directed_flag_complex_in_memory_computer_t>(
		    this, current_dimension, coboundary_matrix[i], index_t(get_index(cell) - coboundary_matrix_offsets[i]),
		    get_coefficient(cell), modulus);
	}

	bool is_top_dimension() { return _is_top_dimension; }

	void computation_result(int, size_t, size_t) {}
	void finished() {}
};

struct store_coboundaries_in_cache_t {
	store_coboundaries_in_cache_t(compressed_sparse_matrix<entry_t>& _coboundary_matrix, int _current_dimension,
	                              const filtered_directed_graph_t& _graph,
	                              const directed_flag_complex_in_memory_t<std::pair<index_t, value_t>>& _complex,
	                              std::vector<size_t>& _cell_index_offsets, size_t _total_cell_number, bool _is_first,
	                              const size_t nb_threads, coefficient_t _modulus = 2)
	    : is_first(_is_first), current_dimension(_current_dimension), coboundary_matrix(_coboundary_matrix),
	      graph(_graph), complex(_complex), cell_index_offsets(_cell_index_offsets),
	      total_cell_number(_total_cell_number), nb_threads(nb_threads), modulus(_modulus) {}

	void done() {
#ifdef INDICATE_PROGRESS
		if (is_first)
			std::cout << "\033[K"
			          << "dimension " << current_dimension << ": computed almost all of the coboundaries" << std::flush
			          << "\r";
#endif
	}
	void operator()(vertex_index_t* first_vertex, int size, std::pair<index_t, value_t>&) {
#ifdef INDICATE_PROGRESS
		if (is_first && (current_index + 1) % 10000 == 0) {
			std::cout << "\033[K"
			          << "dimension " << current_dimension << ": computed ca. " << nb_threads * (current_index + 1);
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
					auto b = __builtin_ctzl(bits);

					// Unset this bit
					bits &= ~(ONE_ << b);

					// Now insert the appropriate vertex at this position
					auto cb = cell.insert_vertex(i, vertex_index_t(vertex_offset + b));
					short thread_index = cb.vertex(0) % nb_threads;
					coboundary_matrix.push_back(make_entry(complex.get_data(current_dimension + 1, cb).first +
					                                           index_t(cell_index_offsets[thread_index]),
					                                       index_t(i & 1 ? -1 + modulus : 1)));
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
	const directed_flag_complex_in_memory_t<std::pair<index_t, value_t>>& complex;
	std::vector<size_t>& cell_index_offsets;
	size_t total_cell_number;
	const size_t nb_threads;
	coefficient_t modulus;
};

void directed_flag_complex_in_memory_computer_t::prepare_next_dimension(int dimension) {
	if (dimension == 0) return;

	assert(dimension == current_dimension + 1);
	current_dimension = dimension;
	cell_count.resize(dimension + 2);

	for (size_t i = 0; i < nb_threads; i++) { coboundary_matrix[i].clear(); }

	if (dimension > max_dimension || _is_top_dimension) return;

	{
		std::vector<size_t> _next_cells_offsets(nb_threads, 0);
		{
			// If we will actually compute coboundaries, then compute the filtration.
			// Also if we need the face filtrations.
			if (dimension + 1 >= min_dimension ||
			    (filtration_algorithm != nullptr && filtration_algorithm->needs_face_filtration())) {
#ifdef INDICATE_PROGRESS
				std::cout << "\033[K"
				          << "preparing dimension " << dimension << ": computing the filtration of all "
				          << (dimension + 1) << "-dimensional cells" << std::flush << "\r";
#endif
				std::vector<compute_filtration_t> compute_filtration;
				for (size_t i = 0; i < nb_threads; i++) {
					compute_filtration.push_back(compute_filtration_t(filtration_algorithm, graph, flag_complex));
				}
				flag_complex.for_each_cell(compute_filtration, dimension + 1);

				size_t _cell_count = 0;
				for (size_t i = 0; i < nb_threads; i++) {
					_next_cells_offsets[i] = _cell_count;
					_cell_count += compute_filtration[i].number_of_cells();
				}
				cell_count[dimension + 1] = _cell_count;
				if (_cell_count == 0) _is_top_dimension = true;

				if (filtration_algorithm != nullptr) {
					// Combine the filtration
					next_filtration.clear();
					next_filtration.reserve(_cell_count);
					for (size_t i = 0; i < nb_threads; i++) {
						for (auto f : compute_filtration[i].filtration()) next_filtration.push_back(f);
					}
				}
			}
		}

		// If we do not want the homology in the next degree, then we can stop
		// here
		if (dimension + 1 < min_dimension) return;

		if (!is_top_dimension()) {
#ifdef INDICATE_PROGRESS
			std::cout << "\033[K"
			          << "preparing dimension " << dimension << ": computing the coboundaries of all " << dimension
			          << "-dimensional cells" << std::flush << "\r";
#endif

			// Now compute the coboundaries
			std::vector<store_coboundaries_in_cache_t> store_coboundaries;
			for (size_t i = 0; i < nb_threads; i++) coboundary_matrix[i] = compressed_sparse_matrix<entry_t>();
			for (size_t i = 0; i < nb_threads; i++) {
				store_coboundaries.push_back(store_coboundaries_in_cache_t(
				    coboundary_matrix[i], dimension, graph, flag_complex, _next_cells_offsets,
				    int(cell_count[dimension]), i == 0, nb_threads, modulus));
			}
			flag_complex.for_each_cell(store_coboundaries, dimension);

			size_t _cell_count = 0;
			for (size_t i = 0; i < nb_threads; i++) {
				coboundary_matrix_offsets[i] = _cell_count;
				_cell_count += coboundary_matrix[i].size();
			}
			cell_count[dimension] = _cell_count;
		}
	}
}
