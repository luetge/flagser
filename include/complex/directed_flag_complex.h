#pragma once

#include <cassert>
#include <thread>
#include <unordered_set>
#include <vector>

#include "../definitions.h"
#include "../directed_graph.h"

class directed_flag_complex_coboundary_cell_t;
class directed_flag_complex_boundary_cell_t;

// Note: To save space we do not store the dimension of the cell.
//       If you want to have that dimension along the cells, uncomment the definition of USE_CELLS_WITH_DIMENSION above
class directed_flag_complex_cell_t {
protected:
	const vertex_index_t* vertices;
#ifndef USE_CELLS_WITHOUT_DIMENSION
	unsigned short dim;
#endif

public:
#ifdef USE_CELLS_WITHOUT_DIMENSION
	directed_flag_complex_cell_t() : vertices(nullptr) {}
	directed_flag_complex_cell_t(const vertex_index_t* vertices) : vertices(vertices) {}
#else
	directed_flag_complex_cell_t(unsigned short _dim) : vertices(nullptr), dim(_dim) {}
	directed_flag_complex_cell_t(unsigned short _dim, const vertex_index_t* vertices) : vertices(vertices), dim(_dim) {}
#endif

	vertex_index_t operator[](size_t index) const { return vertex(index); }
	const vertex_index_t* get_vertex_array() const { return vertices; }
	virtual vertex_index_t vertex(size_t index) const { return vertices[index]; }
	directed_flag_complex_coboundary_cell_t insert_vertex(size_t position, vertex_index_t vertex);
	directed_flag_complex_boundary_cell_t boundary(size_t i);
	void set_vertices(vertex_index_t* _vertices) { vertices = _vertices; }
#ifndef USE_CELLS_WITHOUT_DIMENSION
	virtual unsigned short dimension() const { return dim; }
#endif

	const std::string to_string(
#ifdef USE_CELLS_WITHOUT_DIMENSION
	    unsigned short dim
#endif
	) const {
		const directed_flag_complex_cell_t* t = this;
#ifndef USE_CELLS_WITHOUT_DIMENSION
		unsigned short dim = dimension();
#endif
		std::string s;
		for (unsigned short index = 0; index < dim + 1; index++) {
			s += std::to_string(t->vertex(index));
			if (index < dim) s += "|";
		}
		return s;
	}
};

class directed_flag_complex_coboundary_cell_t : public directed_flag_complex_cell_t {
	size_t _insert_position = -1;
	vertex_index_t _insert_vertex = 0;

public:
#ifdef USE_CELLS_WITHOUT_DIMENSION
	directed_flag_complex_coboundary_cell_t(const vertex_index_t* vertices, vertex_index_t insert_vertex,
	                                        unsigned short insert_position)
	    : directed_flag_complex_cell_t(vertices), _insert_position(insert_position), _insert_vertex(insert_vertex) {}
#else
	directed_flag_complex_coboundary_cell_t(unsigned short dimension, const vertex_index_t* vertices,
	                                        vertex_index_t insert_vertex, unsigned short insert_position)
	    : directed_flag_complex_cell_t(dimension, vertices), _insert_position(insert_position),
	      _insert_vertex(insert_vertex) {}
#endif

	virtual vertex_index_t vertex(size_t index) const override {
		if (index == _insert_position) return _insert_vertex;
		return vertices[index - (index < _insert_position ? 0 : 1)];
	}

#ifndef USE_CELLS_WITHOUT_DIMENSION
	virtual unsigned short dimension() const override { return dim + 1; }
#endif
};

directed_flag_complex_coboundary_cell_t directed_flag_complex_cell_t::insert_vertex(size_t position,
                                                                                    vertex_index_t vertex) {
#ifdef USE_CELLS_WITHOUT_DIMENSION
	return directed_flag_complex_coboundary_cell_t(vertices, vertex, (unsigned short)position);
#else
	return directed_flag_complex_coboundary_cell_t(dim, vertices, vertex, (unsigned short)position);
#endif
}

class directed_flag_complex_boundary_cell_t : public directed_flag_complex_cell_t {
	size_t i;

public:
#ifdef USE_CELLS_WITHOUT_DIMENSION
	directed_flag_complex_boundary_cell_t(const vertex_index_t* vertices, size_t _i)
	    : directed_flag_complex_cell_t(vertices), i(_i) {}
#else
	directed_flag_complex_boundary_cell_t(unsigned short dimension, const vertex_index_t* vertices, size_t _i)
	    : directed_flag_complex_cell_t(dimension, vertices), i(_i) {}
#endif

	virtual vertex_index_t vertex(size_t index) const override { return vertices[index + (index < i ? 0 : 1)]; }

#ifndef USE_CELLS_WITHOUT_DIMENSION
	virtual unsigned short dimension() const override { return dim - 1; }
#endif
};

directed_flag_complex_boundary_cell_t directed_flag_complex_cell_t::boundary(size_t i) {
#ifdef USE_CELLS_WITHOUT_DIMENSION
	return directed_flag_complex_boundary_cell_t(vertices, i);
#else
	return directed_flag_complex_boundary_cell_t(dim, vertices, i);
#endif
}

// This hasher hashes arrays of unsigend ints with the boost algorithm
// Note: If the cells do not contain the dimension you have to call cell_hasher_t::set_current_cell_dimension
//       before using both cell_hasher_t and cell_comparer_t.
struct cell_hasher_t {
	size_t operator()(const directed_flag_complex_cell_t& c) const {
		size_t hash = 0;
		const directed_flag_complex_cell_t* cell = &c;
#ifndef USE_CELLS_WITHOUT_DIMENSION
		unsigned short dimension = c.dimension();
#endif

		for (unsigned short index = 0; index < dimension + 1; index++) {
			hash ^= hasher(cell->vertex(index)) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		}

		return hash;
	}

#ifdef USE_CELLS_WITHOUT_DIMENSION
	static void set_current_cell_dimension(unsigned short _dimension) { dimension = _dimension; }
	static unsigned short dimension;
#endif
private:
	std::hash<vertex_index_t> hasher;
};
#ifdef USE_CELLS_WITHOUT_DIMENSION
unsigned short cell_hasher_t::dimension = 0;
#endif

struct cell_comparer_t {
	bool operator()(const directed_flag_complex_cell_t& lhs, const directed_flag_complex_cell_t& rhs) const {
		// Always at most one of the two is a coboundary cell
		const directed_flag_complex_cell_t* l = &lhs;
		const directed_flag_complex_cell_t* r = &rhs;
#ifdef USE_CELLS_WITHOUT_DIMENSION
		unsigned short dimension = cell_hasher_t::dimension;
#else
		assert(l->dimension() == r->dimension());
		unsigned short dimension = l->dimension();
#endif
		for (unsigned short index = 0; index < dimension + 1; index++) {
			if (l->vertex(index) != r->vertex(index)) return false;
		}

		return true;
	}
};

class directed_flag_complex_t {
public:
	const filtered_directed_graph_t& graph;
	directed_flag_complex_t(const filtered_directed_graph_t& _graph) : graph(_graph) {}

public:
	template <typename Func> void for_each_cell(Func& f, int min_dimension, int max_dimension = -1) {
		std::vector<Func> fs{f};
		for_each_cell(fs, min_dimension, max_dimension);
		// Hack because I cannot wrap the reference in the vector, well technically I could use `std::reference_wrapper`
		// But in this case it doesn't help
		// What I do is storing the computed value again after computation
		f = fs[0];
	}

	template <typename Func> void for_each_cell(std::vector<Func>& fs, int min_dimension, int max_dimension = -1) {
		if (max_dimension == -1) max_dimension = min_dimension;
		auto number_of_threads = int(fs.size());
		std::vector<std::thread> t(number_of_threads);

		for (int index = 0; index < number_of_threads - 1; ++index)
			t[index] = std::thread(&directed_flag_complex_t::worker_thread<Func>, this, number_of_threads, index,
			                       std::ref(fs[index]), min_dimension, max_dimension);

		// Also do work in this thread, namely the last bit
		worker_thread(number_of_threads, number_of_threads - 1, fs[number_of_threads - 1], min_dimension,
		              max_dimension);

		// Wait until all threads stopped
		for (int i = 0; i < number_of_threads - 1; ++i) t[i].join();
	}

private:
	template <typename Func>
	void worker_thread(size_t number_of_threads, int thread_id, Func& f, int min_dimension, int max_dimension) {
		const size_t number_of_vertices = graph.vertex_number();

		std::vector<vertex_index_t> first_position_vertices;
		for (size_t index = thread_id; index < number_of_vertices; index += number_of_threads)
			first_position_vertices.push_back(vertex_index_t(index));

		std::vector<vertex_index_t> prefix(max_dimension + 1);

		do_for_each_cell(f, min_dimension, max_dimension, first_position_vertices, prefix.data(), 0);

		f.done();
	}

	template <typename Func>
	void do_for_each_cell(Func& f, int min_dimension, int max_dimension,
	                      const std::vector<vertex_index_t>& possible_next_vertices, vertex_index_t* prefix,
	                      unsigned short prefix_size = 0) {
		// As soon as we have the correct dimension, execute f
		if (prefix_size >= min_dimension + 1) { f(prefix, prefix_size); }

		// If this is the last dimension we are interested in, exit this branch
		if (prefix_size == max_dimension + 1) return;

		for (auto vertex : possible_next_vertices) {
			// We can write the cell given by taking the current vertex as the maximal element
			prefix[prefix_size] = vertex;

			// And compute the next elements
			std::vector<vertex_index_t> new_possible_vertices;
			if (prefix_size > 0) {
				for (auto v : possible_next_vertices) {
					if (vertex != v && graph.is_connected_by_an_edge(vertex, v)) new_possible_vertices.push_back(v);
				}
			} else {
				// Get outgoing vertices of v in chunks of 64
				for (size_t offset = 0; offset < graph.incidence_row_length; offset++) {
					size_t bits = graph.get_outgoing_chunk(vertex, offset);

					size_t vertex_offset = offset << 6;
					while (bits > 0) {
						// Get the least significant non-zero bit
						auto b = __builtin_ctzl(bits);

						// Unset this bit
						bits &= ~(ONE_ << b);

						new_possible_vertices.push_back(vertex_index_t(vertex_offset + b));
					}
				}
			}

			do_for_each_cell(f, min_dimension, max_dimension, new_possible_vertices, prefix, prefix_size + 1);
		}
	}
};
