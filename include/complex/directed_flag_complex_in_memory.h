#pragma once

#include <array>
#include <cassert>
#include <map>
#include <thread>
#include <unordered_set>
#include <vector>

#include "../definitions.h"
#include "../directed_graph.h"
#include "../parameters.h"

//
// Loading the whole complex into memory, trading memory for computation speed
//
template <typename ExtraData> class directed_flag_complex_cell_in_memory_t {
public:
	vertex_index_t vertex;
	ExtraData data;
	std::map<vertex_index_t, directed_flag_complex_cell_in_memory_t<ExtraData>>* children = nullptr;

	directed_flag_complex_cell_in_memory_t(vertex_index_t _vertex) : vertex(_vertex) {}

	// TODO: Figure out how to do this in the destructor.
	//       The problem was that the classes are copied around, so the
	//       destructor was called before the actual destruction of the complex
	//       took place.
	void free_memory() {
		if (children != nullptr) {
			for (auto c : *children) c.second.free_memory();
			delete children;
		}
	}

	directed_flag_complex_cell_in_memory_t& add_child(vertex_index_t vertex) {
		if (children == nullptr)
			children = new std::map<vertex_index_t, directed_flag_complex_cell_in_memory_t<ExtraData>>();
		directed_flag_complex_cell_in_memory_t<ExtraData> new_cell(vertex);

		// Insert and return the cell
		return children->insert(std::make_pair(new_cell.vertex, new_cell)).first->second;
	}

	template <typename Cell> void set_data(int dimension, const Cell vertices, ExtraData _data, int offset = 0) {
		if (dimension == 0) {
			data = _data;
			return;
		}

		if (children == nullptr) {
			throw std::runtime_error("A cell could not be found in the directed flag complex.");
		}

		offset++;
		auto pair = children->find(vertices[offset]);
		if (pair == children->end()) {
			throw std::runtime_error("A cell could not be found in the directed flag complex.");
		}

		pair->second.set_data(dimension - 1, vertices, _data, offset);
	}

	template <typename Cell> const ExtraData& get_data(int dimension, const Cell vertices, int offset = 0) const {
		if (dimension == 0) return data;

		if (children == nullptr) {
			throw std::runtime_error("A cell could not be found in the directed flag complex.");
		}

		offset++;
		auto pair = children->find(vertices[offset]);
		if (pair == children->end()) {
			throw std::runtime_error("A cell could not be found in the directed flag complex.");
		}

		return pair->second.get_data(dimension - 1, vertices, offset);
	}

	template <typename Func>
	void for_each_cell(Func& f, int min_dimension, int max_dimension, vertex_index_t* prefix, int prefix_size = 0) {
		prefix[prefix_size++] = vertex;
		if (prefix_size >= min_dimension + 1) f(prefix, prefix_size, data);
		if (prefix_size == max_dimension + 1 || children == nullptr) return;

		for (auto child : *children) child.second.for_each_cell(f, min_dimension, max_dimension, prefix, prefix_size);
	}
};

template <typename ExtraData> struct euler_characteristic_computer_t {
	void done() {}
	void operator()(vertex_index_t*, int size, ExtraData) {
		// Add (-1)^size to the Euler characteristic
		if (size & 1)
			ec++;
		else
			ec--;
	}

	index_t euler_characteristic() const { return ec; }

private:
	index_t ec = 0;
};

template <typename ExtraData> class directed_flag_complex_in_memory_t {
public:
	std::vector<directed_flag_complex_cell_in_memory_t<ExtraData>> vertex_cells;
	const size_t nb_threads;

	directed_flag_complex_in_memory_t(const directed_graph_t& graph, const size_t nb_threads, int max_dimension = -1);
	~directed_flag_complex_in_memory_t() {
		for (auto p : vertex_cells) p.free_memory();
	}

	template <typename Cell> void set_data(int dimension, const Cell vertices, ExtraData _data) {
		vertex_cells[*vertices].set_data(dimension, vertices, _data);
	}

	template <typename Cell> const ExtraData& get_data(int dimension, const Cell vertices) const {
		return vertex_cells[vertices[0]].get_data(dimension, vertices);
	}

	template <typename Func> void for_each_cell(Func& f, int min_dimension, int max_dimension = -1) {
		std::vector<Func> fs{{f}};
		for_each_cell(fs, min_dimension, max_dimension);
	}

	index_t euler_characteristic() {
		std::vector<euler_characteristic_computer_t<ExtraData>> compute_euler_characteristic;
		for (auto i = 0ul; i < nb_threads; i++)
			compute_euler_characteristic.push_back(euler_characteristic_computer_t<ExtraData>());

		// TODO: Remove the 10000-hack
		for_each_cell(compute_euler_characteristic, 0, 10000);
		index_t euler_characteristic = 0;
		for (auto i = 0ul; i < nb_threads; i++)
			euler_characteristic += compute_euler_characteristic[i]->euler_characteristic();

		return euler_characteristic;
	}

	template <typename Func> void for_each_cell(std::vector<Func>& fs, int min_dimension, int max_dimension = -1) {
		if (max_dimension == -1) max_dimension = min_dimension;

		auto number_of_threads = int(fs.size());
		std::vector<std::thread> t(number_of_threads);

		for (auto index = 0; index < number_of_threads - 1; ++index)
			t[index] = std::thread(&directed_flag_complex_in_memory_t<ExtraData>::worker_thread<Func>, this,
			                       number_of_threads, index, std::ref(fs[index]), min_dimension, max_dimension);

		// Also do work in this thread, namely the last bit
		worker_thread(number_of_threads, number_of_threads - 1, fs[number_of_threads - 1], min_dimension,
		              max_dimension);

		// Wait until all threads stopped
		for (auto i = 0; i < number_of_threads - 1; ++i) t[i].join();
	}

	template <typename Func>
	void worker_thread(int number_of_threads, int thread_id, Func& f, int min_dimension, int max_dimension) {
		const size_t number_of_vertices = vertex_cells.size();

		std::vector<vertex_index_t> prefix(max_dimension + 1);

		for (vertex_index_t index = thread_id; index < number_of_vertices; index += number_of_threads) {
			vertex_cells[index].for_each_cell(f, min_dimension, max_dimension, prefix.data());
		}

		f.done();
	}
};

template <typename ExtraData>
void construct_children(directed_flag_complex_cell_in_memory_t<ExtraData>& current_cell, const directed_graph_t& graph,
                        const std::vector<vertex_index_t>& possible_next_vertices, int max_dimension,
                        int current_dimension = 0) {
	if (current_dimension == max_dimension) return;

	for (auto vertex : possible_next_vertices) {
		if (vertex == current_cell.vertex) continue;

		// Compute the next elements
		std::vector<vertex_index_t> new_possible_vertices;
		for (auto v : possible_next_vertices)
			if (v != vertex && v != current_cell.vertex && graph.is_connected_by_an_edge(vertex, v))
				new_possible_vertices.push_back(v);

		directed_flag_complex_cell_in_memory_t<ExtraData>& new_cell = current_cell.add_child(vertex);
		construct_children(new_cell, graph, new_possible_vertices, max_dimension, current_dimension + 1);
	}
}

template <typename ExtraData>
void construction_worker_thread(int number_of_threads, int thread_id,
                                directed_flag_complex_in_memory_t<ExtraData>* result, const directed_graph_t* graph,
                                int max_dimension) {
	const size_t number_of_vertices = graph->vertex_number();

	for (size_t index = thread_id; index < number_of_vertices; index += number_of_threads) {
		// Compute possible vertices
		std::vector<vertex_index_t> possible_vertices;
		for (size_t offset = 0; offset < graph->incidence_row_length; offset++) {
			size_t bits = graph->get_outgoing_chunk(vertex_index_t(index), vertex_index_t(offset));
			size_t vertex_offset = offset << 6;

			while (bits > 0) {
				// Get the least significant non-zero bit
				auto b = __builtin_ctzl(bits);

				// Unset this bit
				bits &= ~(ONE_ << b);
				possible_vertices.push_back(vertex_index_t(vertex_offset + b));
			}
		}
		construct_children(result->vertex_cells[index], *graph, possible_vertices, max_dimension);
	}
}

template <typename ExtraData>
directed_flag_complex_in_memory_t<ExtraData>::directed_flag_complex_in_memory_t(const directed_graph_t& graph,
                                                                                const size_t nb_threads,
                                                                                int max_dimension)
    : nb_threads(nb_threads) {
#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	          << "constructing the directed flag complex" << std::flush << "\r";
#endif

	// First we add the vertices of the directed flag complex
	vertex_cells.reserve(graph.number_of_vertices);
	for (vertex_index_t index = 0; index < graph.number_of_vertices; index++)
		vertex_cells.push_back(directed_flag_complex_cell_in_memory_t<ExtraData>(index));

	// Now we start a few threads to construct the flag complex
	std::vector<std::thread> t(nb_threads - 1);

	for (auto index = 0ul; index < nb_threads - 1; ++index)
		t[index] = std::thread(&construction_worker_thread<ExtraData>, nb_threads, index, this, &graph, max_dimension);

	// Also do work in this thread, namely the last bit
	// For this last thread, take all the remaining vertices
	construction_worker_thread(nb_threads, nb_threads - 1, this, &graph, max_dimension);

	// Wait until all threads stopped
	for (auto i = 0ul; i < nb_threads - 1; ++i) t[i].join();
}
