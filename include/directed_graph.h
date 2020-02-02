#pragma once

#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "definitions.h"
#include "persistence.h"

class directed_graph_t {
public:
	// The filtration values of the vertices
	vertex_index_t number_of_vertices;
	std::vector<vertex_index_t> edges;
	mutable std::vector<value_t> edge_filtration;
	std::vector<size_t> outdegrees;
	std::vector<size_t> indegrees;
  bool directed = true;

	// These are the incidences as a matrix of 64-bit masks
	std::vector<size_t> incidence_incoming;
	std::vector<size_t> incidence_outgoing;
	size_t incidence_row_length;

	// Assume by default that the edge density will be roughly one percent
	directed_graph_t(vertex_index_t _number_of_vertices, bool directed = true, float density_hint = 0.01)
	    : number_of_vertices(_number_of_vertices), directed(directed), incidence_row_length((_number_of_vertices >> 6) + 1) {
		outdegrees.resize(_number_of_vertices, 0);
		indegrees.resize(_number_of_vertices, 0);
		incidence_incoming.resize(incidence_row_length * _number_of_vertices, 0);
		incidence_outgoing.resize(incidence_row_length * _number_of_vertices, 0);

		size_t vertex_number = (size_t)_number_of_vertices;
		edges.reserve(size_t(vertex_number * density_hint * vertex_number * 2));
	}

	vertex_index_t vertex_number() const { return number_of_vertices; }
	size_t edge_number() const { return edges.size() / 2; }

	bool add_edge(vertex_index_t v, vertex_index_t w) {
    if (!directed && v > w) return add_edge(w, v);

		const size_t vv = v >> 6;
		const size_t ww = w >> 6;

    // Prevent multiple insertions
		if (incidence_outgoing[v * incidence_row_length + ww] & (1UL << ((w - (ww << 6))))) return false;

		outdegrees[v]++;
		indegrees[w]++;
		edges.push_back(v);
		edges.push_back(w);

		incidence_outgoing[v * incidence_row_length + ww] |= 1UL << ((w - (ww << 6)));

		incidence_incoming[w * incidence_row_length + vv] |= 1UL << (v - (vv << 6));
    return true;
	}

	bool is_connected_by_an_edge(vertex_index_t from, vertex_index_t to) const {
		const auto t = to >> 6;
		return incidence_outgoing[incidence_row_length * from + t] & (1UL << (to - (t << 6)));
	}

	size_t get_outgoing_chunk(vertex_index_t from, size_t chunk_number) const {
		return incidence_outgoing[incidence_row_length * from + chunk_number];
	}

	size_t get_incoming_chunk(vertex_index_t from, size_t chunk_number) const {
		return incidence_incoming[incidence_row_length * from + chunk_number];
	}
};

class filtered_directed_graph_t : public directed_graph_t {
public:
	// The filtration values of the vertices
	std::vector<value_t> vertex_filtration;
	std::vector<value_t> edge_filtration;

	filtered_directed_graph_t(const std::vector<value_t> _vertex_filtration, bool directed)
	    : directed_graph_t(_vertex_filtration.size(), directed), vertex_filtration(_vertex_filtration) {}

	// WARNING: This does not take the filtration into account!
	// TODO: Think about how to do this efficiently.
	filtered_directed_graph_t(filtered_directed_graph_t* big_graph, std::unordered_set<vertex_index_t> subset)
	    : filtered_directed_graph_t(std::vector<value_t>(subset.size(), 0), big_graph->directed) {
		// Add the edges
		std::unordered_map<vertex_index_t, vertex_index_t> vertex_indices;
		vertex_index_t index = 0;
		for (auto v : subset) vertex_indices.insert(std::make_pair(v, index++));

		for (auto v : subset) {
			// Check intersections in chunks of 64
			for (size_t offset = 0; offset < big_graph->incidence_row_length; offset++) {
				auto bits = big_graph->get_outgoing_chunk(v, offset);
				size_t vertex_offset = offset << 6;

				while (bits > 0) {
					// Get the least significant non-zero bit
					int b = __builtin_ctzl(bits);

					// Unset this bit
					bits &= ~(1UL << b);

					if (subset.find(vertex_offset + b) != subset.end())
						add_filtered_edge(vertex_indices[v], vertex_indices[vertex_offset + b], 0);
				}
			}
		}
	}

	void add_filtered_edge(vertex_index_t v, vertex_index_t w, value_t filtration) {
		if (directed_graph_t::add_edge(v, w) && filtration != std::numeric_limits<value_t>::lowest())
      edge_filtration.push_back(filtration);
	}

	std::vector<filtered_directed_graph_t> get_connected_subgraphs(vertex_index_t minimal_number_of_vertices) {
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "computing connected subgraphs" << std::flush << "\r";
#endif
		std::vector<value_t> weights(number_of_vertices, 0);
		filtered_union_find dset(weights);
		const auto number_of_edges = edge_number();

		for (size_t index = 0; index < number_of_edges; index++) {
			index_t u = dset.find(edges[2 * index]), v = dset.find(edges[2 * index + 1]);

			if (u != v) dset.link(u, v);
		}

		std::vector<filtered_directed_graph_t> subgraphs;
		{
			std::unordered_map<vertex_index_t, std::unordered_set<vertex_index_t>> connected_components;
			for (vertex_index_t index = 0; index < number_of_vertices; ++index) {
				vertex_index_t component = dset.find(index);
				if (connected_components.find(component) == connected_components.end()) {
					connected_components.insert(std::make_pair(component, std::unordered_set<vertex_index_t>()));
				}

				connected_components[component].insert(index);
			}

			for (auto pair : connected_components) {
				auto vertex_set = pair.second;
				if (vertex_set.size() < minimal_number_of_vertices) continue;
				subgraphs.push_back(filtered_directed_graph_t(this, vertex_set));
			}
		}

		return subgraphs;
	}
};