// WARNING: This code is generated, all changes will be lost after the next compilation.
#ifndef FLAGSER_ALGORITHMS_H
#define FLAGSER_ALGORITHMS_H
#include <string>
#include <vector>

#include "filtration_algorithms.h"

struct dimension_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&, const value_t*) const {
		return (dimension);
	}
	virtual inline bool needs_face_filtration() const { return false; }
	virtual inline bool overwrite_vertex_filtration() const { return true; }
	virtual inline bool overwrite_edge_filtration() const { return true; }
};

struct zero_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&, const value_t*) const {
		return (0.0f);
	}
	virtual inline bool needs_face_filtration() const { return false; }
	virtual inline bool overwrite_vertex_filtration() const { return true; }
	virtual inline bool overwrite_edge_filtration() const { return true; }
};

struct max_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {
		return (max(boundary_filtration, 0, dimension));
	}
};

struct max3_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {
		return (max(boundary_filtration, 0, dimension));
	}
};

struct max_plus_one_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {
		return (max(boundary_filtration, 0, dimension) + 1.0f);
	}
};

struct product_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {
		return (product(boundary_filtration, 0, dimension));
	}
};

struct sum_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {
		return (sum(boundary_filtration, 0, dimension));
	}
};

struct pmean_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {
		return (max(boundary_filtration, 0, dimension) +
		        pow(
		            [&]() {
			            value_t EXPRESSION_0 = 0.0f;
			            for (int COUNTER_A = (int)0; COUNTER_A <= (int)dimension; COUNTER_A++) {
				            EXPRESSION_0 += (pow(boundary_filtration[COUNTER_A], 2.0f));
			            }
			            return EXPRESSION_0;
		            }(),
		            (1.0f / 2.0f)) /
		            (dimension + 1.0f));
	}
};

struct pmoment_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {
		return (max(boundary_filtration, 0, dimension) +
		        pow(
		            [&]() {
			            value_t EXPRESSION_0 = 0.0f;
			            for (int COUNTER_A = (int)0; COUNTER_A <= (int)dimension; COUNTER_A++) {
				            EXPRESSION_0 += (pow((boundary_filtration[COUNTER_A] -
				                                  sum(boundary_filtration, 0, dimension) / (dimension + 1.0f)),
				                                 2.0f));
			            }
			            return EXPRESSION_0;
		            }(),
		            (1.0f / 2.0f)) /
		            (dimension + 1.0f));
	}
};

struct remove_edges_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t&,
	                                          const filtered_directed_graph_t&,
	                                          const value_t* boundary_filtration) const {

		if (dimension == 0) { return (0.0f); }

		if (dimension == 1) { throw std::runtime_error("Please provide a graph with random weights on the edges."); }
		return (max(boundary_filtration, 0, dimension));
	}
};

struct vertex_degree_filtration : public filtration_algorithm_t {
	virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t& cell,
	                                          const filtered_directed_graph_t& graph,
	                                          const value_t* boundary_filtration) const {

		if (dimension == 0) {
			return (-1.0f *
			        (graph.outdegrees[(int)((cell[(int)((0.0f))]))] + graph.indegrees[(int)((cell[(int)((0.0f))]))]));
		}
		return (max(boundary_filtration, 0, dimension));
	}
	virtual inline bool overwrite_vertex_filtration() const { return true; }
	virtual inline bool overwrite_edge_filtration() const { return true; }
};

filtration_algorithm_t* get_custom_filtration_computer(std::string algorithm) {
	if (algorithm == "dimension") return new dimension_filtration();
	if (algorithm == "zero") return new zero_filtration();
	if (algorithm == "max") return new max_filtration();
	if (algorithm == "max3") return new max3_filtration();
	if (algorithm == "max_plus_one") return new max_plus_one_filtration();
	if (algorithm == "product") return new product_filtration();
	if (algorithm == "sum") return new sum_filtration();
	if (algorithm == "pmean") return new pmean_filtration();
	if (algorithm == "pmoment") return new pmoment_filtration();
	if (algorithm == "remove_edges") return new remove_edges_filtration();
	if (algorithm == "vertex_degree") return new vertex_degree_filtration();

	return nullptr;
}
std::vector<std::string> custom_filtration_computer{"dimension",    "zero",         "max",          "max3",
                                                    "max_plus_one", "product",      "sum",          "pmean",
                                                    "pmoment",      "remove_edges", "vertex_degree"};
#endif // FLAGSER_ALGORITHMS_H
