/*

Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

Copyright 2015-2016 Ulrich Bauer.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

//#define ASSEMBLE_REDUCTION_MATRIX
//#define USE_COEFFICIENTS

//#define INDICATE_PROGRESS
#define PRINT_PERSISTENCE_PAIRS

//#define USE_GOOGLE_HASHMAP

#include "../include/argparser.h"
#include "../include/output/output_classes.h"
#include "../include/usage/ripser.h"
#include "../include/persistence.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>

class binomial_coeff_table {
	std::vector<std::vector<index_t>> B;
	index_t n_max, k_max;

public:
	binomial_coeff_table(index_t n, index_t k) {
		n_max = n;
		k_max = k;

		B.resize(n + 1);
		for (index_t i = 0; i <= n; i++) {
			B[i].resize(k + 1);
			for (index_t j = 0; j <= std::min(i, k); j++) {
				if (j == 0 || j == i)
					B[i][j] = 1;
				else
					B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
			}
		}
	}

	index_t operator()(index_t n, index_t k) const {
		assert(n <= n_max);
		assert(k <= k_max);
		return B[n][k];
	}
};

bool is_prime(const coefficient_t n) {
	if (!(n & 1) || n < 2) return n == 2;
	for (coefficient_t p = 3, q = n / p, r = n % p; p <= q; p += 2, q = n / p, r = n % p)
		if (!r) return false;
	return true;
}

index_t get_next_vertex(index_t& v, const index_t idx, const index_t k, const binomial_coeff_table& binomial_coeff) {
	if (binomial_coeff(v, k) > idx) {
		index_t count = v;
		while (count > 0) {
			index_t i = v;
			index_t step = count >> 1;
			i -= step;
			if (binomial_coeff(i, k) > idx) {
				v = --i;
				count -= step + 1;
			} else
				count = step;
		}
	}
	assert(binomial_coeff(v, k) <= idx);
	assert(binomial_coeff(v + 1, k) > idx);
	return v;
}

template <typename OutputIterator>
OutputIterator get_simplex_vertices(index_t idx, const index_t dim, index_t v,
                                    const binomial_coeff_table& binomial_coeff, OutputIterator out) {
	--v;
	for (index_t k = dim + 1; k > 0; --k) {
		get_next_vertex(v, idx, k, binomial_coeff);
		*out++ = v;
		idx -= binomial_coeff(v, k);
	}
	return out;
}

std::vector<index_t> vertices_of_simplex(const index_t simplex_index, const index_t dim, const index_t n,
                                         const binomial_coeff_table& binomial_coeff) {
	std::vector<index_t> vertices;
	get_simplex_vertices(simplex_index, dim, n, binomial_coeff, std::back_inserter(vertices));
	return vertices;
}

template <typename DistanceMatrix> class rips_filtration_comparator {
public:
	const DistanceMatrix& dist;
	const index_t dim;

private:
	mutable std::vector<index_t> vertices;
	const binomial_coeff_table& binomial_coeff;

public:
	rips_filtration_comparator(const DistanceMatrix& _dist, const index_t _dim,
	                           const binomial_coeff_table& _binomial_coeff)
	    : dist(_dist), dim(_dim), vertices(_dim + 1), binomial_coeff(_binomial_coeff){};

	value_t filtration(const index_t index) const {
		value_t filtr = 0;
		get_simplex_vertices(index, dim, index_t(dist.size()), binomial_coeff, vertices.begin());

		for (index_t i = 0; i <= dim; ++i)
			for (index_t j = 0; j < i; ++j) { filtr = std::max(filtr, dist(vertices[i], vertices[j])); }
		return filtr;
	}

	bool operator()(const index_t a, const index_t b) const {
		assert(a < binomial_coeff(dist.size(), dim + 1));
		assert(b < binomial_coeff(dist.size(), dim + 1));

		return greater_filtration_or_smaller_index<filtration_index_t>()(filtration_index_t(filtration(a), a),
		                                                                 filtration_index_t(filtration(b), b));
	}

	template <typename Entry> bool operator()(const Entry& a, const Entry& b) const {
		return operator()(get_index(a), get_index(b));
	}
};

template <class DistanceMatrix> class simplex_coboundary_enumerator {
private:
	const filtration_entry_t simplex;
	index_t idx_below, idx_above, v, k;
	const coefficient_t modulus;
	const binomial_coeff_table& binomial_coeff;
	const DistanceMatrix& dist;
	std::vector<index_t> vertices;

public:
	simplex_coboundary_enumerator(const filtration_entry_t _simplex, index_t _dim, index_t _n,
	                              const coefficient_t _modulus, const DistanceMatrix& _dist,
	                              const binomial_coeff_table& _binomial_coeff)
	    : simplex(_simplex), idx_below(get_index(_simplex)), idx_above(0), v(_n - 1), k(_dim + 1), modulus(_modulus),
	      binomial_coeff(_binomial_coeff), dist(_dist), vertices(_dim + 1) {
		get_simplex_vertices(get_index(_simplex), _dim, _n, binomial_coeff, vertices.begin());
	}

	bool has_next() {
		while ((v != -1) && (binomial_coeff(v, k) <= idx_below)) {
			idx_below -= binomial_coeff(v, k);
			idx_above += binomial_coeff(v, k + 1);

			--v;
			--k;
			assert(k != -1);
		}
		return v != -1;
	}

	index_t next_index() { return idx_above + binomial_coeff(v--, k + 1) + idx_below; }

	filtration_entry_t next() {
		value_t coface_filtration = get_filtration(simplex);
		for (index_t w : vertices) coface_filtration = std::max(coface_filtration, dist(v, w));
		coefficient_t coface_coefficient = (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;
		return filtration_entry_t(coface_filtration, idx_above + binomial_coeff(v--, k + 1) + idx_below,
		                          coface_coefficient);
	}
};

template <compressed_matrix_layout Layout> class compressed_distance_matrix {
public:
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	void init_rows();

	compressed_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(_distances), rows(size_t(1 + std::sqrt(value_t(1 + 8 * distances.size()))) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	template <typename DistanceMatrix>
	compressed_distance_matrix(const DistanceMatrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();

		for (size_t i = 1ul; i < size(); ++i)
			for (size_t j = 0ul; j < i; ++j) rows[i][j] = mat(index_t(i), index_t(j));
	}

	value_t operator()(const index_t i, const index_t j) const;

	size_t size() const { return rows.size(); }
};

template <> void compressed_distance_matrix<LOWER_TRIANGULAR>::init_rows() {
	value_t* pointer = &distances[0];
	for (auto i = 1ul; i < size(); ++i) {
		rows[i] = pointer;
		pointer += i;
	}
}

template <> void compressed_distance_matrix<UPPER_TRIANGULAR>::init_rows() {
	value_t* pointer = &distances[0] - 1;
	for (auto i = 0ul; i < size() - 1; ++i) {
		rows[i] = pointer;
		pointer += size() - i - 2;
	}
}

template <> value_t compressed_distance_matrix<UPPER_TRIANGULAR>::operator()(index_t i, index_t j) const {
	if (i > j) std::swap(i, j);
	return i == j ? 0 : rows[i][j];
}

template <> value_t compressed_distance_matrix<LOWER_TRIANGULAR>::operator()(index_t i, index_t j) const {
	if (i > j) std::swap(i, j);
	return i == j ? 0 : rows[j][i];
}

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;

class euclidean_distance_matrix {
public:
	std::vector<std::vector<value_t>> points;

	euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points) : points(_points) {}

	value_t operator()(const index_t i, const index_t j) const {
		return std::sqrt(std::inner_product(points[i].begin(), points[i].end(), points[j].begin(), value_t(),
		                                    std::plus<value_t>(),
		                                    [](value_t u, value_t v) { return (u - v) * (u - v); }));
	}

	size_t size() const { return points.size(); }
};

class union_find {
	std::vector<index_t> parent;
	std::vector<uint8_t> rank;

public:
	union_find(index_t n) : parent(n), rank(n, 0) {
		for (index_t i = 0; i < n; ++i) parent[i] = i;
	}

	index_t find(index_t x) {
		index_t y = x, z = parent[y];
		while (z != y) {
			y = z;
			z = parent[y];
		}
		y = parent[x];
		while (z != y) {
			parent[x] = z;
			x = y;
			y = parent[x];
		}
		return z;
	}
	void link(index_t x, index_t y) {
		x = find(x);
		y = find(y);
		if (x == y) return;
		if (rank[x] > rank[y])
			parent[y] = x;
		else {
			parent[x] = y;
			if (rank[x] == rank[y]) ++rank[y];
		}
	}
};

template <typename DistanceMatrix> class vietoris_rips_complex_t {
	DistanceMatrix& distance_matrix;
	size_t n;
	unsigned short max_dimension;
	int current_dimension = 0;
	rips_filtration_comparator<DistanceMatrix>* comparator = nullptr;
	rips_filtration_comparator<DistanceMatrix>* next_comparator = nullptr;
	binomial_coeff_table binomial_coeff;
	mutable std::vector<index_t> _vertices_of_edge;
  coefficient_t modulus;

public:
	vietoris_rips_complex_t(DistanceMatrix& _distance_matrix, unsigned short _max_dimension, coefficient_t _modulus)
	    : distance_matrix(_distance_matrix), n(_distance_matrix.size()), max_dimension(_max_dimension),
	      binomial_coeff(index_t(n), _max_dimension + 2), _vertices_of_edge(2, 0), modulus(_modulus) {}

	size_t number_of_cells(int dimension) const { return binomial_coeff(index_t(n), dimension + 1); }

	const std::vector<value_t> vertex_filtration() const { return std::vector<value_t>(n, 0); }

	bool top_dimension() const { return false; }

	inline value_t filtration(int dimension, index_t index) const {
		if (dimension == current_dimension) return comparator->filtration(index);
		if (dimension == current_dimension + 1) return next_comparator->filtration(index);

		std::cerr << "Called filtration(dim, index) with wrong dimension: " << dimension << " (current dimension is "
		          << current_dimension << ")." << std::endl;
		exit(-1);
	}

	std::pair<vertex_index_t, vertex_index_t> vertices_of_edge(index_t edge) const {
		_vertices_of_edge.clear();
		get_simplex_vertices(edge, 1, index_t(n), binomial_coeff, std::back_inserter(_vertices_of_edge));
		return std::make_pair(_vertices_of_edge[0], _vertices_of_edge[1]);
	}

	// Note: Gets called with consecutive dimensions, starting with zero
	void prepare_next_dimension(int dimension) {
		if (comparator != nullptr) delete comparator;
		comparator = new rips_filtration_comparator<DistanceMatrix>(distance_matrix, dimension, binomial_coeff);

		if (next_comparator != nullptr) delete next_comparator;
		next_comparator =
		    new rips_filtration_comparator<DistanceMatrix>(distance_matrix, dimension + 1, binomial_coeff);

		if (dimension > 0) current_dimension++;
	}

	inline simplex_coboundary_enumerator<DistanceMatrix> coboundary(filtration_entry_t cell) {
		return simplex_coboundary_enumerator<DistanceMatrix>(cell, current_dimension, index_t(n), modulus, distance_matrix, binomial_coeff);
	}

	bool is_top_dimension() { return current_dimension >= max_dimension; }

  size_t top_dimension() { return current_dimension; }

	// Partial result output
	void computation_result(int, long long, long long) {}
	void finished() {}

	// compute_pairs(columns_to_reduce, pivot_column_index, dim, n, threshold,
	// modulus, multiplicative_inverse, dist, comp, comp_prev, binomial_coeff);

	// assemble_columns_to_reduce(columns_to_reduce, pivot_column_index, comp,
	// dim, n, threshold, binomial_coeff);
};

// template <typename DistanceMatrix, typename ComparatorCofaces, typename
// Comparator> void compute_pairs(std::vector<filtration_index_t>&
// columns_to_reduce, hash_map<index_t, index_t>& pivot_column_index,
//  index_t dim, index_t n, value_t threshold, coefficient_t modulus,
//  const std::vector<coefficient_t>& multiplicative_inverse, const
//  DistanceMatrix& dist, const ComparatorCofaces& comp, const Comparator&
//  comp_prev, const binomial_coeff_table& binomial_coeff) {

enum file_format { LOWER_DISTANCE_MATRIX, UPPER_DISTANCE_MATRIX, DISTANCE_MATRIX, POINT_CLOUD, DIPHA };

template <typename T> T read(std::istream& s) {
	T result;
	s.read(reinterpret_cast<char*>(&result), sizeof(T));
	return result; // on little endian: boost::endian::little_to_native(result);
}

compressed_lower_distance_matrix read_point_cloud(std::istream& input_stream) {
	std::vector<std::vector<value_t>> points;

	std::string line;
	value_t value;
	while (std::getline(input_stream, line)) {
		std::vector<value_t> point;
		std::istringstream s(line);
		while (s >> value) {
			point.push_back(value);
			s.ignore();
		}
		if (!point.empty()) points.push_back(point);
		assert(point.size() == points.front().size());
	}

	euclidean_distance_matrix eucl_dist(std::move(points));

	index_t n = index_t(eucl_dist.size());

	std::cout << "point cloud with " << n << " points in dimension " << eucl_dist.points.front().size() << std::endl;

	std::vector<value_t> distances;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_upper_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
}

compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;

	std::string line;
	value_t value;
	for (int i = 0; std::getline(input_stream, line); ++i) {
		std::istringstream s(line);
		for (int j = 0; j < i && s >> value; ++j) {
			distances.push_back(value);
			s.ignore();
		}
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_dipha(std::istream& input_stream) {
	if (read<int64_t>(input_stream) != 8067171840) {
		std::cerr << "input is not a Dipha file (magic number: 8067171840)" << std::endl;
		exit(-1);
	}

	if (read<int64_t>(input_stream) != 7) {
		std::cerr << "input is not a Dipha distance matrix (file type: 7)" << std::endl;
		exit(-1);
	}

	index_t n = index_t(read<int64_t>(input_stream));

	std::vector<value_t> distances;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (i > j)
				distances.push_back(value_t(read<double>(input_stream)));
			else
				read<double>(input_stream);

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream, file_format format) {
	switch (format) {
	case LOWER_DISTANCE_MATRIX:
		return read_lower_distance_matrix(input_stream);
	case UPPER_DISTANCE_MATRIX:
		return read_upper_distance_matrix(input_stream);
	case DISTANCE_MATRIX:
		return read_distance_matrix(input_stream);
	case POINT_CLOUD:
		return read_point_cloud(input_stream);
	case DIPHA:
		return read_dipha(input_stream);
  default:
    throw std::logic_error("Format not known.");
	}
}

int main(int argc, char** argv) {
	file_format format = DISTANCE_MATRIX;

	index_t dim_max = 1;
	index_t dim_min = 0;
	value_t threshold = std::numeric_limits<value_t>::max();

#ifdef USE_COEFFICIENTS
	coefficient_t modulus = 2;
#else
	const coefficient_t modulus = 2;
#endif

	auto arguments = parse_arguments(argc, argv);
	auto positional_arguments = get_positional_arguments(arguments);
	auto named_arguments = get_named_arguments(arguments);

	if (named_arguments.find("help") != named_arguments.end()) { print_usage_and_exit(-1); }
  if (positional_arguments.size() == 0) print_usage_and_exit(-1);

	named_arguments_t::const_iterator it;
	if ((it = named_arguments.find("format")) != named_arguments.end()) {
			if (it->second == std::string("lower-distance"))
				format = LOWER_DISTANCE_MATRIX;
			else if (it->second == std::string("upper-distance"))
				format = UPPER_DISTANCE_MATRIX;
			else if (it->second == std::string("distance"))
				format = DISTANCE_MATRIX;
			else if (it->second == std::string("point-cloud"))
				format = POINT_CLOUD;
			else if (it->second == std::string("dipha"))
				format = DIPHA;
      else {
        std::cerr << "The input format " << format << " is not supported." << std::endl;
        print_usage_and_exit(-1);
      }
  }

	if ((it = named_arguments.find("max-dim")) != named_arguments.end()) {
			std::string parameter = std::string(it->second);
			size_t next_pos;
			dim_max = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
  }

	if ((it = named_arguments.find("threshold")) != named_arguments.end()) {
			std::string parameter = std::string(it->second);
			size_t next_pos;
			threshold = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
  }

#ifdef USE_COEFFICIENTS
	if ((it = named_arguments.find("modulus")) != named_arguments.end()) {
			std::string parameter = std::string(it->second);
      size_t next_pos;
      modulus = std::stol(parameter, &next_pos);
      if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
  }
#endif

	size_t max_entries = std::numeric_limits<size_t>::max();
	if ((it = named_arguments.find("approximate")) != named_arguments.end()) { max_entries = atoi(it->second); }


	std::ifstream file_stream(positional_arguments[0]);
	if (positional_arguments[0] && file_stream.fail()) {
		std::cerr << "couldn't open file " << positional_arguments[0] << std::endl;
		exit(-1);
	}

	compressed_lower_distance_matrix dist = read_file(positional_arguments[0] ? file_stream : std::cin, format);

	index_t n = index_t(dist.size());

	std::cout << "distance matrix with " << n << " points" << std::endl;

	auto value_range = std::minmax_element(dist.distances.begin(), dist.distances.end());
	std::cout << "value range: [" << *value_range.first << "," << *value_range.second << "]" << std::endl;

	dim_max = std::min(dim_max, n - 2);


	vietoris_rips_complex_t<decltype(dist)> vietoris_rips_complex(dist, dim_max, modulus);
	auto output = get_output<decltype(vietoris_rips_complex)>(named_arguments);
  output.set_complex(&vietoris_rips_complex);

	persistence_computer_t<decltype(vietoris_rips_complex)> persistence_computer(vietoris_rips_complex, output, max_entries, modulus, threshold);
	persistence_computer.compute_persistence(dim_min, dim_max, false);
	output.print_aggregated_results();
}
