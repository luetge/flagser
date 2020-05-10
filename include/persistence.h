#pragma once

// #define INDICATE_PROGRESS
// #define USE_COEFFICIENTS
// #define USE_GOOGLE_HASHMAP

#include <cassert>
#include <deque>
#include <iostream>
#include <queue>

#include "definitions.h"
#include "output/base.h"

#ifdef USE_ARRAY_HASHMAP
typedef std::deque<index_t> pivot_column_index_t;
const index_t INVALID_INDEX = std::numeric_limits<index_t>::max();
#else
typedef fast_hash_map<index_t, index_t> pivot_column_index_t;
#endif

std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m) {
	std::vector<coefficient_t> inverse(m);
	inverse[1] = 1;
	// m = a * (m / a) + m % a
	// Multipying with inverse(a) * inverse(m % a):
	// 0 = inverse(m % a) * (m / a)  + inverse(a)  (mod m)
	for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
	return inverse;
}

#ifdef USE_COEFFICIENTS
// TODO: Make the packed attribute work again
struct entry_t {
	// index_t index : 8 * (sizeof(index_t) - sizeof(coefficient_t));
	index_t index;
	coefficient_t coefficient;
	entry_t(index_t _index, coefficient_t _coefficient) : index(_index), coefficient(_coefficient) {}
	entry_t(index_t _index) : index(_index), coefficient(1) {}
	entry_t() : index(0), coefficient(1) {}
};
// } __attribute__((packed));

// static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

entry_t make_entry(index_t _index, coefficient_t _coefficient) { return entry_t(_index, _coefficient); }
index_t get_index(entry_t e) { return e.index; }
index_t get_coefficient(entry_t e) { return e.coefficient; }
void set_coefficient(entry_t& e, const coefficient_t c) { e.coefficient = c; }

bool operator==(const entry_t& e1, const entry_t& e2) {
	return get_index(e1) == get_index(e2) && get_coefficient(e1) == get_coefficient(e2);
}

std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
	stream << get_index(e) << ":" << get_coefficient(e);
	return stream;
}

#else

typedef index_t entry_t;
index_t get_index(entry_t i) { return i; }
index_t get_coefficient(entry_t) { return 1; }
entry_t make_entry(index_t _index, coefficient_t) { return entry_t(_index); }
void set_coefficient(index_t&, const coefficient_t) {}

#endif

const entry_t& get_entry(const entry_t& e) { return e; }

template <typename Entry> struct smaller_index {
	bool operator()(const Entry& a, const Entry& b) { return get_index(a) < get_index(b); }
};

class filtration_index_t : public std::pair<value_t, index_t> {
public:
	filtration_index_t() : std::pair<value_t, index_t>() {}
	filtration_index_t(std::pair<value_t, index_t> p) : std::pair<value_t, index_t>(p) {}
};
value_t get_filtration(filtration_index_t i) { return i.first; }
index_t get_index(filtration_index_t i) { return i.second; }

class filtration_entry_t : public std::pair<value_t, entry_t> {
public:
	filtration_entry_t(std::pair<value_t, entry_t> p) : std::pair<value_t, entry_t>(p) {}
	filtration_entry_t(entry_t e) : std::pair<value_t, entry_t>(value_t(0), e) {}
	filtration_entry_t() : filtration_entry_t(0) {}
	filtration_entry_t(value_t _filtration, index_t _index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(_filtration, make_entry(_index, _coefficient)) {}
	filtration_entry_t(filtration_index_t _filtration_index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(get_filtration(_filtration_index),
	                                  make_entry(get_index(_filtration_index), _coefficient)) {}
	filtration_entry_t(filtration_index_t _filtration_index) : filtration_entry_t(_filtration_index, 1) {}
};

const entry_t& get_entry(const filtration_entry_t& p) { return p.second; }
entry_t& get_entry(filtration_entry_t& p) { return p.second; }
index_t get_index(const filtration_entry_t& p) { return get_index(get_entry(p)); }
coefficient_t get_coefficient(const filtration_entry_t& p) { return get_coefficient(get_entry(p)); }
const value_t& get_filtration(const filtration_entry_t& p) { return p.first; }
void set_coefficient(filtration_entry_t& p, const coefficient_t c) { set_coefficient(get_entry(p), c); }

template <typename Entry> struct greater_filtration_or_smaller_index {
	bool operator()(const Entry& a, const Entry& b) {
		return (get_filtration(a) > get_filtration(b)) ||
		       ((get_filtration(a) == get_filtration(b)) && (get_index(a) < get_index(b)));
	}
};

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

class filtered_union_find {
	std::vector<index_t> parent;
	std::vector<uint8_t> rank;
	const std::vector<value_t> filtration;

public:
	filtered_union_find(const std::vector<value_t>& _filtration)
	    : parent(_filtration.size()), rank(_filtration.size(), 0), filtration(_filtration) {
		for (size_t i = 0ul; i < _filtration.size(); ++i) parent[i] = index_t(i);
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
		if (filtration[x] < filtration[y] || (filtration[x] == filtration[y] && rank[x] > rank[y]))
			parent[y] = x;
		else {
			parent[x] = y;
			if (rank[x] == rank[y]) ++rank[y];
		}
	}
};

template <typename Heap>
filtration_entry_t pop_pivot(Heap& column
#ifdef USE_COEFFICIENTS
                             ,
                             coefficient_t modulus
#endif
) {
	if (column.empty())
		return filtration_entry_t(-1);
	else {
		auto pivot = column.top();

#ifdef USE_COEFFICIENTS
		coefficient_t coefficient = 0;
		do {
			coefficient = (coefficient + get_coefficient(column.top())) % modulus;
			column.pop();

			if (coefficient == 0) {
				if (column.empty())
					return filtration_entry_t(-1);
				else
					pivot = column.top();
			}
		} while (!column.empty() && get_index(column.top()) == get_index(pivot));
		if (get_index(pivot) != -1) { set_coefficient(pivot, coefficient); }
#else
		column.pop();
		while (!column.empty() && get_index(column.top()) == get_index(pivot)) {
			column.pop();
			if (column.empty())
				return filtration_entry_t(-1);
			else {
				pivot = column.top();
				column.pop();
			}
		}
#endif
		return pivot;
	}
}

template <typename Heap>
filtration_entry_t get_pivot(Heap& column
#ifdef USE_COEFFICIENTS
                             ,
                             coefficient_t modulus
#endif
) {
	filtration_entry_t result = pop_pivot(column
#ifdef USE_COEFFICIENTS
	                                      ,
	                                      modulus
#endif
	);
	if (get_index(result) != -1) column.push(result);
	return result;
}

template <typename ValueType> class compressed_sparse_matrix {
	std::deque<size_t> bounds;
	std::deque<ValueType> entries;

public:
	size_t size() const { return bounds.size(); }

	void clear() {
		bounds.clear();
		bounds.shrink_to_fit();
		entries.clear();
		entries.shrink_to_fit();
	}

	typename std::deque<ValueType>::const_iterator cbegin(size_t index) const {
		assert(index < size());
		return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
	}

	typename std::deque<ValueType>::const_iterator cend(size_t index) const {
		assert(index < size());
		return entries.cbegin() + bounds[index];
	}

	template <typename Iterator> void append_column(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) { entries.push_back(*it); }
		bounds.push_back(entries.size());
	}

	void append_column() { bounds.push_back(entries.size()); }

	void push_back(ValueType e) {
		assert(0 < size());
		entries.push_back(e);
		++bounds.back();
	}

	void pop_back() {
		assert(0 < size());
		entries.pop_back();
		--bounds.back();
	}

	template <typename Collection> void append_column(const Collection collection) {
		append_column(collection.cbegin(), collection.cend());
	}
};

template <typename Heap> void push_entry(Heap& column, index_t i, coefficient_t c, value_t filtration) {
	entry_t e = make_entry(i, c);
	column.push(std::make_pair(filtration, e));
}

// This class is just an ordinary priority queue, but once the
// queue gets too long (because a lot of faces are inserted multiple
// times) it starts collecting the coefficients and only inserting each
// new face once
template <class Container, class Comparator>
class priority_queue_t : public std::priority_queue<filtration_entry_t, Container, Comparator> {
#ifdef USE_COEFFICIENTS
	std::unordered_map<index_t, coefficient_t> coefficients;
#else
	std::unordered_map<index_t, bool> coefficients;
#endif
	static const filtration_entry_t dummy;
	// TODO: Enable dynamic switching to the dense version for coefficients
#ifdef USE_COEFFICIENTS
	bool use_dense_version = true;
#else
	bool use_dense_version = false;
#endif
	coefficient_t modulus;
	size_t dense_threshold;

public:
	priority_queue_t(coefficient_t _modulus, size_t _dense_threshold)
	    : modulus(_modulus), dense_threshold(_dense_threshold) {}

	void push(const filtration_entry_t& value) {
		if (use_dense_version) {
			// If we already have this value: update the count and don't push it again
			auto p = coefficients.find(get_index(value));
			if (p != coefficients.end()) {
#ifdef USE_COEFFICIENTS
				p->second = (p->second + get_coefficient(value)) % modulus;
#else
				p->second = !p->second;
#endif
				return;
			}
		}

		std::priority_queue<filtration_entry_t, Container, Comparator>::push(value);

		if (use_dense_version) coefficients.insert(std::make_pair(get_index(value), get_coefficient(value)));

#ifndef USE_COEFFICIENTS
		if (!use_dense_version &&
		    std::priority_queue<filtration_entry_t, Container, Comparator>::size() >= dense_threshold)
			use_dense_version = true;
#endif
	}

	void pop() {
		// Don't use this, only allow get_pivot
		throw std::exception();
	}

	filtration_entry_t pop_pivot() {
		remove_trivial_coefficient_entries();
		if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
			return dummy;
		else {
			auto pivot = get_top();

#ifdef USE_COEFFICIENTS
			coefficient_t coefficient = 0;
			do {
				coefficient = (coefficient + get_coefficient(get_top())) % modulus;
				safe_pop();
				remove_trivial_coefficient_entries();

				if (coefficient == 0) {
					if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
						return dummy;
					else
						pivot = get_top();
				}
			} while (!std::priority_queue<filtration_entry_t, Container, Comparator>::empty() &&
			         get_index(get_top()) == get_index(pivot));
			if (get_index(pivot) != -1) { set_coefficient(pivot, coefficient); }
#else
			safe_pop();
			while (!std::priority_queue<filtration_entry_t, Container, Comparator>::empty() &&
			       get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()) ==
			           get_index(pivot)) {
				safe_pop();
				remove_trivial_coefficient_entries();

				if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
					return dummy;
				else {
					pivot = get_top();
					safe_pop();
				}
			}
#endif
			return pivot;
		}
	}

	filtration_entry_t get_pivot() {
		filtration_entry_t result = pop_pivot();
		if (get_index(result) != -1) { push(result); }
		return result;
	}

private:
	inline filtration_entry_t get_top() {
		auto pivot = std::priority_queue<filtration_entry_t, Container, Comparator>::top();

#ifdef USE_COEFFICIENTS
		if (use_dense_version) {
			auto e = coefficients.find(get_index(pivot));
			if (e != coefficients.end()) set_coefficient(pivot, e->second);
		}
#endif

		return pivot;
	}

	inline void safe_pop() {
		if (use_dense_version) {
			auto e =
			    coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			if (e != coefficients.end()) coefficients.erase(e);
		}
		std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
	}

	inline void remove_trivial_coefficient_entries() {
		if (use_dense_version) {
			if (std::priority_queue<filtration_entry_t, Container, Comparator>::size() == 0) return;
			auto p =
			    coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
#ifdef USE_COEFFICIENTS
			while (p != coefficients.end() && p->second % modulus == 0) {
#else
			while (p != coefficients.end() && p->second == false) {
#endif
				coefficients.erase(p);
				std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
				p = coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			}
		}
	}
};
template <class Container, class Comparator>
const filtration_entry_t
    priority_queue_t<Container, Comparator>::dummy(filtration_entry_t(std::make_pair(value_t(0.0f), entry_t(-1))));

#ifdef SORT_COLUMNS_BY_PIVOT
template <typename Entry, typename Complex> struct greater_filtration_or_better_pivot_or_smaller_index {
	greater_filtration_or_better_pivot_or_smaller_index(Complex& _complex) : complex(_complex) {}
	bool operator()(Entry a, Entry b) const {
		// First order by the filtration value
		if (get_filtration(a) > get_filtration(b)) return true;
		if (get_filtration(a) < get_filtration(b)) return false;

		const auto& ta = get_coboundary_size_and_gap_and_pivot(a);
		const auto& tb = get_coboundary_size_and_gap_and_pivot(b);

		// Then the number of non-trivial coboundary entries
		if (std::get<0>(ta) < std::get<0>(tb)) return true;
		if (std::get<0>(ta) > std::get<0>(tb)) return false;

		// Then order by the better pivoting
		if (std::get<2>(ta) < std::get<2>(tb)) return true;
		if (std::get<2>(ta) > std::get<2>(tb)) return false;

		if (std::get<1>(ta) > std::get<1>(tb)) return true;
		if (std::get<1>(ta) < std::get<1>(tb)) return false;

		// Finally, order by their indices
		return get_index(a) < get_index(b);
	}

private:
	Complex& complex;

	// A column is considered to be a better pivot if the jump from pivot to the next
	// non-trivial element is as big as possible. This prevents accidentally inserting
	// non-trivial elements just below the pivot, which sometimes creates very long
	// reduction chains.
	// The second sort criterium is for it to be small because the small pivots will be
	// used the most.
	std::tuple<size_t, size_t, index_t> get_coboundary_size_and_gap_and_pivot(filtration_entry_t a) const {
		// Look at the first two gaps of the pivot and the next element
		index_t pivot = 0;
		size_t gap_after_pivot = 0;
		auto iterator = complex.coboundary(a);
		size_t coboundary_size = 0;
		while (iterator.has_next()) {
			coboundary_size++;
			index_t next_index = get_index(iterator.next().second);
			if (next_index > pivot) {
				gap_after_pivot = next_index - pivot;
				pivot = next_index;
			}
		}

		return std::make_tuple(coboundary_size, gap_after_pivot, pivot);
	}
};
#endif

template <typename Complex> class persistence_computer_t {
private:
	Complex& complex;
	output_t<Complex>* output;
	value_t max_filtration;
	size_t max_entries;
	index_t euler_characteristic = 0;
	bool print_betti_numbers_to_console = true;

#ifdef USE_COEFFICIENTS
	coefficient_t modulus = 2;
#else
	const coefficient_t modulus = 2;
#endif
	std::vector<coefficient_t> multiplicative_inverse;
	std::deque<filtration_index_t> columns_to_reduce;
#ifdef RETRIEVE_PERSISTENCE
	std::vector<size_t> betti_numbers;
	std::vector<std::vector<std::pair<value_t, value_t>>> birth_deaths_by_dim;
	std::vector<size_t> cell_count;
#endif

public:
	persistence_computer_t(Complex& _complex, output_t<Complex>* _output,
	                       size_t _max_entries = std::numeric_limits<size_t>::max(), int _modulus = 2,
	                       value_t _max_filtration = std::numeric_limits<value_t>::max())
	    : complex(_complex), output(_output), max_filtration(_max_filtration), max_entries(_max_entries),
#ifdef USE_COEFFICIENTS
	      modulus(_modulus),
#endif
	      multiplicative_inverse(multiplicative_inverse_vector(modulus)) {
#ifndef USE_COEFFICIENTS
		if (_modulus != 2) {
			throw std::logic_error("If you want to use modulus != 2, please compile with coefficients.");
		}
#endif
	}

	void set_print_betti_numbers(bool print_betti_numbers) { print_betti_numbers_to_console = print_betti_numbers; }

	void compute_persistence(unsigned short min_dimension = 0,
	                         unsigned short max_dimension = std::numeric_limits<unsigned short>::max(),
	                         bool check_euler_characteristic = true) {
		compute_zeroth_persistence(min_dimension, max_dimension);
		compute_higher_persistence(min_dimension, max_dimension);
		complex.finished();
		output->finished(check_euler_characteristic);

		// Sanity check whether there were any problems computing the homology
		bool computed_full_homology = min_dimension == 0 && max_dimension == std::numeric_limits<unsigned short>::max();
		if (check_euler_characteristic && computed_full_homology && max_entries == std::numeric_limits<size_t>::max()) {
			index_t cell_euler_characteristic = 0;
			for (size_t i = 0; i <= complex.top_dimension(); i++) {
				cell_euler_characteristic += (i % 2 == 1 ? -1 : 1) * index_t(complex.number_of_cells(index_t(i)));
			}

			if (print_betti_numbers_to_console)
				std::cout << "The Euler characteristic is given by: " << cell_euler_characteristic << std::endl;

			if (cell_euler_characteristic != euler_characteristic) {
				std::cerr << "ERROR: The homological Euler characteristic (which is " << euler_characteristic
				          << ") differs from the cellular Euler characteristic (which is " << cell_euler_characteristic
				          << "), apparently there is an error in the program." << std::endl;
			}
		}

#ifdef RETRIEVE_PERSISTENCE
		for (size_t i = min_dimension; i < std::min(complex.top_dimension(), (size_t)max_dimension + 1); i++) {
			cell_count.push_back(index_t(complex.number_of_cells(index_t(i))));
		}
#endif
	}

#ifdef RETRIEVE_PERSISTENCE
	index_t get_euler_characteristic() { return euler_characteristic; }

	std::vector<size_t> get_betti_numbers() { return betti_numbers; }

	size_t get_betti_numbers(size_t dimension) { return betti_numbers[dimension]; }

	std::vector<std::vector<std::pair<value_t, value_t>>> get_persistence_diagram() { return birth_deaths_by_dim; }

	std::vector<std::pair<value_t, value_t>> get_persistence_diagram(size_t dimension) {
		return birth_deaths_by_dim[dimension];
	}

	std::vector<size_t> get_cell_count() { return cell_count; }

	size_t get_cell_count(size_t dimension) { return cell_count[dimension]; }
#endif

protected:
	void compute_zeroth_persistence(unsigned short min_dimension, unsigned short) {
		complex.prepare_next_dimension(0);

		// Only compute this if we actually need it
		if (min_dimension > 1) return;

#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "computing persistent homology in dimension 0" << std::flush << "\r";
#endif

		long long betti_number = 0;
#ifdef RETRIEVE_PERSISTENCE
		std::vector<std::pair<value_t, value_t>> birth_death;
#endif
		size_t n = complex.number_of_cells(0);
		filtered_union_find dset(complex.vertex_filtration());
		std::vector<filtration_index_t> edges;
		index_t number_of_edges = index_t(complex.number_of_cells(1));
		for (index_t index = 0; index < number_of_edges; index++) {
			value_t filtration = complex.filtration(1, index);
			if (filtration <= max_filtration) edges.push_back(std::make_pair(filtration, index));
		}
		std_algorithms::sort(edges.rbegin(), edges.rend(), greater_filtration_or_smaller_index<filtration_index_t>());

		// Let the output class know that we are now computing zeroth degree barcodes
		output->computing_barcodes_in_dimension(0);

		for (auto e : edges) {
			const auto vertices = complex.vertices_of_edge(get_index(e));
			index_t u = dset.find(vertices.first), v = dset.find(vertices.second);

			if (u != v) {
				// Only output bars if we are interested in zeroth homology
				const auto filtration_u = complex.filtration(0, u);
				const auto filtration_v = complex.filtration(0, v);
				dset.link(u, v);
				if (min_dimension == 0 && get_filtration(e) > std::max(filtration_u, filtration_v)) {
					// Check which vertex is merged into which other vertex.
					const auto f = dset.find(u) == u ? filtration_v : filtration_u;
					output->new_barcode(f, get_filtration(e));
#ifdef RETRIEVE_PERSISTENCE
					birth_death.push_back(std::make_pair(f, get_filtration(e)));
#endif
				}
			} else {
				columns_to_reduce.push_back(e);
			}
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

		// If we don't care about zeroth homology, then we can stop here
		if (min_dimension == 1) return;

		for (index_t index = 0; index < index_t(n); ++index) {
			if (dset.find(index) == index) {
				output->new_infinite_barcode(complex.filtration(0, index));
				betti_number++;
#ifdef RETRIEVE_PERSISTENCE
				birth_death.push_back(
				    std::make_pair(complex.filtration(0, index), std::numeric_limits<value_t>::infinity()));
#endif
			}
		}

#ifdef RETRIEVE_PERSISTENCE
		betti_numbers.push_back(betti_number);
		birth_deaths_by_dim.push_back(birth_death);
#endif
		// Report the betti number back to the complex and the output
		complex.computation_result(0, betti_number, 0);
		output->betti_number(betti_number, 0);

		if (print_betti_numbers_to_console) {
			std::cout << "\033[K"
			          << "The dimensions of the mod-" << modulus << " homology of the full complex are:" << std::endl
			          << std::endl
			          << "dim H_0 = " << betti_number << std::endl;
		}
		euler_characteristic += index_t(betti_number);
	}

	void compute_higher_persistence(unsigned short min_dimension, unsigned short max_dimension) {
		for (auto dimension = 1u; dimension <= max_dimension; ++dimension) {
			complex.prepare_next_dimension(dimension);

			if (dimension + 1 == min_dimension) {
				// Here we need to reduce *all* cells because we did not compute anything of smaller dimension
				index_t number_of_cells = index_t(complex.number_of_cells(dimension));
				for (index_t index = 0; index < number_of_cells; index++) {
					columns_to_reduce.push_back(std::make_pair(complex.filtration(dimension, index), index));
				}
			}

			if (dimension + 1 < min_dimension) continue;

			output->computing_barcodes_in_dimension(dimension);

			sort_columns();

#ifdef INDICATE_PROGRESS
			std::cout << "\033[K"
			          << "computing persistent homology in dimension " << dimension << std::flush << "\r";
#endif
#ifdef USE_ARRAY_HASHMAP
			pivot_column_index_t pivot_column_index(complex.number_of_cells(dimension + 1), INVALID_INDEX);
#else
			pivot_column_index_t pivot_column_index;
			pivot_column_index.reserve(complex.number_of_cells(dimension + 1));
#endif

			auto betti = compute_pairs(dimension, pivot_column_index, dimension >= min_dimension);
			if (dimension >= min_dimension) {
				complex.computation_result(dimension, betti.first, betti.second);
#ifdef RETRIEVE_PERSISTENCE
				betti_numbers.push_back(betti.first);
#endif
				output->betti_number(betti.first, betti.second);
				euler_characteristic += (dimension & 1 ? -1 : 1) * betti.first;

				if (print_betti_numbers_to_console) {
					std::cout << "\033[K"
					          << "dim H_" << dimension << " = " << betti.first;
					if (betti.second > 0) { std::cout << " (skipped " << betti.second << ")"; }
					std::cout << std::endl;
				}
			} else if (int(dimension) == min_dimension - 1 && print_betti_numbers_to_console &&
			           max_entries < std::numeric_limits<size_t>::max()) {
				std::cout << "\033[K"
				          << "# Skipped columns in dimension " << dimension << ": " << betti.second << std::endl;
			}
			if (dimension < max_dimension) assemble_columns_to_reduce(dimension, pivot_column_index);

			// Stop early
			if (complex.is_top_dimension()) {
				output->remaining_homology_is_trivial();
				break;
			}
		}
	}

	void assemble_columns_to_reduce(index_t dimension, pivot_column_index_t& pivot_column_index) {
		index_t num_cells = index_t(complex.number_of_cells(dimension + 1));

		columns_to_reduce.clear();

#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "assembling " << num_cells << " columns" << std::flush << "\r";
#endif

		for (index_t index = 0; index < num_cells; ++index) {
			if (
#ifdef USE_ARRAY_HASHMAP
			    pivot_column_index[index] == INVALID_INDEX
#else
			    pivot_column_index.find(index) == pivot_column_index.end()
#endif
			) {
				value_t filtration = complex.filtration(dimension + 1, index);
				if (filtration <= max_filtration) { columns_to_reduce.push_back(std::make_pair(filtration, index)); }
#ifdef INDICATE_PROGRESS
				if ((index + 1) % 100000 == 0)
					std::cout << "\033[K"
					          << "assembled " << columns_to_reduce.size() << " out of " << (index + 1) << "/"
					          << num_cells << " columns" << std::flush << "\r";
#endif
			}
		}
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif
	}

	void sort_columns() {
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "sorting " << columns_to_reduce.size() << " columns" << std::flush << "\r";
#endif

#ifdef SORT_COLUMNS_BY_PIVOT
		std_algorithms::sort(
		    columns_to_reduce.begin(), columns_to_reduce.end(),
		    greater_filtration_or_better_pivot_or_smaller_index<filtration_index_t, decltype(complex)>(complex));
#else
		std_algorithms::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
		                     greater_filtration_or_smaller_index<filtration_index_t>());
#endif
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif
	}

	std::pair<index_t, index_t> compute_pairs(index_t, pivot_column_index_t& pivot_column_index,
	                                          bool generate_output = true) {
		index_t betti = 0;
		index_t betti_error = 0;
#ifdef INDICATE_PROGRESS
		auto verbose_logging_threshold = size_t(columns_to_reduce.size() * 0.90);
#endif

#ifdef ASSEMBLE_REDUCTION_MATRIX
		compressed_sparse_matrix<filtration_entry_t> reduction_coefficients;
#else
#ifdef USE_COEFFICIENTS
		std::vector<filtration_entry_t> reduction_coefficients;
#endif
#endif
#ifdef RETRIEVE_PERSISTENCE
		std::vector<std::pair<value_t, value_t>> birth_death;
#endif

		std::vector<filtration_entry_t> coface_entries;

		for (auto i = 0ul; i < columns_to_reduce.size(); ++i) {
			auto column_to_reduce = columns_to_reduce[i];

#ifdef ASSEMBLE_REDUCTION_MATRIX
			std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>>
			    reduction_column;
#endif

			priority_queue_t<std::deque<filtration_entry_t>, greater_filtration_or_smaller_index<filtration_entry_t>>
			    working_coboundary(modulus, columns_to_reduce.size());

			value_t filtration = get_filtration(column_to_reduce);

#ifdef INDICATE_PROGRESS
			if ((i + 1) % 10000 == 0 || (i >= verbose_logging_threshold && (i + 1) % 1000 == 0)) {
				std::cout << "\033[K"
				          << "reducing column " << i + 1 << "/" << columns_to_reduce.size() << " (filtration "
				          << filtration << ", infinite bars: " << betti;
				if (betti_error > 0) std::cout << " (skipped " << betti_error << ")";
				std::cout << ")" << std::flush << "\r";
			}
#endif

			index_t j = i;

			// start with a dummy pivot entry with coefficient -1 in order to initialize
			// working_coboundary with the coboundary of the simplex with index column_to_reduce
			filtration_entry_t pivot(0, -1, -1 + modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
			// initialize reduction_coefficients as identity matrix
			reduction_coefficients.append_column();
			reduction_coefficients.push_back(filtration_entry_t(column_to_reduce, 1));
#else
#ifdef USE_COEFFICIENTS
			reduction_coefficients.push_back(filtration_entry_t(column_to_reduce, 1));
#endif
#endif

#ifndef SKIP_APPARENT_PAIRS
			bool might_be_apparent_pair = (i == size_t(j));
#endif

			size_t iterations = 0;
			do {
				const coefficient_t factor = modulus - get_coefficient(pivot);

#ifdef ASSEMBLE_REDUCTION_MATRIX
				auto coeffs_begin = reduction_coefficients.cbegin(j), coeffs_end = reduction_coefficients.cend(j);
#else
#ifdef USE_COEFFICIENTS
				auto coeffs_begin = &reduction_coefficients[j], coeffs_end = &reduction_coefficients[j] + 1;
#else
				auto coeffs_begin = &columns_to_reduce[j], coeffs_end = &columns_to_reduce[j] + 1;
#endif
#endif
				for (auto it = coeffs_begin; it != coeffs_end; ++it) {
					filtration_entry_t cell = *it;
					set_coefficient(cell, get_coefficient(cell) * factor % modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
					reduction_column.push(cell);
#endif
					coface_entries.clear();

					auto coboundary = complex.coboundary(cell);

					while (coboundary.has_next()) {
						filtration_entry_t coface = coboundary.next();

						if (get_filtration(coface) <= max_filtration) {
#ifndef SKIP_APPARENT_PAIRS
							coface_entries.push_back(coface);
							if (might_be_apparent_pair && (get_filtration(cell) == get_filtration(coface))) {
#ifdef USE_ARRAY_HASHMAP
								if (pivot_column_index[get_index(coface)] == INVALID_INDEX)
#else
								if (pivot_column_index.find(get_index(coface)) == pivot_column_index.end())
#endif
								{
									pivot = coface;
									goto found_persistence_pair;
								}
								might_be_apparent_pair = false;
							}
#else
							iterations++;
							working_coboundary.push(coface);
#endif
						}
					}

#ifndef SKIP_APPARENT_PAIRS
					for (auto e : coface_entries) {
						iterations++;
						working_coboundary.push(e);
					}
#endif
				}

				if (iterations > max_entries) {
					// Abort, this is too expensive
					if (generate_output) output->skipped_column(filtration);
#ifdef RETRIEVE_PERSISTENCE
					birth_death.push_back(std::make_pair(filtration, std::numeric_limits<value_t>::signaling_NaN()));
#endif
					betti_error++;
					break;
				}

				pivot = working_coboundary.get_pivot();

				if (get_index(pivot) != -1) {
#ifdef USE_ARRAY_HASHMAP
					auto pivot_column_idx = pivot_column_index[get_index(pivot)];

					if (pivot_column_idx != INVALID_INDEX) {
						j = pivot_column_idx;
						continue;
					}
#else
					auto pair = pivot_column_index.find(get_index(pivot));

					if (pair != pivot_column_index.end()) {
						j = pair->second;
						continue;
					}
#endif
				} else {
					if (generate_output) {
						output->new_infinite_barcode(filtration);
						betti++;
#ifdef RETRIEVE_PERSISTENCE
						birth_death.push_back(std::make_pair(filtration, std::numeric_limits<value_t>::infinity()));
#endif
					}
					break;
				}

#ifndef SKIP_APPARENT_PAIRS
			found_persistence_pair:
#endif
				value_t death = get_filtration(pivot);
				if (generate_output && filtration != death) {
					output->new_barcode(filtration, death);
#ifdef RETRIEVE_PERSISTENCE
					birth_death.push_back(std::make_pair(filtration, death));
#endif
				}

#ifdef USE_ARRAY_HASHMAP
				pivot_column_index[get_index(pivot)] = i;
#else
				pivot_column_index.insert(std::make_pair(get_index(pivot), i));
#endif

#ifdef USE_COEFFICIENTS
				const coefficient_t inverse = multiplicative_inverse[get_coefficient(pivot)];
#endif

#ifdef ASSEMBLE_REDUCTION_MATRIX
				// replace current column of reduction_coefficients (with a single diagonal 1 entry)
				// by reduction_column (possibly with a different entry on the diagonal)
				reduction_coefficients.pop_back();
				while (true) {
#ifdef USE_COEFFICIENTS
					filtration_entry_t e = pop_pivot(reduction_column, modulus);
#else
					filtration_entry_t e = pop_pivot(reduction_column);
#endif
					if (get_index(e) == -1) break;
#ifdef USE_COEFFICIENTS
					set_coefficient(e, inverse * get_coefficient(e) % modulus);
					assert(get_coefficient(e) > 0);
#endif
					reduction_coefficients.push_back(e);
				}
#else
#ifdef USE_COEFFICIENTS
				reduction_coefficients.pop_back();
				reduction_coefficients.push_back(filtration_entry_t(column_to_reduce, inverse));
#endif
#endif
				break;
			} while (true);
		}

#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif
#ifdef RETRIEVE_PERSISTENCE
		birth_deaths_by_dim.push_back(birth_death);
#endif
		return std::make_pair(betti, betti_error);
	}
};
