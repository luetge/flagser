#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>

// #define USE_GOOGLE_HASHMAP

#ifndef MANY_VERTICES
// Assume that we have at most 65k vertices, and that there are at most ~2 billion cells
typedef unsigned short vertex_index_t;
typedef int32_t index_t;
#else
// Assume that we have at most 65k vertices, and that there are at most ~2 billion cells
typedef uint32_t vertex_index_t;
typedef int64_t index_t;
#endif

typedef float value_t;
typedef int16_t coefficient_t;

// An issue was discover when porting the code on windows
// 1UL on windows is ALWAYS 32 bits, when on unix systems is pointer size
// See https://stackoverflow.com/a/57213178
// Portable solutio not have a value of ptr size
#define ONE_ uintptr_t(1)

namespace std_algorithms = std;

#ifdef USE_GOOGLE_HASHMAP
#include <sparsehash/sparse_hash_map>
template <class Key, class T, class Hash = std::hash<Key>, class Pred = std::equal_to<Key>>
class hash_map : public google::sparse_hash_map<Key, T, Hash, Pred> {
public:
	inline void reserve(size_t hint) { this->resize(hint); }
};
#else
template <class Key, class T, class Hash = std::hash<Key>, class Pred = std::equal_to<Key>>
class hash_map : public std::unordered_map<Key, T, Hash, Pred> {};
#endif
// If need arises, replace this by a faster hash map.
// This is only used for moderately sized maps, therefore it is separate from
// the other hash map.
template <class Key, class T> class fast_hash_map : public std::unordered_map<Key, T> {};

// RAM reporting
#ifdef __APPLE__
#include <mach/mach.h>
void print_ram_usage(std::string msg = "") {
	std::cout << "\033[K";
	struct task_basic_info t_info;
	mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) {
		std::cout << "Could not determine used ram." << std::endl;
		return;
	}

	std::cout << "RAM usage";
	if (msg.size() > 0) std::cout << " [" << msg << "]";
	std::cout << ": " << printf("%.2f", t_info.resident_size / 1000000000.0) << " GB" << std::endl;
}
#else
void print_ram_usage(std::string) { std::cerr << "Sorry, can only print ram on MacOS for now." << std::endl; }
#endif

#ifndef USE_COEFFICIENTS
const coefficient_t modulus = 2;
#endif
