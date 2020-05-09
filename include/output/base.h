#pragma once

#include <fstream>
#include <string>

#include "../argparser.h"
#include "../definitions.h"

template <typename Complex> class output_t {
public:
	// Due to polymorphism, it's necessary to implement virtual destructor
	virtual ~output_t(){};
	virtual void set_complex(Complex*){};
	virtual void print(std::string){};
	virtual void finished(bool){};
	virtual void computing_barcodes_in_dimension(unsigned short){};
	virtual void new_barcode(value_t, value_t){};
	virtual void new_infinite_barcode(value_t){};
	virtual void skipped_column(value_t){};
	virtual void betti_number(size_t, size_t){};
	virtual void remaining_homology_is_trivial(){};
	virtual void print_aggregated_results(){};
};

template <typename Complex> class file_output_t : public output_t<Complex> {
protected:
	std::ofstream outstream;

public:
	file_output_t(const named_arguments_t& named_arguments) {
		auto filename =
		    get_argument_or_fail(named_arguments, "out", "Please provide an output filename via \"--out filename\" .");

		std::ifstream f(filename);
		if (f.good()) {
			std::cerr << "The output file already exists, aborting." << std::endl;
			exit(-1);
		}

		outstream.open(filename);
		if (outstream.fail()) {
			std::cerr << "couldn't open file " << filename << std::endl;
			exit(-1);
		}
	}

	~file_output_t() {
		if (outstream.is_open()) outstream.close();
	}

	virtual void print(std::string s) { outstream << s; }
};
