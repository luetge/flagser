
SHELL      := bash
GCC        := g++
hdf5_error := $(shell $(GCC) -lhdf5 2>&1 >/dev/null | grep -c "hdf5")
openmp_error := $(shell $(GCC) -fopenmp src/flagser.cpp -o ____flagser.o 2>&1 >/dev/null | grep -c "fopenmp")
_          := $(shell rm -f ____flagser.o)
COMPILE    := $(GCC) -std=c++11 -pthread
DEBUG      := -g
PRODUCTION := -Ofast -D NDEBUG
DETECTED_OS := $(shell sh -c 'uname -s 2>/dev/null || echo not')
ifeq (0, ${hdf5_error})
	COMPILE := ${COMPILE} -lhdf5 -D WITH_HDF5
	HDF5_MSG := "Found the HDF5-library, enabled reading .h5-files."
else
	HDF5_MSG := "\n\n\x1b[31;01m=========================================================================\n| I could not find the HDF5-library on your computer, so I disabled the |\n| functionality to read .h5-files. To re-enable it, please install the  |\n| HDF5-library (https://support.hdfgroup.org/HDF5/) and re-run \"make\".  |\n=========================================================================\x1b[0m\n\n"
ifeq ($(DETECTED_OS),Darwin)  # MacOS
	HDF5_MSG := "\n\n\x1b[31;01m=========================================================================\n| I could not find the HDF5-library on your computer, so I disabled the |\n| functionality to read .h5-files. To re-enable it, please install the  |\n| HDF5-library (https://support.hdfgroup.org/HDF5/) and re-run \"make\".  |\n|                                                                       |\n| On MacOS, the easiest way to do this is installing Homebrew           |\n| (https://brew.sh) and running \"brew install hdf5\" in the terminal.    |\n=========================================================================\x1b[0m\n\n"
endif
endif

ifeq (0, ${openmp_error})
	COMPILE := ${COMPILE} -fopenmp -D OPENMP_SUPPORT
	OPENMP_MSG := "Found OpenMP support, enabled advanced multithreading."
else
	OPENMP_MSG := "\n\n\x1b[31;01m=========================================================================\n| This compiler does not support OpenMP, to increase performance        |\n| install e.g. g++6. To select the compiler, set the Makefile variable  |\n| GCC appropriately, e.g. \"make GCC=g++-6\"                              |\n=========================================================================\x1b[0m\n\n"
endif

CONSENT_WARNING="\n\n\x1b[31;01mWARNING: Please make sure to install all python libraries as listed\nin requirements.txt, otherwise the compilation might fail.\nIf you change your mind about automatic installation,\nrun \"make reset_consent\".\x1b[0m\n\n"

build: er flagser flagser-memory flagser-count pathhom pathhom-count

all: er flagser flagser-memory flagser-debug flagser-count ripser flagser-coefficients flagser-coefficients-memory pathhom pathhom-count

reset_consent:
	@rm -f .install-consent

install_consent:
	@if [ ! -f .install-consent ]; then read -n 1 -p "Are you fine with me installing the required python dependencies for you (y/n)? " answer; if echo "$$answer" | grep -iq "^y" ;then echo "yes" > .install-consent && echo ""; else echo "no" > .install-consent; echo -e ${CONSENT_WARNING}; fi fi

install_libs: install_consent notify_hdf5 notify_openmp
	@if [[ $$(< .install-consent) == yes* ]]; then pip install --user -r requirements.txt --disable-pip-version-check | (grep -v 'already satisfied' || true); fi

notify_hdf5:
	@echo -e ${HDF5_MSG}

notify_openmp:
	@echo -e ${OPENMP_MSG}

flagser: src/flagser.cpp $(wildcard include/*) algorithms.math install_libs
	@echo "Compiling \"flagser\"." && python parser/math_parser.py && ${COMPILE} ${PRODUCTION} src/flagser.cpp -o flagser

flagser-coefficients: src/flagser.cpp $(wildcard include/*) algorithms.math install_libs
	@echo "Compiling \"flagser-coefficients\"." && python parser/math_parser.py && ${COMPILE} ${PRODUCTION} src/flagser.cpp -o flagser-coefficients -D USE_COEFFICIENTS

# Force h5 compilation
flagser-h5: src/flagser.cpp $(wildcard include/*) algorithms.math
	@echo "Compiling \"flagser\"." && python parser/math_parser.py && ${COMPILE} ${PRODUCTION} src/flagser.cpp -o flagser -lhdf5 -D WITH_HDF5

flagser-memory: src/flagser.cpp $(wildcard include/*) algorithms.math install_libs
	@echo "Compiling \"flagser-memory\"." && python parser/math_parser.py && ${COMPILE} ${PRODUCTION} src/flagser.cpp -o flagser-memory -D KEEP_FLAG_COMPLEX_IN_MEMORY

# Force h5 compilation
flagser-memory-h5: src/flagser.cpp $(wildcard include/*) algorithms.math
	@echo "Compiling \"flagser-memory\"." && python parser/math_parser.py && ${COMPILE} ${PRODUCTION} src/flagser.cpp -o flagser -lhdf5 -D WITH_HDF5 -D KEEP_FLAG_COMPLEX_IN_MEMORY
	
flagser-coefficients-memory: src/flagser.cpp $(wildcard include/*) algorithms.math install_libs
	@echo "Compiling \"flagser-coefficients-memory\"." && python parser/math_parser.py && ${COMPILE} ${PRODUCTION} src/flagser.cpp -o flagser-coefficients-memory -D KEEP_FLAG_COM    PLEX_IN_MEMORY -D USE_COEFFICIENTS

flagser-debug: src/flagser.cpp $(wildcard include/*) algorithms.math install_libs
	@echo "Compiling \"flagser-debug\"." && python parser/math_parser.py && ${COMPILE} ${DEBUG} src/flagser.cpp -o flagser-debug

flagser-count: src/flagser-count.cpp $(wildcard include/*) install_libs
	@echo "Compiling \"flagser-count\"." && ${COMPILE} ${PRODUCTION} src/flagser-count.cpp -o flagser-count

pathhom: src/flagser.cpp $(wildcard include/*) install_libs
	@echo "Compiling \"pathhom\"." && ${COMPILE} ${PRODUCTION} src/flagser.cpp -o pathhom -D PATH_HOMOLOGY

pathhom-count: src/flagser-count.cpp $(wildcard include/*) install_libs
	@echo "Compiling \"pathhom-count\"." && ${COMPILE} ${PRODUCTION} src/flagser-count.cpp -o pathhom-count -D PATH_HOMOLOGY

er: src/er.cpp
	@echo "Compiling \"er\"." && ${COMPILE} ${PRODUCTION} src/er.cpp -o tools/er

ripser: src/ripser.cpp $(wildcard include/*)
	@echo "Compiling \"ripser\"." && ${COMPILE} ${PRODUCTION} src/ripser.cpp -o ripser

clean:
	rm -f flagser flagser-memory flagser-debug tools/er ripser flagser-count pathhom pathhom-count flagser-coefficients flagser-coefficients-memory
