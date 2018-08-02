# flagser

### Description

flagser computes the homology of directed flag complexes. For a more detailed
description consult the documentation under `docs/documentation_flagser.pdf`.

### Building

Flagser requires a C++11 compiler and python. Here is how to obtain, build, and
run flagser:

```sh
git clone git@gitlab.com:luetge/flagser.git
cd flagser
(mkdir -p build && cd build && cmake .. && make -j)
./test/run
```

### Running

To call flagser you run

```sh
./flagser --out example.homology test/medium-test-data.flag
```

For more detailed instructions, see `docs/documentation_flagser.pdf`. The
program `flagser-memory` is a variant of flagser storing the full directed flag
complex in memory, speeding up parts of the computation but requiring more
memory.

### Euler characteristic and cell counts

To only compute the Euler characteristic and cell counts, run

```sh
./flagser-count example.flag
```

### Caveats

Persistent homology can at the moment only be computed modulo two. We plan to
support the other finite fields in the future.

### License

flagser is licensed under the LGPL 3.0. Please contact the author if you want to use flagser in your software under a different license.
