# flagser

Copyright © 2017–2018 Daniel Lütgehetmann.

### Description

flagser computes the homology of directed flag complexes. For a more detailed
description consult the documentation under `docs/documentation_flagser.pdf`.

### Building

Flagser requires a C++14 compiler and cmake. Here is how to obtain, build, and
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

Included in this package is also the original program `ripser`, modified only
in that the features of `flagser` are supported.

### Euler characteristic and cell counts

To only compute the Euler characteristic and cell counts, run

```sh
./flagser-count example.flag
```

### License

flagser is licensed under the MIT license (COPYING.txt), with an extra clause (CONTRIBUTING.txt) clarifying the license for modifications released without an explicit written license agreement. Please contact the author if you want to use Ripser in your software under a different license.
