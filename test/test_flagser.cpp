#define ASSEMBLE_REDUCTION_MATRIX
#define INDICATE_PROGRESS
#define SKIP_APPARENT_PAIRS
#define USE_ARRAY_HASHMAP
#define USE_CELLS_WITHOUT_DIMENSION
#define SORT_COLUMNS_BY_PIVOT
#define RETRIEVE_PERSISTENCE
// #define WITH_HDF5
// #define USE_COEFFICIENTS
// #define MANY_VERTICES


#include "./base.h"

int main(int, char**) {
  run_all(true);
  return 0;
}
