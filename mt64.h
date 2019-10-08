#ifndef __SIMPLEMT_HEADER
#define __SIMPLEMT_HEADER

#include <stdint.h>

// Fast Mersenne Twister Random Number Generator optimized for 64-bits CPU
#define __MT19937_N 321
#define MT64_RANDOM_MAX 0xffffffffffffffff

typedef struct {
  // State of the Mersenne Twister Random Generator
  int pos;
  uint64_t state[__MT19937_N];
} MT64_S;

// Initialise SIMPLEMT with a 64-bits seed
void MT64_init(MT64_S *mt, uint64_t seed);

// Raw 64-bits integer number generation
uint64_t MT64_random64(MT64_S *mt);

// Double precision number generation 0 <= r < 1
double MT64_random_double(MT64_S *mt);

// Generate a vector of random variable using N(0,1)
void MT64_random_normal(MT64_S *mt, double *vect, int size_vect);


#endif

