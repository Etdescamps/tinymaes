#include <math.h>
#include "mt64.h"

/* This code uses formulae from mt19937-64.c of Takuji Nishimura and Makoto Matsumoto
 * References:
 *  * T. Nishimura, "Tables of 64-bit Mersenne Twisters"
 *    ACM Transactions on Modeling and Computer Simulation 10. (2000) 348--357
 *  * M. Matsumoto and T. Nishimura, "Mersenne Twister: a 623-dimensionally equidistributed 
 *         uniform pseudorandom number generator"
 *         ACM Transactions on Modeling and  Computer Simulation 8. (Jan. 1998) 3--30.
 */

void MT64_init(MT64_S *mt, uint64_t seed) {
  int i;
  // Fill the state with a linear congruential generator
  for(i = 0; i < __MT19937_N; i++) {
    // From Knuth MMIX
    seed = 6364136223846793005ULL*seed + 1442695040888963407ULL;
    mt->state[i] = seed;
  }
  mt->pos = __MT19937_N; // Force a first twist before using MT
}

inline static uint64_t twist3(uint64_t s1, uint64_t s2, uint64_t sM) {
  const uint64_t matA[2] = {0, 0xB5026F5AA96619E9ULL};
  const uint64_t mostS = 0xFFFFFFFF80000000ULL; // Most significant 33 bits
  const uint64_t lessS = 0x7FFFFFFFULL; // Least significant 31 bits
  uint64_t x = (s1 & mostS) | (s2 & lessS);
  return sM^(x>>1)^matA[((int) x) & 0x1];
}

static void _twist(MT64_S *mt) {
  const int N_M = 165, M = 156, N = __MT19937_N;
  int i;
  for(i = 0; i < N_M; i++)
    mt->state[i] = twist3(mt->state[i], mt->state[i+1], mt->state[i+M]);
  for(i = N_M; i < N-1; i++) 
    mt->state[i] = twist3(mt->state[i], mt->state[i+1], mt->state[i-N_M]);
  mt->state[N-1] = twist3(mt->state[N-1], mt->state[0], mt->state[M-1]);
  mt->pos = 0;
}

inline static uint64_t _nextI64(MT64_S *mt) {
  uint64_t x;
  if(mt->pos >= __MT19937_N)
    _twist(mt);
  x = mt->state[mt->pos++];
  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);
  return x;
}

uint64_t MT64_random64(MT64_S *mt) {
  return _nextI64(mt);
}

inline static double _nextD(MT64_S *mt) {
  union {
    uint64_t u;
    double f;
  } y = {0x3ff0000000000000ULL | (_nextI64(mt) >> 10)};
  return y.f - 1.0;
}

double MT64_random_double(MT64_S *mt) {
  return _nextD(mt);
}

// Generate random normal variables using Box-Muller's method
void MT64_random_normal(MT64_S *mt, double *vect, int size_vect) {
  double dist, f1, f2, w;
  int i;
  for(i = 0; i < size_vect; i+=2) {
    do {
      f1 = 2*_nextD(mt)-1;
      f2 = 2*_nextD(mt)-1;
      dist = f1*f1+f2*f2;
    } while(dist >= 1.0); // Do until the point (f1,f2) is located within the unit circle
    w = sqrt(-2*log(dist)/dist);
    vect[i] = f1*w;
    if(i + 1 < size_vect)
      vect[i+1] = f2*w;
  }
}

