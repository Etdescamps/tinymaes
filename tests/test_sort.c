#include <stdio.h>
#include "tinymaes.h"
#include "heapsort.h"

#define DIM 512
#define MU 256
#define SEED 5555

int main(int argc, char **argv) {
  double v[DIM];
  int idx[MU], i;
  MAES_RANDOM_S rnd;
  MAES_RANDOM_INIT(&rnd, SEED);
  MAES_RANDOM_NORMAL(&rnd, v, DIM);
  heapsort_mu(v, 1, DIM, idx, MU);
  for(i = 0; i < MU; i++)
    printf("%f\n", v[idx[i]]);
}

