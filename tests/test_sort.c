#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tinymaes.h"
#include "heapsort.h"

#define LAMBDA 1635
#define MU 812
#define SEED 5555

// Compare heapsort_mu and qsort

static int _comp(const void *p1, const void *p2) {
  double *f1 = (double *) p1, *f2 = (double *) p2;
  if(*f1 < *f2)
    return -1;
  else if(*f1 > *f2)
    return 1;
  return 0;
}

int main(int argc, char **argv) {
  double v[LAMBDA], v2[LAMBDA];
  int idx[MU], i;
  MAES_RANDOM_S rnd;
  MAES_RANDOM_INIT(&rnd, SEED);
  MAES_RANDOM_NORMAL(&rnd, v, LAMBDA);
  memcpy((void*) v2, (void*) v, sizeof(double)*LAMBDA);
  qsort(v2, LAMBDA, sizeof(double), _comp);
  heapsort_mu(v, 1, LAMBDA, idx, MU);
  for(i = 0; i < MU; i++) {
    printf("%f %f\n", v[idx[i]], v2[i]);
    if(v[idx[i]] != v2[i]) {
      printf("Sort: Error\n");
      return -1;
    }
  }
  printf("Sort: Test OK\n");
  return 0;
}

