#include <stdio.h>
#include "tinymaes.h"
#include "heapsort.h"

const int N = 30, SEED = 12345, nGMax = 1024;
const int lambda = 4*N*N, mu = 2*N*N;

static double frosenbrock(double *X) {
  double s, v;
  int i;
  v = X[0] - 1;
  s = v*v;
  for(i = 1; i < N; i++) {
    v = X[i-1]*X[i-1] - X[i];
    s += 100*v*v;
  }
}

int main(int argc, char **argv) {
  TINYMAES_S *maes = TINYMAES_Create(N, lambda, mu, MAES_W_SUPERLINEAR, SEED);
  int nGen, i, k, idx[mu];
  double *X = TINYMAES_NextStep(maes, NULL);
  double F[lambda];
  for(nGen = 0; nGen < nGMax; nGen++) {
    for(i = 0; i < lambda; i++)
      F[i] = frosenbrock(&X[i*N]);
    heapsort_mu(F, 1, lambda, idx, mu);
    k = idx[0];
    printf("%f: ", F[k]);
    for(i = 0; i < N; i++)
      printf("%f ", X[k*N+i]);
    printf("\n");
    X = TINYMAES_NextStep(maes, idx);
  }
  TINYMAES_Free(maes);
}



