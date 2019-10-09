/* Copyright (c) 2019 Etienne Descamps <etdescdev@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * This test try to reproduce the results found in the article:
 * "Simplify Your Covariance Matrix Adaptation Evolution Strategy"
 * by Hans-Georg Beyer and Bernhard Sendhoff
 * published in IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 21 no. 5, OCTOBER 2017
 *
 * It tries to reproduce the results Figure 7 for the function Rosenbrock with N=30
 * and lambda = 4*N, mu = 2*N, with an initial X0 of (1,...1)
 */

#include <stdio.h>
#include "tinymaes.h"
#include "heapsort.h"

const int N = 30, SEED = 1234567, nGMax = 3000;
const int lambda = 4*N, mu = 2*N;
const double target = 1e-12;

static double frosenbrock(const double *X) {
  double s = 0, v;
  int i;
  for(i = 1; i < N; i++) {
    v = X[i-1]-1;
    s += v*v;
    v = X[i-1]*X[i-1] - X[i];
    s += 100*v*v;
  }
  return s;
}

int main(int argc, char **argv) {
  TINYMAES_S *maes = TINYMAES_Create(N, lambda, mu, MAES_W_SUPERLINEAR, SEED);
  int nGen, i, k, idx[mu];
  const double *X = TINYMAES_NextStep(maes, NULL);
  double X0[N];
  double F[lambda];
  for(i = 0; i < N; i++)
    X0[i] = 1;
  TINYMAES_SetX0(maes, X0);
  for(nGen = 0; nGen < nGMax; nGen++) {
    for(i = 0; i < lambda; i++)
      F[i] = frosenbrock(&X[i*N]);
    heapsort_mu(F, 1, lambda, idx, mu);
    k = idx[0];
    printf("%d %e: ", nGen, F[k]);
    for(i = 0; i < N; i++)
      printf("%f ", X[k*N+i]);
    printf("\n");
    if(F[k] < target) {
      TINYMAES_Free(maes);
      printf("Optimization with MA-ES: Test OK\n");
      return 0;
    }
    X = TINYMAES_NextStep(maes, idx);
  }
  TINYMAES_Free(maes);
  printf("Optimization with MA-ES: Does not converge in reasonable number of steps\n");
  return -1;
}



