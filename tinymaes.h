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

#ifndef __TINYMAES_HEADER
#define __TINYMAES_HEADER
#include <stdint.h>

enum TINYMAES_WEIGHTS {MAES_W_EQUALS = 0, MAES_W_LINEAR, MAES_W_SUPERLINEAR};

// TODO: Add an alternative for 32-bits systems
#include "mt64.h"
#define MAES_RANDOM_S MT64_S
#define MAES_RANDOM_INIT MT64_init
#define MAES_RANDOM_NORMAL MT64_random_normal
#define MAES_RANDOM_REAL MT64_random_double
#define MAES_RANDOM_INT MT64_random64
#define MAES_RANDOM_MAX MT64_RANDOM_MAX

typedef struct {
  // Dimension of the optimiser
  int nDim, mu, lambda, nStep;
  // Optimiser parameters
  double mueff, chiN, cs, cw, c1, dSigma;
  // Standart deviation
  double sigma;
  // State of the random number generator
  MAES_RANDOM_S rnd;
  // State of the optimizer
  double *M; // Square of the covariance matrix [nDim][nDim]
             // Stored line wise
  double *X; // Matrix of individuals [lambda][nDim]
             // Stored column wise
  double *Z; // Vector generated by the gaussian law such as X = M*Z 
             // [lambda][nDim]
  double *X0; // Mean of the multivariate gaussian random law [nDim]
  double *ps; //  Isotropically path (that comes from the best mu element of Z) [nDim]
  double *weights; // weights puts on the best mu element [mu]
} TINYMAES_S;

// Allocate and init the optimiser structure
TINYMAES_S *TINYMAES_Create(int nDim, int lambda, int mu, int weights, uint64_t seed);

// Free the structure
void TINYMAES_Free(TINYMAES_S *maes);

// Reset the covariance matrix, ps and X0 to their initial values
void TINYMAES_Reset(TINYMAES_S *maes);

// Update population
//   idx[mu] -> best individiual ordered from the best to the worst
//              (if idx == NULL, it regenerates X using the same mean and covariance)
// Return the vectors of X
const double *TINYMAES_NextStep(TINYMAES_S *maes, int *orderIdx);

// Resample a selected interval of individuals (for handling a rare constraint)
const double *TINYMAES_Resample(TINYMAES_S *maes, int id0, int nIds);

// Set initial start point
void TINYMAES_SetX0(TINYMAES_S *maes, double *x0);

#endif

