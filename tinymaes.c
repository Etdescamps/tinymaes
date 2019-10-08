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

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tinymaes.h"

#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define MIN(a,b) ((a) > (b) ? (b) : (a))

const double _alphaCov = 2.0;

TINYMAES_S *TINYMAES_Create(int nDim, int lambda, int mu, int weights, uint64_t seed) {
  int i, j;
  TINYMAES_S *maes;
  double sumW; // Sum of the weights
  double val;
  size_t sX = nDim*MAX(nDim,lambda);
  // Compute the total amount of memory to be allocated
  size_t size = sizeof(double)*(nDim*(nDim + 2) + mu + 2*sX);
  maes = (TINYMAES_S *) malloc(size + sizeof(TINYMAES_S));
  // Set size parameters
  maes->nDim = nDim; maes -> mu = mu; maes->lambda = lambda;
  // Set the pointers to the right positions
  maes->M = (double *) (&maes[1]);
  maes->X = &maes->M[nDim*nDim];
  maes->Z = &maes->X[sX];
  maes->X0 = &maes->Z[sX];
  maes->ps = &maes->X0[nDim];
  maes->weights = &maes->ps[nDim];
  // Init the M, X0, X, Z and ps to 0
  memset(maes->M, 0, size);
  // Set M diagonal to 1 (M = Id)
  for(i = 0; i < nDim; i++)
    maes->M[i+nDim*i] = 1;
  switch(weights) {
    case MAES_W_EQUALS:
      for(i = 0; i < mu; i++)
        maes->weights[i] = 1;
      break;
    case MAES_W_LINEAR:
      for(i = 0; i < mu; i++)
        maes->weights[i] = 0.5 + i;
      break;
    case MAES_W_SUPERLINEAR:
      val = log(MAX(mu, lambda/2) + 0.5);
      for(i = 0; i < mu; i++)
        maes->weights[i] = val - log(mu-i);
      break;
    default:
      free(maes);
      return NULL;
  }
  sumW = 0; 
  for(i = 0; i < mu; i++)
    sumW += maes->weights[i];
  for(i = 0; i < mu; i++)
    maes->weights[i] = maes->weights[i]/sumW;
  sumW = 0;
  for(i = 0; i < mu; i++)
    sumW += maes->weights[i]*maes->weights[i];
  maes->mueff = 1/sumW;
  maes->chiN = sqrt(nDim)*(1-1/(4*nDim)+1/(21*nDim*nDim));
  maes->c1 = _alphaCov/((nDim+1.3)*(nDim+1.3) + maes->mueff);
  maes->cs = (maes->mueff + 2)/(maes->mueff + nDim + 5);
  maes->cw = MIN(1 - maes->c1, _alphaCov*(maes->mueff
        + 1/maes->mueff-2)/((nDim+2)*(nDim+2) + _alphaCov*maes->mueff*0.5));
  maes->dSigma = 1 + maes->cs + 2*MAX(0, sqrt((maes->mueff - 1)/(nDim + 1))-1);
  maes->sigma = 1;
  maes->nStep = 0;
  MAES_RANDOM_INIT(&maes->rnd, seed);
  return maes;
}

void TINYMAES_Free(TINYMAES_S *maes) {
  free(maes);
}

void TINYMAES_SetX0(TINYMAES_S *maes, double *x0) {
  int i;
  for(i = 0; i < maes->nDim; i++)
    maes->X0[i] = x0[i];
}

double *TINYMAES_NextStep(TINYMAES_S *maes, int *idx) {
  int i, j, k;
  double s, a, b;
  if(idx) {
    // Compute the new mean X0 = <X>_w = X0 + sigma*<M*Z>_w
    for(j = 0; j < maes->nDim; j++)
      maes->X0[j] = 0;
    for(i = 0; i < maes->mu; i++) {
      k = idx[maes->mu - i - 1];
      s = maes->weights[i];
      for(j = 0; j < maes->nDim; j++)
        maes->X0[j] += s*maes->X[j+k*maes->nDim];
    }
    // X <- Z[idx][:]^t
    for(i = 0; i < maes->mu; i++) {
      k = idx[maes->mu - i - 1];
      for(j = 0; j < maes->nDim; j++)
        maes->X[i+j*maes->mu] = maes->Z[j+k*maes->nDim];
    }
    a = sqrt(maes->mueff*maes->cs*(2-maes->cs));
    b = 1 - maes->cs;
    // Compute the new value of ps using previous value and <Z>_w
    for(j = 0; j < maes->nDim; j++) {
      s = 0;
      for(i = 0; i < maes->mu; i++)
        s += maes->weights[i]*maes->X[i+j*maes->mu];
      maes->ps[j] = b*maes->ps[j] + a*s;
    }
    // Compute the new covariance matrix
    // Determine symmetric matrix DM = Id + c_1/2*(s*s^t - Id)
    //                               + c_w/2*(<z*z^t>_w - Id)
    // Determine only the upper part
    for(i = 0; i < maes->nDim; i++) {
      double d = 0, *vi = &maes->X[i*maes->mu];
      for(k = 0; k < maes->mu; k++)
        d += maes->weights[k]*vi[k]*vi[k];
      maes->Z[i+i*maes->nDim] = 1 + 0.5*maes->cw*(d-1)
                              + 0.5*maes->c1*(maes->ps[i]*maes->ps[i]-1);
      for(j = i+1; j < maes->nDim; j++) {
        double *vj = &maes->X[j*maes->mu];
        d = 0;
        for(k = 0; k < maes->mu; k++)
          d += maes->weights[k]*vi[k]*vj[k];
        maes->Z[j+i*maes->nDim] = 0.5*maes->cw*d + 0.5*maes->c1*maes->ps[i]*maes->ps[j];
      }
    }
    // copy the upper part of DM to the lower part
    for(i = 1; i < maes->nDim; i++)
      for(j = i+1; j < maes->nDim; j++)
        maes->Z[i+j*maes->nDim] = maes->Z[j+i*maes->nDim];
    // Compute D = M*DM
    for(i = 0; i < maes->nDim; i++)
      for(j = 0; j < maes->nDim; j++) {
        double d = 0;
        for(k = 0; k < maes->nDim; k++) // Z symmetrical so we use the most efficient option
          d += maes->M[k+i*maes->nDim]*maes->Z[k+j*maes->nDim];
        maes->X[j+i*maes->nDim] = d;
      }
    // Copy the resulting matrix into M
    memcpy((void*) maes->M, (void*) maes->X, sizeof(double)*maes->nDim*maes->nDim);
    s = 0;
    for(j = 0; j < maes->nDim; j++)
      s += maes->ps[j]*maes->ps[j];
    maes->sigma = maes->sigma*exp(maes->cs/maes->dSigma*(sqrt(s)/maes->chiN-1));
    maes->nStep++;
  }
  MAES_RANDOM_NORMAL(&maes->rnd, maes->Z, maes->nDim*maes->lambda);
  for(i = 0; i < maes->lambda; i++) {
    for(j = 0; j < maes->nDim; j++) {
      s = 0;
      for(k = 0; k < maes->nDim; k++) {
        s += maes->M[k+j*maes->nDim]*maes->Z[k+i*maes->nDim];
      }
      maes->X[j+i*maes->nDim] = maes->sigma*s + maes->X0[j];
    }
  }
  return maes->X;
}


