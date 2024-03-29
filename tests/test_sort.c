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

