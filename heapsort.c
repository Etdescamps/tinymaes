#include "heapsort.h"

inline static int _comp(const double *u, const double *v, int nD) {
  int i;
  for(i = 0; i < nD; i++) {
    if(u[i] == v[i])
      continue;
    if(u[i] < v[i])
      return -1;
    return 1;
  }
  return 0;
}

static void _insert_bottom(const double *X, int nD, int id, int *heap, int nHeap) {
  const double *v = &X[id*nD];
  int i, j, k;
  for(i = nHeap; i; i = j) {
    j = i >> 1;
    k = heap[j-1];
    if(_comp(v, &X[k*nD], nD) <= 0)
      break;
    heap[i-1] = k;
  }
  heap[i-1] = id;
}

static void _put_top(const double *X, int nD, int id, int *heap, int nHeap) {
  const double *v = &X[id*nD];
  int i, j, k1, k2, l;
  for(i = 1;;) {
    j = i << 1;
    if(j > nHeap)
      break;
    if(j == nHeap) {
      k1 = heap[j-1] & 0x7fffffff;
      if(_comp(v, &X[k1*nD], nD) < 0) {
        heap[i-1] = k1;
        i = j;
      }
      break;
    }
    k1 = heap[j-1]; k2 = heap[j];
    if(k1 >> 31 == 0 && k2 >> 31 == 0) {
      l = _comp(&X[k1*nD], &X[k2*nD], nD);
      // If they are equals, both max bit are set
      if(l >= 0) {
        k1 = k1 | 0x80000000;
        heap[j-1] = k1;
      }
      if(l <= 0) {
        k2 = k2 | 0x80000000;
        heap[j] = k2;
      }
    }
    if(k2 >> 31) {
      k2 = k2 & 0x7fffffff;
      if(_comp(v, &X[k2*nD], nD) >= 0)
        break;
      heap[i-1] = k2;
      i = j+1;
    }
    else {
      k1 = k1 & 0x7fffffff;
      if(_comp(v, &X[k1*nD], nD) >= 0)
        break;
      heap[i-1] = k1;
      i = j;
    }
  }
  heap[i-1] = id;
}

void heapsort_mu(const double *X, int nD, int lambda, int *idx, int mu) {
  int i, idL;
  const double *vTop;
  // Insert the first mu element into the binary heap tree
  for(i = 0; i < mu; i++) {
    _insert_bottom(X, nD, i, idx, i+1);
  }
  vTop = &X[idx[0]*nD];
  for(i = mu; i < lambda; i++) {
    if(_comp(vTop, &X[i*nD], nD) > 0) {
      _put_top(X, nD, i, idx, mu);
      vTop = &X[idx[0]*nD];
    }
  }
  for(i = mu-1; i > 0; i--) {
    idL = idx[i] & 0x7fffffff;
    idx[i] = idx[0];
    _put_top(X, nD, idL, idx, i);
  }
}


