
#include "vec.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void dvec_zero(dvec_t* v) {
  for (int i = 0; i < v->n; ++i) {
    v->elems[i] = 0.0;
  }
}

int dvec_init(dvec_t* v, int n, double init_val) {
  v->n = n;
  v->elems = (double*)malloc(n * sizeof(double));
  if (v->elems) {
    for (int i = 0; i < n; ++i) {
      v->elems[i] = init_val;
    }
    return 0;
  }
  return -1;
}

void dvec_muls(dvec_t* v, double s) {
  for (int i = 0; i < v->n; ++i) {
    v->elems[i] *= s;
  }
}

void dvec_copy(dvec_t* dest, const dvec_t* src) {
  if (dest->n != src->n)
    return;
  for (int i = 0; i < src->n; ++i) {
    dest->elems[i] = src->elems[i];
  }
}

void dvec_add(dvec_t* r, const dvec_t* u) {
  if (r->n != u->n)
    return;
  for (int i = 0; i < r->n; ++i) {
    r->elems[i] += u->elems[i];
  }
}

void dvec_sub(dvec_t* r, const dvec_t* u) {
  if (r->n != u->n)
    return;
  for (int i = 0; i < r->n; ++i) {
    r->elems[i] -= u->elems[i];
  }
}

void dvec_make(dvec_t* r, const dvec_t* u, const dvec_t* v) {
  if (u->n != v->n)
    return;
  for (int i = 0; i < r->n; ++i) {
    r->elems[i] = v->elems[i] - u->elems[i];
  }
}

double dvec_dot(const dvec_t* u, const dvec_t* v) {
  if (u->n != v->n)
    return 0.0;
  double sum = 0.0;
  for (int i = 0; i < u->n; ++i) {
    sum += u->elems[i] * v->elems[i];
  }
  return sum;
}

double dvec_lensq(const dvec_t* v) {
  double sum = 0.0;
  for (int i = 0; i < v->n; ++i) {
    sum += v->elems[i] * v->elems[i];
  }
  return sum;
}

double dvec_norm(dvec_t* v) {
  double mag_sq = dvec_lensq(v);
  if (mag_sq > 0.0) {
    double mag = sqrt(mag_sq);
    dvec_muls(v, 1.0 / mag);
    return mag;
  }
  dvec_zero(v);
  return 0.0;
}

void dvec_free(dvec_t* v) {
  if (v) {
    if (v->elems)
      free(v->elems);
    v->elems = NULL;
    v->n = 0;
  }
}

void dvec_outer(const dvec_t* u, const dvec_t* v, dmat_t* m) {
  // m should be a pointer to a n x n array (row-major)
  for (int i = 0; i < u->n; ++i) {
    for (int j = 0; j < v->n; ++j) {
      m->elems[i * v->n + j] = u->elems[i] * v->elems[j];
    }
  }
}

int dmat_init(dmat_t* mat, int m, int n, double init_val) {
  mat->m = m;
  mat->n = n;
  mat->elems = (double*)malloc(m * n * sizeof(double));
  if (mat->elems) {
    for (int i = 0; i < m * n; ++i) {
      mat->elems[i] = init_val;
    }
    return 0;
  }
  return -1;
}

void dmat_zero(dmat_t* mat) {
  for (int i = 0; i < mat->m * mat->n; ++i) {
    mat->elems[i] = 0.0;
  }
}

void dmat_muls(dmat_t* mat, double s) {
  for (int i = 0; i < mat->m * mat->n; ++i) {
    mat->elems[i] *= s;
  }
}

void dmat_add(dmat_t* mat, const dmat_t* other) {
  if (mat->m != other->m || mat->n != other->n)
    return;
  for (int i = 0; i < mat->m * mat->n; ++i) {
    mat->elems[i] += other->elems[i];
  }
}

void dmat_free(dmat_t *mat) {
  if (mat) {
    if (mat->elems)
      free(mat->elems);
    mat->elems = NULL;
    mat->m = 0;
    mat->n = 0;
  }
}