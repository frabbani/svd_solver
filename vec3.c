#include "defs.h"
#include "vec3.h"

#include <math.h>

void dvec3_zero(dvec3 r) { r[0] = r[1] = r[2] = 0.0; }

void dvec3_copy(dvec3 r, const dvec3 u) {
  r[0] = u[0];
  r[1] = u[1];
  r[2] = u[2];
}

void dvec3_make(dvec3 r, const dvec3 u, const dvec3 v) {
  r[0] = v[0] - u[0];
  r[1] = v[1] - u[1];
  r[2] = v[2] - u[2];
}

void dvec3_muls(dvec3 dv3, double s) {
  dv3[0] *= s;
  dv3[1] *= s;
  dv3[2] *= s;
}

void dvec3_add(dvec3 r, const dvec3 u) {
  r[0] += u[0];
  r[1] += u[1];
  r[2] += u[2];
}

void dvec3_sub(dvec3 r, const dvec3 u) {
  r[0] -= u[0];
  r[1] -= u[1];
  r[2] -= u[2];
}

double dvec3_dot(const dvec3 u, const dvec3 v) {
  return (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}

double dvec3_lensq(const dvec3 dv3) { return (dv3[0] * dv3[0] + dv3[1] * dv3[1] + dv3[2] * dv3[2]); }

void dvec3_cross(const dvec3 u, const dvec3 v, dvec3 n) {
  n[0] = (u[1] * v[2] - u[2] * v[1]);
  n[1] = -(u[0] * v[2] - u[2] * v[0]);
  n[2] = (u[0] * v[1] - u[1] * v[0]);
}

double dvec3_norm(dvec3 r) {
  double mag_sq = dvec3_dot(r, r);
  if (mag_sq > 0.0) {
    double mag = sqrt(mag_sq);
    dvec3_muls(r, 1.0 / mag);
    return mag;
  }
  dvec3_zero(r);
  return 0.0;
}

void dvec3_swap(dvec3 u, dvec3 v) {
  dvec3 t;
  dvec3_copy(t, u);
  dvec3_copy(u, v);
  dvec3_copy(v, t);
}

void dvec3_sort(dvec3 u, dvec3 v, dvec3 w) {
  if (dvec3_dot(u, u) < dvec3_dot(v, v)) {
    dvec3_swap(u, v);
  }
  if (dvec3_dot(v, v) < dvec3_dot(w, w)) {
    dvec3_swap(v, w);
  }
  if (dvec3_dot(u, u) < dvec3_dot(v, v)) {
    dvec3_swap(u, v);
  }
}


void dvec3_outer(const dvec3 u, const dvec3 v, dmat3 m) {
  m[0][0] = u[0] * v[0];
  m[0][1] = u[0] * v[1];
  m[0][2] = u[0] * v[2];
  m[1][0] = u[1] * v[0];
  m[1][1] = u[1] * v[1];
  m[1][2] = u[1] * v[2];
  m[2][0] = u[2] * v[0];
  m[2][1] = u[2] * v[1];
  m[2][2] = u[2] * v[2];
}

void dmat3_add(dmat3 a, const dmat3 b) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      a[i][j] += b[i][j];
    }
  }
}

void dmat3_zero(dmat3_t m) {
  for(int i = 0; i < 3; i++)
  for(int j = 0; j < 3; j++){
    m[i][j] = 0.0;
  }
}

void dmat3_muls(dmat3_t m, double s) {
  for(int i = 0; i < 3; i++)
  for(int j = 0; j < 3; j++){
    m[i][j] *= s;
  }
}
