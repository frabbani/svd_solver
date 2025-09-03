#include "defs.h"
#include "floatn.h"

#include <math.h>

void f3zero(float3 r) { r[0] = r[1] = r[2] = 0.0f; }

void f3copy(float3 r, const float3 u) {
  r[0] = u[0];
  r[1] = u[1];
  r[2] = u[2];
}

void f3make(float3 r, const float3 u, const float3 v) {
  r[0] = v[0] - u[0];
  r[1] = v[1] - u[1];
  r[2] = v[2] - u[2];
}

void f3muls(float3 r, float s) {
  r[0] *= s;
  r[1] *= s;
  r[2] *= s;
}

void f3add(float3 r, const float3 u) {
  r[0] += u[0];
  r[1] += u[1];
  r[2] += u[2];
}

void f3sub(float3 r, const float3 u) {
  r[0] -= u[0];
  r[1] -= u[1];
  r[2] -= u[2];
}

float f3dot(const float3 u, const float3 v) {
  return (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}

float f3lensq(const float3 dv3) { return (dv3[0] * dv3[0] + dv3[1] * dv3[1] + dv3[2] * dv3[2]); }

void f3cross(const float3 u, const float3 v, float3 n) {
  n[0] = (u[1] * v[2] - u[2] * v[1]);
  n[1] = -(u[0] * v[2] - u[2] * v[0]);
  n[2] = (u[0] * v[1] - u[1] * v[0]);
}

float f3norm(float3 r) {
  float mag_sq = f3dot(r, r);
  if (mag_sq > 0.0f) {
    float mag = sqrtf(mag_sq);
    f3muls(r, 1.0f / mag);
    return mag;
  }
  f3zero(r);
  return 0.0f;
}

