#pragma once

typedef float float3[3];
typedef float float3x3[3][3];

#define KF3X3(m) ((const float (*)[3])(m))


void f3zero(float3 r);
void f3muls(float3 r, float s);
void f3copy(float3 r, const float3 u);
void f3make(float3 r, const float3 u, const float3 v);
void f3add(float3 r, const float3 u);
void f3sub(float3 r, const float3 u);
float f3dot(const float3 u, const float3 v);
float f3lensq(const float3 r);
float f3norm(float3 r);
void f3cross(const float3 u, const float3 v, float3 n);

void f3x3zero(float3x3 m);
void f3x3ident(float3x3 m);
void f3transf(const float3x3 m, const float3 v, float3 r);