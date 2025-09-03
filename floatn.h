#pragma once

typedef float float3[3];

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
