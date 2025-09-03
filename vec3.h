#pragma once

typedef double dvec3_t[3];
typedef dvec3_t dvec3;

typedef double dmat3_t[3][3];
typedef dmat3_t dmat3;

void dvec3_zero(dvec3 r);
void dvec3_muls(dvec3 r, double s);
void dvec3_copy(dvec3 r, const dvec3 u);
void dvec3_make(dvec3 r, const dvec3 u, const dvec3 v);
void dvec3_add(dvec3 r, const dvec3 u);
void dvec3_sub(dvec3 r, const dvec3 u);
double dvec3_dot(const dvec3 u, const dvec3 v);
double dvec3_lensq(const dvec3 r);
void dvec3_cross(const dvec3 u, const dvec3 v, dvec3 n);
double dvec3_norm(dvec3 r);
void dvec3_outer(const dvec3 u, const dvec3 v, dmat3 m);
void dvec3_swap(dvec3 u, dvec3 v);

void dmat3_add(dmat3 a, const dmat3 b);
void dmat3_zero(dmat3 m);
void dmat3_muls(dmat3 m, double s);
