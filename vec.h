#pragma once

typedef struct {
  int n;
  double* elems;
} dvec_t;

typedef struct {
  int m, n;
  double* elems;
} dmat_t;

int dvec_init(dvec_t* v, int n, double init_val);
void dvec_zero(dvec_t* v);
void dvec_muls(dvec_t* v, double s);
void dvec_copy(dvec_t* dest, const dvec_t* src);
void dvec_add(dvec_t* r, const dvec_t* u);
void dvec_sub(dvec_t* r, const dvec_t* u);
void dvec_point(dvec_t* r, const dvec_t* u, const dvec_t* v);
double dvec_dot(const dvec_t* u, const dvec_t* v);
double dvec_lensq(const dvec_t* v);
double dvec_norm(dvec_t* v);
void dvec_outer(const dvec_t* u, const dvec_t* v, dmat_t* m);
void dvec_free(dvec_t* v);

int dmat_init(dmat_t* mat, int m, int n, double init_val);
void dmat_zero(dmat_t* mat);
void dmat_muls(dmat_t* mat, double s);
void dmat_add(dmat_t* mat, const dmat_t* other);
void dmat_free(dmat_t* mat);
void dmat_transf(const dmat_t* mat, dvec_t* r);
void dmat_transp(const dmat_t* mat, dmat_t* tmat);
void dmat_mul(const dmat_t* mat1, const dmat_t* mat2, dmat_t* mat);