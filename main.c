#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "objfile.h"
#include "vec.h"

void compute_covariance(dmat_t* cov, const dvec_t* pts, size_t num_pts) {
  dvec_t mean;
  dvec_t r;
  dmat_t outer;

  dvec_init(&mean, pts[0].n, 0.0);
  dvec_init(&r, pts[0].n, 0.0);
  dmat_init(&outer, pts[0].n, pts[0].n, 0.0);

  for (size_t i = 0; i < num_pts; i++) {
    dvec_add(&mean, &pts[i]);
  }
  dvec_muls(&mean, 1.0 / num_pts);

  // Compute the covariance matrix
  dmat_free(cov);
  dmat_init(cov, mean.n, mean.n, 0.0);

  for (size_t i = 0; i < num_pts; i++) {
    dvec_make(&r, &mean, &pts[i]);
    dvec_outer(&r, &r, &outer);  // r . r^T
    dmat_add(cov, &outer);
  }
  dmat_muls(cov, 1.0 / num_pts);

  dvec_free(&r);
  dmat_free(&outer);
  dvec_free(&mean);
}

typedef struct {
  double lamda;
  int no;
} kvp_t;

int cmp(const void* a, const void* b) {
  kvp_t* kvp_a = (kvp_t*)a;
  kvp_t* kvp_b = (kvp_t*)b;
  if (kvp_a->lamda < kvp_b->lamda)
    return 1;
  if (kvp_a->lamda > kvp_b->lamda)
    return -1;
  return 0;
}

void eigen_sort(const dvec_t* eig_vs, int n, const double* eig_lamdas, dvec_t** eig_vs_sorted,
                double** eig_lamdas_sorted) {
  kvp_t* kvps = malloc(n * sizeof(kvp_t));
  for (int i = 0; i < n; i++) {
    kvps[i].lamda = eig_lamdas[i];
    kvps[i].no = i;
  }
  qsort(kvps, n, sizeof(kvp_t), cmp);

  *eig_lamdas_sorted = malloc(n * sizeof(double));
  *eig_vs_sorted = malloc(n * sizeof(dvec_t));
  for (int i = 0; i < n; i++) {
    (*eig_lamdas_sorted)[i] = kvps[i].lamda;
    dvec_init(&(*eig_vs_sorted)[i], eig_vs[kvps[i].no].n, 0.0);
    dvec_copy(&(*eig_vs_sorted)[i], &eig_vs[kvps[i].no]);
  }

  free(kvps);
}

// void eigen_symmetric_decompose(dmat3 m, dvec3 eig_vecs[3], double eig_lamdas[3]) {
// }

void simple_test() {
  size_t npts = 3;
  dvec_t pts[3];

  // allocate and set points
  dvec_init(&pts[0], 2, 0.0);
  dvec_init(&pts[1], 2, 0.0);
  dvec_init(&pts[2], 2, 0.0);

  pts[0].elems[0] = 1.0;
  pts[0].elems[1] = 2.0;  // (1,2)
  pts[1].elems[0] = 2.0;
  pts[1].elems[1] = 3.0;  // (2,3)
  pts[2].elems[0] = 3.0;
  pts[2].elems[1] = 5.0;  // (3,5)

  // covariance result
  dmat_t cov;
  dmat_init(&cov, 2, 2, 0.0);

  compute_covariance(&cov, pts, npts);

  // print covariance
  printf("COV MAT:\n");
  for (int i = 0; i < cov.m; i++) {
    printf("|");
    for (int j = 0; j < cov.n; j++) {
      printf(" %8.4f", cov.elems[i * cov.n + j]);
    }
    printf(" |\n");
  }

  // cleanup
  for (size_t i = 0; i < npts; i++) {
    dvec_free(&pts[i]);
  }
  dmat_free(&cov);
}

int main(/*int argc, char* argv[]*/) {
  printf("hello world!\n");

  simple_test();

  printf("goodbye!\n");
  return 0;
}