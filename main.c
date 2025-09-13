#include <math.h>
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
    dvec_point(&r, &mean, &pts[i]);
    dvec_outer(&r, &r, &outer);  // r . r^T
    dmat_add(cov, &outer);
  }
  dmat_muls(cov, 1.0 / num_pts);

  dvec_free(&r);
  dmat_free(&outer);
  dvec_free(&mean);
}

/*
SVD form:
    M = U Σ Vᵗ

    Shapes:
      M: n×d  (n samples, d features)
      U: n×n  (orthonormal basis in sample space)
      Σ: n×d  (diagonal matrix of singular values σ = √λ)
      V: d×d  (orthonormal basis in feature space)

    V: eigenvectors (orthonormal basis, PCA directions) of Mᵗ M (feature space, d×d)
    U: eigenvectors (orthonormal basis) of M Mᵗ (sample space, n×n)
    Σ: diagonal matrix of singular values σ (σ = √λ, not λ directly)

    Vᵗ appears because we are projecting data into the feature PCA basis (world → local).
*/

typedef struct {
  double lamda;
  int no;
} kvp_t;

void compute_pca(const dmat_t* cov, dvec_t** axes) {
  if (cov->m != cov->n)
    return;

  dvec_t r;
  int num_pts = cov->n * 5000;

  dvec_init(&r, cov->n, 0.0);
  dvec_t* elipse_pts = malloc(sizeof(dvec_t) * num_pts);
  for (int i = 0; i < num_pts; i++) {
    dvec_init(&elipse_pts[i], cov->n, 0.0);
    for (int j = 0; j < elipse_pts[i].n; j++) {
      elipse_pts[i].elems[j] = ((double)rand() / RAND_MAX - 0.5) * 2.0;
    }
    dvec_norm(&elipse_pts[i]);
  }

  for (int i = 0; i < num_pts; i++) dmat_transf(cov, &elipse_pts[i]);

  *axes = malloc(sizeof(dvec_t) * cov->n);
  for (int i = 0; i < cov->n; i++) {
    dvec_init(&(*axes)[i], cov->n, 0.0);
  }

  for (int k = 0; k < cov->n; k++) {
    dvec_t* axis = &(*axes)[k];
    double dot_max = 0.0;
    for (int i = 0; i < num_pts; i++) {
      double dot = dvec_dot(&elipse_pts[i], &elipse_pts[i]);
      if (dot > dot_max) {
        dvec_copy(axis, &elipse_pts[i]);
        dot_max = dot;
      }
    }
    dvec_norm(axis);
    for (int i = 0; i < num_pts; i++) {
      dvec_copy(&r, axis);
      double dot = dvec_dot(&elipse_pts[i], axis);
      dvec_muls(&r, dot);
      dvec_sub(&elipse_pts[i], &r);
    }
  }

  for (int i = 0; i < num_pts; i++) {
    dvec_free(&elipse_pts[i]);
  }
  free(elipse_pts);

  dvec_free(&r);
}

void compute_lambdas(const dvec_t* pts, size_t num_pts, const dvec_t* axes, double** lambdas) {
  if (!pts || num_pts == 0 || !axes)
    return;

  int n = pts[0].n;
  *lambdas = malloc(n * sizeof(double));

  dvec_t mean;
  dvec_t r;

  dvec_init(&mean, n, 0.0);
  dvec_init(&r, n, 0.0);

  for (size_t i = 0; i < num_pts; i++) {
    dvec_add(&mean, &pts[i]);
  }
  dvec_muls(&mean, 1.0 / num_pts);

  for (int i = 0; i < n; i++) {
    double* lambda = &(*lambdas)[i];
    *lambda = 0.0;
    for (unsigned int j = 0; j < num_pts; j++) {
      dvec_point(&r, &mean, &pts[j]);
      double dot = dvec_dot(&r, &axes[i]);
      *lambda += (dot * dot);
    }
    *lambda = *lambda / (double)num_pts;
  }

  dvec_free(&mean);
  dvec_free(&r);
}

typedef struct {
  dmat_t V;
  dmat_t Vt;
  dmat_t E;
  dmat_t U;
} svd_t;

void compute_svd(const dvec_t* pts, size_t num_pts, svd_t* svd) {
  dmat_t M = {0};
  dmat_t Mt = {0};
  dmat_t MtM = {0};
  dmat_init(&M, num_pts, pts[0].n, 0.0);
  for (int i = 0; i < num_pts; i++) {
    for (int j = 0; j < M.n; j++) {
      M.elems[i * M.n + j] = pts[i].elems[j];
    }
  }
  dmat_transp(&M, &Mt);
  dmat_mul(&Mt, &M, &MtM);

  dvec_t* eig_axes = NULL;
  compute_pca(&MtM, &eig_axes);

  int num_axes = MtM.n;

  // Rayleigh Quotient
  // M^tM is  in canonical basis. Transform it to feature space
  // compute its projected value along the feature space axes (eigen vectors)
  // (M^T * M) * v_i ​= lambda_i * v_i
  double* eig_lambdas = malloc(num_axes * sizeof(double));
  for (int i = 0; i < num_axes; i++) {
    eig_lambdas[i] = 0.0;
    dvec_t* axis = &eig_axes[i];
    dvec_t v_lambda = {0};
    dvec_copy(&v_lambda, axis);
    // transform eigen vectors by M^tM and find its projected value on eigen vectors (PCA axes)
    dvec_transf(&MtM, &v_lambda);
    eig_lambdas[i] += dvec_dot(axis, &v_lambda);
  }

  // d x d
  for (int i = 0; i < num_axes; i++) dvec_free(&eig_axes[i]);
  free(eig_axes);
  free(eig_lambdas);
  dmat_mul(&Mt, &M, &MtM);
  dmat_free(&Mt);
  dmat_free(&M);
}

/*
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
*/

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

void shape_test() {
  obj_t sphere, shape;

  obj_load(&sphere, "sphere.obj", 0);
  obj_load(&shape, "test.obj", 0);

  printf("point cloud size: %d\n", shape.num_vs);

  dvec_t* pts = NULL;
  int npts = shape.num_vs;
  pts = malloc(sizeof(dvec_t) * npts);

  for (int i = 0; i < npts; i++) {
    dvec_init(&pts[i], 3, 0.0);
    pts[i].elems[0] = shape.vs[i].x;
    pts[i].elems[1] = shape.vs[i].y;
    pts[i].elems[2] = shape.vs[i].z;
  }

  dmat_t cov = {0};
  compute_covariance(&cov, pts, npts);
  printf("COV MAT:\n");
  for (int i = 0; i < cov.m; i++) {
    printf("|");
    for (int j = 0; j < cov.n; j++) {
      printf(" %8.4f", cov.elems[i * cov.n + j]);
    }
    printf(" |\n");
  }

  float3x3 m;
  m[0][0] = cov.elems[0 * cov.n + 0];
  m[0][1] = cov.elems[0 * cov.n + 1];
  m[0][2] = cov.elems[0 * cov.n + 2];
  m[1][0] = cov.elems[1 * cov.n + 0];
  m[1][1] = cov.elems[1 * cov.n + 1];
  m[1][2] = cov.elems[1 * cov.n + 2];
  m[2][0] = cov.elems[2 * cov.n + 0];
  m[2][1] = cov.elems[2 * cov.n + 1];
  m[2][2] = cov.elems[2 * cov.n + 2];

  for (uint32_t i = 0; i < sphere.num_vs; i++) {
    float3 r;
    float3 v = {sphere.vs[i].x, sphere.vs[i].y, sphere.vs[i].z};
    f3transf(KF3X3(m), v, r);
    f3copy(sphere.vs[i].p, r);
  }

  for (int i = 0; i < npts; i++) {
    dvec_free(&pts[i]);
  }
  dmat_free(&cov);

  obj_write(&sphere, "elipsoid.obj");

  free(pts);
  obj_term(&sphere);
  obj_term(&shape);
}

void pca_test() {
  obj_t box, shape;

  obj_load(&box, "box.obj", 0);
  obj_load(&shape, "test.obj", 0);

  printf("point cloud size: %d\n", shape.num_vs);

  dvec_t* pts = NULL;
  int num_pts = shape.num_vs;
  pts = malloc(sizeof(dvec_t) * num_pts);

  for (int i = 0; i < num_pts; i++) {
    dvec_init(&pts[i], 3, 0.0);
    pts[i].elems[0] = shape.vs[i].x;
    pts[i].elems[1] = shape.vs[i].y;
    pts[i].elems[2] = shape.vs[i].z;
  }

  dmat_t cov = {0};
  compute_covariance(&cov, pts, num_pts);
  printf("SHAPE COV MAT:\n");
  for (int i = 0; i < cov.m; i++) {
    printf("|");
    for (int j = 0; j < cov.n; j++) {
      printf(" %8.4f", cov.elems[i * cov.n + j]);
    }
    printf(" |\n");
  }

  dvec_t* axes = NULL;
  compute_pca(&cov, &axes);
  printf("PCA AXES:\n");
  for (int i = 0; i < 3; i++) {
    printf("  <%8.4f, %8.4f, %8.4f>\n", axes[i].elems[0], axes[i].elems[1], axes[i].elems[2]);
  }

  double* lambdas = NULL;
  compute_lambdas(pts, num_pts, axes, &lambdas);
  printf("EIGEN VALUES:\n");
  for (int i = 0; i < 3; i++) {
    printf("  %8.4f\n", lambdas[i]);
  }
  printf("SINGULAR VALUES:\n");
  for (int i = 0; i < 2; i++) {
    printf("  %8.4f\n", sqrt(num_pts * lambdas[i]));
  }

  float3 r, l, u;
  r[0] = axes[0].elems[0];
  r[1] = axes[0].elems[1];
  r[2] = axes[0].elems[2];
  l[0] = axes[1].elems[0];
  l[1] = axes[1].elems[1];
  l[2] = axes[1].elems[2];
  u[0] = axes[2].elems[0];
  u[1] = axes[2].elems[1];
  u[2] = axes[2].elems[2];
  for (unsigned int i = 0; i < box.num_vs; i++) {
    float3 a, b, c;
    f3copy(a, r);
    f3muls(a, box.vs[i].p[0]);
    f3copy(b, l);
    f3muls(b, box.vs[i].p[1]);
    f3copy(c, u);
    f3muls(c, box.vs[i].p[2]);
    f3zero(box.vs[i].p);
    f3add(box.vs[i].p, a);
    f3add(box.vs[i].p, b);
    f3add(box.vs[i].p, c);
  }
  obj_write(&box, "box_transformed.obj");

  dmat_free(&cov);
  for (int i = 0; i < num_pts; i++) {
    dvec_free(&pts[i]);
  }
  free(pts);
  obj_term(&box);
  obj_term(&shape);
}

void sanity_test() {
  //(2,0), (0,1), (4,2), (3,1)
  dvec_t pts[4];
  for (int i = 0; i < 4; i++) {
    dvec_init(&pts[i], 2, 0.0);
  }
  pts[0].elems[0] = 2;
  pts[0].elems[1] = 0;

  pts[1].elems[0] = 0;
  pts[1].elems[1] = 1;

  pts[2].elems[0] = 4;
  pts[2].elems[1] = 2;

  pts[3].elems[0] = 3;
  pts[3].elems[1] = 1;

  dmat_t cov = {0};
  compute_covariance(&cov, pts, 4);
  printf("COV MAT:\n");
  for (int i = 0; i < cov.m; i++) {
    printf("|");
    for (int j = 0; j < cov.n; j++) {
      printf(" %8.4f", cov.elems[i * cov.n + j]);
    }
    printf(" |\n");
  }

  dvec_t* axes;
  compute_pca(&cov, &axes);
  printf("PCA AXES:\n");
  printf("  <%8.4f, %8.4f>\n", axes[0].elems[0], axes[0].elems[1]);
  printf("  <%8.4f, %8.4f>\n", axes[1].elems[0], axes[1].elems[1]);

  double* lambdas = NULL;
  compute_lambdas(pts, 4, axes, &lambdas);
  printf("EIGEN VALUES:\n");
  for (int i = 0; i < 2; i++) {
    printf("  %8.4f\n", lambdas[i]);
  }
  printf("SINGULAR VALUES:\n");
  for (int i = 0; i < 2; i++) {
    printf("  %8.4f\n", sqrt(4 * lambdas[i]));
  }

  for (int i = 0; i < 2; i++) {
    dvec_free(&axes[i]);
  }
  dvec_free(&pts[0]);
  dvec_free(&pts[1]);
  dvec_free(&pts[2]);
  dvec_free(&pts[3]);
  dmat_free(&cov);
}

int main(/*int argc, char* argv[]*/) {
  printf("hello world!\n");

  sanity_test();
  // pca_test();
  // simple_test();
  // shape_test();

  printf("goodbye!\n");
  return 0;
}