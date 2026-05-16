/*
    RSVD.c -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    RSVD.c implements the randomised truncated singular value
    decomposition (Halko, Martinsson & Tropp 2011) used by `CA.rsvd`.
    Performs subspace iteration with configurable oversampling and
    power-iteration count, computes the QR / SVD pair via LAPACK
    primitives, and returns left singular vectors, singular values,
    and right singular vectors through OCaml-side Bigarrays.

    This program was designed and developed by the author(s),
    with the assistance of the following AI tool(s):
      2026 Claude (Anthropic).
    The final logic and implementation were reviewed and verified in
    their entirety by the author(s).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

/*
 * Randomised truncated SVD for KPop.
 * Algorithm: subspace iteration variant of Halko, Martinsson & Tropp (2011),
 * "Finding structure with randomness", SIAM Review 53(2):217-288.
 *
 * Row-major / LAPACK column-major convention throughout:
 * a row-major r×c matrix X at pointer p is seen by LAPACK (with dimensions
 * swapped to c×r) as X^T. All DGEMM calls are derived from this identity.
 *
 * DGEMM helpers (verified derivations):
 *  dgemm_rm_nn: C[rC×cC] = A[rC×k] * B[k×cC]
 *            C_cm[cC×rC] = B_cm[cC×k] * A_cm[k×rC]
 *   → dgemm_("N","N",&cC,&rC,&k,&one,B,&cC,A,&k,&zero,C,&cC)
 *  dgemm_rm_tn: C[rC×cC] = A^T * B[rA×cC] (A stored as rA×rC)
 *            C_cm[cC×rC] = B_cm[cC×rA] * A_cm^T[rA×rC]
 *   → dgemm_("N","T",&cC,&rC,&rA,&one,B,&cC,A,&rC,&zero,C,&cC)
 *
 * LQ orthogonalisation (thin QR via LAPACK LQ on the transposed view):
 *   Y_rm[big×l] → LAPACK sees Y_ptr as l×big col-major → DGELQF gives
 *   L[l×l]*Q_lq[l×big]; after DORGLQ, Q_lq is in Y_ptr and reading it
 *   back as big×l row-major gives Q_Y (orthonormal columns)
 */

#include <stdlib.h>
#include <string.h>
#define CAML_NAME_SPACE
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

extern void openblas_set_num_threads(int num_threads);

/* BLAS level-3: general matrix multiply */
extern void dgemm_(
  const char *transa, const char *transb,
  const int *m, const int *n, const int *k,
  const double *alpha,
  const double *a, const int *lda,
  const double *b, const int *ldb,
  const double *beta,
  double *c, const int *ldc
);
/* LAPACK: LQ factorisation */
extern void dgelqf_(
  const int *m, const int *n,
  double *a, const int *lda,
  double *tau,
  double *work, int *lwork,
  int *info
);
/* LAPACK: generate Q matrix from LQ factorisation */
extern void dorglq_(
  const int *m, const int *n, const int *k,
  double *a, const int *lda,
  const double *tau,
  double *work, int *lwork,
  int *info
);
/* LAPACK: random number generation (idist=3 → N(0,1)) */
extern void dlarnv_(
  const int *idist,
  int *iseed,
  const int *n,
  double *x
);
/* LAPACK: divide-and-conquer SVD */
extern void dgesdd_(
  const char *jobz,
  const int *m, const int *n,
  double *a, const int *lda,
  double *s,
  double *u, const int *ldu,
  double *vt, const int *ldvt,
  double *work, int *lwork,
  int *iwork,
  int *info
);

/*
 * Row-major DGEMM wrappers
 */

/* C[rC×cC] = A[rC×k] * B[k×cC] */
static void mm_nn(int rC, int cC, int k,
                  const double *A, const double *B, double *C) {
  const double one = 1.0, zero = 0.0;
  dgemm_("N","N", &cC, &rC, &k, &one, B, &cC, A, &k, &zero, C, &cC);
}
/* C[rC×cC] = A^T * B[rA×cC]  (A stored as rA×rC row-major) */
static void mm_tn(int rA, int rC, int cC,
                  const double *A, const double *B, double *C) {
  const double one = 1.0, zero = 0.0;
  dgemm_("N","T", &cC, &rC, &rA, &one, B, &cC, A, &rC, &zero, C, &cC);
}

/*
 * Thin QR in-place via LQ of the transposed LAPACK view.
 * On entry: buf holds a big×l row-major matrix  (big ≥ l).
 * On exit: buf holds Q_Y (big×l row-major, orthonormal columns).
 * tau must be pre-allocated with l elements.
 * Returns LAPACK info (0 = success)
 */
static int lq_orthogonalise(int l, int big,
                            double *buf, double *tau,
                            double *work, int lwork) {
  int info;
  /* LAPACK sees buf as l×big col-major */
  dgelqf_(&l, &big, buf, &l, tau, work, &lwork, &info);
  if (info != 0) return info;
  dorglq_(&l, &big, &l, buf, &l, tau, work, &lwork, &info);
  return info;
}

/*
 * OpenblasRSVD
 *
 * Computes a rank-K approximation of the m×n chi-matrix A (row-major,
 * not overwritten) using subspace iteration:
 *  For iter = 0..q-1:
 *   Y = A · Ω (m×l sketch)
 *   Q_Y = thin_QR(Y) (orthonormal basis for col-space of Y)
 *   Ω = A^T · Q_Y (n×l new sketch)
 *   Ω = thin_QR(Ω) (re-orthogonalise)
 *  .
 *  Y = A · Ω (final sketch)
 *  Q_Y = thin_QR(Y)
 *  B = Q_Y^T · A (l×n small matrix)
 *  [Û, Σ, V̂^T] = SVD(B)
 *  U = Q_Y · Û[:, :K] (m×K)
 *  sv = Σ[:K]
 *  VT = V̂^T[:K, :] (K×n)
 * where l = K + n_oversampling is the sketch dimension.
 *
 * OCaml signature (in CA.ml):
 *  external rsvd :
 *   (float, float64_elt, c_layout) Array1.t -> (* a: chi-matrix m*n *)
 *   int -> int -> (* m, n *)
 *   int -> int -> int -> (* k, n_oversampling, n_power_iter *)
 *   (float, float64_elt, c_layout) Array1.t -> (* u: m*k row-major *)
 *   (float, float64_elt, c_layout) Array1.t -> (* sv: k singular values *)
 *   (float, float64_elt, c_layout) Array1.t -> (* vt: k*n row-major *)
 *   int -> (* threads *)
 *   int
 *   = "OpenblasRSVD_bytecode" "OpenblasRSVD"
 *
 * Returns 0 on success, -1 on allocation failure, LAPACK info otherwise
 */
CAMLprim value OpenblasRSVD(
  value o_a,
  value o_m,
  value o_n,
  value o_k,
  value o_p,
  value o_q,
  value o_u,
  value o_sv,
  value o_vt,
  value o_threads
) {
  CAMLparam5(o_a, o_m, o_n, o_k, o_p);
  CAMLxparam5(o_q, o_u, o_sv, o_vt, o_threads);
  const int m = Int_val(o_m);
  const int n = Int_val(o_n);
  const int k = Int_val(o_k);
  const int p = Int_val(o_p);
  const int q = Int_val(o_q);
  const int threads = Int_val(o_threads);
  openblas_set_num_threads(threads);
  const double *a = (const double *)Caml_ba_data_val(o_a);
  double *u_out = (double *)Caml_ba_data_val(o_u);
  double *sv_out = (double *)Caml_ba_data_val(o_sv);
  double *vt_out = (double *)Caml_ba_data_val(o_vt);
  const int l = k + p;
  /* thin SVD dim of B */
  const int l_svd = (l < n) ? l : n;
  const int k_out = (k < l_svd) ? k : l_svd;
  int info = 0;
  /* Allocate temporaries */
  double *omega = malloc((size_t)n * l * sizeof(double));
  double *Y = malloc((size_t)m * l * sizeof(double));
  double *B = malloc((size_t)l * n * sizeof(double));
  double *U_B = malloc((size_t)l * l_svd * sizeof(double));
  double *sv_B = malloc((size_t)l_svd * sizeof(double));
  double *VT_B = malloc((size_t)l_svd * n * sizeof(double));
  double *tau = malloc((size_t)l * sizeof(double));
  int *iwork = malloc((size_t)8 * l_svd * sizeof(int));
  if (!omega || !Y || !B || !U_B || !sv_B || !VT_B || !tau || !iwork) {
    free(omega); free(Y); free(B); free(U_B);
    free(sv_B); free(VT_B); free(tau); free(iwork);
    caml_failwith("OpenblasRSVD: malloc failed");
  }
  /* Workspace queries */
  double wq;
  int lw_q = -1, lwork_gelqf, lwork_orglq, lwork_gesdd, lwork;
  const int big = (m > n) ? m : n;
  /* dgelqf_ query: M=l, N=big (conservative: covers both Y and omega cases) */
  dgelqf_(&l, &big, omega, &l, tau, &wq, &lw_q, &info);
  lwork_gelqf = (info == 0) ? (int)wq : l * 64;
  /* dorglq_ query: M=l, N=big, K=l */
  lw_q = -1;
  dorglq_(&l, &big, &l, omega, &l, tau, &wq, &lw_q, &info);
  lwork_orglq = (info == 0) ? (int)wq : l * 64;
  /* dgesdd_ query: B passed as n×l col-major (M_lap=n, N_lap=l) */
  lw_q = -1;
  {
    const char jobz = 'S';
    dgesdd_(&jobz, &n, &l, B, &n, sv_B,
            VT_B, &n, U_B, &l_svd,
            &wq, &lw_q, iwork, &info);
  }
  lwork_gesdd = (info == 0) ? (int)wq : n * 64;
  lwork = lwork_gelqf;
  if (lwork_orglq > lwork) lwork = lwork_orglq;
  if (lwork_gesdd > lwork) lwork = lwork_gesdd;
  double *work = malloc((size_t)lwork * sizeof(double));
  if (!work) {
    free(omega); free(Y); free(B); free(U_B);
    free(sv_B); free(VT_B); free(tau); free(iwork);
    caml_failwith("OpenblasRSVD: malloc failed for workspace");
  }
  /* Generate initial random sketch Omega[n×l] ~ N(0,1) */
  {
    /* standard normal */
    int idist = 3;
    /* iseed[3] must be odd */
    int iseed[4] = {3, 1, 4, 5};
    int nl = n * l;
    dlarnv_(&idist, iseed, &nl, omega);
  }
  /* Power iterations */
  for (int qi = 0; qi < q; qi++) {
    /* Y = A * Omega  (m×n × n×l → m×l) */
    mm_nn(m, l, n, a, omega, Y);
    /* QR(Y) → Q_Y in Y (big=m) */
    info = lq_orthogonalise(l, m, Y, tau, work, lwork);
    if (info != 0) goto cleanup;
    /* Omega = A^T * Q_Y (m×n^T × m×l → n×l, A stored as m×n) */
    mm_tn(m, n, l, a, Y, omega);
    /* QR(Omega) → Q_Z in omega (big=n) */
    info = lq_orthogonalise(l, n, omega, tau, work, lwork);
    if (info != 0) goto cleanup;
  }
  /* Final sketch: Y = A * Omega, QR(Y) → Q_Y */
  mm_nn(m, l, n, a, omega, Y);
  info = lq_orthogonalise(l, m, Y, tau, work, lwork);
  if (info != 0) goto cleanup;
  /* B = Q_Y^T * A (l×m × m×n → l×n, Q_Y stored as m×l) */
  mm_tn(m, l, n, Y, a, B);
  /*
   * Thin SVD of B[l×n] via dgesdd_ with transpose trick.
   * Pass B as n×l col-major (M_lap=n, N_lap=l, k_svd=l_svd):
   *  LAPACK U (n×l_svd col-major) → our VT_B (l_svd×n row-major)
   *  LAPACK VT (l_svd×l col-major) → our U_B (l×l_svd row-major)
   */
  {
    const char jobz = 'S';
    dgesdd_(&jobz, &n, &l, B, &n, sv_B,
            VT_B, &n, U_B, &l_svd,
            work, &lwork, iwork, &info);
  }
  if (info != 0) goto cleanup;
  /* Truncate singular values and VT rows */
  memcpy(sv_out, sv_B, (size_t)k_out * sizeof(double));
  for (int d = 0; d < k_out; d++)
    memcpy(vt_out + (size_t)d * n,
           VT_B + (size_t)d * n,
           (size_t)n * sizeof(double));
  /* Build contiguous Û[l×k_out] from first k_out cols of U_B[l×l_svd] */
  {
    double *uhat = malloc((size_t)l * k_out * sizeof(double));
    if (!uhat) { info = -1; goto cleanup; }
    for (int i = 0; i < l; i++)
      memcpy(uhat + (size_t)i * k_out,
             U_B + (size_t)i * l_svd,
             (size_t)k_out * sizeof(double));
    /* U_out = Q_Y * Û (m×l × l×k_out → m×k_out) */
    mm_nn(m, k_out, l, Y, uhat, u_out);
    free(uhat);
  }
cleanup:
  free(omega); free(Y); free(B); free(U_B);
  free(sv_B); free(VT_B); free(tau); free(iwork); free(work);
  CAMLreturn(Val_int(info));
}

/* Bytecode wrapper required for OCaml functions with more than 5 arguments */
CAMLprim value OpenblasRSVD_bytecode(value *argv, int argn) {
  (void)argn;
  return OpenblasRSVD(
    argv[0], argv[1], argv[2], argv[3], argv[4],
    argv[5], argv[6], argv[7], argv[8], argv[9]
  );
}

