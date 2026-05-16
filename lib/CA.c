/*
  CA.c -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

  This file is part of KPop, a scalable method for comparative analysis
  of microbial genomes and environmental samples based on full k-mer
  spectra and correspondence analysis (CA).

  CA.c provides the C bindings to LAPACK's `dgesdd_` (thin SVD)
  consumed by `CA.twist`. Wraps the FORTRAN call into an OCaml-callable
  stub that accepts the chi matrix as a row-major Bigarray and writes
  left singular vectors, singular values, and right singular vectors
  back into caller-provided Bigarrays.

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

#include <stdlib.h>
#define CAML_NAME_SPACE
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

/* OpenBLAS thread-count control */
extern void openblas_set_num_threads(int num_threads);

/* LAPACK divide-and-conquer SVD (Fortran calling convention) */
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
 * OpenblasSVD
 *
 * Computes the thin SVD of an m x n matrix stored row-major as a flat
 * float64 Bigarray. The input array is overwritten during computation.
 *
 * We exploit the transpose trick: LAPACK is column-major (Fortran order),
 * so passing our row-major m x n data as column-major makes LAPACK see
 * the n x m transposed matrix. Swapping M/N and the U/VT output buffers
 * gives us exactly what we need, because:
 *  LAPACK's U output (n x k, column-major) == our Vt (k x n, row-major)
 *  LAPACK's VT output (k x m, column-major) == our U (m x k, row-major)
 * where k = min(m, n).
 *
 * After the call:
 *        u[i * k + d] = left singular vector entry: k-mer i, dimension d
 *                s[d] = d-th singular value (non-increasing order)
 *  vt[d * n_cols + j] = right singular vector entry: dimension d, sample j
 *
 * Returns the LAPACK info code (0 = success, <0 = bad argument, >0 = no convergence).
 *
 * OCaml signature (declared in CA.ml):
 *  external dgesdd :
 *   (float, float64_elt, c_layout) Bigarray.Array1.t -> (* a: input m*n, overwritten *)
 *   int -> (* m: rows = k-mers *)
 *   int -> (* n: cols = samples *)
 *   (float, float64_elt, c_layout) Bigarray.Array1.t -> (* u: output m*k row-major *)
 *   (float, float64_elt, c_layout) Bigarray.Array1.t -> (* s: output k singular values *)
 *   (float, float64_elt, c_layout) Bigarray.Array1.t -> (* vt: output k*n row-major *)
 *   int -> (* threads: OpenBLAS thread count *)
 *   int
 *   = "OpenblasSVD_bytecode" "OpenblasSVD"
 */
CAMLprim value OpenblasSVD(
  value o_a,
  value o_m,
  value o_n,
  value o_u,
  value o_s,
  value o_vt,
  value o_threads
) {
  CAMLparam5(o_a, o_m, o_n, o_u, o_s);
  CAMLxparam2(o_vt, o_threads);
  /* rows = k-mers */
  int m = Int_val(o_m);
  /* cols = samples */
  int n = Int_val(o_n);
  int k = (m < n) ? m : n;
  openblas_set_num_threads(Int_val(o_threads));
  double *a = (double *)Caml_ba_data_val(o_a);
  double *u = (double *)Caml_ba_data_val(o_u);
  double *s = (double *)Caml_ba_data_val(o_s);
  double *vt = (double *)Caml_ba_data_val(o_vt);
  const char jobz = 'S';
  int info = 0;
  /* iwork must be allocated before the workspace query */
  int *iwork = (int *)malloc((size_t)(8 * k) * sizeof(int));
  if (!iwork)
    caml_failwith("OpenblasSVD: malloc failed for iwork");
  /*
   * Workspace query: pass lwork = -1 to obtain the optimal workspace size.
   * After the query, work[0] holds the required size as a double.
   *
   * We call with M_lapack = n (samples), N_lapack = m (k-mers):
   *  - LAPACK sees our row-major m x n data as column-major n x m
   *  - U_lapack (n x k, column-major) receives our Vt (k x n, row-major)
   *  - VT_lapack (k x m, column-major) receives our U (m x k, row-major)
   */
  double work_query = 0.0;
  int lwork = -1;
  dgesdd_(&jobz, &n, &m, a, &n, s, vt, &n, u, &k,
          &work_query, &lwork, iwork, &info);
  if (info != 0) {
    free(iwork);
    caml_failwith("OpenblasSVD: workspace query failed");
  }
  lwork = (int)work_query;
  double *work = (double *)malloc((size_t)lwork * sizeof(double));
  if (!work) {
    free(iwork);
    caml_failwith("OpenblasSVD: malloc failed for work array");
  }
  /* Actual SVD computation */
  dgesdd_(&jobz, &n, &m, a, &n, s, vt, &n, u, &k,
          work, &lwork, iwork, &info);
  free(work);
  free(iwork);
  CAMLreturn(Val_int(info));
}

/* Bytecode wrapper required for OCaml functions with more than 5 arguments */
CAMLprim value OpenblasSVD_bytecode(value *argv, int argn) {
  (void)argn;
  return OpenblasSVD(
    argv[0], argv[1], argv[2],
    argv[3], argv[4], argv[5],
    argv[6]
  );
}

