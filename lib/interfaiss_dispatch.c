/*
  interfaiss_dispatch.c -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

  This file is part of KPop, a scalable method for comparative analysis
  of microbial genomes and environmental samples based on full k-mer
  spectra and correspondence analysis (CA).

  interfaiss_dispatch.c provides the plain interfaiss_* C entry points that
  Interfaiss.ml binds to, routing each call to the generic / avx2 / avx512
  build of the FAISS shim selected once, at start-up, from the CPU's actual
  instruction-set support.  Plain function pointers are used (NOT an ifunc
  resolver) so that this works in a statically-linked musl binary.

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

#include "interfaiss.h"

#define VARIANT(suf) \
  index_t* interfaiss_create_flat_index##suf(int); \
  index_t* interfaiss_create_PQ_index##suf(int, int, int); \
  index_t* interfaiss_create_HNSW_index##suf(int, int, int); \
  void interfaiss_query_index##suf(index_t*, int, idx_t, const dim_t*, idx_t*, float**, idx_t**); \
  void interfaiss_add_data_to_index##suf(index_t*, int, idx_t, const dim_t*); \
  void interfaiss_train_index##suf(index_t*, int, idx_t, const dim_t*); \
  void interfaiss_free_index##suf(index_t*);
VARIANT(_generic)
VARIANT(_avx2)
VARIANT(_avx512)

typedef index_t* (*create_flat_t)(int);
typedef index_t* (*create_pq_t)(int, int, int);
typedef index_t* (*create_hnsw_t)(int, int, int);
typedef void (*query_t)(index_t*, int, idx_t, const dim_t*, idx_t*, float**, idx_t**);
typedef void (*add_t)(index_t*, int, idx_t, const dim_t*);
typedef void (*train_t)(index_t*, int, idx_t, const dim_t*);
typedef void (*free_t)(index_t*);

static create_flat_t p_create_flat;
static create_pq_t p_create_pq;
static create_hnsw_t p_create_hnsw;
static query_t p_query;
static add_t p_add;
static train_t p_train;
static free_t p_free;

/* Selected once, before main(), so there is no first-call race. */
__attribute__((constructor))
static void interfaiss_dispatch_init(void) {
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx512f")) {
    p_create_flat = interfaiss_create_flat_index_avx512;
    p_create_pq = interfaiss_create_PQ_index_avx512;
    p_create_hnsw = interfaiss_create_HNSW_index_avx512;
    p_query = interfaiss_query_index_avx512;
    p_add = interfaiss_add_data_to_index_avx512;
    p_train = interfaiss_train_index_avx512;
    p_free = interfaiss_free_index_avx512;
  } else if (__builtin_cpu_supports("avx2")) {
    p_create_flat = interfaiss_create_flat_index_avx2;
    p_create_pq = interfaiss_create_PQ_index_avx2;
    p_create_hnsw = interfaiss_create_HNSW_index_avx2;
    p_query = interfaiss_query_index_avx2;
    p_add = interfaiss_add_data_to_index_avx2;
    p_train = interfaiss_train_index_avx2;
    p_free = interfaiss_free_index_avx2;
  } else {
    p_create_flat = interfaiss_create_flat_index_generic;
    p_create_pq = interfaiss_create_PQ_index_generic;
    p_create_hnsw = interfaiss_create_HNSW_index_generic;
    p_query = interfaiss_query_index_generic;
    p_add = interfaiss_add_data_to_index_generic;
    p_train = interfaiss_train_index_generic;
    p_free = interfaiss_free_index_generic;
  }
}

index_t* interfaiss_create_flat_index(int d) { return p_create_flat(d); }
index_t* interfaiss_create_PQ_index(int d, int m, int n_bits) { return p_create_pq(d, m, n_bits); }
index_t* interfaiss_create_HNSW_index(int d, int m, int ef) { return p_create_hnsw(d, m, ef); }
void interfaiss_query_index(index_t* idx, int d, idx_t n, const dim_t* q, idx_t* k, float** dist, idx_t** ind) {
  p_query(idx, d, n, q, k, dist, ind);
}
void interfaiss_add_data_to_index(index_t* idx, int d, idx_t n, const dim_t* data) { p_add(idx, d, n, data); }
void interfaiss_train_index(index_t* idx, int d, idx_t n, const dim_t* data) { p_train(idx, d, n, data); }
void interfaiss_free_index(index_t* idx) { p_free(idx); }

