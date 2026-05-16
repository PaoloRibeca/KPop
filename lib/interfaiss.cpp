/*
  interfaiss.cpp -- (c) 2024      Ünsal Öztürk, <uensal.oeztuerk@gmail.com>
                    (c) 2024-2025 Paolo Ribeca, <paolo.ribeca@gmail.com>

  This file is part of KPop, a scalable method for comparative analysis
  of microbial genomes and environmental samples based on full k-mer
  spectra and correspondence analysis (CA).

  interfaiss.cpp is the C++ wrapper around the FAISS approximate
  nearest-neighbour library exposed to OCaml via `Interfaiss.ml`.  It
  implements the index taxonomy (`flat`, `hnsw(M)`, `pq(...)`) and
  the `create` / `train` / `add` / `query` lifecycle with a common
  C-callable surface over the `float32 Bigarray.Array2` data
  interchange format.

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
#include <faiss/IndexFlat.h>
#include <faiss/IndexPQ.h>
#include <faiss/IndexHNSW.h>
#include <stdexcept>

extern "C" {

index_t* interfaiss_create_flat_index(int d) {
    index_t* idx = new index_t;
    idx->type = INDEX_FLAT;
    idx->index = new faiss::IndexFlatL2(d);
    idx->d = d;
    return idx;
}
index_t* interfaiss_create_PQ_index(int d, int m, int n_bits) {
    index_t* idx = new index_t;
    idx->type = INDEX_PQ;
    // Check that m * nbits <= d to ensure each dimension is quantized
    if (m * n_bits > d) {
      fprintf(stderr, "Error: Invalid parameters for PQ index. Ensure that m * n_bits <= dimension d.\n");
      delete idx;
      return NULL;
    }
    idx->index = new faiss::IndexPQ(d, m, n_bits);
    idx->d = d;
    return idx;
}
index_t* interfaiss_create_HNSW_index(int d, int m, int) {
    index_t* idx = new index_t;
    idx->type = INDEX_HNSW;
    idx->index = new faiss::IndexHNSWFlat(d, m);
    idx->d = d;
    return idx;
}

void interfaiss_query_index(index_t* idx, int d, idx_t n, const dim_t* queries, idx_t* k, float** distances, idx_t** indices) {
    if (idx->d != d) {
        throw std::invalid_argument("Dimension mismatch between index and queries");
    }
    //
    faiss::idx_t k_idx = reinterpret_cast<faiss::Index*>(idx->index)->ntotal;
    *k = *k < k_idx ? *k : k_idx;
    //
    dim_t *flat_distances = (dim_t*)malloc(n*(*k)*sizeof(dim_t));
    idx_t *flat_indices   = (idx_t*)malloc(n*(*k)*sizeof(idx_t));
    //
    reinterpret_cast<faiss::Index*>(idx->index)->search(n, queries, *k, flat_distances, flat_indices);
    //
    *distances = flat_distances;
    *indices   = flat_indices;
    return;
}

void interfaiss_add_data_to_index(index_t* idx, int d, idx_t n, const dim_t* data) {
    if (idx->d != d) {
        throw std::invalid_argument("Dimension mismatch between index and data");
    }
    reinterpret_cast<faiss::Index*>(idx->index)->add(n, data);
    return;
}

void interfaiss_train_index(index_t* idx, int d, idx_t n, const dim_t* data) {
    if (idx->d != d) {
        throw std::invalid_argument("Dimension mismatch between index and queries");
    }
    if (!idx || !data) {
        fprintf(stderr,"Error: Invalid index or data pointer provided for training.\n");
        return;
    }
    //
    switch (idx->type) {
        case INDEX_PQ: {
            // Training a PQ index
            auto pq_index = dynamic_cast<faiss::IndexPQ*>(reinterpret_cast<faiss::Index*>(idx->index));
            if (!pq_index->is_trained) {
                pq_index->train(n, data);
            }
            break;
        }
        case INDEX_FLAT:
        case INDEX_HNSW:
        default:
            //fprintf(stderr,"Training not required or applicable for this index type.\n");
            break;
    }
    return;
}

void interfaiss_free_index(index_t* idx) {
    delete reinterpret_cast<faiss::Index*>(idx->index);
    delete idx;
    return;
}

}

