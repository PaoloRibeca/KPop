/*
  interfaiss_variant.h -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

  This file is part of KPop, a scalable method for comparative analysis
  of microbial genomes and environmental samples based on full k-mer
  spectra and correspondence analysis (CA).

  interfaiss_variant.h is force-included (via the compiler's -include)
  when interfaiss.cpp is compiled once per ISA variant.  With -DIFSUF=_avx2
  (etc.) it renames the exported C entry points to suffixed symbols, so the
  generic / avx2 / avx512 builds of the FAISS shim can coexist in a single
  statically-linked binary and be selected at run time by interfaiss_dispatch.c.

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

#ifndef INTERFAISS_VARIANT_H
#define INTERFAISS_VARIANT_H
#ifdef IFSUF
#define IFJOIN_(a, b) a##b
#define IFJOIN(a, b) IFJOIN_(a, b)
#define interfaiss_create_flat_index IFJOIN(interfaiss_create_flat_index, IFSUF)
#define interfaiss_create_PQ_index   IFJOIN(interfaiss_create_PQ_index, IFSUF)
#define interfaiss_create_HNSW_index IFJOIN(interfaiss_create_HNSW_index, IFSUF)
#define interfaiss_query_index       IFJOIN(interfaiss_query_index, IFSUF)
#define interfaiss_add_data_to_index IFJOIN(interfaiss_add_data_to_index, IFSUF)
#define interfaiss_train_index       IFJOIN(interfaiss_train_index, IFSUF)
#define interfaiss_free_index        IFJOIN(interfaiss_free_index, IFSUF)
#endif
#endif

