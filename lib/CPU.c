/*
  CPU.c -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

  This file is part of KPop, a scalable method for comparative analysis
  of microbial genomes and environmental samples based on full k-mer
  spectra and correspondence analysis (CA).

  CPU.c exposes to OCaml (via CPU.ml) the SIMD instruction-set tier
  actually available on the running CPU, so the numeric binaries can warn
  when they have been pushed onto a slow, sub-AVX2 fallback path (as
  happens, for instance, inside a virtual machine advertising a pre-AVX
  CPU model).  It uses only compiler built-ins, so it carries no external
  library dependency and is safe to link into every binary.

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

#include <caml/mlvalues.h>

/* 0 = below AVX2 (slow fallback), 1 = AVX2, 2 = AVX-512. This mirrors exactly
   what OpenBLAS's DYNAMIC_ARCH dispatch and interfaiss_dispatch.c select. */
CAMLprim value kpop_isa_tier(value unit) {
  (void)unit;
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx512f"))
    return Val_int(2);
  if (__builtin_cpu_supports("avx2"))
    return Val_int(1);
  return Val_int(0);
}

