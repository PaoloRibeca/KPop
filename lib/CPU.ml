(*
    CPU.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    CPU.ml reports the SIMD instruction-set tier available on the running CPU
    (below-AVX2 / AVX2 / AVX-512) and warns when the numeric binaries have been
    pushed onto a slow, sub-AVX2 fallback path -- as happens, for instance,
    inside a virtual machine that advertises a pre-AVX CPU model even on newer
    hardware. The tier mirrors what OpenBLAS's DYNAMIC_ARCH dispatch and the
    FAISS runtime dispatcher (interfaiss_dispatch.c) actually select.

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
*)

open BiOCamLib.Better

(* 0 = below AVX2 (slow fallback), 1 = AVX2, 2 = AVX-512 *)
external isa_tier: unit -> int = "kpop_isa_tier"

let tier_name = function
  | 0 -> "sub-AVX2 (slow fallback)"
  | 1 -> "AVX2"
  | _ -> "AVX-512"

(* Always warn when there is no AVX2 (numeric routines are on a slow fallback);
   under verbose, also report the SIMD tier actually in use *)
let warn_if_slow ?(verbose = false) () =
  match isa_tier () with
  | 0 ->
    Printf.eprintf
      "(%s): WARNING: this CPU/VM has no AVX2, so OpenBLAS and FAISS are running on a slow fallback path (performance will be substantially reduced)\n%!"
      __FUNCTION__
  | t ->
    if verbose then
      Printf.eprintf "(%s): SIMD tier in use: %s\n%!" __FUNCTION__ (tier_name t)

