(*
    Interfaiss.ml -- (c) 2024-2025 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Interfaiss.ml exposes the FAISS approximate nearest-neighbour library
    to OCaml via the `interfaiss.cpp` C++ wrapper.  It defines the
    `Interfaiss.Type.t` index taxonomy (`flat`, `hnsw(M)`, `pq(...)`),
    the `create` / `train` / `add` / `query` lifecycle, and the
    `float32 Bigarray.Array2` data interchange format used by every
    KPop module that needs k-NN search.

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

include (
  struct
    module Type =
      struct
        type t =
          | Flat
          | PQ of int * int
          | HNSW of int
        let to_string = function
          | Flat -> "flat"
          | PQ (m, bits) -> Printf.sprintf "pq(%d,%d)" m bits
          | HNSW m -> "hnsw(" ^ string_of_int m ^ ")"
        let of_string_re = Str.regexp "[(,)]"
        let of_string s =
          let raise kind =
            BiOCamLib.Better.Exception.raise __FUNCTION__ Initialize (Printf.sprintf "%s index '%s'" kind s) in
          match Str.full_split of_string_re s with
          | [ Text "flat" ] ->
            Flat
          | [ Text "pq"; Delim "("; Text m; Delim ","; Text bits; Delim ")" ]
          | [ Text "PQ"; Delim "("; Text m; Delim ","; Text bits; Delim ")" ] ->
            let m, bits =
              try
                int_of_string m, int_of_string bits
              with _ ->
                raise "Unknown" in
            if m < 1 || bits < 1 then
              raise "Invalid";
            PQ (m, bits)
          | [ Text "hnsw"; Delim "("; Text m; Delim ")" ]
          | [ Text "HNSW"; Delim "("; Text m; Delim ")" ] ->
            let m =
              try
                int_of_string m
              with _ ->
                raise "Unknown" in
            if m < 0 then
              raise "Invalid";
            HNSW m
          | _ ->
            raise "Unknown"
      end
    type t
    external create: Type.t -> int -> t = "InterfaissCreate"
    let create ?(index_type = Type.Flat) n =
      try
        create index_type n
      with _ ->
        (* This might happen if the definition of index_t get out of synchro with the C interface *)
        assert false
    type vectors_t = (float, Bigarray.float32_elt, Bigarray.c_layout) Bigarray.Array2.t
    type offsets_t = (int64, Bigarray.int64_elt, Bigarray.c_layout) Bigarray.Array2.t
    type distances_t = (float, Bigarray.float32_elt, Bigarray.c_layout) Bigarray.Array2.t
    external add: t -> vectors_t -> unit = "InterfaissAdd"
    external train: t -> vectors_t -> unit = "InterfaissTrain"
    external query: t -> vectors_t -> int -> offsets_t * distances_t = "InterfaissQuery"
    external delete: t -> unit = "InterfaissDelete"
  end: sig
    module Type:
      sig
        type t = private
          (* We make the type private to implement constraints *)
          | Flat
          | PQ of int * int
          | HNSW of int
        val to_string: t -> string
        val of_string: string -> t (* Can fail due to wrong format *)
      end
    type t
    val create: ?index_type:Type.t -> int -> t
    type vectors_t = (float, Bigarray.float32_elt, Bigarray.c_layout) Bigarray.Array2.t
    type offsets_t = (int64, Bigarray.int64_elt, Bigarray.c_layout) Bigarray.Array2.t
    type distances_t = (float, Bigarray.float32_elt, Bigarray.c_layout) Bigarray.Array2.t
    val add: t -> vectors_t -> unit
    val train: t -> vectors_t -> unit
    val query: t -> vectors_t -> int -> offsets_t * distances_t
    val delete: t -> unit
  end
)

