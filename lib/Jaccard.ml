(*
    Jaccard.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Jaccard.ml computes mash-like pairwise distances directly from the
    sparse k-mer count rows of a KMerDB, with no MinHash sketching: the
    binary or weighted (Ruzicka) Jaccard between two samples is read off
    their k-mer count columns, and the mash transform
        d(i, j) = -(1 / k) * ln( 2 J / (1 + J) )
    yields a JC-style approximation of substitution rate.  These
    distances feed the optional branch-length refit (see Refit.ml):
    they are computed only on the O(n K) K-NN edges the refit selects,
    never as a full pairwise matrix.

    Cost per pair is O(m) over the k-mer rows m (the count vectors are
    stored dense per sample, column-major: db.core.data.(sample)).

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

include (
  struct
    (* Jaccard flavour fed to the mash transform.  [Binary]: presence/
       absence Jaccard on the k-mer supports.  [Weighted]: Ruzicka
       (weighted Jaccard) keeping the count multiplicities. *)
    module Kind =
      struct
        type t =
          | Binary
          | Weighted
        let of_string = function
          | "binary" | "jaccard" -> Binary
          | "weighted" | "ruzicka" | "jaccard-weighted" -> Weighted
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "Jaccard kind" s
        let to_string = function
          | Binary -> "binary"
          | Weighted -> "weighted"
      end
    (* Mash transform of a Jaccard value: d = -(1/k) ln(2J/(1+J)).
       J = 0 -> capped at a large distance; J = 1 -> 0. *)
    let mash_of_jaccard ~k j =
      if j <= 0. then 5.
      else
        let x = 2. *. j /. (1. +. j) in
        if x <= 0. then 5. else -. (log x) /. float_of_int k
    (* Build a distance callback  name -> name -> mash distance  closing
       over a KMerDB.  Sample names are resolved to column indices once;
       each call scans the m k-mer rows of the two columns. *)
    let make_dist ?(kind = Kind.Binary) ?(k = 12) (db: KMerDB.t) =
      let m = db.core.n_rows in
      let data = db.core.data in
      let read vec r = KMerDB.CountBAVector.(vec.@(r)) in
      let col name =
        match StringHashtbl.find_opt db.col_names_to_idx name with
        | Some c -> c
        | None ->
          Exception.raise __FUNCTION__ IO_Format
            (Printf.sprintf "Sample '%s' not present in the counts register" name) in
      fun a b ->
        if a = b then 0.
        else begin
          let ca = col a and cb = col b in
          let va = data.(ca) and vb = data.(cb) in
          let j =
            match kind with
            | Kind.Binary ->
              let inter = ref 0 and union = ref 0 in
              for r = 0 to m - 1 do
                let pa = read va r > 0. and pb = read vb r > 0. in
                if pa || pb then begin
                  incr union;
                  if pa && pb then incr inter
                end
              done;
              if !union = 0 then 0.
              else float_of_int !inter /. float_of_int !union
            | Kind.Weighted ->
              let mins = ref 0. and maxs = ref 0. in
              for r = 0 to m - 1 do
                let x = read va r and y = read vb r in
                if x > 0. || y > 0. then begin
                  mins := !mins +. (if x < y then x else y);
                  maxs := !maxs +. (if x > y then x else y)
                end
              done;
              if !maxs <= 0. then 0. else !mins /. !maxs in
          mash_of_jaccard ~k j
        end
  end: sig
    module Kind:
      sig
        type t =
          | Binary
          | Weighted
        val of_string: string -> t
        val to_string: t -> string
      end
    val mash_of_jaccard: k:int -> float -> float
    (* [make_dist ?kind ?k db] returns a [name -> name -> distance]
       callback computing the mash-like distance between two samples of
       the counts register [db].  Suitable as Refit.refit's [~dist]. *)
    val make_dist:
      ?kind:Kind.t -> ?k:int -> KMerDB.t -> (string -> string -> float)
  end
)
