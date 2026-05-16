(*
    Matrix.ml -- (c) 2022-2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Matrix.ml extends the BiOCamLib `Matrix` module with the distance
    machinery used across the KPop suite.  It defines `Matrix.Base`
    (row normalisations, embeddings, distance-matrix construction,
    rowwise pairwise distances, FAISS-compatible bigarray packing) and
    the `Matrix.t` wrapper that tags each matrix with its role
    (`Twister`, `Twisted`, `Inertia`, `Vectors`, `DMatrix`, ...) along
    with the tabular and binary I/O formats specific to each.

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

open BiOCamLib
open Better

(* Extends BiOCamLib matrix class with distance machinery.
   Encapsulation checks are not performed at this level *)
module Base:
  sig
    include module type of Matrix
    (* Compute row normalisations *)
    val get_normalizations: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                            Space.Distance.t -> Float.Array.t -> t -> Float.Array.t
    (* Get embeddings (principal coordinates) from standard coordinates *)
    val get_embeddings: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                        Space.Distance.t -> Float.Array.t -> t -> t
    val get_distance_matrix: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                             Space.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                              Space.Distance.t -> Float.Array.t -> t -> t -> t
    (* Copy a rowwise [Float.Array.t array] of length n, each row of
       length d, into a preallocated n x d float32 C-layout Bigarray,
       as FAISS expects.  Shared by every code path that hands data
       to Interfaiss.  Caller is responsible for the allocation so
       that the same destination can be reused (e.g. for both index
       and query data) and so that the precise float32 layout choice
       stays at the call site.
       TODO: parallelise when switching to OCaml 5 *)
    val to_bigarray:
      Float.Array.t array ->
      (float, Bigarray.float32_elt, Bigarray.c_layout) Bigarray.Array2.t -> unit
  end
= struct
    let ( .@() ) = Float.Array.( .@() )
    let ( .@()<- ) = Float.Array.( .@()<- )
    include Matrix
    let to_bigarray rows ba =
      Array.iteri
        (fun i row -> Float.Array.iteri (fun j x -> ba.{i, j} <- x) row)
        rows
    (* Compute normalisations for rows *)
    let get_normalizations ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) distance metric m =
      let n_rows = Array.length m.row_names and n_cols = Array.length m.col_names in
      let res = Float.Array.create n_rows
      and rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0 in
      (* Generate points to be computed by the parallel processs *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_rows then
            let to_do = min rows_per_step (n_rows - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          let res = ref [] in
          (* We iterate backwards so as to avoid to have to reverse the list in the end *)
          for i = hi_row downto lo_row do
            Space.Distance.compute_norm distance metric m.data.(i) |> List.accum res
          done;
          lo_row, !res)
        (fun (lo_row, norms) ->
          List.iteri
            (fun offs_i norm_i ->
              res.@(lo_row + offs_i) <- if norm_i = 0. then 1. else norm_i;
              if verbose && !processed_rows mod elements_per_step = 0 then
                Printf.eprintf "%s\r(%s): Done %d/%d rows%!"
                  String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
              incr processed_rows)
            norms)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d rows.\n%!" String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
      res
    (* Compute embeddings (principal coordinates) from standard coordinates *)
    let get_embeddings ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      let d = Float.Array.length metric in
      if Array.length m.col_names <> d then
        Exception.raise_incompatible_geometries __FUNCTION__ (Array.make d "") m.col_names;
      let inv_power =
        match distance with
        | Space.Distance.Euclidean | Cosine | Angle -> 0.5
        | Minkowski p -> 1. /. p in
      let normalized_metric = Float.Array.map (fun x -> x ** inv_power) metric
      and rows_per_step = max 1 (elements_per_step / d) and processed_rows = ref 0
      and n = Array.length m.row_names in
      let data = Array.make n (Float.Array.create 0) in
      (* Generate points to be computed by the parallel process *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n then
            let to_do = min rows_per_step (n - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          let res = ref [] in
          (* We iterate backwards so as to avoid to have to reverse the list in the end *)
          for i = hi_row downto lo_row do
            let data_row = m.data.(i) in
            let v = Float.Array.init d (fun col -> data_row.@(col) *. normalized_metric.@(col)) in
            if normalize then begin
              let norm = Space.Distance.compute_norm distance metric v in
              if norm <> 0. then
                Float.Array.iteri (fun i x -> v.@(i) <- (x /. norm)) v
            end;
            List.accum res v
          done;
          lo_row, !res)
        (fun (lo_row, rows) ->
          List.iteri
            (fun offs_i row_i ->
              data.(lo_row + offs_i) <- row_i;
              if verbose && !processed_rows mod rows_per_step = 0 then
                Printf.eprintf "%s\r(%s): Done %d/%d rows%!"
                  String.TermIO.clear __FUNCTION__ !processed_rows n;
              incr processed_rows)
            rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d rows.\n%!" String.TermIO.clear __FUNCTION__ !processed_rows n;
      { m with data = data }
    (* Compute rowwise distance *)
    let get_distance_matrix ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      let d = Array.length m.row_names in
      (* We compute normalisations *)
      let norms =
        if normalize then
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m
        else
          Float.Array.make d 1. in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let data = Array.init d (fun _ -> Float.Array.create d) in
      (* Generate points to be computed by the parallel processs *)
      let total = (d * (d + 1)) / 2 and i = ref 0 and j = ref 0 and elts_done = ref 0 and end_reached = ref (d = 0) in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file;
          (* We only compute 1/2 of the matrix, and symmetrise it at the end of the computation *)
          let res = ref [] in
          begin try
            let cntr = ref 0 in
            while !cntr < elements_per_step do
              List.accum res (!i, !j);
              incr j;
              if !j > !i then begin
                incr i;
                if !i = d then begin
                  end_reached := true;
                  raise_notrace Exit (* This one is OK as it will be caught *)
                end;
                j := 0
              end;
              incr cntr
            done
          with Exit -> ()
          end;
          List.rev !res)
        (List.map
          (* We decorate each matrix element coordinate with the respective distance *)
          (fun (i, j) ->
            i, j, begin
              Space.Distance.compute
                ~adaptor_a:(fun a -> a /. norms.@(i)) ~adaptor_b:(fun b -> b /. norms.@(j))
                distance metric m.data.(i) m.data.(j)
            end))
        (List.iter
          (fun (i, j, dist) ->
            (* Only here do we actually fill out the memory for the result *)
            data.(i).@(j) <- dist;
            (* We symmetrise the matrix *)
            data.(j).@(i) <- dist;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "%s\r(%s): Done %d/%d elements%!"
                String.TermIO.clear __FUNCTION__ !elts_done total;
            incr elts_done))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements.\n%!" String.TermIO.clear __FUNCTION__ !elts_done total;
      { col_names = m.row_names;
        row_names = m.row_names;
        data = data }
    let get_distance_rowwise ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m1 m2 =
      if m1.col_names <> m2.col_names then
        Exception.raise_incompatible_geometries __FUNCTION__ m1.col_names m2.col_names;
      let r1 = Array.length m1.row_names and r2 = Array.length m2.row_names in
      (* We compute normalisations *)
      let n1, n2 =
        if normalize then
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m1,
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m2
        else
          Float.Array.make r1 1., Float.Array.make r2 1. in
      (*
      to_file ~verbose (transpose ~verbose {
        col_names = m1.row_names;
        row_names = [| "Normalizations" |];
        data = [| n1 |]
      }) "N1.txt";
      to_file ~verbose (transpose ~verbose {
        col_names = m2.row_names;
        row_names = [| "Normalizations" |];
        data = [| n2 |]
      }) "N2.txt";
      *)
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let data = Array.init r2 (fun _ -> Float.Array.create r1) in
      (* Generate points to be computed by the parallel processs *)
      let prod = r1 * r2 in
      let i = ref 0 and j = ref 0 and elts_done = ref 0 and end_reached = ref (prod = 0) in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file;
          let res = ref [] in
          begin try
            let cntr = ref 0 in
            while !cntr < elements_per_step do
              List.accum res (!i, !j);
              incr j;
              if !j = r2 then begin
                incr i;
                if !i = r1 then begin
                  end_reached := true;
                  raise_notrace Exit (* This one is OK as it will be caught *)
                end;
                j := 0
              end;
              incr cntr
            done
          with Exit -> ()
          end;
          List.rev !res)
        (List.map
          (* We decorate each matrix element coordinate with the respective distance *)
          (fun (i, j) ->
            i, j, begin
              Space.Distance.compute
                ~adaptor_a:(fun a -> a /. n1.@(i)) ~adaptor_b:(fun b -> b /. n2.@(j))
                distance metric m1.data.(i) m2.data.(j)
            end))
        (List.iter
          (fun (i, j, dist) ->
            data.(j).@(i) <- dist;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "%s\r(%s): Done %d/%d elements=%.3g%%%!"
                String.TermIO.clear __FUNCTION__
                !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod);
            incr elts_done))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements=%.3g%%.\n%!"
          String.TermIO.clear __FUNCTION__
          !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod);
      { col_names = m1.row_names;
        row_names = m2.row_names;
        data = data }
  end

(* KPop-specialised matrices (encapsulated).
   We include in order not to have a repeated module prefix *)
include (
  struct
    module Type =
      struct
        type t =
          | Distill
          | Twister
          | Inertia
          | Twisted
          | Vectors
          | DMatrix
        let to_string = function
          | Distill -> "KPopDistill"
          | Twister -> "KPopTwister"
          | Inertia -> "KPopInertia"
          | Twisted -> "KPopTwisted"
          | Vectors -> "KPopVectors"
          | DMatrix -> "KPopDMatrix"
        let of_string = function
          | "KPopDistill" -> Distill
          | "KPopTwister" -> Twister
          | "KPopInertia" -> Inertia
          | "KPopTwisted" -> Twisted
          | "KPopVectors" -> Vectors
          | "KPopDMatrix" -> DMatrix
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "archive type" s
      end
    type t = {
      which: Type.t;
      matrix: Base.t
    }
    module Exception =
      struct
        include Base.Exception
        let raise_unexpected_type __FUNCTION__ expected found =
          Exception.raise __FUNCTION__ IO_Format
            (Printf.sprintf "Unexpected matrix types (expected %s, found %s)"
              (Type.to_string expected) (Type.to_string found))
        let raise_unexpected_columns_in_inertia_file __FUNCTION__ found =
          Exception.raise __FUNCTION__ IO_Format begin
            let buf = Buffer.create 1024 in
            Printf.bprintf buf "Unrecognized column names in inertia file (expected [ \"inertia\" ], found [";
            Array.iter (Printf.bprintf buf " \"%s\"") found;
            Printf.bprintf buf " ])";
            Buffer.contents buf
          end
      end
    let empty which = { which; matrix = Base.empty } (* Immutable *)
    let transpose ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m =
      { m with matrix = Base.transpose ~threads ~elements_per_step ~verbose m.matrix }
    let merge_rowwise ?(verbose = false) m1 m2 =
      if m1.which <> m2.which then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Incompatible matrix types (found %s, %s)"
            (Type.to_string m1.which) (Type.to_string m2.which));
      { which = m1.which; matrix = Base.merge_rowwise ~verbose m1.matrix m2.matrix }
    let multiply_matrix_vector_single_threaded ?(verbose = false) m =
      Base.multiply_matrix_vector_single_threaded ~verbose m.matrix
    let multiply_matrix_sparse_vector_single_threaded ?(verbose = false) m =
      Base.multiply_matrix_sparse_vector_single_threaded ~verbose m.matrix
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m v =
      Base.multiply_matrix_vector ~threads ~elements_per_step ~verbose m.matrix v
    let multiply_matrix_matrix ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) which m1 m2 =
      { which; matrix = Base.multiply_matrix_matrix ~threads ~elements_per_step ~verbose m1.matrix m2.matrix }
    (* The following function implements automatic file naming *)
    let make_filename_table which = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ "." ^ Type.to_string which ^ ".txt"
    (* We redefine the implementation for Matrix in order to set the correct KPop types *)
    let of_file ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) which prefix =
      { which; matrix = make_filename_table which prefix |> Base.of_file ~threads ~bytes_per_step ~verbose }
    let to_file ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) m prefix =
      make_filename_table m.which prefix |> Base.to_file ~precision ~threads ~elements_per_step ~verbose m.matrix
    (* *)
    let archive_version = "2022-04-03"
    (* *)
    let to_channel output m =
      Type.to_string m.which |> output_value output;
      archive_version |> output_value output;
      output_value output m.matrix
    let of_channel input =
      let which = Type.of_string (input_value input: string) in
      let version = (input_value input: string) in
      if version <> archive_version then
        Exception.raise_incompatible_archive_version __FUNCTION__ version archive_version;
      { which; matrix = (input_value input: Base.t) }
  end: sig
    module Type:
      sig
        type t =
          | Distill
          | Twister
          | Inertia
          | Twisted
          | Vectors (* Only used internally *)
          | DMatrix (* Only used internally *)
        val to_string: t -> string
        val of_string: string -> t
      end
    type t = {
      which: Type.t;
      matrix: Base.t
    }
    module Exception:
      sig
        include module type of Base.Exception
        val raise_unexpected_type: string -> Type.t -> Type.t -> unit
        val raise_unexpected_columns_in_inertia_file: string -> string array -> unit
      end
    val empty: Type.t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    (* Merge two matrices - the type of the two inputs must be the same *)
    val merge_rowwise: ?verbose:bool -> t -> t -> t
    (* TODO: No type checks are performed (yet) when multiplying matrices *)
    val multiply_matrix_vector_single_threaded: ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_sparse_vector_single_threaded: ?verbose:bool -> t -> Base.sparse_vector_t -> Float.Array.t
    val multiply_matrix_vector: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                                t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Type.t -> t -> t -> t
    (* All file name arguments are in fact _prefixes_ *)
    val of_file: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> Type.t -> string -> t
    (* This one discards type information - use at your own risk *)
    val to_file: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    (* Binary marshalling of the matrix *)
    val to_channel: out_channel -> t -> unit
    val of_channel: in_channel -> t
  end
)

