(*
    Twister.ml -- (c) 2024-2025 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Twister.ml implements the Twister (CA transformation) database:
    the scaled row coordinates
    `Twister[d,i] = U[i,d] / (sqrt(r[i]) * sv[d])` together with the
    per-dimension inertia.  Exposes the projection routine
    `add_twisted_from_database` that turns a k-mer spectra database
    into Twisted sample coordinates, plus the binary `.KPopTwister`
    file format.

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

(* We cannot open BiOCamLib here due to the ambiguity between BiOCamLib.Matrix and KPop.Matrix *)
module Numbers = BiOCamLib.Numbers
module Processes = BiOCamLib.Processes
module Tools = BiOCamLib.Tools
open BiOCamLib.Better

(* Twister objects are a combination of KPop matrices.
   We include in order not to have a repeated module prefix *)
include (
  struct
    type t = {
      twister: Matrix.t;
      inertia: Matrix.t
    }
    let empty = { twister = Matrix.empty Twister; inertia = Matrix.empty Inertia }
    (* Strictly speaking, we return the _transposed_ of the matrix product here *)
    let add_twisted_from_database
        ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) twister twisted prefix =
      (* Perform some compatibility checks.
         Names of dimensions must be the same for both matrices *)
      let twisted_col_names =
        if twisted.Twisted.twisted.matrix = Matrix.Base.empty then
          twister.twister.matrix.row_names
        else
          twisted.twisted.matrix.col_names in
      if twister.twister.matrix.row_names <> twisted_col_names then
        Matrix.Exception.raise_incompatible_geometries __FUNCTION__
          twister.twister.matrix.row_names twisted_col_names;
      (* We invert the table for k-mer hashes *)
      let num_twister_cols = Array.length twister.twister.matrix.col_names in
      let twister_col_names_to_idx = Hashtbl.create num_twister_cols in
      Array.iteri
        (fun i name ->
          Hashtbl.add twister_col_names_to_idx name i)
        twister.twister.matrix.col_names;
      let db = KMerDB.of_binary ~verbose prefix in
      (* As we process spectra from the database, we have to conform its k-mers
          to the ones in the twister. This can be done once and for all here,
          for all spectra.
         As a bonus, we'll learn the size of the resulting vector *)
      (* As the database has just been loaded there shouldn't be any buffering space,
          but we just double check that this is indeed the case *)
      assert begin
        db.core.n_rows = Array.length db.core.idx_to_row_names &&
        db.core.n_cols = Array.length db.core.idx_to_col_names
      end;
      let db_row_idx_to_twister_col_idx =
        Array.map
          (fun name ->
            match Hashtbl.find_opt twister_col_names_to_idx name with
            | Some idx ->
              idx
            | None ->
              (* We just discard the k-mer *)
              -1)
          db.core.idx_to_row_names in
      (* We make sure there are no duplicate labels *)
      Array.iter
        (fun label ->
          if StringHashtbl.mem db.col_names_to_idx label then
            Exception.raise __FUNCTION__ IO_Format
              (Printf.sprintf "Spectrum '%s' is already present in the destination database" label))
        twisted.twisted.matrix.row_names;
      (*  *)
      let n_cols = db.core.n_cols in
      let cols_per_step = max 1 (elements_per_step / n_cols) and processed_cols = ref 0
      and res = Array.make n_cols (Float.Array.create 0) in
      (* Parallel section *)
      BiOCamLib.Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_cols < n_cols then
            let to_do = min cols_per_step (n_cols - !processed_cols) in
            let new_processed_cols = !processed_cols + to_do in
            let res = !processed_cols, new_processed_cols - 1 in
            processed_cols := new_processed_cols;
            res
          else
            raise End_of_file)
        (fun (lo_col, hi_col) ->
          let res = ref [] in
          for i = lo_col to hi_col do
            let s_v = ref IntMap.empty and acc = ref 0. in
            Array.iteri
              (fun db_row_idx twister_col_idx ->
                if twister_col_idx <> -1 then begin
                  let v = KMerDB.CountBAVector.(db.core.data.(i).@(db_row_idx)) in
                  acc := !acc +. v;
                  s_v := IntMap.add twister_col_idx v !s_v
                end)
              db_row_idx_to_twister_col_idx;
            (* We first normalise and then twist the spectrum *)
            let acc = !acc in
            let s_v = {
              Matrix.Base.length = num_twister_cols;
              elements =
                if acc <> 0. then
                  IntMap.map (fun el -> el /. acc) !s_v
                else
                  !s_v
            } in
            (i, Matrix.multiply_matrix_sparse_vector_single_threaded ~verbose:false twister.twister s_v)
              |> List.accum res
          done;
          !res)
        (List.iter
          (fun (i, row) ->
            (* The transformed column vector becomes a row *)
            res.(i) <- row;
            let new_processed_cols = !processed_cols + 1 in
            if verbose && new_processed_cols / cols_per_step > !processed_cols / cols_per_step then
              Printf.eprintf "%s\r(%s): Twisting spectra: done %d/%d%!"
                String.TermIO.clear __FUNCTION__ new_processed_cols n_cols;
            processed_cols := new_processed_cols))
        threads;
      let row_names = Array.append twisted.twisted.matrix.row_names db.core.idx_to_col_names
      and data = Array.append twisted.twisted.matrix.data res in
      if verbose then
        Printf.eprintf "%s\r(%s): Twisting spectra: done %d/%d.\n%!"
          String.TermIO.clear __FUNCTION__ n_cols n_cols;
      {
        Twisted.inertia = twister.inertia;
        twisted = {
          which = Twisted;
          matrix = { col_names = twisted_col_names; row_names; data }
        }
      }
    (* *)
    let to_files ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) tr prefix =
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.twister prefix;
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.inertia prefix
    let of_files ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) prefix =
      let twister = Matrix.of_file ~threads ~bytes_per_step ~verbose Twister prefix
      and inertia = Matrix.of_file ~threads ~bytes_per_step ~verbose Inertia prefix in
      (* Let's run at least some checks *)
      if inertia.matrix.row_names <> [| "inertia" |] then
        Matrix.Exception.raise_unexpected_columns_in_inertia_file __FUNCTION__ inertia.matrix.row_names;
      if inertia.matrix.col_names <> twister.matrix.row_names then
        Matrix.Exception.raise_incompatible_geometries __FUNCTION__ inertia.matrix.col_names twister.matrix.row_names;
      { twister; inertia }
    (* *)
    let archive_version = "2025-10-08"
    (* *)
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopTwister"
    let to_binary ?(verbose = false) t prefix =
      let path = make_filename_binary prefix in
      let output = open_out path in
      if verbose then
        Printf.eprintf "(%s): Outputting twister to file '%s'...%!" __FUNCTION__ path;
      output_value output "KPopTwister";
      output_value output archive_version;
      Matrix.to_channel output t.twister;
      Matrix.to_channel output t.inertia;
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) prefix =
      let path = make_filename_binary prefix in
      let input = open_in path in
      if verbose then
        Printf.eprintf "(%s): Reading twister from file '%s'...%!" __FUNCTION__ path;
      let which = (input_value input: string) in
      let version = (input_value input: string) in
      if which <> "KPopTwister" || version <> archive_version then
        Matrix.Exception.raise_incompatible_archive_version __FUNCTION__ which version;
      let twister = Matrix.of_channel input in
      let inertia = Matrix.of_channel input in
      close_in input;
      if Matrix.Type.Twister <> twister.which then
        Matrix.Exception.raise_unexpected_type __FUNCTION__ Twister twister.which;
      if Matrix.Type.Inertia <> inertia.which then
        Matrix.Exception.raise_unexpected_type __FUNCTION__ Inertia inertia.which;
      if verbose then
        Printf.eprintf " done.\n%!";
      { twister; inertia }
  end: sig
    type t = {
      (* The matrix describing spectrum-to-coordinates transformation.
        It has type Twister *)
      twister: Matrix.t;
      (* Inertia along each dimension, i.e. the vector of squared singular values.
        The metrics for principal coordinates will be a normalised version of this.
        It is a matrix of type Inertia *)
      inertia: Matrix.t
    }
    val empty: t
    (* This one can fail due to a number of reasons *)
    val add_twisted_from_database: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                                   t -> Twisted.t -> string -> Twisted.t
    (* Input/Output *)
    val to_files: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val of_files: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t
  end
)

