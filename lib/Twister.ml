(*
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

open BiOCamLib.Better (* We cannot open BiOCamLib here due to the ambiguity with Matrix *)

(* Twister objects are a combination of KPop matrices.
   We include in order not to have a repeated module prefix *)
include (
  struct
    type t = {
      twister: Matrix.t;
      inertia: Matrix.t
    }
    let empty = { twister = Matrix.empty Twister; inertia = Matrix.empty Inertia }
    (* *)
    let to_files ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) tr prefix =
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.twister prefix;
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.inertia prefix
    exception Mismatched_twister_files of string array * string array * string array
    let of_files ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) prefix =
      let twister = Matrix.of_file ~threads ~bytes_per_step ~verbose Twister prefix
      and inertia = Matrix.of_file ~threads ~bytes_per_step ~verbose Inertia prefix in
      (* Let's run at least some checks *)
      if begin
        inertia.matrix.row_names <> [| "inertia" |] ||
        twister.matrix.row_names <> inertia.matrix.col_names
      end then begin
        Printf.eprintf "ERROR: twister.row_names:";
        Array.iter (fun el -> Printf.eprintf "\t\"%s\"" el) twister.matrix.row_names;
        Printf.eprintf "\nERROR: inertia.col_names:";
        Array.iter (fun el -> Printf.eprintf "\t\"%s\"" el) inertia.matrix.col_names;
        Printf.eprintf "\nERROR: inertia.row_names:";
        Array.iter (fun el -> Printf.eprintf "\t\"%s\"" el) inertia.matrix.row_names;
        Printf.eprintf "\n%!";
        Mismatched_twister_files
            (twister.matrix.row_names, inertia.matrix.col_names, inertia.matrix.row_names)
          |> raise
      end;
      { twister; inertia }
    (* Strictly speaking, we return the _transposed_ of the matrix product here *)
    exception Incompatible_twister_and_twisted
    exception Wrong_number_of_columns of int * int * int
    exception Header_expected of string
    exception Float_expected of string
    exception Duplicate_label of string
    let add_twisted_from_files
        ?(normalize = true) ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) ?(debug = false)
        twister twisted fnames =
      ignore elements_per_step; (* Not used at the moment *)
      if twisted.Matrix.which <> Twisted then
        Matrix.Unexpected_type (twisted.Matrix.which, Twisted) |> raise;
      let twisted_col_names =
        if twisted.matrix = Matrix.Base.empty then
          twister.twister.matrix.row_names
        else twisted.matrix.col_names in
      if twister.twister.matrix.row_names <> twisted_col_names then
        raise Incompatible_twister_and_twisted;
      (* We invert the table *)
      let num_twister_cols = Array.length twister.twister.matrix.col_names in
      let twister_col_names_to_idx = Hashtbl.create num_twister_cols in
      Array.iteri
        (fun i name ->
          Hashtbl.add twister_col_names_to_idx name i)
        twister.twister.matrix.col_names;
      (* We decompose the existing twisted matrix *)
      let res = ref StringMap.empty in
      Array.iteri
        (fun i name ->
          res := StringMap.add name twisted.matrix.data.(i) !res)
        twisted.matrix.row_names;
      (* First we read spectra from the files.
         We have to conform the k-mers to the ones in the twister.
         As a bonus, we already know the size of the resulting vector *)
      let fnames = Array.of_list fnames in
      let n = Array.length fnames and file_idx = ref 0 and file = open_in fnames.(0) |> ref and line_num = ref 0
      and labels = ref ("", "") and num_spectra = ref 0 in
      (* Parallel section *)
      BiOCamLib.Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !file_idx = n then
            raise End_of_file
          else begin
            let buf = ref [] in
            (* Each file can contain one or more spectra - we read the next one *)
            begin try
              while true do
                let line_s = input_line !file in
                incr line_num;
                let line = String.Split.on_char_as_array '\t' line_s in
                let l = Array.length line in
                if l <> 2 then
                  Wrong_number_of_columns (!line_num, l, 2) |> raise;
                (* Each file must begin with a header *)
                if !line_num = 1 && line.(0) <> "" then
                  Header_expected line_s |> raise;
                if line.(0) = "" then begin
                  (* New header *)
                  labels := snd !labels, Matrix.Base.strip_external_quotes_and_check line.(1);
                  if !line_num > 1 then begin
                    incr num_spectra;
                    raise Exit
                  end
                end else
                  (* A regular line. The first element is the hash, the second one the count *)
                  List.accum buf (line.(0), line.(1))
              done
            with
            | End_of_file ->
              close_in !file;
              labels := snd !labels, "";
              incr num_spectra;
              if verbose then
                Printf.eprintf "%s\r(%s): [%d/%d] File '%s': Read %d %s on %d %s%s%!"
                  String.TermIO.clear __FUNCTION__ (!file_idx + 1) n fnames.(!file_idx)
                  !num_spectra (String.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                  !line_num (String.pluralize_int "line" !line_num) (if !file_idx + 1 = n then ".\n" else "");
              incr file_idx;
              if !file_idx < n then begin
                file := open_in fnames.(!file_idx);
                line_num := 0
              end
            | Exit ->
              ()
            end;
            if verbose && !file_idx < n then
              Printf.eprintf "%s\r(%s): [%d/%d] File '%s': Read %d %s on %d %s%!"
                String.TermIO.clear __FUNCTION__ (!file_idx + 1) n fnames.(!file_idx)
                !num_spectra (String.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                !line_num (String.pluralize_int "line" !line_num);
            (* The lines are passed in reverse order, but that does not really matter much
                as the final order is determined by the twister *)
            fst !labels, !buf
          end)
        (fun (label, rev_lines) ->
          let t0 = if debug then Sys.time () else 0. in
          let s_v = ref IntMap.empty and acc = ref 0. in
          List.iter
            (fun (name, v) ->
              match Hashtbl.find_opt twister_col_names_to_idx name with
              | Some idx ->
                let v =
                  try
                    float_of_string v
                  with _ ->
                    Float_expected v |> raise in
                acc := !acc +. v;
                s_v := begin
                  match IntMap.find_opt idx !s_v with
                  | Some vv ->
                    (* If there are repeated k-mers, we accumulate them *)
                    IntMap.add idx (vv +. v) !s_v
                  | None ->
                    IntMap.add idx v !s_v
                end
              | None ->
                (* In this case, we just discard the k-mer *)
                ())
            rev_lines;
          let t1 = if debug then Sys.time () else 0. in
          (* We first normalise and then transform the spectrum *)
          let acc = !acc in
          let s_v = {
            Matrix.Base.length = num_twister_cols;
            elements =
              if normalize && acc <> 0. then
                IntMap.map (fun el -> el /. acc) !s_v
              else
                !s_v
          } in
          let t2 = if debug then Sys.time () else 0. in
          let res = Matrix.multiply_matrix_sparse_vector_single_threaded ~verbose:false twister.twister s_v in
          let t3 = if debug then Sys.time () else 0. in
          if debug then
            Printf.eprintf "DEBUG=(lines=%d/%d/%d,%.3g,%.3g,%.3g)\n%!"
              (List.length rev_lines) num_twister_cols (Float.Array.length res) (t1 -. t0) (t2 -. t1) (t3 -. t2);
          label, res)
        (fun (label, row) ->
          (* The transformed column vector becomes a row *)
          match StringMap.find_opt label !res with
          | None ->
            res := StringMap.add label row !res
          | Some _ ->
            Duplicate_label label |> raise)
        threads;
      let n = StringMap.cardinal !res in
      let row_names = Array.make n ""
      and data = Array.make n (Float.Array.create 0) in
      StringMap.iteri
        (fun i label row ->
          row_names.(i) <- label;
          data.(i) <- row)
        !res;
      { Matrix.which = Twisted;
        matrix = { Matrix.Base.col_names = twisted_col_names; row_names; data } }
    (* *)
    let get_metrics_vector m t =
      Space.Distance.Metric.compute m t.inertia.matrix.data.(0)
    let get_metrics_matrix m t = {
      Matrix.which = Metrics;
      matrix = {
        row_names = [| "metrics" |];
        col_names = t.inertia.matrix.col_names;
        data = [| get_metrics_vector m t |]
      }
    }
    (* *)
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopTwister"
    let to_binary ?(verbose = false) t prefix =
      let fname = make_filename_binary prefix in
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      Matrix.to_channel output t.twister;
      Matrix.to_channel output t.inertia;
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) prefix =
      let fname = make_filename_binary prefix in
      let input = open_in fname in
      if verbose then
        Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let twister = Matrix.of_channel input in
      let inertia = Matrix.of_channel input in
      close_in input;
      if Matrix.Type.Twister <> twister.which then
        Matrix.Unexpected_type (Twister, twister.which) |> raise;
      if Matrix.Type.Inertia <> inertia.which then
        Matrix.Unexpected_type (Inertia, inertia.which) |> raise;
      if verbose then
        Printf.eprintf " done.\n%!";
      { twister; inertia }
  end: sig
    type t = {
      twister: Matrix.t; (* Coordinate transformation *)
      inertia: Matrix.t  (* Variance per coordinate *)
    }
    val empty: t
    val to_files: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    exception Mismatched_twister_files of string array * string array * string array
    val of_files: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    (* *)
    exception Incompatible_twister_and_twisted
    exception Wrong_number_of_columns of int * int * int
    exception Header_expected of string
    exception Float_expected of string
    exception Duplicate_label of string
    val add_twisted_from_files:
      ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> ?debug:bool ->
      t -> Matrix.t -> string list -> Matrix.t
    (* *)
    val get_metrics_vector: Space.Distance.Metric.t -> t -> Float.Array.t
    val get_metrics_matrix: Space.Distance.Metric.t -> t -> Matrix.t
    (* *)
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t
  end
)

