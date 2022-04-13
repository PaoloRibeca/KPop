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

open BiOCamLib

module IntSet = Tools.IntSet
module StringSet = Tools.StringSet
module StringMap = Tools.StringMap

(* For counts. We assume each count is < 2^31 *)
module IBAVector =
Tools.BA.Vector (
  struct
    include Int32
    type elt_t = Bigarray.int32_elt
    let elt = Bigarray.Int32
  end
)

(* For normalisations and combined vectors. We assume normalisations might be > 2^31 *)
module FBAVector =
Tools.BA.Vector (
  struct
    include Float
    let of_float x = x
    let to_float x = x
    type elt_t = Bigarray.float64_elt
    let elt = Bigarray.Float64
  end
)

module [@warning "-32"] rec Transformation:
  sig
    (* Transformation function *)
    type t =
      | None
      | Normalize of int * float
      | CLR of int * float
      | Pseudo of int * float
    exception Invalid_transformation of string * int * float
    val compute: which:t -> col_num:int -> col_stats:Statistics.t -> row_num:int -> row_stats:Statistics.t -> int -> float
    type parameters_t = {
      which: string;
      threshold: int;
      power: float
    }
    exception Unknown_transformation of string
    val of_parameters: parameters_t -> t
    val to_parameters: t -> parameters_t
  end
= struct
    type t =
      | None
      | Normalize of int * float
      | CLR of int * float
      | Pseudo of int * float
    exception Invalid_transformation of string * int * float
    let epsilon = 0.1
    let [@warning "-27"] compute ~which ~col_num ~col_stats ~row_num ~row_stats counts =
      match which with
      | None ->
        float_of_int counts
      | Normalize (threshold, power) ->
        if counts >= threshold then
          (float_of_int counts ** power) /. col_stats.Statistics.sum
        else
          0.
      | CLR (threshold, power) ->
        let v =
          if counts >= threshold then
            float_of_int counts
          else
            0. in
        let v = max v epsilon in
        (log v *. power) -. (col_stats.Statistics.sum_log /. float_of_int col_stats.non_zero)
      | Pseudo (threshold, power) ->
        if power < 0. then
          Invalid_transformation ("pseudocounts", threshold, power) |> raise;
        let counts = float_of_int counts in
        let v =
          if power = 0. then
            (float_of_int col_stats.Statistics.max) *. log ((counts +. 1.) /. float_of_int threshold)
          else begin
            let red_threshold = threshold - 1 |> float_of_int |> max 0. in
            let c_p = red_threshold ** power in
            if power < 1. then
              ((counts ** power) -. c_p) *. ((float_of_int col_stats.Statistics.max) ** (1. -. power)) /. power
            else
              ((counts ** power) -. c_p) /. ((float_of_int threshold ** power) -. c_p)
          end in
        floor v /. col_stats.sum |> max 0.
    type parameters_t = {
      which: string;
      threshold: int;
      power: float
    }
    exception Unknown_transformation of string
    let of_parameters { which; threshold; power } =
      match which with
      | "none" ->
        None
      | "normalize" | "normalise" | "norm" ->
        Normalize (threshold, power)
      | "clr" | "CLR" ->
        CLR (threshold, power)
      | "pseudocounts" | "pseudo" ->
        Pseudo (threshold, power)
      | s ->
        Unknown_transformation s |> raise
    let to_parameters = function
      | None -> { which = "none"; threshold = 1; power = 1. }
      | Normalize (threshold, power) -> { which = "normalize"; threshold; power }
      | CLR (threshold, power) -> { which = "clr"; threshold; power }
      | Pseudo (threshold, power) -> { which = "pseudocounts"; threshold; power }
  end
and [@warning "-32"] Statistics:
  sig
    type t = {
      non_zero: int;
      min: int;
      max: int;
      sum: float;
      sum_log: float
      (*
      mean: float;
      variance: float
      (* Could be extended with more moments if needed *)
      *)
      }
    val empty: t
    type table_t = {
      col_stats: t array;
      row_stats: t array
    }
    val table_of_db: ?threads:int -> ?verbose:bool -> Transformation.t -> KMerDB.t -> table_t
  end
= struct
    type t = {
      non_zero: int;
      min: int;
      max: int;
      sum: float;
      sum_log: float
      (*
      mean: float;
      variance: float
      *)
    }
    let empty = {
      non_zero = 0;
      min = 0;
      max = 0;
      sum = 0.;
      sum_log = 0.
      (*
      mean = 0.;
      variance = 0.
      *)
    }
    type table_t = {
      col_stats: t array;
      row_stats: t array
    }
    (*
    let resize_statistics_array ?(exact = false) n a =
      resize_array ~exact n empty_statistics a
    *)
    type col_or_row_t =
    | Col
    | Row
    let col_or_row_to_string = function
    | Col -> "column"
    | Row -> "row"
    let table_of_db ?(threads = 1) ?(verbose = false) transf db =
      let { Transformation.threshold; power; _ } = Transformation.to_parameters transf in
      let core = db.KMerDB.core in
      let compute_one what n =
        let non_zero = ref 0 and min = ref 0 and max = ref 0 and sum = ref 0. and sum_log = ref 0.
        and red_len =
          match what with
          | Col ->
            core.n_rows - 1
          | Row ->
            core.n_cols - 1 in
        for i = 0 to red_len do
          let v =
            match what with
            | Col ->
              IBAVector.N.to_int core.storage.(n).IBAVector.=(i)
            | Row ->
              IBAVector.N.to_int core.storage.(i).IBAVector.=(n) in
          let v =
            if v >= threshold then
              v
            else
              0 in
          if v <> 0 then begin
            incr non_zero;
            let f_v = float_of_int v in
            min := Stdlib.min !min v;
            max := Stdlib.max !max v;
            sum := !sum +. (f_v ** power);
            sum_log := !sum_log +. (log f_v *. power)
          end
        done;
        { non_zero = !non_zero;
          min = !min;
          max = !max;
          sum = !sum;
          sum_log = !sum_log } in
      let compute_all what =
        let n =
          match what with
          | Col ->
            core.n_cols
          | Row ->
            core.n_rows in
        let step = n / threads / 5 |> max 1 and processed = ref 0 and res = ref [] in
        Tools.Parallel.process_stream_chunkwise
          (fun () ->
            if verbose then
              Printf.eprintf "\rComputing %s statistics [%d/%d]%!" (col_or_row_to_string what) !processed n;
            let to_do = n - !processed in
            if to_do > 0 then begin
              let to_do = min to_do step in
              let res = !processed, to_do in
              processed := !processed + to_do;
              res
            end else
              raise End_of_file)
          (fun (processed, to_do) ->
            let res = ref [] and red_to_do = to_do - 1 in
            for i = 0 to red_to_do do
              processed + i |> compute_one what |> Tools.List.accum res
            done;
            Tools.Array.of_rlist !res)
          (Tools.List.accum res)
          threads;
        let rec binary_merge_arrays processed to_do =
          match processed, to_do with
          | [], [] ->
            [||]
          | [], [a] ->
            a
          | [res], [] ->
            res
          | [res], [a] ->
            Array.append res a
          | _, [] ->
            List.rev processed |> binary_merge_arrays []
          | _, [a] ->
            a :: processed |> List.rev |> binary_merge_arrays []
          | _, a1 :: a2 :: tl ->
            binary_merge_arrays ((Array.append a1 a2) :: processed) tl in
        let res = List.rev !res |> binary_merge_arrays [] in
        if verbose then
          Printf.eprintf "\rComputing %s statistics [%d/%d]\n%!" (col_or_row_to_string what) n n;
        res in
      { col_stats = compute_all Col;
        row_stats = compute_all Row }
  end
and KMerDB:
(*
  K-mers are stored as columns, labels as rows
*)
  sig
    type marshalled_t = {
      n_cols: int; (* The number of vectors *)
      n_rows: int; (* The number of k-mers *)
      n_meta: int; (* The number of metadata fields *)
      (* We number rows, columns and metadata fields starting from 0 *)
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      idx_to_meta_names: string array;
      (* *)
      meta: string array array; (* Dims = n_cols * n_meta *)
      storage: IBAVector.t array (* Frequencies are stored as integers. Dims = n_cols * n_rows *)
    }
    type t = {
      core: marshalled_t;
      (* Inverted hashes for parsing *)
      col_names_to_idx: (string, int) Hashtbl.t; (* Labels *)
      row_names_to_idx: (string, int) Hashtbl.t; (* Hashes *)
      meta_names_to_idx: (string, int) Hashtbl.t (* Metadata fields *)
    }
    val make_empty: unit -> t
    (* Adds metadata - the first field must be the label *)
    exception Wrong_number_of_columns of int * int * int
    val add_meta: ?verbose:bool -> t -> string -> t
    (* Adds text files containing k-mers. The first line must contain the label.
       Multiple files separated by "\t\n" can be chained in the same input *)
    exception Header_expected of string
    exception Wrong_format of int * string
    (* Add files contaning k-mers. Multiple files can be chained *)
    val add_files: ?verbose:bool -> t -> string list -> t
    (* Select column names identified by regexps on metadata fields *)
    val selected_from_regexps: ?verbose:bool -> t -> (string * Str.regexp) list -> StringSet.t
    val selected_negate: t -> StringSet.t -> StringSet.t
    (* Add a linear combination of the vectors having the given labels,
        and name it as directed *)
    val add_sum_selected: ?threads:int -> ?verbose:bool -> t -> string -> StringSet.t -> t
    (* Remove the vectors with the given labels *)
    val remove_selected: t -> StringSet.t -> t
    (* Output information about the contents *)
    val output_summary: ?verbose:bool -> t -> unit
    (* Binary marshalling of the database *)
    exception Incompatible_archive_version of string * string
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t
    val make_filename_binary: string -> string
    module TableFilter:
      sig
        type t = {
          print_row_names: bool;
          print_col_names: bool;
          print_metadata: bool;
          transpose: bool;
          transform: Transformation.t;
          print_zero_rows: bool;
          filter_columns: StringSet.t;
          precision: int
        }
        val default: t
      end
    (* Readable output *)
    val to_table:
      ?filter:TableFilter.t -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val make_filename_table: string -> string

  end
= struct
    type marshalled_t = {
      n_cols: int;
      n_rows: int;
      n_meta: int;
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      idx_to_meta_names: string array;
      meta: string array array;
      storage: IBAVector.t array
    }
    type t = {
      core: marshalled_t;
      col_names_to_idx: (string, int) Hashtbl.t;
      row_names_to_idx: (string, int) Hashtbl.t;
      meta_names_to_idx: (string, int) Hashtbl.t
    }
    (* *)
    let make_empty () = {
      core = {
        n_cols = 0;
        n_rows = 0;
        n_meta = 0;
        idx_to_col_names = [||];
        idx_to_row_names = [||];
        idx_to_meta_names = [||];
        meta = [||];
        storage = [||]
      };
      col_names_to_idx = Hashtbl.create 16;
      row_names_to_idx = Hashtbl.create 16;
      meta_names_to_idx = Hashtbl.create 16
    }
    let output_summary ?(verbose = false) db =
      Printf.eprintf "[Vector labels (%d)]:" (Array.length db.core.idx_to_col_names);
      Array.iter (Printf.eprintf " '%s'") db.core.idx_to_col_names;
      Printf.eprintf "\n%!";
      if verbose then begin
        Printf.eprintf "[K-mer hashes (%d)]:" (Array.length db.core.idx_to_row_names);
        Array.iter (Printf.eprintf " '%s'") db.core.idx_to_row_names;
        Printf.eprintf "\n%!"
      end;
      Printf.eprintf "[Meta-data fields (%d)]:" (Array.length db.core.idx_to_meta_names);
      Array.iter (Printf.eprintf " '%s'") db.core.idx_to_meta_names;
      Printf.eprintf "\n%!"
    (* *)
    let resize_array ?(exact = false) n null a =
      let l = Array.length a in
      if n > l then begin
        let res =
          Array.make begin
            if exact then
              n
            else
              max n (l * 14 / 10)
          end null in
        Array.blit a 0 res 0 l;
        res
      end else if n < l && exact then
        (* We have to resize in order to honour ~exact *)
        Array.sub a 0 n
      else
        a
    let resize_string_array ?(exact = false) n a =
      resize_array ~exact n "" a
    (* *)
    let _resize_t_array_ ?(exact = false) length resize create_null nx ny a =
      let lx = Array.length a in
      let eff_nx =
        if exact then
          nx
        else if lx < nx then
          max nx (lx * 14 / 10)
        else
          lx
      and eff_ny =
        if exact then
          ny
        else if lx > 0 then begin
          (* We assume all bigarrays to have the same size *)
          let ly = length a.(0) in
          if ly < ny then
            max ny (ly * 14 / 10)
          else
            ly
        end else
          ny in
      (*Printf.eprintf "Resizing to (%d,%d) - asked (%d,%d)...\n%!" eff_nx eff_ny nx ny;*)
      if eff_nx > lx then
        Array.append
          (Array.map (resize ?exact:(Some true) eff_ny) a)
          (Array.init (eff_nx - lx) (fun _ -> create_null eff_ny))
      else if eff_nx < lx then
        Array.map (resize ?exact:(Some true) eff_ny) (Array.sub a 0 eff_nx)
      else (* eff_nx = lx *)
        if lx > 0 && eff_ny = length a.(0) then
          a
        else
          Array.map (resize ?exact:(Some true) eff_ny) a
    let resize_string_array_array ?(exact = false) =
      _resize_t_array_ ~exact Array.length resize_string_array (fun l -> Array.make l "")
    module BAVectorMisc (M: Tools.BA.VectorType) =
      struct
        let resize ?(exact = false) n a =
          let l = M.length a in
          if n > l then begin
            let res =
              M.init begin
                if exact then
                  n
                else
                  max n (l * 14 / 10)
              end M.N.zero in
            M.blit a (M.sub res 0 l);
            res
          end else if n < l && exact then
            (* We have to resize in order to honour ~exact *)
            let res = M.init n M.N.zero in
            M.blit (M.sub a 0 n) res;
            res
          else
            a
        let resize_array ?(exact = false) =
          _resize_t_array_ ~exact M.length resize (fun l -> M.init l M.N.zero)
        (*
        let to_float_array a =
          let l = M.length a in
          let red_l = l - 1
          and res = Float.Array.make l 0. in
          for i = 0 to red_l do
            M.N.to_float a.M.=(i) |> Float.Array.set res i
          done;
          res
        *) 
      end
    module IBAVectorMisc = BAVectorMisc (IBAVector)
    module FBAVectorMisc = BAVectorMisc (FBAVector)
    (* Utility functions *)
    let invert_table a =
      let res = Hashtbl.create (Array.length a) in
      Array.iteri (fun i name -> Hashtbl.add res name i) a;
      res
    let add_empty_column_if_needed db label =
      let n_cols = !db.core.n_cols in
      let aug_n_cols = n_cols + 1 in
      if Hashtbl.mem !db.col_names_to_idx label |> not then begin
        Hashtbl.add !db.col_names_to_idx label n_cols; (* THIS ONE CHANGES !db *)
        db :=
          { !db with
            core =
            { !db.core with
              n_cols = aug_n_cols;
              (* We have to resize all the relevant containers *)
              idx_to_col_names = Array.append !db.core.idx_to_col_names [| label |];
              meta = resize_string_array_array aug_n_cols !db.core.n_meta !db.core.meta;
              storage = IBAVectorMisc.resize_array aug_n_cols !db.core.n_rows !db.core.storage } }
      end
    (* *)
    let archive_version = "2022-04-03"
    (* *)
    exception Incompatible_archive_version of string * string
    let to_binary ?(verbose = false) db fname =
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      output_value output "KPopCounter";
      output_value output archive_version;
      output_value output {
        db.core with
        (* We have to truncate all the containers *)
        idx_to_col_names = resize_string_array ~exact:true db.core.n_cols db.core.idx_to_col_names;
        idx_to_row_names = resize_string_array ~exact:true db.core.n_rows db.core.idx_to_row_names;
        idx_to_meta_names = resize_string_array ~exact:true db.core.n_meta db.core.idx_to_meta_names;
        meta = resize_string_array_array ~exact:true db.core.n_cols db.core.n_meta db.core.meta;
        storage = IBAVectorMisc.resize_array ~exact:true db.core.n_cols db.core.n_rows db.core.storage
      };
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) fname =
      let input = open_in fname in
      if verbose then
        Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let which = (input_value input: string) in
      let version = (input_value input: string) in
      if which <> "KPopCounter" || version <> archive_version then
        Incompatible_archive_version (which, version) |> raise;
      let core = (input_value input: marshalled_t) in
      close_in input;
      if verbose then
        Printf.eprintf " done.\n%!";
      { core = core;
        col_names_to_idx = invert_table core.idx_to_col_names;
        row_names_to_idx = invert_table core.idx_to_row_names;
        meta_names_to_idx = invert_table core.idx_to_meta_names }
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopCounter"
    (* *)
    exception Wrong_number_of_columns of int * int * int
    let add_meta ?(verbose = false) db fname =
      let input = open_in fname and line_num = ref 0 in
      let header = input_line input |> Tools.Split.on_char_as_array '\t' in
      incr line_num;
      (* We add the names *)
      let missing = ref [] in
      Array.iteri
        (fun i name ->
          if i > 0 && Hashtbl.mem db.meta_names_to_idx name |> not then
            Tools.List.accum missing name)
        header;
      let missing = Tools.Array.of_rlist !missing in
      let db = ref db
      and missing_len = Array.length missing in
      if missing_len > 0 then begin
        Array.iteri
          (fun i name -> !db.core.n_meta + i |> Hashtbl.add !db.meta_names_to_idx name)
          missing;
        let n_meta = !db.core.n_meta + missing_len in
        db :=
          { !db with
            core =
              { !db.core with
                n_meta;
                (* We have to resize all the relevant containers *)
                idx_to_meta_names = Array.append !db.core.idx_to_meta_names missing;
                meta = resize_string_array_array !db.core.n_cols n_meta !db.core.meta } }
      end;
      let num_header_fields = Array.length header
      and meta_indices =
        Array.mapi
          (fun i name ->
            if i = 0 then
              -1
            else
              Hashtbl.find !db.meta_names_to_idx name)
          header in
      begin try
        while true do
          let line = input_line input |> Tools.Split.on_char_as_array '\t' in
          incr line_num;
          (* A regular line. The first element is the vector name, the others the values of meta-data fields *)
          let l = Array.length line in
          if l <> num_header_fields then
            Wrong_number_of_columns (!line_num, l, num_header_fields) |> raise;
          add_empty_column_if_needed db line.(0);
          let col_idx = Hashtbl.find !db.col_names_to_idx line.(0) in
          Array.iteri
            (fun i name_idx ->
              if i > 0 then
                !db.core.meta.(col_idx).(name_idx) <- line.(i))
            meta_indices;
          if verbose && !line_num mod 10 = 0 then
            Printf.eprintf "\rFile '%s': Read %d lines%!" fname !line_num
        done
      with End_of_file ->
        close_in input;
        if verbose then
          Printf.eprintf "\rFile '%s': Read %d lines\n%!" fname !line_num
      end;
      !db
    exception Header_expected of string
    exception Wrong_format of int * string
    let add_files ?(verbose = false) db fnames =
      let db = ref db and n = List.length fnames in
      List.iteri
        (fun i fname ->
          let input = open_in fname and line_num = ref 0 and col_idx = ref 0 in
          (* Each file can contain one or more spectra *)
          begin try
            while true do
              let line_s = input_line input in
              incr line_num;
              let line = Tools.Split.on_char_as_array '\t' line_s in
              let l = Array.length line in
              if l <> 2 then
                Wrong_number_of_columns (!line_num, l, 2) |> raise;
              (* Each file must begin with a header *)
              if !line_num = 1 && line.(0) <> "" then
                Header_expected line_s |> raise;
              if line.(0) = "" then begin
                (* Header *)
                add_empty_column_if_needed db line.(1);
                col_idx := Hashtbl.find !db.col_names_to_idx line.(1)
              end else begin
                (* A regular line. The first element is the hash, the second one the count *)
                if Hashtbl.mem !db.row_names_to_idx line.(0) |> not then begin
                  let n_rows = !db.core.n_rows in
                  let aug_n_rows = n_rows + 1 in
                  Hashtbl.add !db.row_names_to_idx line.(0) n_rows;
                  db :=
                    { !db with
                      core =
                        { !db.core with
                          n_rows = aug_n_rows;
                          (* We have to resize all the relevant containers *)
                          idx_to_row_names = begin
                            let res = resize_string_array aug_n_rows !db.core.idx_to_row_names in
                            res.(n_rows) <- line.(0);
                            res
                          end;
                          storage = IBAVectorMisc.resize_array !db.core.n_cols aug_n_rows !db.core.storage } }
                end;
                let row_idx = Hashtbl.find !db.row_names_to_idx line.(0) in
                let v =
                  try
                    IBAVector.N.of_string line.(1)
                  with _ ->
                    Wrong_format (!line_num, line.(1)) |> raise in
                (* If there are repeated k-mers, we just accumulate them *)
                !db.core.storage.(!col_idx).IBAVector.+(row_idx) <- v
              end;
              if verbose && !line_num mod 10000 = 0 then
                Printf.eprintf "\r[%d/%d] File '%s': Read %d lines%!" (i + 1) n fname !line_num
            done
          with End_of_file ->
            close_in input;
            if verbose then
              Printf.eprintf "\r[%d/%d] File '%s': Read %d lines\n%!" (i + 1) n fname !line_num;
          end)
        fnames;
      !db
    (* *)
    let selected_from_regexps ?(verbose = false) db regexps =
      (* We iterate over the columns *)
      if verbose then
        Printf.eprintf "Selecting columns... [%!";
      List.iter
        (fun (what, _) ->
          if verbose && what <> "" && Hashtbl.find_opt db.meta_names_to_idx what = None then
            Printf.eprintf " (WARNING: Metadata field '%s' not found, no column will match)%!" what)
        regexps;
      let res = ref StringSet.empty in
      Array.iteri
        (fun n_col col_name ->
          if begin
            List.fold_left
              (fun start (what, regexp) ->
                start &&
                if what = "" then
                  (* Case of the label *)
                  Str.string_match regexp col_name 0
                else
                  match Hashtbl.find_opt db.meta_names_to_idx what with
                  | None ->
                    false
                  | Some found ->
                    assert (db.core.idx_to_meta_names.(found) = what);
                    Str.string_match regexp db.core.meta.(n_col).(found) 0)
              true regexps
          end then
            res := StringSet.add col_name !res)
        db.core.idx_to_col_names;
      if verbose then
        StringSet.iter (Printf.eprintf " '%s'%!") !res;
      if verbose then
        Printf.eprintf " ] done.\n%!";
      !res
    let selected_negate db selection =
      StringSet.diff (Array.to_list db.core.idx_to_col_names |> StringSet.of_list) selection
    (* It should be OK to have the same label on both LHS and RHS, as a temporary vector is used *)
    let add_sum_selected ?(threads = 1) ?(verbose = false) db new_label selection =
      (* Here we want no thresholding and linear statistics *)
      let transf = Transformation.of_parameters { which = "normalize"; threshold = 1; power = 1. } in
      let stats = Statistics.table_of_db ~threads ~verbose transf db and db = ref db in
      if verbose then
        Printf.eprintf "Adding/replacing vector '%s': [%!" new_label;
      add_empty_column_if_needed db new_label;
      let num_cols = ref 0 and max_norm = ref 0. in
      StringSet.iter
        (fun label ->
          if verbose then
            Printf.eprintf " '%s'%!" label;
          match Hashtbl.find_opt !db.col_names_to_idx label with
          | Some col_idx ->
            incr num_cols;
            max_norm := max !max_norm stats.col_stats.(col_idx).sum
          | None ->
            if verbose then
              Printf.eprintf "(NOT FOUND)%!")
        selection;
      let num_cols = !num_cols and max_norm = !max_norm in
      if verbose then
        Printf.eprintf " ] n=%d max_norm=%.16g [%!" num_cols max_norm;
      let red_n_rows = !db.core.n_rows - 1 and res = FBAVector.init !db.core.n_rows 0. in
      (* We normalise columns *separately* before adding them *)
      StringSet.iter
        (fun label ->
          match Hashtbl.find_opt !db.col_names_to_idx label with
          | Some col_idx ->
            if verbose then
              Printf.eprintf ".%!";
            let col = !db.core.storage.(col_idx)
            and norm = stats.col_stats.(col_idx).sum in
            for i = 0 to red_n_rows do
              res.FBAVector.+(i) <- IBAVector.N.to_float col.IBAVector.=(i) *. max_norm /. norm
            done
          | None -> ())
        selection;
      let new_col_idx = Hashtbl.find !db.col_names_to_idx new_label in
      let new_col = !db.core.storage.(new_col_idx) and norm = ref 0. in
      for i = 0 to red_n_rows do
        let res_i = IBAVector.N.of_float res.FBAVector.=(i) in
        (* Actual copy to storage *)
        new_col.IBAVector.=(i) <- res_i;
        norm := !norm +. IBAVector.N.to_float res_i
      done;
      if verbose then
        Printf.eprintf "] norm=%.16g\n%!" !norm;
      (* If metadata is present in the database, we generate some for the new column too *)
      if !db.core.n_meta > 0 then begin
        (* For each metadata field, we compute the intersection of the values across all selected columns *)
        let res = Array.make !db.core.n_meta StringSet.empty in
        StringSet.iter
          (fun label ->
            match Hashtbl.find_opt !db.col_names_to_idx label with
            | Some col_idx ->
              let col = !db.core.meta.(col_idx)
              and red_n_meta = !db.core.n_meta - 1 in
              for i = 0 to red_n_meta do
                res.(i) <- StringSet.add col.(i) res.(i)
              done
            | None -> ())
          selection;
        !db.core.meta.(new_col_idx) <-
          Array.map
            (fun set ->
              if StringSet.cardinal set = 1 then
                StringSet.min_elt set
              else
                "")
            res
      end;
      !db
    let remove_selected db selected =
      (* First, we compute the indices of the columns to be kept *)
      let idxs = ref IntSet.empty in
      Array.iteri
        (fun col_idx col_name ->
          if StringSet.mem col_name selected |> not then
            idxs := IntSet.add col_idx !idxs)
        db.core.idx_to_col_names;
      let idxs = !idxs in
      let n = IntSet.cardinal idxs in
      let filter_array a null =
        let filtered = Array.make n null and i = ref 0 in
        IntSet.iter
          (fun col_idx ->
            filtered.(!i) <- a.(col_idx);
            incr i)
          idxs;
        filtered in
      let core =
        { db.core with
          n_cols = n;
          idx_to_col_names = filter_array db.core.idx_to_col_names "";
          meta = filter_array db.core.meta [||];
          storage = filter_array db.core.storage IBAVector.empty } in
      { core = core;
        col_names_to_idx = invert_table core.idx_to_col_names;
        row_names_to_idx = invert_table core.idx_to_row_names;
        meta_names_to_idx = invert_table core.idx_to_meta_names }
    (* *)
    module TableFilter =
      struct
        type t = {
          print_row_names: bool;
          print_col_names: bool;
          print_metadata: bool;
          transpose: bool;
          transform: Transformation.t;
          print_zero_rows: bool;
          filter_columns: StringSet.t;
          precision: int
        }
        let default = {
          print_row_names = true;
          print_col_names = true;
          print_metadata = false;
          transpose = false;
          transform = Transformation.of_parameters { which = "normalize"; threshold = 1; power = 1. };
          print_zero_rows = false;
          filter_columns = StringSet.empty;
          precision = 15
        }      
      end
    let to_table
        ?(filter = TableFilter.default) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) db fname =
      let transform = Transformation.compute ~which:filter.transform
      and stats = Statistics.table_of_db ~threads ~verbose filter.transform db
      and output = open_out fname and meta = ref [] and rows = ref [] and cols = ref [] in
      (* We determine which rows and colunms should be output after all filters have been applied *)
      (*  Rows: metadata and k-mers *)
      if filter.print_metadata then
        Array.iteri
          (fun i meta_name ->
            (* There might be additional storage *)
            if i < db.core.n_meta then
              Tools.List.accum meta (meta_name, i))
          db.core.idx_to_meta_names;
      Array.iteri
        (fun i row_name ->
          (* There might be additional storage.
              Also, we only print the row if it has non-zero elements or if we are explicitly requested to do so *)
          if i < db.core.n_rows && (stats.row_stats.(i).sum > 0. || filter.print_zero_rows) then
            Tools.List.accum rows (row_name, i))
        db.core.idx_to_row_names;
      (* Columns *)
      Array.iteri
        (fun i col_name ->
          (* There might be additional storage, or we might need to remove some columns *)
          if i < db.core.n_cols && StringSet.mem col_name filter.filter_columns |> not then
            Tools.List.accum cols (col_name ,i))
        db.core.idx_to_col_names;
      let meta = Tools.Array.of_rlist !meta and rows = Tools.Array.of_rlist !rows
      and cols = Tools.Array.of_rlist !cols in
      let n_rows = Array.length rows and n_cols = Array.length cols in
      (* There must be at least one row and one column to print *)
      if (Array.length meta + n_rows) > 0 && n_cols > 0 then begin
        if filter.transpose then begin
          (* There is at least one row *)
          if filter.print_col_names then begin
            (* We print the row names *)
            if filter.print_row_names then
              Printf.fprintf output "\t";
            let first_done = ref false in
            Array.iter
              (fun (meta_name, _) ->
                Printf.fprintf output "%s%s" (if !first_done then "\t" else "") meta_name;
                first_done := true)
              meta;
            Array.iter
              (fun (row_name, _) ->
                Printf.fprintf output "%s%s" (if !first_done then "\t" else "") row_name;
                first_done := true)
              rows;
            Printf.fprintf output "\n"
          end;
          Printf.fprintf output "%!";
          let processed_cols = ref 0 and buf = Buffer.create 1048576 in
          Tools.Parallel.process_stream_chunkwise
            (fun () ->
              if !processed_cols < n_cols then
                let to_do = max 1 (elements_per_step / n_rows) |> min (n_cols - !processed_cols) in
                let new_processed_cols = !processed_cols + to_do in
                let res = !processed_cols, new_processed_cols - 1 in
                processed_cols := new_processed_cols;
                res
              else
                raise End_of_file)
            (fun (lo_col, hi_col) ->
              Buffer.clear buf;
              for i = lo_col to hi_col do
                let col_name, col_idx = cols.(i) in
                if filter.print_row_names then
                  Printf.bprintf buf "%s\t" col_name;
                let first_done = ref false in
                Array.iter
                  (fun (_, meta_idx) ->
                    Printf.bprintf buf "%s\"%s\"" (if !first_done then "\t" else "") db.core.meta.(col_idx).(meta_idx);
                    first_done := true)
                  meta;
                Array.iter
                  (fun (_, row_idx) ->
                    Printf.bprintf buf "%s%.*g" (if !first_done then "\t" else "") filter.precision begin
                      IBAVector.N.to_int db.core.storage.(col_idx).IBAVector.=(row_idx) |>
                        transform
                          ~col_num:db.core.n_cols ~col_stats:stats.col_stats.(col_idx)
                          ~row_num:db.core.n_rows ~row_stats:stats.row_stats.(row_idx)
                    end;
                    first_done := true)
                  rows;
                Printf.bprintf buf "\n"
              done;
              hi_col - lo_col + 1, Buffer.contents buf)
            (fun (n_processed, block) ->
              Printf.fprintf output "%s" block;
              let old_processed_cols = !processed_cols in
              processed_cols := !processed_cols + n_processed;
              if verbose && !processed_cols / 2 > old_processed_cols / 2 then (* We write one column at the time *)
                Printf.eprintf "\rWriting table to file '%s': done %d/%d lines%!"
                  fname !processed_cols n_cols)
            threads
        end else begin
          if filter.print_col_names then begin
            (* We print column names - there is at least one column *)
            if filter.print_row_names then
              Printf.fprintf output "\t";
            Array.iteri
              (fun i (col_name, _) ->
                Printf.fprintf output "%s%s" (if i = 0 then "" else "\t") col_name)
              cols;
            Printf.fprintf output "\n"
          end;
          (* We print metadata as lines *)
          Array.iter
            (fun (meta_name, meta_idx) ->
              if filter.print_row_names then
                Printf.fprintf output "%s\t" meta_name;
              Array.iteri
                (fun i (_, col_idx) ->
                  Printf.fprintf output "%s\"%s\"" (if i = 0 then "" else "\t") db.core.meta.(col_idx).(meta_idx))
                cols;
              Printf.fprintf output "\n")
            meta;
          Printf.fprintf output "%!";
          let processed_rows = ref 0 and buf = Buffer.create 1048576 in
          Tools.Parallel.process_stream_chunkwise
            (fun () ->
              if !processed_rows < n_rows then
                let to_do = max 1 (elements_per_step / n_cols) |> min (n_rows - !processed_rows) in
                let new_processed_rows = !processed_rows + to_do in
                let res = !processed_rows, new_processed_rows - 1 in
                processed_rows := new_processed_rows;
                res
              else
                raise End_of_file)
            (fun (lo_row, hi_row) ->
              Buffer.clear buf;
              for i = lo_row to hi_row do
                let row_name, row_idx = rows.(i) in
                if filter.print_row_names then
                  Printf.bprintf buf "%s\t" row_name;
                Array.iteri
                  (fun j (_, col_idx) ->
                    Printf.bprintf buf "%s%.*g" (if j = 0 then "" else "\t") filter.precision begin
                      IBAVector.N.to_int db.core.storage.(col_idx).IBAVector.=(row_idx) |>
                        transform
                          ~col_num:db.core.n_cols ~col_stats:stats.col_stats.(col_idx)
                          ~row_num:db.core.n_rows ~row_stats:stats.row_stats.(row_idx)
                    end)
                  cols;
                Printf.bprintf buf "\n"
              done;
              hi_row - lo_row + 1, Buffer.contents buf)
            (fun (n_processed, block) ->
              Printf.fprintf output "%s" block;
              let old_processed_rows = !processed_rows in
              processed_rows := !processed_rows + n_processed;
              if verbose && !processed_rows / 10000 > old_processed_rows / 10000 then
                Printf.eprintf "\rWriting table to file '%s': done %d/%d lines%!" fname !processed_rows n_rows)
            threads
        end;
        let n_done =
          if filter.transpose then
            n_cols
          else
            n_rows in
        if verbose then
          Printf.eprintf "\rWriting table to file '%s': done %d/%d lines.\n%!" fname n_done n_done
      end;
      close_out output
    let make_filename_table = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopCounter.txt"

  end

type to_do_t =
  | Empty
  | Of_file of string
  | Add_meta of string
  | Add_files of string list
  | Add_sum_selected of string (* The new label *)
  | Remove_selected
  | Summary
  | Selected_from_labels of StringSet.t
  | Selected_from_regexps of (string * Str.regexp) list
  | Selected_negate
  | Selected_print
  | Selected_clear
  | Selected_to_filter
  | Table_emit_row_names of bool
  | Table_emit_col_names of bool
  | Table_emit_metadata of bool
  | Table_transpose of bool
  | Table_transform of Transformation.parameters_t
  | Table_emit_zero_rows of bool
  | Table_precision of int
  | To_table of string
  | To_file of string
  
module Defaults =
  struct
    let filter = KMerDB.TableFilter.default
    let threads = Tools.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    (* This is just to correctly propagate values in the program *)
    let transform = Transformation.to_parameters Defaults.filter.transform |> ref
    let program = ref []
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let version = "0.27"

let header =
  Printf.sprintf begin
    "This is the KPopCountDB program (version %s)\n%!" ^^
    " (c) 2020-2022 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!"
  end version

let _ =
  let module TA = Tools.Argv in
  let module TS = Tools.Split in
  let module TL = Tools.List in
  TA.set_header header;
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    [], None, [ "Actions (executed delayed and in order of specification):" ], TA.Optional, (fun _ -> ());
    [ "-e"; "-E"; "--empty" ],
      None,
      [ "put an empty database into the register" ],
      TA.Optional,
      (fun _ -> Empty |> TL.accum Parameters.program);
    [ "-i"; "-I"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load to the register the database present in the specified file";
        " (which must have extension .KPopCounter)" ],
      TA.Optional,
      (fun _ -> Of_file (TA.get_parameter () |> KMerDB.make_filename_binary) |> TL.accum Parameters.program);
    [ "-m"; "-M"; "--metadata"; "--add-metadata" ],
      Some "<metadata_table_file_name>",
      [ "add to the database present in the register metadata from the specified file" ],
      TA.Optional,
      (fun _ -> Add_meta (TA.get_parameter ()) |> TL.accum Parameters.program);
    [ "-f"; "-F"; "--files"; "--add-files" ],
      Some "<k-mer_table_file_name>[','...','<k-mer_table_file_name>]",
      [ "add to the database present in the register k-mers from the specified files" ],
      TA.Optional,
      (fun _ -> Add_files (TA.get_parameter () |> TS.on_char_as_list ',') |> TL.accum Parameters.program);
    [ "--summary" ],
      None,
      [ "print a summary of the database present in the register" ],
      TA.Optional,
      (fun _ -> Summary |> TL.accum Parameters.program);
    [ "-l"; "-L"; "--labels"; "--selection-from-labels" ],
      Some "<vector_label>[','...','<vector_label>]",
      [ "put into the selection register the specified labels" ],
      TA.Optional,
      (fun _ ->
        let labels = TA.get_parameter () in
        if labels <> "" then
        Selected_from_labels (labels |> TS.on_char_as_list ',' |> StringSet.of_list)
          |> TL.accum Parameters.program);
    [ "-r"; "-R"; "--regexps"; "--selection-from-regexps" ],
      Some "<metadata_field>'~'<regexp>[','...','<metadata_field>'~'<regexp>]",
      [ "put into the selection register the labels of the vectors";
        "whose metadata fields match the specified regexps.";
        "An empty metadata field matches the labels" ],
      TA.Optional,
      (fun _ ->
        Selected_from_regexps begin
          List.map
            (fun l ->
              let res = TS.on_char_as_list '~' l in
              if List.length res <> 2 then begin
                TA.usage ();
                List.length res |>
                  Printf.sprintf "Option '-r': Wrong number of fields in list (expected 2, found %d)" |>
                  TA.parse_error (* parse_error exits the program *)
              end;
              List.nth res 0, List.nth res 1 |> Str.regexp)
            (TA.get_parameter () |> TS.on_char_as_list ',')
        end |> TL.accum Parameters.program);
    [ "-a"; "-A"; "--add-sum-selection"; "--selection-add-sum" ],
      Some "<new_vector_label>",
      [ "add to the database present in the register (or replace if the new label exists)";
        "a linear combination of the vectors whose labels are in the selection register" ],
      TA.Optional,
      (fun _ -> Add_sum_selected (TA.get_parameter ()) |> TL.accum Parameters.program);
    [ "-d"; "-D"; "--delete"; "--selection-remove" ],
      None,
      [ "remove from the table the vectors whose labels are in the selection register" ],
      TA.Optional,
      (fun _ -> Remove_selected |> TL.accum Parameters.program);
    [ "-n"; "-N"; "--selection-negate" ],
      None,
      [ "negate the labels that are present in the selection register" ],
      TA.Optional,
      (fun _ -> Selected_negate |> TL.accum Parameters.program);
    [ "-p"; "-P"; "--selection-print" ],
      None,
      [ "print the labels that are present in the selection register" ],
      TA.Optional,
      (fun _ -> Selected_print |> TL.accum Parameters.program);
    [ "-c"; "-C"; "--selection-clear" ],
      None,
      [ "purges the selection register" ],
      TA.Optional,
      (fun _ -> Selected_clear |> TL.accum Parameters.program);
    [ "-s"; "-S"; "--selection-to-table-filter" ],
      None,
      [ "filters out vectors whose labels are present in the selection register";
        "when writing the database as a tab-separated file" ],
      TA.Optional,
      (fun _ -> Selected_to_filter |> TL.accum Parameters.program);
    [ "--table-emit-row-names" ],
      Some "'true'|'false'",
      [ "whether to emit row names for the database present in the register";
        "when writing it as a tab-separated file" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_row_names),
      (fun _ -> Table_emit_row_names (TA.get_parameter_boolean ()) |> TL.accum Parameters.program);
    [ "--table-emit-col-names" ],
      Some "'true'|'false'",
      [ "whether to emit column names for the database present in the register";
        "when writing it as a tab-separated file" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_col_names),
      (fun _ -> Table_emit_col_names (TA.get_parameter_boolean ()) |> TL.accum Parameters.program);
    [ "--table-emit-metadata" ],
      Some "'true'|'false'",
      [ "whether to emit metadata for the database present in the register";
        "when writing it as a tab-separated file" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_metadata),
      (fun _ -> Table_emit_metadata (TA.get_parameter_boolean ()) |> TL.accum Parameters.program);
    [ "--table-transpose" ],
      Some "'true'|'false'",
      [ "whether to transpose the database present in the register";
        "before writing it as a tab-separated file";
        " (if 'true' : rows are vector names, columns are (metadata and) k-mer names;";
        "  if 'false': rows are (metadata and) k-mer names, columns are vector names)" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.transpose),
      (fun _ -> Table_transpose (TA.get_parameter_boolean ()) |> TL.accum Parameters.program);
    [ "--table-threshold" ],
      Some "<non_negative_integer>",
      [ "set to zero all counts that are less than this threshold";
        "before transforming and outputting them" ],
      TA.Default (fun () -> (Transformation.to_parameters Defaults.filter.transform).threshold |> string_of_int),
      (fun _ ->
         Table_transform { !Parameters.transform with threshold = TA.get_parameter_int_non_neg () }
           |> TL.accum Parameters.program);
    [ "--table-power" ],
      Some "<non_negative_float>",
      [ "raise counts to this power before transforming and outputting them.";
        "A power of 0 when the 'pseudocount' method is used";
        "performs a logarithmic transformation" ],
      TA.Default (fun () -> (Transformation.to_parameters Defaults.filter.transform).power |> string_of_float),
      (fun _ ->
         Table_transform { !Parameters.transform with power = TA.get_parameter_float_non_neg () }
           |> TL.accum Parameters.program);
    [ "--table-transform"; "--table-transformation" ],
      Some "'none'|'normalize'|'pseudocount'|'clr'",
      [ "transformation to apply to table elements before outputting them" ],
      TA.Default (fun () -> (Transformation.to_parameters Defaults.filter.transform).which),
      (fun _ ->
         Table_transform { !Parameters.transform with which = TA.get_parameter () } |> TL.accum Parameters.program);
    [ "--table-emit-zero-rows" ],
      Some "'true'|'false'",
      [ "whether to emit rows whose elements are all zero";
        "when writing the database as a tab-separated file" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_zero_rows),
      (fun _ -> Table_emit_zero_rows (TA.get_parameter_boolean ()) |> TL.accum Parameters.program);
    [ "--table-set-precision" ],
      Some "<positive_integer>",
      [ "set the number of precision digits to be used when outputting counts" ],
      TA.Default (fun () -> string_of_int Defaults.filter.precision),
      (fun _ -> Table_precision (TA.get_parameter_int_pos ()) |> TL.accum Parameters.program);
    [ "-t"; "--table" ],
      Some "<file_prefix>",
      [ "write the database present in the register as a tab-separated file";
        " (rows are k-mer names, columns are vector names;";
        "  the file will be given extension .KPopCounter.txt)" ],
      TA.Optional,
      (fun _ -> To_table (TA.get_parameter () |> KMerDB.make_filename_table) |> TL.accum Parameters.program);
    [ "-o"; "-O"; "--output" ],
      Some "<binary_file_prefix>",
      [ "dump the database present in the register to the specified file";
        " (which will be given extension .KPopCounter)" ],
      TA.Optional,
      (fun _ -> To_file (TA.get_parameter () |> KMerDB.make_filename_binary) |> TL.accum Parameters.program);
    [], None, [ "Miscellaneous (executed immediately):" ], TA.Optional, (fun _ -> ());
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (fun () -> string_of_bool !Parameters.verbose),
      (fun _ -> Parameters.verbose := true);
    (* Hidden option to emit help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 0)
  ];
  let program = List.rev !Parameters.program in
  if program = [] then begin
    TA.usage ();
    exit 0
  end;
  if !Parameters.verbose then
    TA.header ();
  (* These are the three registers available to the program *)
  let current = KMerDB.make_empty () |> ref and selected = ref StringSet.empty
  and filter = ref KMerDB.TableFilter.default in
  
  (* The addition of an exception handler would be nice *)
  
  List.iter
    (function
      | Empty ->
        current := KMerDB.make_empty ()
      | Of_file fname ->
        current := KMerDB.of_binary ~verbose:!Parameters.verbose fname
      | Add_meta fname ->
        current := KMerDB.add_meta ~verbose:!Parameters.verbose !current fname
      | Add_files fnames ->
        current := KMerDB.add_files ~verbose:!Parameters.verbose !current fnames
      | Add_sum_selected new_label ->
        current :=
          KMerDB.add_sum_selected ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !current new_label !selected
      | Remove_selected ->
        current := KMerDB.remove_selected !current !selected
      | Summary ->
        KMerDB.output_summary ~verbose:!Parameters.verbose !current
      | Selected_from_labels labels ->
        selected := labels
      | Selected_from_regexps regexps ->
        selected := KMerDB.selected_from_regexps ~verbose:!Parameters.verbose !current regexps
      | Selected_negate ->
        selected := KMerDB.selected_negate !current !selected
      | Selected_print ->
        Printf.eprintf "Currently selected vectors = [";
        StringSet.iter (Printf.eprintf " '%s'%!") !selected;
        Printf.eprintf " ].\n%!"
      | Selected_clear ->
        selected := StringSet.empty
      | Selected_to_filter ->
        filter := { !filter with filter_columns = !selected }
      | Table_emit_row_names print_row_names ->
        filter := { !filter with print_row_names }
      | Table_emit_col_names print_col_names ->
        filter := { !filter with print_col_names }
      | Table_emit_metadata print_metadata ->
        filter := { !filter with print_metadata }
      | Table_transpose transpose ->
        filter := { !filter with transpose }
      | Table_transform transform ->
        filter := { !filter with transform = Transformation.of_parameters transform }
      | Table_emit_zero_rows print_zero_rows ->
        filter := { !filter with print_zero_rows }
      | Table_precision precision ->
        filter := { !filter with precision }
      | To_table fname ->
        KMerDB.to_table ~filter:!filter ~threads:!Parameters.threads ~verbose:!Parameters.verbose !current fname
      | To_file fname ->
        KMerDB.to_binary ~verbose:!Parameters.verbose !current fname)
    program

