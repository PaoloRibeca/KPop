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

module Numbers = BiOCamLib.Numbers
module Processes = BiOCamLib.Processes

(* We put this one here for lack of a better place *)
module Spectra =
  struct
    let make_filename = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopSpectra.txt"
  end

(* For counts. We assume each count to be < 2^31 *)
module IBAVector = Numbers.Bigarray.Vector (
  struct
    include Numbers.Int32
    type elt_t = Bigarray.int32_elt
    let elt = Bigarray.Int32
  end
)
(* For normalisations and combined spectra. We assume normalisations might be > 2^31 *)
module FBAVector = Numbers.Bigarray.Vector (
  struct
    include Numbers.Float
    type elt_t = Bigarray.float64_elt
    let elt = Bigarray.Float64
  end
)

include (
  struct
    let ( .@() ) = IBAVector.( .@() )
    let ( .@()<- ) = IBAVector.( .@()<- )
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
    module ColOrRow =
      struct
        type t =
        | Col
        | Row
        let to_string = function
        | Col -> "column"
        | Row -> "row"
      end
    module Transformation =
      struct
        type statistics_t = {
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
        type statistics_table_t = {
          col_stats: statistics_t array;
          row_stats: statistics_t array
        }
        type t =
          | Binary of float
          | Power of float * float
          | CLR of float * float
          | Pseudo of float * float
        exception Invalid_transformation of string * float * float
        let epsilon = 0.1
        let [@warning "-27"] compute ~which ~col_num ~col_stats ~row_num ~row_stats counts =
          let counts = float_of_int counts
          and threshold =
            match which with
            | Binary threshold | Power (threshold, _) | CLR (threshold, _) | Pseudo (threshold, _) ->
              if threshold < 1. then
                threshold *. col_stats.sum
              else
                threshold in
          match which with
          | Binary _ ->
            if counts >= threshold then
              1.
            else
              0.
          | Power (_, 1.) -> (* Optimisation *)
            if counts >= threshold then
              counts
            else
              0.
          | Power (_, power) ->
            if counts >= threshold then
              counts ** power
            else
              0.
          | CLR (_, power) ->
            let v =
              if counts >= threshold then
                counts
              else
                0. in
            let v = max v epsilon in
            (log v *. power) -. (col_stats.sum_log /. float_of_int col_stats.non_zero)
          | Pseudo (thr, power) ->
            if power < 0. then
              Invalid_transformation ("pseudocounts", thr, power) |> raise;
            let v =
              if power = 0. then
                (float_of_int col_stats.max) *. log ((counts +. 1.) /. threshold)
              else begin
                let red_threshold = threshold -. 1. |> max 0. in
                let c_p = red_threshold ** power in
                if power < 1. then
                  ((counts ** power) -. c_p) *. ((float_of_int col_stats.max) ** (1. -. power)) /. power
                else
                  ((counts ** power) -. c_p) /. ((threshold ** power) -. c_p)
              end in
            floor v /. col_stats.sum |> max 0.
        type parameters_t = {
          which: string;
          threshold: float;
          power: float
        }
        exception Unknown_transformation of string
        let of_parameters { which; threshold; power } =
          match which with
          | "binary" ->
            Binary threshold
          | "power" | "pow" ->
            Power (threshold, power)
          | "clr" | "CLR" ->
            CLR (threshold, power)
          | "pseudocounts" | "pseudo" ->
            Pseudo (threshold, power)
          | s ->
            Unknown_transformation s |> raise
        let to_parameters = function
          | Binary threshold -> { which = "binary"; threshold; power = 1. }
          | Power (threshold, power) -> { which = "power"; threshold; power }
          | CLR (threshold, power) -> { which = "clr"; threshold; power }
          | Pseudo (threshold, power) -> { which = "pseudocounts"; threshold; power }
      end
    (* Implementation function *)
    let stats_table_of_core_db ?(threads = 1) ?(verbose = false) transf core =
      let { Transformation.threshold; power; _ } = Transformation.to_parameters transf in
      let compute_one what n =
        let red_len =
          match what with
          | ColOrRow.Col ->
            core.n_rows - 1
          | Row ->
            core.n_cols - 1 in
        (* Here in order to compute the threshold we have to pre-compute the normalisation *)
        let sum = ref 0. in
        for i = 0 to red_len do
          let v =
            match what with
            | Col ->
              IBAVector.N.to_int core.storage.(n).@(i)
            | Row ->
              IBAVector.N.to_int core.storage.(i).@(n) in
          sum := !sum +. (float_of_int v ** power)
        done;
        let threshold =
          if threshold < 1. then
            threshold *. !sum
          else
            threshold in
        let non_zero = ref 0 and min = ref 0 and max = ref 0 and sum = ref 0. and sum_log = ref 0. in
        for i = 0 to red_len do
          let v =
            match what with
            | Col ->
              IBAVector.N.to_int core.storage.(n).@(i)
            | Row ->
              IBAVector.N.to_int core.storage.(i).@(n) in
          let f_v = float_of_int v in
          if f_v >= threshold then begin
            incr non_zero;
            min := Stdlib.min !min v;
            max := Stdlib.max !max v;
            sum := !sum +. (f_v ** power);
            sum_log := !sum_log +. (log f_v *. power)
          end
        done;
        { Transformation.non_zero = !non_zero;
          min = !min;
          max = !max;
          sum = !sum;
          sum_log = !sum_log } in
      let compute_all what =
        let n =
          match what with
          | ColOrRow.Col ->
            core.n_cols
          | Row ->
            core.n_rows in
        let step = n / threads / 5 |> max 1 and processed = ref 0 and res = ref [] in
        Processes.Parallel.process_stream_chunkwise
          (fun () ->
            if verbose then
              Printf.eprintf "%s\r(%s): Computing %s statistics [%d/%d]%!"
                String.TermIO.clear __FUNCTION__ (ColOrRow.to_string what) !processed n;
            let to_do = n - !processed in
            if to_do > 0 then begin
              let to_do = min to_do step in
              let res = !processed, to_do in
              processed := !processed + to_do;
              res
            end else
              raise End_of_file)
          (fun (processed, to_do) ->
            let res = ref [] in
            for i = 0 to to_do - 1 do
              processed + i |> compute_one what |> List.accum res
            done;
            Array.of_rlist !res)
          (List.accum res)
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
          Printf.eprintf "%s\r(%s): Computing %s statistics [%d/%d]\n%!"
            String.TermIO.clear __FUNCTION__ (ColOrRow.to_string what) n n;
        res in
      { Transformation.col_stats = compute_all Col;
        row_stats = compute_all Row }
    type t = {
      core: marshalled_t;
      col_names_to_idx: int StringHashtbl.t;
      row_names_to_idx: int StringHashtbl.t;
      meta_names_to_idx: int StringHashtbl.t
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
      col_names_to_idx = StringHashtbl.create 16;
      row_names_to_idx = StringHashtbl.create 16;
      meta_names_to_idx = StringHashtbl.create 16
    }
    let output_summary ?(verbose = false) db =
      Printf.eprintf "[Spectrum labels (%d)]:" db.core.n_cols;
      Array.iteri
        (fun i s ->
          if i < db.core.n_cols then
            Printf.eprintf " '%s'" s)
        db.core.idx_to_col_names;
      Printf.eprintf "\n%!";
      if verbose then begin
        Printf.eprintf "[K-mer hashes (%d)]:" db.core.n_rows;
        Array.iteri
          (fun i s ->
            if i < db.core.n_rows then
              Printf.eprintf " '%s'" s)
          db.core.idx_to_row_names;
        Printf.eprintf "\n%!"
      end;
      Printf.eprintf "[Meta-data fields (%d)]:" db.core.n_meta;
      Array.iteri
        (fun i s ->
          if i < db.core.n_meta then
            Printf.eprintf " '%s'" s)
        db.core.idx_to_meta_names;
      Printf.eprintf "\n%!"
    (* *)
    let resize_string_array ?(is_buffer = true) n a =
      Array.resize ~is_buffer n "" a
    (* *)
    let _resize_t_array_ ?(is_buffer = true) length resize create_null nx ny a =
      let lx = Array.length a in
      let eff_nx =
        if is_buffer then begin
          if lx < nx then
            max nx (lx * 14 / 10)
          else
            lx
        end else
          nx
      and eff_ny =
        if is_buffer then begin
          if lx > 0 then begin
            (* We assume all bigarrays to have the same size *)
            let ly = length a.(0) in
            if ly < ny then
              max ny (ly * 14 / 10)
            else
              ly
          end else
            ny
        end else
          ny in
      (*Printf.eprintf "(%s): Resizing to (%d,%d) - asked (%d,%d)...\n%!" __FUNCTION__ eff_nx eff_ny nx ny;*)
      if eff_nx > lx then
        Array.append
          (* We need to provide the is_buffer argument like this because of type resolution *)
          (Array.map (resize ?is_buffer:(Some false) eff_ny) a)
          (Array.init (eff_nx - lx) (fun _ -> create_null eff_ny))
      else if eff_nx < lx then
        Array.map (resize ?is_buffer:(Some false) eff_ny) (Array.sub a 0 eff_nx)
      else (* eff_nx = lx *)
        if lx > 0 && eff_ny = length a.(0) then
          a
        else
          Array.map (resize ?is_buffer:(Some false) eff_ny) a
    let resize_string_array_array ?(is_buffer = true) =
      _resize_t_array_ ~is_buffer Array.length resize_string_array (fun l -> Array.make l "")
    module BAVectorMisc (M: Numbers.Vector_t) =
      struct
        let resize ?(is_buffer = true) n =
          M.resize ~is_buffer ~fill_with:M.N.zero n
        let resize_array ?(is_buffer = true) =
          _resize_t_array_ ~is_buffer M.length resize (fun l -> M.make l M.N.zero)
      end
    module IBAVectorMisc = BAVectorMisc (IBAVector)
    module FBAVectorMisc = BAVectorMisc (FBAVector)
    (* Utility functions *)
    let invert_table a =
      let res = StringHashtbl.create (Array.length a) in
      Array.iteri (fun i name -> StringHashtbl.add res name i) a;
      res
    let add_empty_column_if_needed db label =
      let n_cols = !db.core.n_cols in
      let aug_n_cols = n_cols + 1 in
      if StringHashtbl.mem !db.col_names_to_idx label |> not then begin
        StringHashtbl.add !db.col_names_to_idx label n_cols; (* THIS ONE CHANGES !db *)
        db := {
          !db with
          core = {
            !db.core with
            n_cols = aug_n_cols;
            (* We have to resize all the relevant containers *)
            idx_to_col_names = Array.append !db.core.idx_to_col_names [| label |];
            meta = resize_string_array_array ~is_buffer:true aug_n_cols !db.core.n_meta !db.core.meta;
            storage = IBAVectorMisc.resize_array ~is_buffer:true aug_n_cols !db.core.n_rows !db.core.storage
          }
        }
      end
    (* *)
    let archive_version = "2022-04-03"
    (* *)
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopCounter"
    exception Incompatible_archive_version of string * string
    let to_binary ?(verbose = false) db prefix =
      let fname = make_filename_binary prefix in
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      output_value output "KPopCounter";
      output_value output archive_version;
      output_value output {
        db.core with
        (* We have to truncate all the containers *)
        idx_to_col_names = resize_string_array ~is_buffer:false db.core.n_cols db.core.idx_to_col_names;
        idx_to_row_names = resize_string_array ~is_buffer:false db.core.n_rows db.core.idx_to_row_names;
        idx_to_meta_names = resize_string_array ~is_buffer:false db.core.n_meta db.core.idx_to_meta_names;
        meta = resize_string_array_array ~is_buffer:false db.core.n_cols db.core.n_meta db.core.meta;
        storage = IBAVectorMisc.resize_array ~is_buffer:false db.core.n_cols db.core.n_rows db.core.storage
      };
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) prefix =
      let fname = make_filename_binary prefix in
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
    (* *)
    exception Wrong_number_of_columns of int * int * int
    let add_meta ?(verbose = false) db fname =
      let input = open_in fname and line_num = ref 0 in
      let header =
        input_line input
          |> String.Split.on_char_as_array '\t' |> Array.map Matrix.Base.strip_external_quotes_and_check in
      incr line_num;
      (* We add the names *)
      let missing = ref [] in
      Array.iteri
        (fun i name ->
          if i > 0 && StringHashtbl.mem db.meta_names_to_idx name |> not then
            List.accum missing name)
        header;
      let missing = Array.of_rlist !missing in
      let db = ref db
      and missing_len = Array.length missing in
      if missing_len > 0 then begin
        Array.iteri
          (fun i name -> !db.core.n_meta + i |> StringHashtbl.add !db.meta_names_to_idx name)
          missing;
        let n_meta = !db.core.n_meta + missing_len in
        db := {
          !db with
          core = {
            !db.core with
            n_meta;
            (* We have to resize all the relevant containers *)
            idx_to_meta_names = Array.append !db.core.idx_to_meta_names missing;
            meta = resize_string_array_array ~is_buffer:true !db.core.n_cols n_meta !db.core.meta
          }
        }
      end;
      let num_header_fields = Array.length header
      and meta_indices =
        Array.mapi
          (fun i name ->
            if i = 0 then
              -1
            else
              StringHashtbl.find !db.meta_names_to_idx name)
          header in
      begin try
        while true do
          let line =
            input_line input
              |> String.Split.on_char_as_array '\t' |> Array.map Matrix.Base.strip_external_quotes_and_check in
          incr line_num;
          (* A regular line. The first element is the spectrum name, the others the values of meta-data fields *)
          let l = Array.length line in
          if l <> num_header_fields then
            Wrong_number_of_columns (!line_num, l, num_header_fields) |> raise;
          add_empty_column_if_needed db line.(0);
          let col_idx = StringHashtbl.find !db.col_names_to_idx line.(0) in
          Array.iteri
            (fun i name_idx ->
              if i > 0 then
                !db.core.meta.(col_idx).(name_idx) <- line.(i))
            meta_indices;
          if verbose && !line_num mod 10 = 0 then
            Printf.eprintf "%s\r(%s): File '%s': Read %d lines%!"
              String.TermIO.clear __FUNCTION__ fname !line_num
        done
      with End_of_file ->
        close_in input;
        if verbose then
          Printf.eprintf "%s\r(%s): File '%s': Read %d lines\n%!"
            String.TermIO.clear __FUNCTION__ fname !line_num
      end;
      !db
    exception Header_expected of string
    exception Wrong_format of int * string
    (* Here we cannot easily parallelise because of DB memory management *)
    let add_files ?(verbose = false) db prefixes =
      let db = ref db and n = List.length prefixes and num_spectra = ref (-1) in
      List.iteri
        (fun i prefix ->
          let fname = Spectra.make_filename prefix in
          let input = open_in fname and line_num = ref 0 and col_idx = ref 0 in
          (* Each file can contain one or more spectra *)
          begin try
            while true do
              let line_s = input_line input in
              incr line_num;
              let line = String.Split.on_char_as_array '\t' line_s in
              let l = Array.length line in
              if l <> 2 then
                Wrong_number_of_columns (!line_num, l, 2) |> raise;
              (* Each file must begin with a header *)
              if !line_num = 1 && line.(0) <> "" then
                Header_expected line_s |> raise;
              if line.(0) = "" then begin
                (* Header *)
                let label = Matrix.Base.strip_external_quotes_and_check line.(1) in
                add_empty_column_if_needed db label;
                col_idx := StringHashtbl.find !db.col_names_to_idx label;
                incr num_spectra;
                if verbose then
                  Printf.eprintf "%s\r(%s): [%d/%d] File '%s': Read %d %s on %d %s%!"
                    String.TermIO.clear __FUNCTION__ (i + 1) n fname
                    !num_spectra (String.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                    !line_num (String.pluralize_int "line" !line_num)
              end else begin
                (* A regular line. The first element is the hash, the second one the count *)
                if StringHashtbl.mem !db.row_names_to_idx line.(0) |> not then begin
                  let n_rows = !db.core.n_rows in
                  let aug_n_rows = n_rows + 1 in
                  StringHashtbl.add !db.row_names_to_idx line.(0) n_rows;
                  db := {
                    !db with
                    core = {
                      !db.core with
                      n_rows = aug_n_rows;
                      (* We have to resize all the relevant containers *)
                      idx_to_row_names = begin
                        let res = resize_string_array ~is_buffer:true aug_n_rows !db.core.idx_to_row_names in
                        res.(n_rows) <- line.(0);
                        res
                      end;
                      storage = IBAVectorMisc.resize_array ~is_buffer:true !db.core.n_cols aug_n_rows !db.core.storage
                    }
                  }
                end;
                let row_idx = StringHashtbl.find !db.row_names_to_idx line.(0) in
                let v =
                  try
                    IBAVector.N.of_string line.(1)
                  with _ ->
                    Wrong_format (!line_num, line.(1)) |> raise in
                (* If there are repeated k-mers, we just accumulate them *)
                !db.core.storage.(!col_idx).IBAVector.+(row_idx) <- v
              end
            done
          with End_of_file ->
            close_in input;
            incr num_spectra;
            if verbose then
              Printf.eprintf "%s\r(%s): [%d/%d] File '%s': Read %d %s on %d %s%s%!"
                String.TermIO.clear __FUNCTION__ (i + 1) n fname
                !num_spectra (String.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                !line_num (String.pluralize_int "line" !line_num) (if i + 1 = n then ".\n" else "")
          end)
        prefixes;
      !db
    (* *)
    let selected_from_regexps ?(verbose = false) db regexps =
      (* We iterate over the columns *)
      if verbose then
        Printf.eprintf "(%s): Selecting columns... [%!" __FUNCTION__;
      List.iter
        (fun (what, _) ->
          if verbose && what <> "" && StringHashtbl.find_opt db.meta_names_to_idx what = None then
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
                  match StringHashtbl.find_opt db.meta_names_to_idx what with
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
    exception Unknown_combination_criterion of string
    module CombinationCriterion =
      struct
        type t = RescaledMean | RescaledMedian
        let of_string = function
          | "mean" -> RescaledMean
          | "median" -> RescaledMedian
          | w -> Unknown_combination_criterion w |> raise
        let to_string = function
          | RescaledMean -> "mean"
          | RescaledMedian -> "median"
      end
    (* It should be OK to have the same label on both LHS and RHS,
        as a temporary vector is used to generate the combination *)
    let add_combined_selected ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        db new_label selection criterion =
      (* Here we want no thresholding and linear statistics *)
      let transf = Transformation.of_parameters { which = "power"; threshold = 1.; power = 1. } in
      let stats = stats_table_of_core_db ~threads ~verbose transf db.core and db = ref db in
      if verbose then
        Printf.eprintf "(%s): Adding/replacing spectrum '%s': [%!" __FUNCTION__ new_label;
      (* We allocate the result *)
      add_empty_column_if_needed db new_label;
      let new_col_idx = StringHashtbl.find !db.col_names_to_idx new_label in
      let new_col = !db.core.storage.(new_col_idx) in
      (* Computing valid labels and maximum normalisation across columns *)
      let found_cols = ref [] and max_norm = ref 0. in
      StringSet.iter
        (fun label ->
          if verbose then
            Printf.eprintf " '%s'%!" label;
          (* Some labels might be invalid *)
          match StringHashtbl.find_opt !db.col_names_to_idx label with
          | Some col_idx ->
            List.accum found_cols col_idx;
            max_norm := max !max_norm stats.col_stats.(col_idx).sum
          | None ->
            if verbose then
              Printf.eprintf "(NOT FOUND)%!")
        selection;
      let found_cols = Array.of_list !found_cols in (* We don't really care about the order here *)
      let num_found_cols = Array.length found_cols and max_norm = !max_norm and norm = ref 0. in
      if verbose then
        Printf.eprintf " ] n_found=%d max_norm=%.16g.\n%!" num_found_cols max_norm;
      let n_rows = !db.core.n_rows and processed_rows = ref 0 in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_rows then
            let to_do = max 1 (elements_per_step / num_found_cols) |> min (n_rows - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          let module FVF = Numbers.Frequencies.Vector(Numbers.Float)(Numbers.MakeComparableNumber) in
          (* We need one histogram per row to be able to compute statistics such as the median *)
          let row_combinators = Array.init (hi_row - lo_row + 1) (fun _ -> FVF.empty ~non_negative:true ()) in
          for i = lo_row to hi_row do
            (* We iterate over valid columns *)
            Array.iter
              (fun col_idx ->
                (* We normalise columns *separately* before combining them *)
                let col = !db.core.storage.(col_idx) and norm = stats.col_stats.(col_idx).sum in
                (* All counts are non-negative *)
                if norm > 0. then
                  (* We add the renormalised sum to the suitable row histogram *)
                  IBAVector.N.to_float col.@(i) *. max_norm /. norm |> FVF.add row_combinators.(i - lo_row))
              found_cols
          done;
          (* For each row histogram in the input range, we now generate a combination and pass it along *)
          lo_row,
          Array.map
            (fun combinator ->
              match criterion with
              | CombinationCriterion.RescaledMean ->
                FVF.sum combinator
              | RescaledMedian ->
                FVF.median combinator *. float_of_int num_found_cols)
            row_combinators)
        (fun (lo_row, block) ->
          let n_processed = Array.length block in
          for i = lo_row to lo_row + n_processed - 1 do
            let res_i = block.(i - lo_row) in
            Numbers.Float.(norm ++ res_i);
            let res_i = Int32.of_float res_i in
            (* Actual copy to storage *)
            new_col.@(i) <- res_i
          done;
          let old_processed_rows = !processed_rows in
          processed_rows := !processed_rows + n_processed;
          if verbose && !processed_rows / 10000 > old_processed_rows / 10000 then
            Printf.eprintf "%s\r(%s): Combining spectra: done %d/%d lines%!"
              String.TermIO.clear __FUNCTION__ !processed_rows n_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Combining spectra: done %d/%d lines. Norm=%.16g\n%!"
          String.TermIO.clear __FUNCTION__ !processed_rows n_rows !norm;
      (* If metadata is present in the database, we generate some for the new column too *)
      if !db.core.n_meta > 0 then begin
        (* For each metadata field, we compute the intersection of the values across all selected columns *)
        let res = Array.make !db.core.n_meta StringSet.empty in
        StringSet.iter
          (fun label ->
            match StringHashtbl.find_opt !db.col_names_to_idx label with
            | Some col_idx ->
              let col = !db.core.meta.(col_idx) in
              for i = 0 to !db.core.n_meta - 1 do
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
      (* First, we compute the indices of the columns to be kept.
         We keep the same column order as in the original matrix *)
      let idxs = ref IntSet.empty in
      Array.iteri
        (fun col_idx col_name ->
          if StringSet.mem col_name selected |> not then
            idxs := IntSet.add col_idx !idxs)
        db.core.idx_to_col_names;
      let idxs = IntSet.elements_array !idxs in
      let n = Array.length idxs in
      let filter_array a =
        Array.init n
          (fun i -> a.(idxs.(i))) in
      let core =
        { db.core with
          n_cols = n;
          idx_to_col_names = filter_array db.core.idx_to_col_names;
          meta = filter_array db.core.meta;
          storage = filter_array db.core.storage } in
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
          transform = Transformation.of_parameters { which = "power"; threshold = 1.; power = 1. };
          print_zero_rows = false;
          filter_columns = StringSet.empty;
          precision = 15
        }
      end
    let make_filename_table = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopCounter.txt"
    let to_table
        ?(filter = TableFilter.default) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) db prefix =
      let transform = Transformation.compute ~which:filter.transform
      and stats = stats_table_of_core_db ~threads ~verbose filter.transform db.core
      and fname = make_filename_table prefix in
      let output = open_out fname and meta = ref [] and rows = ref [] and cols = ref [] in
      (* We determine which rows and colunms should be output after all filters have been applied *)
      (*  Rows: metadata and k-mers *)
      if filter.print_metadata then
        Array.iteri
          (fun i meta_name ->
            (* There might be additional storage *)
            if i < db.core.n_meta then
              List.accum meta (meta_name, i))
          db.core.idx_to_meta_names;
      Array.iteri
        (fun i row_name ->
          (* There might be additional storage.
             Also, we only print the row if it has non-zero elements or if we are explicitly requested to do so *)
          if i < db.core.n_rows && (stats.row_stats.(i).sum > 0. || filter.print_zero_rows) then
            List.accum rows (row_name, i))
        db.core.idx_to_row_names;
      (* Columns *)
      Array.iteri
        (fun i col_name ->
          (* There might be additional storage, or we might need to remove some columns *)
          if i < db.core.n_cols && StringSet.mem col_name filter.filter_columns |> not then
            List.accum cols (col_name, i))
        db.core.idx_to_col_names;
      let meta = Array.of_rlist !meta and rows = Array.of_rlist !rows and cols = Array.of_rlist !cols in
      let n_rows = Array.length rows and n_cols = Array.length cols in
      (* There must be at least one row to print.
         If there are no columns, we just print metadata/row names  *)
      if (Array.length meta + n_rows) > 0 then begin
        if filter.transpose then begin
          if filter.print_col_names then begin
            (* We print row names - there is always at least one row *)
            let first_done = ref false in
            Array.iter
              (fun (meta_name, _) ->
                Printf.fprintf output "%s%s" (if !first_done || filter.print_row_names then "\t" else "") meta_name;
                first_done := true)
              meta;
            Array.iter
              (fun (row_name, _) ->
                Printf.fprintf output "%s%s" (if !first_done || filter.print_row_names then "\t" else "") row_name;
                first_done := true)
              rows;
            Printf.fprintf output "\n"
          end;
          Printf.fprintf output "%!";
          let processed_cols = ref 0 and buf = Buffer.create 1048576 in
          Processes.Parallel.process_stream_chunkwise
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
                  Printf.bprintf buf "%s" col_name;
                let first_done = ref false in
                Array.iter
                  (fun (_, meta_idx) ->
                    Printf.bprintf buf "%s%s"
                      (if !first_done || filter.print_row_names then "\t" else "") db.core.meta.(col_idx).(meta_idx);
                    first_done := true)
                  meta;
                Array.iter
                  (fun (_, row_idx) ->
                    Printf.bprintf buf "%s%.*g"
                      (if !first_done || filter.print_row_names then "\t" else "") filter.precision begin
                        IBAVector.N.to_int db.core.storage.(col_idx).@(row_idx) |>
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
                Printf.eprintf "%s\r(%s): Writing table to file '%s': done %d/%d lines%!"
                  String.TermIO.clear __FUNCTION__ fname !processed_cols n_cols)
            threads
        end else begin
          if filter.print_col_names then begin
            (* We print column names *)
            Array.iteri
              (fun i (col_name, _) ->
                Printf.fprintf output "%s%s" (if i > 0 || filter.print_row_names then "\t" else "") col_name)
              cols;
            Printf.fprintf output "\n"
          end;
          (* We print metadata as lines *)
          Array.iter
            (fun (meta_name, meta_idx) ->
              if filter.print_row_names then
                Printf.fprintf output "%s" meta_name;
              Array.iteri
                (fun i (_, col_idx) ->
                  Printf.fprintf output "%s%s"
                    (if i > 0 || filter.print_row_names then "\t" else "") db.core.meta.(col_idx).(meta_idx))
                cols;
              Printf.fprintf output "\n")
            meta;
          Printf.fprintf output "%!";
          let processed_rows = ref 0 and buf = Buffer.create 1048576 in
          Processes.Parallel.process_stream_chunkwise
            (fun () ->
              if !processed_rows < n_rows then
                let to_do = max 1 (elements_per_step / (max 1 n_cols)) |> min (n_rows - !processed_rows) in
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
                  Printf.bprintf buf "%s" row_name;
                Array.iteri
                  (fun i (_, col_idx) ->
                    Printf.bprintf buf "%s%.*g"
                      (if i > 0 || filter.print_row_names then "\t" else "") filter.precision begin
                        IBAVector.N.to_int db.core.storage.(col_idx).@(row_idx) |>
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
                Printf.eprintf "%s\r(%s): Writing table to file '%s': done %d/%d lines%!"
                  String.TermIO.clear __FUNCTION__ fname !processed_rows n_rows)
            threads
        end;
        let n_done =
          if filter.transpose then
            n_cols
          else
            n_rows in
        if verbose then
          Printf.eprintf "%s\r(%s): Writing table to file '%s': done %d/%d lines.\n%!"
            String.TermIO.clear __FUNCTION__ fname n_done n_done
      end;
      close_out output
    let to_spectra
        ?(filter = TableFilter.default) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) db prefix =
      let transform = Transformation.compute ~which:filter.transform
      and stats = stats_table_of_core_db ~threads ~verbose filter.transform db.core
      and fname = Spectra.make_filename prefix in
      let output = open_out fname and rows = ref [] and cols = ref [] in
      (* We determine which rows and colunms should be output after all filters have been applied *)
      (*  Rows: k-mers *)
      Array.iteri
        (fun i row_name ->
          (* There might be additional storage.
            Also, we only print the row if it has non-zero elements or if we are explicitly requested to do so *)
          if i < db.core.n_rows && (stats.row_stats.(i).sum > 0. || filter.print_zero_rows) then
            List.accum rows (row_name, i))
        db.core.idx_to_row_names;
      (* Columns *)
      Array.iteri
        (fun i col_name ->
          (* There might be additional storage, or we might need to remove some columns *)
          if i < db.core.n_cols && StringSet.mem col_name filter.filter_columns |> not then
            List.accum cols (col_name, i))
        db.core.idx_to_col_names;
      let rows = Array.of_rlist !rows and cols = Array.of_rlist !cols in
      let n_rows = Array.length rows and n_cols = Array.length cols in
      (* There must be at least one column to print.
         If there are no rows, we just print column names  *)
      if n_cols > 0 then begin
        let processed_cols = ref 0 and buf = Buffer.create 1048576 in
        Processes.Parallel.process_stream_chunkwise
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
              Printf.bprintf buf "\t%s\n" col_name;
              Array.iter
                (fun (row_name, row_idx) ->
                  let value =
                    IBAVector.N.to_int db.core.storage.(col_idx).@(row_idx) |>
                      transform
                        ~col_num:db.core.n_cols ~col_stats:stats.col_stats.(col_idx)
                        ~row_num:db.core.n_rows ~row_stats:stats.row_stats.(row_idx) in
                  if value > 0. then
                    Printf.bprintf buf "%s\t%.*g\n" row_name filter.precision value)
                rows
            done;
            hi_col - lo_col + 1, Buffer.contents buf)
          (fun (n_processed, block) ->
            Printf.fprintf output "%s" block;
            let old_processed_cols = !processed_cols in
            processed_cols := !processed_cols + n_processed;
            if verbose && !processed_cols / 2 > old_processed_cols / 2 then (* We write one column at the time *)
              Printf.eprintf "%s\r(%s): Writing spectra to file '%s': done %d/%d spectra%!"
                String.TermIO.clear __FUNCTION__ fname !processed_cols n_cols)
          threads;
        if verbose then
          Printf.eprintf "%s\r(%s): Writing spectra to file '%s': done %d/%d spectra.\n%!"
            String.TermIO.clear __FUNCTION__ fname n_cols n_cols
      end;
      close_out output
    let to_distances
        ?(normalise = true) ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false)
        distance db selection_1 selection_2 prefix =
      let transf = Transformation.of_parameters { which = "power"; threshold = 1.; power = 1. } in
      let stats = stats_table_of_core_db ~threads ~verbose transf db.core
      and n_r = db.core.n_rows in (* Does not change *)
      let make_submatrix selection =
        let idxs = ref IntSet.empty in
        Array.iteri
          (fun col_idx col_name ->
            if StringSet.mem col_name selection then
              idxs := IntSet.add col_idx !idxs)
          db.core.idx_to_col_names;
        let idxs = IntSet.elements_array !idxs in
        let n_c = Array.length idxs in
        (* The distance matrix is computed rowwise, and k-mers are physically stored as rows in db:
            we need to transpose *)
        { Matrix.Base.col_names = Array.sub db.core.idx_to_row_names 0 n_r;
          row_names = Array.init n_c (fun i -> db.core.idx_to_col_names.(idxs.(i)));
          data =
            (* Here we just need to convert the counts to floats, and possibly normalise *)
            Array.init n_c
              (fun i ->
                let idx = idxs.(i) in
                let norm = stats.col_stats.(idx).sum in
                let norm =
                  if not normalise || norm = 0. then
                    1.
                  else
                    norm in
                Float.Array.init n_r
                  (fun j -> IBAVector.N.to_float db.core.storage.(idx).@(j) /. norm)) } in
      let metric = Float.Array.make n_r 1.
      and matrix_1 = make_submatrix selection_1 and matrix_2 = make_submatrix selection_2 in
      Matrix.to_binary ~verbose {
        which = DMatrix;
        matrix =
          Matrix.Base.get_distance_rowwise ~threads ~elements_per_step ~verbose distance metric matrix_1 matrix_2
      } prefix
  end: sig
    (* Conceptually, each k-mer spectrum is stored as a column, even though in practice we store the transposed matrix *)
    type marshalled_t = {
      n_cols: int; (* The number of spectra *)
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
    module Transformation:
      sig
        type statistics_t = {
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
        (* Transformation function *)
        type t =
          | Binary of float
          | Power of float * float
          | CLR of float * float
          | Pseudo of float * float
        exception Invalid_transformation of string * float * float
        (* Only used internally *)
        val [@warning "-32"] compute:
          which:t -> col_num:int -> col_stats:statistics_t -> row_num:int -> row_stats:statistics_t -> int -> float
        type parameters_t = {
          which: string;
          (* A value in the interval (0.,1.) is taken as a relative threshold *)
          threshold: float;
          power: float
        }
        exception Unknown_transformation of string
        val of_parameters: parameters_t -> t
        val to_parameters: t -> parameters_t
      end
    type t = {
      core: marshalled_t;
      (* Inverted hashes for parsing *)
      col_names_to_idx: int StringHashtbl.t; (* Labels *)
      row_names_to_idx: int StringHashtbl.t; (* Hashes *)
      meta_names_to_idx: int StringHashtbl.t (* Metadata fields *)
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
    exception Unknown_combination_criterion of string
    module CombinationCriterion:
      sig
        type t = RescaledMean | RescaledMedian
        val of_string: string -> t
        val to_string: t -> string
      end
    (* Generate according to the specified criterion a combination of the spectra having the given labels,
        name the combination as directed, and add it to the database *)
    val add_combined_selected: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                               t -> string -> StringSet.t -> CombinationCriterion.t -> t
    (* Remove spectra with the given labels *)
    val remove_selected: t -> StringSet.t -> t
    (* Output information about the contents *)
    val output_summary: ?verbose:bool -> t -> unit
    (* Binary marshalling of the database *)
    exception Incompatible_archive_version of string * string
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t
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
    (* Spectral output *)
    val to_spectra:
      ?filter:TableFilter.t -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    (* Spectral distance matrix *)
    val to_distances:
      ?normalise:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
      Space.Distance.t -> t -> StringSet.t -> StringSet.t -> string -> unit
  end
)

