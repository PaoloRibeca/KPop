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

(* General matrix class *)
module Base:
  sig
    type t = {
      (* We number rows and columns starting from 0 *)
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      (* Stored row-wise *)
      storage: Float.Array.t array
    }
    val empty: t
    (* We read in a matrix which has conditions as row names
        and a (large) number of tags (genes, k-mers, etc.) as column names.
       In keeping with the convention accepted by R, the first row would be a header,
        and the first column the row names.
       Names might be quoted *)
    exception Quotes_in_name of string
    val strip_external_quotes_and_check: string -> string
    exception Wrong_number_of_columns of int * int * int
    val of_file: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    val to_file: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val transpose_single_threaded: ?verbose:bool -> t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    (* Compute row normalisations *)
    val get_normalizations: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                            Space.Distance.t -> Float.Array.t -> t -> Float.Array.t
    exception Incompatible_geometries of string array * string array
    exception Duplicate_label of string
    val merge_rowwise: ?verbose:bool -> t -> t -> t
    val multiply_matrix_vector:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_vector_single_threaded: ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    type sparse_vector_t = {
      length: int;
      elements: float Tools.IntMap.t
    }
    val multiply_matrix_sparse_vector_single_threaded: ?verbose:bool -> t -> sparse_vector_t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t -> t
    val get_distance_matrix: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                             Space.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                              Space.Distance.t -> Float.Array.t -> t -> t -> t
    module NearestNeighbor:
      sig
        type tt = t
        module Preconditioner:
          sig
            type t

          end
        type t
        val make_rowwise: Space.Distance.t -> Float.Array.t -> Preconditioner.t -> tt -> t
        val find_and_replace:
          ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
          t -> ((string * Float.Array.t * string * Float.Array.t * float) option -> (string * Float.Array.t) option) ->
          unit
      end
  end
= struct
    type t = {
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      storage: Float.Array.t array
    }
    let empty =
      { idx_to_col_names = [||]; idx_to_row_names = [||]; storage = [||] }
    let to_file ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) m fname =
      let n_cols = Array.length m.idx_to_col_names and n_rows = Array.length m.idx_to_row_names
      and output = open_out fname in
      if n_rows > 0 && n_cols > 0 then begin
        (* We output column names *)
        Printf.fprintf output "";
        Array.iter
          (fun name ->
            Printf.fprintf output "\t%s" name)
          m.idx_to_col_names;
        Printf.fprintf output "\n%!";
        let rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0
        and buf = Buffer.create 1048576 in
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
            Buffer.clear buf;
            (* We output rows *)
            for i = lo_row to hi_row do
              (* We output the row name *)
              m.idx_to_row_names.(i) |> Printf.bprintf buf "%s";
              Float.Array.iter (Printf.bprintf buf "\t%.*g" precision) m.storage.(i);
              Printf.bprintf buf "\n"
            done;
            hi_row - lo_row + 1, Buffer.contents buf)
          (fun (n_processed, block) ->
            Printf.fprintf output "%s" block;
            let new_processed_rows = !processed_rows + n_processed in
            if verbose && new_processed_rows / rows_per_step > !processed_rows / rows_per_step then
              Printf.eprintf "%s\r(%s): Writing table to file '%s': done %d/%d rows%!"
                Tools.String.TermIO.clear __FUNCTION__ fname new_processed_rows n_rows;
            processed_rows := new_processed_rows)
          threads
      end;
      if verbose then
        Printf.eprintf "%s\r(%s): Writing table to file '%s': done %d/%d rows.\n%!"
          Tools.String.TermIO.clear __FUNCTION__ fname n_rows n_rows;
      close_out output
    exception Quotes_in_name of string
    let re_quote = Str.regexp "\""
    let strip_external_quotes_and_check s =
      match String.length s, s with
      | 0, _ -> ""
      | 1, "\"" -> Quotes_in_name s |> raise
      | l, _ ->
        let s =
          if s.[0] = '"' && s.[l - 1] = '"' then
            String.sub s 1 (l - 2)
          else
            s in
        try
          Str.search_forward re_quote s 0 |> ignore;
          Quotes_in_name s |> raise
        with Not_found ->
          s
    exception Wrong_number_of_columns of int * int * int
    let of_file ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) filename =
      let input = open_in filename and line_num = ref 0
      and idx_to_col_names = ref [||] and idx_to_row_names = ref [] and storage = ref [] in
      begin try
        (* We process the header *)
        let line = input_line input |> Tools.Split.on_char_as_array '\t' in
        incr line_num;
        let l = Array.length line in
        (* We assume the matrix always to have row names, and ignore the first name in the header if present *)
        let num_cols = l - 1 in
        idx_to_col_names := Array.make num_cols "";
        Array.iteri
          (fun i name ->
            if i > 0 then
              !idx_to_col_names.(i - 1) <- strip_external_quotes_and_check name)
          line;
        (* We process the rest of the lines in parallel. The first element will be the name *)
        let end_reached = ref false and elts_read = ref 0 in
        Processes.Parallel.process_stream_chunkwise
          (fun () ->
            if !end_reached then
              raise End_of_file;
            let res = ref [] in
            begin try
              let cntr = ref 0 in
              while !cntr < bytes_per_step do
                let line = input_line input in
                incr line_num;
                Tools.List.accum res (!line_num, line);
                cntr := !cntr + String.length line
              done
            with End_of_file ->
              end_reached := true;
              if !res = [] then
                raise End_of_file
            end;
            List.rev !res)
          (List.map
            (fun (line_num, line) ->
              (* We decorate the line number with the results of parsing the line *)
              let line = Tools.Split.on_char_as_array '\t' line in
              let l = Array.length line in
              if l <> num_cols + 1 then
                Wrong_number_of_columns (line_num, l, num_cols + 1) |> raise;
              let array = Float.Array.create num_cols in
              Array.iteri
                (fun i el ->
                  if i > 0 then
                    (* The first element is the name *)
                    float_of_string el |> Float.Array.set array (i - 1))
                line;
              line_num, strip_external_quotes_and_check line.(0), array))
          (List.iter
            (fun (obs_line_num, name, numbers) ->
              incr line_num;
              assert (obs_line_num = !line_num);
              (* Only here do we actually fill out the memory for the result *)
              Tools.List.accum idx_to_row_names name;
              Tools.List.accum storage numbers;
              let new_elts_read = !elts_read + num_cols in
              if verbose && new_elts_read / 100000 > !elts_read / 100000 then
                Printf.eprintf "%s\r(%s): On line %d of file '%s': Read %d elements%!"
                  Tools.String.TermIO.clear __FUNCTION__ !line_num filename new_elts_read;
              elts_read := new_elts_read))
          threads;
        close_in input;
        if verbose then
          Printf.eprintf "%s\r(%s): On line %d of file '%s': Read %d elements.\n%!"
            Tools.String.TermIO.clear __FUNCTION__ !line_num filename !elts_read
      with End_of_file ->
        (* Empty file *)
        close_in input
      end;
      { idx_to_col_names = !idx_to_col_names;
        idx_to_row_names = Tools.Array.of_rlist !idx_to_row_names;
        storage = Tools.Array.of_rlist !storage }
    let [@warning "-27"] transpose_single_threaded ?(verbose = false) m =
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_col_names;
        storage =
          Array.init (Array.length m.idx_to_col_names)
            (fun old_col ->
              Float.Array.init (Array.length m.idx_to_row_names)
                (fun old_row -> Float.Array.get m.storage.(old_row) old_col)) }
    let transpose ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m =
      let n_rows = Array.length m.idx_to_col_names and n_cols = Array.length m.idx_to_row_names in
      let storage = Array.init n_rows (fun _ -> Float.Array.create 0)
      and rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0 in
      (* Generate points to be computed by the parallel processs *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_rows then (* The original columns *)
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
            (* The new row is the original column *)
            Float.Array.init n_cols (fun col -> Float.Array.get m.storage.(col) i) |> Tools.List.accum res
          done;
          lo_row, !res)
        (fun (lo_row, rows) ->
          List.iteri
            (fun offs_i row_i ->
              storage.(lo_row + offs_i) <- row_i;
              if verbose && !processed_rows mod rows_per_step = 0 then
                Printf.eprintf "%s\r(%s): Done %d/%d rows%!"
                  Tools.String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
              incr processed_rows)
            rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d rows.\n%!" Tools.String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_col_names;
        storage = storage }
    exception Incompatible_geometries of string array * string array
    exception Duplicate_label of string
    let merge_rowwise ?(verbose = false) m1 m2 =
      let merged_idx_to_col_names =
        if m1 = empty then
          m2.idx_to_col_names
        else m1.idx_to_col_names in
      if merged_idx_to_col_names <> m2.idx_to_col_names then
        Incompatible_geometries (m1.idx_to_col_names, m2.idx_to_col_names) |> raise;
      if verbose then
        Printf.eprintf "(%s): Merging matrices (%d+%d rows)...%!"
          __FUNCTION__ (Array.length m1.idx_to_row_names) (Array.length m2.idx_to_row_names);
      let merged_rows = ref Tools.StringMap.empty in
      Array.iteri
        (fun i name ->
          (* There ought to be no repeated names here *)
          merged_rows := Tools.StringMap.add name m1.storage.(i) !merged_rows)
        m1.idx_to_row_names;
      Array.iteri
        (fun i name ->
          match Tools.StringMap.find_opt name !merged_rows with
          | Some _ ->
            Duplicate_label name |> raise
          | None ->
            merged_rows := Tools.StringMap.add name m2.storage.(i) !merged_rows)
        m2.idx_to_row_names;
      let row_num = Tools.StringMap.cardinal !merged_rows in
      let merged_storage = Array.init row_num (fun _ -> Float.Array.create 0)
      and merged_idx_to_row_names = Array.make row_num "" in
      Tools.StringMap.iteri
        (fun i name arr ->
          merged_storage.(i) <- arr;
          merged_idx_to_row_names.(i) <- name)
        !merged_rows;
      if verbose then
        Printf.eprintf " done.\n%!";
      { idx_to_col_names = merged_idx_to_col_names;
        idx_to_row_names = merged_idx_to_row_names;
        storage = merged_storage }
    let multiply_matrix_vector_single_threaded ?(verbose = false) m v =
      if Array.length m.idx_to_col_names <> Float.Array.length v then
        Incompatible_geometries (m.idx_to_col_names, Array.make (Float.Array.length v) "") |> raise;
      let d = Array.length m.idx_to_row_names in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let res = Float.Array.create d and elts_done = ref 0 in
      (* We decorate each vector element coordinate with the respective value *)
      Array.iteri
        (fun i row ->
          let acc = ref 0. in
          Float.Array.iter2
            (fun el_1 el_2 ->
              acc := !acc +. (el_1 *. el_2))
            row v;
          Float.Array.set res i !acc;
          incr elts_done;
          if verbose && !elts_done mod 100 = 0 then
            Printf.eprintf "%s\r(%s): Done %d/%d elements%!" Tools.String.TermIO.clear __FUNCTION__ !elts_done d)
        m.storage;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements.\n%!" Tools.String.TermIO.clear __FUNCTION__ !elts_done d;
      res
    type sparse_vector_t = {
      length: int;
      elements: float Tools.IntMap.t
    }
    let multiply_matrix_sparse_vector_single_threaded ?(verbose = false) m s_v =
      if Array.length m.idx_to_col_names <> s_v.length then
        Incompatible_geometries (m.idx_to_col_names, Array.make (s_v.length) "") |> raise;
      let d = Array.length m.idx_to_row_names in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let res = Float.Array.make d 0. and elts_done = ref 0 in
      (* We decorate each vector element coordinate with the respective value *)
      Array.iteri
        (fun i row ->
          let acc = ref 0. in
          Tools.IntMap.iter
            (fun j el ->
              acc := !acc +. (Float.Array.get row j *. el))
            s_v.elements;
          Float.Array.set res i !acc;
          incr elts_done;
          if verbose && !elts_done mod 100 = 0 then
            Printf.eprintf "%s\r(%s): Done %d/%d elements%!" Tools.String.TermIO.clear __FUNCTION__ !elts_done d)
        m.storage;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements.\n%!" Tools.String.TermIO.clear __FUNCTION__ !elts_done d;
      res
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m v =
      let n_rows = Array.length m.idx_to_row_names and n_cols = Array.length m.idx_to_col_names in
      if n_cols <> Float.Array.length v then
        Incompatible_geometries (m.idx_to_col_names, Array.make (Float.Array.length v) "") |> raise;
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
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
            (* We decorate each vector element coordinate with the respective value *)
            let acc = ref 0. in
            Float.Array.iter2
              (fun el_1 el_2 ->
                acc := !acc +. (el_1 *. el_2))
              m.storage.(i) v;
            Tools.List.accum res !acc
          done;
          lo_row, !res)
        (fun (lo_row, v) ->
          List.iteri
            (fun offs_i el ->
              (* Only here do we actually fill out the memory for the result *)
              Float.Array.set res (lo_row + offs_i) el;
              if verbose && !processed_rows mod rows_per_step = 0 then
                Printf.eprintf "%s\r(%s): Done %d/%d rows%!"
                  Tools.String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
              incr processed_rows)
            v)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d rows.\n%!" Tools.String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
      res
    let multiply_matrix_matrix ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m1 m2 =
      if m1.idx_to_col_names <> m2.idx_to_row_names then
        Incompatible_geometries (m1.idx_to_col_names, m2.idx_to_row_names) |> raise;
      let row_num = Array.length m1.idx_to_row_names and col_num = Array.length m2.idx_to_col_names in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let storage = Array.init row_num (fun _ -> Float.Array.create col_num) in
      (* Generate points to be computed by the parallel processs *)
      let prod = row_num * col_num in
      let i = ref 0 and j = ref 0 and elts_done = ref 0 and end_reached = ref (prod = 0) in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file;
          let res = ref [] in
          begin try
            let cntr = ref 0 in
            while !cntr < elements_per_step do
              Tools.List.accum res (!i, !j);
              incr j;
              if !j = col_num then begin
                incr i;
                if !i = row_num then begin
                  end_reached := true;
                  raise Exit
                end;
                j := 0
              end;
              incr cntr
            done
          with Exit -> ()
          end;
          List.rev !res)
        (List.map
          (* We decorate each matrix element coordinate with the respective value *)
          (fun (i, j) ->
            let acc = ref 0. in
            Float.Array.iteri
              (fun k el ->
                acc := !acc +. (el *. Float.Array.get m2.storage.(k) j))
              m1.storage.(i);
            i, j, !acc))
        (List.iter
          (fun (i, j, el) ->
            (* Only here do we actually fill out the memory for the result *)
            Float.Array.set storage.(i) j el;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "%s\r(%s): Done %d/%d elements%!" Tools.String.TermIO.clear __FUNCTION__ !elts_done prod;
            incr elts_done))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements.\n%!" Tools.String.TermIO.clear __FUNCTION__ !elts_done prod;
      { idx_to_col_names = m2.idx_to_col_names;
        idx_to_row_names = m1.idx_to_row_names;
        storage = storage }
    (* Compute normalisations for rows *)
    let get_normalizations ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) distance metric m =
      let n_rows = Array.length m.idx_to_row_names and n_cols = Array.length m.idx_to_col_names in
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
            Space.Distance.compute_norm distance metric m.storage.(i) |> Tools.List.accum res
          done;
          lo_row, !res)
        (fun (lo_row, norms) ->
          List.iteri
            (fun offs_i norm_i ->
              Float.Array.set res (lo_row + offs_i) (if norm_i = 0. then 1. else norm_i);
              if verbose && !processed_rows mod elements_per_step = 0 then
                Printf.eprintf "%s\r(%s): Done %d/%d rows%!"
                  Tools.String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
              incr processed_rows)
            norms)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d rows.\n%!" Tools.String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
      res
    (* Compute rowwise distance *)
    let get_distance_matrix ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      let d = Array.length m.idx_to_row_names in
      (* We compute normalisations *)
      let norms =
        if normalize then
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m
        else
          Float.Array.make d 1. in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let storage = Array.init d (fun _ -> Float.Array.create d) in
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
              Tools.List.accum res (!i, !j);
              incr j;
              if !j > !i then begin
                incr i;
                if !i = d then begin
                  end_reached := true;
                  raise Exit
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
                ~adaptor_a:(fun a -> a /. Float.Array.get norms i) ~adaptor_b:(fun b -> b /. Float.Array.get norms j)
                distance metric m.storage.(i) m.storage.(j)
            end))
        (List.iter
          (fun (i, j, dist) ->
            (* Only here do we actually fill out the memory for the result *)
            Float.Array.set storage.(i) j dist;
            (* We symmetrise the matrix *)
            Float.Array.set storage.(j) i dist;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "%s\r(%s): Done %d/%d elements%!"
                Tools.String.TermIO.clear __FUNCTION__ !elts_done total;
            incr elts_done))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements.\n%!" Tools.String.TermIO.clear __FUNCTION__ !elts_done total;
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_row_names;
        storage = storage }
    let get_distance_rowwise ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m1 m2 =
      if m1.idx_to_col_names <> m2.idx_to_col_names then
        Incompatible_geometries (m1.idx_to_col_names, m2.idx_to_col_names) |> raise;
      let r1 = Array.length m1.idx_to_row_names and r2 = Array.length m2.idx_to_row_names in
      (* We compute normalisations *)
      let n1, n2 =
        if normalize then
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m1,
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m2
        else
          Float.Array.make r1 1., Float.Array.make r2 1. in
      (*
      to_file ~verbose (transpose ~verbose {
        idx_to_col_names = m1.idx_to_row_names;
        idx_to_row_names = [| "Normalizations" |];
        storage = [| n1 |]
      }) "N1.txt";
      to_file ~verbose (transpose ~verbose {
        idx_to_col_names = m2.idx_to_row_names;
        idx_to_row_names = [| "Normalizations" |];
        storage = [| n2 |]
      }) "N2.txt";
      *)
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let storage = Array.init r2 (fun _ -> Float.Array.create r1) in
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
              Tools.List.accum res (!i, !j);
              incr j;
              if !j = r2 then begin
                incr i;
                if !i = r1 then begin
                  end_reached := true;
                  raise Exit
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
                ~adaptor_a:(fun a -> a /. Float.Array.get n1 i) ~adaptor_b:(fun b -> b /. Float.Array.get n2 j)
                distance metric m1.storage.(i) m2.storage.(j)
            end))
        (List.iter
          (fun (i, j, dist) ->
            Float.Array.set storage.(j) i dist;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "%s\r(%s): Done %d/%d elements=%.3g%%%!"
                Tools.String.TermIO.clear __FUNCTION__
                !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod);
            incr elts_done))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements=%.3g%%.\n%!"
          Tools.String.TermIO.clear __FUNCTION__
          !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod);
      { idx_to_col_names = m1.idx_to_row_names;
        idx_to_row_names = m2.idx_to_row_names;
        storage = storage }
    module NearestNeighbor =
      struct
        (* These are column statistics *)
        module Statistics =
          struct
            type t = {
              min: float;
              max: float;
              sum: float
            }
            let empty = {
              min = 0.;
              max = 0.;
              sum = 0.
            }
            let [@warning "-32"] compute ?(threads = 1) ?(verbose = false) m =
              (* Column n *)
              let compute_one n =
                let min = ref 0. and max = ref 0. and sum = ref 0. in
                Array.iter
                  (fun row ->
                    let v = Float.Array.get row n in
                    min := Stdlib.min !min v;
                    max := Stdlib.max !max v;
                    sum := !sum +. v)
                  m.storage;
                { min = !min;
                  max = !max;
                  sum = !sum } in
              let n = Array.length m.idx_to_col_names in
              let step = n / threads / 5 |> max 1 and processed = ref 0 and res = Array.make n empty in
              Processes.Parallel.process_stream_chunkwise
                (fun () ->
                  if verbose then
                    Printf.eprintf "\rComputing column statistics [%d/%d]%!" !processed n;
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
                    processed + i |> compute_one |> Tools.List.accum res
                  done;
                  processed, List.rev !res)
                (fun (base, stats) ->
                  List.iteri
                    (fun i s ->
                      res.(base + i) <- s)
                    stats)
                threads;
              res
          end
        module Preconditioner =
          struct
            type t

          end
        type tt = t
        module StringMap = Tools.StringMap
        type t = {
          vectors: Float.Array.t StringMap.t;


        }
        let [@warning "-27"] make_rowwise dist metr precond rowwise =
          (* Input has sample names as rows, dimensions as columns.
             Output is one specific dimension selected according to the preconditioner,
              and an interator built upon it *)

          { vectors = StringMap.empty }

        let [@warning "-27"] find_and_replace ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) nn f =


          ()

      end
  end

(* KPop-specialised matrices.
   We include in order not to have a repeated module prefix *)
include [@warning "-32"] (
  struct
    module Type =
      struct
        type t =
          | Twister
          | Inertia
          | Twisted
          | DMatrix
          | Metrics
        let to_string = function
          | Twister -> "KPopTwister"
          | Inertia -> "KPopInertia"
          | Twisted -> "KPopTwisted"
          | DMatrix -> "KPopDMatrix"
          | Metrics -> "KPopMetrics"
        exception Unknown_type of string
        let of_string = function
          | "KPopTwister" -> Twister
          | "KPopInertia" -> Inertia
          | "KPopTwisted" -> Twisted
          | "KPopDMatrix" -> DMatrix
          | "KPopMetrics" -> Metrics
          | w ->
            Unknown_type w |> raise
      end
    type t = {
      which: Type.t;
      matrix: Base.t
    }
    let empty which =
      { which; matrix = Base.empty }
    (* The three following functions implement automatic file naming *)
    let make_filename_table which = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ "." ^ Type.to_string which ^ ".txt"
    let make_filename_binary which name =
      match which, name with
      | _, _ when String.length name >= 5 && String.sub name 0 5 = "/dev/" -> name
      | Type.Twisted, prefix | DMatrix, prefix | Metrics, prefix -> prefix ^ "." ^ Type.to_string which
      | Twister, _ | Inertia, _ -> assert false (* Should always be done through KPopTwister *)
    let make_filename_summary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopSummary.txt"
    (* We redefine the implementation for Matrix in order to set the correct KPop types *)
    let of_file ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) which prefix =
      { which; matrix = make_filename_table which prefix |> Base.of_file ~threads ~bytes_per_step ~verbose }
    let to_file ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) m prefix =
      make_filename_table m.which prefix |> Base.to_file ~precision ~threads ~elements_per_step ~verbose m.matrix
    let transpose_single_threaded ?(verbose = false) m =
      { m with matrix = Base.transpose_single_threaded ~verbose m.matrix }
    let transpose ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m =
      { m with matrix = Base.transpose ~threads ~elements_per_step ~verbose m.matrix }
    exception Incompatible_matrices of Type.t * Type.t
    let merge_rowwise ?(verbose = false) m1 m2 =
      if m1.which <> m2.which then
        Incompatible_matrices (m1.which, m2.which) |> raise;
      { which = m1.which; matrix = Base.merge_rowwise ~verbose m1.matrix m2.matrix }
    let multiply_matrix_vector_single_threaded ?(verbose = false) m =
      Base.multiply_matrix_vector_single_threaded ~verbose m.matrix
    let multiply_matrix_sparse_vector_single_threaded ?(verbose = false) m =
      Base.multiply_matrix_sparse_vector_single_threaded ~verbose m.matrix
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m v =
      Base.multiply_matrix_vector ~threads ~elements_per_step ~verbose m.matrix v
    let multiply_matrix_matrix ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) which m1 m2 =
      { which; matrix = Base.multiply_matrix_matrix ~threads ~elements_per_step ~verbose m1.matrix m2.matrix }
    let get_distance_matrix ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      { which = DMatrix;
        matrix = Base.get_distance_matrix ~normalize ~threads ~elements_per_step ~verbose distance metric m.matrix }
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    let get_distance_rowwise ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m1 m2 =
      { which = DMatrix;
        matrix =
          Base.get_distance_rowwise ~normalize ~threads ~elements_per_step ~verbose
            distance metric m1.matrix m2.matrix }
    module FloatIntMultimap = Tools.Multimap (Tools.ComparableFloat) (Tools.ComparableInt)
    let summarize_distance_matrix_row req_len row_name row col_names buf =
      let n_cols = Float.Array.length row in
      let f_n_cols = float_of_int n_cols and distr = ref FloatIntMultimap.empty in
      (* We find the median and filter the result *)
      Float.Array.iteri
        (fun col_idx dist ->
          distr := FloatIntMultimap.add dist col_idx !distr)
        row;
      let eff_len = ref 0 and median_pos = n_cols / 2 |> ref and median = ref 0. and acc = ref 0. in
      FloatIntMultimap.iter_set
        (fun dist set ->
          let set_len = FloatIntMultimap.ValSet.cardinal set in
          acc := !acc +. (float_of_int set_len *. dist);
          if !median_pos >= 0 && !median_pos - set_len < 0 then
            median := dist;
          median_pos := !median_pos - set_len;
          if !eff_len < req_len then
            eff_len := !eff_len + set_len)
        !distr;
      let eff_len = !eff_len and median = !median and mean =
        if n_cols > 0 then
          !acc /. f_n_cols
        else
          0. in
      (* We compute standard deviation e MAD *)
      acc := 0.;
      let ddistr = ref Tools.FloatMap.empty in
      Float.Array.iteri
        (fun _ dist ->
          let d = dist -. mean in
          acc := !acc +. (d *. d);
          let d = (dist -. median) |> abs_float in
          ddistr :=
            match Tools.FloatMap.find_opt d !ddistr with
            | None ->
              Tools.FloatMap.add d 1 !ddistr
            | Some n ->
              Tools.FloatMap.add d (n + 1) !ddistr)
        row;
      median_pos := n_cols / 2;
      let mad = ref 0. in
      Tools.FloatMap.iter
        (fun d occs ->
          if !median_pos >= 0 && !median_pos - occs < 0 then
            mad := d;
          median_pos := !median_pos - occs)
        !ddistr;
      let mad = !mad and stddev =
        if n_cols > 1 then
          !acc /. (f_n_cols -. 1.) |> sqrt
        else
          0. in
      Printf.bprintf buf "%s\t%.15g\t%.15g\t%.15g\t%.15g" row_name mean stddev median mad;
      FloatIntMultimap.iteri
        (fun i dist col_idx ->
          if i < eff_len then
            Printf.bprintf buf "\t%s\t%.15g\t%.15g" col_names.(col_idx) dist ((dist -. mean) /. stddev))
        !distr;
      Printf.bprintf buf "\n"
    exception Unexpected_type of Type.t * Type.t
    let summarize_rowwise
        ?(normalize = true) ?(keep_at_most = Some 2) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m1 m2 prefix =
      if m1.which <> Twisted then
        Unexpected_type (m1.which, Twisted) |> raise;
      if m2.which <> Twisted then
        Unexpected_type (m2.which, Twisted) |> raise;
      if m1.matrix.idx_to_col_names <> m2.matrix.idx_to_col_names then
        Base.Incompatible_geometries (m1.matrix.idx_to_col_names, m2.matrix.idx_to_col_names) |> raise;
      let r1 = Array.length m1.matrix.idx_to_row_names and r2 = Array.length m2.matrix.idx_to_row_names in
      (* We compute normalisations *)
      let n1, n2 =
        if normalize then
          Base.get_normalizations ~threads ~elements_per_step ~verbose distance metric m1.matrix,
          Base.get_normalizations ~threads ~elements_per_step ~verbose distance metric m2.matrix
        else
          Float.Array.make r1 1., Float.Array.make r2 1. in
      (*
      Base.to_file ~verbose (Base.transpose ~verbose {
        idx_to_col_names = m1.matrix.idx_to_row_names;
        idx_to_row_names = [| "Normalizations" |];
        storage = [| n1 |]
      }) "N1.txt";
      Base.to_file ~verbose (Base.transpose ~verbose {
        idx_to_col_names = m2.matrix.idx_to_row_names;
        idx_to_row_names = [| "Normalizations" |];
        storage = [| n2 |]
      }) "N2.txt";
      *)
      let fname = make_filename_summary prefix in
      let output = open_out fname
      and n_cols = Array.length m1.matrix.idx_to_col_names in
      let req_len =
        match keep_at_most with
        | None -> r1
        | Some at_most -> at_most in
      let rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0
      and buf = Buffer.create 1048576 in
      (* Parallel section *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < r2 then
            let to_do = min rows_per_step (r2 - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          Buffer.clear buf;
          for j = lo_row to hi_row do
            (* For each row number of m2, we compute the respective distances from the rows of m1... *)
            let distances =
              Float.Array.init r1
                (fun i ->
                  Space.Distance.compute
                    ~adaptor_a:(fun a -> a /. Float.Array.get n1 i) ~adaptor_b:(fun b -> b /. Float.Array.get n2 j)
                    distance metric m1.matrix.storage.(i) m2.matrix.storage.(j)) in
            (* ...and summarise them *)
            summarize_distance_matrix_row
              req_len m2.matrix.idx_to_row_names.(j) distances m1.matrix.idx_to_row_names buf
          done;
          hi_row - lo_row + 1, Buffer.contents buf)
        (fun (n_processed, block) ->
          Printf.fprintf output "%s" block;
          let new_processed_rows = !processed_rows + n_processed in
          if verbose && new_processed_rows / rows_per_step > !processed_rows / rows_per_step then
            Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows%!"
              Tools.String.TermIO.clear __FUNCTION__ fname new_processed_rows r2;
          processed_rows := new_processed_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows.\n%!"
          Tools.String.TermIO.clear __FUNCTION__ fname r2 r2;
      close_out output
    let summarize_distance ?(keep_at_most = Some 2) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        m prefix =
      if m.which <> DMatrix then
        Unexpected_type (m.which, DMatrix) |> raise;
      let fname = make_filename_summary prefix in
      let output = open_out fname
      and n_cols = Array.length m.matrix.idx_to_col_names in
      let req_len =
        match keep_at_most with
        | None -> n_cols
        | Some at_most -> at_most
      and n_rows = Array.length m.matrix.idx_to_row_names in
      (* Parallel section *)
      let rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0
      and buf = Buffer.create 1048576 in
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
          Buffer.clear buf;
          for i = lo_row to hi_row do
            summarize_distance_matrix_row
              req_len m.matrix.idx_to_row_names.(i) m.matrix.storage.(i) m.matrix.idx_to_col_names buf
          done;
          hi_row - lo_row + 1, Buffer.contents buf)
        (fun (n_processed, block) ->
          Printf.fprintf output "%s" block;
          let new_processed_rows = !processed_rows + n_processed in
          if verbose && new_processed_rows / rows_per_step > !processed_rows / rows_per_step then
            Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows%!"
              Tools.String.TermIO.clear __FUNCTION__ fname new_processed_rows n_rows;
          processed_rows := new_processed_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows.\n%!"
          Tools.String.TermIO.clear __FUNCTION__ fname n_rows n_rows;
      close_out output
    (* *)
    let archive_version = "2022-04-03"
    (* *)
    let to_channel output m =
      Type.to_string m.which |> output_value output;
      archive_version |> output_value output;
      output_value output m.matrix
    let to_binary ?(verbose = false) m prefix =
      let fname = make_filename_binary m.which prefix in
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      to_channel output m;
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    exception Incompatible_archive_version of string * string
    let of_channel input =
      let which = (input_value input: string) in
      let version = (input_value input: string) in
      if version <> archive_version then
        Incompatible_archive_version (which, version) |> raise;
      { which = Type.of_string which; matrix = (input_value input: Base.t) }
    let of_binary ?(verbose = false) which prefix =
      let fname = make_filename_binary which prefix in
      let input = open_in fname in
      if verbose then
        Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let res = of_channel input in
      close_in input;
      if which <> res.which then
        Unexpected_type (which, res.which) |> raise;
      if verbose then
        Printf.eprintf " done.\n%!";
      res
  end: sig
    module Type:
      sig
        type t =
          | Twister
          | Inertia
          | Twisted
          | DMatrix
          | Metrics
        val to_string: t -> string
        exception Unknown_type of string
        val of_string: string -> t
      end
    type t = {
      which: Type.t;
      matrix: Base.t
    }
    val empty: Type.t -> t
    (* All file name arguments are in fact _prefixes_ *)
    val of_file: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> Type.t -> string -> t
    (* This one discards type information - use at your own risk *)
    val to_file: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val transpose_single_threaded: ?verbose:bool -> t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    exception Incompatible_matrices of Type.t * Type.t
    val merge_rowwise: ?verbose:bool -> t -> t -> t
    val multiply_matrix_vector_single_threaded: ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_sparse_vector_single_threaded: ?verbose:bool -> t -> Base.sparse_vector_t -> Float.Array.t
    val multiply_matrix_vector: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                                t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Type.t -> t -> t -> t
    (* Compute distances between the rows of a matrix *)
    val get_distance_matrix: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                             Space.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                              Space.Distance.t -> Float.Array.t -> t -> t -> t
    exception Unexpected_type of Type.t * Type.t
    val summarize_rowwise: ?normalize:bool -> ?keep_at_most:int option -> ?threads:int -> ?elements_per_step:int ->
                           ?verbose:bool ->
                           Space.Distance.t -> Float.Array.t -> t -> t -> string -> unit
    val summarize_distance: ?keep_at_most:int option -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                            t -> string -> unit
    (* Binary marshalling of the matrix *)
    val to_channel: out_channel -> t -> unit
    exception Incompatible_archive_version of string * string
    val of_channel: in_channel -> t
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> Type.t -> string -> t
  end
)

