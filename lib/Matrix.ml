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

module Misc =
  struct
    let _resize_array_ length make blit a idx =
      let l = length a and aug_idx = idx + 1 in
      if l < aug_idx then begin
        let res = make (max aug_idx (l * 14 / 10)) in
        blit a 0 res 0 l;
        res
      end else
        a
    let resize_array make = _resize_array_ Array.length make Array.blit
  end

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
    exception Wrong_number_of_columns of int * int * int
    val of_file: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    val to_file: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val transpose_single_threaded: ?verbose:bool -> t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
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
    val multiply_matrix_matrix:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t -> t
    val get_distance_matrix:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Space.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Space.Distance.t -> Float.Array.t -> t -> t -> t
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
        Printf.fprintf output "\"\"";
        Array.iter
          (fun name ->
            Printf.fprintf output "\t\"%s\"" name)
          m.idx_to_col_names;
        Printf.fprintf output "\n%!";
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
            (* We output rows *)
            for i = lo_row to hi_row do
              (* We output the row name *)
              m.idx_to_row_names.(i) |> Printf.bprintf buf "\"%s\"";
              Float.Array.iter (Printf.bprintf buf "\t%.*g" precision) m.storage.(i);
              Printf.bprintf buf "\n"
            done;
            hi_row - lo_row + 1, Buffer.contents buf)
          (fun (n_processed, block) ->
            Printf.fprintf output "%s" block;
            let old_processed_rows = !processed_rows in
            processed_rows := !processed_rows + n_processed;
            if verbose && !processed_rows / 10000 > old_processed_rows / 10000 then
              Printf.eprintf "\rWriting table to file '%s': done %d/%d lines%!"
                fname !processed_rows n_rows)
          threads
      end;
      if verbose then
        Printf.eprintf "\rWriting table to file '%s': done %d/%d lines.\n%!" fname n_rows n_rows;
      close_out output
    let strip_quotes s =
      let l = String.length s in
      if l = 0 then
        ""
      else
        let s =
          if s.[0] = '"' then
            String.sub s 1 (l - 1)
          else
            s in
        let l = String.length s in
        if l = 0 then
          ""
        else
          if s.[l - 1] = '"' then
            String.sub s 0 (l - 1)
          else
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
              !idx_to_col_names.(i - 1) <- strip_quotes name)
          line;
        (* We process the rest of the lines in parallel. The first element will be the name *)
        let end_reached = ref false and elts_read = ref 0 in
        Tools.Parallel.process_stream_chunkwise
          (fun () ->
            if !end_reached then
              raise End_of_file
            else begin
              let res = ref [] and cntr = ref 0 in
              begin try
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
              List.rev !res
            end)
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
              line_num, strip_quotes line.(0), array))
          (List.iter
            (fun (obs_line_num, name, numbers) ->
              incr line_num;
              assert (obs_line_num = !line_num);
              (* Only here do we actually fill out the memory for the result *)
              Tools.List.accum idx_to_row_names name;
              Tools.List.accum storage numbers;
              let new_elts_read = !elts_read + num_cols in
              if verbose && new_elts_read / 100000 > !elts_read / 100000 then
                Printf.eprintf "\r(%s): On line %d of file '%s': Read %d elements%!            \r"
                  __FUNCTION__ !line_num filename !elts_read;
              elts_read := new_elts_read))
          threads;
        close_in input;
        if verbose then
          Printf.eprintf "\r(%s): On line %d of file '%s': Read %d elements.            \n%!"
            __FUNCTION__ !line_num filename !elts_read
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
    let transpose ?(threads = 1) ?(elements_per_step = 1) ?(verbose = false) m =
      let row_num = Array.length m.idx_to_col_names and col_num = Array.length m.idx_to_row_names in
      let storage = Array.init row_num (fun _ -> Float.Array.create 0) in
      (* Generate points to be computed by the parallel processs *)
      let i = ref 0 and end_reached = ref false and elts_done = ref 0 in
      Tools.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file
          else begin
            let res = ref [] and cntr = ref 0 in
            begin try
              while !cntr < elements_per_step do
                Tools.List.accum res !i;
                incr i;
                if !i = row_num then begin (* The original columns *)
                  end_reached := true;
                  raise Exit
                end;
                incr cntr
              done
            with Exit -> ()
            end;
            List.rev !res
          end)
        (List.map
          (* The new row is the original column *)
          (fun i ->
            i, Float.Array.init col_num (fun col -> Float.Array.get m.storage.(col) i)))
        (List.iter
          (fun (i, row_i) ->
            storage.(i) <- row_i;
            incr elts_done;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "\r(%s): Done %d/%d rows%!            \r" __FUNCTION__ !elts_done row_num))
        threads;
      if verbose then
        Printf.eprintf "\r(%s): Done %d/%d rows.            \n%!" __FUNCTION__ !elts_done row_num;
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
            Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done d)
        m.storage;
      if verbose then
        Printf.eprintf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done d;
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
            Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done d)
        m.storage;
      if verbose then
        Printf.eprintf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done d;
      res
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m v =
      if Array.length m.idx_to_col_names <> Float.Array.length v then
        Incompatible_geometries (m.idx_to_col_names, Array.make (Float.Array.length v) "") |> raise;
      let d = Array.length m.idx_to_row_names in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let res = Float.Array.create d and i = ref 0 and end_reached = ref false and elts_done = ref 0 in
      Tools.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file
          else begin
            let res = ref [] and cntr = ref 0 in
            begin try
              while !cntr < elements_per_step do
                Tools.List.accum res !i;
                incr i;
                if !i = d then begin
                  end_reached := true;
                  raise Exit
                end;
                incr cntr
              done
            with Exit -> ()
            end;
            List.rev !res
          end)
        (List.map
          (* We decorate each vector element coordinate with the respective value *)
          (fun i ->
            let acc = ref 0. in
            Float.Array.iter2
              (fun el_1 el_2 ->
                acc := !acc +. (el_1 *. el_2))
              m.storage.(i) v;
            i, !acc))
        (List.iter
          (fun (i, el) ->
            (* Only here do we actually fill out the memory for the result *)
            Float.Array.set res i el;
            incr elts_done;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done d))
        threads;
      if verbose then
        Printf.eprintf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done d;
      res
    let multiply_matrix_matrix ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m1 m2 =
      if m1.idx_to_col_names <> m2.idx_to_row_names then
        Incompatible_geometries (m1.idx_to_col_names, m2.idx_to_row_names) |> raise;
      let row_num = Array.length m1.idx_to_row_names and col_num = Array.length m2.idx_to_col_names in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let storage = Array.init row_num (fun _ -> Float.Array.create col_num) in
      (* Generate points to be computed by the parallel processs *)
      let i = ref 0 and j = ref 0 and end_reached = ref false
      and prod = row_num * col_num and elts_done = ref 0 in
      Tools.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file
          else begin
            let res = ref [] and cntr = ref 0 in
            begin try
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
            List.rev !res
          end)
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
            incr elts_done;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done prod))
        threads;
      if verbose then
        Printf.eprintf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done prod;
      { idx_to_col_names = m2.idx_to_col_names;
        idx_to_row_names = m1.idx_to_row_names;
        storage = storage }
    (* Compute rowwise distance *)
    let get_distance_matrix ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) distance metric m =
      let distance = Space.Distance.compute distance metric and d = Array.length m.idx_to_row_names in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let storage = Array.init d (fun _ -> Float.Array.create d) in
      (* Generate points to be computed by the parallel processs *)
      let i = ref 0 and j = ref 0 and end_reached = ref false
      and total = (d * (d + 1)) / 2 and elts_done = ref 0 in
      Tools.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file
          else begin
            (* We only compute 1/2 of the matrix, and symmetrise it at the end of the computation *)
            let res = ref [] and cntr = ref 0 in
            begin try
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
            List.rev !res
          end)
        (List.map
          (* We decorate each matrix element coordinate with the respective distance *)
          (fun (i, j) ->
            i, j, distance m.storage.(i) m.storage.(j)))
        (List.iter
          (fun (i, j, dist) ->
            (* Only here do we actually fill out the memory for the result *)
            Float.Array.set storage.(i) j dist;
            incr elts_done;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done total))
        threads;
      if verbose then begin
        Printf.eprintf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done total;
        Printf.eprintf "(%s): Symmetrizing...%!" __FUNCTION__
      end;
      (* We symmetrise the matrix *)
      let red_d = d - 1 in
      for i = 0 to red_d do
        for j = i + 1 to red_d do
          Float.Array.get storage.(j) i |> Float.Array.set storage.(i) j
        done
      done;
      if verbose then
        Printf.eprintf " done.\n%!";
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_row_names;
        storage = storage }
    let get_distance_rowwise ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) distance metric m1 m2 =
      let distance = Space.Distance.compute distance metric in
      if m1.idx_to_col_names <> m2.idx_to_col_names then
        Incompatible_geometries (m1.idx_to_col_names, m2.idx_to_col_names) |> raise;
      let r1 = Array.length m1.idx_to_row_names and r2 = Array.length m2.idx_to_row_names in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let storage = Array.init r1 (fun _ -> Float.Array.create r2) in
      (* Generate points to be computed by the parallel processs *)
      let i = ref 0 and j = ref 0 and end_reached = ref false
      and prod = r1 * r2 and elts_done = ref 0 in
      Tools.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file
          else begin
            let res = ref [] and cntr = ref 0 in
            begin try
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
            List.rev !res
          end)
        (List.map
          (* We decorate each matrix element coordinate with the respective distance *)
          (fun (i, j) ->
            i, j, distance m1.storage.(i) m2.storage.(j)))
        (List.iter
          (fun (i, j, dist) ->
            Float.Array.set storage.(i) j dist;
            incr elts_done;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "\r(%s): Done %d/%d elements=%.3g%%%!            \r"
                __FUNCTION__ !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod)))
        threads;
      if verbose then
        Printf.eprintf "\r(%s): Done %d/%d elements=%.3g%%.            \n%!"
          __FUNCTION__ !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod);
      { idx_to_col_names = m2.idx_to_row_names;
        idx_to_row_names = m1.idx_to_row_names;
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
              Tools.Parallel.process_stream_chunkwise
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
                  let res = ref [] and red_to_do = to_do - 1 in
                  for i = 0 to red_to_do do
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

        let [@warning "-27"] find_and_replace ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) nn f =


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
    let transpose ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m =
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
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m v =
      Base.multiply_matrix_vector ~threads ~elements_per_step ~verbose m.matrix v
    let multiply_matrix_matrix ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) which m1 m2 =
      { which; matrix = Base.multiply_matrix_matrix ~threads ~elements_per_step ~verbose m1.matrix m2.matrix }
    let get_distance_matrix ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) distance metric m =
      { which = DMatrix;
        matrix = Base.get_distance_matrix ~threads ~elements_per_step ~verbose distance metric m.matrix }
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    let get_distance_rowwise ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) distance metric m1 m2 =
      { which = DMatrix;
        matrix = Base.get_distance_rowwise ~threads ~elements_per_step ~verbose distance metric m1.matrix m2.matrix }
    exception Unexpected_type of Type.t * Type.t
    module FloatIntMultimap = Tools.Multimap (Tools.ComparableFloat) (Tools.ComparableInt)
    let summarize_distance
        ?(keep_at_most = Some 2) ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m prefix =
      if m.which <> DMatrix then
        Unexpected_type (m.which, DMatrix) |> raise;
      let n_cols = Array.length m.matrix.idx_to_col_names in
      let req_len =
        match keep_at_most with
        | None -> n_cols
        | Some at_most -> at_most
      and f_n_cols = float_of_int n_cols in
      let process_row buf row_idx =
        let name = m.matrix.idx_to_row_names.(row_idx) and row = m.matrix.storage.(row_idx)
        and distr = ref FloatIntMultimap.empty in
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
        Printf.bprintf buf "\"%s\"\t%.15g\t%.15g\t%.15g\t%.15g" name mean stddev median mad;
        FloatIntMultimap.iteri
          (fun i dist col_idx ->
            if i < eff_len then
              Printf.bprintf buf "\t\"%s\"\t%.15g\t%.15g"
                m.matrix.idx_to_col_names.(col_idx) dist ((dist -. mean) /. stddev))
          !distr;
        Printf.bprintf buf "\n" in
      let n_rows = Array.length m.matrix.idx_to_row_names and processed_rows = ref 0 and buf = Buffer.create 1048576
      and fname = make_filename_summary prefix in
      let output = open_out fname in
      (* Parallel section *)
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
            process_row buf i
          done;
          hi_row - lo_row + 1, Buffer.contents buf)
        (fun (n_processed, block) ->
          Printf.fprintf output "%s" block;
          let old_processed_rows = !processed_rows in
          processed_rows := !processed_rows + n_processed;
          if verbose && !processed_rows / 10000 > old_processed_rows / 10000 then
            Printf.eprintf "\r(%s): Writing distance digest to file '%s': done %d/%d lines%!"
              __FUNCTION__ fname !processed_rows n_rows)
        threads;
      if verbose then
        Printf.eprintf "\r(%s): Writing distance digest to file '%s': done %d/%d lines.\n%!"
          __FUNCTION__ fname n_rows n_rows;
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
    val multiply_matrix_vector:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Type.t -> t -> t -> t
    (* Compute distances between the rows of a matrix *)
    val get_distance_matrix:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Space.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Space.Distance.t -> Float.Array.t -> t -> t -> t
    exception Unexpected_type of Type.t * Type.t
    val summarize_distance:
      ?keep_at_most:int option -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    (* Binary marshalling of the matrix *)
    val to_channel: out_channel -> t -> unit
    exception Incompatible_archive_version of string * string
    val of_channel: in_channel -> t
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> Type.t -> string -> t
  end
)

