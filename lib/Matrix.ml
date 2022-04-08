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

(* A number of distance functions.
   They all have signature: float array -> float array -> float array -> float *)
module Distance:
  sig
    module Metric:
      (* Functions to deduce a metric from a vector *)
      sig
        type t =
          | Flat
        val compute: t -> Float.Array.t -> Float.Array.t
        exception Unknown_metric of string
        val of_string: string -> t
        val to_string: t -> string
      end
    type t =
      | Euclidean
      | Minkowski of float (* Theoretically speaking, the parameter should be an integer *)
    exception IncompatibleLengths of int * int * int
    (* What happens when the vectors have incompatible lengths *)
    type mode_t =
      | Fail
      | Infinity
    val set_mode: mode_t -> unit
    val compute: t -> Float.Array.t -> Float.Array.t -> Float.Array.t -> float
    type parameters_t = {
      which: string;
      power: float
    }
    exception Unknown_distance of string
    val of_parameters: parameters_t -> t
    val to_parameters: t -> parameters_t
  end
= struct
    module Metric =
      struct
        type t =
          | Flat
        let compute = function
          | Flat ->
            Float.Array.map
             (fun _ ->
               1.)
        exception Unknown_metric of string
        let of_string name =
          match name with
          | "flat" ->
            Flat
          | w ->
            Unknown_metric w |> raise
        let to_string = function
          | Flat -> "flat"
      end
    type t =
      | Euclidean
      | Minkowski of float
    exception IncompatibleLengths of int * int * int
    type mode_t =
      | Fail
      | Infinity
    let mode = ref Fail
    let set_mode new_mode = mode := new_mode
    let compute f m a b =
      let length_a = Float.Array.length a and length_m = Float.Array.length m and length_b = Float.Array.length b in
      if length_a <> length_m || length_m <> length_b then begin
        match !mode with
        | Fail -> IncompatibleLengths (length_a, length_m, length_b) |> raise
        | Infinity -> infinity
      end else
        match f with
        | Euclidean ->
          let acc = ref 0. in
          Float.Array.iteri
            (fun i el_a ->
              let diff = el_a -. Float.Array.get b i in
              acc := !acc +. (diff *. diff *. Float.Array.get m i))
            a;
          sqrt !acc
        | Minkowski power ->
          let acc = ref 0. in
          Float.Array.iteri
            (fun i el_a ->
              let diff = el_a -. Float.Array.get b i |> abs_float in
              acc := !acc +. ((diff ** power) *. Float.Array.get m i))
            a;
          !acc ** (1. /. power)
    type parameters_t = {
      which: string;
      power: float
    }
    exception Unknown_distance of string
    let of_parameters parameters =
      match parameters.which with
      | "euclidean" ->
        Euclidean
      | "minkowski" ->
        Minkowski parameters.power
      | w ->
        Unknown_distance w |> raise
    let to_parameters = function
      | Euclidean -> { which = "euclidean"; power = 2. }
      | Minkowski power -> { which = "minkowski"; power }
  end

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
include (
  struct
    type t = {
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      storage: Float.Array.t array
    }
    let empty =
      { idx_to_col_names = [||]; idx_to_row_names = [||]; storage = [||] }
    let [@warning "-27"] to_file
        ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) m filename =
      let output = open_out filename in
      if Array.length m.storage > 0 then begin
        (* We output the column names *)
        Printf.fprintf output "\"\"";
        Array.iter
          (fun name ->
            Printf.fprintf output "\t\"%s\"" name)
          m.idx_to_col_names;
        Printf.fprintf output "\n";
        (* We output the rows *)
        Array.iteri
          (fun i row ->
            (* We output the row name *)
            m.idx_to_row_names.(i) |> Printf.fprintf output "\"%s\"";
            Float.Array.iter (Printf.fprintf output "\t%.*g" precision) row;
            Printf.fprintf output "\n")
          m.storage
      end;
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
    exception WrongNumberOfColumns of int * int * int
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
                  Tools.Misc.accum res (!line_num, line);
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
                WrongNumberOfColumns (line_num, l, num_cols + 1) |> raise;
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
              Tools.Misc.accum idx_to_row_names name;
              Tools.Misc.accum storage numbers;
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
        idx_to_row_names = Tools.Misc.array_of_rlist !idx_to_row_names;
        storage = Tools.Misc.array_of_rlist !storage }
    let [@warning "-27"] transpose_single_threaded ?(threads = 1) ?(elements_per_step = 1) ?(verbose = false) m =
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
                Tools.Misc.accum res !i;
                incr i;
                if !i = row_num then begin (* The original columns *)
                  end_reached := true;
                  raise Exit
                end
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
    exception IncompatibleGeometries of string array * string array
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m v =
      if Array.length m.idx_to_col_names <> Float.Array.length v then
        IncompatibleGeometries (m.idx_to_col_names, Array.make (Float.Array.length v) "") |> raise;
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
                Tools.Misc.accum res !i;
                incr i;
                if !i = d then begin
                  end_reached := true;
                  raise Exit
                end;
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
        IncompatibleGeometries (m1.idx_to_col_names, m2.idx_to_row_names) |> raise;
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
                Tools.Misc.accum res (!i, !j);
                incr j;
                if !j = col_num then begin
                  incr i;
                  if !i = row_num then begin
                    end_reached := true;
                    raise Exit
                  end;
                  j := 0
                end
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
    let get_distance_matrix ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) distance metric m =
      let distance = Distance.compute distance metric and d = Array.length m.idx_to_row_names in
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
                Tools.Misc.accum res (!i, !j);
                incr j;
                if !j > !i then begin
                  incr i;
                  if !i = d then begin
                    end_reached := true;
                    raise Exit
                  end;
                  j := 0
                end
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
      let distance = Distance.compute distance metric in
      if m1.idx_to_col_names <> m2.idx_to_col_names then
        IncompatibleGeometries (m1.idx_to_col_names, m2.idx_to_col_names) |> raise;
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
                Tools.Misc.accum res (!i, !j);
                incr j;
                if !j = r2 then begin
                  incr i;
                  if !i = r1 then begin
                    end_reached := true;
                    raise Exit
                  end;
                  j := 0
                end
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
              Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ prod !elts_done))
        threads;
      if verbose then
        Printf.eprintf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ prod !elts_done;
      { idx_to_col_names = m2.idx_to_row_names;
        idx_to_row_names = m1.idx_to_row_names;
        storage = storage }

  end: sig
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
      Keeping with the convention accepted by R, the first row would be a header,
        and the first column the row names.
      Names might be quoted *)
    exception WrongNumberOfColumns of int * int * int
    val of_file: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    val to_file: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val transpose_single_threaded: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    exception IncompatibleGeometries of string array * string array
    val multiply_matrix_vector:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t -> t
    val get_distance_matrix:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Distance.t -> Float.Array.t -> t -> t -> t

  end
)

