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
        type t (* We hide the type to implement constraints *)
        exception Negative_metric_element of float
        exception Unsorted_metric_vector of Float.Array.t
        val compute: t -> Float.Array.t -> Float.Array.t
        exception Unknown_metric of string
        exception Negative_power of float
        exception Negative_tightness of float
        val of_string: string -> t
        val to_string: t -> string
      end
    type t (* We hide the type to implement constraints *)
    exception Incompatible_lengths of int * int * int
    (* What happens when the vectors have incompatible lengths *)
    type mode_t =
      | Fail
      | Infinity
    val set_mode: mode_t -> unit
    val compute: t -> Float.Array.t -> Float.Array.t -> Float.Array.t -> float
    exception Unknown_distance of string
    exception Negative_power of float
    val of_string: string -> t
    val to_string: t -> string
  end
= struct
    module Metric =
      struct
        type t =
          | Flat
          | Power of float
          (* Parameters are: threshold (1/n/t), tightness left, tightness right *)
          | Sigmoid of float * float * float * float
        exception Negative_metric_element of float
        exception Unsorted_metric_vector of Float.Array.t
        let compute = function
          | Flat ->
            Float.Array.map
             (fun el ->
               if el < 0. then
                 Negative_metric_element el |> raise;
               1.)
          | Power power ->
            (fun v ->
              (* We normalise with respect to the transformed max *)
              let m = ref 0. in
              Float.Array.iter
                (fun el ->
                  if el < 0. then
                   Negative_metric_element el |> raise;
                  m := el ** power |> max !m)
                v;
              let m = !m in
              if m > 0. then
                Float.Array.map
                  (fun el ->
                    (el ** power) /. m)
                  v
              else
                (* All zeros *)
                v)
          | Sigmoid (power, threshold, l_tightness, r_tightness) ->
            (fun v ->
              (* We normalise to one *)
              let acc = ref 0. in
              Float.Array.iteri
                (fun i el ->
                  if el < 0. then
                    Negative_metric_element el |> raise;
                  if i > 0 && el > Float.Array.get v (i - 1) then
                    Unsorted_metric_vector v |> raise;
                  acc := !acc +. el ** power)
                v;
              let acc = !acc in
              if acc > 0. then begin
                (* Determine inflection point *)
                let f_n = float_of_int (Float.Array.length v) in
                let threshold = 1. /. threshold /. f_n and t = ref (-1) in
                Float.Array.iteri
                  (fun i el ->
                    let el = (el ** power) /. acc in
                    if el < threshold && !t = (-1) then
                      t := i)
                  v;
                if !t <> (-1) then begin
                  let t = (float_of_int !t +. 0.5) /. f_n in
                  Float.Array.mapi
(*
f<-function(x,t=0.5,kl=10,kr=100){a<-ifelse(x<t,x/t,(x-t)/(1-t)); y<-ifelse(x<t,(-(2*a-3)*a*a-1)+1/2/(1-t)*((a-1)*a*a),1/t/2*((a-2)*a*a+a)-(2*a-3)*a*a); (ifelse(y>0,-y*((1+kr)/kr)/((1+y*kr)/kr),-y*((1+kl)/kl)/((1-y*kl)/kl))+1)/2}
*)
                    (fun i _ ->
                      let x = (float_of_int i +. 0.5) /. f_n in
                      let y =
                        if x < t then begin
                          let a = x /. t in
                          let a_a = a *. a in
                          ((-2. *. a +. 3.) *. a_a -. 1.) +. 0.5 /. (1. -. t) *. ((a -. 1.) *. a_a)
                        end else begin
                          let a = (x -. t) /. (1. -. t) in
                          let a_a = a *. a in
                          0.5 /. t *. ((a -. 2.) *. a_a +. a) -. (2. *. a -. 3.) *. a_a
                        end in
                      let res =
                        if y > 0. then
                          -. y *. ((1. +. r_tightness) /. r_tightness) /. ((1. +. y *. r_tightness) /. r_tightness)
                        else
                          -. y *. ((1. +. l_tightness) /. l_tightness) /. ((1. -. y *. l_tightness) /. l_tightness) in
                      (1. +. res) /. 2. |> max 0. |> min 1.)
                    v
                end else
                  v
              end else
                (* All zeros *)
                v)
        exception Unknown_metric of string
        exception Negative_power of float
        exception Negative_tightness of float
        let of_string_re = Str.regexp "[(,)]"
        let of_string name =
          match name with
          | "flat" ->
            Flat
          | s ->
            match Str.full_split of_string_re s with
            | [ Text "power"; Delim "("; Text power; Delim ")" ] ->
              let power =
                try
                  float_of_string power
                with _ ->
                  Unknown_metric s |> raise in
              if power < 0. then
                Negative_power power |> raise;
              Power power
            | [ Text "sigmoid"; Delim "(";
                Text power; Delim ","; Text thresh; Delim ","; Text l_tight; Delim ","; Text r_tight;
                Delim ")" ] ->
              let power, thresh, l_tight, r_tight =
                try
                  float_of_string power, float_of_string thresh, float_of_string l_tight, float_of_string r_tight
                with _ ->
                  Unknown_metric s |> raise in
              if power < 0. then
                Negative_power power |> raise;
              Sigmoid (power, thresh, l_tight, r_tight)
            | _ ->
              Unknown_metric s |> raise
        let to_string = function
          | Flat ->
            "flat"
          | Power power ->
            Printf.sprintf "power(%.15g)" power
          | Sigmoid (power, thresh, l_tight, r_tight) ->
            Printf.sprintf "sigmoid(%.15g,%.15g,%.15g,%.15g)" power thresh l_tight r_tight
      end
    type t =
      | Euclidean
      | Minkowski of float (* Theoretically speaking, the parameter should be an integer *)
    exception Incompatible_lengths of int * int * int
    type mode_t =
      | Fail
      | Infinity
    let mode = ref Fail
    let set_mode new_mode = mode := new_mode
    let compute f m a b =
      let length_a = Float.Array.length a and length_m = Float.Array.length m and length_b = Float.Array.length b in
      if length_a <> length_m || length_m <> length_b then begin
        match !mode with
        | Fail -> Incompatible_lengths (length_a, length_m, length_b) |> raise
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
    exception Unknown_distance of string
    exception Negative_power of float
    let of_string_re = Str.regexp "[()]"
    let of_string = function
      | "euclidean" ->
        Euclidean
      | s ->
        match Str.full_split of_string_re s with
        | [ Text "minkowski"; Delim "("; Text power; Delim ")" ] ->
          let power =
            try
              float_of_string power
            with _ ->
              Unknown_distance s |> raise in
          if power < 0. then
            Negative_power power |> raise;
          Minkowski power
        | _ ->
          Unknown_distance s |> raise
    let to_string = function
      | Euclidean -> "euclidean"
      | Minkowski power -> Printf.sprintf "minkowski(%.15g)" power
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
      let red_d = d - 1 and res = Float.Array.create d and elts_done = ref 0 in
      (* We decorate each vector element coordinate with the respective value *)
      for i = 0 to red_d do
        let acc = ref 0. in
        Float.Array.iter2
          (fun el_1 el_2 ->
            acc := !acc +. (el_1 *. el_2))
          m.storage.(i) v;
        Float.Array.set res i !acc;
        incr elts_done;
        if verbose && !elts_done mod 100 = 0 then
          Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done d
      done;
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
      let red_d = d - 1 and res = Float.Array.make d 0. and elts_done = ref 0 in
      (* We decorate each vector element coordinate with the respective value *)
      for i = 0 to red_d do
        let m_v = m.storage.(i) and acc = ref 0. in
        Tools.IntMap.iter
          (fun j el ->
            acc := !acc +. (Float.Array.get m_v j *. el))
          s_v.elements;
        Float.Array.set res i !acc;
        incr elts_done;
        if verbose && !elts_done mod 100 = 0 then
          Printf.eprintf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done d
      done;
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
      let distance = Distance.compute distance metric in
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
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Distance.t -> Float.Array.t -> t -> t -> t
  end
)

