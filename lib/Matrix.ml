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
   They all have signature: float array -> float array -> float *)
module Distance:
  sig
    type t = Float.Array.t -> Float.Array.t -> float
    exception IncompatibleLengths of int * int
    (* What happens when the vectors have incompatible lengths *)
    type mode_t =
      | Fail
      | Infinity
    val set_mode: mode_t -> unit
    val euclidean: t
  
  end
= struct
    type t = Float.Array.t -> Float.Array.t -> float
    exception IncompatibleLengths of int * int
    type mode_t =
      | Fail
      | Infinity
    let mode = ref Fail
    let set_mode new_mode = mode := new_mode
    let _proto_ f a b =
      let length_a = Float.Array.length a and length_b = Float.Array.length b in
      if length_a <> length_b then begin
        match !mode with
        | Fail -> IncompatibleLengths (length_a, length_b) |> raise
        | Infinity -> infinity
      end else
        f a b
    let euclidean =
      _proto_
        (fun a b ->
          let acc = ref 0. in
          Float.Array.iter2
            (fun el_a el_b ->
              let diff = el_a -. el_b in
              acc := !acc +. (diff *. diff))
            a b;
          sqrt !acc)
  
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
    let to_file m filename =
      let output = open_out filename in
      if Array.length m.storage > 0 then begin
        (* We output the column names *)
        Array.iteri
          (fun i name ->
            if i > 0 then
              Printf.fprintf output "\t";
            Printf.fprintf output "\"%s\"" name)
          m.idx_to_col_names;
        Printf.fprintf output "\n";
        (* We output the rows *)
        Array.iteri
          (fun i row ->
            (* We output the row name *)
            m.idx_to_row_names.(i) |> Printf.fprintf output "\"%s\"";
            Float.Array.iter (Printf.fprintf output "\t%g") row;
            Printf.fprintf output "\n")
          m.storage
      end;
      close_out output
    let add_row_name names row_idx name =
      names := Misc._resize_array_ Array.length (fun n -> Array.make n "") Array.blit !names row_idx;
      !names.(row_idx) <- name
    let add_row storage row_idx dim =
      storage :=
        Misc._resize_array_ Array.length (fun n -> Array.make n (Float.Array.create 0)) Array.blit !storage row_idx;
      !storage.(row_idx) <- Float.Array.create dim;
      !storage.(row_idx)
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
    exception WrongNumberOfColumns of int * int
    let of_file ?(verbose = false) filename =
      let input = open_in filename and line_num = ref 0
      and num_cols = ref 0 and idx_to_col_names = ref [||] and idx_to_row_names = ref [||]
      and storage = ref (Float.Array.create 0 |> Array.make 16) and elts_read = ref 0 in
      begin try
        while true do
          let line = input_line input |> Tools.Split.on_char_as_array '\t' in
          incr line_num;
          if !line_num = 1 then begin
            (* We process the header *)
            let l = Array.length line in
            if l > 0 then begin
              (* We assume the matrix always to have row names, and ignore the first name in the header if present *)
              num_cols := l - 1;
              idx_to_col_names := Array.make !num_cols "";
              Array.iteri
                (fun i name ->
                  if i > 0 then
                    !idx_to_col_names.(i - 1) <- strip_quotes name)
                line
            end
          end else begin
            (* A regular line.
               The first element is the name *)
            let l = Array.length line in
            if l <> !num_cols + 1 then
              WrongNumberOfColumns (l, !num_cols + 1) |> raise;
            let line_idx = !line_num - 1 in
            let array = add_row storage line_idx !num_cols in
            Array.iteri
              (fun i el ->
                if i = 0 then
                  (* The first element is the name *)
                  add_row_name idx_to_row_names line_idx (strip_quotes el)
                else begin
                  float_of_string el |> Float.Array.set array (i - 1);
                  incr elts_read;
                  if !elts_read mod 100000 = 0 && verbose then
                    Printf.printf "\r(%s): At row %d of file '%s': Read %d elements%!            \r"
                      __FUNCTION__ !line_num filename !elts_read
                end)
              line
          end
        done
      with End_of_file ->
        close_in input;
        if verbose then
          Printf.printf "\r(%s): At row %d of file '%s': Read %d elements.            \n%!"
            __FUNCTION__ !line_num filename !elts_read
      end;
      { (* At this point idx_to_row_names and storage might be longer than needed - we resize them *)
        idx_to_col_names = !idx_to_col_names;
        idx_to_row_names = Array.sub !idx_to_row_names 0 !line_num;
        storage = Array.sub !storage 0 !line_num }
    let [@warning "-32"] transpose_single_threaded m =
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_col_names;
        storage =
          Array.init (Array.length m.idx_to_col_names)
            (fun old_col ->
              Float.Array.init (Array.length m.idx_to_row_names)
                (fun old_row -> Float.Array.get m.storage.(old_row) old_col)) }
    let [@warning "-32"] transpose ?(threads = 64) ?(elements_per_step = 1) m =
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
            if !elts_done mod elements_per_step = 0 then
              Printf.printf "\r(%s): Done %d/%d rows%!            \r" __FUNCTION__ !elts_done row_num))
        threads;
      Printf.printf "\r(%s): Done %d/%d rows.            \n%!" __FUNCTION__ !elts_done row_num;
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_col_names;
        storage = storage }
    exception IncompatibleGeometries of string array * string array
    let multiply_matrix_vector ?(threads = 64) ?(elements_per_step = 100) m v =
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
            if !elts_done mod elements_per_step = 0 then
              Printf.printf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done d))
        threads;
      Printf.printf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done d;
      res
    let multiply_matrix_matrix ?(threads = 64) ?(elements_per_step = 100) m1 m2 =
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
            if !elts_done mod elements_per_step = 0 then
              Printf.printf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done prod))
        threads;
      Printf.printf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done prod;
      { idx_to_col_names = m2.idx_to_col_names;
        idx_to_row_names = m1.idx_to_row_names;
        storage = storage }
    let get_distance_matrix ?(threads = 64) ?(elements_per_step = 100) distance m =
      let d = Array.length m.idx_to_row_names in
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
            if !elts_done mod elements_per_step = 0 then
              Printf.printf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ !elts_done total))
        threads;
      Printf.printf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ !elts_done total;
      Printf.printf "(%s): Symmetrizing...%!" __FUNCTION__;
      (* We symmetrise the matrix *)
      let red_d = d - 1 in
      for i = 0 to red_d do
        for j = i + 1 to red_d do
          Float.Array.get storage.(j) i |> Float.Array.set storage.(i) j
        done
      done;
      Printf.printf " done.\n%!";
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_row_names;
        storage = storage }
    let get_distance_rowwise ?(threads = 64) ?(elements_per_step = 100) distance m1 m2 =
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
            if !elts_done mod elements_per_step = 0 then
              Printf.printf "\r(%s): Done %d/%d elements%!            \r" __FUNCTION__ prod !elts_done))
        threads;
      Printf.printf "\r(%s): Done %d/%d elements.            \n%!" __FUNCTION__ prod !elts_done;
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
    exception WrongNumberOfColumns of int * int
    val of_file: ?verbose:bool -> string -> t
    val to_file: t -> string -> unit
    val transpose_single_threaded: t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> t -> t
    exception IncompatibleGeometries of string array * string array
    val multiply_matrix_vector: ?threads:int -> ?elements_per_step:int -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> t -> t -> t
    val get_distance_matrix: ?threads:int -> ?elements_per_step:int -> Distance.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise: ?threads:int -> ?elements_per_step:int -> Distance.t -> t -> t -> t

  end
)

