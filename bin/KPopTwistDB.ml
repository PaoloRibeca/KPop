(* We read in a matrix which has conditions as row names
    and a (large) number of tags (genes, k-mers, etc.) as column names.
   Keeping with the convention accepted by R, the first row would be a header,
    and the first column the row names.
   Names might be quoted *)

(* A number of distance functions.
   They all have signature: float array -> float array -> float *)
module DistanceFunction:
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

module IntMap = Tools.IntMap

module Matrix:
  sig
    type t = {
      (* We number rows and columns starting from 0 *)
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      (* *)
      storage: Float.Array.t array
    }
    exception WrongNumberOfColumns of int * int
    val of_file: string -> t
    (*val transpose: t -> t*)
    val get_distance: ?threads:int -> ?elements_per_step:int -> DistanceFunction.t -> t -> t
    val to_file: t -> string -> unit

  end
= struct
    type t = {
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      storage: Float.Array.t array
    }
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
    let _resize_array_ length make blit a idx =
      let l = length a in
      if l < idx + 1 then begin
        let res = make (max (idx + 1) (l * 14 / 10)) in
        blit a 0 res 0 l;
        res
      end else
        a
    let add_row_name names row_idx name =
      names := _resize_array_ Array.length (fun n -> Array.make n "") Array.blit !names row_idx;
      !names.(row_idx) <- name
    let add_row storage row_idx dim =
      storage :=
        _resize_array_ Array.length (fun n -> Array.make n (Float.Array.create 0)) Array.blit !storage row_idx;
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
    let of_file filename =
      let input = open_in filename and line_num = ref 0
      and num_cols = ref 0 and idx_to_col_names = ref [||] and idx_to_row_names = ref [||]
      and storage = ref (Float.Array.create 0 |> Array.make 16) and elts_read = ref 0 in
      begin try
        while true do
          let line = input_line input |> Tools.Split.on_char_as_array '\t' in
          if !line_num = 0 then begin
            (* We process the header *)
            let l = Array.length line in
            if l > 0 then
              if line.(0) = "" then begin
                (* The first name is empty: as many columns as l - 1 *)
                num_cols := l - 1;
                idx_to_col_names := Array.make !num_cols "";
                Array.iteri
                  (fun i name ->
                    if i > 0 then
                      !idx_to_col_names.(i - 1) <- strip_quotes name)
                  line
              end else begin
                (* The first name is not empty: l columns *)
                num_cols := l;
                idx_to_col_names := Array.map strip_quotes line
              end
          end else begin
            (* A regular line.
               The first element is the name *)
            let l = Array.length line in
            if l <> !num_cols + 1 then
              WrongNumberOfColumns (l, !num_cols + 1) |> raise;
            let red_line_num = !line_num - 1 in
            let array = add_row storage red_line_num !num_cols in
            Array.iteri
              (fun i el ->
                if i = 0 then
                  (* The first element is the name *)
                  add_row_name idx_to_row_names red_line_num (strip_quotes el)
                else begin
                  float_of_string el |> Float.Array.set array (i - 1);
                  incr elts_read;
                  if !elts_read mod 100000 = 0 then
                    Printf.printf "\r                                        \rAt row %d: Read %d elements%!"
                      (!line_num + 1) !elts_read
                end)
              line
          end;
          incr line_num
        done
      with End_of_file ->
        close_in input;
        Printf.printf "\r%!"
      end;
      let red_line_num = !line_num - 1 in
      { (* At this point idx_to_row_names and storage might be longer than needed - we resize them *)
        idx_to_col_names = !idx_to_col_names;
        idx_to_row_names = Array.sub !idx_to_row_names 0 red_line_num;
        storage = Array.sub !storage 0 red_line_num }
    let get_distance ?(threads = 64) ?(elements_per_step = 100) distance m =
      let d = Array.length m.storage in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let storage = Array.init d (fun _ -> Float.Array.create d) in
      (* Generate points to be computed by the parallel processs *)
      let i = ref 0 and j = ref 0 and end_reached = ref false and elts_done = ref 0 in
      Tools.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file
          else begin
            (* We only compute 1/2 of the matrix, and symmetrise at the end of the computation *)
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
            (* The actual memory allocation for the result is only done here *)
            Float.Array.set storage.(i) j dist;
            incr elts_done;
            if !elts_done mod 100 = 0 then
              Printf.printf "\r                                        \rAt (%d,%d): Read %d elements%!"
                i j !elts_done))
        threads;
      Printf.printf "\r                                        \rSymmetrizing...%!";
      (* We symmetrise the matrix *)
      let red_d = d - 1 in
      for i = 0 to red_d do
        for j = i + 1 to red_d do
          Float.Array.get storage.(j) i |> Float.Array.set storage.(i) j
        done
      done;
      Printf.printf "\r                                        \r%!";
      { idx_to_col_names = m.idx_to_row_names;
        idx_to_row_names = m.idx_to_row_names;
        storage = storage }
  
  end

let () =
  let matrix = Matrix.of_file "/dev/stdin" in
  let distance = Matrix.get_distance DistanceFunction.euclidean matrix in
  Matrix.to_file distance Sys.argv.(1)
