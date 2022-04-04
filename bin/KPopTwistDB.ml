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

module [@warning "-32"] KPopMatrix:
  sig
    type matrix_t =
      | Twister
      | Inertia
      | Twisted
      | DMatrix
    type t = {
      which: matrix_t;
      matrix: Matrix.t
    }
    val of_file: matrix_t -> string -> t
    (* This one discards type information - use at your own risk *)
    val to_file: t -> string -> unit
    (* TRANSPOSITION IS CURRENTLY UNTESTED AND HENCE DISABLED *)
    (* val transpose_single_threaded: t -> t
       val transpose: ?threads:int -> ?elements_per_step:int -> t -> t *)
    val multiply_matrix_vector: ?threads:int -> ?elements_per_step:int -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> matrix_t -> t -> t -> t
    val get_distance_matrix: ?threads:int -> ?elements_per_step:int -> Matrix.DistanceFunction.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise: ?threads:int -> ?elements_per_step:int -> Matrix.DistanceFunction.t -> t -> t -> t
    (* Binary marshalling of the matrix *)
    val to_binary: t -> string -> unit
    val of_binary: string -> t
  end
= struct
    type matrix_t =
      | Twister
      | Inertia
      | Twisted
      | DMatrix
    type t = {
      which: matrix_t;
      matrix: Matrix.t
    }
    (* We redefine the implementation for Matrix in order to set the correct KPop types *)
    let of_file which fname =
      { which; matrix = Matrix.of_file fname }
    let to_file m =
      Matrix.to_file m.matrix
    (* TRANSPOSITION IS CURRENTLY UNTESTED AND HENCE DISABLED *)
    (*
       let transpose_single_threaded m =
         { m with matrix = Matrix.transpose_single_threaded m.matrix }
       let transpose ?(threads = 64) ?(elements_per_step = 100) m =
         { m with matrix = Matrix.transpose ~threads ~elements_per_step m.matrix }
    *)
    let multiply_matrix_vector ?(threads = 64) ?(elements_per_step = 100) m v =
      Matrix.multiply_matrix_vector ~threads ~elements_per_step m.matrix v
    let multiply_matrix_matrix ?(threads = 64) ?(elements_per_step = 100) which m1 m2 =
      { which; matrix = Matrix.multiply_matrix_matrix ~threads ~elements_per_step m1.matrix m2.matrix }
    let get_distance_matrix ?(threads = 64) ?(elements_per_step = 100) distance m =
      { which = DMatrix; matrix = Matrix.get_distance_matrix ~threads ~elements_per_step distance m.matrix }
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    let get_distance_rowwise ?(threads = 64) ?(elements_per_step = 100) distance m1 m2 =
      { which = DMatrix; matrix = Matrix.get_distance_rowwise ~threads ~elements_per_step distance m1.matrix m2.matrix }


    let archive_version = "2022-04-03"

    exception Incompatible_archive_version of string
    let to_binary m fname =
      let output = open_out fname in
      Printf.eprintf "%s: Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      "KPopTwistDB_" ^ archive_version |> output_value output;
      output_value output m;
      close_out output;
      Printf.eprintf " done.\n%!"
    let of_binary fname =
      let input = open_in fname in
      Printf.eprintf "%s: Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let version = (input_value input: string) in
      if version <> "KPopTwistDB_" ^ archive_version then
        Incompatible_archive_version version |> raise;
      let res = (input_value input: t) in
      close_in input;
      Printf.eprintf " done.\n%!";
      res

  end


(*
    module FloatStringMultimap = Tools.OrderedMultimap (Tools.ComparableFloat) (Tools.ComparableString)
      let res = ref FloatStringMultimap.empty in
      (* We find the median and filter the result *)


            res := FloatStringMultimap.add dist v.idx_to_col_names.(i) !res; 


      let res = !res and req_len =
        match keep_at_most with
        | None -> d
        | Some at_most -> at_most
      and eff_len = ref 0 and median_pos = d / 2 |> ref and median = ref 0. in
      FloatStringMultimap.iter_set
        (fun dist set ->
          let set_len = FloatStringMultimap.ValSet.cardinal set in
          if !median_pos >= 0 && !median_pos - set_len < 0 then
            median := dist;
          median_pos := !median_pos - set_len;
          if !eff_len < req_len then
            eff_len := !eff_len + set_len)
        res;
      let eff_len = !eff_len and median = !median in
      (* At this point, we can allocate storage... *)
      let storage = Array.init eff_len (fun _ -> Float.Array.create 2) and row_names = Array.make eff_len "" in
      (* ...and fill it *)
      FloatStringMultimap.iteri
        (fun i dist key ->
          row_names.(i) <- key;
          let arr = storage.(i) in
          Float.Array.set arr 0 dist;
          (dist -. median) /. median |> Float.Array.set arr 1)
        res;
      Printf.printf "\r                                                             \r%!";
      { idx_to_col_names = [| "distance"; "z_score" |];
        idx_to_row_names = row_names;
        storage = storage }
*)


module [@warning "-32"] Twister:
  sig
    type t = {
      twister: Matrix.t; (* Coordinate transformation *)
      inertia: Matrix.t  (* Variance per coordinate *)
    }
    val to_files: t -> string -> unit
    val of_files: string -> t

    exception WrongNumberOfColumns of int * int * int
    exception HeaderExpected of string
    exception WrongFormat of int * string
    val twisted_from_files: ?threads:int -> ?elements_per_step:int -> t -> string list -> Matrix.t

    exception Incompatible_archive_version of string
    val to_binary: t -> string -> unit
    val of_binary: string -> t

  end
= struct
    type t = {
      twister: Matrix.t;
      inertia: Matrix.t
    }

    let to_files tr prefix =
      prefix ^ ".KPopTwister.twister.txt" |> Matrix.to_file tr.twister;
      prefix ^ ".KPopTwister.inertia.txt" |> Matrix.to_file tr.inertia
    let of_files prefix =
      { twister = prefix ^ ".KPopTwister.twister.txt" |> Matrix.of_file;
        inertia = prefix ^ ".KPopTwister.inertia.txt" |> Matrix.of_file }

    (* When everything is properly separated, this should go into a common library file *)

    exception WrongNumberOfColumns of int * int * int
    exception HeaderExpected of string
    exception WrongFormat of int * string
    let add_files f_header f_line f_end fnames =
      let n = List.length fnames in
      List.iteri
        (fun i fname ->
          let input = open_in fname and line_num = ref 1 in
          (* Each file can contain one or more spectra *)
          begin try
            while true do
              let line_s = input_line input in
              let line = Tools.Split.on_char_as_array '\t' line_s in
              let l = Array.length line in
              if l <> 2 then
                WrongNumberOfColumns (!line_num, l, 2) |> raise;
              (* Each file must begin with a header *)
              if !line_num = 1 && line.(0) <> "" then
                HeaderExpected line_s |> raise;
              if line.(0) = "" then begin
                if !line_num > 1 then
                  f_end !line_num;
                (* Header *)
                f_header !line_num line.(1)
              end else
                (* A regular line. The first element is the hash, the second one the count *)
                f_line !line_num line.(0) line.(1);
              incr line_num;
              if !line_num mod 10000 = 0 then
                Printf.eprintf "\r[%d/%d] File '%s': Read %d lines%!" (i + 1) n fname !line_num
            done
          with End_of_file ->
            close_in input;
            f_end !line_num;
            Printf.eprintf "\r[%d/%d] File '%s': Read %d lines\n%!" (i + 1) n fname !line_num;
          end)
        fnames

    (* Strictly speaking, we return the _transposed_ of the matrix product here *)
    let twisted_from_files ?(threads = 64) ?(elements_per_step = 100) tw fnames =
      (* We invert the table *)
      let col_names_to_idx = Hashtbl.create 1024 in
      Array.iteri
        (fun i name ->
          Hashtbl.add col_names_to_idx name i)
        tw.twister.idx_to_col_names;
      (* First we read spectra from the files.
         We have to conform the k-mers to the ones in the twister.
         As a bonus, we already know the size of the resulting vector *)
      let curr_label = ref "" and curr_alloc = Float.Array.create 0 |> ref
      and res_labels = ref [] and res_allocs = ref [] in
      add_files
        (fun _ label ->
          curr_label := label;
          curr_alloc := Float.Array.make (Array.length tw.twister.idx_to_col_names) 0.)
        (fun line name v ->
          let v =
            try
              float_of_string v
            with _ ->
              WrongFormat (line, v) |> raise in
          match Hashtbl.find_opt col_names_to_idx name with
          | Some idx ->
            (* If there are repeated k-mers, we just accumulate them *)
            Float.Array.get !curr_alloc idx +. v |> Float.Array.set !curr_alloc idx
          | None ->
            (* In this case, we just discard the k-mer *)
            ())
        (fun _ ->
          (* Then, we transform the spectra *)
          curr_alloc :=
            (*KPop*)Matrix.multiply_matrix_vector ~threads ~elements_per_step tw.twister !curr_alloc;
          Tools.Misc.accum res_labels !curr_label;
          Tools.Misc.accum res_allocs !curr_alloc)
        fnames;
      { Matrix.idx_to_col_names = tw.twister.idx_to_row_names;
        idx_to_row_names = Tools.Misc.array_of_rlist !res_labels;
        storage = Tools.Misc.array_of_rlist !res_allocs }

    let archive_version = "2022-04-03"

    exception Incompatible_archive_version of string
    let to_binary t fname =
      let output = open_out fname in
      Printf.eprintf "%s: Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      "KPopTwister_" ^ archive_version |> output_value output;
      output_value output t;
      close_out output;
      Printf.eprintf " done.\n%!"
    let of_binary fname =
      let input = open_in fname in
      Printf.eprintf "%s: Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let version = (input_value input: string) in
      if version <> "KPopTwister_" ^ archive_version then
        Incompatible_archive_version version |> raise;
      let res = (input_value input: t) in
      close_in input;
      Printf.eprintf " done.\n%!";
      res

  end


let () =
  let matrix = Matrix.of_file "/dev/stdin" in
  let distance = Matrix.get_distance_matrix Matrix.DistanceFunction.euclidean matrix in
  Matrix.to_file distance Sys.argv.(1)
