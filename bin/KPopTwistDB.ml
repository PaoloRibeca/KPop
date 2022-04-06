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
    module Type:
      sig
        type t =
          | Twister
          | Inertia
          | Twisted
          | DMatrix
        val to_string: t -> string
        exception Unknown_type of string
        val of_string: string -> t
      end
    type t = {
      which: Type.t;
      matrix: Matrix.t
    }
    val empty: Type.t -> t
    val of_file: ?verbose:bool -> Type.t -> string -> t
    (* This one discards type information - use at your own risk *)
    val to_file: t -> string -> unit
    val transpose_single_threaded: t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> t -> t
    val multiply_matrix_vector: ?threads:int -> ?elements_per_step:int -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> Type.t -> t -> t -> t
    val get_distance_matrix: ?threads:int -> ?elements_per_step:int -> Matrix.Distance.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise: ?threads:int -> ?elements_per_step:int -> Matrix.Distance.t -> t -> t -> t
    (* Binary marshalling of the matrix *)
    val to_channel: out_channel -> t -> unit
    val of_channel: in_channel -> t
    val to_binary: t -> string -> unit
    val of_binary: string -> t
    (* To be used when needed *)
    exception TwisterExpected
    exception InertiaExpected
    exception TwistedExpected
    exception DMatrixExpected
  end
= struct
    module Type =
      struct
        type t =
          | Twister
          | Inertia
          | Twisted
          | DMatrix
        let to_string = function
          | Twister -> "KPopTwister"
          | Inertia -> "KPopInertia"
          | Twisted -> "KPopTwisted"
          | DMatrix -> "KPopDMatrix"
        exception Unknown_type of string
        let of_string = function
          | "KPopTwister" -> Twister
          | "KPopInertia" -> Inertia
          | "KPopTwisted" -> Twisted
          | "KPopDMatrix" -> DMatrix
          | w ->
            Unknown_type w |> raise
      end
    type t = {
      which: Type.t;
      matrix: Matrix.t
    }
    let empty which =
      { which; matrix = Matrix.empty }
    (* We redefine the implementation for Matrix in order to set the correct KPop types *)
    let of_file ?(verbose = false) which fname =
      { which; matrix = Matrix.of_file ~verbose fname }
    let to_file m =
      Matrix.to_file m.matrix
    let transpose_single_threaded m =
      { m with matrix = Matrix.transpose_single_threaded m.matrix }
    let transpose ?(threads = 64) ?(elements_per_step = 100) m =
      { m with matrix = Matrix.transpose ~threads ~elements_per_step m.matrix }
    let multiply_matrix_vector ?(threads = 64) ?(elements_per_step = 100) m v =
      Matrix.multiply_matrix_vector ~threads ~elements_per_step m.matrix v
    let multiply_matrix_matrix ?(threads = 64) ?(elements_per_step = 100) which m1 m2 =
      { which; matrix = Matrix.multiply_matrix_matrix ~threads ~elements_per_step m1.matrix m2.matrix }
    let get_distance_matrix ?(threads = 64) ?(elements_per_step = 100) distance m =
      { which = DMatrix; matrix = Matrix.get_distance_matrix ~threads ~elements_per_step distance m.matrix }
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    let get_distance_rowwise ?(threads = 64) ?(elements_per_step = 100) distance m1 m2 =
      { which = DMatrix; matrix = Matrix.get_distance_rowwise ~threads ~elements_per_step distance m1.matrix m2.matrix }
    (* *)
    let archive_version = "2022-04-03"
    (* *)
    exception Incompatible_archive_version of string * string
    let to_channel output m =
      Type.to_string m.which |> output_value output;
      archive_version |> output_value output;
      output_value output m.matrix
    let to_binary m fname =
      let output = open_out fname in
      Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      to_channel output m;
      close_out output;
      Printf.eprintf " done.\n%!"
    let of_channel input =
      let which = (input_value input: string) in
      let version = (input_value input: string) in
      if version <> archive_version then
        Incompatible_archive_version (which, version) |> raise;
      { which = Type.of_string which; matrix = (input_value input: Matrix.t) }
    let of_binary fname =
      let input = open_in fname in
      Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let res = of_channel input in
      close_in input;
      Printf.eprintf " done.\n%!";
      res
    exception TwisterExpected
    exception InertiaExpected
    exception TwistedExpected
    exception DMatrixExpected

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


module [@warning "-32"] KPopTwister:
  sig
    type t = {
      twister: KPopMatrix.t; (* Coordinate transformation *)
      inertia: KPopMatrix.t  (* Variance per coordinate *)
    }
    val empty: t
    val to_files: t -> string -> unit
    exception MismatchedTwisterFiles of string array * string array * string array
    val of_files: ?verbose:bool -> string -> t
    (* *)
    exception IncompatibleTwisterAndTwisted
    exception WrongNumberOfColumns of int * int * int
    exception HeaderExpected of string
    exception WrongFormat of int * string
    val add_twisted_from_files:
      ?threads:int -> ?elements_per_step:int -> t -> KPopMatrix.t -> string list -> KPopMatrix.t
    (* *)
    val to_binary: t -> string -> unit
    val of_binary: string -> t

  end
= struct
    type t = {
      twister: KPopMatrix.t;
      inertia: KPopMatrix.t
    }
    let empty = { twister = KPopMatrix.empty Twister; inertia = KPopMatrix.empty Inertia }
    (* *)
    let to_files tr prefix =
      prefix ^ ".KPopTwister.txt" |> KPopMatrix.to_file tr.twister;
      prefix ^ ".KPopInertia.txt" |> KPopMatrix.to_file tr.inertia
    exception MismatchedTwisterFiles of string array * string array * string array
    let of_files ?(verbose = false) prefix =
      let twister = prefix ^ ".KPopTwister.txt" |> KPopMatrix.of_file ~verbose Twister
      and inertia = prefix ^ ".KPopInertia.txt" |> KPopMatrix.of_file ~verbose Inertia in
      (* Let's run at least some checks *)
      if begin
        inertia.matrix.idx_to_row_names <> [| "inertia" |] ||
        twister.matrix.idx_to_row_names <> inertia.matrix.idx_to_col_names
      end then begin
        Printf.eprintf "ERROR: twister.idx_to_row_names:";
        Array.iter (fun el -> Printf.eprintf "\t\"%s\"" el) twister.matrix.idx_to_row_names;
        Printf.eprintf "\nERROR: inertia.idx_to_col_names:";
        Array.iter (fun el -> Printf.eprintf "\t\"%s\"" el) inertia.matrix.idx_to_col_names;
        Printf.eprintf "\nERROR: inertia.idx_to_row_names:";
        Array.iter (fun el -> Printf.eprintf "\t\"%s\"" el) inertia.matrix.idx_to_row_names;
        Printf.eprintf "\n%!";
        MismatchedTwisterFiles
            (twister.matrix.idx_to_row_names, inertia.matrix.idx_to_col_names, inertia.matrix.idx_to_row_names)
          |> raise
      end;
      { twister; inertia }

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
    exception IncompatibleTwisterAndTwisted
    let add_twisted_from_files ?(threads = 64) ?(elements_per_step = 100) twister twisted fnames =
      if twisted.KPopMatrix.which <> Twisted then
        raise KPopMatrix.TwistedExpected;
      if twister.twister.matrix.idx_to_row_names <> twisted.matrix.idx_to_col_names then
        raise IncompatibleTwisterAndTwisted;
      (* We invert the table *)
      let col_names_to_idx = Hashtbl.create 1024 in
      Array.iteri
        (fun i name ->
          Hashtbl.add col_names_to_idx name i)
        twister.twister.matrix.idx_to_col_names;
      (* First we read spectra from the files.
         We have to conform the k-mers to the ones in the twister.
         As a bonus, we already know the size of the resulting vector *)
      let curr_label = ref "" and curr_alloc = Float.Array.create 0 |> ref
      and res_labels = ref [] and res_allocs = ref [] in
      add_files
        (fun _ label ->
          curr_label := label;
          curr_alloc := Float.Array.make (Array.length twister.twister.matrix.idx_to_col_names) 0.)
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
          (* We first normalise and then transform the spectrum *)
          let acc = ref 0. in
          Float.Array.iter
            (fun el ->
              acc := !acc +. el)
            !curr_alloc;
          let acc = !acc in
          if acc <> 0. then
            Float.Array.iteri
              (fun i el ->
                el /. acc |> Float.Array.set !curr_alloc i)
              !curr_alloc;
          curr_alloc :=
            KPopMatrix.multiply_matrix_vector ~threads ~elements_per_step twister.twister !curr_alloc;
          Tools.Misc.accum res_labels !curr_label;
          Tools.Misc.accum res_allocs !curr_alloc)
        fnames;
      { KPopMatrix.which = Twisted;
        matrix = {
          Matrix.idx_to_col_names = twisted.matrix.idx_to_col_names;
          idx_to_row_names = Tools.Misc.array_of_rlist !res_labels |> Array.append twisted.matrix.idx_to_row_names;
          storage = Tools.Misc.array_of_rlist !res_allocs |> Array.append twisted.matrix.storage } }
    (* *)
    let to_binary t fname =
      let output = open_out fname in
      Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      KPopMatrix.to_channel output t.twister;
      KPopMatrix.to_channel output t.inertia;
      close_out output;
      Printf.eprintf " done.\n%!"
    let of_binary fname =
      let input = open_in fname in
      Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let twister = KPopMatrix.of_channel input in
      let inertia = KPopMatrix.of_channel input in
      close_in input;
      Printf.eprintf " done.\n%!";
      { twister; inertia }

  end

(* There are two registers here *)
type to_do_t =
  | Empty_twister
  | Empty_twisted
  | Binary_to_twister of string
  | Tables_to_twister of string (* Prefix *)
  | Binary_to_twisted of string
  | Table_to_twisted of string
  | Twister_to_binary of string
  | Twister_to_tables of string (* Prefix *)
  | Twisted_to_binary of string
  | Twisted_to_table of string
  | Add_files_to_twisted of string list
  | Twisted_distances_binary of string * string
  | Twisted_distances_table of string * string
  | Distances_to_table of string * string
  | Table_to_distances of string * string
  (*| Summary*)

module Defaults =
  struct
    let distance = Matrix.Distance.euclidean
    let threads = 1
    let verbose = true (*false*)
  end

module Parameters =
  struct
    let program = ref []
    let distance = ref Defaults.distance
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let version = "0.3"

let _ =
  Printf.eprintf "This is the KPopTwistDB program (version %s)\n%!" version;
  Printf.eprintf " (c) 2022 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!";
  let module TA = Tools.Argv in
  let module TS = Tools.Split in
  let module TM = Tools.Misc in
  TA.parse [
    [], None, [ "Actions (executed delayed and in order of specification):" ], TA.Optional, (fun _ -> ());
    [ "-E"; "--Empty" ],
      None,
      [ "put an empty twister database into the twister register" ],
      TA.Optional,
      (fun _ -> Empty_twister |> TM.accum Parameters.program);
    [ "-e"; "--empty" ],
      None,
      [ "put an empty twisted database into the twisted register" ],
      TA.Optional,
      (fun _ -> Empty_twisted |> TM.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "<binary_file_prefix>",
      [ "load to the twister register the database present in the specified file";
        " (which must have extension .KPopTwister)" ],
      TA.Optional,
      (fun _ -> Binary_to_twister (TA.get_parameter () ^ ".KPopTwister") |> TM.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load to the twisted register the database present in the specified file";
        " (which must have extension .KPopTwisted)" ],
      TA.Optional,
      (fun _ -> Binary_to_twisted (TA.get_parameter () ^ ".KPopTwisted") |> TM.accum Parameters.program);
    [ "--twister-from-tables" ],
      Some "<table_file_prefix>",
      [ "load to the twister register the database present in the specified files";
        " (which must have extensions .KPopTwister.txt and .KPopInertia.txt)" ],
      TA.Optional,
      (fun _ -> Tables_to_twister (TA.get_parameter ()) |> TM.accum Parameters.program);
    [ "--twisted-from-table" ],
      Some "<table_file_prefix>",
      [ "load to the twisted register the database present in the specified file";
        " (which must have extension .KPopTwisted.txt)" ],
      TA.Optional,
      (fun _ -> Table_to_twisted (TA.get_parameter () ^ ".KPopTwisted.txt") |> TM.accum Parameters.program);
    [ "-f"; "-F"; "-k"; "-K"; "-a"; "-A";
      "--add"; "--files"; "--kmers";
      "--add-files"; "--add-kmers"; "--twisted-add-files"; "--twisted-add-kmers" ],
      Some "<k-mer_table_file_name>[','...','<k-mer_table_file_name>]",
      [ "twist k-mers from the specified files through the transformation";
        "present in the twister register, and add the result";
        "to the database present in the twisted register" ],
      TA.Optional,
      (fun _ -> Add_files_to_twisted (TA.get_parameter () |> TS.on_char_as_list ',') |> TM.accum Parameters.program);


(* NEED TO ADD DISTANCE SETTING *)


    [ "-d"; "--distances-binary"; "--compute-distances-binary"; "--twisted-compute-distances-binary" ],
      Some "<twisted_binary_file_prefix> <distance_binary_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension .KPopTwisted),";
        "and output them to the specified distance binary file";
        " (which will be given extension .KPopDMatrix)" ],
      TA.Optional,
      (fun _ ->
        let twisted = TA.get_parameter () ^ ".KPopTwisted" in
        let distances = TA.get_parameter () ^ ".KPopDMatrix" in
        Twisted_distances_binary (twisted, distances) |> TM.accum Parameters.program);
    [ "-D"; "--distances-table"; "--compute-distances-table"; "--twisted-compute-distances-table" ],
      Some "<twisted_binary_file_prefix> <distance_table_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension .KPopTwisted),";
        "and output them to the specified distance text file";
        " (which will be given extension .KPopDMatrix.txt)" ],
      TA.Optional,
      (fun _ ->
        let twisted = TA.get_parameter () ^ ".KPopTwisted" in
        let distances = TA.get_parameter () ^ ".KPopDMatrix.txt" in
        Twisted_distances_table (twisted, distances) |> TM.accum Parameters.program);
    [ "--distances-binary-to-table" ],
      Some "<distance_binary_file_prefix> <distance_table_file_prefix>",
      [ "output distances contained in the specified binary file";
        " (which must have extension .KPopDMatrix)";
        "to the specified text file";
        " (which will be given extension .KPopDMatrix.txt)" ],
      TA.Optional,
      (fun _ ->
        let binary = TA.get_parameter () ^ ".KPopDMatrix" in
        let table = TA.get_parameter () ^ ".KPopDMatrix.txt" in
        Distances_to_table (binary, table) |> TM.accum Parameters.program);
    [ "--distances-table-to-binary" ],
      Some "<distance_table_file_prefix> <distance_binary_file_prefix>",
      [ "output distances contained in the specified text file";
        " (which must have extension .KPopDMatrix.txt)";
        "to the specified binary file";
        " (which will be given extension .KPopDMatrix)" ],
      TA.Optional,
      (fun _ ->
        let table = TA.get_parameter () ^ ".KPopDMatrix.txt" in
        let binary = TA.get_parameter () ^ ".KPopDMatrix" in
        Table_to_distances (table, binary) |> TM.accum Parameters.program);

(* SOMETHING TO SUMMARISE DISTANCES? *)

    [ "-O"; "--Output" ],
      Some "<binary_file_prefix>",
      [ "dump the database present in the twister register to the specified file";
        " (which will be given extension .KPopTwister)" ],
      TA.Optional,
      (fun _ -> Twister_to_binary (TA.get_parameter () ^ ".KPopTwister") |> TM.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "dump the database present in the twisted register to the specified file";
        " (which will be given extension .KPopTwisted)" ],
      TA.Optional,
      (fun _ -> Twisted_to_binary (TA.get_parameter () ^ ".KPopTwisted") |> TM.accum Parameters.program);
    [ "--twister-to-tables" ],
      Some "<text_file_prefix>",
      [ "dump the database present in the twister register to the specified files";
        " (which will be given extension .KPopTwister.txt and .KPopInertia.txt)" ],
      TA.Optional,
      (fun _ -> Twister_to_tables (TA.get_parameter ()) |> TM.accum Parameters.program);
    [ "--twisted-to-table" ],
      Some "<text_file_prefix>",
      [ "dump the database present in the twisted register to the specified file";
        " (which will be given extension .KPopTwisted.txt)" ],
      TA.Optional,
      (fun _ -> Twisted_to_table (TA.get_parameter () ^ ".KPopTwisted.txt") |> TM.accum Parameters.program);
    [], None, [ "Miscellaneous (executed immediately):" ], TA.Optional, (fun _ -> ());
    [ "-t"; "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (fun () -> string_of_bool !Parameters.verbose),
      (fun _ -> Parameters.verbose := true);
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
  let current_twister = ref KPopTwister.empty and current_twisted = KPopMatrix.empty Twisted |> ref in

  (* The addition of an exception handler would be nice *)

  List.iter
    (function
      | Empty_twister ->
        current_twister := KPopTwister.empty
      | Empty_twisted ->
        current_twisted := KPopMatrix.empty Twisted
      | Binary_to_twister fname ->
        current_twister := KPopTwister.of_binary fname
      | Tables_to_twister prefix ->
        current_twister := KPopTwister.of_files ~verbose:!Parameters.verbose prefix
      | Binary_to_twisted fname ->
        current_twisted := KPopMatrix.of_binary fname
      | Table_to_twisted fname ->
        current_twisted := KPopMatrix.of_file ~verbose:!Parameters.verbose Twisted fname
      | Twister_to_binary fname ->
        KPopTwister.to_binary !current_twister fname
      | Twister_to_tables prefix ->
        KPopTwister.to_files !current_twister prefix
      | Twisted_to_binary fname ->
        KPopMatrix.to_binary !current_twisted fname
      | Twisted_to_table fname ->
        KPopMatrix.to_file !current_twisted fname
      | Add_files_to_twisted fnames ->
        current_twisted :=
          KPopTwister.add_twisted_from_files ~threads:!Parameters.threads !current_twister !current_twisted fnames
      | Twisted_distances_binary (twisted_db_fname, fname) | Twisted_distances_table (twisted_db_fname, fname) as w ->
        let twisted_db = KPopMatrix.of_binary twisted_db_fname in
        if twisted_db.which <> Twisted then
          raise KPopMatrix.TwistedExpected;
        let distances =
          KPopMatrix.get_distance_rowwise ~threads:!Parameters.threads !Parameters.distance !current_twisted twisted_db in
        begin match w with
        | Twisted_distances_binary _ ->
          KPopMatrix.to_binary distances fname
        | Twisted_distances_table _ ->
          KPopMatrix.to_file distances fname
        | _ -> assert false
        end
      | Distances_to_table (binary, fname) ->
        let distances = KPopMatrix.of_binary binary in
        if distances.which <> DMatrix then
          raise KPopMatrix.DMatrixExpected;
        KPopMatrix.to_file distances fname
      | Table_to_distances (table, fname) ->
        KPopMatrix.to_binary (KPopMatrix.of_file ~verbose:!Parameters.verbose DMatrix table) fname)
    program

