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
    (* All file name arguments are in fact _prefixes_ *)
    val of_file: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> Type.t -> string -> t
    (* This one discards type information - use at your own risk *)
    val to_file: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val transpose_single_threaded: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    val multiply_matrix_vector:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Type.t -> t -> t -> t
    val get_distance_matrix:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Matrix.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Matrix.Distance.t -> Float.Array.t -> t -> t -> t
    (* Binary marshalling of the matrix *)
    val to_channel: out_channel -> t -> unit
    exception Incompatible_archive_version of string * string
    val of_channel: in_channel -> t
    val to_binary: ?verbose:bool -> t -> string -> unit
    exception Unexpected_type of Type.t * Type.t
    val of_binary: ?verbose:bool -> Type.t -> string -> t
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
    (* The two following functions implement automatic file naming *)
    let make_filename_table which = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ "." ^ Type.to_string which ^ ".txt"
    let make_filename_binary which name =
      match which, name with
      | _, _ when String.length name >= 5 && String.sub name 0 5 = "/dev/" -> name
      | Type.Twisted, prefix | DMatrix, prefix -> prefix ^ "." ^ Type.to_string which
      | Twister, _ | Inertia, _ -> assert false (* Should always be done through KPopTwister *)
    (* We redefine the implementation for Matrix in order to set the correct KPop types *)
    let of_file ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) which prefix =
      { which; matrix = make_filename_table which prefix |> Matrix.of_file ~threads ~bytes_per_step ~verbose }
    let to_file ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) m prefix =
      make_filename_table m.which prefix |> Matrix.to_file ~precision ~threads ~elements_per_step ~verbose m.matrix
    let transpose_single_threaded ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) m =
      { m with matrix = Matrix.transpose_single_threaded ~threads ~elements_per_step ~verbose m.matrix }
    let transpose ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m =
      { m with matrix = Matrix.transpose ~threads ~elements_per_step ~verbose m.matrix }
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) m v =
      Matrix.multiply_matrix_vector ~threads ~elements_per_step ~verbose m.matrix v
    let multiply_matrix_matrix ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) which m1 m2 =
      { which; matrix = Matrix.multiply_matrix_matrix ~threads ~elements_per_step ~verbose m1.matrix m2.matrix }
    let get_distance_matrix ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) distance metric m =
      { which = DMatrix;
        matrix = Matrix.get_distance_matrix ~threads ~elements_per_step ~verbose distance metric m.matrix }
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    let get_distance_rowwise ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) distance metric m1 m2 =
      { which = DMatrix;
        matrix = Matrix.get_distance_rowwise ~threads ~elements_per_step ~verbose distance metric m1.matrix m2.matrix }
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
      { which = Type.of_string which; matrix = (input_value input: Matrix.t) }
    exception Unexpected_type of Type.t * Type.t
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
      Printf.eprintf "\r                                                             \r%!";
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
    val to_files: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    exception Mismatched_twister_files of string array * string array * string array
    val of_files: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    (* *)
    exception Incompatible_twister_and_twisted
    exception Wrong_number_of_columns of int * int * int
    exception Header_expected of string
    exception Wrong_format of int * string
    exception Duplicate_label of string
    val add_twisted_from_files:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> KPopMatrix.t -> string list -> KPopMatrix.t
    (* *)
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t
  end
= struct
    type t = {
      twister: KPopMatrix.t;
      inertia: KPopMatrix.t
    }
    let empty = { twister = KPopMatrix.empty Twister; inertia = KPopMatrix.empty Inertia }
    (* *)
    let to_files ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) tr prefix =
      KPopMatrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.twister prefix;
      KPopMatrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.inertia prefix
    exception Mismatched_twister_files of string array * string array * string array
    let of_files ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) prefix =
      let twister = KPopMatrix.of_file ~threads ~bytes_per_step ~verbose Twister prefix
      and inertia = KPopMatrix.of_file ~threads ~bytes_per_step ~verbose Inertia prefix in
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
        Mismatched_twister_files
            (twister.matrix.idx_to_row_names, inertia.matrix.idx_to_col_names, inertia.matrix.idx_to_row_names)
          |> raise
      end;
      { twister; inertia }
    (* Strictly speaking, we return the _transposed_ of the matrix product here *)
    exception Incompatible_twister_and_twisted
    exception Wrong_number_of_columns of int * int * int
    exception Header_expected of string
    exception Wrong_format of int * string
    exception Duplicate_label of string
    let [@warning "-27"] add_twisted_from_files
        ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) twister twisted fnames =
      if twisted.KPopMatrix.which <> Twisted then
        KPopMatrix.Unexpected_type (twisted.KPopMatrix.which, Twisted) |> raise;
      let twisted_idx_to_col_names =
        if twisted.matrix = Matrix.empty then
          twister.twister.matrix.idx_to_row_names
        else twisted.matrix.idx_to_col_names in
      if twister.twister.matrix.idx_to_row_names <> twisted_idx_to_col_names then
        raise Incompatible_twister_and_twisted;
      (* We invert the table *)
      let num_twister_cols = Array.length twister.twister.matrix.idx_to_col_names in
      let twister_col_names_to_idx = Hashtbl.create num_twister_cols in
      Array.iteri
        (fun i name ->
          Hashtbl.add twister_col_names_to_idx name i)
        twister.twister.matrix.idx_to_col_names;
      (* We decompose the existing twisted matrix *)
      let res = ref Tools.StringMap.empty in
      Array.iteri
        (fun i name ->
          res := Tools.StringMap.add name twisted.matrix.storage.(i) !res)
        twisted.matrix.idx_to_row_names;
      (* First we read spectra from the files.
         We have to conform the k-mers to the ones in the twister.
         As a bonus, we already know the size of the resulting vector *)
      let fnames = Array.of_list fnames in
      let n = Array.length fnames and file_idx = ref 0 and file = open_in fnames.(0) |> ref and line_num = ref 0
      and labels = ref ("", "") and num_spectra = ref 0 in
      (* Parallel section *)
      Tools.Parallel.process_stream_chunkwise
        (fun () ->
          if !file_idx = n then
            raise End_of_file
          else begin
            let buf = ref [] in
            (* Each file can contain one or more spectra - we read the next one *)
            begin try
              while true do
                let line_s = input_line !file in
                incr line_num;
                let line = Tools.Split.on_char_as_array '\t' line_s in
                let l = Array.length line in
                if l <> 2 then
                  Wrong_number_of_columns (!line_num, l, 2) |> raise;
                (* Each file must begin with a header *)
                if !line_num = 1 && line.(0) <> "" then
                  Header_expected line_s |> raise;
                if line.(0) = "" then begin
                  (* New header *)
                  labels := snd !labels, line.(1);
                  if !line_num > 1 then begin
                    incr num_spectra;
                    raise Exit
                  end
                end else
                  (* A regular line. The first element is the hash, the second one the count *)
                  Tools.Misc.accum buf (line.(0), line.(1))
              done
            with
            | End_of_file ->
              close_in !file;
              labels := snd !labels, "";
              incr num_spectra;
              Printf.eprintf "\r[%d/%d] File '%s': Read %d %s on %d %s.\n%!"
                (!file_idx + 1) n fnames.(!file_idx)
                !num_spectra (Tools.Misc.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                !line_num (Tools.Misc.pluralize_int "line" !line_num);
              incr file_idx;
              if !file_idx < n then begin
                file := open_in fnames.(!file_idx);
                line_num := 0
              end
            | Exit ->
              ()
            end;
            if verbose && !file_idx < n then
              Printf.eprintf "\r[%d/%d] File '%s': Read %d %s on %d %s%!"
                (!file_idx + 1) n fnames.(!file_idx)
                !num_spectra (Tools.Misc.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                !line_num (Tools.Misc.pluralize_int "line" !line_num);
            (* The lines are passed in reverse order, but that does not really matter much
                as the final order is determined by the twister *)
            fst !labels, !buf
          end)
        (fun (label, rev_lines) ->
          let n = List.length rev_lines
          (* We allocate the column vector to be transformed *)
          and col = Float.Array.make num_twister_cols 0. in
          List.iteri
            (fun i (name, v) ->
              let v =
                try
                  float_of_string v
                with _ ->
                  Wrong_format (n - i, v) |> raise in
              match Hashtbl.find_opt twister_col_names_to_idx name with
              | Some idx ->
                (* If there are repeated k-mers, we just accumulate them *)
                Float.Array.get col idx +. v |> Float.Array.set col idx
              | None ->
                (* In this case, we just discard the k-mer *)
                ())
            rev_lines;
          (* We first normalise and then transform the spectrum *)
          let acc = ref 0. in
          Float.Array.iter
            (fun el ->
              acc := !acc +. el)
            col;
          let acc = !acc in
          if acc <> 0. then
            Float.Array.iteri
              (fun i el ->
                el /. acc |> Float.Array.set col i)
              col;
          label, KPopMatrix.multiply_matrix_vector ~threads:1 ~elements_per_step:n ~verbose:false twister.twister col)
        (fun (label, row) ->
          (* The transformed column vector becomes a row *)
          match Tools.StringMap.find_opt label !res with
          | None ->
            res := Tools.StringMap.add label row !res
          | Some _ ->
            Duplicate_label label |> raise)
        threads;
      let n = Tools.StringMap.cardinal !res in
      let idx_to_row_names = Array.make n ""
      and storage = Array.make n (Float.Array.create 0) in
      let module StringMap = Tools.Misc.Map(Tools.ComparableString) in
      StringMap.iteri
        (fun i label row ->
          idx_to_row_names.(i) <- label;
          storage.(i) <- row)
        !res;
      { KPopMatrix.which = Twisted; 
        matrix = { Matrix.idx_to_col_names = twisted_idx_to_col_names; idx_to_row_names; storage } }
    (* *)
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopTwister"
    let to_binary ?(verbose = false) t prefix =
      let fname = make_filename_binary prefix in
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      KPopMatrix.to_channel output t.twister;
      KPopMatrix.to_channel output t.inertia;
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) prefix =
      let fname = make_filename_binary prefix in
      let input = open_in fname in
      if verbose then
        Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let twister = KPopMatrix.of_channel input in
      let inertia = KPopMatrix.of_channel input in
      close_in input;
      if KPopMatrix.Type.Twister <> twister.which then
        KPopMatrix.Unexpected_type (Twister, twister.which) |> raise;
      if KPopMatrix.Type.Inertia <> inertia.which then
        KPopMatrix.Unexpected_type (Inertia, inertia.which) |> raise;
      if verbose then
        Printf.eprintf " done.\n%!";
      { twister; inertia }

  end

module RegisterType =
  struct
    type t =
      | Twister
      | Twisted
      | Distances
    exception Invalid_register_type
    let of_string = function
      | "T" -> Twister
      | "t" -> Twisted
      | "d" -> Distances
      | _ ->
        Tools.Argv.usage ();
        raise Invalid_register_type
  end

(* There are three registers in this program, one per DB type *)
type to_do_t =
  | Empty of RegisterType.t
  | Binary_to_register of RegisterType.t * string
  | Tables_to_register of RegisterType.t * string
  | Add_files_to_twisted of string list
  | Register_to_binary of RegisterType.t * string
  | Set_precision of int
  | Register_to_tables of RegisterType.t * string
  | Set_distance of Matrix.Distance.parameters_t
  | Set_metric of Matrix.Distance.Metric.t
  | Distances_from_twisted_binary of string
  (*| Summary*)

module Defaults =
  struct
    let distance = { Matrix.Distance.which = "euclidean"; power = 2. }
    let metric = Matrix.Distance.Metric.Flat
    let precision = 15
    let threads = Tools.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    (* This is just to correctly propagate values in the program *)
    let distance = ref Defaults.distance
    let program = ref []
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let version = "0.9"

let header =
  Printf.sprintf begin
    "This is the KPopTwistDB program (version %s)\n%!" ^^
    " (c) 2022 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!"
  end version

let _ =
  let module TA = Tools.Argv in
  let module TS = Tools.Split in
  let module TM = Tools.Misc in
  TA.set_header header;
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    [], None, [ "Actions (executed delayed and in order of specification):" ], TA.Optional, (fun _ -> ());
    [ "-e"; "--empty" ],
      Some "T|t|d",
      [ "load an empty twisted database into the specified register";
        " (T=twister; t=twisted; d=distance)" ],
      TA.Optional,
      (fun _ -> Empty (TA.get_parameter () |> RegisterType.of_string) |> TM.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "T|t|d <binary_file_prefix>",
      [ "load the specified binary database into the specified register";
        " (T=twister; t=twisted; d=distance).";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwister; .KPopTwisted; or .KPopDMatrix, respectively)" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Binary_to_register (register_type, TA.get_parameter ()) |> TM.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "T|t|d <table_file_prefix>",
      [ "load the specified tabular database(s) into the specified register";
        " (T=twister; t=twisted; d=distance).";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwister.txt and .KPopInertia.txt; .KPopTwisted.txt;";
        "  or .KPopDMatrix, respectively)" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Tables_to_register (register_type, TA.get_parameter ()) |> TM.accum Parameters.program);
    [ "-f"; "-F"; "-k"; "-K"; "-a"; "-A";
      "--add"; "--files"; "--kmers"; "--add-files"; "--add-kmers" ],
      Some "<k-mer_table_file_name>[','...','<k-mer_table_file_name>]",
      [ "twist k-mers from the specified files through the transformation";
        "present in the twister register, and add the results";
        "to the database present in the twisted register" ],
      TA.Optional,
      (fun _ -> Add_files_to_twisted (TA.get_parameter () |> TS.on_char_as_list ',') |> TM.accum Parameters.program);
    [ "--distance"; "--distance-function"; "--set-distance"; "--set-distance-function" ],
      Some "'euclidean'|'minkowski'",
      [ "set the function to be used when computing distances" ],
      TA.Default (fun () -> Defaults.distance.which),
      (fun _ -> Set_distance { !Parameters.distance with which = TA.get_parameter () } |> TM.accum Parameters.program);
    [ "--distance-power"; "--set-distance-power" ],
      Some "<non_negative_float>",
      [ "set the power to be used when computing (Minkowski) distances" ],
      TA.Default (fun () -> string_of_float Defaults.distance.power),
      (fun _ ->
        Set_distance { !Parameters.distance with power = TA.get_parameter_float_non_neg () }
          |> TM.accum Parameters.program);
    [ "--metric"; "--metric-function"; "--set-metric"; "--set-metric-function" ],
      Some "'flat'",
      [ "set the metric function to be used when computing distances" ],
      TA.Default (fun () -> Matrix.Distance.Metric.to_string Defaults.metric),
      (fun _ -> Set_metric (TA.get_parameter () |> Matrix.Distance.Metric.of_string) |> TM.accum Parameters.program);
    [ "-d"; "--distances"; "--compute-distances"; "--compute-distances-twisted" ],
      Some "<twisted_binary_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension .KPopTwisted)" ],
      TA.Optional,
      (fun _ -> Distances_from_twisted_binary (TA.get_parameter ()) |> TM.accum Parameters.program);

(* SOMETHING TO SUMMARISE DISTANCES? *)

    [ "-o"; "--output" ],
      Some "T|t|d <binary_file_prefix>",
      [ "dump the database present in the specified register";
        " (T=twister; t=twisted; d=distance) to the specified binary file.";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwister; .KPopTwisted; or .KPopDMatrix, respectively)" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_binary (register_type, TA.get_parameter ()) |> TM.accum Parameters.program);
    [ "--precision"; "--set-precision"; "--set-table-precision" ],
      Some "<positive_integer>",
      [ "set the number of precision digits to be used when outputting numbers" ],
      TA.Default (fun () -> string_of_int Defaults.precision),
      (fun _ -> Set_precision (TA.get_parameter_int_pos ()) |> TM.accum Parameters.program);
    [ "-O"; "--Output" ],
      Some "T|t|d <table_file_prefix>",
      [ "dump the database present in the specified register";
        " (T=twister; t=twisted; d=distance) to the specified tabular file(s).";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwister.txt and .KPopInertia.txt; .KPopTwisted.txt;";
        "  or .KPopDMatrix, respectively)" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_tables (register_type, TA.get_parameter ()) |> TM.accum Parameters.program);
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
  (* These are the two registers available to the program *)
  let twister = ref KPopTwister.empty and twisted = KPopMatrix.empty Twisted |> ref
  and distance = Matrix.Distance.of_parameters Defaults.distance |> ref
  and metric = Matrix.Distance.Metric.of_string "flat" |> Matrix.Distance.Metric.compute |> ref
  and distances = KPopMatrix.empty DMatrix |> ref and precision = ref Defaults.precision in

  (* The addition of an exception handler would be nice *)

  List.iter
    (function
      | Empty RegisterType.Twister ->
        twister := KPopTwister.empty
      | Empty RegisterType.Twisted ->
        twisted := KPopMatrix.empty Twisted
      | Empty RegisterType.Distances ->
        distances := KPopMatrix.empty DMatrix
      | Binary_to_register (RegisterType.Twister, prefix) ->
        twister := KPopTwister.of_binary ~verbose:!Parameters.verbose prefix
      | Binary_to_register (RegisterType.Twisted, prefix) ->
        twisted := KPopMatrix.of_binary ~verbose:!Parameters.verbose Twisted prefix;
        if !twisted.which <> Twisted then
          KPopMatrix.Unexpected_type (!twisted.which, Twisted) |> raise
      | Binary_to_register (RegisterType.Distances, prefix) ->
        distances := KPopMatrix.of_binary ~verbose:!Parameters.verbose DMatrix prefix;
        if !distances.which <> DMatrix then
          KPopMatrix.Unexpected_type (!distances.which, DMatrix) |> raise
      | Tables_to_register (RegisterType.Twister, prefix) ->
        twister := KPopTwister.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose prefix
      | Tables_to_register (RegisterType.Twisted, prefix) ->
        twisted := KPopMatrix.of_file ~threads:!Parameters.threads ~verbose:!Parameters.verbose Twisted prefix;
        if !twisted.which <> Twisted then
          KPopMatrix.Unexpected_type (!twisted.which, Twisted) |> raise
      | Tables_to_register (RegisterType.Distances, prefix) ->
        distances := KPopMatrix.of_file ~threads:!Parameters.threads ~verbose:!Parameters.verbose DMatrix prefix;
        if !distances.which <> DMatrix then
          KPopMatrix.Unexpected_type (!distances.which, DMatrix) |> raise
      | Add_files_to_twisted fnames ->
        twisted :=
          KPopTwister.add_twisted_from_files
            ~threads:!Parameters.threads ~verbose:!Parameters.verbose !twister !twisted fnames
      | Register_to_binary (RegisterType.Twister, prefix) ->
        KPopTwister.to_binary ~verbose:!Parameters.verbose !twister prefix
      | Register_to_binary (RegisterType.Twisted, prefix) ->
        KPopMatrix.to_binary ~verbose:!Parameters.verbose !twisted prefix
      | Register_to_binary (RegisterType.Distances, prefix) ->
        KPopMatrix.to_binary ~verbose:!Parameters.verbose !distances prefix
      | Set_precision prec ->
        precision := prec
      | Register_to_tables (RegisterType.Twister, prefix) ->
        KPopTwister.to_files ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose
          !twister prefix
      | Register_to_tables (RegisterType.Twisted, prefix) ->
        KPopMatrix.to_file
          ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose !twisted prefix
      | Register_to_tables (RegisterType.Distances, prefix) ->
        KPopMatrix.to_file
          ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose !distances prefix
      | Set_distance dist ->
        distance := Matrix.Distance.of_parameters dist
      | Set_metric metr ->
        metric := Matrix.Distance.Metric.compute metr
      | Distances_from_twisted_binary prefix ->
        let twisted_db = KPopMatrix.of_binary ~verbose:!Parameters.verbose Twisted prefix in
        if !twisted.which <> Twisted then
          KPopMatrix.Unexpected_type (!twisted.which, Twisted) |> raise;
        distances :=
          KPopMatrix.get_distance_rowwise ~verbose:!Parameters.verbose ~threads:!Parameters.threads
            !distance begin
              (* We compute the metric vector *)
              if Array.length !twister.inertia.matrix.idx_to_row_names > 0 then
                (* We use the metric induced by inertia *)
                !metric !twister.inertia.matrix.storage.(0)
              else begin
                (* We assume a flat inertia *)
                let len = Array.length !twisted.matrix.idx_to_col_names in
                if len = 0 then
                  Float.Array.create 0
                else
                  1. /. float_of_int len |> Float.Array.make len |> !metric
              end 
            end !twisted twisted_db)
    program

