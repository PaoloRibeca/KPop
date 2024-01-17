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
open KPop

module [@warning "-32"] Twister:
  sig
    type t = {
      twister: Matrix.t; (* Coordinate transformation *)
      inertia: Matrix.t  (* Variance per coordinate *)
    }
    val empty: t
    val to_files: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    exception Mismatched_twister_files of string array * string array * string array
    val of_files: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    (* *)
    exception Incompatible_twister_and_twisted
    exception Wrong_number_of_columns of int * int * int
    exception Header_expected of string
    exception Float_expected of string
    exception Duplicate_label of string
    val add_twisted_from_files:
      ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> ?debug:bool ->
      t -> Matrix.t -> string list -> Matrix.t
    (* *)
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t
  end
= struct
    type t = {
      twister: Matrix.t;
      inertia: Matrix.t
    }
    let empty = { twister = Matrix.empty Twister; inertia = Matrix.empty Inertia }
    (* *)
    let to_files ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) tr prefix =
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.twister prefix;
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose tr.inertia prefix
    exception Mismatched_twister_files of string array * string array * string array
    let of_files ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) prefix =
      let twister = Matrix.of_file ~threads ~bytes_per_step ~verbose Twister prefix
      and inertia = Matrix.of_file ~threads ~bytes_per_step ~verbose Inertia prefix in
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
    exception Float_expected of string
    exception Duplicate_label of string
    let [@warning "-27"] add_twisted_from_files
        ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false) ?(debug = false) twister twisted fnames =
      if twisted.Matrix.which <> Twisted then
        Matrix.Unexpected_type (twisted.Matrix.which, Twisted) |> raise;
      let twisted_idx_to_col_names =
        if twisted.matrix = Matrix.Base.empty then
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
      Processes.Parallel.process_stream_chunkwise
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
                  Tools.List.accum buf (line.(0), line.(1))
              done
            with
            | End_of_file ->
              close_in !file;
              labels := snd !labels, "";
              incr num_spectra;
              if verbose then
                Printf.eprintf "%s\r(%s): [%d/%d] File '%s': Read %d %s on %d %s%s%!"
                  Tools.String.TermIO.clear __FUNCTION__ (!file_idx + 1) n fnames.(!file_idx)
                  !num_spectra (Tools.String.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                  !line_num (Tools.String.pluralize_int "line" !line_num) (if !file_idx + 1 = n then ".\n" else "");
              incr file_idx;
              if !file_idx < n then begin
                file := open_in fnames.(!file_idx);
                line_num := 0
              end
            | Exit ->
              ()
            end;
            if verbose && !file_idx < n then
              Printf.eprintf "%s\r(%s): [%d/%d] File '%s': Read %d %s on %d %s%!"
                Tools.String.TermIO.clear __FUNCTION__ (!file_idx + 1) n fnames.(!file_idx)
                !num_spectra (Tools.String.pluralize_int ~plural:"spectra" "spectrum" !num_spectra)
                !line_num (Tools.String.pluralize_int "line" !line_num);
            (* The lines are passed in reverse order, but that does not really matter much
                as the final order is determined by the twister *)
            fst !labels, !buf
          end)
        (fun (label, rev_lines) ->
          let t0 = Sys.time () in
          let s_v = ref Tools.IntMap.empty and acc = ref 0. in
          List.iter
            (fun (name, v) ->
              match Hashtbl.find_opt twister_col_names_to_idx name with
              | Some idx ->
                let v =
                  try
                    float_of_string v
                  with _ ->
                    Float_expected v |> raise in
                acc := !acc +. v;
                s_v := begin
                  match Tools.IntMap.find_opt idx !s_v with
                  | Some vv ->
                    (* If there are repeated k-mers, we accumulate them *)
                    Tools.IntMap.add idx (vv +. v) !s_v
                  | None ->
                    Tools.IntMap.add idx v !s_v
                end
              | None ->
                (* In this case, we just discard the k-mer *)
                ())
            rev_lines;
          let t1 = Sys.time () in
          (* We first normalise and then transform the spectrum *)
          let acc = !acc in
          let s_v = {
            Matrix.Base.length = num_twister_cols;
            elements =
              if acc <> 0. then
                Tools.IntMap.map
                  (fun el ->
                    el /. acc)
                  !s_v
              else
                !s_v
          } in
          let t2 = Sys.time () in
          let res = Matrix.multiply_matrix_sparse_vector_single_threaded ~verbose:false twister.twister s_v in
          let t3 = Sys.time () in
          if debug then
            Printf.eprintf "DEBUG=(lines=%d/%d/%d,%.3g,%.3g,%.3g)\n%!"
              (List.length rev_lines) num_twister_cols (Float.Array.length res) (t1 -. t0) (t2 -. t1) (t3 -. t2);
          label, res)
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
      Tools.StringMap.iteri
        (fun i label row ->
          idx_to_row_names.(i) <- label;
          storage.(i) <- row)
        !res;
      { Matrix.which = Twisted;
        matrix = { Matrix.Base.idx_to_col_names = twisted_idx_to_col_names; idx_to_row_names; storage } }
    (* *)
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopTwister"
    let to_binary ?(verbose = false) t prefix =
      let fname = make_filename_binary prefix in
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      Matrix.to_channel output t.twister;
      Matrix.to_channel output t.inertia;
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) prefix =
      let fname = make_filename_binary prefix in
      let input = open_in fname in
      if verbose then
        Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let twister = Matrix.of_channel input in
      let inertia = Matrix.of_channel input in
      close_in input;
      if Matrix.Type.Twister <> twister.which then
        Matrix.Unexpected_type (Twister, twister.which) |> raise;
      if Matrix.Type.Inertia <> inertia.which then
        Matrix.Unexpected_type (Inertia, inertia.which) |> raise;
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
      | Metrics
    exception Invalid_register_type of string
    let of_string = function
      | "T" -> Twister
      | "t" -> Twisted
      | "d" -> Distances
      | "m" -> Metrics
      | w ->
        Tools.Argv.usage ();
        Invalid_register_type w |> raise
  end

module KeepAtMost =
  struct
    type t = int option
    exception Invalid_keep_at_most of string
    let of_string = function
      | "all" -> None
      | w ->
        try
          let res = int_of_string w in
          if res <= 0 then
            raise Exit;
          Some res
        with _ ->
          Tools.Argv.usage ();
          Invalid_keep_at_most w |> raise
    let to_string = function
      | None -> "all"
      | Some n -> string_of_int n
  end

(* There are three registers in this program, one per DB type *)
type to_do_t =
  | Empty of RegisterType.t
  | Binary_to_register of RegisterType.t * string
  | Tables_to_register of RegisterType.t * string
  | Add_binary_to_register of RegisterType.t * string
  | Add_tables_to_register of RegisterType.t * string
  | Add_kmer_files_to_twisted of string list
  | Register_to_binary of RegisterType.t * string
  | Set_precision of int
  | Register_to_tables of RegisterType.t * string
  | Set_distance of Space.Distance.t
  | Set_distance_normalize of bool
  | Set_metric of Space.Distance.Metric.t
  | Distances_from_twisted_binary of string
  | Set_keep_at_most of KeepAtMost.t
  | Summary_from_twisted_binary of string * string
  | Summary_from_distances of string

module Defaults =
  struct
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalize = true
    let metric = Space.Distance.Metric.of_string "powers(1,1.,2)"
    let precision = 15
    let keep_at_most = Some 2
  end

module Parameters =
  struct
    let program = ref []
    let threads = Processes.Parallel.get_nproc () |> ref
    let verbose = ref false
  end

let info = {
  Tools.Info.name = "KPopTwistDB";
  version = "27";
  date = "17-Jan-2024"
} and authors = [
  "2022-2024", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.make_header info authors [ BiOCamLib.Info.info; KPop.Info.info ] |> TA.set_header;
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    [ "-e"; "--empty" ],
      Some "T|t|d",
      [ "load an empty database into the specified register";
        " (T=twister; t=twisted; d=distance)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics ->
          TA.parse_error "You cannot load content into the metrics register"
        | Twister | Twisted | Distances as register_type ->
          Empty register_type |> Tools.List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "T|t|d <binary_file_prefix>",
      [ "load the specified binary database into the specified register";
        " (T=twister; t=twisted; d=distance).";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwister; .KPopTwisted; or .KPopDMatrix, respectively)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics ->
          TA.parse_error "You cannot load content into the metrics register"
        | Twister | Twisted | Distances as register_type ->
          Binary_to_register (register_type, TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "T|t|d <table_file_prefix>",
      [ "load the specified tabular database(s) into the specified register";
        " (T=twister; t=twisted; d=distance).";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwister.txt and .KPopInertia.txt; .KPopTwisted.txt;";
        "  or .KPopDMatrix.txt, respectively)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics ->
          TA.parse_error "You cannot load content into the metrics register"
        | Twister | Twisted | Distances as register_type ->
          Tables_to_register (register_type, TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "-a"; "--add" ],
      Some "t|d <binary_file_prefix>",
      [ "add the contents of the specified binary database to the specified register";
        " (t=twisted; d=distance).";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwisted; or .KPopDMatrix, respectively)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Metrics ->
          TA.parse_error "You cannot add content to the twister or metrics registers"
        | Twisted | Distances as register_type ->
          Add_binary_to_register (register_type, TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "-A"; "--Add" ],
      Some "t|d <table_file_prefix>",
      [ "add the contents of the specified tabular database to the specified register";
        " (t=twisted; d=distance).";
        "File extension is automatically determined depending on database type";
        " (will be: .KPopTwisted.txt; or .KPopDMatrix.txt, respectively)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Metrics ->
          TA.parse_error "You cannot add content to the twister or metrics registers"
        | Twisted | Distances as register_type ->
          Add_tables_to_register (register_type, TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "-k"; "-K"; "--kmers"; "--add-kmers"; "--add-kmer-files" ],
      Some "<k-mer_table_file_name>[','...','<k-mer_table_file_name>]",
      [ "twist k-mers from the specified files through the transformation";
        "present in the twister register, and add the results";
        "to the database present in the twisted register" ],
      TA.Optional,
      (fun _ ->
        Add_kmer_files_to_twisted
          (TA.get_parameter () |> Tools.Split.on_char_as_list ',') |> Tools.List.accum Parameters.program);
    [ "--distance"; "--distance-function"; "--set-distance"; "--set-distance-function" ],
      Some "'euclidean'|'cosine'|'minkowski(<non_negative_float>)'",
      [ "set the function to be used when computing distances.";
        "The parameter for Minkowski is the power.";
        "Euclidean is the same as Minkowski(2); Cosine is the same as (Euclidean^2)/2" ],
      TA.Default (fun () -> Space.Distance.to_string Defaults.distance),
      (fun _ -> Set_distance (TA.get_parameter () |> Space.Distance.of_string) |> Tools.List.accum Parameters.program);
    [ "--distance-normalization"; "--set-distance-normalization" ],
      Some "'true'|'false'",
      [ "set whether twisted vectors should be normalized prior to computing distances" ],
      TA.Default (fun () -> string_of_bool Defaults.distance_normalize),
      (fun _ -> Set_distance_normalize (TA.get_parameter_boolean ()) |> Tools.List.accum Parameters.program);
    [ "-m"; "--metric"; "--metric-function"; "--set-metric"; "--set-metric-function" ],
      Some "'flat'|'powers('POWERS_PARAMETERS')'",
      [ "where POWERS_PARAMETERS :=";
        " <non_negative_float>','<fractional_float>','<non_negative_float> :";
        "set the metric function to be used when computing distances.";
        "Parameters are:";
        " internal power; fractional accumulative threshold; external power." ],
      TA.Default (fun () -> Space.Distance.Metric.to_string Defaults.metric),
      (fun _ ->
        Set_metric (TA.get_parameter () |> Space.Distance.Metric.of_string) |> Tools.List.accum Parameters.program);
    [ "-d"; "--distances"; "--compute-distances"; "--compute-twisted-distances" ],
      Some "<twisted_binary_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension .KPopTwisted)";
        "using the metric provided by the twister present in the twister register.";
        "The result will be placed in the distance register" ],
      TA.Optional,
      (fun _ -> Distances_from_twisted_binary (TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "T|t|d <binary_file_prefix>",
      [ "dump the database present in the specified register";
        " (T=twister; t=twisted; d=distance)";
        "to the specified binary file.";
        "File extension is automatically assigned depending on database type";
        " (will be: .KPopTwister; .KPopTwisted; or .KPopDMatrix, respectively)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics ->
          TA.parse_error "You cannot output binary content from the metrics registers"
        | Twister | Twisted | Distances as register_type ->
          Register_to_binary (register_type, TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "--precision"; "--set-precision"; "--set-table-precision" ],
      Some "<positive_integer>",
      [ "set the number of precision digits to be used when outputting numbers" ],
      TA.Default (fun () -> string_of_int Defaults.precision),
      (fun _ -> Set_precision (TA.get_parameter_int_pos ()) |> Tools.List.accum Parameters.program);
    [ "-O"; "--Output" ],
      Some "T|t|d|m <table_file_prefix>",
      [ "dump the database present in the specified register";
        " (T=twister; t=twisted; d=distance; m=metric)";
        "to the specified tabular file(s).";
        "File extension is automatically assigned depending on database type";
        " (will be: .KPopTwister.txt and .KPopInertia.txt; .KPopTwisted.txt;";
        "  .KPopDMatrix.txt; or .KPopMetrics.txt, respectively)" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_tables (register_type, TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "--keep-at-most"; "--set-keep-at-most"; "--summary-keep-at-most" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of closest target sequences";
        "to be kept when summarizing distances.";
        "Note that more might be printed anyway in case of ties" ],
      TA.Default (fun () -> KeepAtMost.to_string Defaults.keep_at_most),
      (fun _ -> Set_keep_at_most (TA.get_parameter () |> KeepAtMost.of_string) |> Tools.List.accum Parameters.program);
    [ "-s"; "--compute-and-summarize-distances"; "--compute-and-summarize-twisted-distances" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension .KPopTwisted)";
        "using the metric provided by the twister present in the twister register;";
        "summarize them and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be .KPopSummary.txt)" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_binary (twisted_prefix, TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    [ "-S"; "--summarize-distances" ],
      Some "<summary_file_prefix>",
      [ "summarize the distances present in the distance register";
        "and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be .KPopSummary.txt)" ],
      TA.Optional,
      (fun _ -> Summary_from_distances (TA.get_parameter ()) |> Tools.List.accum Parameters.program);
    TA.make_separator_multiline [ "Miscellaneous options."; "They are set immediately." ];
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
    [ "-V"; "--version" ],
      None,
      [ "print version and exit" ],
      TA.Optional,
      (fun _ -> Printf.printf "%s\n%!" info.version; exit 0);
    (* Hidden option to emit help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
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
  let twister_loaded = ref false in
  List.iter
    (function
      | Binary_to_register (Twister, _) | Tables_to_register (Twister, _) ->
        twister_loaded := true
      | Register_to_tables (Metrics, _)
      | Distances_from_twisted_binary _ | Summary_from_twisted_binary _ | Summary_from_distances _ ->
        (* A twister must have been loaded to provide the metric induced by inertia *)
        if not !twister_loaded then
          TA.parse_error
            "Options '-O m', '-d', '-s', and '-S' require a twister in the twister register to provide a metric!"
      | Register_to_binary (Metrics, _) ->
        (* This is not really an option *)
        assert false
      | Set_distance _ | Set_distance_normalize _ | Set_metric _ ->
        (* If we really wanted to, we might add these too *)
        ()
      | Empty _
      | Binary_to_register (Metrics, _) | Binary_to_register (Twisted, _) | Binary_to_register (Distances, _)
      | Tables_to_register (Metrics, _) | Tables_to_register (Twisted, _) | Tables_to_register (Distances, _)
      | Add_binary_to_register _ | Add_tables_to_register _ | Add_kmer_files_to_twisted _
      | Register_to_tables (Twister, _) | Register_to_tables (Twisted, _) | Register_to_tables (Distances, _)
      | Register_to_binary (Twister, _) | Register_to_binary (Twisted, _) | Register_to_binary (Distances, _)
      | Set_precision _ | Set_keep_at_most _ ->
        ())
    program;
  (* These are the registers available to the program *)
  let twister = ref Twister.empty and twisted = Matrix.empty Twisted |> ref
  and distance = ref Defaults.distance
  and distance_normalize = ref Defaults.distance_normalize
  and metric = Space.Distance.Metric.compute Defaults.metric |> ref
  and distances = Matrix.empty DMatrix |> ref and precision = ref Defaults.precision
  and keep_at_most = ref Defaults.keep_at_most in
  let compute_metrics () =
    { Matrix.Base.idx_to_row_names = [| "metrics" |];
      idx_to_col_names = !twister.inertia.matrix.idx_to_col_names;
      storage = [| !metric !twister.inertia.matrix.storage.(0) |] } in
  try
    List.iter
      (function
        | Empty RegisterType.Twister ->
          twister := Twister.empty
        | Empty RegisterType.Twisted ->
          twisted := Matrix.empty Twisted
        | Empty RegisterType.Distances ->
          distances := Matrix.empty DMatrix
        | Empty RegisterType.Metrics ->
          assert false
        | Binary_to_register (RegisterType.Twister, prefix) ->
          twister := Twister.of_binary ~verbose:!Parameters.verbose prefix
        | Binary_to_register (RegisterType.Twisted, prefix) ->
          twisted := Matrix.of_binary ~verbose:!Parameters.verbose Twisted prefix;
          if !twisted.which <> Twisted then
            Matrix.Unexpected_type (!twisted.which, Twisted) |> raise
        | Binary_to_register (RegisterType.Distances, prefix) ->
          distances := Matrix.of_binary ~verbose:!Parameters.verbose DMatrix prefix;
          if !distances.which <> DMatrix then
            Matrix.Unexpected_type (!distances.which, DMatrix) |> raise
        | Binary_to_register (RegisterType.Metrics, _) ->
          assert false
        | Add_binary_to_register (RegisterType.Twister, _) ->
          assert false
        | Add_binary_to_register (RegisterType.Twisted, prefix) ->
          let additional = Matrix.of_binary ~verbose:!Parameters.verbose Twisted prefix in
          if additional.which <> Twisted then
            Matrix.Unexpected_type (additional.which, Twisted) |> raise;
          twisted := Matrix.merge_rowwise ~verbose:!Parameters.verbose !twisted additional
        | Add_binary_to_register (RegisterType.Distances, prefix) ->
          let additional = Matrix.of_binary ~verbose:!Parameters.verbose DMatrix prefix in
          if additional.which <> DMatrix then
            Matrix.Unexpected_type (additional.which, DMatrix) |> raise;
          distances := Matrix.merge_rowwise ~verbose:!Parameters.verbose !distances additional
        | Add_binary_to_register (RegisterType.Metrics, _) ->
          assert false
        | Tables_to_register (RegisterType.Twister, prefix) ->
          twister := Twister.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose prefix
        | Tables_to_register (RegisterType.Twisted, prefix) ->
          twisted := Matrix.of_file ~threads:!Parameters.threads ~verbose:!Parameters.verbose Twisted prefix;
          if !twisted.which <> Twisted then
            Matrix.Unexpected_type (!twisted.which, Twisted) |> raise
        | Tables_to_register (RegisterType.Distances, prefix) ->
          distances := Matrix.of_file ~threads:!Parameters.threads ~verbose:!Parameters.verbose DMatrix prefix;
          if !distances.which <> DMatrix then
            Matrix.Unexpected_type (!distances.which, DMatrix) |> raise
        | Tables_to_register (RegisterType.Metrics, _) ->
          assert false
        | Add_tables_to_register (RegisterType.Twister, _) ->
          assert false
        | Add_tables_to_register (RegisterType.Twisted, prefix) ->
          let additional = Matrix.of_file ~threads:!Parameters.threads ~verbose:!Parameters.verbose Twisted prefix in
          if additional.which <> Twisted then
            Matrix.Unexpected_type (additional.which, Twisted) |> raise;
          twisted := Matrix.merge_rowwise ~verbose:!Parameters.verbose !twisted additional
        | Add_tables_to_register (RegisterType.Distances, prefix) ->
          let additional = Matrix.of_file ~threads:!Parameters.threads ~verbose:!Parameters.verbose DMatrix prefix in
          if additional.which <> DMatrix then
            Matrix.Unexpected_type (additional.which, DMatrix) |> raise;
          distances := Matrix.merge_rowwise ~verbose:!Parameters.verbose !distances additional
        | Add_tables_to_register (RegisterType.Metrics, _) ->
          assert false
        | Add_kmer_files_to_twisted fnames ->
          twisted :=
            Twister.add_twisted_from_files
              ~threads:!Parameters.threads ~verbose:!Parameters.verbose !twister !twisted fnames
        | Register_to_binary (RegisterType.Twister, prefix) ->
          Twister.to_binary ~verbose:!Parameters.verbose !twister prefix
        | Register_to_binary (RegisterType.Twisted, prefix) ->
          Matrix.to_binary ~verbose:!Parameters.verbose !twisted prefix
        | Register_to_binary (RegisterType.Distances, prefix) ->
          Matrix.to_binary ~verbose:!Parameters.verbose !distances prefix
        | Register_to_binary (RegisterType.Metrics, _) ->
          assert false
        | Set_precision prec ->
          precision := prec
        | Register_to_tables (RegisterType.Twister, prefix) ->
          Twister.to_files ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !twister prefix
        | Register_to_tables (RegisterType.Twisted, prefix) ->
          Matrix.to_file
            ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose !twisted prefix
        | Register_to_tables (RegisterType.Distances, prefix) ->
          Matrix.to_file
            ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose !distances prefix
        | Register_to_tables (RegisterType.Metrics, prefix) ->
          Matrix.to_file
            ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose {
              which = Metrics;
              matrix = compute_metrics ()
            } prefix
        | Set_distance dist ->
          distance := dist
        | Set_distance_normalize norm ->
          distance_normalize := norm
        | Set_metric metr ->
          metric := Space.Distance.Metric.compute metr
        | Distances_from_twisted_binary prefix ->
          let twisted_db = Matrix.of_binary ~verbose:!Parameters.verbose Twisted prefix in
          if !twisted.which <> Twisted then
            Matrix.Unexpected_type (!twisted.which, Twisted) |> raise;
          distances :=
            Matrix.get_distance_rowwise
              ~normalize:!distance_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !distance (compute_metrics ()).storage.(0) !twisted twisted_db
        | Set_keep_at_most kam ->
          keep_at_most := kam
        | Summary_from_twisted_binary (prefix_in, prefix_out) ->
          let twisted_db = Matrix.of_binary ~verbose:!Parameters.verbose Twisted prefix_in in
          if !twisted.which <> Twisted then
            Matrix.Unexpected_type (!twisted.which, Twisted) |> raise;
          Matrix.summarize_rowwise
            ~keep_at_most:!keep_at_most ~normalize:!distance_normalize
            ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !distance (compute_metrics ()).storage.(0) !twisted twisted_db prefix_out
        | Summary_from_distances prefix ->
          Matrix.summarize_distance
            ~keep_at_most:!keep_at_most ~threads:!Parameters.threads ~verbose:!Parameters.verbose !distances prefix)
      program
  with exc ->
    TA.usage ();
    raise exc

