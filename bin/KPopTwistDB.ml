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
open Better
open KPop

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
    let metric = Space.Distance.Metric.of_string "powers(1,1,2)"
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
  Tools.Argv.name = "KPopTwistDB";
  version = "31";
  date = "16-Jun-2024"
} and authors = [
  "2022-2024", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    [ "-e"; "--empty" ],
      Some "'T'|'t'|'d'",
      [ "load an empty database into the specified register";
        " (T=twister; t=twisted; d=distance)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics ->
          TA.parse_error "You cannot load content into the metrics register"
        | Twister | Twisted | Distances as register_type ->
          Empty register_type |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "'T'|'t'|'d' <binary_file_prefix>",
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
          Binary_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "'T'|'t'|'d' <table_file_prefix>",
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
          Tables_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-a"; "--add" ],
      Some "'t'|'d' <binary_file_prefix>",
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
          Add_binary_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-A"; "--Add" ],
      Some "'t'|'d' <table_file_prefix>",
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
          Add_tables_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-k"; "--kmers"; "--add-kmers"; "--add-kmer-files" ],
      Some "<k-mer_table_file_name>[','...','<k-mer_table_file_name>]",
      [ "twist k-mers from the specified files through the transformation";
        "present in the twister register, and add the results";
        "to the database present in the twisted register" ],
      TA.Optional,
      (fun _ ->
        Add_kmer_files_to_twisted
          (TA.get_parameter () |> String.Split.on_char_as_list ',') |> List.accum Parameters.program);
    [ "--distance"; "--distance-function"; "--set-distance"; "--set-distance-function" ],
      Some "'euclidean'|'cosine'|'minkowski(<non_negative_float>)'",
      [ "set the function to be used when computing distances.";
        "The parameter for 'minkowski' is the power.";
        "Note that 'euclidean' is the same as 'minkowski(2)',";
        "and 'cosine' is the same as ('euclidean'^2)/2" ],
      TA.Default (fun () -> Space.Distance.to_string Defaults.distance),
      (fun _ -> Set_distance (TA.get_parameter () |> Space.Distance.of_string) |> List.accum Parameters.program);
    [ "--distance-normalization"; "--set-distance-normalization" ],
      Some "'true'|'false'",
      [ "set whether twisted vectors should be normalized prior to computing distances" ],
      TA.Default (fun () -> string_of_bool Defaults.distance_normalize),
      (fun _ -> Set_distance_normalize (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "-m"; "--metric"; "--metric-function"; "--set-metric"; "--set-metric-function" ],
      Some "'flat'|'powers('POWERS_PARAMETERS')'",
      [ "where POWERS_PARAMETERS :=";
        " <non_negative_float>','<fractional_float>','<non_negative_float> :";
        "set the metric function to be used when computing distances.";
        "Parameters are:";
        " internal power; fractional accumulative threshold; external power." ],
      TA.Default (fun () -> Space.Distance.Metric.to_string Defaults.metric),
      (fun _ ->
        Set_metric (TA.get_parameter () |> Space.Distance.Metric.of_string) |> List.accum Parameters.program);
    [ "-d"; "--distances"; "--compute-distances"; "--compute-twisted-distances" ],
      Some "<twisted_binary_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension .KPopTwisted)";
        "using the metric provided by the twister present in the twister register.";
        "The result will be placed in the distance register" ],
      TA.Optional,
      (fun _ -> Distances_from_twisted_binary (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "'T'|'t'|'d' <binary_file_prefix>",
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
          Register_to_binary (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--precision"; "--set-precision"; "--set-table-precision" ],
      Some "<positive_integer>",
      [ "set the number of precision digits to be used when outputting numbers" ],
      TA.Default (fun () -> string_of_int Defaults.precision),
      (fun _ -> Set_precision (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-O"; "--Output" ],
      Some "'T'|'t'|'d'|'m' <table_file_prefix>",
      [ "dump the database present in the specified register";
        " (T=twister; t=twisted; d=distance; m=metric)";
        "to the specified tabular file(s).";
        "File extension is automatically assigned depending on database type";
        " (will be: .KPopTwister.txt and .KPopInertia.txt; .KPopTwisted.txt;";
        "  .KPopDMatrix.txt; or .KPopMetrics.txt, respectively)" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_tables (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-K"; "--keep-at-most"; "--set-keep-at-most"; "--summary-keep-at-most" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of closest target sequences";
        "to be kept when summarizing distances.";
        "Note that more might be printed anyway in case of ties" ],
      TA.Default (fun () -> KeepAtMost.to_string Defaults.keep_at_most),
      (fun _ -> Set_keep_at_most (TA.get_parameter () |> KeepAtMost.of_string) |> List.accum Parameters.program);
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
        Summary_from_twisted_binary (twisted_prefix, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-S"; "--summarize-distances"; "--summarize-twisted-distances" ],
      Some "<summary_file_prefix>",
      [ "summarize the distances present in the distance register";
        "and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be .KPopSummary.txt)" ],
      TA.Optional,
      (fun _ -> Summary_from_distances (TA.get_parameter ()) |> List.accum Parameters.program);
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
      | Register_to_tables (Metrics, _) | Distances_from_twisted_binary _ | Summary_from_twisted_binary _ ->
        (* A twister must have been loaded to provide the metric induced by inertia *)
        if not !twister_loaded then
          TA.parse_error
            "Options '-O m', '-d', and '-s' require a twister in the twister register to provide a metric!"
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
      | Set_precision _ | Set_keep_at_most _ | Summary_from_distances _ ->
        ())
    program;
  (* These are the registers available to the program *)
  let twister = ref Twister.empty and twisted = Matrix.empty Twisted |> ref
  and distance = ref Defaults.distance and distance_normalize = ref Defaults.distance_normalize
  and metric = Defaults.metric |> ref and distances = Matrix.empty DMatrix |> ref
  and precision = ref Defaults.precision and keep_at_most = ref Defaults.keep_at_most in
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
            ~precision:!precision ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            (Twister.get_metrics_matrix !metric !twister) prefix
        | Set_distance dist ->
          distance := dist
        | Set_distance_normalize norm ->
          distance_normalize := norm
        | Set_metric metr ->
          metric := metr
        | Distances_from_twisted_binary prefix ->
          let twisted_db = Matrix.of_binary ~verbose:!Parameters.verbose Twisted prefix in
          if !twisted.which <> Twisted then
            Matrix.Unexpected_type (!twisted.which, Twisted) |> raise;
          distances :=
            Matrix.get_distance_rowwise
              ~normalize:!distance_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !distance (Twister.get_metrics_vector !metric !twister) !twisted twisted_db
        | Set_keep_at_most kam ->
          keep_at_most := kam
        | Summary_from_twisted_binary (prefix_in, prefix_out) ->
          let twisted_db = Matrix.of_binary ~verbose:!Parameters.verbose Twisted prefix_in in
          if !twisted.which <> Twisted then
            Matrix.Unexpected_type (!twisted.which, Twisted) |> raise;
          Matrix.summarize_rowwise
            ~keep_at_most:!keep_at_most ~normalize:!distance_normalize
            ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !distance (Twister.get_metrics_vector !metric !twister) !twisted twisted_db prefix_out
        | Summary_from_distances prefix ->
          Matrix.summarize_distance
            ~keep_at_most:!keep_at_most ~threads:!Parameters.threads ~verbose:!Parameters.verbose !distances prefix)
      program
  with exc ->
    Printf.eprintf "[#%s]: (%s): %s\n%!" (Unix.getpid () |> string_of_int |> String.TermIO.blue) __FUNCTION__
      ("FATAL: Uncaught exception: " ^ Printexc.to_string exc |> String.TermIO.red)

