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
      | Metrics
      | Twister
      | Twisted
      | Embeddings
      | Distances
      | Splits
    exception Invalid_register_type of string
    let of_string = function
      | "m" -> Metrics
      | "T" -> Twister
      | "t" -> Twisted
      | "e" -> Embeddings
      | "d" -> Distances
      | "s" -> Splits
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

type to_do_t =
  | Empty of RegisterType.t
  | Binary_to_register of RegisterType.t * string
  | Tables_to_register of RegisterType.t * string
  | Add_binary_to_register of RegisterType.t * string
  | Add_tables_to_register of RegisterType.t * string
  | Set_kmers_normalize of bool
  | Add_kmers_files_to_twisted of string list
(* | Add_kmers_binary_to_twisted of string *)
  | Register_to_binary of RegisterType.t * string
  | Set_precision_tables of int
  | Set_precision_splits of int
  | Register_to_tables of RegisterType.t * string
  | Set_distance of Space.Distance.t
  | Set_distance_normalize of bool
  | Set_metric of Space.Distance.Metric.t
  | Embeddings_from_twisted
  | Set_splits_algorithm of Matrix.SplitsAlgorithm.t
  | Set_splits_keep_at_most of int
  | Splits_from_embeddings
  | Distances_from_twisted_binary of string
  | Set_summary_keep_at_most of KeepAtMost.t
  | Summary_from_twisted_binary of string * string
  | Summary_from_distances of string

module Defaults =
  struct
    let kmers_normalize = true
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalize = true
    let metric = Space.Distance.Metric.of_string "powers(1,1,2)"
    let precision_tables = 15
    let precision_splits = 10
    let splits_algorithm = Matrix.SplitsAlgorithm.of_string "gaps"
    let splits_keep_at_most = 10000
    let summary_keep_at_most = Some 2
  end

module Parameters =
  struct
    let program = ref []
    let threads = Processes.Parallel.get_nproc () |> ref
    let debug_twisting = ref false
    let verbose = ref false
  end

let info = {
  Tools.Argv.name = "KPopTwistDB";
  version = "38";
  date = "10-Feb-2025"
} and authors = [
  "2022-2025", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    [ "-z"; "--zero"; "--empty" ],
      Some "'T'|'t'|'e'|'d'",
      [ "load an empty database into the specified register";
        " ('T'=twister; 't'=twisted; 'e'=embeddings; 'd'=distances)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics | Splits ->
          TA.parse_error "You cannot load content into the metric or splits registers"
        | Twister | Twisted | Embeddings | Distances as register_type ->
          Empty register_type |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "'T'|'t'|'e'|'d' <binary_file_prefix>",
      [ "load the specified binary database into the specified register";
        " ('T'=twister; 't'=twisted; 'e'=embeddings; 'd'=distances).";
        "File extension is automatically determined depending on database type";
        " (will be: '.KPopTwister'; '.KPopTwisted'; '.KPopVectors';";
        "  or '.KPopDMatrix', respectively, unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics | Splits ->
          TA.parse_error "You cannot load content into the metric or splits registers"
        | Twister | Twisted | Embeddings | Distances as register_type ->
          Binary_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "'T'|'t'|'e'|'d' <table_file_prefix>",
      [ "load the specified tabular database(s) into the specified register";
        " ('T'=twister; 't'=twisted; 'e'=embeddings; 'd'=distances).";
        "File extension is automatically determined depending on database type";
        " (will be: '.KPopTwister.txt' and '.KPopInertia.txt'; '.KPopTwisted.txt';";
        "  '.KPopVectors.txt'; or '.KPopDMatrix.txt', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics | Splits ->
          TA.parse_error "You cannot load content into the metric or splits registers"
        | Twister | Twisted | Embeddings | Distances as register_type ->
          Tables_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-a"; "--add" ],
      Some "'t'|'e'|'d' <binary_file_prefix>",
      [ "add the contents of the specified binary database to the specified register";
        " ('t'=twisted; 'e'=embeddings; 'd'=distances).";
        "File extension is automatically determined depending on database type";
        " (will be: '.KPopTwisted'; '.KPopVectors'; or '.KPopDMatrix', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Metrics | Splits ->
          TA.parse_error "You cannot add content to the twister, metric or splits registers"
        | Twisted | Embeddings | Distances as register_type ->
          Add_binary_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-A"; "--Add" ],
      Some "'t'|'e'|'d'|'s' <table_file_prefix>",
      [ "add the contents of the specified tabular database to the specified register";
        " ('t'=twisted; 'e'=embeddings; 'd'=distances).";
        "File extension is automatically determined depending on database type";
        " (will be: '.KPopTwisted.txt'; '.KPopVectors'; or '.KPopDMatrix.txt',";
        "  respectively, unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Metrics | Splits ->
          TA.parse_error "You cannot add content to the twister, metric or splits registers"
        | Twisted | Embeddings | Distances as register_type ->
          Add_tables_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--counts-normalize"; "--counts-normalization" ],
      Some "'true'|'false'",
      [ "whether to normalize spectra before twisting" ],
      TA.Default (fun () -> string_of_bool Defaults.kmers_normalize),
      (fun _ -> Set_kmers_normalize (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "-k"; "--kmers"; "--add-kmers"; "--add-kmer-files" ],
      Some "<k-mer_table_file_name>[','...','<k-mer_table_file_name>]",
      [ "twist k-mers from the specified files through the transformation";
        "present in the twister register, and add the results";
        "to the database present in the twisted register" ],
      TA.Optional,
      (fun _ ->
        Add_kmers_files_to_twisted
          (TA.get_parameter () |> String.Split.on_char_as_list ',') |> List.accum Parameters.program);
    [ "--distance"; "--distance-function" ],
      Some "'euclidean'|'cosine'|'minkowski(<non_negative_float>)'",
      [ "set the function to be used when computing distances.";
        "The parameter for 'minkowski' is the power.";
        "Note that 'euclidean' is the same as 'minkowski(2)',";
        "and 'cosine' is the same as ('euclidean'^2)/2" ],
      TA.Default (fun () -> Space.Distance.to_string Defaults.distance),
      (fun _ -> Set_distance (TA.get_parameter () |> Space.Distance.of_string) |> List.accum Parameters.program);
    [ "--distance-normalize"; "--distance-normalization" ],
      Some "'true'|'false'",
      [ "whether to normalize twisted vectors before computing distances" ],
      TA.Default (fun () -> string_of_bool Defaults.distance_normalize),
      (fun _ -> Set_distance_normalize (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "-m"; "--metric"; "--metric-function" ],
      Some "'flat'|'powers('POWERS_PARAMETERS')'",
      [ "where POWERS_PARAMETERS :=";
        " <non_negative_float>','<fractional_float>','<non_negative_float> :";
        "set the metric function to be used when computing distances.";
        "Parameters are:";
        " internal power; fractional accumulative threshold; external power." ],
      TA.Default (fun () -> Space.Distance.Metric.to_string Defaults.metric),
      (fun _ ->
        Set_metric (TA.get_parameter () |> Space.Distance.Metric.of_string) |> List.accum Parameters.program);
    [ "-e"; "--embeddings"; "--compute-embeddings"; "--twisted-to-embeddings" ],
      None,
      [ "compute embeddings from the vectors present in the twisted register";
        "using the metric provided by the twister present in the twister register.";
        "The result will be placed in the embeddings register" ],
      TA.Optional,
      (fun _ -> List.accum Parameters.program Embeddings_from_twisted);
    [ "--splits-algorithm" ],
      Some "'gaps'|'centroids'",
      [ "algorithm to use when computing splits from embeddings" ],
      TA.Default (fun () -> Matrix.SplitsAlgorithm.to_string Defaults.splits_algorithm),
      (fun _ ->
        Set_splits_algorithm (TA.get_parameter () |> Matrix.SplitsAlgorithm.of_string)
        |> List.accum Parameters.program);
    [ "--splits-at-most"; "--splits-keep-at-most" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of phylogenetic splits to be kept";
        "when generating them from embeddings" ],
      TA.Default (fun () -> string_of_int Defaults.splits_keep_at_most),
      (fun _ -> Set_splits_keep_at_most (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-p"; "--splits"; "--compute-splits"; "--embeddings-to-splits" ],
      None,
      [ "compute phylogenetic splits from the vectors present in the embeddings";
        "register. The result will be placed in the splits register" ],
      TA.Optional,
      (fun _ -> List.accum Parameters.program Splits_from_embeddings);
    [ "-d"; "--distances"; "--compute-distances"; "--compute-twisted-distances" ],
      Some "<twisted_binary_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension '.KPopTwisted' unless file is '/dev/*')";
        "using the metric provided by the twister present in the twister register.";
        "The result will be placed in the distance register" ],
      TA.Optional,
      (fun _ -> Distances_from_twisted_binary (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "'T'|'t'|'e'|'d'|'s' <binary_file_prefix>",
      [ "dump the database present in the specified register";
        " ('T'=twister; 't'=twisted; 'e'=embeddings; 'd'=distances; 's'=splits)";
        "to the specified binary file.";
        "File extension is automatically assigned depending on database type";
        " (will be: '.KPopTwister'; '.KPopTwisted'; '.KPopVectors'; .KPopDMatrix;";
        "  or '.PhyloSplits', respectively, unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Metrics ->
          TA.parse_error "You cannot output binary content from the metrics registers"
        | Twister | Twisted | Embeddings | Distances | Splits as register_type ->
          Register_to_binary (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--precision-for-tables" ],
      Some "<positive_integer>",
      [ "set how many precision digits should be used when outputting numbers";
        "in tabular formats" ],
      TA.Default (fun () -> string_of_int Defaults.precision_tables),
      (fun _ -> Set_precision_tables (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "--precision-for-splits" ],
      Some "<positive_integer>",
      [ "set how many precision digits should be used when outputting splits";
        "in plain-text format" ],
      TA.Default (fun () -> string_of_int Defaults.precision_splits),
      (fun _ -> Set_precision_splits (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-O"; "--Output" ],
      Some "'T'|'t'|'e'|'d'|'m'|'s' <table_file_prefix>",
      [ "dump the database present in the specified register";
        " ('T'=twister; 't'=twisted; 'e'=embeddings; 'd'=distances; 'm'=metric;";
        "  's'=splits)";
        "to the specified tabular file(s).";
        "File extension is automatically assigned depending on database type";
        " (will be: '.KPopTwister.txt' and '.KPopInertia.txt'; '.KPopTwisted.txt';";
        "  '.KPopVectors.txt'; '.KPopDMatrix.txt'; '.KPopMetrics.txt';";
        "  or '.PhyloSplits.txt', respectively, unless files is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_tables (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--summary-at-most"; "--summary-keep-at-most" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of closest target sequences to be kept";
        "when summarizing distances. Note that more might be printed anyway";
        "in case of ties" ],
      TA.Default (fun () -> KeepAtMost.to_string Defaults.summary_keep_at_most),
      (fun _ ->
        Set_summary_keep_at_most (TA.get_parameter () |> KeepAtMost.of_string) |> List.accum Parameters.program);
    [ "-s"; "--compute-and-summarize-distances"; "--compute-and-summarize-twisted-distances" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "compute distances between all the vectors present in the twisted register";
        "and all the vectors present in the specified twisted binary file";
        " (which must have extension '.KPopTwisted' unless file is '/dev/*')";
        "using the metric provided by the twister present in the twister register;";
        "summarize them and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopSummary.txt' unless files is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_binary (twisted_prefix, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-S"; "--summarize-distances"; "--summarize-twisted-distances" ],
      Some "<summary_file_prefix>",
      [ "summarize the distances present in the distance register";
        "and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopSummary.txt' unless file is '/dev/*')" ],
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
    (* Hidden option to profile twisting *)
    [ "--debug-twisting" ], None, [], TA.Optional, (fun _ -> Parameters.debug_twisting := true);
    (* Hidden option to emit help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    (* Hidden option to print exception backtrace *)
    [ "-x"; "--print-exception-backtrace" ], None, [], TA.Optional, (fun _ -> Printexc.record_backtrace true);
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
      | Add_kmers_files_to_twisted _ ->
        (* A twister must have been loaded to twist spectra *)
        if not !twister_loaded then
          TA.parse_error
            "Option '-k' requires a twister in the twister register!"
      | Register_to_tables (Metrics, _) | Embeddings_from_twisted
      | Distances_from_twisted_binary _ | Summary_from_twisted_binary _ ->
        (* A twister must have been loaded to provide the metric induced by inertia *)
        if not !twister_loaded then
          TA.parse_error
            "Options '-O m', '-e', '-d', and '-s' require a twister in the twister register to provide a metric!"
      | Register_to_binary (Metrics, _) ->
        (* This is not really an option *)
        assert false
      | Set_distance _ | Set_distance_normalize _ | Set_metric _ ->
        (* If we really wanted to, we might add these too *)
        ()
      | Empty _
      | Binary_to_register (Metrics, _) | Binary_to_register (Twisted, _)
      | Binary_to_register (Embeddings, _) | Binary_to_register (Distances, _)
      | Binary_to_register (Splits, _)
      | Tables_to_register (Metrics, _) | Tables_to_register (Twisted, _)
      | Tables_to_register (Embeddings, _) | Tables_to_register (Distances, _)
      | Tables_to_register (Splits, _)
      | Add_binary_to_register _ | Add_tables_to_register _ | Set_kmers_normalize _
      | Register_to_tables (Twister, _) | Register_to_tables (Twisted, _)
      | Register_to_tables (Embeddings, _) | Register_to_tables (Distances, _)
      | Register_to_tables (Splits, _)
      | Register_to_binary (Twister, _) | Register_to_binary (Twisted, _)
      | Register_to_binary (Embeddings, _) | Register_to_binary (Distances, _)
      | Register_to_binary (Splits, _)
      | Set_precision_tables _ | Set_precision_splits _
      | Set_splits_algorithm _ | Set_splits_keep_at_most _ | Splits_from_embeddings
      | Set_summary_keep_at_most _ | Summary_from_distances _ ->
        ())
    program;
  (* These are the registers available to the program *)
  let twister = ref Twister.empty and twisted = Matrix.empty Twisted |> ref
  and embeddings = Matrix.empty Vectors |> ref and metric = Defaults.metric |> ref
  and kmers_normalize = ref Defaults.kmers_normalize and distance = ref Defaults.distance
  and distance_normalize = ref Defaults.distance_normalize and distances = Matrix.empty DMatrix |> ref
  and splits_keep_at_most = ref Defaults.splits_keep_at_most
  and splits_algorithm = ref Defaults.splits_algorithm and splits = Trees.Splits.create [||] |> ref
  and summary_keep_at_most = ref Defaults.summary_keep_at_most
  and precision_tables = ref Defaults.precision_tables and precision_splits = ref Defaults.precision_splits in
  let open_and_check f ty prefix =
    let res = f ?verbose:(Some !Parameters.verbose) ty prefix in
    if res.Matrix.which <> ty then
      Matrix.Unexpected_type (res.which, ty) |> raise;
    res
  and open_and_check_and_merge f ty prefix current =
    let additional = f ?verbose:(Some !Parameters.verbose) ty prefix in
    if additional.Matrix.which <> ty then
      Matrix.Unexpected_type (additional.which, ty) |> raise;
    Matrix.merge_rowwise ~verbose:!Parameters.verbose current additional in
  let open_binary_and_check = open_and_check Matrix.of_binary
  and open_binary_and_check_and_merge = open_and_check_and_merge Matrix.of_binary
  and open_file_and_check = open_and_check (Matrix.of_file ~threads:!Parameters.threads ~bytes_per_step:4194304)
  and open_file_and_check_and_merge =
    open_and_check_and_merge (Matrix.of_file ~threads:!Parameters.threads ~bytes_per_step:4194304)
  and matrix_to_binary = Matrix.to_binary ~verbose:!Parameters.verbose
  and matrix_to_file =
    Matrix.to_file ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose in
  try
    List.iter
      (function
        | Empty Metrics | Empty Splits ->
          assert false
        | Empty Twister ->
          twister := Twister.empty
        | Empty Twisted ->
          twisted := Matrix.empty Twisted
        | Empty Embeddings ->
          embeddings := Matrix.empty Vectors
        | Empty Distances ->
          distances := Matrix.empty DMatrix
        | Binary_to_register (Twister, prefix) ->
          twister := Twister.of_binary ~verbose:!Parameters.verbose prefix
        | Binary_to_register (Twisted, prefix) ->
          twisted := open_binary_and_check Twisted prefix
        | Binary_to_register (Embeddings, prefix) ->
          embeddings := open_binary_and_check Vectors prefix
        | Binary_to_register (Distances, prefix) ->
          distances := open_binary_and_check DMatrix prefix
        | Binary_to_register (Metrics, _) | Binary_to_register (Splits, _) ->
          assert false
        | Add_binary_to_register (Twister, _) | Add_binary_to_register (Metrics, _)
        | Add_binary_to_register (Splits, _) ->
          assert false
        | Add_binary_to_register (Twisted, prefix) ->
          twisted := open_binary_and_check_and_merge Twisted prefix !twisted
        | Add_binary_to_register (Embeddings, prefix) ->
          embeddings := open_binary_and_check_and_merge Vectors prefix !embeddings
        | Add_binary_to_register (Distances, prefix) ->
          distances := open_binary_and_check_and_merge DMatrix prefix !distances
        | Tables_to_register (Metrics, _) | Tables_to_register (Splits, _) ->
          assert false
        | Tables_to_register (Twister, prefix) ->
          twister := Twister.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose prefix
        | Tables_to_register (Twisted, prefix) ->
          twisted := open_file_and_check Twisted prefix
        | Tables_to_register (Embeddings, prefix) ->
          twisted := open_file_and_check Vectors prefix
        | Tables_to_register (Distances, prefix) ->
          distances := open_file_and_check DMatrix prefix
        | Add_tables_to_register (Twister, _) | Add_tables_to_register (Metrics, _)
        | Add_tables_to_register (Splits, _) ->
          assert false
        | Add_tables_to_register (Twisted, prefix) ->
          twisted := open_file_and_check_and_merge Twisted prefix !twisted
        | Add_tables_to_register (Embeddings, prefix) ->
          embeddings := open_file_and_check_and_merge Vectors prefix !embeddings
        | Add_tables_to_register (Distances, prefix) ->
          distances := open_file_and_check_and_merge DMatrix prefix !distances
        | Set_kmers_normalize norm ->
          kmers_normalize := norm
        | Add_kmers_files_to_twisted fnames ->
          twisted :=
            Twister.add_twisted_from_files
              ~normalize:!kmers_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              ~debug:!Parameters.debug_twisting !twister !twisted fnames
        | Embeddings_from_twisted ->
          embeddings :=
            Matrix.get_embeddings
              ~normalize:!distance_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !distance (Twister.get_metrics_vector !metric !twister) !twisted
        | Set_splits_algorithm algo ->
          splits_algorithm := algo
        | Set_splits_keep_at_most kam ->
          splits_keep_at_most := kam
        | Splits_from_embeddings ->
          splits :=
            Matrix.get_splits ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !splits_algorithm !splits_keep_at_most !embeddings
        | Register_to_binary (Metrics, _) ->
          assert false
        | Register_to_binary (Twister, prefix) ->
          Twister.to_binary ~verbose:!Parameters.verbose !twister prefix
        | Register_to_binary (Twisted, prefix) ->
          matrix_to_binary !twisted prefix
        | Register_to_binary (Embeddings, prefix) ->
          matrix_to_binary !embeddings prefix
        | Register_to_binary (Distances, prefix) ->
          matrix_to_binary !distances prefix
        | Register_to_binary (Splits, prefix) ->
          Trees.Splits.to_binary ~verbose:!Parameters.verbose !splits prefix
        | Set_precision_tables prec ->
          precision_tables := prec
        | Set_precision_splits prec ->
          precision_splits := prec
        | Register_to_tables (Twister, prefix) ->
          Twister.to_files ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !twister prefix
        | Register_to_tables (Twisted, prefix) ->
          matrix_to_file !twisted prefix
        | Register_to_tables (Embeddings, prefix) ->
          matrix_to_file !embeddings prefix
        | Register_to_tables (Distances, prefix) ->
          matrix_to_file !distances prefix
        | Register_to_tables (Metrics, prefix) ->
          matrix_to_file (Twister.get_metrics_matrix !metric !twister) prefix
        | Register_to_tables (Splits, prefix) ->
          Trees.Splits.to_file ~precision:!precision_splits !splits prefix
        | Set_distance dist ->
          distance := dist
        | Set_distance_normalize norm ->
          distance_normalize := norm
        | Set_metric metr ->
          metric := metr
        | Distances_from_twisted_binary prefix ->
          distances :=
            Matrix.get_distance_rowwise
              ~normalize:!distance_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !distance (Twister.get_metrics_vector !metric !twister) !twisted (open_binary_and_check Twisted prefix)
        | Set_summary_keep_at_most kam ->
          summary_keep_at_most := kam
        | Summary_from_twisted_binary (prefix_in, prefix_out) ->
          Matrix.summarize_rowwise
            ~keep_at_most:!summary_keep_at_most ~normalize:!distance_normalize
            ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !distance (Twister.get_metrics_vector !metric !twister)
            !twisted (open_binary_and_check Twisted prefix_in) prefix_out
        | Summary_from_distances prefix ->
          Matrix.summarize_distance
            ~keep_at_most:!summary_keep_at_most ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !distances prefix)
      program
  with exc ->
    Printf.peprintf "(%s): %s\n%!" __FUNCTION__
      ("FATAL: Uncaught exception: " ^ Printexc.to_string exc |> String.TermIO.red);
    Printf.peprintf "(%s): This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!" __FUNCTION__;
    Printf.peprintf "(%s): You might also wish to rerun me with option -x to get a full backtrace.\n%!" __FUNCTION__;
    Printexc.print_backtrace stderr
