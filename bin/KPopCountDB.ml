(*
    KPopCountDB.ml -- (c) 2022-2025 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    KPopCountDB is the binary for manipulating, annotating, and
    inspecting binary `.KPopSpectra` databases.  Operations include
    merging spectra, filtering by metadata, attaching class labels,
    tabular import/export, and k-mer distillation (sorting k-mers based on
    their information content against a class label set).

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

type to_do_t =
  | Empty
  | Of_binary of string
  | Of_tables of string
  | Set_merge_criterion of KMerDB.MergeCriterion.t
  | Merge_binary of string
  | Add_metadata of string
  | Set_combination_criterion of KMerDB.CombinationCriterion.t
  | Split_spectra of string
  | Add_combined_selected of string (* The new label *)
  | Remove_selected
  | Transform of KMerDB.Transformation.t
  | Distill_kmers of string * string
  | Summary
  | Selected_from_labels of StringSet.t
  | Selected_from_regexps of regexps_t
  | Selected_negate
  | Selected_print
  | Selected_clear
  | To_binary of string
  | Set_output_zero_kmers of bool
  | Set_precision of int
  | To_tables of string
  | Distance_set of Space.Distance.t
  | Distance_normalisation_set of bool
  | To_distances of regexps_t * regexps_t * string
and regexps_t = (string * Str.regexp) list

module Defaults =
  struct
    let merge_criterion = KMerDB.MergeCriterion.of_string "union"
    let combination_criterion = KMerDB.CombinationCriterion.of_string "mean"
    let transformation = KMerDB.Transformation.of_string "power(1)"
    let output_zero_kmers = true
    let precision = 15
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalise = true
    let threads = Processes.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let program = ref []
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let info = {
  Tools.Argv.name = "KPopCountDB";
  version = "55";
  date = "10-Dec-2025"
} and authors = [
  "2020-2025", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[ACTIONS]";
  let parse_regexp_selector option s =
    List.map
    (fun l ->
      let res = String.Split.on_char_as_list '~' l in
      if List.length res <> 2 then begin
        TA.usage ();
        List.length res |>
          Printf.sprintf "Option '%s': Wrong number of fields in list (expected 2, found %d)" option |>
          TA.parse_error (* parse_error exits the program *)
      end;
      List.nth res 0, List.nth res 1 |> Str.regexp)
    (String.Split.on_char_as_list ',' s) in
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    TA.make_separator_multiline [ ""; "Actions on the database register - Input/Output operations:" ];
    [ "-0"; "--zero"; "--empty" ],
      None,
      [ "load an empty database into the register" ],
      TA.Optional,
      (fun _ -> Empty |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load into the register the database present in the specified binary file";
        " (which must have extension '.KPopSpectra' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Of_binary (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "<tabular_file_prefix>",
      [ "load into the register the database present in the specified tabular files";
        " (which must have extensions '.KPopKMatrix.txt' and '.KPopMMatrix.txt'";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Of_tables (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--addition-criterion"; "--database-addition-criterion" ],
      Some "'first'|'second'|'union'|'intersect'",
      [ "set the criterion used to combine databases of spectra.";
        "Possibilities are:";
        " 'first'";
        "  The resulting database will have the same k-mer and metadata labels";
        "  as the first database";
        " 'second'";
        "  The resulting database will have the same k-mer and metadata labels";
        "  as the second database";
        " 'union'";
        "  The resulting database will have as k-mer and metadata labels the";
        "  union of the k-mer and metadata labels of the two databases";
        " 'intersect'";
        "  The resulting database will have as k-mer and metadata labels the";
        "  intersection of the k-mer and metadata labels of the two databases.";
        "For criteria 'first', 'second', and 'union', missing data will be set";
        "to zero in the case of k-mer counts, and to the empty string in the case";
        "of metadata entries" ],
      TA.Default (KMerDB.MergeCriterion.to_string Defaults.merge_criterion |> Fun.const),
      (fun _ ->
        Set_merge_criterion (TA.get_parameter () |> KMerDB.MergeCriterion.of_string)
          |> List.accum Parameters.program);
    [ "-a"; "--add"; "--add-database" ],
      Some "<binary_file_prefix>",
      [ "add to the register the database present in the specified binary file";
        " (which must have extension '.KPopSpectra' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Merge_binary (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-m"; "--metadata"; "--add-metadata" ],
      Some "<metadata_table_file_name>",
      [ "add to the register metadata from the specified tabular file.";
        "Metadata should be presented as a tab-separated text table, with a header";
        "containing spectrum labels and with row names being metadata fields labels.";
        "Spectrum labels, metadata field names and metadata values must not contain";
        "double quote '\"' characters" ],
      TA.Optional,
      (fun _ -> Add_metadata (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--summary" ],
      None,
      [ "print a summary of the database present in the register" ],
      TA.Optional,
      (fun _ -> Summary |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "save the database present in the register to the specified file";
        " (which will be given extension '.KPopSpectra' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> To_binary (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--precision" ],
      Some "<positive_integer>",
      [ "set the number of precision digits to be used for tabular output" ],
      TA.Default (string_of_int Defaults.precision |> Fun.const),
      (fun _ -> Set_precision (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "--Output-zero-kmers"; "--Output-zero-k-mers" ],
      Some "'true'|'false'",
      [ "whether to output k-mers whose frequencies are all zero";
        "when writing the database as tabular files" ],
      TA.Default (string_of_bool Defaults.output_zero_kmers |> Fun.const),
      (fun _ -> Set_output_zero_kmers (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "-O"; "--Output" ],
      Some "<tabular_file_prefix>",
      [ "write the database present in the register as tab-separated files.";
        "File extensions are automatically assigned";
        " (will be '.KPopKMatrix.txt' and '.KPopMMatrix.txt',";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> To_tables (TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Actions on the database register - Other operations:" ];
    [ "--combination-criterion"; "--spectrum-combination-criterion" ],
      Some "'mean'|'median'",
      [ "set the criterion used to combine the k-mer frequencies of spectra.";
        "To avoid rounding issues, each k-mer frequency is also rescaled";
        "by the largest normalization across spectra";
        " ('mean' averages frequencies across spectra;";
        "  'median' computes the median across spectra)" ],
      TA.Default (KMerDB.CombinationCriterion.to_string Defaults.combination_criterion |> Fun.const),
      (fun _ ->
        Set_combination_criterion (TA.get_parameter () |> KMerDB.CombinationCriterion.of_string)
          |> List.accum Parameters.program);
    [ "-c"; "--combine"; "--combine-by-class"; "--combine-spectra-by-class" ],
      Some "<classes_metadata_field_name>",
      [ "split the database into classes according to the labels contained in the";
        "specified metadata field and combine the spectra belonging to each class";
        "into a separate vector named as the class label. Delete original spectra.";
        "Class label cannot be the same as the name of an existing spectrum" ],
      TA.Optional,
      (fun _ -> Split_spectra (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--transform" ],
      Some "TRANSFORMATION",
      [ "replace the database with the one obtained from the specified transformation.";
        "Transformations are defined as follows:";
        " TRANSFORMATION := 'threshold(' <non-negative_float>')'";
        "                 | 'power(<float>)'";
        "                 | 'binary'";
        "                 | 'clr'";
        "                 | 'pseudocounts('POWER','QUANTIZE')'";
        " POWER := <non-negative_float>";
        " QUANTIZE := 'false'|'true'";
        "A value such that 0. <= THRESHOLD < 1. is interpreted as a fraction relative";
        "to the sum of all the counts in the spectrum; values such that THRESHOLD >= 1.";
        "are considered absolute thresholds; 'binary' is an alias for 'power(0)'.";
        "For the exact definition of transformations 'clr' and 'pseudocounts', see";
        " https://github.com/PaoloRibeca/KPop";
        "or";
        " https://doi.org/10.1186/s13059-025-03585-8" ],
      TA.Default (KMerDB.Transformation.to_string Defaults.transformation |> Fun.const),
      (fun _ ->
        Transform (TA.get_parameter () |> KMerDB.Transformation.of_string) |> List.accum Parameters.program);
    [ "-d"; "--distill"; "--distill-kmers" ],
      Some "<classes_metadata_field_name> <summary_file_prefix>",
      [ "optimize k-mers by identifying which ones are most informative";
        "according to the labels contained in the specified metadata field";
        "and by re-sorting k-mers in decreasing order accordingly.";
        "The labels must identify at least two equivalence classes, and fewer classes";
        "than the number of k-mers.";
        "Details of the procedure will be written to the specified summary file";
        " (which will be given extension '.KPopDistill.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let classes_label = TA.get_parameter () in
        Distill_kmers (classes_label, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--distance"; "--distance-function" ],
      Some "'euclidean'|'minkowski(<non-negative_float>)'",
      [ "set the function to be used when computing distances.";
        "The parameter for 'minkowski()' is the power" ],
      TA.Default (Space.Distance.to_string Defaults.distance |> Fun.const),
      (fun _ -> Distance_set (TA.get_parameter () |> Space.Distance.of_string) |> List.accum Parameters.program);
    [ "--distance-normalize"; "--distance-normalization" ],
      Some "'true'|'false'",
      [ "whether spectra should be normalized prior to computing distances" ],
      TA.Default (string_of_bool Defaults.distance_normalise |> Fun.const),
      (fun _ -> Distance_normalisation_set (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "--distances"; "--compute-distances"; "--compute-spectral-distances" ],
      Some "REGEXP_SELECTOR REGEXP_SELECTOR <binary_file_prefix>",
      [ "where REGEXP_SELECTOR :=";
        " <metadata_field>'~'<regexp>[','...','<metadata_field>'~'<regexp>]";
        "and regexps are defined according to <https://ocaml.org/api/Str.html>:";
        "select two sets of spectra from the register";
        "and compute and output distances between all possible pairs";
        " (metadata fields must match the regexps specified in the selector;";
        "  an empty metadata field makes the regexp match labels.";
        "  The result will be given extension '.KPopDMatrix' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let regexps_1 = TA.get_parameter () |> parse_regexp_selector "-d" in
        let regexps_2 = TA.get_parameter () |> parse_regexp_selector "-d" in
        To_distances (regexps_1, regexps_2, TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Actions involving the selection register:" ];
    [ "-L"; "--labels"; "--selection-from-labels" ],
      Some "<spectrum_label>[','...','<spectrum_label>]",
      [ "put into the selection register the specified labels" ],
      TA.Optional,
      (fun _ ->
        let labels = TA.get_parameter () in
        if labels <> "" then
        Selected_from_labels (labels |> String.Split.on_char_as_list ',' |> StringSet.of_list)
          |> List.accum Parameters.program);
    [ "-R"; "--regexps"; "--selection-from-regexps" ],
      Some "<metadata_field>'~'<regexp>[','...','<metadata_field>'~'<regexp>]",
      [ "put into the selection register the labels of the spectra";
        "whose metadata fields match the specified regexps";
        "and where regexps are defined according to <https://ocaml.org/api/Str.html>.";
        "An empty metadata field makes the regexp match labels" ],
      TA.Optional,
      (fun _ ->
        Selected_from_regexps (TA.get_parameter () |> parse_regexp_selector "-R")
          |> List.accum Parameters.program);
    [ "-A"; "--add-combined-selection"; "--selection-combine-and-add" ],
      Some "<spectrum_label>",
      [ "combine the spectra whose labels are in the selection register ";
        "and add the result (or replace it if a spectrum named <spectrum_label>";
        "already exists) to the database present in the database register" ],
      TA.Optional,
      (fun _ -> Add_combined_selected (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-D"; "--delete"; "--selection-delete" ],
      None,
      [ "drop the spectra whose labels are in the selection register";
        "from the database present in the register" ],
      TA.Optional,
      (fun _ -> Remove_selected |> List.accum Parameters.program);
    [ "-N"; "--selection-negate" ],
      None,
      [ "negate the labels that are present in the selection register" ],
      TA.Optional,
      (fun _ -> Selected_negate |> List.accum Parameters.program);
    [ "-P"; "--selection-print" ],
      None,
      [ "print the labels that are present in the selection register" ],
      TA.Optional,
      (fun _ -> Selected_print |> List.accum Parameters.program);
    [ "-C"; "--selection-clear" ],
      None,
      [ "purge the selection register" ],
      TA.Optional,
      (fun _ -> Selected_clear |> List.accum Parameters.program);
    TA.make_separator_multiline [ "Miscellaneous options."; "They are set immediately" ];
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (string_of_int Defaults.threads |> Fun.const),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (Fun.const "quiet execution"),
      (fun _ -> Parameters.verbose := true);
    [ "-V"; "--version" ],
      None,
      [ "print version and exit" ],
      TA.Optional,
      (fun _ -> Printf.printf "%s\n%!" info.version; exit 0);
    (* Hidden option to output help in markdown format *)
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
  (* These are the registers available to the program *)
  let current = KMerDB.empty () |> ref and selected = ref StringSet.empty
  and merge_criterion = ref Defaults.merge_criterion and combination_criterion = ref Defaults.combination_criterion
  and output_zero_kmers = ref Defaults.output_zero_kmers and precision = ref Defaults.precision
  and distance = ref Defaults.distance and distance_normalise = ref Defaults.distance_normalise in
  try
    List.iter
      (function
        | Empty ->
          current := KMerDB.empty ()
        | Of_binary prefix ->
          current := KMerDB.of_binary ~verbose:!Parameters.verbose prefix
        | Of_tables prefix ->
          current := KMerDB.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose prefix
        | Set_merge_criterion criterion ->
          merge_criterion := criterion
        | Merge_binary prefix ->
          current :=
            KMerDB.of_binary ~verbose:!Parameters.verbose prefix |>
            KMerDB.merge ~verbose:!Parameters.verbose !merge_criterion !current
        | Add_metadata path ->
          current := KMerDB.add_metadata_file ~threads:!Parameters.threads ~verbose:!Parameters.verbose !current path
        | Set_combination_criterion criterion ->
          combination_criterion := criterion
        | Split_spectra classes_label ->
          current :=
            KMerDB.split_spectra ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !current classes_label !combination_criterion
        | Add_combined_selected new_label ->
          current :=
            KMerDB.add_combined_selected ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !current new_label !selected !combination_criterion
        | Remove_selected ->
          current := KMerDB.remove_selected !current !selected
        | Transform transformation ->
          current :=
            KMerDB.transform ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            transformation !current
        | Distill_kmers (classes_label, summary_prefix) ->
          KMerDB.distill_kmers ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !current classes_label summary_prefix
        | Summary ->
          KMerDB.output_summary ~verbose:!Parameters.verbose !current
        | Selected_from_labels labels ->
          selected := labels
        | Selected_from_regexps regexps ->
          selected := KMerDB.selected_from_regexps ~verbose:!Parameters.verbose !current regexps
        | Selected_negate ->
          selected := KMerDB.selected_negate !current !selected
        | Selected_print ->
          Printf.eprintf "Currently selected spectra = [";
          StringSet.iter (Printf.eprintf " '%s'%!") !selected;
          Printf.eprintf " ].\n%!"
        | Selected_clear ->
          selected := StringSet.empty
        | Set_output_zero_kmers ozr ->
          output_zero_kmers := ozr
        | Set_precision p ->
          precision := p
        | To_tables prefix ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              KMerDB.to_files ~precision:!precision ~output_zero_kmers:!output_zero_kmers
                              ~threads:!Parameters.threads ~verbose:!Parameters.verbose
                              !current prefix)
        | To_binary prefix ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () -> KMerDB.to_binary ~verbose:!Parameters.verbose !current prefix)
        | Distance_set dist ->
          distance := dist
        | Distance_normalisation_set normalise ->
          distance_normalise := normalise
        | To_distances (regexps_1, regexps_2, prefix) ->
          let selected_1 = KMerDB.selected_from_regexps ~verbose:!Parameters.verbose !current regexps_1
          and selected_2 = KMerDB.selected_from_regexps ~verbose:!Parameters.verbose !current regexps_2 in
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              KMerDB.to_distances
                ~precision:!precision ~normalise:!distance_normalise
                ~threads:!Parameters.threads ~verbose:!Parameters.verbose
                !distance !current selected_1 selected_2 prefix))
      program
  with e ->
    Exception.handle __FUNCTION__ TA.usage (fun () ->
      Printf.peprintf "(%s): This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!" __FUNCTION__;
      Printf.peprintf "(%s): You might also wish to rerun me with option -x to get a full backtrace.\n%!" __FUNCTION__
    ) e

