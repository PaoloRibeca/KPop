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

type to_do_t =
  | Empty
  | Of_file of string
  | Add_meta of string
  | Add_files of string list
  | Combination_criterion_set of KMerDB.CombinationCriterion.t
  | Add_combined_selected of string (* The new label *)
  | Remove_selected
  | Summary
  | Selected_from_labels of StringSet.t
  | Selected_from_regexps of regexps_t
  | Selected_negate
  | Selected_print
  | Selected_clear
  | Selected_to_filter
  | To_file of string
  | Table_output_row_names of bool
  | Table_output_col_names of bool
  | Table_output_metadata of bool
  | Table_transpose of bool
  | Table_transform_threshold of float
  | Table_transform_power of float
  | Table_transform_which of string
  | Table_output_zero_rows of bool
  | Table_precision of int
  | To_table of string
  | To_spectra of string
  | Distance_set of Space.Distance.t
  | Distance_normalisation_set of bool
  | To_distances of regexps_t * regexps_t * string
and regexps_t = (string * Str.regexp) list

module Defaults =
  struct
    let combination_criterion = KMerDB.CombinationCriterion.of_string "mean"
    let filter = KMerDB.TableFilter.default
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalise = true
  end

module Parameters =
  struct
    let program = ref []
    let threads = Processes.Parallel.get_nproc () |> ref
    let verbose = ref false
  end

let info = {
  Tools.Argv.name = "KPopCountDB";
  version = "47";
  date = "05-Feb-2025"
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
    TA.make_separator_multiline [ ""; "Actions on the database register:" ];
    [ "-e"; "--empty" ],
      None,
      [ "put an empty database into the register" ],
      TA.Optional,
      (fun _ -> Empty |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load into the register the database present in the specified file";
        " (which must have extension '.KPopCounter' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Of_file (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-m"; "--metadata"; "--add-metadata" ],
      Some "<metadata_table_file_name>",
      [ "add to the database present in the register metadata from the specified file.";
        "Metadata field names and values must not contain double quote '\"' characters" ],
      TA.Optional,
      (fun _ -> Add_meta (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-k"; "--kmers"; "--add-kmers"; "--add-kmer-files" ],
      Some "<k-mer_table_file_name>[','...','<k-mer_table_file_name>]",
      [ "add to the database present in the register k-mers from the specified files";
        " (which must have extension '.KPopSpectra.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        Add_files (TA.get_parameter () |> String.Split.on_char_as_list ',') |> List.accum Parameters.program);
    [ "--summary" ],
      None,
      [ "print a summary of the database present in the register" ],
      TA.Optional,
      (fun _ -> Summary |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "dump the database present in the register to the specified file";
        " (which will be given extension '.KPopCounter' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> To_file (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--distance"; "--distance-function" ],
      Some "'euclidean'|'minkowski(<non_negative_float>)'",
      [ "set the function to be used when computing distances.";
        "The parameter for 'minkowski()' is the power" ],
      TA.Default (fun () -> Space.Distance.to_string Defaults.distance),
      (fun _ -> Distance_set (TA.get_parameter () |> Space.Distance.of_string) |> List.accum Parameters.program);
    [ "--distance-normalize"; "--distance-normalization" ],
      Some "'true'|'false'",
      [ "whether spectra should be normalized prior to computing distances" ],
      TA.Optional,
      (fun _ -> Distance_normalisation_set (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "-d"; "--distances"; "--compute-distances"; "--compute-spectral-distances" ],
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
    [ "--table-output-row-names" ],
      Some "'true'|'false'",
      [ "whether to output row names for the database present in the register";
        "when writing it as a tab-separated file" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_row_names),
      (fun _ -> Table_output_row_names (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "--table-output-col-names" ],
      Some "'true'|'false'",
      [ "whether to output column names for the database present in the register";
        "when writing it as a tab-separated file" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_col_names),
      (fun _ -> Table_output_col_names (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "--table-output-metadata" ],
      Some "'true'|'false'",
      [ "whether to output metadata for the database present in the register";
        "when writing it as a tab-separated file" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_metadata),
      (fun _ -> Table_output_metadata (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "--table-transpose" ],
      Some "'true'|'false'",
      [ "whether to transpose the database present in the register";
        "before writing it as a tab-separated file";
        " (if 'true': rows are spectrum names, columns [metadata and] k-mer names;";
        "  if 'false': rows are [metadata and] k-mer names, columns spectrum names)" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.transpose),
      (fun _ -> Table_transpose (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "--counts-threshold" ],
      Some "<non_negative_integer>",
      [ "set to zero all counts that are less than this threshold";
        "before transforming and outputting them as table or spectra.";
        "A fractional threshold between 0. and 1. is taken as a relative one";
        "with respect to the sum of all counts in the spectrum" ],
      TA.Default
        (fun () -> (KMerDB.Transformation.to_parameters Defaults.filter.transform).threshold |> string_of_float),
      (fun _ ->
        Table_transform_threshold (TA.get_parameter_float_non_neg ()) |> List.accum Parameters.program);
    [ "--counts-power" ],
      Some "<non_negative_float>",
      [ "raise counts to this power before transforming and outputting them";
        "as table or spectra.";
        "A power of 0 when the 'pseudocounts' method is used";
        "performs a logarithmic transformation" ],
      TA.Default (fun () -> (KMerDB.Transformation.to_parameters Defaults.filter.transform).power |> string_of_float),
      (fun _ ->
        Table_transform_power (TA.get_parameter_float_non_neg ()) |> List.accum Parameters.program);
    [ "--counts-transform"; "--counts-transformation" ],
      Some "'binary'|'power'|'pseudocounts'|'clr'",
      [ "transformation to apply to counts before outputting them as table or spectra" ],
      TA.Default (fun () -> (KMerDB.Transformation.to_parameters Defaults.filter.transform).which),
      (fun _ ->
        Table_transform_which (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--counts-output-zero-kmers"; "--counts-output-zero-k-mers" ],
      Some "'true'|'false'",
      [ "whether to output k-mers whose frequencies are all zero";
        "when writing the database as table or spectra" ],
      TA.Default (fun () -> string_of_bool Defaults.filter.print_zero_rows),
      (fun _ -> Table_output_zero_rows (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "--counts-precision" ],
      Some "<positive_integer>",
      [ "set the number of precision digits to be used when outputting counts";
        "as table or spectra" ],
      TA.Default (fun () -> string_of_int Defaults.filter.precision),
      (fun _ -> Table_precision (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-t"; "--table" ],
      Some "<file_prefix>",
      [ "write the database present in the register as a tab-separated file";
        " (rows are k-mer names, columns are spectrum names.";
        "  The result will be given extension '.KPopCounter.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> To_table (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-s"; "--spectra" ],
      Some "<file_prefix>",
      [ "write the database present in the register as k-mer spectra";
        " (the result will be given extension '.KPopSpectra.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> To_spectra (TA.get_parameter ()) |> List.accum Parameters.program);
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
    [ "--selection-combination-criterion"; "--combination-criterion" ],
      Some "'mean'|'median'",
      [ "set the criterion used to combine the k-mer frequencies of selected spectra.";
        "To avoid rounding issues, each k-mer frequency is also rescaled";
        "by the largest normalization across spectra";
        " ('mean' averages frequencies across spectra;";
        "  'median' computes the median across spectra)" ],
      TA.Default (fun () -> KMerDB.CombinationCriterion.to_string Defaults.combination_criterion),
      (fun _ ->
        Combination_criterion_set (TA.get_parameter () |> KMerDB.CombinationCriterion.of_string)
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
    [ "-F"; "--selection-to-table-filter" ],
      None,
      [ "filter out spectra whose labels are present in the selection register";
        "when writing the database as a tab-separated file" ],
      TA.Optional,
      (fun _ -> Selected_to_filter |> List.accum Parameters.program);
    TA.make_separator_multiline [ "Miscellaneous options."; "They are set immediately" ];
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
    (* Hidden option to output help in markdown format *)
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
  (* These are the registers available to the program *)
  let current = KMerDB.make_empty () |> ref and selected = ref StringSet.empty
  and combination_criterion = ref Defaults.combination_criterion
  and transform = KMerDB.Transformation.to_parameters Defaults.filter.transform |> ref
  and filter = ref Defaults.filter
  and distance = ref Defaults.distance and distance_normalise = ref Defaults.distance_normalise in
  try
    List.iter
      (function
        | Empty ->
          current := KMerDB.make_empty ()
        | Of_file prefix ->
          current := KMerDB.of_binary ~verbose:!Parameters.verbose prefix
        | Add_meta fname ->
          current := KMerDB.add_meta ~verbose:!Parameters.verbose !current fname
        | Add_files prefixes ->
          current := KMerDB.add_files ~verbose:!Parameters.verbose !current prefixes
        | Combination_criterion_set criterion ->
          combination_criterion := criterion
        | Add_combined_selected new_label ->
          current :=
            KMerDB.add_combined_selected ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !current new_label !selected !combination_criterion
        | Remove_selected ->
          current := KMerDB.remove_selected !current !selected
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
        | Selected_to_filter ->
          filter := { !filter with filter_columns = !selected }
        | Table_output_row_names print_row_names ->
          filter := { !filter with print_row_names }
        | Table_output_col_names print_col_names ->
          filter := { !filter with print_col_names }
        | Table_output_metadata print_metadata ->
          filter := { !filter with print_metadata }
        | Table_transpose transpose ->
          filter := { !filter with transpose }
        | Table_transform_threshold threshold ->
          transform := { !transform with threshold };
          filter := { !filter with transform = KMerDB.Transformation.of_parameters !transform }
        | Table_transform_power power ->
          transform := { !transform with power };
          filter := { !filter with transform = KMerDB.Transformation.of_parameters !transform }
        | Table_transform_which which ->
          transform := { !transform with which };
          filter := { !filter with transform = KMerDB.Transformation.of_parameters !transform }
        | Table_output_zero_rows print_zero_rows ->
          filter := { !filter with print_zero_rows }
        | Table_precision precision ->
          filter := { !filter with precision }
        | To_table prefix ->
          KMerDB.to_table ~filter:!filter ~threads:!Parameters.threads ~verbose:!Parameters.verbose !current prefix
        | To_spectra prefix ->
          KMerDB.to_spectra ~filter:!filter ~threads:!Parameters.threads ~verbose:!Parameters.verbose !current prefix
        | To_file prefix ->
          KMerDB.to_binary ~verbose:!Parameters.verbose !current prefix
        | Distance_set dist ->
          distance := dist
        | Distance_normalisation_set normalise ->
          distance_normalise := normalise
        | To_distances (regexps_1, regexps_2, prefix) ->
          let selected_1 = KMerDB.selected_from_regexps ~verbose:!Parameters.verbose !current regexps_1
          and selected_2 = KMerDB.selected_from_regexps ~verbose:!Parameters.verbose !current regexps_2 in
          KMerDB.to_distances ~normalise:!distance_normalise ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !distance !current selected_1 selected_2 prefix)
      program
  with exc ->
    Printf.peprintf "(%s): %s\n%!" __FUNCTION__
      ("FATAL: Uncaught exception: " ^ Printexc.to_string exc |> String.TermIO.red)

