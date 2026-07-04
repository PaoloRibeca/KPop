(*
    KPopTwist.ml -- (c) 2022-2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    KPopTwist is the binary that runs correspondence analysis on a
    binary `.KPopSpectra` database and emits the `.KPopTwister`
    (CA transformation) and `.KPopTwisted` (sample standard
    coordinates) outputs.  Selects between the full LAPACK SVD
    (`CA.twist`) and the randomised truncated SVD (`CA.rsvd`)
    according to the `--dimensions` flag, and exposes the k-mer
    filtering control parameters (keep-list, random subsampling,
    row-sum threshold, condition-number filter).

    This program was designed and developed by the author(s),
    with the assistance of the following AI tool(s):
      2026 Claude (Anthropic).
    The final logic and implementation were reviewed and verified in
    their entirety by the author(s).

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

module Defaults =
  struct
    let kmers_keep = ""
    let kmers_sample = 1.
    let threshold_kmers = CA.Filter.Auto
    let condition_number = CA.Filter.Off
    let dimensions = 0 (* 0 means: use full SVD via twist *)
    let threads = Processes.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let input = ref ""
    let output = ref ""
    let output_kmers = ref ""
    let kmers_keep = ref Defaults.kmers_keep
    let kmers_sample = ref Defaults.kmers_sample
    let threshold_kmers = ref Defaults.threshold_kmers
    let condition_number = ref Defaults.condition_number
    let dimensions = ref Defaults.dimensions
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let info = {
  Tools.Argv.name = "KPopTwist";
  version = "33";
  date = "07-Apr-2026"
} and authors = [
  "2022-2026", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "-i|--input <binary_input_prefix> -o|--output <binary_output_prefix> [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "-d"; "--dimensions" ],
      Some "<positive_integer>",
      [ "number of CA dimensions to compute using the randomised SVD";
        "(Halko, Martinsson & Tropp 2011).";
        "When not set, all min(k-mers, samples) - 1 dimensions are computed";
        "using the full LAPACK SVD, which is more accurate but slower" ],
      TA.Default (Fun.const "all dimensions (full SVD)"),
      (fun _ -> Parameters.dimensions := TA.get_parameter_int_pos ());
    [ "--keep"; "--keep-kmers"; "--kmers-keep" ],
      Some "<kmer_list_file>",
      [ "discard the k-mers not listed in this file before twisting the table.";
        "The file must contain one k-mer label per line and no header" ],
      TA.Default (Fun.const "keep all"),
      (fun _ -> Parameters.kmers_keep := TA.get_parameter ());
    [ "--sample"; "--sample-kmers"; "--kmers-sample" ],
      Some "<fractional_float>",
      [ "fraction of the k-mers to be randomly resampled and kept";
        "after parameter -k has been applied and before twisting" ],
      TA.Default (string_of_float Defaults.kmers_sample |> Fun.const),
      (fun _ -> Parameters.kmers_sample := TA.get_parameter_float_fraction ());
    [ "--kmers-threshold" ],
      Some "'off'|'auto'|<non-negative_float>",
      [ "compute the sum of all counts for each k-mer, and eliminate k-mers";
        "such that the corresponding sum is less than a cutoff.";
        "'off' (or '0') disables the filter; <float> sets the cutoff to the";
        "largest row sum rescaled by that fraction (legacy semantics);";
        "'auto' picks the cutoff at the Kneedle elbow of the sorted-ascending";
        "row-sum distribution, removing the noise tail of rare and singleton";
        "k-mers without a user-supplied magic number.";
        "This filters out k-mers having low frequencies across all spectra" ],
      TA.Default (CA.Filter.to_string Defaults.threshold_kmers |> Fun.const),
      (fun _ -> Parameters.threshold_kmers := CA.Filter.of_string (TA.get_parameter ()));
    [ "--kmers-condition-number" ],
      Some "'off'|'auto'|<positive_float>",
      [ "compute the row contribution to total inertia CTR_i = ||S[i,:]||^2 for";
        "each k-mer, and eliminate k-mers whose CTR_i is below a cutoff.";
        "'off' (or '0') disables the filter; <float> sets the cutoff to";
        "max(CTR) / parameter (legacy semantics: a larger value retains more";
        "k-mers); 'auto' picks the cutoff at the Kneedle elbow of the";
        "sorted-ascending CTR distribution, removing nearly-uniform k-mers in";
        "the noise tail.";
        "This filters out k-mers that are nearly uniform across all spectra" ],
      TA.Default (CA.Filter.to_string Defaults.condition_number |> Fun.const),
      (fun _ -> Parameters.condition_number := CA.Filter.of_string (TA.get_parameter ()));
    TA.make_separator "Input/Output";
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load the specified k-mer database in binary format and twist it.";
        "File extension is automatically determined";
        " (will be '.KPopSpectra' unless file is '/dev/*')" ],
      TA.Mandatory,
      (fun _ -> Parameters.input := TA.get_parameter ());
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "use this prefix when saving the generated twister and twisted sequences.";
        "File extensions are automatically determined";
        " (will be '.KPopTwister' and '.KPopTwisted' unless file is '/dev/*')" ],
      TA.Mandatory,
      (fun _ -> Parameters.output := TA.get_parameter ());
    [ "-k"; "--output-kmers"; "--output-twisted-kmers" ],
      Some "<binary_file_prefix>",
      [ "use this prefix when saving the generated twisted k-mers.";
        "File extension is automatically determined";
        " (will be '.KPopTwisted' unless file is '/dev/*')" ],
      TA.Default (Fun.const "do not output"),
      (fun _ -> Parameters.output_kmers := TA.get_parameter ());
    TA.make_separator "Miscellaneous";
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
    (* Hidden option to emit help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 0)
  ];
  if !Parameters.verbose then
    TA.header ();
  CPU.warn_if_slow ~verbose:!Parameters.verbose ();
  (* Seed the random number generator for k-mer subsampling *)
  Random.self_init ();
  (* Load the k-mer database *)
  let db = KMerDB.of_binary ~verbose:!Parameters.verbose !Parameters.input in
  (* Perform correspondence analysis *)
  let twister, twisted_samples, twisted_kmers =
    if !Parameters.dimensions = 0 then
      CA.twist
        ~kmers_keep:!Parameters.kmers_keep
        ~kmers_sample:!Parameters.kmers_sample
        ~threshold_kmers:!Parameters.threshold_kmers
        ~condition_number:!Parameters.condition_number
        ~threads:!Parameters.threads
        ~verbose:!Parameters.verbose
        db
    else
      CA.rsvd
        ~kmers_keep:!Parameters.kmers_keep
        ~kmers_sample:!Parameters.kmers_sample
        ~threshold_kmers:!Parameters.threshold_kmers
        ~condition_number:!Parameters.condition_number
        ~threads:!Parameters.threads
        ~verbose:!Parameters.verbose
        db
        !Parameters.dimensions in
  (* Write twister (transformation matrix + inertia) *)
  Twister.to_binary ~verbose:!Parameters.verbose twister !Parameters.output;
  (* Write twisted sequences (column standard coordinates + inertia) *)
  Twisted.to_binary ~verbose:!Parameters.verbose twisted_samples !Parameters.output;
  (* Optionally write twisted k-mers (row standard coordinates + inertia) *)
  if !Parameters.output_kmers <> "" then
    Twisted.to_binary ~verbose:!Parameters.verbose twisted_kmers !Parameters.output_kmers

