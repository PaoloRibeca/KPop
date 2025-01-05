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

module Parameters =
  struct
    let input = ref ""
    let output = ref ""
    let output_kmers = ref ""
    (* The following three are KPopCountDB.TableFilter.default *)
    let sampling = ref 1.
    let threshold_counts = ref 1.
    let power = ref 1.
    (*let precision = 15*)
    let transformation = ref "power"
    let normalize = ref true
    let threshold_kmers = ref 0.
    let threads = Processes.Parallel.get_nproc () |> ref
    let temporaries = ref false
    let verbose = ref false
  end

let info = {
  Tools.Argv.name = "KPopTwist";
  version = "21";
  date = "29-Oct-2024"
} and authors = [
  "2022-2024", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "-i|--input <binary_input_prefix> -o|--output <binary_output_prefix> [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "-f"; "-F"; "-s"; "-S"; "--fraction"; "--sampling"; "--sampling-fraction" ],
      Some "<fractional_float>",
      [ "fraction of the k-mers to be considered and resampled before twisting" ],
      TA.Default (fun () -> string_of_float !Parameters.sampling),
      (fun _ -> Parameters.sampling := TA.get_parameter_float_fraction ());
    [ "--threshold-counts" ],
      Some "<non_negative_integer>",
      [ "set to zero all counts that are less than this threshold";
        "before transforming them.";
        "A fractional threshold between 0. and 1. is taken as a relative one";
        "with respect to the sum of all counts in the spectrum" ],
      TA.Default (fun () -> string_of_float !Parameters.threshold_counts),
      (fun _ -> Parameters.threshold_counts := TA.get_parameter_float_non_neg ());
    [ "--power" ],
      Some "<non_negative_float>",
      [ "raise counts to this power before transforming them.";
        "A power of 0 when the 'pseudocounts' method is used";
        "performs a logarithmic transformation" ],
      TA.Default (fun () -> string_of_float !Parameters.power),
      (fun _ -> Parameters.power := TA.get_parameter_float_non_neg ());
    [ "--transform"; "--transformation" ],
      Some "'binary'|'power'|'pseudocounts'|'clr'",
      [ "transformation to apply to table elements" ],
      TA.Default (fun () -> !Parameters.transformation),
      (fun _ -> Parameters.transformation := TA.get_parameter ());
    [ "-n"; "--normalize"; "--normalize-counts" ],
      Some "'true'|'false'",
      [ "whether to normalize spectra after transformation and before twisting" ],
      TA.Default (fun () -> string_of_bool !Parameters.normalize),
      (fun _ -> Parameters.normalize := TA.get_parameter_boolean ());
    [ "--threshold-kmers" ],
      Some "<non_negative_integer>",
      [ "compute the sum of all transformed (and possibly normalized) counts";
        "for each k-mer, and eliminate k-mers such that the corresponding sum";
        "is less than the largest sum rescaled by this threshold.";
        "This filters out k-mers having low frequencies across all spectra" ],
      TA.Default (fun () -> string_of_float !Parameters.threshold_kmers),
      (fun _ -> Parameters.threshold_kmers := TA.get_parameter_float_non_neg ());
    TA.make_separator "Input/Output";
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load the specified k-mer database in the register and twist it.";
        "File extension is automatically determined";
        " (will be .KPopCounter).";
        "The prefix is then re-used for output";
        " (and the output files will be given extensions";
        "  .KPopTwister and .KPopTwisted)" ],
      TA.Mandatory,
      (fun _ -> Parameters.input := TA.get_parameter ());
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "use this prefix when dumping generated twister and twisted sequences.";
        "File extensions are automatically determined";
        " (will be .KPopTwister and .KPopTwisted)" ],
      TA.Mandatory,
      (fun _ -> Parameters.output := TA.get_parameter ());
    [ "-k"; "--output-kmers" ],
      Some "<binary_file_prefix>",
      [ "use this prefix when dumping generated twisted k-mers.";
        "File extension is automatically determined";
        " (will be .KPopTwisted)" ],
      TA.Default (fun () -> !Parameters.output_kmers),
      (fun _ -> Parameters.output_kmers := TA.get_parameter ());
    TA.make_separator "Miscellaneous";
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
    [ "--keep-temporaries" ],
      None,
      [ "keep temporary files rather than deleting them in the end" ],
      TA.Default (fun () -> string_of_bool !Parameters.temporaries),
      (fun _ -> Parameters.temporaries := true);
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
  (*let program = List.rev !Parameters.program in
  if program = [] then begin
    TA.usage ();
    exit 0
  end;*)
  if !Parameters.verbose then
    TA.header ();

  (* For the moment, we just check if the file is there, and repeat input parameters *)

  begin try
    !Parameters.input ^ ".KPopCounter" |> open_in |> close_in
  with e ->
    TA.header ();
    raise e
  end;
  Printf.printf "%s\001%.12g\001%.12g\001%.12g\001%s\001%b\001%.12g\001%s\001%s\001%d\001%b\001%b\n%!"
    !Parameters.input !Parameters.sampling !Parameters.threshold_counts !Parameters.power !Parameters.transformation
    !Parameters.normalize !Parameters.threshold_kmers !Parameters.output !Parameters.output_kmers !Parameters.threads
    !Parameters.temporaries !Parameters.verbose

