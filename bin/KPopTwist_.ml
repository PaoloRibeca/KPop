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
    (* The following three are KPopCountDB.TableFilter.default *)
    let transformation = ref "normalize"
    let threshold = ref 1
    let power = ref 1.
    let sampling = ref 1.
    (*let precision = 15*)
    let threads = Processes.Parallel.get_nproc () |> ref
    let temporaries = ref false
    let verbose = ref false
  end

let info = {
  Tools.Argv.name = "KPopTwist";
  version = "17";
  date = "02-Jan-2024"
} and authors = [
  "2022-2024", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "-i|--input <input_table_prefix> [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "-f"; "-F"; "-s"; "-S"; "--fraction"; "--sampling"; "--sampling-fraction" ],
      Some "<fractional_float>",
      [ "fraction of the rows to be considered and resampled before twisting" ],
      TA.Default (fun () -> string_of_float !Parameters.sampling),
      (fun _ -> Parameters.sampling := TA.get_parameter_float_fraction ());
    [ "--threshold" ],
      Some "<non_negative_integer>",
      [ "set to zero all counts that are less than this threshold";
        "before transforming them" ],
      TA.Default (fun () -> string_of_int !Parameters.threshold),
      (fun _ -> Parameters.threshold := TA.get_parameter_int_non_neg ());
    [ "--power" ],
      Some "<non_negative_float>",
      [ "raise counts to this power before transforming them.";
        "A power of 0 when the 'pseudocount' method is used";
        "performs a logarithmic transformation" ],
      TA.Default (fun () -> string_of_float !Parameters.power),
      (fun _ -> Parameters.power := TA.get_parameter_float_non_neg ());
    [ "--transform"; "--transformation" ],
      Some "'none'|'normalize'|'pseudocount'|'clr'",
      [ "transformation to apply to table elements" ],
      TA.Default (fun () -> !Parameters.transformation),
      (fun _ -> Parameters.transformation := TA.get_parameter ());
    TA.make_separator "Input/Output";
    [ "-i"; "--input" ],
      Some "<input_table_prefix>",
      [ "load the specified k-mer database in the register and twist it.";
        "File extension is automatically determined";
        " (will be .KPopCounter).";
        "The prefix is then re-used for output";
        " (and the output files will be given extensions";
        "  .KPopTwister and .KPopTwisted)" ],
      TA.Mandatory,
      (fun _ -> Parameters.input := TA.get_parameter ());
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
  Printf.printf "%s\001%s\001%d\001%.12g\001%.12g\001%d\001%b\001%b\n%!"
    !Parameters.input !Parameters.transformation !Parameters.threshold !Parameters.power
    !Parameters.sampling !Parameters.threads !Parameters.temporaries !Parameters.verbose

