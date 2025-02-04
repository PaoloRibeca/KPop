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

module Defaults =
  struct
    let algorithm = Matrix.SplitsAlgorithm.of_string "gaps"
    let splits = 10000
    let threads = Processes.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let algorithm = ref Defaults.algorithm
    let splits = ref Defaults.splits
    let input = ref ""
    let output = ref ""
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let info = {
  Tools.Argv.name = "KPopPhylo";
  version = "7";
  date = "04-Feb-2025"
} and authors = [
  "2024-2025", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "-i|--input <binary_input_prefix> -o|--output <newick_output_file> [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "--algorithm"; "--splits-algorithm" ],
      Some "'gaps'",
      [ "algorithm to use when computing splits from embeddings" ],
      TA.Default (fun () -> Matrix.SplitsAlgorithm.to_string Defaults.algorithm),
      (fun _ -> Parameters.algorithm := TA.get_parameter () |> Matrix.SplitsAlgorithm.of_string);
    [ "-m"; "--max-splits" ],
      Some "<positive_integer>",
      [ "keep at most that many splits when generating the tree" ],
      TA.Default (fun () -> string_of_int Defaults.splits),
      (fun _ -> Parameters.splits := TA.get_parameter_int_pos ());
    TA.make_separator "Input/Output";
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load the specified embedding database and generate a tree from it.";
        "File extension is automatically determined";
        " (will be '.KPopVectors' unless file is '/dev/*')" ],
      TA.Mandatory,
      (fun _ -> Parameters.input := TA.get_parameter ());
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "use this prefix when outputting the generated tree and filtered splits.";
        "File extensions are automatically determined";
        " (will be '.nwk' for the tree; '.KPopPSplits.txt' for the splits)" ],
      TA.Mandatory,
      (fun _ -> Parameters.output := TA.get_parameter ());
    TA.make_separator "Miscellaneous";
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (fun () -> string_of_int Defaults.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (fun () -> string_of_bool Defaults.verbose),
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
  begin try
    let input =
      let res = Matrix.of_binary ?verbose:(Some !Parameters.verbose) Vectors !Parameters.input in
      if res.Matrix.which <> Vectors then
        Matrix.Unexpected_type (res.which, Vectors) |> raise;
      res in
    let splits =
      Matrix.get_splits ~threads:!Parameters.threads ~verbose:!Parameters.verbose
      !Parameters.algorithm !Parameters.splits input in
    let output_tree = !Parameters.output ^ ".nwk"
    and output_splits_ok = !Parameters.output ^ ".ok.KPopPSplits.txt"
    and output_splits_ko = !Parameters.output ^ ".ko.KPopPSplits.txt" in
    let splits_ok, tree, splits_ko = Trees.Splits.to_tree ~verbose:!Parameters.verbose splits in
    Trees.Newick.to_file tree output_tree;
    Trees.Splits.to_file splits_ok output_splits_ok;
    Trees.Splits.to_file splits_ko output_splits_ko
  with e ->
    TA.header ();
    raise e
  end;

