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

module KMerCounter (KMH: KMer.KMerHash with type t = int):
  sig
    val compute: ?verbose:bool -> linter:(string -> string) -> KMer.ReadFiles.t -> int -> string -> string -> unit
  end
= struct
    let compute ?(verbose = false) ~linter store max_results_size label fname =
      let output =
        if fname = "" then
          stdout
        else
          open_out fname in
      (* Header with label *)
      Printf.fprintf output "\t%s\n" label;
      let output_format = Scanf.format_from_string (Printf.sprintf "%%0%dx\t%%d\n" ((KMH.k + 1) / 2)) "%d%d" in
      let reads_cntr = ref 0 and res = KMer.IntHashtbl.create max_results_size in
      KMer.ReadFiles.iter ~linter ~verbose:false
        (fun _ segm_id read ->
          KMH.iter
            (fun hash occs ->
              KMer.add_to_kmer_counter res hash occs;
              if KMer.IntHashtbl.length res > max_results_size then begin
                let min_binding =
                  KMer.IntHashtbl.fold (fun _ occs old_min -> min old_min !occs) res max_int in
                if verbose then
                  Printf.eprintf "\rKMerCounter.compute: Outputting and removing hashes having #%d occurrences...%!"
                    min_binding;
                KMer.IntHashtbl.filter_map_inplace
                  (fun hash occs ->
                    if !occs = min_binding then begin
                      Printf.fprintf output output_format hash !occs;
                      (* Discard the binding *)
                      None
                    end else
                      (* Keep the binding *)
                      Some occs)
                  res;
                if verbose then
                  Printf.eprintf " done.\n%!"
              end)
            read.seq;
          if verbose && !reads_cntr mod 10000 = 0 then
            Printf.eprintf "\rKMerCounter.compute: Added %d reads%!" !reads_cntr;
          if segm_id = 0 then
            incr reads_cntr)
        store;
      if verbose then begin
        Printf.eprintf "\rKMerCounter.compute: Added %d reads\n%!" !reads_cntr;
        Printf.eprintf "KMerCounter.compute: Outputting hashes...%!";
      end;
      KMer.IntHashtbl.iter (fun hash occs -> Printf.fprintf output output_format hash !occs) res;
      if verbose then
        Printf.eprintf " done.%!\n";
      close_out output
  end

module Content =
  struct
    type t =
      | DNA
      | Protein
    exception Invalid_content of string
    let of_string = function
      | "DNA" | "dna" -> DNA
      | "protein" | "prot" -> Protein
      | w -> Invalid_content w |> raise
    let to_string = function
      | DNA -> "DNA"
      | Protein -> "protein"
  end

module Defaults =
  struct
    let content = Content.DNA
    let k = 12
    let max_results_size = 16777216 (* Or: 4^12 *)
    let output = ""
(*
    let threads = Tools.Parallel.get_nproc ()
*)
    let verbose = false
  end

module Parameters =
  struct
    let content = ref Defaults.content
    let k = ref Defaults.k
    let max_results_size = ref Defaults.max_results_size
    let inputs = ref []
    let label = ref ""
    let output = ref Defaults.output
    (*let threads = ref Defaults.threads*)
    let verbose = ref Defaults.verbose
  end

let version = "0.5"

let header =
  Printf.sprintf begin
    "This is the KPopCount program (version %s)\n%!" ^^
    " (c) 2017-2022 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!"
  end version

let _ =
  let module TA = Tools.Argv in
  TA.set_header header;
  TA.set_synopsis "-l|--label <output_vector_label> [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "-k"; "-K"; "--k-mer-size"; "--k-mer-length" ],
      Some "<k_mer_length>",
      [ "k-mer length"; "(must be positive, and <= 30 for DNA or <= 12 for protein)" ],
      TA.Default (fun () -> string_of_int !Parameters.k),
      (fun _ -> Parameters.k := TA.get_parameter_int_pos ());
    [ "-m"; "-M"; "--max-results-size" ],
      Some "<positive_integer>",
      [ "maximum number of k-mer signatures to be kept in memory at any given time.";
        "If more are present, the ones corresponding to the lowest cardinality";
        "will be removed from memory and printed out, and there will be";
        "repeated signatures in the output" ],
      TA.Default (fun () -> string_of_int !Parameters.max_results_size),
      (fun _ -> Parameters.max_results_size := TA.get_parameter_int_pos ());
    TA.make_separator "Input/Output";
    [ "-c"; "-C"; "--content"; "--mode" ],
      Some "'DNA'|'protein'",
      [ "how file contents should be interpreted" ],
      TA.Default (fun () -> Content.to_string !Parameters.content),
      (fun _ -> Parameters.content := TA.get_parameter () |> Content.of_string);
    [ "-f"; "-F"; "--fasta" ],
      Some "<fasta_file_name>",
      [ "FASTA input file containing sequences" ],
      TA.Optional,
      (fun _ -> KMer.ReadFiles.FASTA (TA.get_parameter ()) |> Tools.List.accum Parameters.inputs);
    [ "-s"; "-S"; "--single-end" ],
      Some "<fastq_file_name>",
      [ "FASTQ input file containing single-end sequencing reads" ],
      TA.Optional,
      (fun _ -> KMer.ReadFiles.SingleEndFASTQ (TA.get_parameter ()) |> Tools.List.accum Parameters.inputs);
    [ "-p"; "-P"; "--paired-end" ],
      Some "<fastq_file_name1> <fastq_file_name2>",
      [ "FASTQ input files containing paired-end sequencing reads" ],
      TA.Optional,
      (fun _ ->
        let name1 = TA.get_parameter () in
        let name2 = TA.get_parameter () in
        KMer.ReadFiles.PairedEndFASTQ (name1, name2) |> Tools.List.accum Parameters.inputs);
    [ "-l"; "--label" ],
      Some "<output_vector_label>",
      [ "label of the k-mer vector in the output file" ],
      TA.Mandatory,
      (fun _ -> Parameters.label := TA.get_parameter() );
    [ "-o"; "--output" ],
      Some "<output_file_name>",
      [ "name of generated output file" ],
      TA.Default (fun () -> if !Parameters.output = "" then "<stdout>" else !Parameters.output),
      (fun _ -> Parameters.output := TA.get_parameter() );
    TA.make_separator "Miscellaneous";
(*
    [ "-t"; "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
*)
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (fun () -> string_of_bool !Parameters.verbose),
      (fun _ -> Parameters.verbose := true);
    (* Hidden option to emit help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 1)
  ];
  let module KMCD = KMerCounter (KMer.DNAEncodingHash (struct let value = !Parameters.k end)) in
  let module KMCP = KMerCounter (KMer.ProteinEncodingHash (struct let value = !Parameters.k end)) in
  Parameters.inputs := List.rev !Parameters.inputs;
  if !Parameters.inputs <> [] then begin
    let store = ref KMer.ReadFiles.empty in
    List.iter
      (fun input -> store := KMer.ReadFiles.add_from_files !store input)
      !Parameters.inputs;
    begin match !Parameters.content with
    | DNA -> KMCD.compute ~linter:(Sequences.Lint.dnaize ~keep_dashes:false)
    | Protein -> KMCP.compute ~linter:(Sequences.Lint.proteinize ~keep_dashes:false)
    end ~verbose:!Parameters.verbose !store !Parameters.max_results_size !Parameters.label !Parameters.output
  end

