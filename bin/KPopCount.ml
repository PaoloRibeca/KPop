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

module KMerCounter (KMH: KMers.Hash_t with type t = int):
  sig
    val compute: ?verbose:bool -> linter:(string -> string) -> Files.ReadsIterate.t -> int -> string -> string -> unit
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
      let reads_cntr = ref 0 and res = KMers.HashFrequencies.Base.create max_results_size in
      Files.ReadsIterate.iter ~linter ~verbose:false
        (fun _ segm_id read ->
          KMH.iterc
            (fun hash occs ->
              KMers.HashFrequencies.add res hash occs;
              if KMers.HashFrequencies.Base.length res > max_results_size then begin
                let min_binding =
                  KMers.HashFrequencies.Base.fold (fun _ occs old_min -> min old_min !occs) res max_int in
                if verbose then
                  Printf.eprintf "%s\r(%s): Outputting and removing hashes having #%d occurrences...%!"
                    Tools.String.TermIO.clear __FUNCTION__ min_binding;
                KMers.HashFrequencies.Base.filter_map_inplace
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
          if verbose && !reads_cntr mod 1_000 = 0 then
            Printf.eprintf "%s\r(%s): Added %d reads%!" Tools.String.TermIO.clear __FUNCTION__ !reads_cntr;
          if segm_id = 0 then
            incr reads_cntr)
        store;
      if verbose then begin
        Printf.eprintf "%s\r(%s): Added %d reads.\n%!" Tools.String.TermIO.clear __FUNCTION__ !reads_cntr;
        Printf.eprintf "(%s): Outputting hashes...%!" __FUNCTION__;
      end;
      KMers.HashFrequencies.iter (Printf.fprintf output output_format) res;
      if verbose then
        Printf.eprintf " done.%!\n";
      close_out output
  end

module Content =
  struct
    type t =
      | DNA_ss
      | DNA_ds
      | Protein
    exception Invalid_content of string
    let of_string = function
      | "DNA-ss" | "DNA-single-stranded" -> DNA_ss
      | "DNA-ds" | "DNA-double-stranded" -> DNA_ds
      | "protein" | "prot" -> Protein
      | w -> Invalid_content w |> raise
  end

module Parameters =
  struct
    let content = ref None
    let k = ref 12
    let max_results_size = ref 16777216 (* Or: 4^12 *)
    let inputs = ref []
    let label = ref ""
    let output = ref ""
    (*let threads = Tools.Parallel.get_nproc () |> ref*)
    let verbose = ref false
  end

let info = {
  Tools.Info.name = "KPopCount";
  version = "10";
  date = "17-Jan-2024"
} and authors = [
  "2017-2024", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.make_header info authors [ BiOCamLib.Info.info; KPop.Info.info ] |> TA.set_header;
  TA.set_synopsis "-l|--label <output_vector_label> [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "-k"; "-K"; "--k-mer-size"; "--k-mer-length" ],
      Some "<k_mer_length>",
      [ "k-mer length"; "(must be positive, and <= 30 for DNA or <= 12 for protein)" ],
      TA.Default (fun () -> string_of_int !Parameters.k),
      (fun _ -> Parameters.k := TA.get_parameter_int_pos ());
    [ "-M"; "--max-results-size" ],
      Some "<positive_integer>",
      [ "maximum number of k-mer signatures to be kept in memory at any given time.";
        "If more are present, the ones corresponding to the lowest cardinality";
        "will be removed from memory and printed out, and there will be";
        "repeated signatures in the output" ],
      TA.Default (fun () -> string_of_int !Parameters.max_results_size),
      (fun _ -> Parameters.max_results_size := TA.get_parameter_int_pos ());
    TA.make_separator "Input/Output";
    [ "-C"; "--content" ],
      Some "'DNA-ss'|'DNA-single-stranded'|'DNA-ds'|'DNA-double-stranded'|'protein'",
      [ "how file contents should be interpreted.";
        "When content is 'DNA-ss' or 'protein', the sequence is hashed;";
        "when content is 'DNA-ds', both sequence and reverse complement are hashed" ],
      TA.Default (fun _ -> "'DNA-ss' for FASTA input, 'DNA-ds' for FASTQ input"),
      (fun _ -> Parameters.content := Some (TA.get_parameter () |> Content.of_string));
    [ "-f"; "--fasta" ],
      Some "<fasta_file_name>",
      [ "FASTA input file containing sequences.";
        "You can specify more than one FASTA input, but not FASTA and FASTQ inputs";
        "at the same time. Contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ -> Files.Type.FASTA (TA.get_parameter ()) |> Tools.List.accum Parameters.inputs);
    [ "-s"; "--single-end" ],
      Some "<fastq_file_name>",
      [ "FASTQ input file containing single-end sequencing reads";
        "You can specify more than one FASTQ input, but not FASTQ and FASTA inputs";
        "at the same time. Contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ -> SingleEndFASTQ (TA.get_parameter ()) |> Tools.List.accum Parameters.inputs);
    [ "-p"; "--paired-end" ],
      Some "<fastq_file_name1> <fastq_file_name2>",
      [ "FASTQ input files containing paired-end sequencing reads";
        "You can specify more than one FASTQ input, but not FASTQ and FASTA inputs";
        "at the same time. Contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let name1 = TA.get_parameter () in
        let name2 = TA.get_parameter () in
        PairedEndFASTQ (name1, name2) |> Tools.List.accum Parameters.inputs);
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
      (fun _ -> TA.usage (); exit 1)
  ];
  if !Parameters.verbose then
    TA.header ();
  Parameters.inputs := List.rev !Parameters.inputs;
  if !Parameters.inputs <> [] then begin
    let is_format_fasta = ref false and store = ref Files.ReadsIterate.empty in
    List.iteri
      (fun i input ->
        if i = 0 then
          is_format_fasta := begin
            match input with
            | Files.Type.FASTA _ -> true
            | SingleEndFASTQ _ | PairedEndFASTQ _ -> false
            | _ -> assert false
          end
        else
          if begin
            match input with
            | FASTA _ -> not !is_format_fasta
            | SingleEndFASTQ _ | PairedEndFASTQ _ -> !is_format_fasta
            | _ -> assert false
          end then
            TA.parse_error "You cannot process FASTA and FASTQ inputs together!";        
        store := Files.ReadsIterate.add_from_files !store input)
      !Parameters.inputs;
    if !Parameters.content = None then
      Parameters.content := begin
        if !is_format_fasta then
          Some DNA_ss
        else
          Some DNA_ds
      end;
    begin match !Parameters.content with
    | None -> assert false
    | Some DNA_ss ->
      let module KMC = KMerCounter (KMers.DNAHashSingleStranded (struct let value = !Parameters.k end)) in
      KMC.compute ~linter:(Sequences.Lint.dnaize ~keep_lowercase:false ~keep_dashes:false)
    | Some DNA_ds ->
      let module KMC = KMerCounter (KMers.DNAHashDoubleStranded (struct let value = !Parameters.k end)) in
      KMC.compute ~linter:(Sequences.Lint.dnaize ~keep_lowercase:false ~keep_dashes:false)
    | Some Protein ->
      let module KMC = KMerCounter (KMers.ProteinHash (struct let value = !Parameters.k end)) in
      KMC.compute ~linter:(Sequences.Lint.proteinize ~keep_lowercase:false ~keep_dashes:false)
    end ~verbose:!Parameters.verbose !store !Parameters.max_results_size !Parameters.label !Parameters.output
  end

