open BiOCamLib

module KMerCounter (KMH: KMer.KMerHash with type t = int):
  sig
    val compute: ?verbose:bool -> ?display_step:int -> KMer.ReadStore.t -> int -> string -> string -> unit
  end
= struct
    let compute ?(verbose = false) ?(display_step = 10000) store max_results_size label fname =
      let output =
        if fname = "" then
          stdout
        else
          open_out fname in
      (* Header with label *)
      Printf.fprintf output "\t%s\n" label;
      let output_format = Scanf.format_from_string (Printf.sprintf "%%0%dx\t%%d\n" ((KMH.k + 1) / 2)) "%d%d" in
      let reads_cntr = ref 0 and res = KMer.IntHashtbl.create max_results_size in
      KMer.ReadStore.iter
        (fun _ segm_id read ->
          KMH.iter
            (fun hash occs ->
              KMer.add_to_kmer_counter res hash occs;
              if KMer.IntHashtbl.length res > max_results_size then begin
                let min_binding =
                  KMer.IntHashtbl.fold (fun _ occs old_min -> min old_min !occs) res max_int in
                if verbose then
                  Printf.eprintf "\rKMerCounter.compute: Outputting and removing hashes having #%d occurrences...%!" min_binding;
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
          if !reads_cntr mod display_step = 0 then
            Printf.eprintf "\rKMerCounter.compute: Added %d reads%!" !reads_cntr;
          if segm_id = 0 then
            incr reads_cntr)
        store;
      Printf.eprintf "\rKMerCounter.compute: Added %d reads\n%!" !reads_cntr;
      Printf.eprintf "KMerCounter.compute: Outputting hashes...%!";
      KMer.IntHashtbl.iter (fun hash occs -> Printf.fprintf output output_format hash !occs) res;
      Printf.eprintf " done.%!\n";
      close_out output

  end


module Defaults =
  struct
    let k = 12
    let max_results_size = 16777216 (* Or: 4^12 *)
    let output = ""
    (*let threads = 1*)
    let verbose = false
  end

module Parameters =
  struct
    let k = ref Defaults.k
    let max_results_size = ref Defaults.max_results_size
    let inputs = ref []
    let label = ref ""
    let output = ref Defaults.output
    (*let threads = ref Defaults.threads*)
    let verbose = ref Defaults.verbose
  end

let version = "0.3"

let _ =
  Printf.eprintf "This is the KPopCount program (version %s)\n%!" version;
  Printf.eprintf " (c) 2017-2022 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!";
  let module TA = Tools.Argv in
  let module RS = KMer.ReadStore in
  TA.parse [
    [], None, [ "Algorithmic parameters:" ], TA.Optional, (fun _ -> ());
    [ "-k"; "-K"; "--k-mer-size"; "--k-mer-length" ],
      Some "<k_mer_length>",
      [ "k-mer length"; "(must be positive and <= 30)" ],
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
    [], None, [ "Input/Output:" ], TA.Optional, (fun _ -> ());
    [ "-f"; "-F"; "--fasta" ],
      Some "<fasta_file_name>",
      [ "FASTA input file containing sequences" ],
      TA.Optional,
      (fun _ -> RS.FASTA (TA.get_parameter ()) |> Tools.Misc.accum Parameters.inputs);
    [ "-s"; "-S"; "--single-end" ],
      Some "<fastq_file_name>",
      [ "FASTQ input file containing single-end sequencing reads" ],
      TA.Optional,
      (fun _ -> RS.SingleEndFASTQ (TA.get_parameter ()) |> Tools.Misc.accum Parameters.inputs);
    [ "-p"; "-P"; "--paired-end" ],
      Some "<fastq_file_name1> <fastq_file_name2>",
      [ "FASTQ input files containing paired-end sequencing reads" ],
      TA.Optional,
      (fun _ ->
        let name1 = TA.get_parameter () in
        let name2 = TA.get_parameter () in
        RS.PairedEndFASTQ (name1, name2) |> Tools.Misc.accum Parameters.inputs);
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
    [], None, [ "Miscellaneous:" ], TA.Optional, (fun _ -> ());
(*
    [ "-t"; "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
*)
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (fun () -> string_of_bool !Parameters.verbose),
      (fun _ -> Parameters.verbose := true);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 1)
  ];
  let module KMC = KMerCounter (KMer.EncodingHash (struct let value = !Parameters.k end)) in
  Parameters.inputs := List.rev !Parameters.inputs;
  if !Parameters.inputs <> [] then begin
    let store = ref RS.empty in
    List.iter
      (fun input -> store := RS.add_from_files !store input)
      !Parameters.inputs;
    let store = !store in
    KMC.compute
      ~verbose:!Parameters.verbose store !Parameters.max_results_size !Parameters.label !Parameters.output
  end
