(*
    KPopCount.ml -- (c) 2022-2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    KPopCount is the binary that extracts k-mer spectra from FASTA,
    FASTQ, or tabular inputs and emits a binary `.KPopSpectra`
    database.  It exposes options for k-mer length, read filtering,
    metadata labelling, and per-sequence vs per-batch accumulation.

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

module CoverageFromName =
  struct
    let regexp = Str.regexp {|[+-]?\(\([0-9]+[.]?[0-9]*\)\|\([.][0-9]+\)\)\([eE][-+]?[0-9]+\)?|}
    (* Field number is 1-based! *)
    let extract n_field s =
      let numbers =
        try
          (* In principle this part should be safe *)
          Str.full_split regexp s
            |> List.filter_map
                 (function
                   | Str.Delim s -> Some (float_of_string s)
                   | Text _ -> None)
            |> Array.of_list
        with _ ->
          assert false in
      let n = Array.length numbers in
      if n_field > n then
        Exception.raise __FUNCTION__ IO_Format (Printf.sprintf "Invalid field %d (found %d)" n_field n)
      else
        numbers.(n_field - 1)
  end

module KMerIterator =
  struct
    include KMers.Iterator
    module Hasher =
      struct
        include KMers.Iterator.Hasher
        (* We overwrite the stock library function to better suit this program *)
        let to_string = function
        | K_mers k -> Printf.sprintf "continuous k-mers of size %d" k
        | Gapped (k, g) -> Printf.sprintf "gapped k-mers of size %d (%d+%d+%d)" (2*k+g) k g k
      end
  end
module KMI = KMerIterator

module Action =
  struct
    type t =
      | Empty
      | Of_file of string
      | Set_content of KMI.Content.t
      | Set_hasher of KMI.Hasher.t
      | Set_weight_extractor of int
      (* The string is the label - if empty, each sequence names is treated as a label *)
      | Add_sequences of (string * Files.Reads.t) list
      | To_file of string
    (* Collects together into a list all consecutive inputs having the same format *)
    let compact_add_sequences program s input =
      program :=
        match !program with
        | Add_sequences l :: tl ->
          Add_sequences ((s, input) :: l) :: tl
        | l ->
          Add_sequences [s, input] :: l
  end

module Defaults =
  struct
    let content = KMI.Content.of_string "ds-DNA"
    let hasher = KMI.Hasher.K_mers 12
    let weight_field = 0
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
  Tools.Argv.name = "KPopCount";
  version = "29";
  date = "10-Jan-2026"
} and authors = [
  "2017-2026", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    TA.make_separator_multiline [ ""; "Input/Output of spectra databases:" ];
    [ "-0"; "--zero"; "--empty" ],
      None,
      [ "load an empty database into the register" ],
      TA.Optional,
      (fun _ -> Action.Empty |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load into the register the database present in the specified file";
        " (which must have extension '.KPopSpectra' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Of_file (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "save the database present in the register to the specified file";
        " (which will be given extension '.KPopSpectra' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> To_file (TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Algorithmic parameters:" ];
    [ "-k"; "--k-mer-size"; "--k-mer-length" ],
      Some "<positive_integer>",
      [ "set the hashing strategy to iteration over regular k-mers";
        "and specify the k-mer length to be used.";
        "Options '-k' and '-g' are mutually exclusive; if multiple are specified,";
        "the last one will take effect" ],
      TA.Default
        ((match Defaults.hasher with
          | K_mers _ -> KMI.Hasher.to_string Defaults.hasher
          | Gapped _ -> "not used") |> Fun.const),
      (fun _ -> Set_hasher (KMI.Hasher.K_mers (TA.get_parameter_int_pos ())) |> List.accum Parameters.program);
    [ "-g"; "--gapped-k-mer-sizes"; "--gapped-k-mer-lengths" ],
      Some "BLOCK_SIZE GAP_SIZE",
      [ "where";
        " BLOCK-SIZE := <positive_integer>";
        " GAP-SIZE := <positive_integer>";
        "Set the hashing strategy to iteration over symmetrical gapped k-mers";
        "(having a BLOCK-GAP-BLOCK structure, with BLOCKs of the same size) and";
        "specify their geometry in terms of BLOCK and GAP sizes, respectively.";
        "For instance, option";
        " '-g 5 1'";
        "will iterate on all existing k-mers of size 11 (5+1+5) and not take the";
        "central nucleotide into account for the purpose of computing the hash.";
        "Options '-k' and '-g' are mutually exclusive; if multiple are specified,";
        "the last one will take effect" ],
      TA.Default
        ((match Defaults.hasher with
          | K_mers _ -> "not_used"
          | Gapped _ -> KMI.Hasher.to_string Defaults.hasher) |> Fun.const),
      (fun _ ->
        let k = TA.get_parameter_int_pos () in
        let g = TA.get_parameter_int_pos () in
        Set_hasher (KMI.Hasher.Gapped (k, g)) |> List.accum Parameters.program);
    [ "-c"; "--content" ],
      Some "'ss-DNA'|'single-stranded-DNA'|'ds-DNA'|'double-stranded-DNA'|'protein'|FULL",
      [ "set how contents of following input files should be interpreted.";
        "When content is 'ss-DNA', 'protein' or 'text', only the sequence is hashed;";
        "when content is 'ds-DNA', both sequence and reverse complement are hashed.";
        "'ss-DNA' prevents automatic matching of reverse-complemented sequences;";
        "use it only when comparing a set of single, homogeneus sequences.";
        "These are shortcuts for the full form of this option, which is defined as";
        " FULL := 'DNA('STRANDEDNESS','CASE_SENSITIVITY','UNKNOWN_CHAR_ACTION')'";
        "       | 'protein('UNKNOWN_CHAR_ACTION')'";
        "       | 'text('CASE_SENSITIVITY','UNKNOWN_CHAR_ACTION',";
        "               '<dictionary_file_name>')'";
        "where";
        " STRANDEDNESS := 'ss'|'single-stranded'|'ds'|'double-stranded'";
        " CASE_SENSITIVITY := 'ci'|'case-insensitive'|'cs'|'case-sensitive'";
        " UNKNOWN_CHAR_ACTION := 'split'|'ignore'|'error'";
        "If 'case-insensitive' is specified, DNA/protein sequences are converted to";
        "uppercase characters, while text sequences (and dictionary entries)";
        "are converted to lowercase characters.";
        "UNKNOWN_CHAR_ACTION decides what happens when an unknown character is found";
        "in the input. Option 'split' (the default for DNA/protein sequences) splits";
        "the input sequence and skips unknown characters (for instance N/X) whenever";
        "they are encountered; option 'ignore' (the default for text) silently skips";
        "unknown characters (for instance whitespace); option 'error' causes the";
        "program to abort.";
        "If a dictionary file is specified, each of its lines is interpreted";
        "as a different dictionary entry/token" ],
      TA.Default (KMI.Content.to_string Defaults.content |> Fun.const),
      (fun _ -> Set_content (TA.get_parameter () |> KMI.Content.of_string) |> List.accum Parameters.program);
    [ "-w"; "--weights"; "--weights-from-sequence-names" ],
      Some "<non_negative_integer>",
      [ "given the index n specified as a parameter, extract the n-th number";
        "from each sequence name and weigh the corresponding sequence accordingly.";
        "Indices are 1-based; a value of 0 disables weighting.";
        "If no such field exists, the program will fail.";
        "If the weight is a float number, the ceiling of such number will be used" ],
      TA.Default
        ((if Defaults.weight_field = 0 then
            "do not weigh"
          else
            string_of_int Defaults.weight_field) |> Fun.const),
      (fun _ -> Set_weight_extractor (TA.get_parameter_int_non_neg ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Input/Output of sequences for processing:" ];
    [ "-f"; "--fasta" ],
      Some "<label> <fasta_file_name>",
      [ "process the sequences contained in the specified FASTA input file.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let label = TA.get_parameter () in
        let path = TA.get_parameter () in
        Action.compact_add_sequences Parameters.program label (Files.Reads.FASTA path));
    [ "-s"; "--single-end" ],
      Some "<label> <fastq_file_name>",
      [ "process the sequences contained in the specified FASTQ input file";
        "containing single-end sequencing reads.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let label = TA.get_parameter () in
        let path = TA.get_parameter () in
        Action.compact_add_sequences Parameters.program label (SingleEndFASTQ path));
    [ "-p"; "--paired-end" ],
      Some "<label> <fastq_file_name1> <fastq_file_name2>",
      [ "process the sequences contained in the specified FASTQ input file";
        "containing paired-end sequencing reads.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let label = TA.get_parameter () in
        let path1 = TA.get_parameter () in
        let path2 = TA.get_parameter () in
        Action.compact_add_sequences Parameters.program label (PairedEndFASTQ (path1, path2)));
    [ "-t"; "--tabular" ],
      Some "<label> <tabular_file_name>",
      [ "process the sequences contained in the specified tabular input file.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let label = TA.get_parameter () in
        let path = TA.get_parameter () in
        Action.compact_add_sequences Parameters.program label (Files.Reads.Tabular path));
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
      (fun _ -> TA.usage (); exit 1)
  ];
  let program = List.rev !Parameters.program in
  if program = [] then begin
    TA.usage ();
    exit 0
  end;
  if !Parameters.verbose then
    TA.header ();
  (* These are the registers available to the program *)
  let db = KMerDB.empty () |> ref and content = ref Defaults.content and hasher = ref Defaults.hasher
  and weight_field = ref Defaults.weight_field in
  try
    List.iter
      (function
        | Action.Empty ->
          db := KMerDB.empty ()
        | Of_file prefix ->
          db := KMerDB.of_binary ~verbose:!Parameters.verbose prefix
        | Set_content c ->
          content := c
        | Set_hasher h ->
          hasher := h
        | Set_weight_extractor wf ->
          weight_field := wf
        | Add_sequences files ->
          (* Files will have been inverted during input compaction *)
          let files = Array.of_rlist files in
          let num_files = Array.length files and num_reads = ref 0 in
          if !Parameters.threads = 1 then begin
            let file_cntr = ref 0 and spectrum_cntr = ref 0 and read_label = ref ""            
            and k_mer_iterator, k_mer_finalizer, install_label =
              (* Not for public consumption *)
              let col_idx = -1 |> ref and hashes_to_rows = ref [||] in
              let k_mer_iterator, k_mer_finalizer =
                KMI.make ~verbose:!Parameters.verbose !content !hasher
                  (fun hashes row_idx n ->
                    (* We adjust the size of the translation table as needed *)
                    hashes_to_rows := Array.resize ~is_buffer:true (Array.length hashes) (-1) !hashes_to_rows;
                    if !hashes_to_rows.(row_idx) = -1 then
                      (* The row corresponding to the k-mer in the data matrix
                          has not been allocated yet.
                        Note that the call only ever happens exactly once per row *)
                      !hashes_to_rows.(row_idx) <- KMerDB.add_empty_row_if_needed db hashes.(row_idx);
                    KMerDB.add_counts_unsafe db !col_idx !hashes_to_rows.(row_idx) n) in
              k_mer_iterator, k_mer_finalizer,
              (fun label ->
                col_idx := KMerDB.add_empty_column_if_needed db label) in
            Array.iteri
              (fun file_idx (file_label, input) ->
                let new_file_cntr = file_idx + 1 in
                (* Note that linting is done automatically at a lower level by KMerIterator
                    depending on the sequence type, so we disable it here *)
                Files.Reads.iter ~linter:Sequences.Lint.none ~verbose:false
                  (fun (_, _, read) ->
                    incr num_reads;
                    (* If no global label has been specified, we output one spectrum per sequence *)
                    let new_read_label =
                      if file_label = "" then
                        read.tag
                      else
                        file_label in
                    let is_label_new = new_read_label <> !read_label in
                    (* Remember that there might be more than one read per file *)
                    if new_file_cntr > !file_cntr || is_label_new then begin
                      k_mer_finalizer ();
                      if is_label_new then begin
                        (* Remember that installing the label is expensive,
                            and we want to do it only when necessary *)
                        install_label new_read_label;
                        read_label := new_read_label
                      end
                    end;
                    let weight =
                      if !weight_field = 0 then
                        1.
                      else
                        CoverageFromName.extract !weight_field read.tag in
                    k_mer_iterator ~weight read.seq;
                    let new_spectrum_cntr = !db.core.n_cols in
                    if !Parameters.verbose && begin
                      new_file_cntr > !file_cntr || new_spectrum_cntr > !spectrum_cntr ||
                      !num_reads < 25 || !num_reads mod 25 = 0
                    end then
                      Printf.eprintf "%s\r(%s): Hashed and added %d %s from %d %s and %d %s%!"
                        String.TermIO.clear __FUNCTION__
                        !spectrum_cntr (String.pluralize_int ~plural:"spectra" "spectrum" !spectrum_cntr)
                        !num_reads (String.pluralize_int "read" !num_reads)
                        (new_file_cntr - 1) (String.pluralize_int "file" (new_file_cntr - 1));
                    file_cntr := new_file_cntr;
                    spectrum_cntr := new_spectrum_cntr)
                  input;
                (* We need to make sure we process all buffers here *)
                k_mer_finalizer ())
              files
          end else begin
            (* The READER's variables.
               Here we use C++-style iterators to parallelise over sequences
                rather than files, which should give a better granularity.
               We also try to balance load by accumulating reads up to the maximum
                sequence length seen by the reader so far *)
            let module Iter = Files.Reads.Iterator in
            let file_idx = -1 |> ref and file_label = ref ""
            and it = Iter.empty () |> ref and max_len = ref 8192
            (* Everyone else's variables *)
            and file_cntr = ref 0 and spectrum_cntr = ref 0 and read_label = ref ""
            and k_mer_iterator, dump_counts =
              (* Not for public consumption.
                 The hash set always grows and never gets passed around - we only transmit differences *)
              let old_num_hashes = ref 0 and hashes = ref [||] in
              (* These are stacks backed by bigarrays so that passing them around is easier *)
              let module IntBAS = Numbers.IntBAVectorStack in
              let module F32BAS = Numbers.F32BAVectorStack in
              let keys = IntBAS.empty () and vals = F32BAS.empty () in
              let k_mer_iterator, k_mer_finalizer =
                KMI.make ~verbose:!Parameters.verbose !content !hasher
                  (fun hs row_idx n ->
                    hashes := hs;
                    IntBAS.push keys row_idx;
                    F32BAS.push vals n) in
              k_mer_iterator,
              (fun file_cntr read_label res ->
                IntBAS.clear keys;
                F32BAS.clear vals;
                k_mer_finalizer ();
                (* We send new hashes and counts associated with the label *)
                let keys = IntBAS.contents keys and vals = F32BAS.contents vals in
                let num_hashes = Array.length !hashes in
                (* Note that it is OK for some labels not to have any associated k-mers,
                    so in such cases we have to send the label anyway *)
                if Numbers.IntBAVector.length keys > 0 || read_label <> "" then
                  Tools.ArrayStack.push res begin
                    file_cntr, read_label, !old_num_hashes,
                    Array.sub !hashes !old_num_hashes (num_hashes - !old_num_hashes),
                    keys, vals
                  end;
                old_num_hashes := num_hashes)
            and hashes_to_rows_table = IntHashtbl.create !Parameters.threads in
            Processes.Parallel.process_stream_chunkwise
              (fun () ->
                assert (!file_idx <= num_files);
                let acc_len = ref 0 and res = ref [] in
                (* Still some files to process *)
                while Iter.is_empty !it && !file_idx < num_files do
                  Iter.delete !it;
                  incr file_idx;
                  if !file_idx < num_files then begin
                    (* If we are done with the current file, let's open the next one *)
                    let new_file_label, input = files.(!file_idx) in
                    file_label := new_file_label;
                    (* Note that linting is done automatically at a lower level by KMerIterator
                        depending on the sequence type, so we disable it here *)
                    it := Iter.create (Sequences.Lint.none ~keep_lowercase:true ~keep_dashes:true, input)
                  end
                done;
                if !file_idx < num_files then begin
                  (* If we got here, the iterator is valid *)
                  while Iter.is_empty !it |> not && !acc_len < !max_len do
                    Iter.get_and_incr !it
                      (fun (_, _, read) ->
                        let len = String.length read.seq in
                        max_len := max !max_len len;
                        Numbers.Int.(acc_len ++ len);
                        (* If no global label has been specified, we output one spectrum per sequence *)
                        let label =
                          if !file_label = "" then
                            read.tag
                          else
                            !file_label in
                        List.accum res (!file_idx + 1, label, read))
                  done;
                  !res
                end else
                  raise End_of_file)
              (fun jobs ->
                (* Keeping the original order is important here, as the encoding of hashes
                    is transmitted incrementally *)
                let num_reads = ref 0 and res = Tools.ArrayStack.empty () in                
                List.iter
                  (fun (new_file_cntr, new_label, read) ->
                    if new_file_cntr > !file_cntr || new_label <> !read_label then begin
                      dump_counts !file_cntr !read_label res;
                      file_cntr := new_file_cntr;
                      read_label := new_label
                    end;
                    let weight =
                      if !weight_field = 0 then
                        1.
                      else
                        CoverageFromName.extract !weight_field read.Files.Base.Read.tag in
                    k_mer_iterator ~weight read.seq;
                    incr num_reads)
                  jobs;
                dump_counts !file_cntr !read_label res;
                Unix.getpid (), !num_reads, Tools.ArrayStack.contents res)
              (fun (pid, n_reads, res) ->
                Array.iter
                  (fun (new_file_cntr, label, old_num_hashes, new_hashes, keys, vals) ->
                    file_cntr := max !file_cntr new_file_cntr;
                    (* We update the hashes for the worker process that sent the result *)
                    let hashes_to_rows =
                      match IntHashtbl.find_opt hashes_to_rows_table pid with
                      | Some hashes_to_rows ->
                        hashes_to_rows
                      | None ->
                        let res = ref [||] in
                        IntHashtbl.add hashes_to_rows_table pid res;
                        res in
                    let num_new_hashes = Array.length new_hashes in
                    (* We resize the vector to the size provided by the subprocess *)
                    hashes_to_rows :=
                      Array.resize ~is_buffer:true (old_num_hashes + num_new_hashes) (-1) !hashes_to_rows;
                    Array.iteri
                      (fun i hash ->
                        let i = old_num_hashes + i in
                        if !hashes_to_rows.(i) = -1 then
                          !hashes_to_rows.(i) <- KMerDB.add_empty_row_if_needed db hash)
                      new_hashes;
                    (* Which spectrum are we updating? *)
                    let new_spectrum_cntr = !db.core.n_cols in
                    if !Parameters.verbose && begin
                      new_spectrum_cntr > !spectrum_cntr || !num_reads < 25 || !num_reads mod 25 = 0
                    end then
                      Printf.eprintf "%s\r(%s): Hashed and added %d %s from %d %s and %d %s%!"
                        String.TermIO.clear __FUNCTION__
                        new_spectrum_cntr (String.pluralize_int ~plural:"spectra" "spectrum" new_spectrum_cntr)
                        !num_reads (String.pluralize_int "read" !num_reads)
                        (!file_cntr - 1) (String.pluralize_int "file" (!file_cntr - 1));
                    spectrum_cntr := new_spectrum_cntr;
                    let col_idx = KMerDB.add_empty_column_if_needed db label in
                    (* We iterate on the list of k-mers *)
                    assert (Numbers.IBAVector.length keys = Numbers.F32BAVector.length vals);
                    Numbers.IntBAVector.iteri
                      (fun i key ->
                        KMerDB.add_counts_unsafe db col_idx !hashes_to_rows.(key) Numbers.F32BAVector.(vals.@(i)))
                      keys)
                  res;
                Numbers.Int.(num_reads ++ n_reads))
              !Parameters.threads;
          end;
          if !Parameters.verbose then
            let spectrum_cntr = !db.core.n_cols in
            Printf.eprintf "%s\r(%s): Hashed and added %d %s from %d %s and %d %s.\n%!"
              String.TermIO.clear __FUNCTION__
              spectrum_cntr (String.pluralize_int ~plural:"spectra" "spectrum" spectrum_cntr)
              !num_reads (String.pluralize_int "read" !num_reads)
              num_files (String.pluralize_int "file" num_files)
        | To_file prefix ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () -> KMerDB.to_binary ~verbose:!Parameters.verbose !db prefix))
      program;
    if !Parameters.verbose && !Parameters.threads = 1 then
      Printf.eprintf "(%s): Timers: (iterate=%g (lint=%g, encode=%g, accumulate=%g), finalize=%g)\n%!" __FUNCTION__
        (Tools.Timer.read "KMers.Iterator.Iterator") (Tools.Timer.read "KMers.Iterator.Linter")
        (Tools.Timer.read "KMers.Iterator.Encoder") (Tools.Timer.read "KMers.Iterator.Accumulator")
        (Tools.Timer.read "KMers.Iterator.Finalizer")
  with e ->
    Exception.handle __FUNCTION__ TA.usage (fun () ->
      Printf.peprintf "(%s): This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!" __FUNCTION__;
      Printf.peprintf "(%s): You might also wish to rerun me with option -x to get a full backtrace.\n%!" __FUNCTION__
    ) e

