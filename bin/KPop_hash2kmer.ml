(*
    KPop-hash2kmer.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    KPop-hash2kmer is a small helper that back-translates the k-mer hashes
    produced by KPopCount (one per line, as found in the row names of a
    '.KPopSpectra' database) into the corresponding sequences.  It accepts
    the subset of KPopCount options that govern hashing, reads one hash per
    line from standard input, and writes the back-translated sequence to
    standard output, flushing after each line.

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

module Defaults =
  struct
    let content = KMI.Content.of_string "ds-DNA"
    let hasher = KMI.Hasher.K_mers 12
    let gap = '-'
    let verbose = false
  end

module Parameters =
  struct
    let content = ref Defaults.content
    let hasher = ref Defaults.hasher
    let gap = ref Defaults.gap
    let verbose = ref Defaults.verbose
  end

let info = {
  Tools.Argv.name = "KPop-hash2kmer";
  version = "1";
  date = "29-Jun-2026"
} and authors = [
  "2026", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[OPTIONS]";
  TA.parse [
    TA.make_separator_multiline
      [ "Back-translate into sequences the k-mer hashes produced by KPopCount.";
        "Hashes are read one per line from standard input (they are the row names";
        "of a '.KPopSpectra' database); the corresponding sequences are written";
        "one per line to standard output, which is flushed after each line.";
        "The hashing parameters below must match those used to produce the hashes." ];
    TA.make_separator_multiline [ "Algorithmic parameters:" ];
    [ "-k"; "--k-mer-size"; "--k-mer-length" ],
      Some "<positive_integer>",
      [ "set the hashing strategy to iteration over regular k-mers";
        "and specify the k-mer length that was used to produce the hashes.";
        "Options '-k' and '-g' are mutually exclusive; if multiple are specified,";
        "the last one will take effect" ],
      TA.Default
        ((match Defaults.hasher with
          | K_mers _ -> KMI.Hasher.to_string Defaults.hasher
          | Gapped _ -> "not used") |> Fun.const),
      (fun _ -> Parameters.hasher := KMI.Hasher.K_mers (TA.get_parameter_int_pos ()));
    [ "-g"; "--gapped-k-mer-sizes"; "--gapped-k-mer-lengths" ],
      Some "BLOCK_SIZE GAP_SIZE",
      [ "where";
        " BLOCK-SIZE := <positive_integer>";
        " GAP-SIZE := <positive_integer>";
        "set the hashing strategy to iteration over symmetrical gapped k-mers";
        "(having a BLOCK-GAP-BLOCK structure, with BLOCKs of the same size) and";
        "specify their geometry in terms of BLOCK and GAP sizes, respectively.";
        "The GAP residues are not hashed; option '-G' sets how they are rendered.";
        "Options '-k' and '-g' are mutually exclusive; if multiple are specified,";
        "the last one will take effect" ],
      TA.Default
        ((match Defaults.hasher with
          | K_mers _ -> "not used"
          | Gapped _ -> KMI.Hasher.to_string Defaults.hasher) |> Fun.const),
      (fun _ ->
        let k = TA.get_parameter_int_pos () in
        let g = TA.get_parameter_int_pos () in
        Parameters.hasher := KMI.Hasher.Gapped (k, g));
    [ "-c"; "--content" ],
      Some "'ss-DNA'|'single-stranded-DNA'|'ds-DNA'|'double-stranded-DNA'|'protein'|FULL",
      [ "set how the hashes should be interpreted, i.e. which alphabet is used";
        "to back-translate them. Note that for 'ds-DNA' the hashes are";
        "strand-canonical, and hence are back-translated to the smaller of each";
        "k-mer and its reverse complement; 'ss-DNA' and 'ds-DNA' are otherwise";
        "equivalent here, as are the strandedness and unknown-character settings";
        "of the full form, which do not affect back-translation.";
        "These are shortcuts for the full form of this option, which is defined as";
        " FULL := 'DNA('STRANDEDNESS','CASE_SENSITIVITY','UNKNOWN_CHAR_ACTION')'";
        "       | 'protein('UNKNOWN_CHAR_ACTION')'";
        "       | 'text('CASE_SENSITIVITY','UNKNOWN_CHAR_ACTION',";
        "               '<dictionary_file_name>')'";
        "where";
        " STRANDEDNESS := 'ss'|'single-stranded'|'ds'|'double-stranded'";
        " CASE_SENSITIVITY := 'ci'|'case-insensitive'|'cs'|'case-sensitive'";
        " UNKNOWN_CHAR_ACTION := 'split'|'ignore'|'error'";
        "If a dictionary file is specified, each of its lines is interpreted";
        "as a different dictionary entry/token" ],
      TA.Default (KMI.Content.to_string Defaults.content |> Fun.const),
      (fun _ -> Parameters.content := TA.get_parameter () |> KMI.Content.of_string);
    [ "-G"; "--gap"; "--gap-character" ],
      Some "<character>",
      [ "set the single character used to render each unhashed gap position";
        "when back-translating gapped k-mers produced with option '-g'; it has";
        "no effect otherwise" ],
      TA.Default (Printf.sprintf "%c" Defaults.gap |> Fun.const),
      (fun _ ->
        let s = TA.get_parameter () in
        if String.length s <> 1 then
          Exception.raise __FUNCTION__ IO_Format
            (Printf.sprintf "The gap specifier must be a single character (found '%s')" s);
        Parameters.gap := s.[0]);
    TA.make_separator_multiline [ "Miscellaneous options:" ];
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
    [ "-x"; "--print-exception-backtrace" ], None, [], TA.Optional,
      (fun _ -> Printexc.record_backtrace true);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 1)
  ];
  if !Parameters.verbose then
    TA.header ();
  try
    let decode =
      KMI.Decoder.make ~verbose:!Parameters.verbose ~gap:!Parameters.gap
        !Parameters.content !Parameters.hasher in
    while true do
      input_line stdin |> decode |> Printf.printf "%s\n%!"
    done
  with
  | End_of_file -> ()
  | e ->
    Exception.handle __FUNCTION__ TA.usage (fun () ->
      Printf.peprintf
        "(%s): This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!"
        __FUNCTION__;
      Printf.peprintf
        "(%s): You might also wish to rerun me with option -x to get a full backtrace.\n%!"
        __FUNCTION__
    ) e

