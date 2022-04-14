open BiOCamLib

module Defaults =
  struct
    let input = ""
    let sampling = 1.
    (*let precision = 15*)
    let threads = Tools.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let input = ref Defaults.input
    let sampling = ref Defaults.sampling
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let version = "0.11"

let header =
  Printf.sprintf begin
    "This is the KPopTwist program (version %s)\n%!" ^^
    " (c) 2022 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!"
  end version

let _ =
  let module TA = Tools.Argv in
  let module TS = Tools.Split in
  let module TL = Tools.List in
  TA.set_header header;
  TA.set_synopsis "-i|--input <input_table_prefix> [OPTIONS]";
  TA.parse [
    [], None, [ "Algorithmic parameters:" ], TA.Optional, (fun _ -> ());
    [ "-f"; "-F"; "-s"; "-S"; "--fraction"; "--sampling"; "--sampling-fraction" ],
      Some "<non_negative_float>",
      [ "fraction of the rows to be considered and resampled before twisting" ],
      TA.Default (fun () -> string_of_float Defaults.sampling),
      (fun _ -> Parameters.sampling := TA.get_parameter_float_non_neg ());
    [], None, [ "Input/Output:" ], TA.Optional, (fun _ -> ());
    [ "-i"; "--input" ],
      Some "<input_table_prefix>",
      [ "load the specified k-mer database in the register and twist it.";
        "File extension is automatically determined";
        " (will be .KPopCounter).";
        "The prefix is then re-used for output";
        " (and the output file will be given prefix .KPopTwisted)" ],
      TA.Mandatory,
      (fun _ -> Parameters.input := TA.get_parameter ());
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
  Printf.printf "%s\001%.12g\001%d\001%b\n%!"
    !Parameters.input !Parameters.sampling !Parameters.threads !Parameters.verbose

