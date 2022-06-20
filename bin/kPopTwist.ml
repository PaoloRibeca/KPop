open BiOCamLib

module Defaults =
  struct
    let input = ""
    (* The following three are KPopCountDB.TableFilter.default *)
    let transformation = "normalize"
    let threshold = 1
    let power = 1.
    let sampling = 1.
    (*let precision = 15*)
    let threads = Tools.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let input = ref Defaults.input
    let transformation = ref Defaults.transformation
    let threshold = ref Defaults.threshold
    let power = ref Defaults.power
    let sampling = ref Defaults.sampling
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let version = "0.14"

let header =
  Printf.sprintf begin
    "This is the KPopTwist program (version %s)\n%!" ^^
    " (c) 2022 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!"
  end version

let _ =
  let module TA = Tools.Argv in
  TA.set_header header;
  TA.set_synopsis "-i|--input <input_table_prefix> [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "-f"; "-F"; "-s"; "-S"; "--fraction"; "--sampling"; "--sampling-fraction" ],
      Some "<non_negative_float>",
      [ "fraction of the rows to be considered and resampled before twisting" ],
      TA.Default (fun () -> string_of_float Defaults.sampling),
      (fun _ -> Parameters.sampling := TA.get_parameter_float_non_neg ());
    [ "--threshold" ],
      Some "<non_negative_integer>",
      [ "set to zero all counts that are less than this threshold";
        "before transforming them" ],
      TA.Default (fun () -> string_of_int Defaults.threshold),
      (fun _ -> Parameters.threshold := TA.get_parameter_int_non_neg ());
    [ "--power" ],
      Some "<non_negative_float>",
      [ "raise counts to this power before transforming them.";
        "A power of 0 when the 'pseudocount' method is used";
        "performs a logarithmic transformation" ],
      TA.Default (fun () -> string_of_float Defaults.power),
      (fun _ -> Parameters.power := TA.get_parameter_float_non_neg ());
    [ "--transform"; "--transformation" ],
      Some "'none'|'normalize'|'pseudocount'|'clr'",
      [ "transformation to apply to table elements" ],
      TA.Default (fun () -> Defaults.transformation),
      (fun _ -> Parameters.transformation := TA.get_parameter ());
    TA.make_separator "Input/Output";
    [ "-i"; "--input" ],
      Some "<input_table_prefix>",
      [ "load the specified k-mer database in the register and twist it.";
        "File extension is automatically determined";
        " (will be .KPopCounter).";
        "The prefix is then re-used for output";
        " (and the output file will be given prefix .KPopTwisted)" ],
      TA.Mandatory,
      (fun _ -> Parameters.input := TA.get_parameter ());
    TA.make_separator "Miscellaneous";
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
  Printf.printf "%s\001%s\001%d\001%.12g\001%.12g\001%d\001%b\n%!"
    !Parameters.input !Parameters.transformation !Parameters.threshold !Parameters.power
    !Parameters.sampling !Parameters.threads !Parameters.verbose

