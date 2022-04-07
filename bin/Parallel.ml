open BiOCamLib

module Misc =
  struct

(* Should be moved to Tools.ml ? *)

    let input_line_num = ref 0
    let line__failwith s =
      Printf.eprintf "On line %d: %s\n%!" !input_line_num s;
      exit 1

  end


module Defaults =
  struct
    let lines_per_block = 10000
    let threads =
      try
        Tools.Subprocess.spawn_and_read_single_line "nproc" |> int_of_string
      with _ ->
        1
    let verbose = false
    let debug = false
  end

module Parameters =
  struct
    let command = ref ""
    let args = ref [||]
    let lines_per_block = ref Defaults.lines_per_block
    let input = ref ""
    let output = ref ""
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
    let debug = ref Defaults.debug
  end

let version = "0.1"

let () =
  let module TA = Tools.Argv in
  TA.set_synopsis "[OPTIONS] -- [COMMAND TO PARALLELIZE AND ITS OPTIONS]";
  TA.parse [
    [], None, [ "Command to parallelize" ], TA.Optional, (fun _ -> ());
    [ "--" ],
      None,
      [ "consider all the subsequent parameters";
        "as the command to be executed in parallel.";
        "At least one command must be specified" ],
      TA.Mandatory,
      (fun _ ->
        Parameters.command := TA.get_parameter ();
        Parameters.args := Array.append [| !Parameters.command |] (TA.get_remaining_parameters ()));
    [], None, [ "Input/Output" ], TA.Optional, (fun _ -> ());
    [ "--lines-per-block" ],
      Some "<positive_integer>",
      [ "number of lines to be processed per block" ],
      TA.Default (fun () -> string_of_int !Parameters.lines_per_block),
      (fun _ -> Parameters.lines_per_block := TA.get_parameter_int_pos ());
    [ "-i"; "--input" ],
      Some "<input_file>",
      [ "name of input file" ],
      TA.Default (fun () -> "stdin"),
      (fun _ -> Parameters.input := TA.get_parameter ());
    [ "-o"; "--output" ],
      Some "<output_file>",
      [ "name of output file" ],
      TA.Default (fun () -> "stdout"),
      (fun _ -> Parameters.output := TA.get_parameter ());
    [], None, [ "Miscellaneous" ], TA.Optional, (fun _ -> ());
    [ "-t"; "--threads" ],
      Some "<positive_integer>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (fun () -> string_of_bool !Parameters.verbose),
      (fun _ -> Parameters.verbose := true);
    [ "-d"; "--debug" ],
      None,
      [ "output debugging information" ],
      TA.Default (fun () -> string_of_bool !Parameters.debug),
      (fun _ -> Parameters.debug := true);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 1)
  ];
  let verbose = !Parameters.verbose and debug = !Parameters.debug in
  if verbose then begin
    Printf.eprintf "This is the Parallel program (version %s)\n%!" version;
    Printf.eprintf " (c) 2019 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!"
  end;
  let input =
    if !Parameters.input = "" then
      stdin
    else
      open_in !Parameters.input
  and output =
    if !Parameters.output = "" then
      stdout
    else
      open_out !Parameters.output in
  let print_num_lines what =
    if verbose then
      Printf.sprintf "%d %s %s" !Misc.input_line_num (Tools.Misc.pluralize_int "line" !Misc.input_line_num) what |>
      Tools.Printf.pteprintf "%s\n%!"
  and wait_and_check pid =
    match let _, status = Unix.waitpid [] pid in status with
    | Unix.WEXITED 0 -> ()
    | _ ->
      Printf.sprintf "[%07d] !!! Subprocess %d terminated with an error" (Unix.getpid ()) pid |>
      Misc.line__failwith in
  let lines_per_block = !Parameters.lines_per_block
  and eof = ref false and processing_buffer = Buffer.create 16777216 in
  Tools.Parallel.process_stream_chunkwise
    (fun () ->
      if not !eof then begin
        if !Misc.input_line_num mod lines_per_block = 0 then
          print_num_lines "read";
        let buf = ref [] and base_input_line_num = !Misc.input_line_num + 1 in
        begin try
          for _ = 1 to lines_per_block do
            incr Misc.input_line_num;
            input_line input |> Tools.Misc.accum buf
          done
        with End_of_file ->
          decr Misc.input_line_num;
          print_num_lines "read";
          eof := true
        end;
        base_input_line_num, !buf
      end else
        raise End_of_file)
    (fun (base_input_line_num, input_lines) ->
      (* We want to keep counters realiable in case of error *)
      Misc.input_line_num := base_input_line_num;
      let input_lines_processed = ref 0 in
      Buffer.clear processing_buffer;
      (* We start the subprocess *)
      let in_pipe_out, in_pipe_in = Unix.pipe ()
      and out_pipe_out, out_pipe_in = Unix.pipe () in
      match Unix.fork () with
      | 0 -> (* Child *)
        Unix.close in_pipe_in;
        Unix.close out_pipe_out;
        (* We connect the first pipe to stdin and the second pipe to stdout *)
        Unix.dup2 in_pipe_out Unix.stdin;
        Unix.dup2 out_pipe_in Unix.stdout;
        Unix.close in_pipe_out;
        Unix.close out_pipe_in;
        Unix.execvp !Parameters.command !Parameters.args
        (* We will never get here *)
      | pid_child -> (* Parent *)
        Unix.close in_pipe_out;
        Unix.close out_pipe_in;
        (* We communicate with the subprocess through the pipes *)
        let process_input = Unix.out_channel_of_descr in_pipe_in
        and process_output = Unix.in_channel_of_descr out_pipe_out in
(* Is this really necessary?
        set_binary_mode_in process_output false;
        set_binary_mode_out process_input false;
*)
        match Unix.fork () with
        | 0 -> (* Grandchild *)
          close_in process_output;
          (* We feed the input lines to the subprocess *)
          List.iter
            (fun line ->
              if debug then
                Tools.Printf.peprintf ">>> Pushing line '%s'\n%!" line;
              output_string process_input line;
              output_char process_input '\n')
            (List.rev input_lines);
            flush process_input;
          close_out process_input;
  (*        Unix.close in_pipe_in; *)
          if debug then
            Tools.Printf.peprintf ">>> Pushing done.\n%!";
          (* Commit suicide *)
          exit 0
        | pid_grandchild -> (* Grandparent *)
          close_out process_input;
          (* We have to update the counter here, otherwise the information will be lost *)
          input_lines_processed := !input_lines_processed + List.length input_lines;
          (* We read the output from the subprocess, and accumulate it into the buffer.
             Note that in principle the number of lines in output
             could differ from the number of lines in input *)
          begin try
            while true do
              let line = input_line process_output in
              if debug then
                Tools.Printf.peprintf "<<< Pulling line '%s'\n%!" line;
              Buffer.add_string processing_buffer line;
              Buffer.add_char processing_buffer '\n'
            done
          with End_of_file -> ()
          end;
          if debug then
            Tools.Printf.peprintf "<<< Pulling done.\n%!";
  (*        Unix.close out_pipe_out; *)
          close_in process_output;
          (* We make sure that the children are not zombying around *)
          wait_and_check pid_child;
          wait_and_check pid_grandchild;
          !input_lines_processed, Buffer.contents processing_buffer)
    (fun (processed, buf) ->
      if !Misc.input_line_num = 0 then
        print_num_lines "processed";
      let old_periods = !Misc.input_line_num / lines_per_block in
      Misc.input_line_num := !Misc.input_line_num + processed;
      if !Misc.input_line_num / lines_per_block > old_periods then
        print_num_lines "processed";
      Printf.fprintf output "%s%!" buf)
    !Parameters.threads;
  print_num_lines "processed";
  (* Cleanup actions *)
  close_in input;
  close_out output
