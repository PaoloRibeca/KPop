(*
    KPop_autotuner.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    Find a partition of a spectrum database without being told what to look for.

    This is the whole de-novo pipeline as one program: a projection to embed the
    spectra in, a Monte-Carlo search for the partition that best fits them, and
    then the loop that feeds the partition back to define the projection and
    searches again, until the two agree.  Each half already existed --
    KPopTwist builds a twister, KPopTwistDB clusters a twisted register -- and
    what did not was the loop between them, which spans three registers and so
    belonged in a program of its own rather than inside either of theirs.

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

module Defaults =
  struct
    let projection_sample = 0
    let iterations = 12
    let report_anyway = false
    let metric = Space.Distance.Metric.of_string "powers(1,1,1)"
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalize = false
    let combination_criterion = KMerDB.CombinationCriterion.of_string "mean"
    let montecarlo_replicas = 4
    let montecarlo_steps = 2000
    let montecarlo_sample = 250
    let montecarlo_temperature = 0.002
    let montecarlo_decades = 3.
    let montecarlo_cooling = 0.999
    let seed = 17
    let threads = Processes.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let input = ref ""
    let output = ref ""
    let projection_sample = ref Defaults.projection_sample
    let iterations = ref Defaults.iterations
    let report_anyway = ref Defaults.report_anyway
    let metric = ref Defaults.metric
    let distance = ref Defaults.distance
    let distance_normalize = ref Defaults.distance_normalize
    let combination_criterion = ref Defaults.combination_criterion
    let montecarlo_replicas = ref Defaults.montecarlo_replicas
    let montecarlo_steps = ref Defaults.montecarlo_steps
    let montecarlo_sample = ref Defaults.montecarlo_sample
    let montecarlo_temperature = ref Defaults.montecarlo_temperature
    let montecarlo_decades = ref Defaults.montecarlo_decades
    let montecarlo_cooling = ref Defaults.montecarlo_cooling
    let seed = ref Defaults.seed
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let info = { Tools.Argv.name = "KPop-autotuner"; version = Info.version; date = Info.date }
and authors = [
  "2026", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

(* THE CANONICAL FORM OF A PARTITION, which is what the loop compares to decide it has
   finished.  Two partitions are the same partition when they induce the same equivalence
   relation, whatever they call their classes, and the labels here are the names of
   representative samples and so change whenever the representative does -- comparing them
   directly would report a difference every round.  Mapping each element to the position of the
   first element of its own class erases the naming and leaves the relation *)
let canonical assign =
  let n = Array.length assign in
  let first = Hashtbl.create 64 and canon = Array.make n 0 in
  Array.iteri
    (fun i a ->
      canon.(i) <-
        (match Hashtbl.find_opt first a with
        | Some j -> j
        | None -> Hashtbl.add first a i; i))
    assign;
  canon

(* How much larger than the sample the pool it is chosen from is.  Four is enough for the choice
   to have something to reject and small enough that analysing the pool stays cheap *)
let pool_factor = 4

(* THE SAMPLE THE FIRST PROJECTION IS BUILT FROM, chosen to SPAN the corpus rather than to
   represent it -- and those are opposite things here.  A corpus whose classes differ greatly in
   size is mostly its big ones, half of the norovirus VP2 set being a single genotype, so a
   uniform sample is mostly that genotype too and the analysis of it spans that genotype's
   internal variation instead of what tells one class from another, which is the structure the
   search is looking for.  Measured on that corpus, a uniform sample of 150 gives a projection
   in which the whole set shows no groups at all, where a spanning sample of about the same size
   gives one in which it falls into 87.
   So a larger pool is drawn at random, analysed on its own, and the sample taken from it by
   farthest-point: the first member is arbitrary and each one after it is whichever of the pool
   lies furthest from everything taken so far.  The pool's own projection is skewed in exactly
   the same way, and does not need not to be: all it has to show is that an unusual sequence
   lies far from the bulk, which a projection dominated by the bulk shows perfectly well.
   Returns the spectra the sample does NOT name, those being what has to go for the database to
   become the sample, [remove_selected] being the way to keep a subset *)
let complement_of_diverse_sample ~threads ~verbose state db n =
  let names = db.KMerDB_Base.core.KMerDB_Base.idx_to_col_names
  and n_cols = db.KMerDB_Base.core.KMerDB_Base.n_cols in
  (* A pass over the whole set taking each with the probability that leaves the selection the
     size it should be, so that no name can be drawn twice and the cost does not depend on how
     much of the set is wanted *)
  let take_random m =
    let picked = Array.make n_cols false and left = ref (min m n_cols) and rest = ref n_cols in
    for i = 0 to n_cols - 1 do
      if Random.State.int state !rest < !left then begin picked.(i) <- true; decr left end;
      decr rest
    done;
    picked in
  let complement_of picked =
    let acc = ref StringSet.empty in
    for i = 0 to n_cols - 1 do
      if not picked.(i) then acc := StringSet.add names.(i) !acc
    done;
    !acc in
  if n * pool_factor >= n_cols then
    (* The pool would be the whole database, and analysing all of it to choose part of it costs
       what analysing all of it costs, which is what taking a sample was avoiding *)
    complement_of (take_random n)
  else begin
    let pool = KMerDB_Base.remove_selected db (complement_of (take_random (n * pool_factor))) in
    let _, pool_twisted, _ = CA.twist ~threads ~verbose pool in
    let mat = pool_twisted.Twisted.twisted.Matrix.matrix in
    let coords = mat.Matrix.Base.data and pool_names = mat.Matrix.Base.row_names in
    let m = Array.length pool_names in
    let d = if m > 0 then Float.Array.length coords.(0) else 0 in
    (* IN PRINCIPAL COORDINATES AND NOT IN THE STORED STANDARD ONES.  A pool of p spectra is
       analysed into p-1 axes of which only the leading few carry any structure, and standard
       coordinates weigh all of them alike -- so a distance taken over them is mostly noise,
       and farthest-point, which seeks whatever is most extreme, then selects whatever is most
       extreme in that noise.  Measured on norovirus VP1, a sample chosen that way yields a
       projection in which the corpus shows no groups at all.  Weighting each axis by its
       inertia puts the structure back in charge of the distance, which is the same reason
       every clustering here is done in principal coordinates *)
    let iv = pool_twisted.Twisted.inertia.Matrix.matrix.Matrix.Base.data.(0) in
    let w = Float.Array.init d (fun k -> sqrt (Float.Array.get iv k)) in
    let dist i j =
      let acc = ref 0. in
      for k = 0 to d - 1 do
        let v =
          (Float.Array.get coords.(i) k -. Float.Array.get coords.(j) k)
            *. Float.Array.get w k in
        acc := !acc +. (v *. v)
      done;
      sqrt !acc in
    (* [near] is each candidate's distance to the nearest member already taken, which is what
       farthest-point maximises, and it is updated rather than recomputed: adding a member can
       only bring a candidate's nearest closer *)
    let chosen = Array.make m false and near = Array.make m infinity in
    let cur = ref (if m > 0 then Random.State.int state m else 0) in
    for _ = 1 to min n m do
      chosen.(!cur) <- true;
      let best = ref (-1) and best_d = ref neg_infinity in
      for i = 0 to m - 1 do
        if not chosen.(i) then begin
          let dd = dist i !cur in
          if dd < near.(i) then near.(i) <- dd;
          if near.(i) > !best_d then begin best_d := near.(i); best := i end
        end
      done;
      if !best >= 0 then cur := !best
    done;
    let keep = ref StringSet.empty in
    Array.iteri (fun i nm -> if chosen.(i) then keep := StringSet.add nm !keep) pool_names;
    let acc = ref StringSet.empty in
    for i = 0 to n_cols - 1 do
      if not (StringSet.mem names.(i) !keep) then acc := StringSet.add names.(i) !acc
    done;
    !acc
  end

let () =
  let module TA = Tools.Argv in
  let prefix = Printf.sprintf "(%s):" __FUNCTION__ in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "-i <binary_file_prefix> -o <output_file_prefix> [OPTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Input/Output."; "" ];
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load the specified binary database of k-mer spectra.";
        "File extension is automatically assigned (will be '.KPopSpectra')" ],
      TA.Mandatory,
      (fun _ -> Parameters.input := TA.get_parameter ());
    [ "-o"; "--output" ],
      Some "<output_file_prefix>",
      [ "prefix for the three results, their extensions being assigned rather than typed:";
        "the twister the loop settled on ('.KPopTwister'), the spectra as it embeds them";
        "('.KPopTwisted'), and the partition ('.KPopClasses.txt'), the last in the two-line";
        "form 'KPopCountDB -m <file> -c CLASS' reads" ],
      TA.Mandatory,
      (fun _ -> Parameters.output := TA.get_parameter ());
    TA.make_separator_multiline [ "Algorithm."; "" ];
    [ "--projection-sample" ],
      Some "<non_negative_integer>",
      [ "number of spectra defining the initial projection, or 0 to use every one of them.";
        "A sample is usually enough, and is much cheaper than a correspondence analysis of";
        "the whole database, because the projection only has to keep the classes apart and";
        "their arrangement occupies far fewer dimensions than the database is stored in.";
        "Note that the Johnson-Lindenstrauss bound does NOT justify this: being";
        "distribution-free it asks for more dimensions than there are samples here, and what";
        "makes a small sample work is the low intrinsic dimension of this particular data";
        "rather than any worst-case guarantee.  It matters only for the first round: from the";
        "second the axes come from the partition rather than from a sample" ],
      TA.Default (string_of_int Defaults.projection_sample |> Fun.const),
      (fun _ -> Parameters.projection_sample := TA.get_parameter_int_non_neg ());
    [ "--iterations" ],
      Some "<positive_integer>",
      [ "maximum number of rounds.  The loop stops earlier, and usually does, when a round";
        "returns the partition it was given -- that being a fixed point of the whole";
        "procedure and the answer it is looking for" ],
      TA.Default (string_of_int Defaults.iterations |> Fun.const),
      (fun _ -> Parameters.iterations := TA.get_parameter_int_pos ());
    [ "--report-anyway" ],
      None,
      [ "write a partition out even when the samples show no group structure.  Off by default,";
        "because a partition of a set that does not cluster says more about the criterion that";
        "produced it than about the data, and the one thing this program should not do is hand";
        "one over without saying so" ],
      TA.Default (fun () -> string_of_bool Defaults.report_anyway),
      (fun _ -> Parameters.report_anyway := true);
    [ "-m"; "--metric" ],
      Some "<metric>",
      [ "metric the search measures distances under" ],
      TA.Default (Space.Distance.Metric.to_string Defaults.metric |> Fun.const),
      (fun _ -> Parameters.metric := TA.get_parameter () |> Space.Distance.Metric.of_string);
    [ "-d"; "--distance" ],
      Some "<distance>",
      [ "distance the search measures distances with" ],
      TA.Default (Space.Distance.to_string Defaults.distance |> Fun.const),
      (fun _ -> Parameters.distance := TA.get_parameter () |> Space.Distance.of_string);
    [ "--distance-normalize" ],
      Some "<bool>",
      [ "whether to normalise vectors before measuring the distance between them" ],
      TA.Default (string_of_bool Defaults.distance_normalize |> Fun.const),
      (fun _ -> Parameters.distance_normalize := TA.get_parameter_boolean ());
    [ "--combination-criterion" ],
      Some "'mean'|'median'",
      [ "criterion combining the spectra of a class into the one representing it when the";
        "partition is fed back to define the next projection" ],
      TA.Default (KMerDB.CombinationCriterion.to_string Defaults.combination_criterion
                  |> Fun.const),
      (fun _ ->
        Parameters.combination_criterion :=
          TA.get_parameter () |> KMerDB.CombinationCriterion.of_string);
    TA.make_separator_multiline
      [ "The search."; "These mirror the '--clusters-montecarlo-*' family of KPopTwistDB." ];
    [ "--montecarlo-replicas" ],
      Some "<positive_integer>",
      [ "number of chains run in parallel at a ladder of temperatures" ],
      TA.Default (string_of_int Defaults.montecarlo_replicas |> Fun.const),
      (fun _ -> Parameters.montecarlo_replicas := TA.get_parameter_int_pos ());
    [ "--montecarlo-steps" ],
      Some "<positive_integer>",
      [ "number of moves each round of the search attempts" ],
      TA.Default (string_of_int Defaults.montecarlo_steps |> Fun.const),
      (fun _ -> Parameters.montecarlo_steps := TA.get_parameter_int_pos ());
    [ "--montecarlo-sample" ],
      Some "<positive_integer>",
      [ "number of points a chain scores a partition on" ],
      TA.Default (string_of_int Defaults.montecarlo_sample |> Fun.const),
      (fun _ -> Parameters.montecarlo_sample := TA.get_parameter_int_pos ());
    [ "--montecarlo-temperature" ],
      Some "<positive_float>",
      [ "coldest rung of the temperature ladder, which retunes itself from there" ],
      TA.Default (string_of_float Defaults.montecarlo_temperature |> Fun.const),
      (fun _ -> Parameters.montecarlo_temperature := TA.get_parameter_float_pos ());
    [ "--montecarlo-decades" ],
      Some "<positive_float>",
      [ "number of decades of temperature the ladder spans" ],
      TA.Default (string_of_float Defaults.montecarlo_decades |> Fun.const),
      (fun _ -> Parameters.montecarlo_decades := TA.get_parameter_float_pos ());
    [ "--montecarlo-cooling" ],
      Some "<fractional_float>",
      [ "factor the temperature is multiplied by after each move" ],
      TA.Default (string_of_float Defaults.montecarlo_cooling |> Fun.const),
      (fun _ -> Parameters.montecarlo_cooling := TA.get_parameter_float_fraction ());
    TA.make_separator_multiline [ "Miscellaneous."; "" ];
    [ "--seed" ],
      Some "<integer>",
      [ "seed for the random number generator" ],
      TA.Default (string_of_int Defaults.seed |> Fun.const),
      (fun _ -> Parameters.seed := TA.get_parameter_int ());
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned" ],
      TA.Default (string_of_int Defaults.threads |> Fun.const),
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
      (fun _ -> Printf.printf "%s\n%!" info.TA.version; exit 0);
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 1)
  ];
  let verbose = !Parameters.verbose and threads = !Parameters.threads in
  let state = Random.State.make [| !Parameters.seed |] in
  let db = KMerDB.of_binary ~verbose !Parameters.input in
  (* ROUND ZERO NEEDS AXES AND HAS NO PARTITION TO TAKE THEM FROM, so it takes them from the
     spectra themselves.  A sample of them will do -- see --projection-sample -- and the
     database becomes that sample by removing everything the sample does not name *)
  let twister_of_sample () =
    let n = !Parameters.projection_sample in
    let sampled =
      if n <= 0 then db
      else
        KMerDB_Base.remove_selected db
          (complement_of_diverse_sample ~threads ~verbose state db n) in
    if verbose then
      Printf.eprintf "%s Building the initial projection...\n%!" prefix;
    let twister, _, _ = CA.twist ~threads ~verbose sampled in
    twister
  (* AND EVERY LATER ROUND TAKES THEM FROM THE PARTITION, combining the spectra of each class
     into the one that stands for it and running the analysis on those, so that the axes become
     the directions telling the classes apart rather than the directions a sample happened to
     span.  Every spectrum is named by [names], those being the rows of the register the whole
     database was projected into, so none is left without a class -- which matters, the
     combination putting everything unnamed into a single class of its own that would then sit
     in the middle of the space and win queries against the real ones *)
  (* DRIVEN OFF THE ASSIGNMENT AND NOT OFF THE NAMES.  A matrix's array of row names can be
     longer than the matrix has rows, being grown in steps and left with room to spare, so
     walking the names and indexing the assignment by the same counter runs off the end of the
     assignment -- and the two lengths agreeing in one round is no promise that they will in the
     next.  The assignment has exactly one entry per row, and every value in it is a row index,
     so both lookups below are in range by construction *)
  and twister_of_partition names assign =
    let acc = ref db in
    Array.iteri
      (fun i a ->
        KMerDB_Base.set_metadata acc names.(i) "CLASS" (Printf.sprintf "C@%s" names.(a)))
      assign;
    let combined =
      KMerDB.split_spectra ~threads ~verbose !acc "CLASS" !Parameters.combination_criterion in
    let twister, _, _ = CA.twist ~threads ~verbose combined in
    twister in
  let cluster twister =
    let twisted =
      Twister.add_twisted_from_database ~threads ~verbose twister Twisted.empty
        !Parameters.input in
    let mat = twisted.Twisted.twisted.Matrix.matrix
    and iv = twisted.Twisted.inertia.Matrix.matrix.Matrix.Base.data.(0) in
    (* ASKED BEFORE THE SEARCH RATHER THAN AFTER IT, because the search cannot answer it: some
       partition always maximises the criterion, so one always comes back, and over a set that
       does not fall into groups that partition is a property of the criterion.  What decides
       is whether the distances are two-moded at all *)
    let structured, _ =
      Clustering.assess_structure ~verbose ~seed:!Parameters.seed ~what_label:"samples"
        ~metric:!Parameters.metric ~distance:!Parameters.distance
        ~distance_normalize:!Parameters.distance_normalize mat.Matrix.Base.data iv in
    (* A SINGLE MODE HAS TWO READINGS AND THEY ARE DIFFERENT CLAIMS, so the message says which
       one applies rather than asserting the stronger.  An embedding built from a sample can
       miss groups that are really there: a sample drawn uniformly from a corpus whose classes
       differ greatly in size is mostly the big ones, and its axes then span their internal
       variation rather than what separates one class from another *)
    if not structured && not !Parameters.report_anyway then
      Exception.raise __FUNCTION__ IO_Format
        (if !Parameters.projection_sample > 0 then
           Printf.sprintf
             "the embedding built from a sample of %d spectra shows no groups: its pairwise \
              distances have a single mode and no valley between two.  The sample may simply \
              not span the corpus -- raise --projection-sample, or set it to 0 to use every \
              spectrum -- or pass --report-anyway to have a partition written out regardless"
             !Parameters.projection_sample
         else
           "the samples do not fall into groups: their pairwise distances have a single mode \
            and no valley between two, so any partition of them would be an artefact of the \
            criterion used to find it rather than a property of the data.  Pass --report-anyway \
            to have one written out regardless");
    let assign =
      Clustering.run_montecarlo ~verbose ~threads ~seed:!Parameters.seed
        ~replicas:!Parameters.montecarlo_replicas ~what_label:"samples"
        ~steps:!Parameters.montecarlo_steps ~sample:!Parameters.montecarlo_sample
        ~temperature:!Parameters.montecarlo_temperature ~decades:!Parameters.montecarlo_decades
        ~cooling:!Parameters.montecarlo_cooling ~metric:!Parameters.metric
        ~distance:!Parameters.distance ~distance_normalize:!Parameters.distance_normalize
        mat.Matrix.Base.data mat.Matrix.Base.row_names iv in
    twisted, mat.Matrix.Base.row_names, assign in
  (* THE LOOP, which is the program.  It stops when a round hands back what it was given,
     because that partition reproduces itself through a projection built from it and is
     therefore a property of the database rather than of the sample round zero began with *)
  let twisted = ref (Twisted.empty: Twisted.t)
  and twister = ref Twister.empty and names = ref [||] and assign = ref [||]
  and settled = ref false and round = ref 0 in
  while not !settled && !round < !Parameters.iterations do
    incr round;
    let tw = if !round = 1 then twister_of_sample () else twister_of_partition !names !assign in
    let tws, nms, asg = cluster tw in
    let n_classes = canonical asg |> Array.to_list |> IntSet.of_list |> IntSet.cardinal in
    settled := !round > 1 && canonical asg = canonical !assign;
    twister := tw; twisted := tws; names := nms; assign := asg;
    if verbose then
      Printf.eprintf "%s Round %d: %d %s%s.\n%!" prefix !round n_classes
        (String.pluralize_int ~plural:"classes" "class" n_classes)
        (if !settled then ", which is the one it was given -- settled" else "")
  done;
  if verbose && not !settled then
    Printf.eprintf "%s Stopped at the %d-round limit without settling.\n%!" prefix !round;
  Twister.to_binary ~verbose !twister !Parameters.output;
  Twisted.to_binary ~verbose !twisted !Parameters.output;
  (* The partition, in the two-line form KPopCountDB reads back *)
  let classes = !Parameters.output ^ ".KPopClasses.txt" in
  Exception.catch_unexpected_end_of_output __FUNCTION__
    (fun () ->
      let output = open_out classes in
      (* As many names as there are entries in the assignment, and no more: the array holding
         them can be longer than the register has rows, and a header naming samples the body
         then has no label for is a class file no reader can line up *)
      Array.iteri
        (fun i _ -> output_char output '\t'; output_string output (!names).(i)) !assign;
      output_char output '\n';
      output_string output "CLASS";
      Array.iter
        (fun a ->
          output_char output '\t';
          output_string output (Printf.sprintf "C@%s" (!names).(a)))
        !assign;
      output_char output '\n';
      close_out output);
  if verbose then
    Printf.eprintf "%s Partition written to '%s'.\n%!" prefix classes

