(*
    KPopTwistDB.ml -- (c) 2022-2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    KPopTwistDB is the binary for projecting new samples through an
    existing Twister and computing downstream artefacts: pairwise
    distances, summaries, nearest neighbours, embeddings and greedy
    leader clustering.  (Phylogenetic-tree construction now lives in
    the separate KPopPhylo binary.)  It implements the CLI surface
    (`Tools.Argv`-based options, `to_do_t` two-pass dry-run /
    execution model) for every operation on the Twister and Twisted
    database registers.

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

module RegisterType =
  struct
    type t =
      | Twister
      | Twisted
    let of_string = function
      | "T" -> Twister
      | "t" -> Twisted
      | s ->
        Exception.raise_unrecognized_initializer __FUNCTION__ "register type" s
  end

module KeepAtMost =
  struct
    type t = int option
    let of_string = function
      | "all" -> None
      | s ->
        try
          let res = int_of_string s in
          (* This raise is OK because it is caught by the surrounding 'with _' *)
          if res <= 0 then
            raise_notrace Exit;
          Some res
        with _ ->
          Exception.raise_unrecognized_initializer __FUNCTION__ "value for --keep-at-most option" s
    let to_string = function
      | None -> "all"
      | Some n -> string_of_int n
  end

type to_do_t =
  | Empty of RegisterType.t
  (* Parameter is input prefix *)
  | Binary_to_register of RegisterType.t * string
  (* Parameter is input prefix *)
  | Tables_to_register of RegisterType.t * string
  (* Parameter is input prefix *)
  | Add_binary_to_twisted of string
  (* Parameter is input prefix *)
  | Twist_database of string
  (* Parameter is output prefix *)
  | Register_to_binary of RegisterType.t * string
  | Set_precision_tables of int
  (* Parameter is output prefix *)
  | Register_to_tables of RegisterType.t * string
  | Set_metric of Space.Distance.Metric.t
  | Set_distance of Space.Distance.t
  | Set_distance_normalize of bool
  (* Computes embeddings from twisted vectors using the current metric/distance.
     Parameter is output prefix *)
  | Embeddings_from_twisted of string
  | Set_summary_keep_at_most of KeepAtMost.t
  (* Parameters are input probes, output prefix for summary/distance matrix,
      and a boolean specifying whether the distance matrix should be output *)
  | Summary_from_twisted_binary of string * string * bool
  | Set_neighbors_keep_at_most of KeepAtMost.t
  | Set_neighbors_guard_policy of Twisted.NeighborsPolicy.t
  | Set_neighbors_index_type of Interfaiss.Type.t
  (* Parameters are input probes and output summary *)
  | Summary_from_twisted_neighbors of string * string
  (* Clustering *)
  | Set_clusters_method of Clustering.Algorithm.t
  | Set_clusters_greedy_epsilon of Clustering.GreedyEpsilon.t
  | Set_clusters_greedy_order of Clustering.Order.t
  | Set_clusters_greedy_density_sample_number of int
  | Set_clusters_greedy_index_type of Interfaiss.Type.t
  | Set_clusters_hdbscan_min_cluster_size of int
  | Set_clusters_hdbscan_min_samples of int option
  | Set_clusters_hdbscan_mst_mode of Clustering.HdbscanMstMode.t
  | Set_clusters_hdbscan_num_neighbors of int option
  | Set_clusters_hdbscan_index_type of Interfaiss.Type.t
  | Clusters_kmers of string
  | Clusters_samples of string

module Defaults =
  struct
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalize = false
    let metric = Space.Distance.Metric.of_string "powers(1,1,1)"
    let precision_tables = 15
    let summary_keep_at_most = Some 2
    let neighbors_keep_at_most = Some 6
    let neighbors_guard_policy = Twisted.NeighborsPolicy.of_string "times(2)"
    let neighbors_index_type = Interfaiss.Type.of_string "hnsw(32)"
    (* Clustering *)
    let clusters_method = Clustering.Algorithm.Greedy
    let clusters_greedy_epsilon = Clustering.GreedyEpsilon.FirstNN
    let clusters_greedy_order = Clustering.Order.Inertia
    let clusters_greedy_density_sample_number = 200
    let clusters_greedy_index_type = Interfaiss.Type.of_string "hnsw(32)"
    let clusters_hdbscan_min_cluster_size = 1
    let clusters_hdbscan_min_samples = (None : int option)
    let clusters_hdbscan_mst_mode = Clustering.HdbscanMstMode.of_string "auto"
    let clusters_hdbscan_num_neighbors = (None : int option)
    let clusters_hdbscan_index_type = Interfaiss.Type.of_string "hnsw(32)"
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
  Tools.Argv.name = "KPopTwistDB";
  version = "48";
  date = "07-Apr-2026"
} and authors = [
  "2022-2026", "Paolo Ribeca", "paolo.ribeca@gmail.com";
  "2024     ", "Ünsal Öztürk", "uensal.oeztuerk@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  let prefix = Printf.sprintf "(%s):" __FUNCTION__ in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    TA.make_separator_multiline [ ""; "Actions on the database registers - Input/Output operations:" ];
    [ "-0"; "--zero"; "--empty" ],
      Some "'T'|'t'",
      [ "load an empty database into the specified register";
        " ('T'=twister; 't'=twisted)" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Empty register_type |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "'T'|'t' <binary_file_prefix>",
      [ "load the specified binary database into the specified register";
        " ('T'=twister; 't'=twisted).";
        "File extension is automatically determined depending on database type";
        " (will be '.KPopTwister' or '.KPopTwisted', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Binary_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "'T'|'t' <tabular_file_prefix>",
      [ "load the specified tabular database(s) into the specified register";
        " ('T'=twister; 't'=twisted).";
        "File extension is automatically determined depending on database type";
        " (will be: '.KPopTwister.txt' and '.KPopInertia.txt;';
                or: '.KPopInertia.txt' and '.KPopTwisted.txt', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Tables_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-a"; "--add"; "--add-to-twisted" ],
      Some "<binary_file_prefix>",
      [ "add the contents of the specified binary database to the twisted register.";
        "File extension is automatically determined";
        " (will be '.KPopTwisted', unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Add_binary_to_twisted (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "'T'|'t' <binary_file_prefix>",
      [ "save the database present in the specified register";
        " ('T'=twister; 't'=twisted)";
        "to the specified binary file.";
        "File extension is automatically assigned depending on database type";
        " (will be '.KPopTwister' or '.KPopTwisted', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_binary (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--precision-for-tables" ],
      Some "<positive_integer>",
      [ "set how many precision digits should be used when outputting numbers";
        "in tabular formats" ],
      TA.Default (string_of_int Defaults.precision_tables |> Fun.const),
      (fun _ -> Set_precision_tables (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-O"; "--Output" ],
      Some "'T'|'t' <tabular_file_prefix>",
      [ "save the database present in the specified register";
        " ('T'=twister; 't'=twisted)";
        "to the specified tabular files.";
        "File extensions are automatically assigned depending on database type";
        " (will be: '.KPopTwister.txt' and '.KPopInertia.txt';";
        "       or: '.KPopInertia.txt' and '.KPopTwisted.txt', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_tables (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Actions on the database registers - Other operations:" ];
    [ "-t"; "--twist"; "--twist-kmers"; "--twist-spectra" ],
      Some "<binary_file_prefix>",
      [ "twist the k-mer spectra contained in the specified binary database";
        "according to the transformation present in the twister register,";
        "and add the results to the database loaded in the twisted register" ],
      TA.Optional,
      (fun _ -> Twist_database (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-m"; "--metric"; "--metric-function" ],
      Some "'flat'|'powers('POWERS_PARAMETERS')'",
      [ "where POWERS_PARAMETERS := ";
        "        INTERNAL_POWER','FRACTIONAL_ACCUMULATIVE_THRESHOLD','EXTERNAL_POWER";
        "      INTERNAL_POWER := <non-negative_float>";
        "      FRACTIONAL_ACCUMULATIVE_THRESHOLD := <fractional_float>";
        "      EXTERNAL_POWER := <non-negative_float>";
        "Set the metric function to be used when computing distances.";
        "The 'power' transformation is computed as follows:";
        " (1) the inertia vector is raised to INTERNAL_POWER and normalized;";
        " (2) elements are summed in order until FRACTIONAL_ACCUMULATIVE_THRESHOLD";
        "     (a number between 0. and 1.) is reached, while the elements";
        "     above the threshold are set to zero";
        " (3) the resulting vector is raised to EXTERNAL_POWER and normalized.";
        "Note that";
        " 'flat'";
        "(which is equivalent to 'power(0,1,1)' or 'power(1,1,0)')";
        "disregards inertia, i.e. it is the same as standard coordinates, while";
        " 'power(1,1,1)'";
        "leaves inertia unchanged, i.e. it is the same as principal coordinates" ],
      TA.Default (Space.Distance.Metric.to_string Defaults.metric |> Fun.const),
      (fun _ ->
        Set_metric (TA.get_parameter () |> Space.Distance.Metric.of_string) |> List.accum Parameters.program);
    [ "--distance"; "--distance-function" ],
      Some "'euclidean'|'cosine'|'angle'|'minkowski('<non-negative_float>')'",
      [ "set the function to be used when computing distances.";
        "The parameter for 'minkowski' is the power.";
        "Note that:";
        " 'euclidean' is the same as 'minkowski(2)';";
        " 'cosine' is the same as ('euclidean'^2)/2, or 1 - cos(theta);";
        " 'angle' is the same as arccos(1 - ('euclidean'^2)/2), or theta,";
        "where theta is the relative angle between the two embeddings" ],
      TA.Default (Space.Distance.to_string Defaults.distance |> Fun.const),
      (fun _ -> Set_distance (TA.get_parameter () |> Space.Distance.of_string) |> List.accum Parameters.program);
    [ "--distance-normalize"; "--distance-normalization" ],
      Some "'true'|'false'",
      [ "whether to normalize twisted vectors before computing distances.";
        "It must be 'true' when the distance function is 'cosine' or 'angle'" ],
      TA.Default (string_of_bool Defaults.distance_normalize |> Fun.const),
      (fun _ -> Set_distance_normalize (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "-e"; "--embeddings"; "--compute-embeddings"; "--twisted-to-embeddings" ],
      Some "<tabular_file_prefix>",
      [ "compute embeddings from the vectors present in the twisted register";
        "using the current metric function, distance function and normalization.";
        "The result will be written to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopVectors.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Embeddings_from_twisted (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--distances-summarize-at-most"; "--distances-in-summary" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of closest sequences to be printed";
        "when summarizing distances.";
        "Note that more might be printed anyway in case of ties.";
        "The statistics in the summary will be computed on all sequences" ],
      TA.Default (KeepAtMost.to_string Defaults.summary_keep_at_most |> Fun.const),
      (fun _ ->
        Set_summary_keep_at_most (TA.get_parameter () |> KeepAtMost.of_string) |> List.accum Parameters.program);
    [ "-d"; "--summarize-distances"; "--compute-and-summarize-distances" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "for each vector present in the twisted register, compute distances";
        "to all vectors present in the specified twisted binary file";
        " (which must have extension '.KPopTwisted' unless file is '/dev/*')";
        "using the current metric function, distance function and normalization;";
        "summarize them and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopSummary.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_binary (twisted_prefix, TA.get_parameter (), false) |> List.accum Parameters.program);
    [ "-D"; "--summarize-and-output-distances"; "--compute-summarize-and-output-distances" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "same as option '-d', but additionally output the distance matrix";
        "in tabular form.";
        "File extensions are automatically assigned";
        " (will be '.KPopSummary.txt' and '.KPopDMatrix.txt',";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_binary (twisted_prefix, TA.get_parameter (), true) |> List.accum Parameters.program);
    [ "--neighbors-index-type"; "--neighbors-faiss-index-type" ],
      Some "'flat'|'pq('PQ_PARAMETERS')'|'hnsw('<positive_integer>')'",
      [ "where PQ_PARAMETERS :=";
        " <positive_integer>','<positive_integer>'";
        "Set the type of Faiss index used to compute nearest neighbors.";
        "Parameters for 'pq' are:";
        " number of subquantizers; bits per subquantizer.";
        "Note that the product of the two must be less than or equal to";
        "the number of dimensions of the twisted vectors.";
        "The parameter for 'hnsw' is";
        " hyperparameter M.";
        "Note that some indices may not be able to return all the existing";
        "nearest neighbors" ],
      TA.Default (Interfaiss.Type.to_string Defaults.neighbors_index_type |> Fun.const),
      (fun _ ->
        Set_neighbors_index_type (TA.get_parameter () |> Interfaiss.Type.of_string) |> List.accum Parameters.program);
    [ "--neighbors-summarize-at-most"; "--neighbors-in-summary" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of closest sequences to be printed";
        "when summarizing nearest neighbors.";
        "Note that more might be printed anyway in case of ties.";
        "The statistics in the summary will be computed on all the neighbors";
        "explored according to the policy specified by option";
        "'--neighbors-guard-policy'" ],
      TA.Default (KeepAtMost.to_string Defaults.neighbors_keep_at_most |> Fun.const),
      (fun _ ->
        Set_neighbors_keep_at_most (TA.get_parameter () |> KeepAtMost.of_string) |> List.accum Parameters.program);
    [ "--neighbors-guard-policy"; "--neighbors-exploration-policy" ],
      Some "'times('<float_no_less_than_one>')'|'plus(<non-negative_integer>)'",
      [ "set the number of nearest neighbors to be explored";
        "when summarizing them.";
        "Note that this is greater than or equal to the number of neighbors";
        "specified with option '--neighbors-summarize-at-most'.";
        "Calling the latter n,";
        " policy 'times('m')'";
        "will explore m*n nearest neighbors, while";
        " policy 'plus('m')'";
        "will explore m+n nearest neighbors.";
        "The additional neighbors explored are not printed,";
        "but used to compute overall statistics" ],
      TA.Default (Twisted.NeighborsPolicy.to_string Defaults.neighbors_guard_policy |> Fun.const),
      (fun _ ->
        Set_neighbors_guard_policy
          (TA.get_parameter () |> Twisted.NeighborsPolicy.of_string) |> List.accum Parameters.program);
    [ "-n"; "--summarize-neighbors"; "--find-and-summarize-neighbors" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "for each vector present in the twisted register, find nearest neighbors";
        "among the vectors present in the specified twisted binary file";
        " (which must have extension '.KPopTwisted' unless file is '/dev/*')";
        "using the current metric function, distance function and normalization;";
        "summarize distances and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopSummary.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_neighbors (twisted_prefix, TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Actions on the database registers - Clustering operations:" ];
    [ "--clusters-method" ],
      Some "'greedy'|'hdbscan'",
      [ "clustering algorithm.";
        "'greedy': greedy-leader clustering (see --clusters-greedy-* knobs).";
        "'hdbscan': HDBSCAN* density-based clustering (see --clusters-hdbscan-*";
        "  knobs).  Same metric/distance/normalisation pre-scaling as 'greedy'." ],
      TA.Default (Clustering.Algorithm.to_string Defaults.clusters_method |> Fun.const),
      (fun _ ->
        Set_clusters_method (TA.get_parameter () |> Clustering.Algorithm.of_string) |> List.accum Parameters.program);
    [ "--clusters-greedy-epsilon" ],
      Some "'firstNN'|'density'|'adaptive'",
      [ "epsilon-estimation strategy for the greedy-leader algorithm:";
        "'firstNN': kneedle elbow in sorted FAISS 1-NN distances";
        "  (O(n log n) with HNSW, O(n^2) with flat);";
        "'density': kneedle elbow in sorted dist_star values, where dist_star";
        "  is the distance maximising k/V(d_k,D) for each point";
        "  (O(n_sample * n), or O(n^2) when --clusters-greedy-order density is also set;";
        "  use --clusters-greedy-order firstNN for O(n log n) ordering instead);";
        "'adaptive': per-point dist_star as the absorption threshold (no global";
        "  epsilon).  Computes dist_star for all n points (O(n^2)) and forces the";
        "  processing order to ascending dist_star, so --clusters-greedy-order";
        "  is ignored in this mode.  Handles multi-scale cluster structure that";
        "  a single global threshold cannot capture.";
        "Ignored unless --clusters-method 'greedy' is in effect." ],
      TA.Default (Clustering.GreedyEpsilon.to_string Defaults.clusters_greedy_epsilon |> Fun.const),
      (fun _ ->
        Set_clusters_greedy_epsilon (TA.get_parameter () |> Clustering.GreedyEpsilon.of_string) |> List.accum Parameters.program);
    [ "--clusters-greedy-order" ],
      Some "'inertia'|'firstNN'|'density'",
      [ "order in which points are processed by the greedy-leader clusterer:";
        "'inertia': decreasing row inertia proxy sum_d(lambda_d * T[i,d]^2),";
        "  so the most informative k-mers / most distinctive samples become";
        "  cluster representatives;";
        "'firstNN': increasing FAISS 1-NN distance (densest regions first, O(n log n));";
        "  when --clusters-greedy-epsilon firstNN is also set, the FAISS distances computed";
        "  for epsilon estimation are reused for ordering at no extra cost;";
        "'density': increasing dist_star (densest regions first, O(n^2));";
        "  when --clusters-greedy-epsilon density is also set, the dist_star values";
        "  computed for epsilon estimation are reused for ordering.";
        "Ignored unless --clusters-method 'greedy' is in effect, or when";
        "--clusters-greedy-epsilon 'adaptive' (which forces ascending-dist_star order)." ],
      TA.Default (Clustering.Order.to_string Defaults.clusters_greedy_order |> Fun.const),
      (fun _ ->
        Set_clusters_greedy_order (TA.get_parameter () |> Clustering.Order.of_string) |> List.accum Parameters.program);
    [ "--clusters-greedy-density-sample-number" ],
      Some "<positive_integer>",
      [ "number of points randomly sampled for dist_star estimation";
        "when --clusters-greedy-epsilon density and --clusters-greedy-order inertia or firstNN are set.";
        "When --clusters-greedy-order density is also set, all n points are used.";
        "Ignored unless --clusters-method 'greedy' is in effect." ],
      TA.Default (string_of_int Defaults.clusters_greedy_density_sample_number |> Fun.const),
      (fun _ ->
        Set_clusters_greedy_density_sample_number (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    [ "--clusters-greedy-index-type" ],
      Some "'flat'|'hnsw('<positive_integer>')'",
      [ "FAISS index type used for 1-NN estimation and greedy-leader clustering.";
        "Ignored unless --clusters-method 'greedy' is in effect." ],
      TA.Default (Interfaiss.Type.to_string Defaults.clusters_greedy_index_type |> Fun.const),
      (fun _ ->
        Set_clusters_greedy_index_type (TA.get_parameter () |> Interfaiss.Type.of_string) |> List.accum Parameters.program);
    [ "--clusters-hdbscan-min-cluster-size" ],
      Some "<positive_integer>",
      [ "minimum cluster size for the 'hdbscan' clustering algorithm.";
        "Same semantics as KPopPhylo's --hdbscan-min-cluster-size, but settable";
        "independently for the clusters consumer of the HDBSCAN core.";
        "Ignored unless --clusters-method 'hdbscan' is in effect." ],
      TA.Default (string_of_int Defaults.clusters_hdbscan_min_cluster_size |> Fun.const),
      (fun _ ->
        Set_clusters_hdbscan_min_cluster_size (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    [ "--clusters-hdbscan-min-samples" ],
      Some "<positive_integer>",
      [ "k for the core-distance neighbourhood of the 'hdbscan' clustering";
        "algorithm.  When unset (default), k is taken equal to";
        "--clusters-hdbscan-min-cluster-size, matching the reference HDBSCAN";
        "one-knob ergonomic.";
        "Ignored unless --clusters-method 'hdbscan' is in effect." ],
      TA.Default (Fun.const "same as --clusters-hdbscan-min-cluster-size"),
      (fun _ ->
        Set_clusters_hdbscan_min_samples (Some (TA.get_parameter_int_pos ()))
        |> List.accum Parameters.program);
    [ "--clusters-hdbscan-mst-mode" ],
      Some "'auto'|'sparse'|'dense'",
      [ "minimum-spanning-tree construction strategy for the 'hdbscan'";
        "clustering algorithm.  Same semantics as KPopPhylo's --hdbscan-mst-mode";
        "but settable independently.";
        "Ignored unless --clusters-method 'hdbscan' is in effect." ],
      TA.Default (Clustering.HdbscanMstMode.to_string Defaults.clusters_hdbscan_mst_mode |> Fun.const),
      (fun _ ->
        Set_clusters_hdbscan_mst_mode (TA.get_parameter () |> Clustering.HdbscanMstMode.of_string)
        |> List.accum Parameters.program);
    [ "--clusters-hdbscan-num-neighbors" ],
      Some "<positive_integer>",
      [ "number of nearest neighbours per point used to build the FAISS k-NN";
        "candidate graph for the sparse 'hdbscan' MST.  When unset, auto-computed";
        "as max(--clusters-hdbscan-min-samples + 1, min(n - 1, 30)).";
        "Ignored unless --clusters-method 'hdbscan' and";
        "--clusters-hdbscan-mst-mode 'sparse' are in effect." ],
      TA.Default (Fun.const "auto (max(min_samples + 1, min(n - 1, 30)))"),
      (fun _ ->
        Set_clusters_hdbscan_num_neighbors (Some (TA.get_parameter_int_pos ()))
        |> List.accum Parameters.program);
    [ "--clusters-hdbscan-index-type" ],
      Some "'flat'|'pq('PQ_PARAMETERS')'|'hnsw('<positive_integer>')'",
      [ "FAISS index type used by the sparse 'hdbscan' MST when used as the";
        "clustering algorithm.  Same syntax as --neighbors-index-type.";
        "Ignored unless --clusters-method 'hdbscan' and";
        "--clusters-hdbscan-mst-mode 'sparse' are in effect." ],
      TA.Default (Interfaiss.Type.to_string Defaults.clusters_hdbscan_index_type |> Fun.const),
      (fun _ ->
        Set_clusters_hdbscan_index_type (TA.get_parameter () |> Interfaiss.Type.of_string)
        |> List.accum Parameters.program);
    [ "-c"; "--clusters" ],
      Some "'T' <kmer_list_file>|'t' <class_file>",
      [ "apply clustering to the contents of the specified register";
        " ('T'=twister, clusters k-mers; 't'=twisted, clusters samples).";
        "The algorithm is selected by --clusters-method (default 'greedy';";
        "  see also 'hdbscan').";
        "Uses the current metric, distance, and normalization settings.";
        "The cluster assignment table is written to stdout.";
        "For 'T': k-mer standard coordinates are recovered from the twister";
        "  as km_std[i][d] = Twister[d,i] * sqrt(inertia[d]);";
        "  the names of representative k-mers are written to <kmer_list_file>,";
        "  one per line and with no header, ready to be passed to KPopTwist --keep.";
        "  Greedy writes the cluster representatives; HDBSCAN writes one";
        "  representative k-mer per cluster plus every noise k-mer.";
        "For 't': a two-line tab-separated class file is written to <class_file>";
        "  (header line of sample names; 'CLASS' line of class labels),";
        "  ready to be passed to KPopCountDB -m <class_file> -c CLASS.";
        "  Greedy labels each sample as 'C@<representative_name>';";
        "  HDBSCAN labels assigned samples as 'C@<integer>' (cluster id) and";
        "  outlier samples as 'noise'." ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister -> Clusters_kmers (TA.get_parameter ()) |> List.accum Parameters.program
        | Twisted -> Clusters_samples (TA.get_parameter ()) |> List.accum Parameters.program);
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
    (* Hidden option to emit help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    (* Hidden option to print exception backtrace *)
    [ "-x"; "--print-exception-backtrace" ], None, [], TA.Optional, (fun _ -> Printexc.record_backtrace true);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 0)
  ];
  let program = List.rev !Parameters.program in
  if program = [] then begin
    TA.usage ();
    exit 0
  end;
  if !Parameters.verbose then
    TA.header ();
  CPU.warn_if_slow ~verbose:!Parameters.verbose ();
  Random.self_init ();
  (* We perform a dry run of the program to detect possible errors *)
  let twister_loaded = ref false
  and distance = ref Defaults.distance and distance_normalize = ref Defaults.distance_normalize in
  List.iter
    (function
      | Empty _ ->
        ()
      | Binary_to_register (Twister, _) | Tables_to_register (Twister, _) ->
        twister_loaded := true
      | Binary_to_register (Twisted, _) | Tables_to_register (Twisted, _)
      | Add_binary_to_twisted _ ->
        ()
      | Twist_database _ ->
        (* A twister must have been loaded to twist spectra *)
        if not !twister_loaded then
          TA.parse_error
            "Option '-k'/'-s' requires a twister in the twister register!"
      | Register_to_binary _
      | Set_precision_tables _ | Register_to_tables _
      | Set_metric _ ->
        ()
      | Set_distance dist ->
        distance := dist
      | Set_distance_normalize norm ->
        distance_normalize := norm
      | Embeddings_from_twisted _
      | Summary_from_twisted_binary _ | Summary_from_twisted_neighbors _ ->
        begin match !distance, !distance_normalize with
        | Space.Distance.Cosine, false | Angle, false ->
          TA.parse_error "Distances 'cosine' and 'angle' require embeddings to be normalized"
        | Cosine, true | Angle, true | Euclidean, _ | Minkowski _, _ ->
          ()
        end
      | Set_summary_keep_at_most _
      | Set_neighbors_keep_at_most _ | Set_neighbors_guard_policy _ | Set_neighbors_index_type _
      | Set_clusters_method _
      | Set_clusters_greedy_epsilon _ | Set_clusters_greedy_order _ | Set_clusters_greedy_density_sample_number _
      | Set_clusters_greedy_index_type _
      | Set_clusters_hdbscan_min_cluster_size _ | Set_clusters_hdbscan_min_samples _
      | Set_clusters_hdbscan_mst_mode _ | Set_clusters_hdbscan_num_neighbors _ | Set_clusters_hdbscan_index_type _ ->
        ()
      | Clusters_kmers _ ->
        if not !twister_loaded then
          TA.parse_error
            "Option '--clusters T' requires a twister to have been loaded first!";
        (match !distance, !distance_normalize with
        | Space.Distance.Cosine, false | Space.Distance.Angle, false ->
          TA.parse_error
            "Distances 'cosine' and 'angle' require embeddings to be normalized \
             (add --distance-normalize true before --clusters)"
        | _ -> ())
      | Clusters_samples _ ->
        (match !distance, !distance_normalize with
        | Space.Distance.Cosine, false | Space.Distance.Angle, false ->
          TA.parse_error
            "Distances 'cosine' and 'angle' require embeddings to be normalized \
             (add --distance-normalize true before --clusters)"
        | _ -> ()))
    program;
  (* These are the registers available to the program *)
  let twister = ref Twister.empty and twisted = ref Twisted.empty and metric = ref Defaults.metric
  and distance = ref Defaults.distance and distance_normalize = ref Defaults.distance_normalize
  and summary_keep_at_most = ref Defaults.summary_keep_at_most
  and neighbors_keep_at_most = ref Defaults.neighbors_keep_at_most
  and neighbors_guard_policy = ref Defaults.neighbors_guard_policy
  and neighbors_index_type = ref Defaults.neighbors_index_type
  and precision_tables = ref Defaults.precision_tables
  and clusters_method = ref Defaults.clusters_method
  and clusters_greedy_epsilon = ref Defaults.clusters_greedy_epsilon and clusters_greedy_order = ref Defaults.clusters_greedy_order
  and clusters_greedy_density_sample_number = ref Defaults.clusters_greedy_density_sample_number
  and clusters_greedy_index_type = ref Defaults.clusters_greedy_index_type
  and clusters_hdbscan_min_cluster_size = ref Defaults.clusters_hdbscan_min_cluster_size
  and clusters_hdbscan_min_samples = ref Defaults.clusters_hdbscan_min_samples
  and clusters_hdbscan_mst_mode = ref Defaults.clusters_hdbscan_mst_mode
  and clusters_hdbscan_num_neighbors = ref Defaults.clusters_hdbscan_num_neighbors
  and clusters_hdbscan_index_type = ref Defaults.clusters_hdbscan_index_type in
  let twisted_of_binary = Twisted.of_binary ~verbose:!Parameters.verbose
  and twisted_of_files = Twisted.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose in
  try
    List.iter
      (function
        | Empty Twister ->
          twister := Twister.empty
        | Empty Twisted ->
          twisted := Twisted.empty
        | Binary_to_register (Twister, prefix) ->
          twister := Twister.of_binary ~verbose:!Parameters.verbose prefix
        | Binary_to_register (Twisted, prefix) ->
          twisted := twisted_of_binary prefix
        | Tables_to_register (Twister, prefix) ->
          twister := Twister.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose prefix
        | Tables_to_register (Twisted, prefix) ->
          twisted := twisted_of_files prefix
        | Add_binary_to_twisted prefix ->
          twisted := twisted_of_binary prefix |> Twisted.merge_rowwise !twisted
        | Twist_database prefix ->
          twisted :=
            Twister.add_twisted_from_database
              ~threads:!Parameters.threads ~verbose:!Parameters.verbose !twister !twisted prefix
        | Register_to_binary (Twister, prefix) ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () -> Twister.to_binary ~verbose:!Parameters.verbose !twister prefix)
        | Register_to_binary (Twisted, prefix) ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () -> Twisted.to_binary ~verbose:!Parameters.verbose !twisted prefix)
        | Set_precision_tables prec ->
          precision_tables := prec
        | Register_to_tables (Twister, prefix) ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              Twister.to_files
                ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
                !twister prefix)
        | Register_to_tables (Twisted, prefix) ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              Twisted.to_files
                ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
                !twisted prefix)
        | Set_metric metr ->
          metric := metr
        | Set_distance dist ->
          distance := dist
        | Set_distance_normalize norm ->
          distance_normalize := norm
        | Embeddings_from_twisted prefix ->
          let res =
            Twisted.to_embeddings
              ~normalize:!distance_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !distance !metric !twisted in
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              Matrix.to_file
                ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
                res prefix)
        | Set_summary_keep_at_most kam ->
          summary_keep_at_most := kam
        | Summary_from_twisted_binary (prefix_in, prefix_out, output_distance_matrix) ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              Twisted.summarize_distances_rowwise
                ~normalize:!distance_normalize ~keep_at_most:!summary_keep_at_most ~output_distance_matrix
                ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
                !distance !metric (twisted_of_binary prefix_in) !twisted prefix_out)
        | Set_neighbors_keep_at_most nhm ->
          neighbors_keep_at_most := nhm
        | Set_neighbors_guard_policy gp ->
          neighbors_guard_policy := gp
        | Set_neighbors_index_type it ->
          neighbors_index_type := it
        | Summary_from_twisted_neighbors (prefix_in, prefix_out) ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              Twisted.summarize_neighbors
                ~normalize:!distance_normalize ~how_many:!neighbors_keep_at_most
                ~policy:!neighbors_guard_policy ~index_type:!neighbors_index_type
                ~threads:!Parameters.threads ~verbose:!Parameters.verbose
                !metric (twisted_of_binary prefix_in) !twisted prefix_out)
        (* Clustering setters *)
        | Set_clusters_method m ->
          clusters_method := m
        | Set_clusters_greedy_epsilon m ->
          clusters_greedy_epsilon := m
        | Set_clusters_greedy_order o ->
          clusters_greedy_order := o
        | Set_clusters_greedy_density_sample_number n ->
          clusters_greedy_density_sample_number := n
        | Set_clusters_greedy_index_type it ->
          clusters_greedy_index_type := it
        | Set_clusters_hdbscan_min_cluster_size n ->
          clusters_hdbscan_min_cluster_size := n
        | Set_clusters_hdbscan_min_samples k ->
          clusters_hdbscan_min_samples := k
        | Set_clusters_hdbscan_mst_mode mode ->
          clusters_hdbscan_mst_mode := mode
        | Set_clusters_hdbscan_num_neighbors k ->
          clusters_hdbscan_num_neighbors := k
        | Set_clusters_hdbscan_index_type idx ->
          clusters_hdbscan_index_type := idx
        (* Clustering actions *)
        | Clusters_kmers output_file ->
          (* Recover standard k-mer coordinates from the Twister register:
             km_std[i][dim] = Twister[dim,i] * sqrt(inertia[dim]) *)
          let tw = !twister in
          let twm = tw.Twister.twister.Matrix.matrix in
          let iv = tw.Twister.inertia.Matrix.matrix.Matrix.Base.data.(0) in
          let mm = Array.length twm.Matrix.Base.col_names in
          if mm = 0 then
            Exception.raise __FUNCTION__ IO_Format
              "twister register is empty, nothing to cluster";
          let kk = Array.length twm.Matrix.Base.row_names in
          let ( .@!() ) = Float.Array.( .@!() ) in
          let kc = Array.init mm (fun i ->
            Float.Array.init kk (fun dim ->
              twm.Matrix.Base.data.(dim).@!(i) *. sqrt iv.@!(dim))) in
          (* Write representative k-mer names to output file, one per line, no header.
             Greedy writes the cluster representatives; HDBSCAN writes one k-mer per
             cluster plus every noise point (each noise k-mer is unique) *)
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              let output = open_out output_file in
              (match !clusters_method with
              | Clustering.Algorithm.Greedy ->
                let rep_orig =
                  Clustering.run_greedy
                    ~verbose:!Parameters.verbose
                    ~what_label:"k-mers"
                    ~epsilon_:!clusters_greedy_epsilon
                    ~order_:!clusters_greedy_order
                    ~density_sample_number:!clusters_greedy_density_sample_number
                    ~index_type:!clusters_greedy_index_type
                    ~metric:!metric
                    ~distance:!distance
                    ~distance_normalize:!distance_normalize
                    kc twm.Matrix.Base.col_names iv in
                Array.iteri
                  (fun i ri ->
                    if ri = i then
                      Printf.fprintf output "%s\n" twm.Matrix.Base.col_names.(i))
                  rep_orig
              | Clustering.Algorithm.Hdbscan ->
                let cluster_of =
                  Clustering.run_hdbscan
                    ~verbose:!Parameters.verbose ~threads:!Parameters.threads
                    ~what_label:"k-mers"
                    ~min_cluster_size:!clusters_hdbscan_min_cluster_size
                    ?min_samples:!clusters_hdbscan_min_samples
                    ~mst_mode:!clusters_hdbscan_mst_mode
                    ?num_neighbors:!clusters_hdbscan_num_neighbors
                    ~index_type:!clusters_hdbscan_index_type
                    ~metric:!metric
                    ~distance:!distance
                    ~distance_normalize:!distance_normalize
                    kc twm.Matrix.Base.col_names iv in
                let seen_clusters = Hashtbl.create 16 in
                Array.iteri
                  (fun i cid ->
                    if cid < 0 || not (Hashtbl.mem seen_clusters cid) then begin
                      if cid >= 0 then Hashtbl.add seen_clusters cid ();
                      Printf.fprintf output "%s\n" twm.Matrix.Base.col_names.(i)
                    end)
                  cluster_of);
              close_out output);
          if !Parameters.verbose then
            Printf.eprintf "%s Representative k-mer list written to '%s'.\n%!"
              prefix output_file
        | Clusters_samples output_file ->
          let mat = (!twisted).Twisted.twisted.Matrix.matrix in
          let iv = (!twisted).Twisted.inertia.Matrix.matrix.Matrix.Base.data.(0) in
          if Array.length mat.Matrix.Base.row_names = 0 then
            Exception.raise __FUNCTION__ IO_Format
              "twisted register is empty, nothing to cluster";
          let n = Array.length mat.Matrix.Base.row_names in
          (* Write a two-line tab-separated class file readable by KPopCountDB -m/-c:
             line 1: sample names (header); line 2: CLASS followed by per-sample labels.
             Greedy labels each sample with its representative's name; HDBSCAN labels with
             a cluster id (C@<integer>) or 'noise' for outlier samples *)
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () ->
              let output = open_out output_file in
              for i = 0 to n - 1 do
                output_char output '\t';
                output_string output mat.Matrix.Base.row_names.(i)
              done;
              output_char output '\n';
              output_string output "CLASS";
              (match !clusters_method with
              | Clustering.Algorithm.Greedy ->
                let rep_orig =
                  Clustering.run_greedy
                    ~verbose:!Parameters.verbose
                    ~what_label:"samples"
                    ~epsilon_:!clusters_greedy_epsilon
                    ~order_:!clusters_greedy_order
                    ~density_sample_number:!clusters_greedy_density_sample_number
                    ~index_type:!clusters_greedy_index_type
                    ~metric:!metric
                    ~distance:!distance
                    ~distance_normalize:!distance_normalize
                    mat.Matrix.Base.data mat.Matrix.Base.row_names iv in
                for i = 0 to n - 1 do
                  output_char output '\t';
                  output_string output ("C@" ^ mat.Matrix.Base.row_names.(rep_orig.(i)))
                done
              | Clustering.Algorithm.Hdbscan ->
                let cluster_of =
                  Clustering.run_hdbscan
                    ~verbose:!Parameters.verbose ~threads:!Parameters.threads
                    ~what_label:"samples"
                    ~min_cluster_size:!clusters_hdbscan_min_cluster_size
                    ?min_samples:!clusters_hdbscan_min_samples
                    ~mst_mode:!clusters_hdbscan_mst_mode
                    ?num_neighbors:!clusters_hdbscan_num_neighbors
                    ~index_type:!clusters_hdbscan_index_type
                    ~metric:!metric
                    ~distance:!distance
                    ~distance_normalize:!distance_normalize
                    mat.Matrix.Base.data mat.Matrix.Base.row_names iv in
                for i = 0 to n - 1 do
                  output_char output '\t';
                  if cluster_of.(i) >= 0 then
                    output_string output ("C@" ^ string_of_int cluster_of.(i))
                  else
                    output_string output "noise"
                done);
              output_char output '\n';
              close_out output);
          if !Parameters.verbose then
            Printf.eprintf "%s Sample class file written to '%s'.\n%!"
              prefix output_file)
      program
  with e ->
    Exception.handle __FUNCTION__ TA.usage (fun () ->
      Printf.peprintf "%s This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!" prefix;
      Printf.peprintf "%s You might also wish to rerun me with option -x to get a full backtrace.\n%!" prefix
    ) e

