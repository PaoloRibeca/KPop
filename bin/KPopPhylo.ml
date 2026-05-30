(*
    KPopPhylo.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    KPopPhylo builds phylogenetic trees from KPop data.  It keeps three
    registers -- a counts register (k-mer spectra), a twister, and a
    twisted register -- and selects what to do from which of them are
    populated:

      * counts only            : twist (CA) -> topology -> refit
      * counts + twisted        : (no twist) topology from twisted -> refit
      * twisted only            : topology only
      * counts + input tree     : refit the given tree only

    Topology is built by sparse-NJ (or the other splits methods); when a
    counts register is present the branch lengths are, by default,
    refitted by sparse OLS against mash-like Jaccard distances read
    directly from the k-mer rows (no MinHash sketching), on the K-NN
    edges of the output tree.  The chosen workflow is announced on
    stderr.

    The CLI mirrors KPopTwistDB's structure (Tools.Argv options,
    two-pass dry-run / execution model over a to_do_t list); the
    tree-construction options are those of KPopTwistDB with the
    'phylo-' prefix stripped.

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
      | Spectra
      | Twister
      | Twisted
    let of_string = function
      | "s" -> Spectra
      | "T" -> Twister
      | "t" -> Twisted
      | s ->
        Exception.raise_unrecognized_initializer __FUNCTION__ "register type" s
  end

(* Branch-length refit source. *)
module RefitWith =
  struct
    type t =
      | NoRefit
      | Mash of Jaccard.Kind.t
    let of_string = function
      | "none" -> NoRefit
      | "jaccard" | "mash" | "binary" -> Mash Jaccard.Kind.Binary
      | "jaccard-weighted" | "mash-weighted" | "weighted" -> Mash Jaccard.Kind.Weighted
      | s ->
        Exception.raise_unrecognized_initializer __FUNCTION__ "refit source" s
    let to_string = function
      | NoRefit -> "none"
      | Mash Jaccard.Kind.Binary -> "jaccard"
      | Mash Jaccard.Kind.Weighted -> "jaccard-weighted"
  end

let info = {
  Tools.Argv.name = "KPopPhylo";
  version = "1";
  date = "30-May-2026"
} and authors = [
  "2026", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

type to_do_t =
  | Binary_to_register of RegisterType.t * string
  | Tree_to_register of string
  | Set_method of Twisted.Phylo.Method.t
  | Set_splits_keep_at_most of int
  | Set_metric of Space.Distance.Metric.t
  | Set_distance of Space.Distance.t
  | Set_distance_normalize of bool
  | Set_centroids_balance_penalty of Twisted.BalancePenalty.t
  | Set_centroids_num_seeds of int
  | Set_centroids_seed of int
  | Set_gaps_prefilter_kneedle of bool
  | Set_hdbscan_min_cluster_size of int
  | Set_hdbscan_min_samples of int
  | Set_hdbscan_mst_mode of Clustering.HdbscanMstMode.t
  | Set_hdbscan_num_neighbors of int
  | Set_hdbscan_index_type of Interfaiss.Type.t
  | Set_hdbscan_lengths_mode of Clustering.Hdbscan.LengthsMode.t
  | Set_snj_index_type of Interfaiss.Type.t
  | Set_snj_mode of SparseNJ.Mode.t
  | Set_snj_num_neighbors of int
  | Set_snj_k_query_factor of int
  | Set_snj_hyp_scale of float
  | Set_snj_distance of SparseNJ.Distance.t
  | Set_snj_row_sum of SparseNJ.RowSum.t
  | Set_snj_symmetry of SparseNJ.Symmetry.t
  | Set_refit_with of RefitWith.t
  | Set_refit_strategy of Refit.Strategy.t
  | Set_refit_num_neighbors of int
  | Set_refit_k of int
  | Compute_tree of string (* Output prefix *)

module Defaults =
  struct
    let threads = Processes.Parallel.get_nproc ()
    let verbose = false
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalize = false
    let metric = Space.Distance.Metric.of_string "powers(1,1,1)"
    let method_ = Twisted.Phylo.Method.of_string "sparse-nj"
    let splits_keep_at_most = 10000
    let centroids_balance_penalty = Twisted.BalancePenalty.default
    let centroids_num_seeds = 10
    let centroids_seed = (None : int option)
    let gaps_prefilter_kneedle = true
    let hdbscan_min_cluster_size = 1
    let hdbscan_min_samples = (None : int option)
    let hdbscan_mst_mode = Clustering.HdbscanMstMode.of_string "auto"
    let hdbscan_num_neighbors = (None : int option)
    let hdbscan_index_type = Interfaiss.Type.of_string "hnsw(32)"
    let hdbscan_lengths_mode = Clustering.Hdbscan.LengthsMode.of_string "mreach"
    let snj_index_type = Interfaiss.Type.of_string "hnsw(32)"
    let snj_mode = SparseNJ.Mode.of_string "dense"
    let snj_num_neighbors = 10
    let snj_k_query_factor = 3
    let snj_hyp_scale = 1.0
    let snj_distance = SparseNJ.Distance.of_string "saitou-nei"
    let snj_row_sum = SparseNJ.RowSum.of_string "knn"
    let snj_symmetry = SparseNJ.Symmetry.of_string "one"
    (* Branch-length refit: when a counts register is present, refit by
       default; the user can disable with --refit-with none. *)
    let refit_with = RefitWith.Mash Jaccard.Kind.Binary
    let refit_strategy = Refit.Strategy.of_string "reuse-knn"
    let refit_num_neighbors = 10
    let refit_k = 12
  end

module Parameters =
  struct
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
    let program = ref []
  end

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[-i|--input s|T|t <binary_prefix>]* [--input-tree <newick>] [OPTIONS] -o|--output <output_prefix>";
  TA.parse [
    TA.make_separator_multiline [ "Registers."; "";
      "KPopPhylo holds a counts register ('s', k-mer spectra), a twister";
      "register ('T') and a twisted register ('t').  The workflow run is";
      "decided from which registers are populated:";
      "  counts only          -> twist, build topology, refit";
      "  counts + twisted     -> build topology from twisted, refit";
      "  twisted only         -> build topology only";
      "  counts + input tree  -> refit the given tree only";
      "(refit happens whenever counts are present, unless --refit-with none)." ];
    [ "-i"; "--input" ],
      Some "s|T|t <binary_file_prefix>",
      [ "load into the specified register the binary KPop database with";
        "the specified prefix ('s'=counts/.KPopSpectra (in fact a count";
        "database), 'T'=.KPopTwister, 't'=.KPopTwisted)" ],
      TA.Mandatory,
      (fun _ ->
        let rt = TA.get_parameter () |> RegisterType.of_string in
        Binary_to_register (rt, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--input-tree"; "--tree" ],
      Some "<newick_file>",
      [ "load a topology in Newick format whose branch lengths are to be";
        "refitted (refit-only workflow; requires a counts register)" ],
      TA.Optional,
      (fun _ -> Tree_to_register (TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Embedding distance (used to build the topology):" ];
    [ "-m"; "--metric" ],
      Some "'flat'|'powers('<p_int>','<p_med>','<p_ext>')'|...",
      [ "metric (inertia weighting) applied to twisted vectors before";
        "computing distances for the topology" ],
      TA.Default (Space.Distance.Metric.to_string Defaults.metric |> Fun.const),
      (fun _ -> Set_metric (TA.get_parameter () |> Space.Distance.Metric.of_string)
        |> List.accum Parameters.program);
    [ "--distance"; "--distance-function" ],
      Some "'euclidean'|'cosine'|'angle'|'minkowski('<p>')'",
      [ "distance function used to build the topology" ],
      TA.Default (Space.Distance.to_string Defaults.distance |> Fun.const),
      (fun _ -> Set_distance (TA.get_parameter () |> Space.Distance.of_string)
        |> List.accum Parameters.program);
    [ "--distance-normalize"; "--distance-normalization" ],
      Some "'true'|'false'",
      [ "whether to normalize twisted vectors before computing distances";
        "(must be 'true' for 'cosine'/'angle')" ],
      TA.Default (string_of_bool Defaults.distance_normalize |> Fun.const),
      (fun _ -> Set_distance_normalize (TA.get_parameter () |> bool_of_string)
        |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Topology construction:" ];
    [ "--method" ],
      Some "'sparse-nj'|'gaps'|'centroids'|'hdbscan'",
      [ "topology-construction method (default 'sparse-nj')" ],
      TA.Default (Twisted.Phylo.Method.to_string Defaults.method_ |> Fun.const),
      (fun _ -> Set_method (TA.get_parameter () |> Twisted.Phylo.Method.of_string)
        |> List.accum Parameters.program);
    [ "--splits-keep-at-most" ],
      Some "<positive_integer>",
      [ "for splits-based methods, the maximum number of splits to keep" ],
      TA.Default (string_of_int Defaults.splits_keep_at_most |> Fun.const),
      (fun _ -> Set_splits_keep_at_most (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    [ "--snj-mode" ],
      Some "'dense'|'subquadratic'|'hyperbolic'|'cover-tree'|'periodic-rebuild'|'rp-forest'",
      [ "sparse-NJ implementation mode (see KPop-PhyloSplits-Subquadratic).";
        "'dense' (default, validated O(n^3)); 'periodic-rebuild' is the";
        "  recommended exact quadratic engine; 'rp-forest' with";
        "  '--snj-distance hyperbolic' is the sub-quadratic-memory option." ],
      TA.Default (SparseNJ.Mode.to_string Defaults.snj_mode |> Fun.const),
      (fun _ -> Set_snj_mode (TA.get_parameter () |> SparseNJ.Mode.of_string)
        |> List.accum Parameters.program);
    [ "--snj-distance" ],
      Some "'saitou-nei'|'centroid'|'hyperbolic'|'hyperbolic-frechet'|'hybrid'",
      [ "distance model fed to the sparse-NJ Q-formula (only honoured by";
        "'--snj-mode rp-forest')" ],
      TA.Default (SparseNJ.Distance.to_string Defaults.snj_distance |> Fun.const),
      (fun _ -> Set_snj_distance (TA.get_parameter () |> SparseNJ.Distance.of_string)
        |> List.accum Parameters.program);
    [ "--snj-num-neighbors" ],
      Some "<positive_integer>",
      [ "number K of nearest active neighbours kept per cluster" ],
      TA.Default (string_of_int Defaults.snj_num_neighbors |> Fun.const),
      (fun _ -> Set_snj_num_neighbors (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    [ "--snj-k-query-factor" ],
      Some "<positive_integer>",
      [ "FAISS/forest expansion factor (K_QUERY = K * factor)" ],
      TA.Default (string_of_int Defaults.snj_k_query_factor |> Fun.const),
      (fun _ -> Set_snj_k_query_factor (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    [ "--snj-hyp-scale" ],
      Some "<positive_float>",
      [ "radial lift scale for hyperbolic modes (best ~0.08 for the";
        "rp-forest hyperbolic distance)" ],
      TA.Default (string_of_float Defaults.snj_hyp_scale |> Fun.const),
      (fun _ -> Set_snj_hyp_scale (TA.get_parameter_float_pos ())
        |> List.accum Parameters.program);
    [ "--snj-rowsum" ],
      Some "'knn'|'topk'|'full'",
      [ "row-sum estimator for the sparse-NJ Q-formula" ],
      TA.Default (SparseNJ.RowSum.to_string Defaults.snj_row_sum |> Fun.const),
      (fun _ -> Set_snj_row_sum (TA.get_parameter () |> SparseNJ.RowSum.of_string)
        |> List.accum Parameters.program);
    [ "--snj-symmetry" ],
      Some "'one'|'both'",
      [ "K-NN membership policy for the sparse-NJ candidate set" ],
      TA.Default (SparseNJ.Symmetry.to_string Defaults.snj_symmetry |> Fun.const),
      (fun _ -> Set_snj_symmetry (TA.get_parameter () |> SparseNJ.Symmetry.of_string)
        |> List.accum Parameters.program);
    [ "--snj-index-type" ],
      Some "'flat'|'hnsw('<M>')'",
      [ "FAISS index type for sparse-NJ K-NN maintenance" ],
      TA.Default (Interfaiss.Type.to_string Defaults.snj_index_type |> Fun.const),
      (fun _ -> Set_snj_index_type (TA.get_parameter () |> Interfaiss.Type.of_string)
        |> List.accum Parameters.program);
    [ "--hdbscan-min-cluster-size" ],
      Some "<positive_integer>",
      [ "HDBSCAN minimum cluster size (for '--method hdbscan')" ],
      TA.Default (string_of_int Defaults.hdbscan_min_cluster_size |> Fun.const),
      (fun _ -> Set_hdbscan_min_cluster_size (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Branch-length refit (when a counts register is present):" ];
    [ "--refit-with" ],
      Some "'none'|'jaccard'|'jaccard-weighted'",
      [ "branch-length refit source.  'jaccard' (default): mash-like";
        "  distance from binary k-mer Jaccard read off the counts register;";
        "  'jaccard-weighted': weighted (Ruzicka) Jaccard; 'none': keep the";
        "  topology stage's native branch lengths.  Refit runs only when a";
        "  counts register is loaded." ],
      TA.Default (RefitWith.to_string Defaults.refit_with |> Fun.const),
      (fun _ -> Set_refit_with (TA.get_parameter () |> RefitWith.of_string)
        |> List.accum Parameters.program);
    [ "--refit-strategy" ],
      Some "'reuse-knn'|'per-branch'|'all'",
      [ "which pair subset constrains the OLS refit ('reuse-knn' default:";
        "topological K-NN edges of the output tree, sub-quadratic)" ],
      TA.Default (Refit.Strategy.to_string Defaults.refit_strategy |> Fun.const),
      (fun _ -> Set_refit_strategy (TA.get_parameter () |> Refit.Strategy.of_string)
        |> List.accum Parameters.program);
    [ "--refit-num-neighbors" ],
      Some "<positive_integer>",
      [ "K for the 'reuse-knn' refit pair selection" ],
      TA.Default (string_of_int Defaults.refit_num_neighbors |> Fun.const),
      (fun _ -> Set_refit_num_neighbors (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    [ "--refit-k"; "--refit-kmer-size" ],
      Some "<positive_integer>",
      [ "k-mer size used in the mash transform d = -(1/k) ln(2J/(1+J))";
        "(a global scale; affects absolute, not relative, branch lengths)" ],
      TA.Default (string_of_int Defaults.refit_k |> Fun.const),
      (fun _ -> Set_refit_k (TA.get_parameter_int_pos ())
        |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Output and miscellaneous:" ];
    [ "-o"; "--output" ],
      Some "<output_file_prefix>",
      [ "run the workflow selected by the loaded registers and write the";
        "resulting tree to '<prefix>.nwk'" ],
      TA.Mandatory,
      (fun _ -> Compute_tree (TA.get_parameter ()) |> List.accum Parameters.program);
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
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    [ "-x"; "--print-exception-backtrace" ], None, [], TA.Optional,
      (fun _ -> Printexc.record_backtrace true);
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
  Random.self_init ();
  (* Dry run: validate that the requested workflow is runnable. *)
  let has_spectra = ref false and has_twisted = ref false and has_tree = ref false
  and refit_with = ref Defaults.refit_with in
  List.iter
    (function
      | Binary_to_register (RegisterType.Spectra, _) -> has_spectra := true
      | Binary_to_register (RegisterType.Twisted, _) -> has_twisted := true
      | Binary_to_register (RegisterType.Twister, _) -> ()
      | Tree_to_register _ -> has_tree := true
      | Set_refit_with rw -> refit_with := rw
      | Compute_tree _ ->
        if not (!has_tree || !has_twisted || !has_spectra) then
          TA.parse_error
            "Option '-o' needs at least one of: a counts register ('-i s'), a \
             twisted register ('-i t'), or an input tree ('--input-tree')!";
        if !has_tree && not !has_spectra then
          TA.parse_error
            "An input tree was given but no counts register ('-i s'): there is \
             nothing to do (refit needs counts)!"
      | _ -> ())
    program;
  ignore !refit_with;
  (* Execution. *)
  let prefix = String.TermIO.(grey (Printf.sprintf "(%s):" "KPopPhylo")) in
  let spectra = ref (KMerDB.empty ()) and spectra_loaded = ref false
  and twister = ref Twister.empty
  and twisted = ref Twisted.empty and twisted_loaded = ref false
  and tree = ref None
  and metric = ref Defaults.metric and distance = ref Defaults.distance
  and distance_normalize = ref Defaults.distance_normalize
  and method_ = ref Defaults.method_
  and splits_keep_at_most = ref Defaults.splits_keep_at_most
  and centroids_balance_penalty = ref Defaults.centroids_balance_penalty
  and centroids_num_seeds = ref Defaults.centroids_num_seeds
  and centroids_seed = ref Defaults.centroids_seed
  and gaps_prefilter_kneedle = ref Defaults.gaps_prefilter_kneedle
  and hdbscan_min_cluster_size = ref Defaults.hdbscan_min_cluster_size
  and hdbscan_min_samples = ref Defaults.hdbscan_min_samples
  and hdbscan_mst_mode = ref Defaults.hdbscan_mst_mode
  and hdbscan_num_neighbors = ref Defaults.hdbscan_num_neighbors
  and hdbscan_index_type = ref Defaults.hdbscan_index_type
  and hdbscan_lengths_mode = ref Defaults.hdbscan_lengths_mode
  and snj_index_type = ref Defaults.snj_index_type
  and snj_mode = ref Defaults.snj_mode
  and snj_num_neighbors = ref Defaults.snj_num_neighbors
  and snj_k_query_factor = ref Defaults.snj_k_query_factor
  and snj_hyp_scale = ref Defaults.snj_hyp_scale
  and snj_distance = ref Defaults.snj_distance
  and snj_row_sum = ref Defaults.snj_row_sum
  and snj_symmetry = ref Defaults.snj_symmetry
  and refit_with = ref Defaults.refit_with
  and refit_strategy = ref Defaults.refit_strategy
  and refit_num_neighbors = ref Defaults.refit_num_neighbors
  and refit_k = ref Defaults.refit_k in
  let verbose = !Parameters.verbose and threads = !Parameters.threads in
  let announce fmt = Printf.eprintf ("%s " ^^ fmt ^^ "\n%!") prefix in
  List.iter
    (function
      | Binary_to_register (RegisterType.Spectra, p) ->
        spectra := KMerDB.of_binary ~verbose p; spectra_loaded := true
      | Binary_to_register (RegisterType.Twister, p) ->
        twister := Twister.of_binary ~verbose p
      | Binary_to_register (RegisterType.Twisted, p) ->
        twisted := Twisted.of_binary ~verbose p; twisted_loaded := true
      | Tree_to_register p ->
        tree := Some (Trees.Newick.of_file
                        ~negative_branches:Trees.Newick.NegativeBranchesPolicy.Zero p)
      | Set_method m -> method_ := m
      | Set_splits_keep_at_most n -> splits_keep_at_most := n
      | Set_metric m -> metric := m
      | Set_distance d -> distance := d
      | Set_distance_normalize b -> distance_normalize := b
      | Set_centroids_balance_penalty p -> centroids_balance_penalty := p
      | Set_centroids_num_seeds n -> centroids_num_seeds := n
      | Set_centroids_seed s -> centroids_seed := Some s
      | Set_gaps_prefilter_kneedle b -> gaps_prefilter_kneedle := b
      | Set_hdbscan_min_cluster_size n -> hdbscan_min_cluster_size := n
      | Set_hdbscan_min_samples k -> hdbscan_min_samples := Some k
      | Set_hdbscan_mst_mode m -> hdbscan_mst_mode := m
      | Set_hdbscan_num_neighbors k -> hdbscan_num_neighbors := Some k
      | Set_hdbscan_index_type i -> hdbscan_index_type := i
      | Set_hdbscan_lengths_mode m -> hdbscan_lengths_mode := m
      | Set_snj_index_type i -> snj_index_type := i
      | Set_snj_mode m -> snj_mode := m
      | Set_snj_num_neighbors k -> snj_num_neighbors := k
      | Set_snj_k_query_factor k -> snj_k_query_factor := k
      | Set_snj_hyp_scale s -> snj_hyp_scale := s
      | Set_snj_distance d -> snj_distance := d
      | Set_snj_row_sum rs -> snj_row_sum := rs
      | Set_snj_symmetry sym -> snj_symmetry := sym
      | Set_refit_with rw -> refit_with := rw
      | Set_refit_strategy s -> refit_strategy := s
      | Set_refit_num_neighbors k -> refit_num_neighbors := k
      | Set_refit_k k -> refit_k := k
      | Compute_tree out_prefix ->
        (* 1. Obtain the topology. *)
        let topology_tree =
          match !tree with
          | Some t ->
            announce "workflow: refit-only (using the provided input tree)";
            t
          | None ->
            let twisted_for_topology =
              if !twisted_loaded then begin
                announce "workflow: topology from the loaded twisted register%s"
                  (if !spectra_loaded then ", then refit from counts" else "");
                !twisted
              end else begin
                announce "workflow: only counts loaded -- twisting (CA), \
                          then topology, then refit";
                let tw, ts, _ =
                  CA.twist ~threads ~verbose !spectra in
                twister := tw; twisted := ts; twisted_loaded := true;
                ts
              end in
            Twisted.get_phylo_tree
              ~normalize:!distance_normalize ~threads ~verbose
              ~balance_penalty:!centroids_balance_penalty
              ~gaps_prefilter_kneedle:!gaps_prefilter_kneedle
              ~num_seeds:!centroids_num_seeds
              ?seed:!centroids_seed
              ~hdbscan_min_cluster_size:!hdbscan_min_cluster_size
              ?hdbscan_min_samples:!hdbscan_min_samples
              ~hdbscan_mst_mode:!hdbscan_mst_mode
              ?hdbscan_num_neighbors:!hdbscan_num_neighbors
              ~hdbscan_index_type:!hdbscan_index_type
              ~hdbscan_lengths_mode:!hdbscan_lengths_mode
              ~sparse_nj_index_type:!snj_index_type
              ~sparse_nj_mode:!snj_mode
              ~sparse_nj_num_neighbors:!snj_num_neighbors
              ~sparse_nj_k_query_factor:!snj_k_query_factor
              ~sparse_nj_hyp_scale:!snj_hyp_scale
              ~sparse_nj_distance:!snj_distance
              ~sparse_nj_row_sum:!snj_row_sum
              ~sparse_nj_symmetry:!snj_symmetry
              !distance !metric !method_ !splits_keep_at_most twisted_for_topology in
        (* 2. Refit branch lengths if counts present and requested. *)
        let final_tree =
          match !refit_with with
          | RefitWith.NoRefit ->
            if !spectra_loaded then
              announce "refit: skipped ('--refit-with none')";
            topology_tree
          | RefitWith.Mash kind ->
            if not !spectra_loaded then begin
              announce "refit: skipped (no counts register loaded)";
              topology_tree
            end else begin
              announce "refit: branch lengths by sparse OLS vs mash-%s Jaccard \
                        (strategy=%s, K=%d, k-mer=%d)"
                (Jaccard.Kind.to_string kind)
                (Refit.Strategy.to_string !refit_strategy)
                !refit_num_neighbors !refit_k;
              let dist = Jaccard.make_dist ~kind ~k:!refit_k !spectra in
              Refit.refit ~verbose ~strategy:!refit_strategy
                ~k:!refit_num_neighbors ~dist topology_tree
            end in
        Exception.catch_unexpected_end_of_output __FUNCTION__
          (fun () -> Trees.Newick.to_file final_tree (out_prefix ^ ".nwk"));
        announce "wrote %s.nwk" out_prefix)
    program
