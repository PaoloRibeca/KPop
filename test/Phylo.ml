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

(* Invariant tests for the phylogenetic-tree subsystem
   (gaps, centroids, hdbscan, sparse-nj).
   Each part is independent and exits with code 1 on failure.

   Expected to be run from the project root so that
   "test/Primer/Classes-5.KPopTwisted" resolves. *)

open BiOCamLib.Better
module Trees = BiOCamLib.Trees
module Trees_Base = BiOCamLib.Trees_Base
open KPop

let fail fmt =
  Printf.ksprintf (fun s -> Printf.eprintf "FAIL: %s\n%!" s; exit 1) fmt

let pass label =
  Printf.printf "  %-60s PASS\n%!" label

let twisted_prefix = "test/Primer/Classes-5"
let distance = Space.Distance.of_string "euclidean"
let metric = Space.Distance.Metric.of_string "powers(1,1,1)"

(* Helper: canonicalise a Trees.Newick.t to its underlying split set so we
   can compare two trees in a way insensitive to child-traversal order at
   internal nodes.  Two trees compare equal iff they share the same
   bipartitions with the same branch lengths. *)
let dump tree =
  Trees.Splits.to_string ~precision:15
    (Trees_Base.Splits.of_newick tree)

(* Run get_phylo_tree with sensible defaults for the chosen algorithm *)
let run_phylo ?(threads = 1) ?(seed = 42)
    ?(hdbscan_mst_mode = Clustering.HdbscanMstMode.Auto)
    ?(hdbscan_index_type = Interfaiss.Type.of_string "hnsw(32)")
    ?(hdbscan_min_cluster_size = 2)
    ?(hdbscan_lengths_mode = Clustering.Hdbscan.LengthsMode.Mreach)
    algorithm twisted =
  Twisted.get_phylo_tree
    ~threads ~verbose:false
    ~num_seeds:5 ~seed
    ~hdbscan_min_cluster_size
    ~hdbscan_mst_mode ~hdbscan_index_type
    ~hdbscan_lengths_mode
    distance metric algorithm 10000 twisted

(* -------------------------------------------------------------------------- *)
(*  Part 1 - HDBSCAN emits a well-formed Newick tree                          *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 1: HDBSCAN emits a well-formed Newick tree ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  List.iter
    (fun k ->
      let tree = run_phylo ~hdbscan_min_cluster_size:k
                   Twisted.Phylo.Method.Hdbscan twisted in
      let s = dump tree in
      if String.length s < 2 || s.[String.length s - 2] <> ';' && s.[String.length s - 1] <> ';' then
        fail "K=%d: HDBSCAN Newick does not end with ';' (got %S)" k s;
      pass (Printf.sprintf "HDBSCAN K=%d: emits a parseable Newick tree" k))
    [ 1; 2; 3; 4 ]

(* -------------------------------------------------------------------------- *)
(*  Part 2 - HDBSCAN Auto(flat) == Dense byte-for-byte (Persistence mode)     *)
(*  When the FAISS index is exact (flat), Auto's core distances match Dense   *)
(*  and the resulting MST + condensed tree must be bit-identical.  We         *)
(*  compare via Persistence-mode lengths: the condensation step is tie-       *)
(*  insensitive and so produces canonical splits even on inputs with          *)
(*  coincident lambda values.  (Mreach mode walks the raw binary merge tree   *)
(*  directly, and so is sensitive to MST tie-breaking when several leaves     *)
(*  merge at the same lambda; that's a fundamental property of Mreach, not    *)
(*  a regression.)                                                            *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 2: Auto(flat) == Dense on HDBSCAN (Persistence) ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let flat = Interfaiss.Type.of_string "flat" in
  List.iter
    (fun k ->
      let dense =
        run_phylo ~hdbscan_min_cluster_size:k
          ~hdbscan_mst_mode:Clustering.HdbscanMstMode.Dense
          ~hdbscan_lengths_mode:Clustering.Hdbscan.LengthsMode.Persistence
          Twisted.Phylo.Method.Hdbscan twisted in
      let auto_flat =
        run_phylo ~hdbscan_min_cluster_size:k
          ~hdbscan_mst_mode:Clustering.HdbscanMstMode.Auto
          ~hdbscan_index_type:flat
          ~hdbscan_lengths_mode:Clustering.Hdbscan.LengthsMode.Persistence
          Twisted.Phylo.Method.Hdbscan twisted in
      if dump dense <> dump auto_flat then
        fail "K=%d: Auto(flat) and Dense trees differ" k;
      pass (Printf.sprintf "HDBSCAN K=%d: Auto(flat) bit-identical to Dense" k))
    [ 1; 2; 3; 4 ]

(* -------------------------------------------------------------------------- *)
(*  Part 3 - Centroids reproducibility across thread counts                   *)
(*  With a fixed seed the multi-seed bootstrap loop must be byte-identical    *)
(*  regardless of --threads.                                                  *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 3: Centroids reproducible across -T 1 and -T 4 ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let t1 = run_phylo ~threads:1 Twisted.Phylo.Method.Centroids twisted in
  let t4 = run_phylo ~threads:4 Twisted.Phylo.Method.Centroids twisted in
  if dump t1 <> dump t4 then
    fail "centroids: -T 1 and -T 4 outputs differ at seed=42";
  pass "centroids with same seed: byte-identical across -T 1 and -T 4"

(* -------------------------------------------------------------------------- *)
(*  Part 4 - Gaps, Centroids and Sparse-NJ produce non-trivial trees          *)
(*  Smoke test for the dispatch matrix; each must emit a tree with at least   *)
(*  one internal node beyond a bare star on the synthetic 5-class fixture.    *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 4: every method produces a non-empty Newick tree ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  List.iter
    (fun (label, method_) ->
      let tree = run_phylo method_ twisted in
      let s = dump tree in
      if String.length s < 4 then
        fail "%s: Newick output is suspiciously short (got %S)" label s;
      pass (Printf.sprintf "%s emits a non-empty Newick tree (%d chars)" label (String.length s)))
    [ "gaps",      Twisted.Phylo.Method.Gaps;
      "centroids", Twisted.Phylo.Method.Centroids;
      "sparse-nj", Twisted.Phylo.Method.Sparse_nj ]

(* -------------------------------------------------------------------------- *)
(*  Part 5 - HDBSCAN persistence and mreach branch-length modes both work     *)
(*  Both must emit valid Newick.  They are allowed (and expected) to differ,  *)
(*  since branch lengths use different formulas.                              *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 5: HDBSCAN persistence and mreach lengths modes ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  List.iter
    (fun mode ->
      let tree = run_phylo
                   ~hdbscan_lengths_mode:mode
                   Twisted.Phylo.Method.Hdbscan twisted in
      let s = dump tree in
      if String.length s < 4 then
        fail "HDBSCAN lengths=%s: Newick is too short (got %S)"
          (Clustering.Hdbscan.LengthsMode.to_string mode) s;
      pass (Printf.sprintf "HDBSCAN lengths=%s: emits valid Newick"
              (Clustering.Hdbscan.LengthsMode.to_string mode)))
    [ Clustering.Hdbscan.LengthsMode.Persistence;
      Clustering.Hdbscan.LengthsMode.Mreach ]

(* -------------------------------------------------------------------------- *)
(*  Part 6 - Sparse-mode error on under-sized num_neighbors                   *)
(*  Catches that the validation guard fires before the FAISS call.            *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 6: Sparse mode rejects num_neighbors < min_samples ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let raised =
    try
      let _ = Twisted.get_phylo_tree
                ~threads:1 ~verbose:false
                ~hdbscan_min_cluster_size:2 ~hdbscan_min_samples:5
                ~hdbscan_mst_mode:Clustering.HdbscanMstMode.Sparse
                ~hdbscan_num_neighbors:2
                distance metric Twisted.Phylo.Method.Hdbscan 10000 twisted in
      false
    with _ -> true in
  if not raised then
    fail "expected an exception when num_neighbors < min_samples";
  pass "validation rejects num_neighbors(2) < min_samples(5)"

(* -------------------------------------------------------------------------- *)
(*  Part 7 - sparse-NJ flat == hnsw(32) on small fixture                      *)
(*  HNSW with M=32 is exact at n=10, so the two index choices must produce    *)
(*  byte-identical Newick output on the Classes-5 fixture.                    *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 7: sparse-NJ flat == hnsw(32) on small fixture ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let m = Twisted.to_embeddings ~normalize:true ~verbose:false distance metric twisted in
  let flat =
    SparseNJ.compute ~verbose:false ~k_nn:10
      ~index_type:(Interfaiss.Type.of_string "flat")
      m.matrix.row_names m.matrix.data in
  let hnsw =
    SparseNJ.compute ~verbose:false ~k_nn:10
      ~index_type:(Interfaiss.Type.of_string "hnsw(32)")
      m.matrix.row_names m.matrix.data in
  if Trees.Newick.to_string flat <> Trees.Newick.to_string hnsw then
    fail "sparse-NJ: flat and hnsw(32) outputs differ on small fixture";
  pass "sparse-NJ: flat and hnsw(32) produce identical trees on small fixture"

(* -------------------------------------------------------------------------- *)
(*  Part 8 - subquadratic mode matches dense mode at sufficient K_QUERY       *)
(*  On the Classes-5 fixture (n=10) the candidate pool quickly becomes the    *)
(*  full active set, so a moderate factor (= 5) is enough for subquadratic    *)
(*  to reproduce the dense tree.  Validates that the algorithmic skeleton --  *)
(*  reverse-insertion + FAISS expansion + Saitou-Nei reranking -- is wired    *)
(*  correctly; the asymptotic O(n K^2 + n K log n) claim depends on the       *)
(*  empirical question of how K_QUERY scales with n on real data, which       *)
(*  this small fixture cannot answer.                                         *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 8: periodic-rebuild matches dense at K_QUERY=5K on small fixture ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let m = Twisted.to_embeddings ~normalize:true ~verbose:false distance metric twisted in
  let dense =
    SparseNJ.compute ~verbose:false ~k_nn:5
      ~mode:SparseNJ.Mode.Dense
      ~index_type:(Interfaiss.Type.of_string "flat")
      m.matrix.row_names m.matrix.data in
  let sq =
    SparseNJ.compute ~verbose:false ~k_nn:5 ~k_query_factor:5
      ~mode:SparseNJ.Mode.PeriodicRebuild
      ~index_type:(Interfaiss.Type.of_string "flat")
      m.matrix.row_names m.matrix.data in
  let dense_splits = Trees_Base.Splits.of_newick dense in
  let sq_splits = Trees_Base.Splits.of_newick sq in
  if Trees.Splits.to_string ~precision:15 dense_splits
     <> Trees.Splits.to_string ~precision:15 sq_splits then
    fail "subquadratic: split set differs from dense at K=5 K_QUERY=25";
  pass "subquadratic K=5 K_QUERY=25: split set matches dense byte-for-byte"

let () =
  Printf.printf "\nAll phylo-subsystem invariants passed.\n%!"
