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

(* Invariant tests for the splits subsystem (gaps, centroids, hdbscan).
   Each part is independent and exits with code 1 on failure.

   Expected to be run from the project root so that
   "test/Primer/Classes-5.KPopTwisted" resolves. *)

open BiOCamLib.Better
module Trees = BiOCamLib.Trees
open KPop

let fail fmt =
  Printf.ksprintf (fun s -> Printf.eprintf "FAIL: %s\n%!" s; exit 1) fmt

let pass label =
  Printf.printf "  %-60s PASS\n%!" label

let twisted_prefix = "test/Primer/Classes-5"
let distance = Space.Distance.of_string "euclidean"
let metric = Space.Distance.Metric.of_string "powers(1,1,1)"

(* Helper: stringify a Trees.Splits.t so we can compare two instances exactly. *)
let dump splits =
  Trees.Splits.to_string ~precision:15 splits

(* Run get_splits with sensible defaults for the chosen algorithm *)
let run_splits ?(threads = 1) ?(seed = 42)
    ?(hdbscan_mst_mode = Clustering.HdbscanMstMode.Auto)
    ?(hdbscan_index_type = Interfaiss.Type.of_string "hnsw(32)")
    ?(hdbscan_min_cluster_size = 2)
    algorithm twisted =
  Twisted.get_splits
    ~threads ~verbose:false
    ~num_seeds:5 ~seed
    ~hdbscan_min_cluster_size
    ~hdbscan_mst_mode ~hdbscan_index_type
    distance metric algorithm 10000 twisted

(* -------------------------------------------------------------------------- *)
(*  Part 1 - Empty Trees.Splits round-trip                                    *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 1: empty Trees.Splits round-trip ===\n%!";
  let names = [| "a"; "b"; "c" |] in
  let empty = Trees.Splits.create names in
  let s = dump empty in
  let parsed = Trees.Splits.of_string s in
  if Trees.Splits.cardinal parsed <> 0 then
    fail "empty round-trip: cardinal = %d (expected 0)" (Trees.Splits.cardinal parsed);
  if Trees.Splits.get_names parsed <> names then
    fail "empty round-trip: names mismatch (got %s)"
      (String.concat "," (Array.to_list (Trees.Splits.get_names parsed)));
  pass "empty splits text round-trip preserves names + cardinal 0"

(* -------------------------------------------------------------------------- *)
(*  Part 2 - HDBSCAN emits jointly-compatible splits                          *)
(*  Buneman: a set of splits is tree-realisable iff every pair is compatible. *)
(*  Trees.Splits.to_tree partitions input into (kept, tree, dropped); for     *)
(*  HDBSCAN output the dropped set must be empty (nested-by-construction).    *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 2: HDBSCAN splits are jointly compatible ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  List.iter
    (fun k ->
      let splits = run_splits ~hdbscan_min_cluster_size:k
                     Twisted.Splits.Method.Hdbscan twisted in
      let _, _, dropped = Trees.Splits.to_tree splits in
      let nd = Trees.Splits.cardinal dropped in
      if nd <> 0 then
        fail "K=%d: %d incompatible splits dropped (expected 0)" k nd;
      pass (Printf.sprintf "HDBSCAN K=%d: 0 dropped on Trees.Splits.to_tree" k))
    [ 1; 2; 3; 4 ]

(* -------------------------------------------------------------------------- *)
(*  Part 3 - HDBSCAN Auto(flat) == Dense byte-for-byte                        *)
(*  When the FAISS index is exact (flat), Auto's core distances match Dense   *)
(*  and the resulting MST + condensed tree + splits are bit-identical.        *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 3: Auto(flat) == Dense on HDBSCAN ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let flat = Interfaiss.Type.of_string "flat" in
  List.iter
    (fun k ->
      let dense =
        run_splits ~hdbscan_min_cluster_size:k
          ~hdbscan_mst_mode:Clustering.HdbscanMstMode.Dense
          Twisted.Splits.Method.Hdbscan twisted in
      let auto_flat =
        run_splits ~hdbscan_min_cluster_size:k
          ~hdbscan_mst_mode:Clustering.HdbscanMstMode.Auto
          ~hdbscan_index_type:flat
          Twisted.Splits.Method.Hdbscan twisted in
      if dump dense <> dump auto_flat then
        fail "K=%d: Auto(flat) and Dense splits differ" k;
      pass (Printf.sprintf "HDBSCAN K=%d: Auto(flat) bit-identical to Dense" k))
    [ 1; 2; 3; 4 ]

(* -------------------------------------------------------------------------- *)
(*  Part 4 - Centroids reproducibility across thread counts                   *)
(*  With a fixed seed the multi-seed bootstrap loop must be byte-identical    *)
(*  regardless of --threads.                                                  *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 4: Centroids reproducible across -T 1 and -T 4 ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let t1 = run_splits ~threads:1 Twisted.Splits.Method.Centroids twisted in
  let t4 = run_splits ~threads:4 Twisted.Splits.Method.Centroids twisted in
  if dump t1 <> dump t4 then
    fail "centroids: -T 1 and -T 4 outputs differ at seed=42";
  pass "centroids with same seed: byte-identical across -T 1 and -T 4"

(* -------------------------------------------------------------------------- *)
(*  Part 5 - Gaps and Centroids produce non-trivial output on Classes-5       *)
(*  Smoke test for the dispatch matrix; both should emit at least one split   *)
(*  given the synthetic 5-class structure of the fixture.                     *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 5: Gaps and Centroids produce splits on Classes-5 ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let gaps = run_splits Twisted.Splits.Method.Gaps twisted in
  let cen = run_splits Twisted.Splits.Method.Centroids twisted in
  if Trees.Splits.cardinal gaps = 0 then
    fail "gaps: zero splits emitted on Classes-5";
  if Trees.Splits.cardinal cen = 0 then
    fail "centroids: zero splits emitted on Classes-5";
  pass (Printf.sprintf "gaps emits %d splits" (Trees.Splits.cardinal gaps));
  pass (Printf.sprintf "centroids emits %d splits" (Trees.Splits.cardinal cen))

(* -------------------------------------------------------------------------- *)
(*  Part 6 - Sparse-mode error on under-sized num_neighbors                   *)
(*  Catches that the validation guard fires before the FAISS call.            *)
(* -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 6: Sparse mode rejects num_neighbors < min_samples ===\n%!";
  let twisted = Twisted.of_binary twisted_prefix in
  let raised =
    try
      let _ = Twisted.get_splits
                ~threads:1 ~verbose:false
                ~hdbscan_min_cluster_size:2 ~hdbscan_min_samples:5
                ~hdbscan_mst_mode:Clustering.HdbscanMstMode.Sparse
                ~hdbscan_num_neighbors:2
                distance metric Twisted.Splits.Method.Hdbscan 10000 twisted in
      false
    with _ -> true in
  if not raised then
    fail "expected an exception when num_neighbors < min_samples";
  pass "validation rejects num_neighbors(2) < min_samples(5)"

let () =
  Printf.printf "\nAll splits-subsystem invariants passed.\n%!"

