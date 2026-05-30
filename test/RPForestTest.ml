(* Smoke test for lib/RPForest.ml.  Validates that the forest's K-NN
   results overlap substantially with brute-force on synthetic
   Euclidean data; RP-forest is approximate, so we measure recall
   rather than exact equality. *)

open BiOCamLib.Better
open KPop

module Pt = struct
  type t = int * Float.Array.t
  let equal (i, _) (j, _) = i = j
  let hash (i, _) = i
  let embedding (_, e) = e
end

module Forest = RPForest.Make (Pt)

let eucl_dist a b =
  let n = Float.Array.length a in
  let acc = ref 0. in
  for k = 0 to n - 1 do
    let d = Float.Array.unsafe_get a k -. Float.Array.unsafe_get b k in
    acc := !acc +. d *. d
  done;
  sqrt !acc

let brute_knn points query k =
  let arr =
    Array.mapi (fun i p -> (eucl_dist (snd p) (snd query), i, p)) points in
  Array.sort (fun (a, _, _) (b, _, _) -> compare a b) arr;
  let m = min k (Array.length arr) in
  Array.sub arr 0 m |> Array.map (fun (_, _, p) -> p) |> Array.to_list

let make_random_points n d seed =
  Random.init seed;
  Array.init n (fun i ->
    let e = Float.Array.init d (fun _ -> Random.float 1.0 -. 0.5) in
    (i, e))

let fail fmt =
  Printf.ksprintf (fun s -> Printf.eprintf "FAIL: %s\n%!" s; exit 1) fmt
let pass label =
  Printf.printf "  %-60s PASS\n%!" label

let recall ct bf =
  let in_bf p =
    List.exists (fun q -> Pt.equal p q) bf in
  let hits = List.fold_left (fun acc p -> if in_bf p then acc + 1 else acc) 0 ct in
  if List.length bf = 0 then 1.0
  else float_of_int hits /. float_of_int (List.length bf)

(* ============================================================ *)
let () =
  Printf.printf "=== Part 1: build + knn recall (n=200, d=20, k=10) ===\n%!";
  let n = 200 and d = 20 and k = 10 in
  let points = make_random_points n d 42 in
  let forest =
    Forest.build ~leaf_size:16 ~n_trees:10 ~seed:31337 points in
  if Forest.size forest <> n then
    fail "size mismatch: %d (expected %d)" (Forest.size forest) n;
  pass (Printf.sprintf "built forest of size %d" (Forest.size forest));
  let total_recall = ref 0. in
  let n_queries = 20 in
  for q_idx = 0 to n_queries - 1 do
    let q = points.(q_idx) in
    let ct = Forest.knn forest q k in
    let bf = brute_knn points q k in
    total_recall := !total_recall +. recall ct bf
  done;
  let mean_recall = !total_recall /. float_of_int n_queries in
  Printf.printf "  mean K-NN recall over %d queries (k=%d): %.3f\n%!"
    n_queries k mean_recall;
  if mean_recall < 0.5 then
    fail "recall too low: %.3f" mean_recall;
  pass "recall >= 0.5"

(* ============================================================ *)
let () =
  Printf.printf "\n=== Part 2: more trees -> higher recall (n=200, d=20, k=10) ===\n%!";
  let n = 200 and d = 20 and k = 10 in
  let points = make_random_points n d 11 in
  let measure n_trees =
    let forest = Forest.build ~leaf_size:16 ~n_trees ~seed:42 points in
    let total = ref 0. in
    let qs = 20 in
    for q_idx = 0 to qs - 1 do
      let q = points.(q_idx) in
      let ct = Forest.knn forest q k in
      let bf = brute_knn points q k in
      total := !total +. recall ct bf
    done;
    !total /. float_of_int qs in
  let r5 = measure 5 in
  let r20 = measure 20 in
  Printf.printf "  recall at n_trees=5:  %.3f\n%!" r5;
  Printf.printf "  recall at n_trees=20: %.3f\n%!" r20;
  if r20 < r5 -. 0.05 then
    fail "more trees should not decrease recall (got %.3f -> %.3f)" r5 r20;
  pass "recall monotone in n_trees"

(* ============================================================ *)
let () =
  Printf.printf "\n=== Part 3: high-d recall depends on n_trees (n=400, d=50) ===\n%!";
  (* RP-forest recall degrades with embedding dimension; this part
     documents the sensitivity rather than asserting a fixed bound.
     The sparse-NJ caller compensates via more trees and an over-fetch
     factor in K-NN queries. *)
  let n = 400 and d = 50 and k = 10 in
  let kq = 3 * k in
  let points = make_random_points n d 7 in
  let measure n_trees =
    let forest = Forest.build ~leaf_size:16 ~n_trees ~seed:99 points in
    let total = ref 0. in
    let qs = 20 in
    for q_idx = 0 to qs - 1 do
      let q = points.(q_idx) in
      let ct = Forest.knn forest q kq in
      let bf = brute_knn points q k in
      total := !total +. recall ct bf
    done;
    !total /. float_of_int qs in
  let r10 = measure 10 in
  let r25 = measure 25 in
  let r50 = measure 50 in
  Printf.printf
    "  over-fetched (3k) recall vs true %d-NN at n_trees=10 / 25 / 50: %.3f / %.3f / %.3f\n%!"
    k r10 r25 r50;
  if r50 < 0.8 then
    fail "recall at n_trees=50 should reach >= 0.8 (got %.3f)" r50;
  pass "high n_trees recovers high recall at d=50"

let () =
  Printf.printf "\nAll RPForest tests passed.\n%!"
