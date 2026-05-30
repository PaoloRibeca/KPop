(* Smoke test for lib/CoverTree.ml.  Validates insert + knn + remove
   against brute-force on random Euclidean points. *)

open BiOCamLib.Better
open KPop

module FloatArray_Metric = struct
  type t = Float.Array.t
  let equal a b =
    let n = Float.Array.length a in
    if Float.Array.length b <> n then false
    else
      let same = ref true in
      let i = ref 0 in
      while !same && !i < n do
        if Float.Array.get a !i <> Float.Array.get b !i then same := false;
        incr i
      done;
      !same
  let dist a b =
    let n = Float.Array.length a in
    let acc = ref 0. in
    for k = 0 to n - 1 do
      let delta = Float.Array.unsafe_get a k -. Float.Array.unsafe_get b k in
      acc := !acc +. delta *. delta
    done;
    sqrt !acc
end

module CT = CoverTree.Make (FloatArray_Metric)

let brute_knn points query k =
  let cand =
    Array.mapi (fun i p -> (FloatArray_Metric.dist p query, i, p)) points in
  Array.sort (fun (a, _, _) (b, _, _) -> compare a b) cand;
  let m = min k (Array.length cand) in
  Array.sub cand 0 m
  |> Array.map (fun (_, _, p) -> p)
  |> Array.to_list

let make_random_points n d seed =
  Random.init seed;
  Array.init n (fun _ ->
    Float.Array.init d (fun _ -> Random.float 1.0 -. 0.5))

let fail fmt =
  Printf.ksprintf (fun s -> Printf.eprintf "FAIL: %s\n%!" s; exit 1) fmt
let pass label =
  Printf.printf "  %-60s PASS\n%!" label

(* Compare two K-NN lists by their point-set (order doesn't matter for
   K-NN equality; we just need the same K points). *)
let knn_set_equal a b =
  let sa = List.map FloatArray_Metric.dist a in
  let _ = sa in
  let in_b p = List.exists (FloatArray_Metric.equal p) b in
  List.length a = List.length b &&
  List.for_all in_b a

(* ============================================================ *)
let () =
  Printf.printf "=== Part 1: insert + knn against brute force (n=50, d=10) ===\n%!";
  let n = 50 and d = 10 and k = 5 in
  let points = make_random_points n d 42 in
  let tree = CT.create () in
  Array.iter (CT.insert tree) points;
  if CT.size tree <> n then
    fail "size mismatch after insert: %d (expected %d)" (CT.size tree) n;
  pass (Printf.sprintf "%d insertions, tree size = %d" n (CT.size tree));
  (* Spot-check K-NN against brute force on a few queries *)
  let ok = ref 0 and ko = ref 0 in
  for q_idx = 0 to 9 do
    let q = points.(q_idx) in
    let ct = CT.knn tree q k in
    let bf = brute_knn points q k in
    if knn_set_equal ct bf then incr ok else begin
      incr ko;
      let ct_dists = List.map (FloatArray_Metric.dist q) ct in
      let bf_dists = List.map (FloatArray_Metric.dist q) bf in
      Printf.eprintf "  query %d: ct dists %s, bf dists %s\n%!" q_idx
        (String.concat ", " (List.map (Printf.sprintf "%.4f") ct_dists))
        (String.concat ", " (List.map (Printf.sprintf "%.4f") bf_dists))
    end
  done;
  if !ko > 0 then fail "K-NN mismatch on %d / 10 queries" !ko;
  pass (Printf.sprintf "K-NN matches brute force on 10/10 queries (k=%d)" k)

(* ============================================================ *)
let () =
  Printf.printf "\n=== Part 2: remove + knn (n=30, d=5, k=3) ===\n%!";
  let n = 30 and d = 5 and k = 3 in
  let points = make_random_points n d 7 in
  let tree = CT.create () in
  Array.iter (CT.insert tree) points;
  (* Remove every other point, verify K-NN on the remaining *)
  let removed = Array.make n false in
  for i = 0 to n - 1 do
    if i mod 2 = 0 then begin
      let r = CT.remove tree points.(i) in
      if not r then fail "remove returned false at i=%d" i;
      removed.(i) <- true
    end
  done;
  let expected_size = (n + 1) / 2 in
  if CT.size tree <> expected_size then
    fail "size after removes: %d (expected %d)" (CT.size tree) expected_size;
  pass (Printf.sprintf "removed %d / %d points; tree size = %d" (n / 2) n (CT.size tree));
  (* Reduced point set *)
  let remaining =
    let acc = ref [] in
    Array.iteri (fun i p -> if not removed.(i) then acc := p :: !acc) points;
    Array.of_list (List.rev !acc) in
  let ok = ref 0 and ko = ref 0 in
  for q_idx = 0 to min 9 (Array.length remaining - 1) do
    let q = remaining.(q_idx) in
    let ct = CT.knn tree q k in
    let bf = brute_knn remaining q k in
    if knn_set_equal ct bf then incr ok else incr ko
  done;
  if !ko > 0 then fail "post-remove K-NN mismatch on %d queries" !ko;
  pass "K-NN matches brute force after removes"

(* ============================================================ *)
let () =
  Printf.printf "\n=== Part 3: stress test (n=200, d=20, k=10, interleaved insert / remove / query) ===\n%!";
  let n = 200 and d = 20 and k = 10 in
  let points = make_random_points n d 17 in
  let tree = CT.create () in
  Array.iter (CT.insert tree) points;
  (* Randomly remove half *)
  Random.init 17;
  let removed = Array.make n false in
  for _ = 0 to n / 2 - 1 do
    let i = Random.int n in
    if not removed.(i) then begin
      let _ = CT.remove tree points.(i) in
      removed.(i) <- true
    end
  done;
  let remaining =
    let acc = ref [] in
    Array.iteri (fun i p -> if not removed.(i) then acc := p :: !acc) points;
    Array.of_list (List.rev !acc) in
  let r_size = Array.length remaining in
  pass (Printf.sprintf "post-remove tree size = %d (expected = %d)" (CT.size tree) r_size);
  if CT.size tree <> r_size then
    fail "size disagrees with expected after random removes";
  (* Validate K-NN over the remaining set *)
  let ok = ref 0 in
  for q_idx = 0 to min 19 (r_size - 1) do
    let q = remaining.(q_idx) in
    let ct = CT.knn tree q k in
    let bf = brute_knn remaining q k in
    if knn_set_equal ct bf then incr ok
  done;
  if !ok < 20 then fail "K-NN mismatch: %d / 20 ok after interleaved ops" !ok;
  pass (Printf.sprintf "K-NN matches brute force on 20/20 queries after random removes")

let () =
  Printf.printf "\nAll CoverTree tests passed.\n%!"
