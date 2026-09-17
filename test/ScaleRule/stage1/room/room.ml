(* room <twisted.txt> <inertia.txt> <labels> <tag> [<queries> [<k>]]

   How much room each rung of the doubling ladder gives the neighbourhoods, measured from the data
   alone.  At 1, 2, 4, ... axes and at all of them, over a sample of labelled query points:
   - the two-nearest-neighbour intrinsic dimension (Facco et al. 2017), computed as
     Clustering.intrinsic_dimension computes it;
   - how many of each query's k nearest neighbours are still among its k nearest when the axes
     double, and how many are among its k nearest in the full space;
   - how many of them carry the query's own label, which judges the other two and takes no part in
     them.
   Coordinates are weighted by the square root of the inertia, the powers(1,1,1) embedding.  A
   truncation rescales every distance alike, so neighbours and distance ratios do not depend on how
   the metric normalises the truncated inertia.  Only labelled rows are used, as in detect3.  One
   line per rung, tab-separated. *)

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4) in
  let nq_req = if Array.length Sys.argv > 5 then int_of_string Sys.argv.(5) else 1000
  and k = if Array.length Sys.argv > 6 then int_of_string Sys.argv.(6) else 10 in
  if k < 2 then failwith "k must be at least 2";
  let labs = Kio.read_labels labf and iv = Kio.read_inertia inf in
  let names, coords = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs names.(i))
    |> Array.of_list in
  let n = Array.length rows in
  let avail = min (Array.length iv) (Array.length coords.(rows.(0))) in
  let rungs =
    let rec up acc r = if r >= avail then List.rev (avail :: acc) else up (r :: acc) (2 * r) in
    Array.of_list (up [] 1) in
  let nr = Array.length rungs in
  let pc = Array.map (fun r -> Array.init avail (fun a -> coords.(r).(a) *. sqrt iv.(a))) rows in
  let lab =
    let tbl = Hashtbl.create 64 in
    Array.map
      (fun r ->
        let s = Hashtbl.find labs names.(r) in
        match Hashtbl.find_opt tbl s with
        | Some id -> id
        | None ->
          let id = Hashtbl.length tbl in
          Hashtbl.add tbl s id;
          id)
      rows in
  (* A deterministic sample of distinct query points *)
  let st = Random.State.make [| 20260915 |] in
  let perm = Array.init n Fun.id in
  for i = n - 1 downto 1 do
    let j = Random.State.int st (i + 1) in
    let t = perm.(i) in
    perm.(i) <- perm.(j);
    perm.(j) <- t
  done;
  let nq = min nq_req n in
  let knn = Array.init nr (fun _ -> Array.make_matrix nq k (-1))
  and r1 = Array.make_matrix nr nq 0. and r2 = Array.make_matrix nr nq 0. in
  let d2 = Array.make n 0. and bd = Array.make k infinity and bj = Array.make k (-1) in
  for q = 0 to nq - 1 do
    let qi = perm.(q) in
    Array.fill d2 0 n 0.;
    let prev = ref 0 in
    Array.iteri
      (fun ri d ->
        (* Squared distances from the query, the axes of this rung added to those of the last *)
        for a = !prev to d - 1 do
          let x = pc.(qi).(a) in
          for j = 0 to n - 1 do
            let y = x -. pc.(j).(a) in
            d2.(j) <- d2.(j) +. (y *. y)
          done
        done;
        prev := d;
        Array.fill bd 0 k infinity;
        Array.fill bj 0 k (-1);
        for j = 0 to n - 1 do
          if j <> qi then begin
            let v = d2.(j) in
            if v < bd.(k - 1) then begin
              let p = ref (k - 1) in
              while !p > 0 && bd.(!p - 1) > v do
                bd.(!p) <- bd.(!p - 1);
                bj.(!p) <- bj.(!p - 1);
                decr p
              done;
              bd.(!p) <- v;
              bj.(!p) <- j
            end
          end
        done;
        Array.blit bj 0 knn.(ri).(q) 0 k;
        r1.(ri).(q) <- sqrt bd.(0);
        r2.(ri).(q) <- sqrt bd.(1))
      rungs
  done;
  let shared a b =
    let c = ref 0 in
    Array.iter (fun x -> if Array.exists (( = ) x) b then incr c) a;
    !c in
  Printf.printf "# %s: %d labelled spectra, %d queries, k = %d, %d axes available\n" tag n nq k avail;
  Printf.printf "tag\trung\tintrinsic_dimension\tkept_at_next_rung\tshared_with_full\tsame_label\n";
  for ri = 0 to nr - 1 do
    let tot = ref 0. and cnt = ref 0 and nx = ref 0 and fu = ref 0 and sl = ref 0 in
    for q = 0 to nq - 1 do
      if r1.(ri).(q) > 0. then begin
        tot := !tot +. log (r2.(ri).(q) /. r1.(ri).(q));
        incr cnt
      end;
      if ri < nr - 1 then nx := !nx + shared knn.(ri).(q) knn.(ri + 1).(q);
      fu := !fu + shared knn.(ri).(q) knn.(nr - 1).(q);
      let lq = lab.(perm.(q)) in
      Array.iter (fun j -> if lab.(j) = lq then incr sl) knn.(ri).(q)
    done;
    let per x = float_of_int x /. float_of_int (nq * k) in
    Printf.printf "%s\t%d\t%s\t%s\t%.3f\t%.3f\n%!" tag rungs.(ri)
      (if !cnt > 0 && !tot > 0. then Printf.sprintf "%.2f" (float_of_int !cnt /. !tot) else "-")
      (if ri < nr - 1 then Printf.sprintf "%.3f" (per !nx) else "-")
      (per !fu) (per !sl)
  done
