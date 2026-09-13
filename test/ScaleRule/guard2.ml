(* guard2 <twisted> <inertia> <labels> <d> <tag> <tsv>
   Guard variants for the silhouettes the Monte-Carlo partition search maximises, scored on
   reference, degenerate and peeled partitions of the labelled spectra.
   Embedding = first d coordinates times sqrt(inertia), as cs.ml; the library's powers(1,1,1)
   weights are d * inertia / sum(inertia), a global factor both silhouettes ignore.
   Silhouettes, averaged over all n points:
     classical   a = mean distance to the other members of the own cluster, b = smallest mean
                 distance to the members of another cluster;
     simplified  a = distance to the own centroid (the point included, as in the library),
                 b = smallest distance to another centroid.  Computed from the distance matrix
                 through |x - mean C|^2 = sum_{j in C} |x - x_j|^2 / |C| - sum_{j<l in C} |x_j - x_l|^2 / |C|^2
   s(a, b) = (b - a) / max(a, b); 0 when b is undefined or a = b = 0.
   Guards at size s (a cluster is ELIGIBLE when it has at least s members):
     G0   no guard: clusters of >= 2 members score against every other cluster, singletons 0;
     G1   as G0, but only eligible clusters can be b (no eligible other cluster: 0);
     G2   only members of eligible clusters score, only eligible clusters can be b, others 0;
     G3a  projection: every member of an ineligible cluster is reassigned to the eligible
          cluster nearest by the silhouette's own dissimilarity (centroid distance for the
          simplified, mean distance for the classical), measured on the partition as given;
          the score is G0 of the resulting partition, all of whose clusters are eligible
          (0 when there are fewer than two eligible clusters);
     G3b  as G2, but members of ineligible clusters score -1 *)
type res = {
  k : int; fp : float; noise : int array; ne : int array;
  g0 : float array; g1 : float array array; g2 : float array array;
  g3a : float array array; g3b : float array array
}

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3)
  and dreq = int_of_string Sys.argv.(4) and tag = Sys.argv.(5) and tsvf = Sys.argv.(6) in
  let log fmt = Printf.ksprintf (fun s -> Printf.eprintf "[%s %6.1fs] %s\n%!" tag (Sys.time ()) s) fmt in
  let labs = Kio.read_labels labf and iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs names_all.(i)) |> Array.of_list in
  let n = Array.length rows in
  let d = min dreq (Array.length iv) in
  let x = Array.map (fun r -> Array.init d (fun k -> coords_all.(r).(k) *. sqrt iv.(k))) rows in
  let ss = [| 4; 8; 16 |] in
  let nt = Array.length ss and nsil = 2 in
  let dist a b =
    let s = ref 0. in
    for k = 0 to d - 1 do let t = a.(k) -. b.(k) in s := !s +. (t *. t) done;
    sqrt !s in
  let np = n * (n - 1) / 2 in
  let base = Array.init n (fun i -> (i * n) - (i * (i + 1) / 2)) in
  let dm = Float.Array.make np 0. in
  (let idx = ref 0 in
   for i = 0 to n - 2 do
     let xi = x.(i) in
     for j = i + 1 to n - 1 do Float.Array.unsafe_set dm !idx (dist xi x.(j)); incr idx done
   done);
  let dd i j =
    if i = j then 0.
    else
      let i, j = if i < j then i, j else j, i in
      Float.Array.unsafe_get dm (base.(i) + (j - i - 1)) in
  log "n=%d d=%d, distance matrix done" n d;
  let canon p =
    let h = Hashtbl.create 64 in
    let q =
      Array.map
        (fun l ->
          match Hashtbl.find_opt h l with
          | Some v -> v
          | None -> let v = Hashtbl.length h in Hashtbl.add h l v; v)
        p in
    q, Hashtbl.length h in
  let sizes p k = let s = Array.make k 0 in Array.iter (fun c -> s.(c) <- s.(c) + 1) p; s in
  let members p k =
    let l = Array.make k [] in
    for i = n - 1 downto 0 do l.(p.(i)) <- i :: l.(p.(i)) done;
    Array.map Array.of_list l in
  let centroids p k =
    let c = Array.make_matrix k d 0. and s = sizes p k in
    Array.iteri (fun i ci -> for t = 0 to d - 1 do c.(ci).(t) <- c.(ci).(t) +. x.(i).(t) done) p;
    Array.iteri (fun ci v -> if s.(ci) > 0 then for t = 0 to d - 1 do v.(t) <- v.(t) /. float s.(ci) done) c;
    c in
  let unit_centroid un =
    let c = Array.make d 0. in
    Array.iter (fun i -> for t = 0 to d - 1 do c.(t) <- c.(t) +. x.(i).(t) done) un;
    let f = float (Array.length un) in
    Array.map (fun v -> v /. f) c in
  let sil a b =
    if b = infinity then 0.
    else let mx = Float.max a b in if mx > 0. then (b -. a) /. mx else 0. in
  (* THE SCORING PASS: one O(n^2) sweep accumulating, for every point, the sums of distances and
     of squared distances to every cluster, from which both silhouettes follow under every guard *)
  let s1 = Array.make n 0. and s2 = Array.make n 0. in
  let bcl = Array.make nt infinity and bsi = Array.make nt infinity
  and rcl = Array.make nt (-1) and rsi = Array.make nt (-1) in
  let base_score ?pts p k want_proj =
    let sz = sizes p k and mem = members p k in
    let w = Array.make k 0. in
    Array.iteri
      (fun c ms ->
        let m = Array.length ms in
        if m >= 2 then begin
          let acc = ref 0. in
          for u = 0 to m - 2 do
            let iu = ms.(u) in
            for v = u + 1 to m - 1 do let t = dd iu ms.(v) in acc := !acc +. (t *. t) done
          done;
          w.(c) <- !acc /. (float m *. float m)
        end)
      mem;
    let ne = Array.map (fun s -> Array.fold_left (fun acc m -> if m >= s then acc + 1 else acc) 0 sz) ss
    and noise = Array.map (fun s -> Array.fold_left (fun acc m -> if m < s then acc + m else acc) 0 sz) ss in
    let proj = Array.init nt (fun t -> want_proj && noise.(t) > 0 && ne.(t) >= 2) in
    let q = Array.init nsil (fun _ -> Array.init nt (fun t -> if proj.(t) then Array.copy p else [||])) in
    let g0 = Array.make nsil 0. and g1 = Array.make_matrix nsil nt 0.
    and g2 = Array.make_matrix nsil nt 0. and g3b = Array.make_matrix nsil nt 0. in
    for i = 0 to n - 1 do
      Array.fill s1 0 k 0.;
      Array.fill s2 0 k 0.;
      for j = 0 to i - 1 do
        let t = Float.Array.unsafe_get dm (base.(j) + i - j - 1) and c = p.(j) in
        s1.(c) <- s1.(c) +. t;
        s2.(c) <- s2.(c) +. (t *. t)
      done;
      let bi = base.(i) - i - 1 in
      for j = i + 1 to n - 1 do
        let t = Float.Array.unsafe_get dm (bi + j) and c = p.(j) in
        s1.(c) <- s1.(c) +. t;
        s2.(c) <- s2.(c) +. (t *. t)
      done;
      let c = p.(i) in
      let m = sz.(c) in
      let acl = if m >= 2 then s1.(c) /. float (m - 1) else 0.
      and asi = if m >= 2 then sqrt (Float.max 0. ((s2.(c) /. float m) -. w.(c))) else 0. in
      let b0cl = ref infinity and b0si = ref infinity in
      Array.fill bcl 0 nt infinity;
      Array.fill bsi 0 nt infinity;
      Array.fill rcl 0 nt (-1);
      Array.fill rsi 0 nt (-1);
      for c' = 0 to k - 1 do
        if c' <> c then begin
          let mc = sz.(c') in
          let fm = float mc in
          let vcl = s1.(c') /. fm and vsi = sqrt (Float.max 0. ((s2.(c') /. fm) -. w.(c'))) in
          if vcl < !b0cl then b0cl := vcl;
          if vsi < !b0si then b0si := vsi;
          for t = 0 to nt - 1 do
            if mc >= ss.(t) then begin
              if vcl < bcl.(t) then begin bcl.(t) <- vcl; rcl.(t) <- c' end;
              if vsi < bsi.(t) then begin bsi.(t) <- vsi; rsi.(t) <- c' end
            end
          done
        end
      done;
      if m >= 2 then begin
        g0.(0) <- g0.(0) +. sil acl !b0cl;
        g0.(1) <- g0.(1) +. sil asi !b0si
      end;
      (match pts with
       | Some pa ->
         pa.(0).(i) <- (if m >= 2 then sil acl !b0cl else 0.);
         pa.(1).(i) <- (if m >= 2 then sil asi !b0si else 0.)
       | None -> ());
      for t = 0 to nt - 1 do
        if m >= 2 then begin
          g1.(0).(t) <- g1.(0).(t) +. sil acl bcl.(t);
          g1.(1).(t) <- g1.(1).(t) +. sil asi bsi.(t)
        end;
        if m >= ss.(t) then begin
          let vc = sil acl bcl.(t) and vs = sil asi bsi.(t) in
          g2.(0).(t) <- g2.(0).(t) +. vc;
          g2.(1).(t) <- g2.(1).(t) +. vs;
          g3b.(0).(t) <- g3b.(0).(t) +. vc;
          g3b.(1).(t) <- g3b.(1).(t) +. vs
        end else begin
          g3b.(0).(t) <- g3b.(0).(t) -. 1.;
          g3b.(1).(t) <- g3b.(1).(t) -. 1.;
          if proj.(t) then begin q.(0).(t).(i) <- rcl.(t); q.(1).(t).(i) <- rsi.(t) end
        end
      done
    done;
    let fn = float n in
    let dv a = Array.map (fun r -> Array.map (fun v -> v /. fn) r) a in
    sz, ne, noise, Array.map (fun v -> v /. fn) g0, dv g1, dv g2, dv g3b, q in
  let cache = ref [] in
  let g0_of q =
    match List.assoc_opt q !cache with
    | Some v -> v
    | None ->
      let qc, kq = canon q in
      let _, _, _, g0, _, _, _, _ = base_score qc kq false in
      cache := (q, g0) :: !cache;
      g0 in
  let fp_of sz =
    Array.fold_left (fun acc v -> acc +. (float v *. float (v - 1))) 0. sz /. (float n *. float (n - 1)) in
  let score_q p k =
    cache := [];
    let sz, ne, noise, g0, g1, g2, g3b, q = base_score p k true in
    let g3a =
      Array.init nsil (fun sl ->
          Array.init nt (fun t ->
              if noise.(t) = 0 then g0.(sl)
              else if ne.(t) <= 1 then 0.
              else (g0_of q.(sl).(t)).(sl))) in
    { k; fp = fp_of sz; noise; ne; g0; g1; g2; g3a; g3b }, q in
  let score p k = fst (score_q p k) in
  let oc = open_out tsvf in
  let results = ref [] in
  let add_res cat name r =
    results := (cat, name, r) :: !results;
    let silname = [| "classical"; "simplified" |] in
    for sl = 0 to nsil - 1 do
      Printf.fprintf oc "%s\t%d\t%s\t%s\t%d\t%.6f\t%d\t%d\t%d\t%s\tG0\t0\t%.6f\n" tag d cat name r.k r.fp
        r.noise.(0) r.noise.(1) r.noise.(2) silname.(sl) r.g0.(sl);
      for t = 0 to nt - 1 do
        List.iter
          (fun (gn, g) ->
            Printf.fprintf oc "%s\t%d\t%s\t%s\t%d\t%.6f\t%d\t%d\t%d\t%s\t%s\t%d\t%.6f\n" tag d cat name r.k r.fp
              r.noise.(0) r.noise.(1) r.noise.(2) silname.(sl) gn ss.(t) g.(sl).(t))
          [ "G1", r.g1; "G2", r.g2; "G3a", r.g3a; "G3b", r.g3b ]
      done
    done;
    flush oc;
    log "%-60s k=%4d G0 cl %.4f si %.4f" name r.k r.g0.(0) r.g0.(1) in
  let add cat name (p, k) = let r = score p k in add_res cat name r; r in
  (* CDC and the references of cs.ml, reproduced step for step *)
  let p0, k0 = canon (Array.map (fun r -> Hashtbl.find labs names_all.(r)) rows) in
  let sz0 = sizes p0 k0 and mem0 = members p0 k0 and cen0 = centroids p0 k0 in
  let r_cdc = add "ref" "CDC" (p0, k0) in
  (let tot = ref 0. in
   for i = 0 to n - 1 do
     let c = p0.(i) in
     if sz0.(c) >= 2 then begin
       let a = dist x.(i) cen0.(c) and b = ref infinity in
       for c' = 0 to k0 - 1 do if c' <> c then begin let v = dist x.(i) cen0.(c') in if v < !b then b := v end done;
       tot := !tot +. sil a !b
     end
   done;
   Printf.printf "# check: CDC simplified G0 from explicit centroids %.6f, from the distance matrix %.6f\n"
     (!tot /. float n) r_cdc.g0.(1));
  let st = Random.State.make [| 4242 |] in
  let big = ref 0 in
  Array.iteri (fun c m -> if Array.length m > Array.length mem0.(!big) then big := c) mem0;
  let p1 = Array.copy p0 in
  (let ms = mem0.(!big) in
   let far_from c =
     Array.fold_left (fun (bj, bv) j -> let v = dist x.(j) c in if v > bv then (j, v) else (bj, bv)) (-1, neg_infinity) ms
     |> fst in
   let s1 = far_from cen0.(!big) in
   let s2 = far_from x.(s1) in
   let m1 = Array.copy x.(s1) and m2 = Array.copy x.(s2) in
   for _ = 1 to 20 do
     let a1 = Array.make d 0. and a2 = Array.make d 0. and n1 = ref 0 and n2 = ref 0 in
     Array.iter
       (fun j ->
         if dist x.(j) m1 <= dist x.(j) m2 then begin
           p1.(j) <- !big; incr n1; for t = 0 to d - 1 do a1.(t) <- a1.(t) +. x.(j).(t) done
         end else begin
           p1.(j) <- k0; incr n2; for t = 0 to d - 1 do a2.(t) <- a2.(t) +. x.(j).(t) done
         end)
       ms;
     if !n1 > 0 && !n2 > 0 then
       for t = 0 to d - 1 do m1.(t) <- a1.(t) /. float !n1; m2.(t) <- a2.(t) /. float !n2 done
   done);
  ignore (add "crit" "CDC largest class split (2-means)" (canon p1));
  (let bc = ref (0, 1) and bv = ref infinity in
   for a = 0 to k0 - 1 do
     for b = a + 1 to k0 - 1 do
       let s = ref 0. in
       for t = 0 to d - 1 do let u = cen0.(a).(t) -. cen0.(b).(t) in s := !s +. (u *. u) done;
       if !s < !bv then begin bv := !s; bc := (a, b) end
     done
   done;
   let ca, cb = !bc in
   ignore (add "crit" "CDC two nearest classes merged" (canon (Array.map (fun c -> if c = cb then ca else c) p0))));
  ignore
    (add "crit" "CDC genogroups (label prefix)"
       (canon
          (Array.map
             (fun r ->
               let l = Hashtbl.find labs names_all.(r) in
               match String.index_opt l '.' with Some z -> String.sub l 0 z | None -> l)
             rows)));
  (let p3 = Array.copy p0 in
   for _ = 1 to n / 100 do
     let i = Random.State.int st n in
     let c = p0.(i) and best = ref (-1) and bvv = ref infinity in
     for c' = 0 to k0 - 1 do
       if c' <> c then begin let v = dist x.(i) cen0.(c') in if v < !bvv then begin bvv := v; best := c' end end
     done;
     p3.(i) <- !best
   done;
   ignore (add "pert" "CDC 1% moved to nearest other class" (canon p3)));
  (* References derived from CDC: the noise bin G2 pays for, one pass of reassignment, and the
     partitions G3a actually scores CDC as *)
  (let pa = Array.make_matrix nsil n 0. in
   ignore (base_score ~pts:pa p0 k0 false);
   Array.iteri
     (fun sl sn ->
       let p = Array.mapi (fun i c -> if pa.(sl).(i) < 0. then k0 + i else c) p0 in
       let neg = Array.fold_left (fun a v -> if v < 0. then a + 1 else a) 0 pa.(sl) in
       ignore (add "noise" (Printf.sprintf "CDC, its %d points of negative %s G0 as singletons" neg sn) (canon p)))
     [| "classical"; "simplified" |]);
  (let pc =
     Array.map
       (fun xi ->
         let best = ref 0 and bv = ref infinity in
         Array.iteri (fun c ce -> let v = dist xi ce in if v < !bv then begin bv := v; best := c end) cen0;
         !best)
       x in
   let moved = ref 0 in
   Array.iteri (fun i c -> if c <> p0.(i) then incr moved) pc;
   ignore (add "crit" (Printf.sprintf "CDC polished: %d points to nearest centroid" !moved) (canon pc)));
  (let pm = Array.make n 0 and acc = Array.make k0 0. in
   for i = 0 to n - 1 do
     Array.fill acc 0 k0 0.;
     for j = 0 to n - 1 do if j <> i then acc.(p0.(j)) <- acc.(p0.(j)) +. dd i j done;
     let best = ref (-1) and bv = ref infinity in
     for c = 0 to k0 - 1 do
       let m = if c = p0.(i) then sz0.(c) - 1 else sz0.(c) in
       if m > 0 then begin let v = acc.(c) /. float m in if v < !bv then begin bv := v; best := c end end
     done;
     pm.(i) <- !best
   done;
   let moved = ref 0 in
   Array.iteri (fun i c -> if c <> p0.(i) then incr moved) pm;
   ignore (add "crit" (Printf.sprintf "CDC polished: %d points to smallest mean distance" !moved) (canon pm)));
  (let _, q = score_q p0 k0 in
   Array.iteri
     (fun t s ->
       Array.iteri
         (fun sl sn ->
           if Array.length q.(sl).(t) > 0 then
             ignore (add "proj" (Printf.sprintf "G3a projection of CDC at s=%d (%s)" s sn) (canon q.(sl).(t))))
         [| "classical"; "simplified" |])
     ss);
  (* Outliers *)
  let g = Array.make d 0. in
  Array.iter (fun v -> for t = 0 to d - 1 do g.(t) <- g.(t) +. (v.(t) /. float n) done) x;
  let dg = Array.map (fun v -> dist v g) x in
  let by_far = Array.init n Fun.id in
  Array.stable_sort (fun a b -> compare dg.(b) dg.(a)) by_far;
  let far = by_far.(0) in
  let kk = 16 in
  let knn_i = Array.make_matrix n kk (-1) and knn_d = Array.make_matrix n kk infinity in
  for i = 0 to n - 1 do
    let bd = knn_d.(i) and bix = knn_i.(i) in
    for j = 0 to n - 1 do
      if j <> i then begin
        let t = dd i j in
        if t < bd.(kk - 1) then begin
          let pos = ref (kk - 1) in
          while !pos > 0 && bd.(!pos - 1) > t do
            bd.(!pos) <- bd.(!pos - 1); bix.(!pos) <- bix.(!pos - 1); decr pos
          done;
          bd.(!pos) <- t; bix.(!pos) <- j
        end
      end
    done
  done;
  log "neighbours done";
  let group_of seed size = Array.append [| seed |] (Array.sub knn_i.(seed) 0 (size - 1)) in
  let singletons_but groups =
    let p = Array.init n (fun i -> i + List.length groups) in
    List.iteri (fun gi gr -> Array.iter (fun i -> p.(i) <- gi) gr) groups;
    canon p in
  let one_vs_rest gr = let p = Array.make n 1 in Array.iter (fun i -> p.(i) <- 0) gr; canon p in
  List.iter
    (fun k ->
      let p = Array.make n 0 in
      for r = 0 to k - 1 do p.(by_far.(r)) <- r + 1 done;
      ignore (add "degen" (Printf.sprintf "all but farthest %d, as singletons" k) (canon p)))
    [ 1; 2; 3; 4; 6; 8; 16 ];
  List.iter
    (fun k ->
      ignore (add "degen" (Printf.sprintf "all but farthest %d, as one group" k) (one_vs_rest (Array.sub by_far 0 k))))
    [ 2; 3; 4; 6; 8; 16 ];
  List.iter
    (fun k ->
      ignore (add "degen" (Printf.sprintf "farthest point + %d NN as one group vs rest" (k - 1)) (one_vs_rest (group_of far k))))
    [ 2; 3; 4; 5; 6; 7; 8; 16 ];
  List.iter
    (fun k ->
      let bsi_ = ref None and bcl_ = ref None in
      for r = 0 to 31 do
        let res = score (fst (one_vs_rest (group_of by_far.(r) k))) 2 in
        (match !bsi_ with Some (_, rb) when rb.g0.(1) >= res.g0.(1) -> () | _ -> bsi_ := Some (r, res));
        (match !bcl_ with Some (_, rb) when rb.g0.(0) >= res.g0.(0) -> () | _ -> bcl_ := Some (r, res))
      done;
      let emit which = function
        | Some (r, res) ->
          add_res "degen" (Printf.sprintf "best seed+%d NN group, 32 farthest seeds, by %s (rank %d)" (k - 1) which r) res
        | None -> () in
      emit "si" !bsi_;
      emit "cl" !bcl_)
    [ 4; 8; 16 ];
  (* Peeled satellites *)
  let by_size = Array.init k0 Fun.id in
  Array.stable_sort (fun a b -> compare sz0.(b) sz0.(a)) by_size;
  let top5 = Array.sub by_size 0 (min 5 k0) in
  let peel keep pp =
    let p = Array.copy p0 and cnt = ref 0 in
    for c = 0 to k0 - 1 do
      if keep c then begin
        incr cnt;
        let ms = Array.copy mem0.(c) in
        let dc = Array.map (fun i -> i, dist x.(i) cen0.(c)) ms in
        Array.stable_sort (fun (_, u) (_, v) -> compare v u) dc;
        for r = 0 to pp - 1 do p.(fst dc.(r)) <- k0 + c done
      end
    done;
    canon p, !cnt in
  List.iter
    (fun pp ->
      let pk, _ = peel (fun c -> Array.mem c top5) pp in
      ignore (add "peel" (Printf.sprintf "peel %d farthest from each of the 5 largest classes" pp) pk);
      let pk, cnt = peel (fun c -> sz0.(c) >= 2 * pp) pp in
      ignore (add "peel" (Printf.sprintf "peel %d farthest from all %d classes >= %d" pp cnt (2 * pp)) pk))
    [ 2; 3; 4; 5; 8; 16 ];
  List.iter
    (fun pp ->
      let tight keep =
        let p = Array.copy p0 and cnt = ref 0 in
        for c = 0 to k0 - 1 do
          if keep c then begin
            incr cnt;
            let ms = mem0.(c) in
            let far_m = Array.fold_left (fun b i -> if dist x.(i) cen0.(c) > dist x.(b) cen0.(c) then i else b) ms.(0) ms in
            let byd = Array.copy ms in
            Array.stable_sort (fun u v -> compare (dd far_m u) (dd far_m v)) byd;
            for r = 0 to pp - 1 do p.(byd.(r)) <- k0 + c done
          end
        done;
        canon p, !cnt in
      let pk, _ = tight (fun c -> Array.mem c top5) in
      ignore (add "peel" (Printf.sprintf "peel farthest member + %d within-class NN, 5 largest classes" (pp - 1)) pk);
      let pk, cnt = tight (fun c -> sz0.(c) >= 2 * pp) in
      ignore (add "peel" (Printf.sprintf "peel farthest member + %d within-class NN, all %d classes >= %d" (pp - 1) cnt (2 * pp)) pk))
    [ 4; 8; 16 ];
  List.iter
    (fun (pp, t) ->
      let pk, cnt = peel (fun c -> sz0.(c) >= 2 * pp) pp in
      let _, q = score_q (fst pk) (snd pk) in
      Array.iteri
        (fun sl sn ->
          if Array.length q.(sl).(t) > 0 then
            ignore
              (add "proj" (Printf.sprintf "G3a projection at s=%d of peel %d from all %d classes (%s)" ss.(t) pp cnt sn)
                 (canon q.(sl).(t))))
        [| "classical"; "simplified" |])
    [ 2, 0; 3, 0; 4, 1; 5, 1; 8, 2 ];
  (* Atomised partitions: greedy matching by distance, the odd unit out joining its nearest match *)
  let decode bs m idx =
    let lo = ref 0 and hi = ref (m - 1) in
    while !lo < !hi do
      let mid = (!lo + !hi + 1) / 2 in
      if bs.(mid) <= idx then lo := mid else hi := mid - 1
    done;
    let u = !lo in
    u, idx - bs.(u) + u + 1 in
  let sorted_pairs =
    lazy
      (let a = Array.init np Fun.id in
       Array.stable_sort (fun u v -> Float.compare (Float.Array.unsafe_get dm u) (Float.Array.unsafe_get dm v)) a;
       log "sorted all pairs";
       a) in
  let attach_leftover units l =
    let best = ref 0 and bv = ref infinity in
    Array.iteri (fun ui un -> let v = dist x.(l) (unit_centroid un) in if v < !bv then begin bv := v; best := ui end) units;
    units.(!best) <- Array.append units.(!best) [| l |] in
  let greedy_pairs excluded =
    let matched = Array.copy excluded in
    let remaining = ref (Array.fold_left (fun acc e -> if e then acc else acc + 1) 0 excluded) in
    let units = ref [] and t = ref 0 in
    let sp = Lazy.force sorted_pairs in
    while !remaining > 1 do
      let u, v = decode base n sp.(!t) in
      incr t;
      if not matched.(u) && not matched.(v) then begin
        matched.(u) <- true; matched.(v) <- true; remaining := !remaining - 2;
        units := [| u; v |] :: !units
      end
    done;
    let units = Array.of_list (List.rev !units) in
    if !remaining = 1 then begin
      let l = ref (-1) in
      Array.iteri (fun i mt -> if not mt then l := i) matched;
      attach_leftover units !l
    end;
    units in
  let merge_level units =
    let m = Array.length units in
    let cu = Array.map unit_centroid units in
    let mp = m * (m - 1) / 2 in
    let bm = Array.init m (fun i -> (i * m) - (i * (i + 1) / 2)) in
    let du = Float.Array.make mp 0. in
    let idx = ref 0 in
    for u = 0 to m - 2 do for v = u + 1 to m - 1 do Float.Array.unsafe_set du !idx (dist cu.(u) cu.(v)); incr idx done done;
    let sp = Array.init mp Fun.id in
    Array.stable_sort (fun a b -> Float.compare (Float.Array.unsafe_get du a) (Float.Array.unsafe_get du b)) sp;
    let matched = Array.make m false and remaining = ref m and res = ref [] and t = ref 0 in
    while !remaining > 1 do
      let u, v = decode bm m sp.(!t) in
      incr t;
      if not matched.(u) && not matched.(v) then begin
        matched.(u) <- true; matched.(v) <- true; remaining := !remaining - 2;
        res := Array.append units.(u) units.(v) :: !res
      end
    done;
    let res = Array.of_list (List.rev !res) in
    if !remaining = 1 then begin
      let l = ref (-1) in
      Array.iteri (fun i mt -> if not mt then l := i) matched;
      let lc = cu.(!l) in
      let best = ref 0 and bv = ref infinity in
      Array.iteri (fun ui un -> let v = dist lc (unit_centroid un) in if v < !bv then begin bv := v; best := ui end) res;
      res.(!best) <- Array.append res.(!best) units.(!l)
    end;
    res in
  let of_units extra units =
    let p = Array.make n (-1) in
    List.iteri (fun gi gr -> Array.iter (fun i -> p.(i) <- gi) gr) extra;
    let off = List.length extra in
    Array.iteri (fun ui un -> Array.iter (fun i -> p.(i) <- off + ui) un) units;
    assert (Array.for_all (fun v -> v >= 0) p);
    canon p in
  let atom1 = greedy_pairs (Array.make n false) in
  ignore (add "degen" "greedy NN pairs (atomised, sizes 2..3)" (of_units [] atom1));
  let atom2 = merge_level atom1 in
  ignore (add "degen" "greedy pairs of pairs (atomised, sizes >= 4)" (of_units [] atom2));
  let atom3 = merge_level atom2 in
  ignore (add "degen" "greedy atomised, sizes >= 8" (of_units [] atom3));
  let atom4 = merge_level atom3 in
  ignore (add "degen" "greedy atomised, sizes >= 16" (of_units [] atom4));
  ignore (add "degen" "all singletons" (Array.init n Fun.id, n));
  let r1 = add "degen" "one cluster" (Array.make n 0, 1) in
  ignore r1;
  (* Gradient from the all-singleton start *)
  let closest = ref 0 in
  for i = 1 to n - 1 do if knn_d.(i).(0) < knn_d.(!closest).(0) then closest := i done;
  ignore (add "grad" "closest pair + singletons" (singletons_but [ group_of !closest 2 ]));
  Array.iter
    (fun s ->
      let by_rad = Array.init n Fun.id in
      Array.stable_sort (fun a b -> compare knn_d.(a).(s - 2) knn_d.(b).(s - 2)) by_rad;
      let g1_ = group_of by_rad.(0) s in
      let r = ref 1 in
      while Array.exists (fun i -> Array.mem i g1_) (group_of by_rad.(!r) s) do incr r done;
      let g2_ = group_of by_rad.(!r) s in
      ignore (add "grad" (Printf.sprintf "tightest %d-group + singletons" s) (singletons_but [ g1_ ]));
      ignore (add "grad" (Printf.sprintf "two tightest disjoint %d-groups + singletons" s) (singletons_but [ g1_; g2_ ]));
      let excl = Array.make n false in
      Array.iter (fun i -> excl.(i) <- true) g1_;
      ignore (add "grad" (Printf.sprintf "greedy NN pairs + tightest %d-group" s) (of_units [ g1_ ] (greedy_pairs excl))))
    ss;
  close_out oc;
  (* Tables *)
  let res = List.rev !results in
  let c0 = r_cdc in
  Printf.printf "#### %s d=%d n=%d classes=%d; CDC classes below s=4/8/16: %d/%d/%d (%d/%d/%d points)\n" tag d n k0
    (Array.fold_left (fun a m -> if m < 4 then a + 1 else a) 0 sz0)
    (Array.fold_left (fun a m -> if m < 8 then a + 1 else a) 0 sz0)
    (Array.fold_left (fun a m -> if m < 16 then a + 1 else a) 0 sz0)
    c0.noise.(0) c0.noise.(1) c0.noise.(2);
  Printf.printf "# one cluster: the library's full_silhouette returns -1 (no point has a b, so none is counted); here 0\n";
  Array.iteri
    (fun sl sn ->
      Printf.printf "## %s silhouette; '*' = above CDC under the same guard and s\n" sn;
      Printf.printf "  %-5s %-62s %5s %6s %8s | %-35s | %-35s | %-35s\n" "cat" "partition" "k" "f_p" "G0"
        "s=4:  G1      G2      G3a     G3b" "s=8:  G1      G2      G3a     G3b" "s=16: G1      G2      G3a     G3b";
      List.iter
        (fun (cat, name, r) ->
          let cell v cv = Printf.sprintf "%7.4f%s" v (if cat <> "ref" && v > cv +. 1e-9 then "*" else " ") in
          Printf.printf "  %-5s %-62s %5d %6.4f %s |" cat name r.k r.fp (cell r.g0.(sl) c0.g0.(sl));
          for t = 0 to nt - 1 do
            Printf.printf " %s%s%s%s |" (cell r.g1.(sl).(t) c0.g1.(sl).(t)) (cell r.g2.(sl).(t) c0.g2.(sl).(t))
              (cell r.g3a.(sl).(t) c0.g3a.(sl).(t)) (cell r.g3b.(sl).(t) c0.g3b.(sl).(t))
          done;
          print_newline ())
        res)
    [| "classical"; "simplified" |];
  log "all done"
