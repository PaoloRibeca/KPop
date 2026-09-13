(* cs <twisted> <inertia> <labels> <d> <tag> <draws> <fv,fv,...>
   Scratch check of the draft's couples silhouette (m partners per cluster, shared by all scored
   points in a draw), of the push-back penalty on a few reference partitions, and of k-means++
   against uniform seeding in the split move (5 Lloyd iterations, as run_chain does) *)
let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3)
  and dreq = int_of_string Sys.argv.(4) and tag = Sys.argv.(5)
  and draws = int_of_string Sys.argv.(6)
  and fvs = String.split_on_char ',' Sys.argv.(7) |> List.map float_of_string in
  let labs = Kio.read_labels labf and iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs names_all.(i)) |> Array.of_list in
  let n = Array.length rows in
  let d = min dreq (Array.length iv) in
  let x = Array.map (fun r -> Array.init d (fun k -> coords_all.(r).(k) *. sqrt iv.(k))) rows in
  let np = n * (n - 1) / 2 in
  let base = Array.init n (fun i -> (i * n) - (i * (i + 1) / 2)) in
  let dm = Float.Array.make np 0. in
  let idx = ref 0 in
  for i = 0 to n - 2 do
    for j = i + 1 to n - 1 do
      let s = ref 0. in
      for k = 0 to d - 1 do let t = x.(i).(k) -. x.(j).(k) in s := !s +. (t *. t) done;
      Float.Array.unsafe_set dm !idx (sqrt !s); incr idx
    done
  done;
  let dd i j =
    if i = j then 0.
    else let i, j = if i < j then i, j else j, i in Float.Array.unsafe_get dm (base.(i) + (j - i - 1)) in
  let dist_to i c = let s = ref 0. in for k = 0 to d - 1 do let t = x.(i).(k) -. c.(k) in s := !s +. (t *. t) done; sqrt !s in
  let canon p =
    let h = Hashtbl.create 64 in
    let q = Array.map (fun l -> match Hashtbl.find_opt h l with Some v -> v | None -> let v = Hashtbl.length h in Hashtbl.add h l v; v) p in
    q, Hashtbl.length h in
  let members p k =
    let l = Array.make k [] in
    for i = n - 1 downto 0 do l.(p.(i)) <- i :: l.(p.(i)) done;
    Array.map Array.of_list l in
  let centroids p k =
    let c = Array.make_matrix k d 0. and s = Array.make k 0 in
    Array.iteri (fun i ci -> s.(ci) <- s.(ci) + 1; for t = 0 to d - 1 do c.(ci).(t) <- c.(ci).(t) +. x.(i).(t) done) p;
    Array.iteri (fun ci v -> if s.(ci) > 0 then for t = 0 to d - 1 do v.(t) <- v.(t) /. float s.(ci) done) c;
    c in
  let fp p k =
    let s = Array.make k 0 in
    Array.iter (fun c -> s.(c) <- s.(c) + 1) p;
    Array.fold_left (fun acc v -> acc +. (float v *. float (v - 1))) 0. s /. (float n *. float (n - 1)) in
  let exact_classical p k =
    let mem = members p k in
    let tot = ref 0. in
    for i = 0 to n - 1 do
      let c = p.(i) in
      if Array.length mem.(c) >= 2 then begin
        let a = ref 0. in
        Array.iter (fun j -> a := !a +. dd i j) mem.(c);
        let a = !a /. float (Array.length mem.(c) - 1) and b = ref infinity in
        for c' = 0 to k - 1 do
          if c' <> c && Array.length mem.(c') > 0 then begin
            let s = ref 0. in
            Array.iter (fun j -> s := !s +. dd i j) mem.(c');
            let v = !s /. float (Array.length mem.(c')) in
            if v < !b then b := v
          end
        done;
        if !b < infinity then tot := !tot +. ((!b -. a) /. Float.max a !b)
      end
    done;
    !tot /. float n in
  let exact_simplified p k =
    let mem = members p k and cen = centroids p k in
    let tot = ref 0. in
    for i = 0 to n - 1 do
      let c = p.(i) in
      if Array.length mem.(c) >= 2 then begin
        let a = dist_to i cen.(c) and b = ref infinity in
        for c' = 0 to k - 1 do
          if c' <> c && Array.length mem.(c') > 0 then begin
            let v = dist_to i cen.(c') in
            if v < !b then b := v
          end
        done;
        if !b < infinity then tot := !tot +. ((!b -. a) /. Float.max a !b)
      end
    done;
    !tot /. float n in
  let st = Random.State.make [| 4242 |] in
  let p0, k0 = canon (Array.map (fun r -> Hashtbl.find labs names_all.(r)) rows) in
  let mem0 = members p0 k0 and cen0 = centroids p0 k0 in
  (* P1: largest class split by 2-means *)
  let big = ref 0 in
  Array.iteri (fun c m -> if Array.length m > Array.length mem0.(!big) then big := c) mem0;
  let p1 = Array.copy p0 in
  let ms = mem0.(!big) in
  let far_from c = Array.fold_left (fun (bj, bv) j -> let v = dist_to j c in if v > bv then (j, v) else (bj, bv)) (-1, neg_infinity) ms |> fst in
  let s1 = far_from cen0.(!big) in
  let s2 = far_from x.(s1) in
  let m1 = Array.copy x.(s1) and m2 = Array.copy x.(s2) in
  for _ = 1 to 20 do
    let a1 = Array.make d 0. and a2 = Array.make d 0. and n1 = ref 0 and n2 = ref 0 in
    Array.iter (fun j ->
        if dist_to j m1 <= dist_to j m2 then begin p1.(j) <- !big; incr n1; for t = 0 to d - 1 do a1.(t) <- a1.(t) +. x.(j).(t) done end
        else begin p1.(j) <- k0; incr n2; for t = 0 to d - 1 do a2.(t) <- a2.(t) +. x.(j).(t) done end) ms;
    if !n1 > 0 && !n2 > 0 then for t = 0 to d - 1 do m1.(t) <- a1.(t) /. float !n1; m2.(t) <- a2.(t) /. float !n2 done
  done;
  let p1, k1 = canon p1 in
  (* P2: two classes with the nearest centroids merged *)
  let bc = ref (0, 1) and bv = ref infinity in
  for a = 0 to k0 - 1 do for b = a + 1 to k0 - 1 do
    let s = ref 0. in for t = 0 to d - 1 do let u = cen0.(a).(t) -. cen0.(b).(t) in s := !s +. (u *. u) done;
    if !s < !bv then begin bv := !s; bc := (a, b) end done done;
  let ca, cb = !bc in
  let p2, k2 = canon (Array.map (fun c -> if c = cb then ca else c) p0) in
  (* P3: 1% of points moved to their nearest other centroid *)
  let p3 = Array.copy p0 in
  for _ = 1 to n / 100 do
    let i = Random.State.int st n in
    let c = p0.(i) and best = ref (-1) and bvv = ref infinity in
    for c' = 0 to k0 - 1 do if c' <> c then begin let v = dist_to i cen0.(c') in if v < !bvv then begin bvv := v; best := c' end end done;
    p3.(i) <- !best
  done;
  let p3, k3 = canon p3 in
  (* P4: all but the point farthest from the global centroid *)
  let g = Array.make d 0. in
  Array.iter (fun v -> for t = 0 to d - 1 do g.(t) <- g.(t) +. (v.(t) /. float n) done) x;
  let far = ref 0 in
  for i = 1 to n - 1 do if dist_to i g > dist_to !far g then far := i done;
  let p4 = Array.init n (fun i -> if i = !far then 1 else 0) and k4 = 2 in
  let parts = [| "CDC", p0, k0; "CDC+largest split", p1, k1; "CDC+nearest merged", p2, k2;
                 "CDC+1% moved", p3, k3; "all but farthest point", p4, k4 |] in
  let np_ = Array.length parts in
  Printf.printf "#### %s d=%d n=%d classes=%d draws=%d\n%!" tag d n k0 draws;
  let ex = Array.map (fun (_, p, k) -> exact_classical p k) parts
  and es = Array.map (fun (_, p, k) -> exact_simplified p k) parts
  and fps = Array.map (fun (_, p, k) -> fp p k) parts in
  Array.iteri (fun t (nm, _, k) ->
      Printf.printf "  %-24s k=%3d f_p=%.4f classical %.4f simplified %.4f\n" nm k fps.(t) ex.(t) es.(t)) parts;
  List.iter (fun fv ->
      Printf.printf "  push-back with f_v=%.3f: simplified - lambda*|ln(f_p/f_v)| for lambda 0.05/0.1/0.2\n" fv;
      Array.iteri (fun t (nm, _, _) ->
          let pen = Float.abs (log (fps.(t) /. fv)) in
          Printf.printf "    %-24s penalty/lambda %.3f -> %.4f %.4f %.4f\n" nm pen
            (es.(t) -. (0.05 *. pen)) (es.(t) -. (0.1 *. pen)) (es.(t) -. (0.2 *. pen))) parts) fvs;
  (* Couples estimator *)
  let couples p k m prio scoring =
    let mem = members p k in
    let partners =
      Array.map (fun ms -> let a = Array.copy ms in Array.sort (fun u v -> compare prio.(u) prio.(v)) a;
                  Array.sub a 0 (min (m + 1) (Array.length a))) mem in
    let tot = ref 0. and cnt = ref 0 in
    Array.iter (fun i ->
        let c = p.(i) in
        incr cnt;
        if Array.length mem.(c) >= 2 then begin
          let s = ref 0. and q = ref 0 in
          Array.iter (fun j -> if j <> i && !q < m then begin s := !s +. dd i j; incr q end) partners.(c);
          let a = !s /. float !q and b = ref infinity in
          for c' = 0 to k - 1 do
            if c' <> c && Array.length partners.(c') > 0 then begin
              let pa = partners.(c') in
              let lim = min m (Array.length pa) in
              let s = ref 0. in
              for t = 0 to lim - 1 do s := !s +. dd i pa.(t) done;
              let v = !s /. float lim in
              if v < !b then b := v
            end
          done;
          if !b < infinity then tot := !tot +. ((!b -. a) /. Float.max a !b)
        end) scoring;
    !tot /. float !cnt in
  let mean_sd a =
    let n = float (Array.length a) in
    let m = Array.fold_left ( +. ) 0. a /. n in
    m, sqrt (Array.fold_left (fun acc v -> acc +. ((v -. m) *. (v -. m))) 0. a /. (n -. 1.)) in
  List.iter (fun (m, nsc) ->
      let est = Array.make_matrix draws np_ 0. and est_ind = Array.make_matrix draws np_ 0. in
      for r = 0 to draws - 1 do
        let prio = Array.init n (fun _ -> Random.State.float st 1.) in
        let scoring = if nsc >= n then Array.init n Fun.id else Array.init nsc (fun _ -> Random.State.int st n) in
        Array.iteri (fun t (_, p, k) -> est.(r).(t) <- couples p k m prio scoring) parts;
        (* independent draws per partition, for comparison *)
        Array.iteri (fun t (_, p, k) ->
            let prio = Array.init n (fun _ -> Random.State.float st 1.) in
            let scoring = if nsc >= n then Array.init n Fun.id else Array.init nsc (fun _ -> Random.State.int st n) in
            est_ind.(r).(t) <- couples p k m prio scoring) parts
      done;
      Printf.printf "  couples m=%d scored on %s:\n" m (if nsc >= n then "all points" else string_of_int nsc ^ " sampled points");
      Array.iteri (fun t (nm, _, _) ->
          let mu, sd = mean_sd (Array.map (fun e -> e.(t)) est) in
          let dmu, dsd = mean_sd (Array.map (fun e -> e.(0) -. e.(t)) est) in
          let _, dsdi = mean_sd (Array.map (fun e -> e.(0) -. e.(t)) est_ind) in
          Printf.printf "    %-24s mean %.4f (bias %+.4f vs classical) sd %.4f | CDC minus this: exact classical %+.4f, estimated %+.4f sd %.4f shared draw, sd %.4f independent draws\n"
            nm mu (mu -. ex.(t)) sd (ex.(0) -. ex.(t)) dmu dsd dsdi) parts)
    [ 5, 250; 5, n; 20, n ];
  (* Split seeding *)
  let wcss ms = let c = Array.make d 0. in
    Array.iter (fun j -> for t = 0 to d - 1 do c.(t) <- c.(t) +. (x.(j).(t) /. float (Array.length ms)) done) ms;
    Array.fold_left (fun acc j -> let v = dist_to j c in acc +. (v *. v)) 0. ms in
  let trial ms d2 =
    let sz = Array.length ms in
    let k1 = ms.(Random.State.int st sz) in
    let k2 =
      if not d2 then ms.(Random.State.int st sz)
      else begin
        let w = Array.map (fun j -> let v = dd k1 j in v *. v) ms in
        let tot = Array.fold_left ( +. ) 0. w in
        let u = Random.State.float st tot and acc = ref 0. and pick = ref ms.(sz - 1) and found = ref false in
        Array.iteri (fun t wv -> if not !found then begin acc := !acc +. wv; if !acc > u then begin pick := ms.(t); found := true end end) w;
        !pick
      end in
    if k1 = k2 then 0, 1.
    else begin
      let m1 = Array.copy x.(k1) and m2 = Array.copy x.(k2) and side = Array.make sz true in
      for _ = 1 to 5 do
        let n1 = ref 0 and n2 = ref 0 in
        Array.iteri (fun t j -> if dist_to j m1 <= dist_to j m2 then begin side.(t) <- true; incr n1 end else begin side.(t) <- false; incr n2 end) ms;
        if !n1 > 0 && !n2 > 0 then begin
          let a1 = Array.make d 0. and a2 = Array.make d 0. in
          Array.iteri (fun t j -> let a = if side.(t) then a1 else a2 in for u = 0 to d - 1 do a.(u) <- a.(u) +. x.(j).(u) done) ms;
          for u = 0 to d - 1 do m1.(u) <- a1.(u) /. float !n1; m2.(u) <- a2.(u) /. float !n2 done
        end
      done;
      let l1 = ref [] and l2 = ref [] in
      Array.iteri (fun t j -> if side.(t) then l1 := j :: !l1 else l2 := j :: !l2) ms;
      let a1 = Array.of_list !l1 and a2 = Array.of_list !l2 in
      let small = min (Array.length a1) (Array.length a2) in
      if small = 0 then 0, 1. else small, (wcss a1 +. wcss a2) /. wcss ms
    end in
  let trials = 100 in
  List.iter (fun d2 ->
      let noop = ref 0 and le2 = ref 0 and le5pct = ref 0 and tot = ref 0 and ratios = ref [] in
      let big5 = Array.to_list mem0 |> List.filter (fun m -> Array.length m >= 6)
                 |> List.sort (fun a b -> compare (Array.length b) (Array.length a)) in
      List.iter (fun ms ->
          for _ = 1 to trials do
            let small, ratio = trial ms d2 in
            incr tot;
            if small = 0 then incr noop
            else begin
              if small <= 2 then incr le2;
              if float small <= 0.05 *. float (Array.length ms) then incr le5pct;
              ratios := ratio :: !ratios
            end
          done) big5;
      let ra = Array.of_list !ratios in
      Array.sort compare ra;
      Printf.printf "  split seeding %-8s over %d classes x %d trials: no-op %.1f%%, smaller side <= 2 points %.1f%%, <= 5%% of class %.1f%%, median WCSS after/before %.3f\n"
        (if d2 then "k-means++" else "uniform") (List.length big5) trials
        (100. *. float !noop /. float !tot) (100. *. float !le2 /. float !tot)
        (100. *. float !le5pct /. float !tot)
        (if Array.length ra > 0 then ra.(Array.length ra / 2) else nan)) [ false; true ]
