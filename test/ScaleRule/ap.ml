(* ap <twisted> <inertia> <labels> <tag> <rungs,> <bphi> <bouter> <filter 0|1>
   Scratch check of the draft valley detector run on ALL pairs: counts rescaled to 100000 pairs,
   variance phi * Poisson, z >= 3.  Reports phi measured two ways (resamples drawing all pairs,
   and resamples drawing 100000 pairs as vboot did), its spread over groups of 10 resamples,
   whether +-1 smoothing and trough-peak differences obey the model, the valleys found, the
   ladder pick, and the pick under a bootstrap over sequences *)
let ref_pairs = 100000.

let smooth_series c =
  let n = Array.length c in
  Array.init n (fun b ->
      let s = ref 0. and k = ref 0 in
      for j = b - 1 to b + 1 do
        if j >= 0 && j < n then begin s := !s +. c.(j); incr k end
      done;
      !s /. float !k)

let detect c ~phi ~zmin =
  let n = Array.length c and w = 3. in
  let sm = smooth_series c in
  let ext = ref [] and b = ref 0 in
  while !b < n do
    let e = ref !b in
    while !e + 1 < n && sm.(!e + 1) = sm.(!b) do incr e done;
    let left = if !b > 0 then sm.(!b - 1) else neg_infinity
    and right = if !e + 1 < n then sm.(!e + 1) else neg_infinity in
    let mid = (!b + !e) / 2 in
    if sm.(!b) > left && sm.(!b) > right then ext := (true, mid, sm.(!b)) :: !ext
    else if sm.(!b) < left && sm.(!b) < right then ext := (false, mid, sm.(!b)) :: !ext;
    b := !e + 1
  done;
  let alt = ref [] in
  List.iter
    (fun ((is_peak, _, v) as x) ->
      match !alt with
      | [] -> if is_peak then alt := [ x ]
      | (lp, _, lv) :: rest ->
        if lp = is_peak then begin
          if (is_peak && v > lv) || ((not is_peak) && v < lv) then alt := x :: rest
        end else alt := x :: !alt)
    (List.rev !ext);
  let alt = match !alt with (false, _, _) :: rest -> rest | l -> l in
  let a = ref (Array.of_list (List.rev alt)) in
  let z arr k =
    let _, _, pl = arr.(k - 1) and _, _, t = arr.(k) and _, _, pr = arr.(k + 1) in
    let p = Float.min pl pr in
    (p -. t) /. sqrt (phi *. Float.max (p +. t) 1. /. w) in
  let go = ref true in
  while !go do
    let arr = !a in
    let worst = ref (-1) and wz = ref infinity and k = ref 1 in
    while !k < Array.length arr - 1 do
      let zk = z arr !k in
      if zk < !wz then begin wz := zk; worst := !k end;
      k := !k + 2
    done;
    if !worst < 0 || !wz >= zmin then go := false
    else begin
      let _, _, pl = arr.(!worst - 1) and _, _, pr = arr.(!worst + 1) in
      let drop = if pl < pr then !worst - 1 else !worst + 1 in
      let lo = min !worst drop in
      a := Array.append (Array.sub arr 0 lo) (Array.sub arr (lo + 2) (Array.length arr - lo - 2))
    end
  done;
  let arr = !a and out = ref [] and k = ref 1 in
  while !k < Array.length arr - 1 do
    let _, tb, _ = arr.(!k) and _, pbl, pl = arr.(!k - 1) and _, pbr, pr = arr.(!k + 1) in
    out := (tb, z arr !k, if pl <= pr then pbl else pbr) :: !out;
    k := !k + 2
  done;
  List.rev !out

(* Guards on share below; c and tot on the same scale *)
let guarded c tot vs =
  List.filter_map
    (fun (tb, z, pb) ->
      let below = ref 0. in
      for b = 0 to tb - 1 do below := !below +. c.(b) done;
      let s = !below /. tot in
      if s >= 0.002 && s <= 0.9 then Some (s, z, tb, pb) else None)
    vs

let agg c400 bins =
  let f = 400 / bins in
  let out = Array.make bins 0. in
  for b = 0 to 399 do out.(b / f) <- out.(b / f) +. c400.(b) done;
  out

let show vs =
  String.concat " " (List.map (fun (s, z, _, _) -> Printf.sprintf "%.1f%%(z%.1f)" (100. *. s) z) vs)

let stats a =
  let n = float (Array.length a) in
  let m = Array.fold_left ( +. ) 0. a /. n in
  let v = Array.fold_left (fun acc x -> acc +. ((x -. m) *. (x -. m))) 0. a /. (n -. 1.) in
  m, sqrt v, Array.fold_left Float.min infinity a, Array.fold_left Float.max neg_infinity a

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4)
  and rungs_req = String.split_on_char ',' Sys.argv.(5) |> List.map int_of_string
  and bphi = int_of_string Sys.argv.(6) and bouter = int_of_string Sys.argv.(7)
  and filter = Sys.argv.(8) = "1" in
  let labs = Kio.read_labels labf and iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> (not filter) || Hashtbl.mem labs names_all.(i))
    |> Array.of_list in
  let n = Array.length rows in
  let avail = min (Array.length iv) (Array.length coords_all.(rows.(0))) in
  let rungs = List.map (fun d -> min d avail) rungs_req |> List.sort_uniq compare in
  let nr = List.length rungs in
  let dmax = List.fold_left max 0 rungs in
  let pc = Array.map (fun r -> Array.init dmax (fun k -> coords_all.(r).(k) *. sqrt iv.(k))) rows in
  let np = n * (n - 1) / 2 in
  Printf.printf "#### %s: n=%d pairs=%d axes available=%d rungs=%s bphi=%d bouter=%d\n%!" tag n np
    avail (String.concat "," (List.map string_of_int rungs)) bphi bouter;
  let base = Array.init n (fun i -> (i * n) - (i * (i + 1) / 2)) in
  let sq = Float.Array.make np 0. in
  let idx400 = Bigarray.Array1.create Bigarray.int16_unsigned Bigarray.c_layout np in
  let fine = 1_000_000 in
  let fh = Float.Array.make fine 0. in
  let iter_pairs_w w f =
    let idx = ref 0 in
    for i = 0 to n - 2 do
      let wi = w.(i) in
      if wi = 0 then idx := !idx + (n - 1 - i)
      else begin
        let fwi = float wi in
        for j = i + 1 to n - 1 do
          let wj = w.(j) in
          if wj > 0 then f !idx (fwi *. float wj);
          incr idx
        done
      end
    done in
  let ones = Array.make n 1 in
  let q99 w maxd =
    Float.Array.fill fh 0 fine 0.;
    let fw = maxd *. (1. +. 1e-9) /. float fine in
    let tot = ref 0. in
    iter_pairs_w w (fun id x ->
        let fb = int_of_float (sqrt (Float.Array.unsafe_get sq id) /. fw) in
        let fb = if fb >= fine then fine - 1 else fb in
        Float.Array.unsafe_set fh fb (Float.Array.unsafe_get fh fb +. x);
        tot := !tot +. x);
    let target = Float.of_int (int_of_float (0.99 *. !tot)) in
    let cum = ref 0. and fb = ref 0 in
    while !fb < fine - 1 && !cum +. Float.Array.get fh !fb <= target do
      cum := !cum +. Float.Array.get fh !fb; incr fb
    done;
    float (!fb + 1) *. fw, !tot in
  let ngroups = bphi / 10 in
  let mk () = Array.make nr 0 in
  let res_b10 = mk () and res_b200 = mk () and res_a10 = mk () and res_a200 = mk ()
  and res_raw = mk () and res_b100bins = mk () and res_b400bins = mk () and res_p5 = mk () in
  let res_groups = Array.make_matrix nr ngroups 0 in
  let res_outer = Array.make_matrix nr bouter 0 and prom_outer = Array.make_matrix nr bouter nan in
  let prev = ref 0 in
  List.iteri
    (fun ri d ->
      let t0 = Sys.time () in
      let idx = ref 0 in
      for i = 0 to n - 2 do
        let a = pc.(i) in
        for j = i + 1 to n - 1 do
          let b = pc.(j) in
          let s = ref (Float.Array.unsafe_get sq !idx) in
          for k = !prev to d - 1 do
            let x = Array.unsafe_get a k -. Array.unsafe_get b k in
            s := !s +. (x *. x)
          done;
          Float.Array.unsafe_set sq !idx !s;
          incr idx
        done
      done;
      prev := d;
      let maxd = ref 0. in
      for id = 0 to np - 1 do
        let v = Float.Array.unsafe_get sq id in
        if v > !maxd then maxd := v
      done;
      let maxd = sqrt !maxd in
      let hi0, _ = q99 ones maxd in
      let w400 = hi0 /. 400. in
      let c400 = Array.make 401 0. in
      for id = 0 to np - 1 do
        let b = int_of_float (sqrt (Float.Array.unsafe_get sq id) /. w400) in
        let b = if b > 400 then 400 else b in
        Bigarray.Array1.unsafe_set idx400 id b;
        c400.(b) <- c400.(b) +. 1.
      done;
      (* Resamples of the sequences, fixed bins *)
      let st = Random.State.make [| 20260913; d |] in
      let all_c = Array.make_matrix bphi 401 0. and k_c = Array.make_matrix bphi 401 0. in
      let all_tot = Array.make bphi 0. and k_tot = Array.make bphi 0. in
      for r = 0 to bphi - 1 do
        let member = Array.init n (fun _ -> Random.State.int st n) in
        let w = Array.make n 0 in
        Array.iter (fun m -> w.(m) <- w.(m) + 1) member;
        let acc = all_c.(r) in
        let idx = ref 0 in
        for i = 0 to n - 2 do
          let wi = w.(i) in
          if wi = 0 then idx := !idx + (n - 1 - i)
          else begin
            let fwi = float wi in
            for j = i + 1 to n - 1 do
              let wj = Array.unsafe_get w j in
              if wj > 0 then begin
                let b = Bigarray.Array1.unsafe_get idx400 !idx in
                Array.unsafe_set acc b (Array.unsafe_get acc b +. (fwi *. float wj))
              end;
              incr idx
            done
          end
        done;
        all_tot.(r) <- Array.fold_left ( +. ) 0. acc;
        let acc = k_c.(r) and m = ref 0 in
        for _ = 1 to 100000 do
          let p = Random.State.int st n and q = Random.State.int st n in
          let a = member.(p) and b = member.(q) in
          if p <> q && a <> b then begin
            let i = min a b and j = max a b in
            let bb = Bigarray.Array1.unsafe_get idx400 (base.(i) + (j - i - 1)) in
            acc.(bb) <- acc.(bb) +. 1.;
            incr m
          end
        done;
        k_tot.(r) <- float !m
      done;
      let scaled c tot bins = Array.map (fun x -> x *. ref_pairs /. tot) (agg c bins) in
      let range a b = List.init (b - a) (fun i -> a + i) in
      let phi_of mats tots reps bins smoothed =
        let series =
          List.map
            (fun r ->
              let s = scaled mats.(r) tots.(r) bins in
              if smoothed then smooth_series s else s)
            reps
          |> Array.of_list in
        let nrep = float (Array.length series) in
        let num = ref 0. and den = ref 0. in
        for b = 1 to bins - 2 do
          let mean = ref 0. in
          Array.iter (fun s -> mean := !mean +. s.(b)) series;
          let mean = !mean /. nrep in
          if mean >= 50. *. 200. /. float bins then begin
            let var = ref 0. in
            Array.iter (fun s -> let x = s.(b) -. mean in var := !var +. (x *. x)) series;
            num := !num +. (!var /. (nrep -. 1.));
            den := !den +. mean
          end
        done;
        if smoothed then !num /. (!den /. 3.) else !num /. !den in
      let allr = range 0 bphi in
      let phi_all200 = phi_of all_c all_tot allr 200 false
      and phi_k200 = phi_of k_c k_tot allr 200 false
      and phi_all100 = phi_of all_c all_tot allr 100 false
      and phi_all400 = phi_of all_c all_tot allr 400 false
      and phis_all = phi_of all_c all_tot allr 200 true
      and phis_k = phi_of k_c k_tot allr 200 true in
      let g_all = Array.init ngroups (fun g -> phi_of all_c all_tot (range (10 * g) ((10 * g) + 10)) 200 false)
      and g_k = Array.init ngroups (fun g -> phi_of k_c k_tot (range (10 * g) ((10 * g) + 10)) 200 false) in
      let orig bins = scaled c400 (float np) bins in
      let det bins phi zmin = let c = orig bins in guarded c ref_pairs (detect c ~phi ~zmin) in
      let v_b10 = det 200 g_all.(0) 3. and v_b200 = det 200 phi_all200 3.
      and v_a10 = det 200 g_k.(0) 3. and v_a200 = det 200 phi_k200 3.
      and v_p5 = det 200 1. 5.
      and v_b100 = det 100 phi_all100 3. and v_b400 = det 400 phi_all400 3. in
      let raw = agg c400 200 in
      let v_raw = guarded raw (float np) (detect raw ~phi:1. ~zmin:5.) in
      res_b10.(ri) <- List.length v_b10; res_b200.(ri) <- List.length v_b200;
      res_a10.(ri) <- List.length v_a10; res_a200.(ri) <- List.length v_a200;
      res_raw.(ri) <- List.length v_raw; res_b100bins.(ri) <- List.length v_b100;
      res_b400bins.(ri) <- List.length v_b400; res_p5.(ri) <- List.length v_p5;
      Array.iteri (fun g ph -> res_groups.(ri).(g) <- List.length (det 200 ph 3.)) g_all;
      let m1, s1, lo1, hi1 = stats g_all and m2, s2, lo2, hi2 = stats g_k in
      Printf.printf "== %s d=%d hi0=%.4g\n" tag d hi0;
      Printf.printf "  phi from %d resamples, 200 bins: all-pairs resamples %.2f | 100k-pair resamples %.2f | all-pairs at 100/400 bins %.2f/%.2f\n"
        bphi phi_all200 phi_k200 phi_all100 phi_all400;
      Printf.printf "  phi from groups of 10 (%d groups): all-pairs mean %.2f sd %.2f range %.2f-%.2f | 100k mean %.2f sd %.2f range %.2f-%.2f\n"
        ngroups m1 s1 lo1 hi1 m2 s2 lo2 hi2;
      Printf.printf "  inflation of +-1-smoothed bins over mean/3: all-pairs %.2f (x%.2f of raw phi) | 100k %.2f (x%.2f)\n"
        phis_all (phis_all /. phi_all200) phis_k (phis_k /. phi_k200);
      Printf.printf "  raw Poisson z>=5 on all %d pairs: %d valleys: %s\n" np (List.length v_raw) (show v_raw);
      Printf.printf "  Poisson z>=5 after rescaling all pairs to 100k (phi=1): %d: %s\n" (List.length v_p5) (show v_p5);
      Printf.printf "  draft A (phi from 10 100k-pair resamples = %.2f, z>=3): %d: %s\n" g_k.(0) (List.length v_a10) (show v_a10);
      Printf.printf "  draft A (phi_200 = %.2f): %d: %s\n" phi_k200 (List.length v_a200) (show v_a200);
      Printf.printf "  draft B (phi from 10 all-pair resamples = %.2f, z>=3): %d: %s\n" g_all.(0) (List.length v_b10) (show v_b10);
      Printf.printf "  draft B (phi_200 = %.2f): %d: %s\n" phi_all200 (List.length v_b200) (show v_b200);
      Printf.printf "  draft B at 100 bins: %d: %s\n  draft B at 400 bins: %d: %s\n" (List.length v_b100) (show v_b100)
        (List.length v_b400) (show v_b400);
      Printf.printf "  draft B valley count under each group-of-10 phi: %s\n"
        (String.concat " " (Array.to_list (Array.map string_of_int res_groups.(ri))));
      (* Does the prominence obey the model? *)
      let sm_all = Array.init bphi (fun r -> smooth_series (scaled all_c.(r) all_tot.(r) 200))
      and sm_k = Array.init bphi (fun r -> smooth_series (scaled k_c.(r) k_tot.(r) 200)) in
      let so = smooth_series (orig 200) in
      List.iter
        (fun (s, z, tb, pb) ->
          let vr sms =
            let ds = Array.map (fun a -> a.(pb) -. a.(tb)) sms
            and ms = Array.map (fun a -> a.(pb) +. a.(tb)) sms in
            let _, sd, _, _ = stats ds and mm, _, _, _ = stats ms in
            sd *. sd, mm in
          let va, ma = vr sm_all and vk, mk = vr sm_k in
          Printf.printf
            "   valley %.1f%% z %.1f (trough bin %d, lower flank bin %d, prominence %.1f of peak %.1f): var(prominence)/model: all-pairs %.2f | 100k %.2f (vs pure Poisson %.2f)\n"
            (100. *. s) z tb pb (so.(pb) -. so.(tb)) so.(pb)
            (va /. (phi_all200 *. ma /. 3.)) (vk /. (phi_k200 *. mk /. 3.)) (vk /. (mk /. 3.)))
        v_b200;
      (* Bootstrap over sequences: whole detector per resample, phi fixed at phi_all200 *)
      for r = 0 to bouter - 1 do
        let st = Random.State.make [| 777; r |] in
        let w = Array.make n 0 in
        for _ = 1 to n do let m = Random.State.int st n in w.(m) <- w.(m) + 1 done;
        let hi, tot = q99 w maxd in
        let wb = hi /. 200. in
        let c = Array.make 200 0. in
        iter_pairs_w w (fun id x ->
            let b = int_of_float (sqrt (Float.Array.unsafe_get sq id) /. wb) in
            if b < 200 then c.(b) <- c.(b) +. x);
        let cs = Array.map (fun x -> x *. ref_pairs /. tot) c in
        let vs = guarded cs ref_pairs (detect cs ~phi:phi_all200 ~zmin:3.) in
        res_outer.(ri).(r) <- List.length vs;
        match vs with
        | [] -> ()
        | h :: _ ->
          let s, _, _, _ =
            List.fold_left (fun ((_, zb, _, _) as best) ((_, z, _, _) as v) -> if z > zb then v else best) h vs in
          prom_outer.(ri).(r) <- s
      done;
      let hist = Hashtbl.create 8 in
      Array.iter (fun c -> Hashtbl.replace hist c (1 + try Hashtbl.find hist c with Not_found -> 0)) res_outer.(ri);
      let keys = Hashtbl.fold (fun k _ acc -> k :: acc) hist [] |> List.sort compare in
      Printf.printf "  bootstrap (%d resamples) valley count distribution: %s\n" bouter
        (String.concat " " (List.map (fun k -> Printf.sprintf "%d:%d" k (Hashtbl.find hist k)) keys));
      let ph = Hashtbl.create 8 in
      Array.iter
        (fun s ->
          let k = if Float.is_nan s then "none" else Printf.sprintf "%.0f" (100. *. s) in
          Hashtbl.replace ph k (1 + try Hashtbl.find ph k with Not_found -> 0))
        prom_outer.(ri);
      let pk = Hashtbl.fold (fun k v acc -> (k, v) :: acc) ph [] |> List.sort compare in
      Printf.printf "  bootstrap most prominent valley, share below rounded to 1%%: %s\n"
        (String.concat " " (List.map (fun (k, v) -> Printf.sprintf "%s%%:%d" k v) pk));
      Printf.printf "  [rung time %.0f s]\n%!" (Sys.time () -. t0))
    rungs;
  let ladder l = List.map (fun d -> min d avail) l |> List.sort_uniq compare in
  let ladders = [ "evidence 5..80,120,all", ladder [ 5; 10; 20; 40; 80; 120; avail ];
                  "draft 5..80,160,all", ladder [ 5; 10; 20; 40; 80; 160; avail ] ] in
  let index d = let rec go i = function [] -> -1 | x :: r -> if x = d then i else go (i + 1) r in go 0 rungs in
  let pick lad count_of =
    let best = ref 0 and bd = ref (-1) in
    List.iter (fun d -> let c = count_of (index d) in if c > 0 && c >= !best then begin best := c; bd := d end) lad;
    !bd in
  List.iter
    (fun (lname, lad) ->
      let row name arr =
        Printf.printf "  %-28s counts %s | pick %d\n" name
          (String.concat " " (List.map (fun d -> string_of_int arr.(index d)) lad)) (pick lad (fun i -> arr.(i))) in
      Printf.printf "-- %s ladder %s (%s)\n" tag lname (String.concat "," (List.map string_of_int lad));
      row "raw Poisson z5, all pairs" res_raw;
      row "Poisson z5 rescaled 100k" res_p5;
      row "draft A phi10" res_a10;
      row "draft A phi200" res_a200;
      row "draft B phi10" res_b10;
      row "draft B phi200" res_b200;
      row "draft B phi200 100 bins" res_b100bins;
      row "draft B phi200 400 bins" res_b400bins;
      let gp = Hashtbl.create 8 in
      for g = 0 to ngroups - 1 do
        let p = pick lad (fun i -> res_groups.(i).(g)) in
        Hashtbl.replace gp p (1 + try Hashtbl.find gp p with Not_found -> 0)
      done;
      Printf.printf "  pick under each group-of-10 phi: %s\n"
        (String.concat " " (Hashtbl.fold (fun k v acc -> Printf.sprintf "%d:%d" k v :: acc) gp [] |> List.sort compare));
      let op = Hashtbl.create 8 in
      for r = 0 to bouter - 1 do
        let p = pick lad (fun i -> res_outer.(i).(r)) in
        Hashtbl.replace op p (1 + try Hashtbl.find op p with Not_found -> 0)
      done;
      Printf.printf "  pick under bootstrap over sequences: %s\n%!"
        (String.concat " " (Hashtbl.fold (fun k v acc -> Printf.sprintf "%d:%d" k v :: acc) op [] |> List.sort compare)))
    ladders
