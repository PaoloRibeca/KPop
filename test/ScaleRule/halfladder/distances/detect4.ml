(* detect4 <twisted.txt> <inertia.txt> <labels> <tag> <rungs,> <resamples> <labelled-only 0|1> <flat|powers> <euclidean|manhattan|angle>
   detect3 with the metric and the distance as arguments, computed as KPop's Space.Distance does (dist.ml);
   under powers and euclidean it is detect3.
   detect3 --synthetic
   The calibrated valley detector of the scale-rule draft, revision 2, with the evidence it needs.
   Per rung, on ALL pairs of the chosen spectra:
   - 200 bins up to the 99th percentile (from a 10^6-bin histogram, not a sort), smoothed over +-1
     bin; the overflow bin counts towards all pairs but takes no part in smoothing or extrema;
   - extrema arranged peak-trough-...-peak, plateaus collapsed to their middle, so a trough needs
     a peak on both sides and the first and last bins are never troughs;
   - persistence simplification by relative size, smallest score first, until every score is at
     least 0.25.  A trough T between peaks L and R scores (min(L,R) - T) / min(L,R) and goes with
     the lower of L and R (the right one on a tie); an interior peak P between troughs A and B
     scores (P - max(A,B)) / min(L',R'), L' and R' being the peaks beyond A and B, and goes with
     the higher of A and B (the right one on a tie).  A bump inside a gap is therefore judged
     against the modes the gap separates, not against the gap's floor;
   - each surviving trough calibrated on its own: its prominence min(L,R) - T at fixed bins,
     recomputed in R resamples of the sequences by weighting every pair with the product of the
     two multiplicities; z = prominence / its resampling SD (infinite when the SD is 0);
     reappearance = share of resamples where its relative depth is still at least 0.25;
   - kept when z >= 3 and 0.2% <= share below <= 90%; kept valleys whose shares below are chained
     by gaps under 1.5 points form one group, which keeps its highest z (the lower share on a tie).
   The axis rules look only at kept valleys with at most 50% of pairs below:
   - count: on 5, 10, 20, 40, 80, 160, max, the rung with most valleys, ties to fewer axes;
   - doubling: on 1, 2, 4, ..., 128, max, walking up from the first rung with valleys, the last
     rung before one without any (and, for comparison, the last rung with valleys anywhere);
   - half powers of two: on 1, 2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64, 91, 128, 181, max -- two to
     the power k/2 rounded, the ladder the autotuner takes -- under the same rules as the doubling
     one, so that the two are read off the same detections.
   The level is the leftmost such valley at the chosen rung.  The same R resamples give each rule's
   pick and level in every resample, a valley being present in a resample when it reappears there.
   Also reports, at every rung, the label-based share of misordered comparisons between same-class
   and different-class pairs, which the rules never see.  Every candidate goes to stderr as a
   CAND line, for the threshold table *)
type valley = { k: int; bin: int; radius: float; share: float; z: float; reapp: float; depth: float }

type rung = { d: int; ladder: valley list; rds: float array array; mis: float }

let smooth_series c =
  let n = Array.length c in
  Array.init n (fun b ->
      let s = ref 0. and k = ref 0 in
      for j = b - 1 to b + 1 do
        if j >= 0 && j < n then begin s := !s +. c.(j); incr k end
      done;
      !s /. float !k)

(* Extrema as an alternating peak-trough-...-peak array of (is_peak, bin, value) *)
let extrema sm =
  let n = Array.length sm in
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
  Array.of_list (List.rev alt)

(* Removes adjacent peak-trough pairs, smallest relative score first, as described above *)
let simplify arr depth_min =
  let a = ref arr and go = ref true in
  while !go do
    let arr = !a in
    let len = Array.length arr in
    let v k = let _, _, x = arr.(k) in x in
    let best = ref infinity and lo = ref (-1) in
    let k = ref 1 in
    while !k < len - 1 do
      let p = Float.min (v (!k - 1)) (v (!k + 1)) in
      let s = if p > 0. then (p -. v !k) /. p else 0. in
      if s < !best then begin
        best := s;
        lo := if v (!k - 1) < v (!k + 1) then !k - 1 else !k
      end;
      k := !k + 2
    done;
    let k = ref 2 in
    while !k < len - 2 do
      let outer = Float.min (v (!k - 2)) (v (!k + 2)) in
      let s = if outer > 0. then (v !k -. Float.max (v (!k - 1)) (v (!k + 1))) /. outer else 0. in
      if s < !best then begin
        best := s;
        lo := if v (!k - 1) > v (!k + 1) then !k - 1 else !k
      end;
      k := !k + 2
    done;
    if !lo < 0 || !best >= depth_min then go := false
    else a := Array.append (Array.sub arr 0 !lo) (Array.sub arr (!lo + 2) (len - !lo - 2))
  done;
  !a

(* Troughs as (trough bin, left peak bin, right peak bin) *)
let troughs arr =
  List.init ((Array.length arr - 1) / 2) (fun i ->
      let k = (2 * i) + 1 in
      let _, tb, _ = arr.(k) and _, lb, _ = arr.(k - 1) and _, rb, _ = arr.(k + 1) in
      (tb, lb, rb))

let synthetic () =
  let bins = 200 in
  let g c h s b = h *. exp (-.((float b -. c) ** 2.) /. (2. *. s *. s)) in
  let base b = 20. +. g 40. 1000. 10. b +. g 140. 2000. 15. b in
  let bump h b = base b +. if b >= 78 && b <= 82 then h else 0. in
  let cases =
    [ "two modes (1000 at 40, 2000 at 140), floor 20", base;
      "bump +5 over bins 78-82", bump 5.;
      "bump +8 over bins 78-82", bump 8.;
      "bump +10 over bins 78-82", bump 10.;
      "bump +30 over bins 78-82", bump 30.;
      "bump +100 over bins 78-82", bump 100.;
      "third mode 100 at bin 85 (sd 5)", (fun b -> base b +. g 85. 100. 5. b);
      "third mode 400 at bin 85 (sd 5)", (fun b -> base b +. g 85. 400. 5. b);
      "ripple +-3 on every bin", (fun b -> base b +. (3. *. sin (float b *. 2.)));
      "10% dip at the top of the first mode", (fun b -> base b -. if b >= 39 && b <= 41 then 100. else 0.);
      "flat floor of 25 over bins 70-100", (fun b -> if b >= 70 && b <= 100 then 25. else base b);
      "rising from bin 0 (minimum at the edge)", (fun b -> 20. +. g 140. 2000. 15. b) ] in
  List.iter
    (fun (name, f) ->
      let sm = smooth_series (Array.init bins f) in
      let arr = extrema sm in
      let simp = simplify arr 0.25 in
      Printf.printf "%-48s | %2d extrema -> %d | %s\n" name (Array.length arr) (Array.length simp)
        (String.concat "; "
           (List.map
              (fun (tb, lb, rb) ->
                let p = Float.min sm.(lb) sm.(rb) in
                Printf.sprintf "trough at %d (peaks %d, %d), prominence %.1f, relative depth %.2f" tb lb rb
                  (p -. sm.(tb)) ((p -. sm.(tb)) /. p))
              (troughs simp))))
    cases

let main () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4)
  and rungs_req = String.split_on_char ',' Sys.argv.(5) |> List.map int_of_string
  and nres = int_of_string Sys.argv.(6) and filter = Sys.argv.(7) = "1" in
  let flat = Sys.argv.(8) = "flat" and kind = Dist.kind_of_string Sys.argv.(9) in
  let zmin = 3. and depth_min = 0.25 and share_lo = 0.002 and share_hi = 0.9 and ladder_hi = 0.5
  and dup = 0.015 in
  let labs = Kio.read_labels labf and iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> (not filter) || Hashtbl.mem labs names_all.(i))
    |> Array.of_list in
  let n = Array.length rows in
  let avail = min (Array.length iv) (Array.length coords_all.(rows.(0))) in
  let rungs = List.map (fun d -> min d avail) rungs_req |> List.sort_uniq compare in
  let dmax = List.fold_left max 0 rungs in
  let pc = Array.map (fun r -> Array.init dmax (fun k -> Dist.weigh kind ~flat iv coords_all.(r).(k) k)) rows in
  let lab =
    let tbl = Hashtbl.create 64 in
    Array.map
      (fun r ->
        match Hashtbl.find_opt labs names_all.(r) with
        | None -> -1
        | Some s -> (
          match Hashtbl.find_opt tbl s with
          | Some id -> id
          | None -> let id = Hashtbl.length tbl in Hashtbl.add tbl s id; id))
      rows in
  let np = n * (n - 1) / 2 and bins = 200 in
  Printf.printf "#### %s: %d spectra (%s), %d pairs, %d axes available, rungs %s, %d resamples\n%!" tag n
    (if filter then "labelled" else "all") np avail (String.concat "," (List.map string_of_int rungs)) nres;
  let acc = Float.Array.make np 0. and dist = Float.Array.make np 0. and norms = Array.make n 0. in
  let idx = Bigarray.Array1.create Bigarray.int16_unsigned Bigarray.c_layout np in
  let st = Random.State.make [| 20260913 |] in
  let weights =
    Array.init nres (fun _ ->
        let w = Array.make n 0 in
        for _ = 1 to n do let m = Random.State.int st n in w.(m) <- w.(m) + 1 done;
        w) in
  let fine = 1_000_000 in
  let fh = Float.Array.make fine 0. in
  let prev = ref 0 and results = ref [] in
  List.iter
    (fun d ->
      let t0 = Sys.time () in
      Dist.accumulate kind pc acc norms !prev d;
      Dist.to_dist kind n acc norms dist;
      prev := d;
      let maxsq = ref 0. in
      Float.Array.iter (fun v -> if v > !maxsq then maxsq := v) dist;
      Float.Array.fill fh 0 fine 0.;
      let fw = !maxsq *. (1. +. 1e-9) /. float fine in
      Float.Array.iter
        (fun v ->
          let fb = int_of_float (v /. fw) in
          let fb = if fb >= fine then fine - 1 else fb in
          Float.Array.unsafe_set fh fb (Float.Array.unsafe_get fh fb +. 1.))
        dist;
      let target = 0.99 *. float np and cum = ref 0. and fb = ref 0 in
      while !fb < fine - 1 && !cum +. Float.Array.get fh !fb <= target do
        cum := !cum +. Float.Array.get fh !fb;
        incr fb
      done;
      let hi = float (!fb + 1) *. fw in
      let w = hi /. float bins in
      let c = Array.make (bins + 1) 0. and same_c = Array.make (bins + 1) 0. and diff_c = Array.make (bins + 1) 0. in
      let id = ref 0 in
      for i = 0 to n - 2 do
        let li = lab.(i) in
        for j = i + 1 to n - 1 do
          let b = int_of_float (Float.Array.unsafe_get dist !id /. w) in
          let b = if b > bins then bins else b in
          Bigarray.Array1.unsafe_set idx !id b;
          c.(b) <- c.(b) +. 1.;
          let lj = lab.(j) in
          if li >= 0 && lj >= 0 then begin
            if li = lj then same_c.(b) <- same_c.(b) +. 1. else diff_c.(b) <- diff_c.(b) +. 1.
          end;
          incr id
        done
      done;
      let s_tot = Array.fold_left ( +. ) 0. same_c and d_tot = Array.fold_left ( +. ) 0. diff_c in
      let above = ref d_tot and num = ref 0. in
      for b = 0 to bins do
        above := !above -. diff_c.(b);
        num := !num +. (same_c.(b) *. (!above +. (0.5 *. diff_c.(b))))
      done;
      let misordered = if s_tot > 0. && d_tot > 0. then 1. -. (!num /. (s_tot *. d_tot)) else nan in
      let sm = smooth_series (Array.sub c 0 bins) in
      let cand = Array.of_list (troughs (simplify (extrema sm) depth_min)) in
      let nc = Array.length cand in
      let prom s (tb, lb, rb) = Float.min s.(lb) s.(rb) -. s.(tb) in
      let rdepth s (tb, lb, rb) = let p = Float.min s.(lb) s.(rb) in if p > 0. then (p -. s.(tb)) /. p else 0. in
      let proms = Array.make_matrix nc nres 0. and rds = Array.make_matrix nc nres 0. in
      let tot_full = float np in
      for r = 0 to nres - 1 do
        let wr = weights.(r) and acc = Array.make (bins + 1) 0. and id = ref 0 in
        for i = 0 to n - 2 do
          let wi = wr.(i) in
          if wi = 0 then id := !id + (n - 1 - i)
          else begin
            let fwi = float wi in
            for j = i + 1 to n - 1 do
              let wj = Array.unsafe_get wr j in
              if wj > 0 then begin
                let b = Bigarray.Array1.unsafe_get idx !id in
                Array.unsafe_set acc b (Array.unsafe_get acc b +. (fwi *. float wj))
              end;
              incr id
            done
          end
        done;
        let scale = tot_full /. Array.fold_left ( +. ) 0. acc in
        let smr = smooth_series (Array.init bins (fun b -> acc.(b) *. scale)) in
        Array.iteri (fun k cd -> proms.(k).(r) <- prom smr cd; rds.(k).(r) <- rdepth smr cd) cand
      done;
      let below = Array.make (bins + 1) 0. in
      for b = 1 to bins do below.(b) <- below.(b - 1) +. c.(b - 1) done;
      let all =
        Array.to_list
          (Array.mapi
             (fun k ((tb, _, _) as cd) ->
               let pf = prom sm cd in
               let mean = Array.fold_left ( +. ) 0. proms.(k) /. float nres in
               let var =
                 Array.fold_left (fun acc x -> acc +. ((x -. mean) *. (x -. mean))) 0. proms.(k) /. float (nres - 1) in
               let z = if var > 0. then pf /. sqrt var else infinity in
               let reapp =
                 float (Array.fold_left (fun acc x -> if x >= depth_min then acc + 1 else acc) 0 rds.(k)) /. float nres in
               { k; bin = tb; radius = (float tb +. 0.5) *. w;
                 share = (below.(tb) +. (0.5 *. c.(tb))) /. tot_full; z; reapp; depth = rdepth sm cd })
             cand) in
      let passing =
        List.filter (fun v -> v.z >= zmin && v.share >= share_lo && v.share <= share_hi) all
        |> List.sort (fun a b -> compare a.share b.share) in
      let groups =
        List.fold_left
          (fun acc v ->
            match acc with
            | ((p :: _) as g) :: rest when v.share -. p.share < dup -> (v :: g) :: rest
            | _ -> [ v ] :: acc)
          [] passing in
      let kept =
        List.rev_map
          (fun g ->
            match List.rev g with
            | first :: rest -> List.fold_left (fun b v -> if v.z > b.z then v else b) first rest
            | [] -> assert false)
          groups in
      let ladder = List.filter (fun v -> v.share <= ladder_hi) kept in
      let status v =
        if List.exists (fun w -> w.k = v.k) kept then "kept"
        else if v.z < zmin then "z"
        else if v.share < share_lo || v.share > share_hi then "share"
        else "duplicate" in
      List.iter
        (fun v ->
          Printf.eprintf "CAND\t%s\t%d\t%.6g\t%.5f\t%.3f\t%.3f\t%.3f\t%s\n" tag d v.radius v.share v.z v.reapp v.depth (status v))
        all;
      let fmt v = Printf.sprintf "%.1f%%(z%.1f,re%.0f%%,dep%.2f)" (100. *. v.share) v.z (100. *. v.reapp) v.depth in
      Printf.printf "  d=%-3d | %d candidates, %d kept | %s | score %d | misordered %.2f%% | %.0f s\n%!" d nc
        (List.length kept)
        (String.concat " "
           (List.map (fun v -> let s = status v in if s = "kept" then fmt v else Printf.sprintf "[%s %s]" s (fmt v)) all))
        (List.length ladder) (100. *. misordered) (Sys.time () -. t0);
      results := { d; ladder; rds; mis = misordered } :: !results)
    rungs;
  let res = List.rev !results in
  let on l = List.filter (fun r -> List.mem r.d l || r.d = avail) res in
  let count_ladder = on [ 5; 10; 20; 40; 80; 160 ]
  and doubling_ladder = on [ 1; 2; 4; 8; 16; 32; 64; 128 ]
  and half_ladder = on [ 1; 2; 3; 4; 6; 8; 11; 16; 23; 32; 45; 64; 91; 128; 181 ] in
  let leftmost = function
    | [] -> None
    | v :: rest -> Some (List.fold_left (fun b v -> if v.share < b.share then v else b) v rest) in
  (* Rungs come in increasing order, so a strict improvement keeps the fewer axes on a tie *)
  let pick_count ladder valleys_of =
    List.fold_left
      (fun (bd, bs) r -> let s = List.length (valleys_of r) in if s > 0 && s > bs then (Some r, s) else (bd, bs))
      (None, 0) ladder
    |> fst in
  (* Ties to MORE axes: a later rung with an equal count replaces the one held, which
     leaves the search the largest space that still resolves as much structure *)
  let pick_count_more ladder valleys_of =
    List.fold_left
      (fun (bd, bs) r -> let s = List.length (valleys_of r) in if s > 0 && s >= bs then (Some r, s) else (bd, bs))
      (None, 0) ladder
    |> fst in
  (* The FIRST maximal plateau, and its right end: the first rung reaching the maximum,
     extended only across rungs adjoining it at the same score, so that a later and
     unconnected rung of equal height does not win.  Within the plateau the rightmost
     rung is taken, leaving the search the largest space that resolves that structure *)
  let pick_first_plateau ladder valleys_of =
    let scored = List.map (fun r -> (r, List.length (valleys_of r))) ladder in
    let best = List.fold_left (fun acc (_, s) -> max acc s) 0 scored in
    if best = 0 then None
    else begin
      let rec find = function
        | [] -> None
        | (r, s) :: rest -> if s = best then Some (extend r rest) else find rest
      and extend r = function
        | (r2, s2) :: rest when s2 = best -> extend r2 rest
        | _ -> r in
      find scored
    end in
  let pick_doubling ladder valleys_of =
    let rec go last = function
      | [] -> last
      | r :: rest ->
        if valleys_of r <> [] then go (Some r) rest
        else (match last with Some _ -> last | None -> go None rest) in
    go None ladder in
  let pick_last ladder valleys_of =
    List.fold_left (fun acc r -> if valleys_of r <> [] then Some r else acc) None ladder in
  let summary l = String.concat " " (List.map (fun r -> Printf.sprintf "%d:%d" r.d (List.length r.ladder)) l) in
  let best l = List.fold_left (fun (bd, bm) r -> if r.mis < bm then (r.d, r.mis) else (bd, bm)) (-1, infinity) l in
  let bc, bcm = best count_ladder and bd, bdm = best doubling_ladder and bh, bhm = best half_ladder in
  Printf.printf "  COUNT LADDER scores %s | best rung by labels %d (%.2f%%)\n" (summary count_ladder) bc (100. *. bcm);
  Printf.printf "  DOUBLING LADDER scores %s | best rung by labels %d (%.2f%%)\n" (summary doubling_ladder) bd
    (100. *. bdm);
  Printf.printf "  HALF-POWERS LADDER scores %s | best rung by labels %d (%.2f%%)\n" (summary half_ladder) bh
    (100. *. bhm);
  let report name ladder rule =
    let pick = rule ladder (fun r -> r.ladder) in
    let tally = Hashtbl.create 8 and levels = ref [] in
    for i = 0 to nres - 1 do
      let present r = List.filter (fun v -> r.rds.(v.k).(i) >= depth_min) r.ladder in
      let key, lv = match rule ladder present with None -> (-1, None) | Some r -> (r.d, leftmost (present r)) in
      Hashtbl.replace tally key (1 + Option.value ~default:0 (Hashtbl.find_opt tally key));
      Option.iter (fun v -> levels := v.share :: !levels) lv
    done;
    let ts = Hashtbl.fold (fun k v acc -> (k, v) :: acc) tally [] |> List.sort (fun (_, a) (_, b) -> compare b a) in
    let picks =
      String.concat " "
        (List.map
           (fun (k, v) -> Printf.sprintf "%s:%.0f%%" (if k < 0 then "none" else string_of_int k) (100. *. float v /. float nres))
           ts) in
    let lvs = Array.of_list (List.sort compare !levels) in
    let nl = Array.length lvs in
    let q p = lvs.(int_of_float (p *. float (nl - 1))) in
    let level_res =
      if nl = 0 then "no level"
      else Printf.sprintf "level median %.1f%%, 5-95%% %.1f-%.1f%%" (100. *. q 0.5) (100. *. q 0.05) (100. *. q 0.95) in
    match pick with
    | None -> Printf.printf "  RULE %-34s pick none | resamples: picks %s; %s\n" name picks level_res
    | Some r ->
      let lv = Option.get (leftmost r.ladder) in
      Printf.printf "  RULE %-34s pick %3d (misordered %.2f%%) | leftmost %.1f%% (z %.1f, re %.0f%%) | resamples: picks %s; %s\n"
        name r.d (100. *. r.mis) (100. *. lv.share) lv.z (100. *. lv.reapp) picks level_res in
  report "count, ties to fewer axes" count_ladder pick_count;
  report "doubling, until valleys vanish" doubling_ladder pick_doubling;
  report "doubling, last rung with valleys" doubling_ladder pick_last;
  report "doubling ladder, most valleys, ties to more" doubling_ladder pick_count_more;
  report "doubling ladder, most valleys, ties to fewer" doubling_ladder pick_count;
  report "doubling ladder, most valleys, first plateau" doubling_ladder pick_first_plateau;
  report "half-powers ladder, most valleys, first plateau" half_ladder pick_first_plateau;
  report "half-powers ladder, most valleys, ties to fewer" half_ladder pick_count;
  print_endline ""

let () = if Array.length Sys.argv > 1 && Sys.argv.(1) = "--synthetic" then synthetic () else main ()
