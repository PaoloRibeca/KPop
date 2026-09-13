(* vtree <twisted.txt> <inertia.txt> <tag> <d,d,...> <level,names> <cutspec> <labels> [<labels>]
   The tree's view of the pair-distance histogram.  For each d, the average-linkage (UPGMA)
   dendrogram of the labelled spectra over the first d principal axes, built by the nearest-
   neighbour chain; the histogram of the heights at which it joins pairs (each join of clusters A
   and B counting |A||B| pairs), split by whether the pairs share a class at each label level, on
   the same linear and log bins as vdump and scaled to the same number of sampled pairs so that
   the two histograms can be drawn on one axis; its valleys; and the tree cut at every valley --
   of the pair-distance histogram (cutspec "d=r,r;d=r", "-" for none) and of its own -- scored
   against each label level, together with the best cut for each level *)
let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and tag = Sys.argv.(3)
  and ds = List.map int_of_string (String.split_on_char ',' Sys.argv.(4))
  and lnames = String.split_on_char ',' Sys.argv.(5) and cutspec = Sys.argv.(6) in
  let nlab = Array.length Sys.argv - 7 in
  let labs = Array.init nlab (fun i -> Kio.read_labels Sys.argv.(7 + i)) in
  let pair_cuts d =
    if cutspec = "-" then []
    else
      String.split_on_char ';' cutspec
      |> List.concat_map (fun s ->
             match String.split_on_char '=' s with
             | [ dd; rs ] when int_of_string dd = d && rs <> "-" && rs <> "" ->
               List.map float_of_string (String.split_on_char ',' rs)
             | _ -> []) in
  let iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs.(0) names_all.(i)) |> Array.of_list in
  let n = Array.length rows and maxd = Array.length iv in
  let ncl = Array.make nlab 0 in
  let lab =
    Array.mapi
      (fun lv h ->
        let tbl = Hashtbl.create 64 in
        let ids =
          Array.map
            (fun r ->
              let s = Hashtbl.find h names_all.(r) in
              match Hashtbl.find_opt tbl s with
              | Some id -> id
              | None -> let id = Hashtbl.length tbl in Hashtbl.add tbl s id; id)
            rows in
        ncl.(lv) <- Hashtbl.length tbl;
        ids)
      labs in
  let npairs = n * (n - 1) / 2 in
  let total_same =
    Array.init nlab (fun lv ->
        let c = Array.make ncl.(lv) 0 in
        Array.iter (fun l -> c.(l) <- c.(l) + 1) lab.(lv);
        Array.fold_left (fun acc x -> acc +. (float x *. float (x - 1) /. 2.)) 0. c) in
  let idx i j = let i, j = if i < j then i, j else j, i in (i * (2 * n - i - 1) / 2) + (j - i - 1) in
  (* The detector's pairs, for the histogram range and the scale of the counts *)
  let state = Random.State.make [| 17 |] in
  let pa = Array.make 100000 0 and pb = Array.make 100000 0 and cnt = ref 0 in
  for _ = 1 to 100000 do
    let a = Random.State.int state n and b = Random.State.int state n in
    if a <> b then begin pa.(!cnt) <- a; pb.(!cnt) <- b; incr cnt end
  done;
  let m = !cnt and bins = 200 in
  let sq = Float.Array.make npairs 0. and dm = Float.Array.make npairs 0. and done_axes = ref 0 in
  let jarr f a = "[" ^ String.concat "," (Array.to_list (Array.map f a)) ^ "]" in
  List.iter
    (fun d ->
      let d = min d maxd in
      for k = !done_axes to d - 1 do
        let wk = iv.(k) in
        let col = Array.map (fun r -> coords_all.(r).(k)) rows in
        for i = 0 to n - 2 do
          let ci = col.(i) and base = (i * (2 * n - i - 1) / 2) - i - 1 in
          for j = i + 1 to n - 1 do
            let x = ci -. col.(j) in
            Float.Array.unsafe_set sq (base + j) (Float.Array.unsafe_get sq (base + j) +. (wk *. x *. x))
          done
        done
      done;
      done_axes := d;
      for p = 0 to npairs - 1 do Float.Array.unsafe_set dm p (sqrt (Float.Array.unsafe_get sq p)) done;
      let samp = Array.init m (fun p -> Float.Array.get dm (idx pa.(p) pb.(p))) in
      Array.sort compare samp;
      let zeros = ref 0 in
      Array.iter (fun x -> if x <= 0. then incr zeros) samp;
      let hi = samp.(int_of_float (float m *. 0.99)) in
      let lo = samp.(!zeros + int_of_float (float (m - !zeros) *. 0.001)) in
      (* Average linkage by the nearest-neighbour chain, with class counts carried per cluster *)
      let active = Array.make n true and size = Array.make n 1 in
      let cc = Array.init nlab (fun lv -> Array.init n (fun i -> let a = Array.make ncl.(lv) 0 in a.(lab.(lv).(i)) <- 1; a)) in
      let chain = Array.make (n + 1) 0 and top = ref 0 and start = ref 0 and nm = ref 0 in
      let mi = Array.make (n - 1) 0 and mj = Array.make (n - 1) 0 and mh = Array.make (n - 1) 0.
      and mp = Array.make (n - 1) 0. and msame = Array.init nlab (fun _ -> Array.make (n - 1) 0.) in
      while !nm < n - 1 do
        if !top = 0 then begin
          while not active.(!start) do incr start done;
          chain.(0) <- !start; top := 1
        end;
        let a = chain.(!top - 1) in
        let prev = if !top >= 2 then chain.(!top - 2) else -1 in
        let best = ref prev and bd = ref (if prev >= 0 then Float.Array.get dm (idx a prev) else infinity) in
        for k = 0 to n - 1 do
          if active.(k) && k <> a then begin
            let dk = Float.Array.get dm (idx a k) in
            if dk < !bd then begin bd := dk; best := k end
          end
        done;
        if prev >= 0 && !best = prev then begin
          top := !top - 2;
          let i = min a prev and j = max a prev in
          let si = float size.(i) and sj = float size.(j) in
          mi.(!nm) <- i; mj.(!nm) <- j; mh.(!nm) <- !bd; mp.(!nm) <- si *. sj;
          for lv = 0 to nlab - 1 do
            let ci = cc.(lv).(i) and cj = cc.(lv).(j) in
            let s = ref 0. in
            for l = 0 to ncl.(lv) - 1 do
              s := !s +. (float ci.(l) *. float cj.(l));
              ci.(l) <- ci.(l) + cj.(l)
            done;
            msame.(lv).(!nm) <- !s
          done;
          incr nm;
          for k = 0 to n - 1 do
            if active.(k) && k <> i && k <> j then
              Float.Array.set dm (idx i k)
                (((si *. Float.Array.get dm (idx i k)) +. (sj *. Float.Array.get dm (idx j k))) /. (si +. sj))
          done;
          size.(i) <- size.(i) + size.(j);
          active.(j) <- false
        end else begin chain.(!top) <- !best; incr top end
      done;
      let nmerge = n - 1 in
      let joined_below r =
        let p = ref 0. and s = Array.make nlab 0. in
        for t = 0 to nmerge - 1 do
          if mh.(t) < r then begin
            p := !p +. mp.(t);
            for lv = 0 to nlab - 1 do s.(lv) <- s.(lv) +. msame.(lv).(t) done
          end
        done;
        !p /. float npairs, Array.mapi (fun lv x -> x /. total_same.(lv)) s in
      (* Cuts: the sweep applies merges in height order through a union-find *)
      let order = Array.init nmerge Fun.id in
      Array.sort (fun a b -> compare mh.(a) mh.(b)) order;
      let score_cuts heights =
        let parent = Array.init n Fun.id in
        let rec find x = if parent.(x) = x then x else begin let r = find parent.(x) in parent.(x) <- r; r end in
        let pos = ref 0 in
        List.map
          (fun h ->
            while !pos < nmerge && mh.(order.(!pos)) <= h do
              let a = find mi.(order.(!pos)) and b = find mj.(order.(!pos)) in
              if a <> b then parent.(a) <- b;
              incr pos
            done;
            let cl = Array.init n find in
            let csize = Hashtbl.create 1024 in
            Array.iter (fun c -> Hashtbl.replace csize c (1 + Option.value ~default:0 (Hashtbl.find_opt csize c))) cl;
            let k = Hashtbl.length csize and largest = Hashtbl.fold (fun _ v acc -> max v acc) csize 0 in
            let c2 x = x *. (x -. 1.) /. 2. and fn = float n in
            let per =
              Array.init nlab (fun lv ->
                  let nij = Hashtbl.create 4096 and bj = Hashtbl.create 64 in
                  Array.iteri
                    (fun p c ->
                      let l = lab.(lv).(p) in
                      let key = (c * 4096) + l in
                      Hashtbl.replace nij key (1 + Option.value ~default:0 (Hashtbl.find_opt nij key));
                      Hashtbl.replace bj l (1 + Option.value ~default:0 (Hashtbl.find_opt bj l)))
                    cl;
                  let sij = Hashtbl.fold (fun _ v acc -> acc +. c2 (float v)) nij 0.
                  and sa = Hashtbl.fold (fun _ v acc -> acc +. c2 (float v)) csize 0.
                  and sb = Hashtbl.fold (fun _ v acc -> acc +. c2 (float v)) bj 0. in
                  let expected = sa *. sb /. c2 fn in
                  let den = ((sa +. sb) /. 2.) -. expected in
                  let ari = if den = 0. then 1. else (sij -. expected) /. den in
                  let hc = Hashtbl.fold (fun _ v acc -> let p = float v /. fn in acc -. (p *. log p)) bj 0.
                  and hk = Hashtbl.fold (fun _ v acc -> let p = float v /. fn in acc -. (p *. log p)) csize 0. in
                  let hck = ref 0. and hkc = ref 0. in
                  Hashtbl.iter
                    (fun key v ->
                      let c = key / 4096 and l = key mod 4096 in
                      let pv = float v /. fn in
                      hck := !hck -. (pv *. log (float v /. float (Hashtbl.find csize c)));
                      hkc := !hkc -. (pv *. log (float v /. float (Hashtbl.find bj l))))
                    nij;
                  ari, (if hc > 0. then 1. -. (!hck /. hc) else 1.), (if hk > 0. then 1. -. (!hkc /. hk) else 1.)) in
            h, k, largest, per)
          heights in
      let cut_json at (h, k, largest, per) =
        Printf.sprintf "{\"at\":%S,\"r\":%.6g,\"k\":%d,\"largest\":%d,\"ari\":%s,\"hom\":%s,\"compl\":%s}" at h k largest
          (jarr (fun (a, _, _) -> Printf.sprintf "%.4f" a) per)
          (jarr (fun (_, x, _) -> Printf.sprintf "%.4f" x) per)
          (jarr (fun (_, _, x) -> Printf.sprintf "%.4f" x) per) in
      let lin_tree_valleys = ref [] in
      List.iter
        (fun scale ->
          let to_axis, of_axis, a0, a1 = if scale = "log" then log, exp, log lo, log hi else Fun.id, Fun.id, 0., hi in
          let w = (a1 -. a0) /. float bins in
          let factor = float m /. float npairs in
          let hist = Array.make bins 0. and hs = Array.init nlab (fun _ -> Array.make bins 0.) in
          for t = 0 to nmerge - 1 do
            let h = mh.(t) in
            if h > 0. || scale <> "log" then begin
              let b = int_of_float ((to_axis h -. a0) /. w) in
              if b >= 0 && b < bins then begin
                hist.(b) <- hist.(b) +. (mp.(t) *. factor);
                for lv = 0 to nlab - 1 do hs.(lv).(b) <- hs.(lv).(b) +. (msame.(lv).(t) *. factor) done
              end
            end
          done;
          let smooth s =
            Array.init bins (fun b ->
                let t = ref 0. and c = ref 0 in
                for j = b - s to b + s do if j >= 0 && j < bins then begin t := !t +. hist.(j); incr c end done;
                !t /. float !c) in
          (* The fixed detector of simpson/detector/detect3.ml, on the join-height histogram:
             extrema arranged peak-trough-...-peak, then adjacent peak-trough pairs removed
             smallest relative score first -- a trough against the lower of its flanking peaks,
             an interior peak against the modes beyond the troughs it sits between -- so a bump
             inside a gap no longer splits it.  What is NOT ported is the calibration: z comes
             from resampling the sequences, and there is no such resampling here, so these
             valleys carry a relative depth and no significance.  The share cap is the same.
             s is kept in the signature so the page*s two smoothing controls still resolve; the
             detector smooths by +-1 and both settings now answer alike. *)
          let valleys _s sm =
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
            let a = ref (Array.of_list (List.rev alt)) and go = ref true in
            while !go do
              let arr = !a in
              let len = Array.length arr in
              let v k = let _, _, x = arr.(k) in x in
              let best = ref infinity and lo = ref (-1) in
              let k = ref 1 in
              while !k < len - 1 do
                let p = Float.min (v (!k - 1)) (v (!k + 1)) in
                let sc = if p > 0. then (p -. v !k) /. p else 0. in
                if sc < !best then begin best := sc; lo := if v (!k - 1) < v (!k + 1) then !k - 1 else !k end;
                k := !k + 2
              done;
              let k = ref 2 in
              while !k < len - 2 do
                let outer = Float.min (v (!k - 2)) (v (!k + 2)) in
                let sc = if outer > 0. then (v !k -. Float.max (v (!k - 1)) (v (!k + 1))) /. outer else 0. in
                if sc < !best then begin best := sc; lo := if v (!k - 1) > v (!k + 1) then !k - 1 else !k end;
                k := !k + 2
              done;
              if !lo < 0 || !best >= 0.25 then go := false
              else a := Array.append (Array.sub arr 0 !lo) (Array.sub arr (!lo + 2) (len - !lo - 2))
            done;
            let arr = !a in
            List.init ((Array.length arr - 1) / 2) (fun i ->
                let _, tb, _ = arr.((2 * i) + 1) in
                (tb, of_axis (a0 +. ((float tb +. 0.5) *. w)))) in
          let vjson vl =
            String.concat ","
              (List.map
                 (fun (b, r) ->
                   let below, caught = joined_below r in
                   Printf.sprintf "{\"bin\":%d,\"r\":%.6g,\"below\":%.4f,\"caught\":%s}" b r below
                     (jarr (Printf.sprintf "%.4f") caught))
                 vl) in
          let sm5 = smooth 5 and sm2 = smooth 2 in
          let v5 = valleys 5 sm5 and v2 = valleys 2 sm2 in
          if scale = "lin" then lin_tree_valleys := List.map snd v5;
          let extra =
            if scale <> "lin" then ""
            else begin
              let pr = pair_cuts d and tr = !lin_tree_valleys in
              let grid = List.init 120 (fun g -> hi *. float (g + 1) /. 120.) in
              let all = List.sort_uniq compare (grid @ pr @ tr) in
              let scored = score_cuts all in
              let find r = List.find (fun (h, _, _, _) -> h = r) scored in
              let cuts = List.map (fun r -> cut_json "pairs" (find r)) pr @ List.map (fun r -> cut_json "tree" (find r)) tr in
              let best =
                List.init nlab (fun lv ->
                    let b =
                      List.fold_left
                        (fun acc ((_, _, _, per) as x) ->
                          match acc with
                          | None -> Some x
                          | Some (_, _, _, pb) -> let a, _, _ = per.(lv) and c, _, _ = pb.(lv) in if a > c then Some x else acc)
                        None scored in
                    cut_json (List.nth lnames lv) (Option.get b)) in
              Printf.sprintf ",\"cuts\":[%s],\"best\":[%s]" (String.concat "," cuts) (String.concat "," best)
            end in
          Printf.printf
            "{\"tag\":%S,\"d\":%d,\"axes\":%d,\"scale\":%S,\"a0\":%.6g,\"w\":%.6g,\"pairs\":%d,\"height\":%.6g,\"levels\":%s,\"base\":%s,\"counts\":%s,\"same\":%s,\"smooth5\":%s,\"smooth2\":%s,\"valleys5\":[%s],\"valleys2\":[%s]%s}\n%!"
            tag d maxd scale a0 w m
            (Array.fold_left Float.max 0. mh)
            (jarr (Printf.sprintf "%S") (Array.of_list lnames))
            (jarr (fun t -> Printf.sprintf "%.4f" (t /. float npairs)) total_same)
            (jarr (Printf.sprintf "%.1f") hist)
            (jarr (fun h -> jarr (Printf.sprintf "%.1f") h) hs)
            (jarr (Printf.sprintf "%.1f") sm5) (jarr (Printf.sprintf "%.1f") sm2)
            (vjson v5) (vjson v2) extra)
        [ "lin"; "log" ])
    ds
