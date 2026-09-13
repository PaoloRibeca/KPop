(* treecut <twisted.txt> <inertia.txt> <d> <tag> <r,r,...|-> <labels> [<labels>]
   Does a valley in the pair-distance histogram mark a level of the tree?  Builds the average-
   linkage (UPGMA) dendrogram of the labelled spectra over the first d principal axes, by the
   nearest-neighbour chain, then (1) sets the histogram of all pairwise distances beside the
   histogram of cophenetic distances -- the height at which the tree joins each pair -- valley for
   valley, and (2) cuts the tree at a grid of heights and at each valley radius given, scoring
   every cut against each label level *)
let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and d = int_of_string Sys.argv.(3) and tag = Sys.argv.(4) in
  let radii = if Sys.argv.(5) = "-" then [] else List.map float_of_string (String.split_on_char ',' Sys.argv.(5)) in
  let nlab = Array.length Sys.argv - 6 in
  let labs = Array.init nlab (fun i -> Kio.read_labels Sys.argv.(6 + i)) in
  let iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs.(0) names_all.(i)) |> Array.of_list in
  let n = Array.length rows and d = min d (Array.length iv) in
  let pc = Array.map (fun r -> Array.init d (fun k -> coords_all.(r).(k) *. sqrt iv.(k))) rows in
  let lab =
    Array.map
      (fun h ->
        let tbl = Hashtbl.create 64 in
        Array.map
          (fun r ->
            let s = Hashtbl.find h names_all.(r) in
            match Hashtbl.find_opt tbl s with
            | Some id -> id
            | None -> let id = Hashtbl.length tbl in Hashtbl.add tbl s id; id)
          rows)
      labs in
  let npairs = n * (n - 1) / 2 in
  let idx i j = let i, j = if i < j then i, j else j, i in i * (2 * n - i - 1) / 2 + (j - i - 1) in
  let dm = Float.Array.make npairs 0. in
  for i = 0 to n - 2 do
    let pi = pc.(i) and base = i * (2 * n - i - 1) / 2 - i - 1 in
    for j = i + 1 to n - 1 do
      let pj = pc.(j) and s = ref 0. in
      for k = 0 to d - 1 do let x = pi.(k) -. pj.(k) in s := !s +. (x *. x) done;
      Float.Array.unsafe_set dm (base + j) (sqrt !s)
    done
  done;
  (* The detector's range: the 99th percentile of its 100000 seed-17 pairs *)
  let state = Random.State.make [| 17 |] in
  let samp = Array.make 100000 0. and m = ref 0 in
  for _ = 1 to 100000 do
    let a = Random.State.int state n and b = Random.State.int state n in
    if a <> b then begin samp.(!m) <- Float.Array.get dm (idx a b); incr m end
  done;
  let samp = Array.sub samp 0 !m in
  Array.sort compare samp;
  let hi = samp.(int_of_float (float !m *. 0.99)) and bins = 200 in
  let w = hi /. float bins in
  let raw = Array.make bins 0. in
  Float.Array.iter (fun x -> let b = int_of_float (x /. w) in if b >= 0 && b < bins then raw.(b) <- raw.(b) +. 1.) dm;
  let below_raw r = let c = ref 0 in Float.Array.iter (fun x -> if x < r then incr c) dm; float !c /. float npairs in
  (* Average linkage by the nearest-neighbour chain, Lance-Williams updates in place *)
  let active = Array.make n true and size = Array.make n 1 in
  let chain = Array.make (n + 1) 0 and top = ref 0 in
  let mi = Array.make (n - 1) 0 and mj = Array.make (n - 1) 0 and mh = Array.make (n - 1) 0.
  and ms = Array.make (n - 1) 0. and nm = ref 0 and start = ref 0 in
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
      mi.(!nm) <- i; mj.(!nm) <- j; mh.(!nm) <- !bd; ms.(!nm) <- si *. sj; incr nm;
      for k = 0 to n - 1 do
        if active.(k) && k <> i && k <> j then
          Float.Array.set dm (idx i k)
            (((si *. Float.Array.get dm (idx i k)) +. (sj *. Float.Array.get dm (idx j k))) /. (si +. sj))
      done;
      size.(i) <- size.(i) + size.(j);
      active.(j) <- false
    end else begin chain.(!top) <- !best; incr top end
  done;
  let coph = Array.make bins 0. in
  for t = 0 to n - 2 do
    let b = int_of_float (mh.(t) /. w) in
    if b >= 0 && b < bins then coph.(b) <- coph.(b) +. ms.(t)
  done;
  let below_coph r = let c = ref 0. in for t = 0 to n - 2 do if mh.(t) < r then c := !c +. ms.(t) done; !c /. float npairs in
  let valleys hist =
    let s = 5 in
    let sm = Array.init bins (fun b -> let t = ref 0. and c = ref 0 in for j = b - s to b + s do if j >= 0 && j < bins then begin t := !t +. hist.(j); incr c end done; !t /. float !c) in
    let minb a b = let v = ref sm.(min a b) in for j = min a b to max a b do if sm.(j) < !v then v := sm.(j) done; !v in
    let argm a b = let v = ref (min a b) in for j = min a b to max a b do if sm.(j) < sm.(!v) then v := j done; !v in
    let is_peak b = sm.(b) > 0. && (let ok = ref true in for j = b - s to b + s do if j >= 0 && j < bins && j <> b && (sm.(j) > sm.(b) || (sm.(j) = sm.(b) && j < b)) then ok := false done; !ok) in
    let rec merge = function
      | p1 :: p2 :: rest -> if minb p1 p2 <= 0.5 *. Float.min sm.(p1) sm.(p2) then p1 :: merge (p2 :: rest) else merge ((if sm.(p1) >= sm.(p2) then p1 else p2) :: rest)
      | l -> l in
    let rec fix ps = let ps' = merge ps in if ps' = ps then ps else fix ps' in
    let rec vs = function p1 :: (p2 :: _ as rest) -> ((float (argm p1 p2) +. 0.5) *. w) :: vs rest | _ -> [] in
    vs (fix (List.filter is_peak (List.init bins Fun.id))) in
  let vraw = valleys raw and vcoph = valleys coph in
  Printf.printf "# %s, %d axes, %d labelled spectra; tree height %.4g; histogram range 0-%.4g\n" tag d n
    (Array.fold_left Float.max 0. mh) hi;
  Printf.printf "  valleys, all pairwise distances:  %s\n"
    (String.concat "  " (List.map (fun r -> Printf.sprintf "%.3g (%.1f%% below)" r (100. *. below_raw r)) vraw));
  Printf.printf "  valleys, tree joining heights:    %s\n"
    (String.concat "  " (List.map (fun r -> Printf.sprintf "%.3g (%.1f%% joined below)" r (100. *. below_coph r)) vcoph));
  (* Cuts *)
  let order = Array.init (n - 1) Fun.id in
  Array.sort (fun a b -> compare mh.(a) mh.(b)) order;
  let parent = Array.init n Fun.id in
  let rec find x = if parent.(x) = x then x else begin let r = find parent.(x) in parent.(x) <- r; r end in
  let score () =
    let cl = Array.init n find in
    let csize = Hashtbl.create 1024 in
    Array.iter (fun c -> Hashtbl.replace csize c (1 + try Hashtbl.find csize c with Not_found -> 0)) cl;
    let k = Hashtbl.length csize and largest = Hashtbl.fold (fun _ v acc -> max v acc) csize 0 in
    let c2 x = x *. (x -. 1.) /. 2. and fn = float n in
    let per =
      Array.init nlab (fun lv ->
        let nij = Hashtbl.create 4096 and bj = Hashtbl.create 64 in
        Array.iteri
          (fun p c ->
            let l = lab.(lv).(p) in
            let key = (c * 4096) + l in
            Hashtbl.replace nij key (1 + try Hashtbl.find nij key with Not_found -> 0);
            Hashtbl.replace bj l (1 + try Hashtbl.find bj l with Not_found -> 0))
          cl;
        let sij = Hashtbl.fold (fun _ v acc -> acc +. c2 (float v)) nij 0.
        and sa = Hashtbl.fold (fun _ v acc -> acc +. c2 (float v)) csize 0.
        and sb = Hashtbl.fold (fun _ v acc -> acc +. c2 (float v)) bj 0. in
        let expected = sa *. sb /. c2 fn in
        let ari = if (sa +. sb) /. 2. -. expected = 0. then 1. else (sij -. expected) /. (((sa +. sb) /. 2.) -. expected) in
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
    k, largest, per in
  let grid = List.sort_uniq compare (List.init 120 (fun g -> hi *. float (g + 1) /. 120.) @ radii) in
  let pos = ref 0 and results = ref [] in
  List.iter
    (fun h ->
      while !pos < n - 1 && mh.(order.(!pos)) <= h do
        let a = find mi.(order.(!pos)) and b = find mj.(order.(!pos)) in
        if a <> b then parent.(a) <- b;
        incr pos
      done;
      results := (h, score ()) :: !results)
    grid;
  let results = List.rev !results in
  let lvname lv = if nlab = 2 && lv = 1 then "finest label" else "class" in
  let show (h, (k, largest, per)) =
    Printf.sprintf "h=%-7.3g %5d clusters, largest %4d |%s" h k largest
      (String.concat " |"
         (Array.to_list (Array.mapi (fun lv (a, hm, cm) -> Printf.sprintf " %s ARI %.3f hom %.3f compl %.3f" (lvname lv) a hm cm) per))) in
  List.iter (fun r -> Printf.printf "  cut at valley %s\n" (show (List.find (fun (h, _) -> h = r) results))) radii;
  for lv = 0 to nlab - 1 do
    let best = List.fold_left (fun acc ((_, (_, _, per)) as x) -> match acc with None -> Some x | Some (_, (_, _, pb)) -> let a, _, _ = per.(lv) and b, _, _ = pb.(lv) in if a > b then Some x else acc) None results in
    Option.iter (fun x -> Printf.printf "  best cut for %s:  %s\n" (lvname lv) (show x)) best
  done;
  print_endline ""
