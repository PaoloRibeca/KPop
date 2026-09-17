(* silhnull <twisted.txt> <inertia.txt> <labels> <tag> <detect3.cand> <replicates> <rungs,>

   Whether each rung's own partition is better separated than the same number of clusters would be
   with no structure at all.  At every rung given:
   - the classical silhouette of the labels, to judge and never part of a rule;
   - the partition leader clustering finds at the rung's finest usable valley (kept, share <= 0.5,
     radius from detect3's candidate file), its number of clusters k and its silhouette;
   - over <replicates> copies of the embedding in which every axis has had its coordinates shuffled
     across the spectra -- which keeps each axis's spread and destroys what the axes share -- the
     silhouette of the partition leader clustering finds there, its radius tuned by bisection to
     give k clusters as nearly as it can.
   The excess of the real silhouette over the mean shuffled one is what the partition owes to
   structure rather than to having few clusters.  Coordinates are weighted by the square root of
   the inertia, as in detect3, whose radii are in the same units; a singleton scores 0. *)

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4)
  and candf = Sys.argv.(5) and replicates = int_of_string Sys.argv.(6)
  and rungs_req = String.split_on_char ',' Sys.argv.(7) |> List.map int_of_string in
  let labs = Kio.read_labels labf and iv = Kio.read_inertia inf in
  let names, coords = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs names.(i))
    |> Array.of_list in
  let n = Array.length rows in
  let avail = min (Array.length iv) (Array.length coords.(rows.(0))) in
  let rungs = List.map (fun d -> min d avail) rungs_req |> List.sort_uniq compare |> Array.of_list in
  let nr = Array.length rungs in
  let pc = Array.map (fun r -> Array.init avail (fun a -> coords.(r).(a) *. sqrt iv.(a))) rows in
  let lab, nlab =
    let tbl = Hashtbl.create 64 in
    let a =
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
    a, Hashtbl.length tbl in
  let finest = Hashtbl.create 16 in
  let ic = open_in candf in
  (try
     while true do
       match String.split_on_char '\t' (input_line ic) with
       | [ "CAND"; t; rung; radius; share; _; _; _; "kept" ] when t = tag ->
         let rung = int_of_string rung and share = float_of_string share
         and radius = float_of_string radius in
         if share <= 0.5 then begin
           match Hashtbl.find_opt finest rung with
           | Some (s, _) when s <= share -> ()
           | _ -> Hashtbl.replace finest rung (share, radius)
         end
       | _ -> ()
     done
   with End_of_file -> ());
  close_in ic;
  (* Leader clustering of a coordinate matrix at a radius in its first d axes, in index order *)
  let leader mat d radius =
    let assign = Array.make n (-1) and leaders = ref [] and nl = ref 0 in
    let r2 = radius *. radius in
    for i = 0 to n - 1 do
      let row = mat.(i) and best = ref (-1) and best_d = ref infinity in
      List.iter
        (fun l ->
          let lrow = mat.(l) and s = ref 0. and a = ref 0 in
          while !a < d && !s <= r2 do
            let y = row.(!a) -. lrow.(!a) in
            s := !s +. (y *. y);
            incr a
          done;
          if !s <= r2 && !s < !best_d then begin
            best_d := !s;
            best := l
          end)
        !leaders;
      if !best >= 0 then assign.(i) <- assign.(!best)
      else begin
        leaders := i :: !leaders;
        assign.(i) <- !nl;
        incr nl
      end
    done;
    assign, !nl in
  (* The radius at which leader clustering gives [target] clusters, by bisection: fewer clusters
     as the radius grows.  The partition whose count is nearest the target is returned *)
  let leader_k mat d target =
    let top =
      let s = ref 0. in
      for a = 0 to d - 1 do
        let lo = ref infinity and hi = ref neg_infinity in
        Array.iter (fun row -> let v = row.(a) in if v < !lo then lo := v; if v > !hi then hi := v) mat;
        s := !s +. ((!hi -. !lo) *. (!hi -. !lo))
      done;
      sqrt !s *. 1.000001 in
    let best = ref (leader mat d top) and lo = ref 0. and hi = ref top and it = ref 0 in
    while !it < 40 && snd !best <> target do
      let mid = (!lo +. !hi) /. 2. in
      let (_, k) as p = leader mat d mid in
      if abs (k - target) < abs (snd !best - target) then best := p;
      if k > target then lo := mid else hi := mid;
      incr it
    done;
    !best in
  let np = n * (n - 1) / 2 in
  let sq = Float.Array.make np 0. and col = Array.make n 0. in
  let accumulate mat from_ to_ =
    for a = from_ to to_ - 1 do
      for i = 0 to n - 1 do col.(i) <- mat.(i).(a) done;
      let p = ref 0 in
      for i = 0 to n - 2 do
        let x = col.(i) in
        for j = i + 1 to n - 1 do
          let y = x -. col.(j) in
          Float.Array.unsafe_set sq !p (Float.Array.unsafe_get sq !p +. (y *. y));
          incr p
        done
      done
    done in
  let silhouette assign k =
    let s = Array.make (n * k) 0. and cnt = Array.make k 0 in
    Array.iter (fun c -> cnt.(c) <- cnt.(c) + 1) assign;
    let p = ref 0 in
    for i = 0 to n - 2 do
      let ci = assign.(i) and base_i = i * k in
      for j = i + 1 to n - 1 do
        let dist = sqrt (Float.Array.unsafe_get sq !p) in
        incr p;
        let cj = assign.(j) in
        s.(base_i + cj) <- s.(base_i + cj) +. dist;
        s.((j * k) + ci) <- s.((j * k) + ci) +. dist
      done
    done;
    let tot = ref 0. in
    for i = 0 to n - 1 do
      let ci = assign.(i) in
      if cnt.(ci) > 1 then begin
        let a = s.((i * k) + ci) /. float_of_int (cnt.(ci) - 1) and b = ref infinity in
        for c = 0 to k - 1 do
          if c <> ci && cnt.(c) > 0 then begin
            let m = s.((i * k) + c) /. float_of_int cnt.(c) in
            if m < !b then b := m
          end
        done;
        let den = Float.max a !b in
        if !b < infinity && den > 0. then tot := !tot +. ((!b -. a) /. den)
      end
    done;
    !tot /. float_of_int n in
  (* The real embedding: the labels' silhouette, and each rung's own partition *)
  let s_lab = Array.make nr nan and own_k = Array.make nr 0 and s_own = Array.make nr nan in
  let prev = ref 0 in
  Array.iteri
    (fun ri d ->
      accumulate pc !prev d;
      prev := d;
      s_lab.(ri) <- silhouette lab nlab;
      match Hashtbl.find_opt finest d with
      | Some (_, radius) ->
        let a, k = leader pc d radius in
        own_k.(ri) <- k;
        s_own.(ri) <- silhouette a k
      | None -> ())
    rungs;
  (* The shuffled copies, one at a time *)
  let null_s = Array.make_matrix nr replicates nan and null_k = Array.make_matrix nr replicates 0 in
  for r = 0 to replicates - 1 do
    let st = Random.State.make [| 20260916; r |] in
    let sh = Array.map Array.copy pc in
    for a = 0 to avail - 1 do
      for i = n - 1 downto 1 do
        let j = Random.State.int st (i + 1) in
        let t = sh.(i).(a) in
        sh.(i).(a) <- sh.(j).(a);
        sh.(j).(a) <- t
      done
    done;
    Float.Array.fill sq 0 np 0.;
    let prev = ref 0 in
    Array.iteri
      (fun ri d ->
        accumulate sh !prev d;
        prev := d;
        if own_k.(ri) > 0 then begin
          let a, k = leader_k sh d own_k.(ri) in
          null_k.(ri).(r) <- k;
          null_s.(ri).(r) <- silhouette a k
        end)
      rungs;
    Printf.eprintf "%s: shuffle %d of %d done\n%!" tag (r + 1) replicates
  done;
  Printf.printf "# %s: %d labelled spectra, %d label classes, %d axes available, %d shuffles\n" tag n nlab
    avail replicates;
  Printf.printf "tag\trung\tsil_labels\town_clusters\tsil_own\tsil_null_mean\tsil_null_sd\tnull_clusters_mean\texcess\n";
  Array.iteri
    (fun ri d ->
      if own_k.(ri) = 0 then
        Printf.printf "%s\t%d\t%.4f\t-\t-\t-\t-\t-\t-\n" tag d s_lab.(ri)
      else begin
        let v = null_s.(ri) and m = float_of_int replicates in
        let mean = Array.fold_left ( +. ) 0. v /. m in
        let sd =
          if replicates > 1 then
            sqrt (Array.fold_left (fun acc x -> acc +. ((x -. mean) *. (x -. mean))) 0. v /. (m -. 1.))
          else nan in
        let km = float_of_int (Array.fold_left ( + ) 0 null_k.(ri)) /. m in
        Printf.printf "%s\t%d\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.1f\t%.4f\n" tag d s_lab.(ri) own_k.(ri) s_own.(ri)
          mean sd km (s_own.(ri) -. mean)
      end)
    rungs
