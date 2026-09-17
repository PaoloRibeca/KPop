(* silhgrid <twisted.txt> <inertia.txt> <labels> <tag> <detect4.cand> <replicates> <rungs,> <flat|powers> <euclidean|manhattan|angle>

   silhnull under a chosen metric and distance (dist.ml).  At every rung given: the classical
   silhouette of the labels, to judge and never part of a rule; the partition leader clustering finds
   at the rung's finest usable valley -- kept, share <= 0.5, radius from detect4's candidate file run
   under the same metric and distance -- its number of clusters k and its silhouette; and, over
   <replicates> copies of the embedding with every axis's coordinates shuffled across the spectra,
   the silhouette of the partition leader clustering finds there with its radius tuned to give k
   clusters.  The excess of the real silhouette over the mean shuffled one is what the partition owes
   to structure rather than to having few clusters. *)

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4)
  and candf = Sys.argv.(5) and replicates = int_of_string Sys.argv.(6)
  and rungs_req = String.split_on_char ',' Sys.argv.(7) |> List.map int_of_string
  and flat = Sys.argv.(8) = "flat" and kind = Dist.kind_of_string Sys.argv.(9) in
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
  let pc = Array.map (fun r -> Array.init avail (fun a -> Dist.weigh kind ~flat iv coords.(r).(a) a)) rows in
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
  let np = n * (n - 1) / 2 in
  let acc = Float.Array.make np 0. and dist = Float.Array.make np 0. and norms = Array.make n 0. in
  (* The radius at which leader clustering gives [target] clusters, by bisection *)
  let leader_k mat d target =
    let top = Dist.top_radius kind mat d in
    let best = ref (Dist.leader kind mat norms d top) and lo = ref 0. and hi = ref top and it = ref 0 in
    while !it < 40 && snd !best <> target do
      let mid = (!lo +. !hi) /. 2. in
      let (_, k) as p = Dist.leader kind mat norms d mid in
      if abs (k - target) < abs (snd !best - target) then best := p;
      if k > target then lo := mid else hi := mid;
      incr it
    done;
    !best in
  let s_lab = Array.make nr nan and own_k = Array.make nr 0 and s_own = Array.make nr nan in
  let prev = ref 0 in
  Array.iteri
    (fun ri d ->
      Dist.accumulate kind pc acc norms !prev d;
      Dist.to_dist kind n acc norms dist;
      prev := d;
      s_lab.(ri) <- Dist.silhouette n dist lab nlab;
      match Hashtbl.find_opt finest d with
      | Some (_, radius) ->
        let a, k = Dist.leader kind pc norms d radius in
        own_k.(ri) <- k;
        s_own.(ri) <- Dist.silhouette n dist a k
      | None -> ())
    rungs;
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
    Float.Array.fill acc 0 np 0.;
    Array.fill norms 0 n 0.;
    let prev = ref 0 in
    Array.iteri
      (fun ri d ->
        Dist.accumulate kind sh acc norms !prev d;
        Dist.to_dist kind n acc norms dist;
        prev := d;
        if own_k.(ri) > 0 then begin
          let a, k = leader_k sh d own_k.(ri) in
          null_k.(ri).(r) <- k;
          null_s.(ri).(r) <- Dist.silhouette n dist a k
        end)
      rungs;
    Printf.eprintf "%s: shuffle %d of %d done\n%!" tag (r + 1) replicates
  done;
  Printf.printf "# %s: %d labelled spectra, %d label classes, %d axes available, %d shuffles, %s %s\n" tag n
    nlab avail replicates (if flat then "flat" else "powers(1,1,1)") (Dist.kind_to_string kind);
  Printf.printf "tag\trung\tsil_labels\town_clusters\tsil_own\tsil_null_mean\tsil_null_sd\tnull_clusters_mean\texcess\n";
  Array.iteri
    (fun ri d ->
      if own_k.(ri) = 0 then Printf.printf "%s\t%d\t%.4f\t-\t-\t-\t-\t-\t-\n" tag d s_lab.(ri)
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
