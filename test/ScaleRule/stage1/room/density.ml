(* density <twisted.txt> <inertia.txt> <labels> <tag> <detect3.cand> [<queries>]

   Whether a cluster density can push the axis ladder back from too many axes.  At 1, 2, 4, ...
   axes and at all of them, for two partitions of the labelled spectra -- the labels themselves,
   which do not change with the rung, and the leader clustering at the rung's finest usable valley
   (kept, share <= 0.5, radius from detect3's candidate file), which is where a search starts --
   it reports:
   - the two-nearest-neighbour intrinsic dimension at that rung, over a sample of queries;
   - R_c / R_all: the RMS distance of a cluster's members to its centroid over that of all points to
     theirs, the median over members;
   - the log relative density of a cluster, log((n_c / N) * (R_all / R_c)^e), averaged over
     members, once with e the number of axes and once with e the intrinsic dimension.  A density
     in d axes has units of length^-d, so only the ratio to the whole set's density in the same
     axes can be compared between rungs;
   - how many clusters are singletons and how many have zero volume (identical members), which
     both have no finite density and are left out of the averages.
   Coordinates are weighted by the square root of the inertia, as in detect3, whose radii are in
   the same units. *)

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4)
  and candf = Sys.argv.(5) in
  let nq_req = if Array.length Sys.argv > 6 then int_of_string Sys.argv.(6) else 1000 in
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
  (* The finest usable kept valley at each rung of this embedding, as (share, radius) *)
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
  (* Intrinsic dimension per rung, over a deterministic sample of queries *)
  let st = Random.State.make [| 20260915 |] in
  let perm = Array.init n Fun.id in
  for i = n - 1 downto 1 do
    let j = Random.State.int st (i + 1) in
    let t = perm.(i) in
    perm.(i) <- perm.(j);
    perm.(j) <- t
  done;
  let nq = min nq_req n in
  let id_tot = Array.make nr 0. and id_cnt = Array.make nr 0 in
  let d2 = Array.make n 0. in
  for q = 0 to nq - 1 do
    let qi = perm.(q) in
    Array.fill d2 0 n 0.;
    let prev = ref 0 in
    Array.iteri
      (fun ri d ->
        for a = !prev to d - 1 do
          let x = pc.(qi).(a) in
          for j = 0 to n - 1 do
            let y = x -. pc.(j).(a) in
            d2.(j) <- d2.(j) +. (y *. y)
          done
        done;
        prev := d;
        let b1 = ref infinity and b2 = ref infinity in
        for j = 0 to n - 1 do
          if j <> qi then begin
            let v = d2.(j) in
            if v < !b1 then begin b2 := !b1; b1 := v end else if v < !b2 then b2 := v
          end
        done;
        if !b1 > 0. then begin
          id_tot.(ri) <- id_tot.(ri) +. (0.5 *. log (!b2 /. !b1));
          id_cnt.(ri) <- id_cnt.(ri) + 1
        end)
      rungs
  done;
  (* Statistics of one partition in the first d axes *)
  let stats d idim assign =
    let ncl = Array.fold_left max (-1) assign + 1 in
    let size = Array.make ncl 0 and cent = Array.make_matrix ncl d 0. and all = Array.make d 0. in
    Array.iteri
      (fun i c ->
        size.(c) <- size.(c) + 1;
        for a = 0 to d - 1 do
          cent.(c).(a) <- cent.(c).(a) +. pc.(i).(a);
          all.(a) <- all.(a) +. pc.(i).(a)
        done)
      assign;
    for c = 0 to ncl - 1 do
      if size.(c) > 0 then
        for a = 0 to d - 1 do cent.(c).(a) <- cent.(c).(a) /. float_of_int size.(c) done
    done;
    for a = 0 to d - 1 do all.(a) <- all.(a) /. float_of_int n done;
    let ss = Array.make ncl 0. and ss_all = ref 0. in
    Array.iteri
      (fun i c ->
        for a = 0 to d - 1 do
          let y = pc.(i).(a) -. cent.(c).(a) and z = pc.(i).(a) -. all.(a) in
          ss.(c) <- ss.(c) +. (y *. y);
          ss_all := !ss_all +. (z *. z)
        done)
      assign;
    let r_all = sqrt (!ss_all /. float_of_int n) in
    let live = ref 0 and singletons = ref 0 and zero = ref 0 and counted = ref 0
    and s_nom = ref 0. and s_id = ref 0. and ratios = ref [] in
    for c = 0 to ncl - 1 do
      if size.(c) > 0 then begin
        incr live;
        if size.(c) = 1 then incr singletons
        else if ss.(c) <= 0. then incr zero
        else begin
          let m = size.(c) in
          let ratio = sqrt (ss.(c) /. float_of_int m) /. r_all in
          let base = log (float_of_int m /. float_of_int n) in
          counted := !counted + m;
          s_nom := !s_nom +. (float_of_int m *. (base -. (float_of_int d *. log ratio)));
          s_id := !s_id +. (float_of_int m *. (base -. (idim *. log ratio)));
          ratios := (ratio, m) :: !ratios
        end
      end
    done;
    let median =
      let sorted = List.sort compare !ratios in
      let half = !counted / 2 and acc = ref 0 and res = ref nan in
      List.iter (fun (r, m) -> if Float.is_nan !res then begin acc := !acc + m; if !acc > half then res := r end) sorted;
      !res in
    let per x = if !counted > 0 then x /. float_of_int !counted else nan in
    !live, !singletons, !zero, !counted, median, per !s_nom, per !s_id in
  Printf.printf "# %s: %d labelled spectra, %d queries, %d axes available\n" tag n nq avail;
  Printf.printf
    "tag\trung\tintrinsic_dimension\tpartition\tvalley_share\tclusters\tsingletons\tzero_volume\tpoints_counted\tmedian_Rc_over_Rall\tlog_rel_density_axes\tlog_rel_density_id\n";
  Array.iteri
    (fun ri d ->
      let idim = if id_cnt.(ri) > 0 && id_tot.(ri) > 0. then float_of_int id_cnt.(ri) /. id_tot.(ri) else nan in
      let line kind share (live, singletons, zero, counted, median, s_nom, s_id) =
        Printf.printf "%s\t%d\t%.2f\t%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.2f\t%.2f\n%!" tag d idim kind share live
          singletons zero counted median s_nom s_id in
      line "labels" "-" (stats d idim lab);
      match Hashtbl.find_opt finest d with
      | None -> Printf.printf "%s\t%d\t%.2f\tvalley\tnone\t-\t-\t-\t-\t-\t-\t-\n%!" tag d idim
      | Some (share, radius) ->
        (* Leader clustering at the valley radius, in index order, as run_montecarlo starts *)
        let assign = Array.make n (-1) and leaders = ref [] and nl = ref 0 in
        for i = 0 to n - 1 do
          let best = ref (-1) and best_d = ref infinity in
          List.iter
            (fun l ->
              let s = ref 0. and a = ref 0 in
              while !a < d && !s <= radius *. radius do
                let y = pc.(i).(!a) -. pc.(l).(!a) in
                s := !s +. (y *. y);
                incr a
              done;
              if !s <= radius *. radius && !s < !best_d then begin best_d := !s; best := l end)
            !leaders;
          if !best >= 0 then assign.(i) <- assign.(!best)
          else begin leaders := i :: !leaders; assign.(i) <- !nl; incr nl end
        done;
        line "valley" (Printf.sprintf "%.4f" share) (stats d idim assign))
    rungs
