(* silh <twisted.txt> <inertia.txt> <labels> <tag> <detect3.cand> <pick> <rungs,>

   How the silhouette changes as the axes grow.  At every rung given, over all pairs of the labelled
   spectra, the classical silhouette (a singleton scoring 0) of three partitions:
   - the labels, to judge the curves and never part of a rule;
   - the partition leader clustering finds at that rung's finest usable valley (kept, share <= 0.5,
     radius from detect3's candidate file), in that rung's own space -- the partition a search there
     starts from;
   - the partition leader clustering finds at the rung <pick>, built once and carried unchanged to
     every rung, so that its silhouette moves only because axes are added or taken away.
   Coordinates are weighted by the square root of the inertia, as in detect3, whose radii are in the
   same units.  Squared distances are accumulated axis by axis, so every rung costs one pass per
   new axis and one per partition. *)

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4)
  and candf = Sys.argv.(5) and pick = int_of_string Sys.argv.(6)
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
  (* The finest usable kept valley at each rung, as (share, radius) *)
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
  (* Leader clustering at the valley radius in the first d axes, in index order *)
  let leader d radius =
    let assign = Array.make n (-1) and leaders = ref [] and nl = ref 0 in
    let r2 = radius *. radius in
    for i = 0 to n - 1 do
      let best = ref (-1) and best_d = ref infinity in
      List.iter
        (fun l ->
          let s = ref 0. and a = ref 0 in
          while !a < d && !s <= r2 do
            let y = pc.(i).(!a) -. pc.(l).(!a) in
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
  let np = n * (n - 1) / 2 in
  let sq = Float.Array.make np 0. in
  (* The classical silhouette of a partition on the distances accumulated so far *)
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
  let carried = Option.map (fun (_, radius) -> leader pick radius) (Hashtbl.find_opt finest pick) in
  Printf.printf "# %s: %d labelled spectra, %d label classes, %d axes available, pick %d%s\n%!" tag n nlab
    avail pick
    (match carried with Some (_, k) -> Printf.sprintf " (%d clusters carried)" k | None -> " (no valley there)");
  Printf.printf "tag\trung\tsil_labels\town_share\town_clusters\tsil_own\tpick\tpick_clusters\tsil_carried\n%!";
  let prev = ref 0 and col = Array.make n 0. in
  Array.iter
    (fun d ->
      for a = !prev to d - 1 do
        for i = 0 to n - 1 do col.(i) <- pc.(i).(a) done;
        let p = ref 0 in
        for i = 0 to n - 2 do
          let x = col.(i) in
          for j = i + 1 to n - 1 do
            let y = x -. col.(j) in
            Float.Array.unsafe_set sq !p (Float.Array.unsafe_get sq !p +. (y *. y));
            incr p
          done
        done
      done;
      prev := d;
      let sl = silhouette lab nlab in
      let own =
        match Hashtbl.find_opt finest d with
        | Some (share, radius) ->
          let a, k = leader d radius in
          Printf.sprintf "%.4f\t%d\t%.4f" share k (silhouette a k)
        | None -> "-\t-\t-" in
      let car =
        match carried with
        | Some (a, k) -> Printf.sprintf "%d\t%d\t%.4f" pick k (silhouette a k)
        | None -> Printf.sprintf "%d\t-\t-" pick in
      Printf.printf "%s\t%d\t%.4f\t%s\t%s\n%!" tag d sl own car)
    rungs
