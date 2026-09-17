(* hc <twisted.txt> <inertia.txt> <labels> <tag> <detect3.cand>

   Homogeneity and completeness against the labels of the leader partition at each rung's finest
   usable kept valley, in index order, as density.ml builds it: whether the partition a rung hands
   the search is already coarser than the labels (low homogeneity) or finer (low completeness). *)

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and labf = Sys.argv.(3) and tag = Sys.argv.(4)
  and candf = Sys.argv.(5) in
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
  (* Homogeneity and completeness of partition a against reference b, as Clustering.compare_partitions *)
  let hc a b =
    let cell = Hashtbl.create 64 and ca = Hashtbl.create 64 and cb = Hashtbl.create 64 in
    let bump t k = Hashtbl.replace t k (1 + try Hashtbl.find t k with Not_found -> 0) in
    for i = 0 to n - 1 do bump cell (a.(i), b.(i)); bump ca a.(i); bump cb b.(i) done;
    let fn = float_of_int n in
    let log2 x = log x /. log 2. in
    let entropy t = Hashtbl.fold (fun _ v acc -> let p = float_of_int v /. fn in acc -. (p *. log2 p)) t 0. in
    let h_a = entropy ca and h_b = entropy cb in
    let conditional on_b =
      Hashtbl.fold
        (fun (ka, kb) v acc ->
          let p = float_of_int v /. fn
          and marginal = float_of_int (Hashtbl.find (if on_b then cb else ca) (if on_b then kb else ka)) in
          acc -. (p *. log2 (float_of_int v /. marginal)))
        cell 0. in
    let homogeneity = if h_b > 0. then 1. -. (conditional false /. h_b) else 1.
    and completeness = if h_a > 0. then 1. -. (conditional true /. h_a) else 1. in
    homogeneity, completeness in
  Printf.printf "tag\trung\tshare\tclusters\thomogeneity\tcompleteness\n";
  Array.iter
    (fun d ->
      match Hashtbl.find_opt finest d with
      | None -> Printf.printf "%s\t%d\tnone\t-\t-\t-\n%!" tag d
      | Some (share, radius) ->
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
              if !s <= r2 && !s < !best_d then begin best_d := !s; best := l end)
            !leaders;
          if !best >= 0 then assign.(i) <- assign.(!best)
          else begin leaders := i :: !leaders; assign.(i) <- !nl; incr nl end
        done;
        let h, c = hc assign lab in
        Printf.printf "%s\t%d\t%.4f\t%d\t%.3f\t%.3f\n%!" tag d share !nl h c)
    rungs
