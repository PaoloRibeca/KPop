(* agree <twisted.txt> <inertia.txt> <labels> <tag> <detect3.cand>

   Whether the partition stops changing as the axes double.  At 1, 2, 4, ... axes and at all of
   them, the labelled spectra are clustered by leaders at the rung's finest usable valley (kept,
   share <= 0.5, radius from detect3's candidate file), the partition a search would start from.
   For each rung it reports that partition's agreement with the partitions at the next rungs up --
   the adjusted Rand index, and homogeneity and completeness of the finer partition against the
   coarser -- and, to judge the criterion and not as part of it, its adjusted Rand index,
   homogeneity and completeness against the labels.  Coordinates are weighted by the square root
   of the inertia, as in detect3, whose radii are in the same units. *)

(* Adjusted Rand index, homogeneity and completeness of partition a judged against reference b, as
   Clustering.compare_partitions computes them *)
let compare_partitions a b =
  let n = Array.length a in
  let cell = Hashtbl.create 256 and ca = Hashtbl.create 64 and cb = Hashtbl.create 64 in
  let bump t k = Hashtbl.replace t k (1 + Option.value ~default:0 (Hashtbl.find_opt t k)) in
  for i = 0 to n - 1 do
    bump cell (a.(i), b.(i));
    bump ca a.(i);
    bump cb b.(i)
  done;
  let choose2 x = let x = float_of_int x in x *. (x -. 1.) /. 2. in
  let sum t f = Hashtbl.fold (fun _ v acc -> acc +. f v) t 0. in
  let index = sum cell choose2 and sa = sum ca choose2 and sb = sum cb choose2
  and total = choose2 n and fn = float_of_int n in
  let expected = if total > 0. then sa *. sb /. total else 0. and maximum = (sa +. sb) /. 2. in
  let ari = if maximum -. expected <> 0. then (index -. expected) /. (maximum -. expected) else 1. in
  let log2 x = log x /. log 2. in
  let entropy t = sum t (fun v -> let p = float_of_int v /. fn in -.(p *. log2 p)) in
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
  ari, homogeneity, completeness

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
  let parts =
    Array.map
      (fun d ->
        match Hashtbl.find_opt finest d with
        | None -> None
        | Some (share, radius) -> let assign, k = leader d radius in Some (share, assign, k))
      rungs in
  Printf.printf "# %s: %d labelled spectra, %d axes available\n" tag n avail;
  Printf.printf
    "tag\trung\tvalley_share\tclusters\tari_next\thom_next\tcomp_next\tari_next2\tari_labels\thom_labels\tcomp_labels\n";
  Array.iteri
    (fun ri d ->
      match parts.(ri) with
      | None -> Printf.printf "%s\t%d\tnone\t-\t-\t-\t-\t-\t-\t-\t-\n%!" tag d
      | Some (share, assign, k) ->
        let against j =
          if j < nr then
            match parts.(j) with
            | Some (_, other, _) -> Some (compare_partitions other assign)
            | None -> None
          else None in
        let show = function Some (x, _, _) -> Printf.sprintf "%.3f" x | None -> "-" in
        let next = against (ri + 1) in
        let hn, cn = match next with Some (_, h, c) -> Printf.sprintf "%.3f" h, Printf.sprintf "%.3f" c | None -> "-", "-" in
        let al, hl, cl = compare_partitions assign lab in
        Printf.printf "%s\t%d\t%.4f\t%d\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\n%!" tag d share k (show next) hn cn
          (show (against (ri + 2))) al hl cl)
    rungs
