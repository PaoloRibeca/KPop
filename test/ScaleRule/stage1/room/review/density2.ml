(* density2 <twisted.txt> <inertia.txt> <labels> <tag> <detect3.cand>

   Review of density.ml.  Reproduces its ID-exponent density (same queries, same seed, same
   leader pass, same finest usable kept valley) and adds, per rung and partition:
   - base = member-weighted mean log(n_c / N) and L = member-weighted mean log(R_c / R_all), over
     the clusters density.ml counts, so that dens = base - e * L for any exponent e;
   - W/T, the pooled within-cluster sum of squares over the total;
   - dens_nn: singletons included, with R_c = half the distance to their nearest neighbour;
   - dens_ge5: clusters with fewer than 5 members left out;
   - the two-nearest-neighbour dimension on random subsamples of n/8 and n/32 points, which
     measures it at larger scales, and the density with the n/32 value as exponent;
   - the valley partition under three shuffled leader orders as well as index order. *)

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
  (* Full-set dimension, exactly as density.ml *)
  let st = Random.State.make [| 20260915 |] in
  let perm = Array.init n Fun.id in
  for i = n - 1 downto 1 do
    let j = Random.State.int st (i + 1) in
    let t = perm.(i) in
    perm.(i) <- perm.(j);
    perm.(j) <- t
  done;
  let nq = min 1000 n in
  let id_tot = Array.make nr 0. and id_cnt = Array.make nr 0 and skipped = Array.make nr 0 in
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
        end else skipped.(ri) <- skipped.(ri) + 1)
      rungs
  done;
  let idim = Array.init nr (fun ri -> float_of_int id_cnt.(ri) /. id_tot.(ri)) in
  (* Dimension on subsamples of m points, all of them queries *)
  let decimated m nsub seed =
    let st = Random.State.make [| seed |] in
    let tot = Array.make nr 0. and cnt = Array.make nr 0 in
    for _ = 1 to nsub do
      let p = Array.init n Fun.id in
      for i = 0 to m - 1 do
        let j = i + Random.State.int st (n - i) in
        let t = p.(i) in
        p.(i) <- p.(j);
        p.(j) <- t
      done;
      let s = Array.sub p 0 m in
      let dm = Array.make_matrix m m 0. in
      let prev = ref 0 in
      Array.iteri
        (fun ri d ->
          for a = !prev to d - 1 do
            for x = 0 to m - 1 do
              let px = pc.(s.(x)).(a) in
              for y = x + 1 to m - 1 do
                let v = px -. pc.(s.(y)).(a) in
                dm.(x).(y) <- dm.(x).(y) +. (v *. v)
              done
            done
          done;
          prev := d;
          for x = 0 to m - 1 do
            let b1 = ref infinity and b2 = ref infinity in
            for y = 0 to m - 1 do
              if y <> x then begin
                let v = if x < y then dm.(x).(y) else dm.(y).(x) in
                if v < !b1 then begin b2 := !b1; b1 := v end else if v < !b2 then b2 := v
              end
            done;
            if !b1 > 0. then begin
              tot.(ri) <- tot.(ri) +. (0.5 *. log (!b2 /. !b1));
              cnt.(ri) <- cnt.(ri) + 1
            end
          done)
        rungs
    done;
    Array.init nr (fun ri -> if tot.(ri) > 0. then float_of_int cnt.(ri) /. tot.(ri) else nan) in
  let id8 = decimated (n / 8) 8 1 and id32 = decimated (n / 32) 32 2 in
  let stats d assign =
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
    let ss = Array.make ncl 0. and ss_all = ref 0. and first = Array.make ncl (-1) in
    Array.iteri
      (fun i c ->
        if first.(c) < 0 then first.(c) <- i;
        for a = 0 to d - 1 do
          let y = pc.(i).(a) -. cent.(c).(a) and z = pc.(i).(a) -. all.(a) in
          ss.(c) <- ss.(c) +. (y *. y);
          ss_all := !ss_all +. (z *. z)
        done)
      assign;
    let r_all = sqrt (!ss_all /. float_of_int n) in
    let nn i =
      let best = ref infinity in
      for j = 0 to n - 1 do
        if j <> i then begin
          let s = ref 0. and a = ref 0 in
          while !a < d && !s < !best do
            let y = pc.(i).(!a) -. pc.(j).(!a) in
            s := !s +. (y *. y);
            incr a
          done;
          if !s < !best then best := !s
        end
      done;
      sqrt !best in
    let live = ref 0 and singletons = ref 0 in
    let c0 = ref 0 and b0 = ref 0. and l0 = ref 0. in
    let c1 = ref 0 and b1 = ref 0. and l1 = ref 0. in
    let c5 = ref 0 and b5 = ref 0. and l5 = ref 0. in
    let ss_within = ref 0. in
    for c = 0 to ncl - 1 do
      if size.(c) > 0 then begin
        incr live;
        ss_within := !ss_within +. ss.(c);
        let m = size.(c) in
        let fm = float_of_int m in
        let base = log (fm /. float_of_int n) in
        if m = 1 then begin
          incr singletons;
          let r = nn first.(c) /. 2. in
          if r > 0. then begin
            incr c1;
            b1 := !b1 +. base;
            l1 := !l1 +. log (r /. r_all)
          end
        end else if ss.(c) > 0. then begin
          let lr = log (sqrt (ss.(c) /. fm) /. r_all) in
          c0 := !c0 + m; b0 := !b0 +. (fm *. base); l0 := !l0 +. (fm *. lr);
          c1 := !c1 + m; b1 := !b1 +. (fm *. base); l1 := !l1 +. (fm *. lr);
          if m >= 5 then begin c5 := !c5 + m; b5 := !b5 +. (fm *. base); l5 := !l5 +. (fm *. lr) end
        end
      end
    done;
    let per s c = if !c > 0 then !s /. float_of_int !c else nan in
    ( !live, !singletons, !ss_within /. !ss_all,
      (per b0 c0, per l0 c0), (per b1 c1, per l1 c1), (per b5 c5, per l5 c5) ) in
  let leader d radius order =
    let assign = Array.make n (-1) and leaders = ref [] and nl = ref 0 in
    let r2 = radius *. radius in
    Array.iter
      (fun i ->
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
        else begin leaders := i :: !leaders; assign.(i) <- !nl; incr nl end)
      order;
    assign in
  let shuffled seed =
    let st = Random.State.make [| seed |] in
    let p = Array.init n Fun.id in
    for i = n - 1 downto 1 do
      let j = Random.State.int st (i + 1) in
      let t = p.(i) in
      p.(i) <- p.(j);
      p.(j) <- t
    done;
    p in
  let orders = [ "index", Array.init n Fun.id; "shuf1", shuffled 101; "shuf2", shuffled 102;
                 "shuf3", shuffled 103 ] in
  Printf.printf "# %s: %d labelled spectra, %d axes available\n" tag n avail;
  Printf.printf
    "tag\trung\tpartition\torder\tshare\tclusters\tsingletons\tid\tid_skipped\tid_n8\tid_n32\tWT\tbase\tL\tdens_id\tdens_id_nn\tdens_id_ge5\tdens_id_n32\tdens_id_n32_nn\n%!";
  Array.iteri
    (fun ri d ->
      let line kind order share (live, singletons, wt, (b0, l0), (b1, l1), (b5, l5)) =
        Printf.printf "%s\t%d\t%s\t%s\t%s\t%d\t%d\t%.2f\t%d\t%.2f\t%.2f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n%!"
          tag d kind order share live singletons idim.(ri) skipped.(ri) id8.(ri) id32.(ri) wt b0 l0
          (b0 -. (idim.(ri) *. l0)) (b1 -. (idim.(ri) *. l1)) (b5 -. (idim.(ri) *. l5))
          (b0 -. (id32.(ri) *. l0)) (b1 -. (id32.(ri) *. l1)) in
      line "labels" "-" "-" (stats d lab);
      match Hashtbl.find_opt finest d with
      | None -> ()
      | Some (share, radius) ->
        List.iter
          (fun (oname, order) ->
            line "valley" oname (Printf.sprintf "%.4f" share) (stats d (leader d radius order)))
          orders)
    rungs
