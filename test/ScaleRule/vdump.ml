(* vdump <twisted.txt> <inertia.txt> <tag> <d,d,...> <level,names> <cutspec> <labels> [<labels>]
   One JSON object per (d, scale) with everything needed to draw the pair-distance histogram the
   structure test reads: raw counts, counts of same-class pairs per label level, the +-5 and +-2
   moving averages for the drawn curve, the stock detector's own pick, and -- at the radii the
   CALIBRATED detector found on ALL pairs, passed in as "d=r,r;d=r" -- the share of pairs below each
   and the share of each label level it catches.  This program detects nothing: one detector per page. *)
let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and tag = Sys.argv.(3)
  and ds = List.map int_of_string (String.split_on_char ',' Sys.argv.(4))
  and lnames = String.split_on_char ',' Sys.argv.(5) and cutspec = Sys.argv.(6) in
  let nlab = Array.length Sys.argv - 7 in
  let labs = Array.init nlab (fun i -> Kio.read_labels Sys.argv.(7 + i)) in
  (* The radii for a rung, as the calibrated detector reported them *)
  let pair_cuts d =
    if cutspec = "-" then []
    else
      String.split_on_char ';' cutspec
      |> List.filter_map (fun part ->
             match String.split_on_char '=' part with
             | [ dd; rs ] when int_of_string dd = d && rs <> "-" && rs <> "" ->
               Some (List.map float_of_string (String.split_on_char ',' rs))
             | _ -> None)
      |> List.concat in
  let iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs.(0) names_all.(i)) |> Array.of_list in
  let n = Array.length rows and maxd = Array.length iv in
  let state = Random.State.make [| 17 |] in
  let pairs = 100000 in
  let pa = Array.make pairs 0 and pb = Array.make pairs 0 and cnt = ref 0 in
  for _ = 1 to pairs do
    let a = Random.State.int state n and b = Random.State.int state n in
    if a <> b then begin pa.(!cnt) <- a; pb.(!cnt) <- b; incr cnt end
  done;
  let m = !cnt in
  let label lv r = try Hashtbl.find labs.(lv) names_all.(rows.(r)) with Not_found -> "\000" in
  let same = Array.init nlab (fun lv -> Array.init m (fun p -> label lv pa.(p) = label lv pb.(p))) in
  let base = Array.map (Array.fold_left (fun acc x -> if x then acc + 1 else acc) 0) same in
  let bins = 200 in
  let sq = Array.make m 0. and done_axes = ref 0 in
  let jarr f a = "[" ^ String.concat "," (Array.to_list (Array.map f a)) ^ "]" in
  List.iter
    (fun d ->
      let d = min d maxd in
      for k = !done_axes to d - 1 do
        let wk = iv.(k) in
        for p = 0 to m - 1 do
          let x = coords_all.(rows.(pa.(p))).(k) -. coords_all.(rows.(pb.(p))).(k) in
          sq.(p) <- sq.(p) +. wk *. x *. x
        done
      done;
      done_axes := d;
      let tot_inertia = Array.fold_left ( +. ) 0. iv and kept = ref 0. in
      for k = 0 to d - 1 do kept := !kept +. iv.(k) done;
      let dist = Array.map sqrt sq in
      let sorted = Array.copy dist in
      Array.sort compare sorted;
      let zeros = ref 0 in
      Array.iter (fun x -> if x <= 0. then incr zeros) sorted;
      let hi = sorted.(int_of_float (float m *. 0.99)) in
      let lo = sorted.(!zeros + int_of_float (float (m - !zeros) *. 0.001)) in
      List.iter
        (fun scale ->
          let to_axis, of_axis, a0, a1 = if scale = "log" then log, exp, log lo, log hi else Fun.id, Fun.id, 0., hi in
          let w = (a1 -. a0) /. float bins in
          let bin_of x = if scale = "log" && x <= 0. then -1 else int_of_float ((to_axis x -. a0) /. w) in
          let hist = Array.make bins 0 and hs = Array.init nlab (fun _ -> Array.make bins 0) in
          Array.iteri
            (fun p x ->
              let b = bin_of x in
              if b >= 0 && b < bins then begin
                hist.(b) <- hist.(b) + 1;
                for lv = 0 to nlab - 1 do if same.(lv).(p) then hs.(lv).(b) <- hs.(lv).(b) + 1 done
              end)
            dist;
          let smooth s =
            Array.init bins (fun b ->
              let t = ref 0 and c = ref 0 in
              for j = b - s to b + s do if j >= 0 && j < bins then begin t := !t + hist.(j); incr c end done;
              float !t /. float !c) in
          let describe_valleys s sm =
            let minb a b = let v = ref sm.(min a b) in for j = min a b to max a b do if sm.(j) < !v then v := sm.(j) done; !v in
            let argm a b = let v = ref (min a b) in for j = min a b to max a b do if sm.(j) < sm.(!v) then v := j done; !v in
            let stats r =
              let nb = ref 0 and sb = Array.make nlab 0 in
              for p = 0 to m - 1 do
                if dist.(p) < r then begin
                  incr nb; for lv = 0 to nlab - 1 do if same.(lv).(p) then sb.(lv) <- sb.(lv) + 1 done
                end
              done;
              Printf.sprintf "\"r\":%.6g,\"below\":%.4f,\"caught\":%s" r (float !nb /. float m)
                (jarr (fun lv -> Printf.sprintf "%.4f" (float sb.(lv) /. float base.(lv))) (Array.init nlab Fun.id)) in
            (* The radii come from the calibrated detector, which read every pair; this only
               says where each one falls on the sampled histogram drawn here *)
            let vs rs =
              List.map
                (fun r ->
                  let b = int_of_float ((to_axis r -. a0) /. w) in
                  Printf.sprintf "{\"bin\":%d,%s}" (if b < 0 then 0 else if b >= bins then bins - 1 else b) (stats r))
                rs in
            let stock =
              if scale <> "lin" || s <> 5 then "null"
              else begin
                let top = ref 0 in
                Array.iteri (fun b v -> if v > sm.(!top) then top := b) sm;
                let second = ref (-1) in
                for b = 0 to bins - 1 do
                  if abs (b - !top) > 2 * s && minb b !top <= 0.5 *. sm.(b) && (!second < 0 || sm.(b) > sm.(!second)) then second := b
                done;
                if !second < 0 then "null"
                else begin
                  let lo_b = min !second !top and hi_b = max !second !top in
                  let vb = argm lo_b hi_b and lower = Float.min sm.(lo_b) sm.(hi_b) in
                  let lo_e = ref vb and hi_e = ref vb in
                  while !lo_e > lo_b && sm.(!lo_e - 1) <= 0.05 *. lower do decr lo_e done;
                  while !hi_e < hi_b && sm.(!hi_e + 1) <= 0.05 *. lower do incr hi_e done;
                  Printf.sprintf "{\"bin\":%d,%s}" vb (stats ((float (!lo_e + !hi_e) /. 2. +. 0.5) *. w))
                end
              end in
            Printf.sprintf "\"modes%d\":%s,\"valleys%d\":[%s]%s" s (jarr string_of_int [||]) s
              (String.concat "," (vs (pair_cuts d))) (if scale = "lin" && s = 5 then ",\"stock\":" ^ stock else "") in
          let sm5 = smooth 5 and sm2 = smooth 2 in
          Printf.printf
            "{\"tag\":%S,\"d\":%d,\"axes\":%d,\"inertia_kept\":%.4f,\"scale\":%S,\"a0\":%.6g,\"w\":%.6g,\"pairs\":%d,\"levels\":%s,\"base\":%s,\"counts\":%s,\"same\":%s,\"smooth5\":%s,\"smooth2\":%s,%s,%s}\n"
            tag d maxd (!kept /. tot_inertia) scale a0 w m
            (jarr (Printf.sprintf "%S") (Array.of_list lnames))
            (jarr (fun b -> Printf.sprintf "%.4f" (float b /. float m)) base)
            (jarr string_of_int hist)
            (jarr (fun h -> jarr string_of_int h) hs)
            (jarr (Printf.sprintf "%.1f") sm5) (jarr (Printf.sprintf "%.1f") sm2)
            (describe_valleys 5 sm5) (describe_valleys 2 sm2))
        [ "lin"; "log" ])
    ds
