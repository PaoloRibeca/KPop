(* splitcmp <twisted.txt> <inertia.txt> <d> <tag> <radii|-> <tree.nwk> <level,names> <labels> [<labels>]
   Are the partitions got by cutting the average-linkage tree at the valleys clades of another
   tree?  Builds the UPGMA dendrogram of the labelled spectra over the first d principal axes,
   cuts it at each radius given, and asks of every resulting cluster (two or more members) how
   closely it matches a clade of the Newick tree -- best Jaccard over both sides of every branch,
   and whether it matches one exactly.  Then the same of each label class, in both trees, and the
   Robinson-Foulds distance between the two trees over the labelled leaves.  A round trip of the
   UPGMA tree through its own Newick checks the parser and the split code: it must give RF 0 and
   every cluster an exact clade *)

(* A tree: children per node, and the labelled index of each leaf (-1 for internal nodes and for
   leaves the labels do not name) *)
type tree = { children: int array array; leaf: int array; root: int }

let postorder t =
  let n = Array.length t.children in
  let out = Array.make n 0 and k = ref 0 in
  let stack = Stack.create () in
  Stack.push (t.root, false) stack;
  while not (Stack.is_empty stack) do
    let v, expanded = Stack.pop stack in
    if expanded || Array.length t.children.(v) = 0 then begin out.(!k) <- v; incr k end
    else begin
      Stack.push (v, true) stack;
      Array.iter (fun c -> Stack.push (c, false) stack) t.children.(v)
    end
  done;
  Array.sub out 0 !k

let parse_newick s names_to_idx =
  let len = String.length s and pos = ref 0 in
  let kids = ref [] and leaves = ref [] and count = ref 0 in
  let rec skip () =
    while !pos < len && (match s.[!pos] with ' ' | '\t' | '\n' | '\r' -> true | _ -> false) do incr pos done;
    if !pos < len && s.[!pos] = '[' then begin
      while !pos < len && s.[!pos] <> ']' do incr pos done;
      incr pos;
      skip ()
    end in
  let read_label () =
    skip ();
    if !pos < len && (s.[!pos] = '\'' || s.[!pos] = '"') then begin
      let q = s.[!pos] in
      incr pos;
      let b = Buffer.create 32 in
      let fin = ref false in
      while not !fin do
        if !pos >= len then fin := true
        else if s.[!pos] = q && !pos + 1 < len && s.[!pos + 1] = q then begin Buffer.add_char b q; pos := !pos + 2 end
        else if s.[!pos] = q then begin incr pos; fin := true end
        else begin Buffer.add_char b s.[!pos]; incr pos end
      done;
      Buffer.contents b
    end else begin
      let st = !pos in
      while !pos < len && not (String.contains ",():;[ \t\n\r" s.[!pos]) do incr pos done;
      String.sub s st (!pos - st)
    end in
  let skip_length () =
    skip ();
    if !pos < len && s.[!pos] = ':' then begin
      incr pos;
      skip ();
      while !pos < len && not (String.contains ",();[ \t\n\r" s.[!pos]) do incr pos done
    end in
  let rec subtree () =
    skip ();
    let id = !count in
    incr count;
    if !pos < len && s.[!pos] = '(' then begin
      incr pos;
      let ch = ref [] and go = ref true in
      while !go do
        ch := subtree () :: !ch;
        skip ();
        if !pos < len && s.[!pos] = ',' then incr pos
        else if !pos < len && s.[!pos] = ')' then begin incr pos; go := false end
        else failwith (Printf.sprintf "Newick: unexpected character at offset %d" !pos)
      done;
      kids := (id, Array.of_list (List.rev !ch)) :: !kids;
      ignore (read_label ())
    end else begin
      let nm = read_label () in
      leaves := (id, Option.value ~default:(-1) (Hashtbl.find_opt names_to_idx nm)) :: !leaves
    end;
    skip_length ();
    id in
  let root = subtree () in
  let children = Array.make !count [||] and leaf = Array.make !count (-1) in
  List.iter (fun (id, ch) -> children.(id) <- ch) !kids;
  List.iter (fun (id, l) -> leaf.(id) <- l) !leaves;
  { children; leaf; root }

let () =
  let tw = Sys.argv.(1) and inf = Sys.argv.(2) and d = int_of_string Sys.argv.(3) and tag = Sys.argv.(4) in
  let radii = if Sys.argv.(5) = "-" then [] else List.map float_of_string (String.split_on_char ',' Sys.argv.(5)) in
  let nwk = Sys.argv.(6) and lnames = Array.of_list (String.split_on_char ',' Sys.argv.(7)) in
  let nlab = Array.length Sys.argv - 8 in
  let labs = Array.init nlab (fun i -> Kio.read_labels Sys.argv.(8 + i)) in
  let iv = Kio.read_inertia inf in
  let names_all, coords_all = Kio.read_twisted tw in
  let rows =
    List.init (Array.length names_all) Fun.id
    |> List.filter (fun i -> Hashtbl.mem labs.(0) names_all.(i)) |> Array.of_list in
  let n = Array.length rows and d = min d (Array.length iv) in
  let names = Array.map (fun r -> names_all.(r)) rows in
  let names_to_idx = Hashtbl.create (2 * n) in
  Array.iteri (fun i nm -> Hashtbl.replace names_to_idx nm i) names;
  let lab =
    Array.map
      (fun h ->
        let tbl = Hashtbl.create 64 in
        Array.map
          (fun nm ->
            let s = Hashtbl.find h nm in
            match Hashtbl.find_opt tbl s with
            | Some id -> id
            | None -> let id = Hashtbl.length tbl in Hashtbl.add tbl s id; id)
          names)
      labs in
  (* UPGMA over the labelled spectra *)
  let pc = Array.map (fun r -> Array.init d (fun k -> coords_all.(r).(k) *. sqrt iv.(k))) rows in
  let npairs = n * (n - 1) / 2 in
  let idx i j = let i, j = if i < j then i, j else j, i in (i * (2 * n - i - 1) / 2) + (j - i - 1) in
  let dm = Float.Array.make npairs 0. in
  for i = 0 to n - 2 do
    let pi = pc.(i) and base = (i * (2 * n - i - 1) / 2) - i - 1 in
    for j = i + 1 to n - 1 do
      let pj = pc.(j) and s = ref 0. in
      for k = 0 to d - 1 do let x = pi.(k) -. pj.(k) in s := !s +. (x *. x) done;
      Float.Array.unsafe_set dm (base + j) (sqrt !s)
    done
  done;
  let active = Array.make n true and size = Array.make n 1 in
  let chain = Array.make (n + 1) 0 and top = ref 0 and start = ref 0 and nm = ref 0 in
  let mi = Array.make (n - 1) 0 and mj = Array.make (n - 1) 0 and mh = Array.make (n - 1) 0. in
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
      mi.(!nm) <- i; mj.(!nm) <- j; mh.(!nm) <- !bd; incr nm;
      for k = 0 to n - 1 do
        if active.(k) && k <> i && k <> j then
          Float.Array.set dm (idx i k)
            (((si *. Float.Array.get dm (idx i k)) +. (sj *. Float.Array.get dm (idx j k))) /. (si +. sj))
      done;
      size.(i) <- size.(i) + size.(j);
      active.(j) <- false
    end else begin chain.(!top) <- !best; incr top end
  done;
  (* The UPGMA tree as a generic tree: leaves 0..n-1, merge t is node n+t *)
  let upgma =
    let slot = Array.init n Fun.id and children = Array.make ((2 * n) - 1) [||] in
    let order = Array.init (n - 1) Fun.id in
    Array.sort (fun a b -> compare mh.(a) mh.(b)) order;
    (* The chain emits merges out of height order but each slot's current node must be the one
       formed last, so nodes are numbered in the order the chain made them *)
    for t = 0 to n - 2 do
      let v = n + t in
      children.(v) <- [| slot.(mi.(t)); slot.(mj.(t)) |];
      slot.(mi.(t)) <- v
    done;
    ignore order;
    { children; leaf = Array.init ((2 * n) - 1) (fun v -> if v < n then v else -1); root = (2 * n) - 2 } in
  let height v = if v < n then 0. else mh.(v - n) in
  (* Clusters from cutting at a height: union-find over merges no higher than it *)
  let cut h =
    let parent = Array.init n Fun.id in
    let rec find x = if parent.(x) = x then x else begin let r = find parent.(x) in parent.(x) <- r; r end in
    for t = 0 to n - 2 do
      if mh.(t) <= h then begin
        let a = find mi.(t) and b = find mj.(t) in
        if a <> b then parent.(a) <- b
      end
    done;
    Array.init n find in
  (* For a partition of the labelled leaves, how its classes sit in a tree *)
  let clade_match t post part =
    let nnodes = Array.length t.children in
    let lc = Array.make nnodes 0 in
    Array.iter (fun v -> if t.leaf.(v) >= 0 then lc.(v) <- 1 else Array.iter (fun c -> lc.(v) <- lc.(v) + lc.(c)) t.children.(v)) post;
    let total = lc.(t.root) in
    let members = Hashtbl.create 256 in
    Array.iteri (fun i c -> Hashtbl.replace members c (i :: Option.value ~default:[] (Hashtbl.find_opt members c))) part;
    let ind = Array.make n false and cnt = Array.make nnodes 0 in
    let nclust = ref 0 and nexact = ref 0 and wexact = ref 0 and wj = ref 0. and wtot = ref 0 in
    Hashtbl.iter
      (fun _ ms ->
        let c = List.length ms in
        if c >= 2 then begin
          List.iter (fun i -> ind.(i) <- true) ms;
          let bestj = ref 0. and exact = ref false in
          Array.iter
            (fun v ->
              cnt.(v) <-
                (if t.leaf.(v) >= 0 then (if ind.(t.leaf.(v)) then 1 else 0)
                 else Array.fold_left (fun acc ch -> acc + cnt.(ch)) 0 t.children.(v));
              if v <> t.root then begin
                let inter = cnt.(v) and sz = lc.(v) in
                let j1 = float inter /. float (c + sz - inter) in
                let inter2 = c - inter and sz2 = total - sz in
                let j2 = if c + sz2 - inter2 > 0 then float inter2 /. float (c + sz2 - inter2) else 0. in
                if j1 > !bestj then bestj := j1;
                if j2 > !bestj then bestj := j2;
                if (inter = c && sz = c) || (inter2 = c && sz2 = c) then exact := true
              end)
            post;
          List.iter (fun i -> ind.(i) <- false) ms;
          incr nclust;
          if !exact then begin incr nexact; wexact := !wexact + c end;
          wj := !wj +. (float c *. !bestj);
          wtot := !wtot + c
        end)
      members;
    !nclust, !nexact, (if !wtot > 0 then float !wexact /. float !wtot else nan), (if !wtot > 0 then !wj /. float !wtot else nan) in
  (* Non-trivial splits over the labelled leaves, as canonical bit strings *)
  let splits t post =
    let nb = (n + 7) / 8 in
    let sets = Array.make (Array.length t.children) Bytes.empty and lc = Array.make (Array.length t.children) 0 in
    let tbl = Hashtbl.create 8192 in
    Array.iter
      (fun v ->
        let b = Bytes.make nb '\000' in
        if t.leaf.(v) >= 0 then begin
          let i = t.leaf.(v) in
          Bytes.set b (i lsr 3) (Char.chr (1 lsl (i land 7)));
          lc.(v) <- 1
        end else
          Array.iter
            (fun c ->
              for k = 0 to nb - 1 do
                Bytes.set b k (Char.chr (Char.code (Bytes.get b k) lor Char.code (Bytes.get sets.(c) k)))
              done;
              lc.(v) <- lc.(v) + lc.(c))
            t.children.(v);
        sets.(v) <- b;
        let small = min lc.(v) (n - lc.(v)) in
        if v <> t.root && small >= 2 then begin
          let key =
            if Char.code (Bytes.get b 0) land 1 = 1 then begin
              let c = Bytes.copy b in
              for k = 0 to nb - 1 do Bytes.set c k (Char.chr (lnot (Char.code (Bytes.get c k)) land 255)) done;
              let extra = (nb * 8) - n in
              if extra > 0 then Bytes.set c (nb - 1) (Char.chr (Char.code (Bytes.get c (nb - 1)) land (255 lsr extra)));
              Bytes.to_string c
            end else Bytes.to_string b in
          Hashtbl.replace tbl key small
        end)
      post;
    Array.iteri (fun i _ -> if i < Array.length sets then sets.(i) <- Bytes.empty) sets;
    tbl in
  let compare_splits a b =
    let shared = ref 0 and deep_a = ref 0 and deep_shared = ref 0 in
    Hashtbl.iter
      (fun k small ->
        let present = Hashtbl.mem b k in
        if present then incr shared;
        if small >= 10 then begin incr deep_a; if present then incr deep_shared end)
      a;
    let na = Hashtbl.length a and nb = Hashtbl.length b in
    na, nb, !shared,
    (if na + nb > 0 then 1. -. (2. *. float !shared /. float (na + nb)) else 0.),
    !deep_a, !deep_shared in
  (* The UPGMA tree through its own Newick, as a check on the parser and the split code *)
  let buf = Buffer.create (64 * n) in
  let rec write v parent_h =
    if v < n then Buffer.add_string buf (Printf.sprintf "'%s':%.6g" names.(v) (parent_h -. 0.))
    else begin
      Buffer.add_char buf '(';
      Array.iteri (fun k c -> if k > 0 then Buffer.add_char buf ','; write c (height v)) upgma.children.(v);
      Buffer.add_string buf (Printf.sprintf "):%.6g" (parent_h -. height v))
    end in
  write upgma.root (height upgma.root);
  Buffer.add_char buf ';';
  let upgma_nwk = Buffer.contents buf in
  let oc = open_out (tag ^ ".d" ^ string_of_int d ^ ".upgma.nwk") in
  output_string oc upgma_nwk;
  close_out oc;
  let back = parse_newick upgma_nwk names_to_idx in
  let post_u = postorder upgma and post_b = postorder back in
  let su = splits upgma post_u and sb = splits back post_b in
  let _, _, _, rf_self, _, _ = compare_splits su sb in
  let nj =
    let ic = open_in_bin nwk in
    let s = really_input_string ic (in_channel_length ic) in
    close_in ic;
    parse_newick s names_to_idx in
  let post_n = postorder nj in
  let nj_leaves = Array.fold_left (fun acc l -> if l >= 0 then acc + 1 else acc) 0 nj.leaf in
  let nj_all_leaves = Array.fold_left (fun acc ch -> if Array.length ch = 0 then acc + 1 else acc) 0 nj.children in
  let sn = splits nj post_n in
  Printf.printf "# %s, %d axes: %d labelled spectra; %s has %d leaves, %d of them labelled\n" tag d n
    (Filename.basename nwk) nj_all_leaves nj_leaves;
  let na, nb, shared, nrf, deep, deep_shared = compare_splits su sn in
  Printf.printf "  check: UPGMA vs its own Newick RF %.4f (must be 0)\n" rf_self;
  Printf.printf "  UPGMA vs tree: %d and %d non-trivial splits, %d shared, normalised RF %.3f; splits with >= 10 leaves on each side: %d of %d shared (%.1f%%)\n"
    na nb shared nrf deep_shared deep (if deep > 0 then 100. *. float deep_shared /. float deep else nan);
  let jrows = ref [] in
  let describe what jhead part =
    let c1, e1, w1, j1 = clade_match nj post_n part in
    let _, e2, w2, j2 = clade_match back post_b part in
    Printf.printf "  %-34s %4d groups | in KPopPhylo tree: %4d exact (%.1f%% of members), best-match Jaccard %.3f | in UPGMA: %4d exact (%.1f%%), Jaccard %.3f\n"
      what c1 e1 (100. *. w1) j1 e2 (100. *. w2) j2;
    jrows :=
      Printf.sprintf "{%s,\"groups\":%d,\"nj_exact\":%d,\"nj_pct\":%.4f,\"nj_j\":%.4f,\"up_exact\":%d,\"up_pct\":%.4f,\"up_j\":%.4f}"
        jhead c1 e1 w1 j1 e2 w2 j2
      :: !jrows in
  List.iter (fun r -> describe (Printf.sprintf "valley cut at %.4g" r) (Printf.sprintf "\"kind\":\"valley\",\"r\":%.6g" r) (cut r)) radii;
  Array.iteri
    (fun lv l -> describe (Printf.sprintf "CDC %s classes" lnames.(lv)) (Printf.sprintf "\"kind\":\"class\",\"level\":%S" lnames.(lv)) l)
    lab;
  Printf.printf
    "JSON {\"tag\":%S,\"d\":%d,\"n\":%d,\"leaves\":%d,\"labelled\":%d,\"check_rf\":%.4f,\"upgma_splits\":%d,\"nj_splits\":%d,\"shared\":%d,\"nrf\":%.4f,\"deep\":%d,\"deep_shared\":%d,\"rows\":[%s]}\n"
    tag d n nj_all_leaves nj_leaves rf_self na nb shared nrf deep deep_shared (String.concat "," (List.rev !jrows));
  print_endline ""
