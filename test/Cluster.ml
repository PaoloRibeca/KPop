(*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*)

(* Greedy leader clustering on KPop twisted (standard) coordinates,
   with optional built-in epsilon estimation.

   Four modes (exactly one required):
     Explicit epsilon (-e <f>): user provides the global absorption threshold
       directly; no estimation is performed.  O(n log n) with HNSW.

     1-NN elbow (-1): sorts the 1-NN distances of all n points (O(n log n)
       with HNSW, O(n^2) with flat),
       identifies the kneedle elbow, outputs the sorted distance table (same
       format as Epsilon Part 1), then clusters with that elbow distance as
       the global epsilon.

     Dist_star elbow (-3): computes the maximum of k/V(d_k,D) for each of
       n_sample randomly chosen starting points (O(n_sample * n)), identifies
       the kneedle elbow in the sorted dist_star values, outputs the table
       (same format as Epsilon Part 3), then clusters with that elbow
       distance as the global epsilon.

     Adaptive (-a): computes dist_star for ALL n points (O(n^2)), outputs the
       full sorted dist_star table (with elbow marked for reference), then
       clusters with each point's own dist_star as its absorption threshold,
       processing points in order of increasing dist_star (densest first).

   Output (TSV, stdout):
     For modes -1, -3, and -a an estimation table is printed first (same
     column layout as Epsilon Parts 1 and 3 respectively), followed by a
     blank line and then the clustering results.  Columns in the clustering
     section:
       name  [kmer]  representative  [rep_kmer]  status
     where status is "rep" (representative) or "abs" (absorbed).

   When -k is given, a DNA sequence column (fwd/revcomp) is added after each
   hash name using the same back-translation as in Epsilon.

   Note: --order is ignored in adaptive mode; points are always processed
   in order of increasing dist_star. *)

open BiOCamLib.Better
open KPop

type mode_t =
  (* -e: user-provided global epsilon *)
  | Explicit of float
  (* -1: use Part-1 kneedle elbow as global epsilon *)
  | Elbow1nn
  (* -3: use Part-3 dist_star elbow as global epsilon *)
  | ElbowDs
  (* -a: per-point dist_star threshold *)
  | Adaptive

module Defaults =
  struct
    let n_sample = 200
    let index_type = "hnsw(32)"
    let order = "random"
    let kmer_len = 0
    let verbose = false
  end

module Parameters =
  struct
    let input = ref ""
    let mode = ref (None: mode_t option)
    let n_sample = ref Defaults.n_sample
    let index_type = ref Defaults.index_type
    let order = ref Defaults.order
    let kmer_len = ref Defaults.kmer_len
    let verbose = ref Defaults.verbose
  end

let set_mode m =
  match !Parameters.mode with
  | Some _ ->
    Printf.eprintf
      "Error: only one of -e, -1, -3, -a may be given\n%!";
    exit 1
  | None ->
    Parameters.mode := Some m

let usage () =
  Printf.eprintf
    "Usage: %s <twisted_prefix> (-e <epsilon> | -1 | -3 | -a) [OPTIONS]\n\
     Modes (exactly one required):\n\
     \  -e, --epsilon <f>        explicit global epsilon\n\
     \  -1, --elbow-1nn          estimate epsilon from 1-NN elbow (O(n log n) with HNSW)\n\
     \  -3, --elbow-ds           estimate epsilon from dist_star elbow (O(n_sample*n))\n\
     \  -a, --adaptive           adaptive per-point dist_star epsilon (O(n^2))\n\
     Options:\n\
     \  -n, --n-sample <n>       random sample size for -3 mode (default: %d)\n\
     \  -I, --index <type>       FAISS index type: flat, hnsw(m) (default: %s)\n\
     \  -r, --order <order>      processing order for global modes: random, file (default: %s)\n\
     \  -k, --kmer <len>         k-mer length for hash back-translation (default: no decoding)\n\
     \  -v, --verbose            verbose output on stderr\n\
     \  -h, --help               print this message and exit\n%!"
    Sys.argv.(0)
    Defaults.n_sample
    Defaults.index_type
    Defaults.order

let () =
  let argc = Array.length Sys.argv in
  if argc < 2 then begin usage (); exit 1 end;
  let i = ref 1 in
  while !i < argc do
    let arg = Sys.argv.(!i) in
    let get_next flag =
      incr i;
      if !i >= argc then begin
        Printf.eprintf "Error: %s requires an argument\n%!" flag;
        usage (); exit 1
      end;
      Sys.argv.(!i) in
    (match arg with
    | "-e" | "--epsilon" ->
      let v = get_next arg in
      (match float_of_string_opt v with
      | Some f when f > 0. -> set_mode (Explicit f)
      | _ ->
        Printf.eprintf "Error: --epsilon must be a positive float (got '%s')\n%!" v;
        usage (); exit 1)
    | "-1" | "--elbow-1nn" ->
      set_mode Elbow1nn
    | "-3" | "--elbow-ds" ->
      set_mode ElbowDs
    | "-a" | "--adaptive" ->
      set_mode Adaptive
    | "-n" | "--n-sample" ->
      let v = get_next arg in
      (match int_of_string_opt v with
      | Some n when n > 0 -> Parameters.n_sample := n
      | _ ->
        Printf.eprintf "Error: --n-sample must be a positive integer (got '%s')\n%!" v;
        usage (); exit 1)
    | "-I" | "--index" ->
      let v = get_next arg in
      (try ignore (Interfaiss.Type.of_string v)
       with _ ->
         Printf.eprintf "Error: unknown FAISS index type '%s'\n%!" v;
         usage (); exit 1);
      Parameters.index_type := v
    | "-r" | "--order" ->
      let v = get_next arg in
      (match v with
      | "random" | "file" -> Parameters.order := v
      | _ ->
        Printf.eprintf "Error: --order must be random or file (got '%s')\n%!" v;
        usage (); exit 1)
    | "-k" | "--kmer" ->
      let v = get_next arg in
      (match int_of_string_opt v with
      | Some k when k > 0 -> Parameters.kmer_len := k
      | _ ->
        Printf.eprintf "Error: --kmer must be a positive integer (got '%s')\n%!" v;
        usage (); exit 1)
    | "-v" | "--verbose" ->
      Parameters.verbose := true
    | "-h" | "--help" ->
      usage (); exit 0
    | s when s.[0] = '-' ->
      Printf.eprintf "Error: unknown option '%s'\n%!" s;
      usage (); exit 1
    | s ->
      if !Parameters.input = "" then
        Parameters.input := s
      else begin
        Printf.eprintf "Error: unexpected argument '%s'\n%!" s;
        usage (); exit 1
      end);
    incr i
  done;
  if !Parameters.input = "" then begin
    Printf.eprintf "Error: no input twisted prefix specified\n%!";
    usage (); exit 1
  end;
  if !Parameters.mode = None then begin
    Printf.eprintf "Error: one of -e, -1, -3, -a is required\n%!";
    usage (); exit 1
  end

let twisted = Twisted.of_binary ~verbose:!Parameters.verbose !Parameters.input

let () = Random.self_init ()

let () =
  let kmer_len = !Parameters.kmer_len in
  let index_type = Interfaiss.Type.of_string !Parameters.index_type in
  let mode = match !Parameters.mode with Some m -> m | None -> assert false in

  let mat = twisted.Twisted.twisted.Matrix.matrix in
  let data = mat.Matrix.Base.data in
  let names = mat.Matrix.Base.row_names in
  let n = Array.length names in
  let d = if n > 0 then Float.Array.length data.(0) else 0 in
  if n < 2 then begin
    Printf.eprintf "Error: at least 2 points required (found %d)\n%!" n;
    exit 1
  end;
  if d = 0 then begin
    Printf.eprintf "Error: dimensionality is 0\n%!";
    exit 1
  end;

  (* Plain Euclidean distance between two rows *)
  let euclidean a b =
    let s = ref 0. in
    for i = 0 to d - 1 do
      let x = Float.Array.unsafe_get a i -. Float.Array.unsafe_get b i in
      s := !s +. x *. x
    done;
    sqrt !s in

  (* log(Gamma(d/2 + 1)) for D-ball volume (needed by ElbowDs and Adaptive) *)
  let log_gamma_d2_plus1 =
    if d mod 2 = 0 then begin
      let s = ref 0. in
      for i = 1 to d / 2 do s := !s +. log (float_of_int i) done;
      !s
    end else begin
      let m = (d + 1) / 2 in
      let s = ref (0.5 *. log Float.pi) in
      for i = 1 to m do s := !s +. log (float_of_int (2 * i - 1)) done;
      s := !s -. float_of_int m *. log 2.;
      !s
    end in
  let log_vol1 = 0.5 *. float_of_int d *. log Float.pi -. log_gamma_d2_plus1 in
  let fd = float_of_int d in

  (* k-mer back-translation (same as Epsilon) *)
  let bases = [| 'A'; 'C'; 'G'; 'T' |] in
  let comp = [| 'T'; 'G'; 'C'; 'A' |] in
  let decode_segment n_bases n =
    let fwd = String.init n_bases (fun i ->
      let shift = 2 * (n_bases - 1 - i) in
      bases.(Int64.to_int (Int64.logand (Int64.shift_right_logical n shift) 3L))) in
    let rc = String.init n_bases (fun i ->
      let shift = 2 * i in
      comp.(Int64.to_int (Int64.logand (Int64.shift_right_logical n shift) 3L))) in
    fwd, rc in
  let decode_hash hex_str =
    match String.split_on_char '_' hex_str with
    | [ s ] ->
      if kmer_len = 0 then "?/?"
      else
        (match Int64.of_string_opt ("0x" ^ s) with
        | None -> "?/?"
        | Some n ->
          let fwd, rc = decode_segment kmer_len n in
          fwd ^ "/" ^ rc)
    | segs ->
      if kmer_len = 0 then "?/?"
      else
        let n_segs = List.length segs in
        let n_bases = kmer_len / n_segs in
        let results = List.map (fun s ->
          match Int64.of_string_opt ("0x" ^ s) with
          | None -> "?", "?"
          | Some n -> decode_segment n_bases n) segs in
        let fwd_segs = List.map fst results in
        let rc_segs = List.rev_map snd results in
        String.concat "_" fwd_segs ^ "/" ^ String.concat "_" rc_segs in
  let kmer_col name =
    if kmer_len = 0 then ""
    else "\t" ^ decode_hash name in
  let kmer_hdr suffix =
    if kmer_len = 0 then ""
    else "\t" ^ suffix in

  (* Blank line between output sections *)
  let parts_done = ref 0 in
  let section_sep () =
    if !parts_done > 0 then print_char '\n';
    incr parts_done in

  (* FAISS buffer: 1-row float32 bigarray, reused for every query and add *)
  let buf = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout 1 d in
  let load_buf i =
    Float.Array.iteri (fun j x -> buf.{0, j} <- x) data.(i) in

  (* Greedy leader clustering via FAISS.
     order:      original indices in the order they are processed.
     epsilon_of: original index -> absorption threshold.
     Returns (rep_orig, n_reps) where rep_orig.(i) is the original index of
     i's representative (rep_orig.(i) = i iff i is itself a representative). *)
  let greedy_leader order epsilon_of =
    let rep_orig     = Array.make n (-1) in
    let rep_of_faiss = Array.make n (-1) in
    let n_reps = ref 0 in
    let index = Interfaiss.create ~index_type d in
    Array.iter (fun i ->
      if !n_reps > 0 then begin
        load_buf i;
        let offsets, sq_dists = Interfaiss.query index buf 1 in
        let faiss_slot = Int64.to_int offsets.{0, 0} in
        let dist = sqrt sq_dists.{0, 0} in
        if dist < epsilon_of i then
          rep_orig.(i) <- rep_of_faiss.(faiss_slot)
        else begin
          Interfaiss.add index buf;
          rep_of_faiss.(!n_reps) <- i;
          incr n_reps;
          rep_orig.(i) <- i
        end
      end else begin
        load_buf i;
        Interfaiss.add index buf;
        rep_of_faiss.(0) <- i;
        n_reps := 1;
        rep_orig.(i) <- i
      end)
    order;
    Interfaiss.delete index;
    rep_orig, !n_reps in

  (* ----------------------------------------------------------------
     Mode-specific epsilon computation and optional estimation table
     ---------------------------------------------------------------- *)

  (* Compute dist_star for a set of starting points given by orig_indices.
     Returns (name, max_density, k_star, dist_star) for each. *)
  let compute_dist_star orig_indices =
    let m = Array.length orig_indices in
    Array.init m (fun si ->
      let i = orig_indices.(si) in
      let dists = Array.init (n - 1) (fun j ->
        let j' = if j < i then j else j + 1 in
        euclidean data.(i) data.(j')) in
      Array.sort compare dists;
      let best = ref 0. and best_k = ref 0 and best_dist = ref infinity in
      Array.iteri
        (fun idx dist ->
          if dist > 0. then begin
            let density =
              float_of_int (idx + 1) *. exp (-. log_vol1 -. fd *. log dist) in
            if density > !best then begin
              best := density; best_k := idx + 1; best_dist := dist
            end
          end)
        dists;
      if !Parameters.verbose && (si + 1) mod 100 = 0 then
        Printf.eprintf "%s\r(%s): dist_star: done %d/%d points%!"
          String.TermIO.clear Sys.argv.(0) (si + 1) m;
      names.(i), !best, !best_k, !best_dist) in

  (* Print a dist_star table (shared by ElbowDs and Adaptive) *)
  let print_dist_star_table section_title max_densities =
    let m = Array.length max_densities in
    let elbow_rank =
      Clustering.kneedle_elbow (fun i -> let _, _, _, r = max_densities.(i) in r) m in
    section_sep ();
    Printf.printf
      "=== %s ===\n\
       # elbow at rank %d, dist_star %.15g\n\
       rank\tname%s\tmax_density\tk_star\tdist_star\telbow\n"
      section_title
      elbow_rank (let _, _, _, r = max_densities.(elbow_rank) in r)
      (kmer_hdr "kmer");
    Array.iteri
      (fun rank (name, max_d, k_star, dist_star) ->
        Printf.printf "%d\t%s%s\t%.15g\t%d\t%.15g\t%s\n"
          rank name (kmer_col name) max_d k_star dist_star
          (if rank = elbow_rank then "*" else ""))
      max_densities;
    elbow_rank in

  (* Determine epsilon and processing order for greedy leader *)
  let epsilon, order, dist_star_all =
    match mode with

    | Explicit epsilon ->
      let order = Array.init n Fun.id in
      if !Parameters.order = "random" then begin
        for i = n - 1 downto 1 do
          let j = Random.int (i + 1) in
          let tmp = order.(i) in order.(i) <- order.(j); order.(j) <- tmp
        done
      end;
      epsilon, order, None

    | Elbow1nn ->
      (* Build a FAISS index over all n points and query each for its 2
         nearest neighbours.  The first result is always the point itself
         (distance ≈ 0); the second is the true (or approximate) 1-NN.
         Complexity: O(n log n) with HNSW, O(n^2) with flat. *)
      let index_str = Interfaiss.Type.to_string index_type in
      if !Parameters.verbose then
        Printf.eprintf
          "(%s): Computing 1-NN distances for %d points in %d dimensions via FAISS %s...\n%!"
          Sys.argv.(0) n d index_str;
      let data_all = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout n d in
      Array.iteri (fun i row ->
        Float.Array.iteri (fun j x -> data_all.{i, j} <- x) row) data;
      let idx1nn = Interfaiss.create ~index_type d in
      Interfaiss.train idx1nn data_all;
      Interfaiss.add idx1nn data_all;
      let offsets, sq_dists = Interfaiss.query idx1nn data_all 2 in
      Interfaiss.delete idx1nn;
      if !Parameters.verbose then
        Printf.eprintf "(%s): done.\n%!" Sys.argv.(0);
      (* offsets.{i,0} should be i (self); take the second result.
         If the index returns something other than self first (e.g. exact
         duplicate), fall back to the first result. *)
      let sorted_nn1 = Array.init n (fun i ->
        let r0 = Int64.to_int offsets.{i, 0} in
        let dist = sqrt (if r0 = i then sq_dists.{i, 1} else sq_dists.{i, 0}) in
        names.(i), dist) in
      Array.sort (fun (_, d1) (_, d2) -> compare d1 d2) sorted_nn1;
      let elbow_rank = Clustering.kneedle_elbow (fun i -> snd sorted_nn1.(i)) n in
      let epsilon = snd sorted_nn1.(elbow_rank) in
      section_sep ();
      Printf.printf
        "=== Part 1: sorted 1-NN distances (n=%d, D=%d, index=%s) ===\n\
         # elbow at rank %d, distance %.15g\n\
         rank\tname%s\t1nn_distance\telbow\n"
        n d index_str elbow_rank epsilon
        (kmer_hdr "kmer");
      Array.iteri
        (fun rank (name, dist) ->
          Printf.printf "%d\t%s%s\t%.15g\t%s\n"
            rank name (kmer_col name) dist
            (if rank = elbow_rank then "*" else ""))
        sorted_nn1;
      let order = Array.init n Fun.id in
      if !Parameters.order = "random" then begin
        for i = n - 1 downto 1 do
          let j = Random.int (i + 1) in
          let tmp = order.(i) in order.(i) <- order.(j); order.(j) <- tmp
        done
      end;
      epsilon, order, None

    | ElbowDs ->
      let n_sample = !Parameters.n_sample in
      let eff_n_sample = min n_sample n in
      let idx = Array.init n Fun.id in
      for i = 0 to eff_n_sample - 1 do
        let j = i + Random.int (n - i) in
        let tmp = idx.(i) in idx.(i) <- idx.(j); idx.(j) <- tmp
      done;
      let sample = Array.sub idx 0 eff_n_sample in
      if !Parameters.verbose then
        Printf.eprintf
          "(%s): Computing dist_star for %d sampled points (out of %d) in %d dimensions...\n%!"
          Sys.argv.(0) eff_n_sample n d;
      let max_densities = compute_dist_star sample in
      if !Parameters.verbose then
        Printf.eprintf "%s\r(%s): dist_star: done %d/%d points.\n%!"
          String.TermIO.clear Sys.argv.(0) eff_n_sample eff_n_sample;
      Array.sort (fun (_, _, _, r1) (_, _, _, r2) -> compare r1 r2) max_densities;
      let title =
        Printf.sprintf
          "Part 3: elbow in dist_star, sampled (n_sample=%d, n=%d, D=%d, log_vol(unit D-ball)=%.6g)"
          eff_n_sample n d log_vol1 in
      let elbow_rank = print_dist_star_table title max_densities in
      let epsilon = let _, _, _, r = max_densities.(elbow_rank) in r in
      let order = Array.init n Fun.id in
      if !Parameters.order = "random" then begin
        for i = n - 1 downto 1 do
          let j = Random.int (i + 1) in
          let tmp = order.(i) in order.(i) <- order.(j); order.(j) <- tmp
        done
      end;
      epsilon, order, None

    | Adaptive ->
      if !Parameters.verbose then
        Printf.eprintf
          "(%s): Computing dist_star for all %d points in %d dimensions (O(n^2))...\n%!"
          Sys.argv.(0) n d;
      let all_indices = Array.init n Fun.id in
      let max_densities = compute_dist_star all_indices in
      if !Parameters.verbose then
        Printf.eprintf "%s\r(%s): dist_star: done %d/%d points.\n%!"
          String.TermIO.clear Sys.argv.(0) n n;
      (* Build dist_star lookup by original index before sorting *)
      let ds_by_orig = Array.make n infinity in
      Array.iteri
        (fun orig_i (_, _, _, ds) -> ds_by_orig.(orig_i) <- ds)
        max_densities;
      (* Sort by increasing dist_star for the table and for processing order *)
      Array.sort (fun (_, _, _, r1) (_, _, _, r2) -> compare r1 r2) max_densities;
      let title =
        Printf.sprintf
          "Adaptive: dist_star for all points (n=%d, D=%d, log_vol(unit D-ball)=%.6g)"
          n d log_vol1 in
      ignore (print_dist_star_table title max_densities);
      (* Processing order: increasing dist_star *)
      let order = Array.init n Fun.id in
      Array.sort (fun i j -> compare ds_by_orig.(i) ds_by_orig.(j)) order;
      (* epsilon is irrelevant in adaptive mode; pass per-point ds_by_orig via dist_star_all *)
      0., order, Some ds_by_orig in

  (* ----------------------------------------------------------------
     Greedy leader clustering
     ---------------------------------------------------------------- *)
  let epsilon_of =
    match dist_star_all with
    | Some ds -> (fun i -> ds.(i))
    | None    -> (fun _ -> epsilon) in

  if !Parameters.verbose then
    Printf.eprintf "(%s): Running greedy leader clustering...\n%!" Sys.argv.(0);

  let rep_orig, n_reps = greedy_leader order epsilon_of in
  let n_absorbed = n - n_reps in

  (* ----------------------------------------------------------------
     Output clustering section
     ---------------------------------------------------------------- *)
  let mode_str =
    match mode with
    | Explicit f -> Printf.sprintf "explicit epsilon=%.15g" f
    | Elbow1nn   -> Printf.sprintf "1-NN elbow epsilon=%.15g" epsilon
    | ElbowDs    -> Printf.sprintf "dist_star elbow epsilon=%.15g" epsilon
    | Adaptive   -> "adaptive (per-point dist_star)" in
  section_sep ();
  Printf.printf
    "=== Clustering: greedy leader (%s, index=%s, D=%d) ===\n\
     # n=%d n_representatives=%d n_absorbed=%d compression=%.1f%%\n\
     name%s\trepresentative%s\tstatus\n"
    mode_str (Interfaiss.Type.to_string index_type) d
    n n_reps n_absorbed
    (100. *. float_of_int n_absorbed /. float_of_int n)
    (kmer_hdr "kmer") (kmer_hdr "rep_kmer");

  for i = 0 to n - 1 do
    let rep = rep_orig.(i) in
    Printf.printf "%s%s\t%s%s\t%s\n"
      names.(i)   (kmer_col names.(i))
      names.(rep) (kmer_col names.(rep))
      (if rep = i then "rep" else "abs")
  done;

  if !Parameters.verbose then
    Printf.eprintf
      "(%s): done. %d representatives, %d absorbed (%.1f%% compression).\n%!"
      Sys.argv.(0) n_reps n_absorbed
      (100. *. float_of_int n_absorbed /. float_of_int n)
