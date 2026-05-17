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

(* Epsilon estimation diagnostics for KPop twisted (standard) coordinates.
   Given a KPopTwisted file, outputs up to three TSV tables to stdout:

   Part 1 - Elbow in sorted 1-NN distances:
     For each point the distance to its nearest neighbour is computed
     (plain Euclidean, no inertia weighting).  The distances are sorted
     and the elbow is identified via the kneedle heuristic.  The elbow
     distance is a data-driven candidate for the epsilon threshold in
     greedy leader clustering.

   Part 2 - Running density estimate per point:
     For each point all n-1 pairwise distances are computed and sorted.
     At every rank k the quantity k / V(d_k, D) is emitted, where
       V(r, D) = pi^(D/2) / Gamma(D/2+1) * r^D
     is the volume of a D-dimensional Euclidean ball of radius r.
     The row with the maximum density for each point is marked with *.

   Part 3 - Elbow in max running density (sampled):
     A random sample of n_sample starting points is chosen.  For each,
     all exact distances to the remaining points are computed and the
     maximum of k/V(d_k,D) over all ranks k is recorded.  The n_sample
     maxima are sorted and the kneedle elbow is found.  This is
     sub-quadratic (O(n_sample * n)) and provides an alternative epsilon
     estimate complementary to Part 1.

   When -k is given, a DNA sequence column is added after each hash name
   by back-translating the hex hash: each hex digit pair encodes two bases
   in base 4 (A=0, C=1, G=2, T=3).

   NOTE: Part 2 is O(n^2) in both time and output volume. *)

open BiOCamLib.Better
open KPop

module Defaults =
  struct
    let output = "all"
    let n_sample = 200
    (* 0 = do not decode hashes *)
    let kmer_len = 0
    let verbose = false
  end

module Parameters =
  struct
    let input = ref ""
    let output = ref Defaults.output
    let n_sample = ref Defaults.n_sample
    let kmer_len = ref Defaults.kmer_len
    let verbose = ref Defaults.verbose
  end

let usage () =
  Printf.eprintf
    "Usage: %s <twisted_prefix> [OPTIONS]\n\
     Options:\n\
     \  -o, --output <parts>    which parts to output: 1, 2, 3, or all (default: %s)\n\
     \  -n, --n-sample <n>      random sample size for Part 3 (default: %d)\n\
     \  -k, --kmer <len>        k-mer length for hash back-translation (default: no decoding)\n\
     \  -v, --verbose           verbose output on stderr\n\
     \  -h, --help              print this message and exit\n%!"
    Sys.argv.(0)
    Defaults.output
    Defaults.n_sample

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
    | "-o" | "--output" ->
      let v = get_next arg in
      (match v with
      | "1" | "2" | "3" | "all" -> Parameters.output := v
      | _ ->
        Printf.eprintf "Error: --output must be 1, 2, 3, or all (got '%s')\n%!" v;
        usage (); exit 1)
    | "-n" | "--n-sample" ->
      let v = get_next arg in
      (match int_of_string_opt v with
      | Some n when n > 0 -> Parameters.n_sample := n
      | _ ->
        Printf.eprintf "Error: --n-sample must be a positive integer (got '%s')\n%!" v;
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
  end

let twisted = Twisted.of_binary ~verbose:!Parameters.verbose !Parameters.input

let () = Random.self_init ()

let () =
  let do_part1 = !Parameters.output = "1" || !Parameters.output = "all" in
  let do_part2 = !Parameters.output = "2" || !Parameters.output = "all" in
  let do_part3 = !Parameters.output = "3" || !Parameters.output = "all" in
  let n_sample = !Parameters.n_sample in
  let kmer_len = !Parameters.kmer_len in

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

  (* Plain Euclidean distance between two standard-coordinate rows *)
  let euclidean a b =
    let s = ref 0. in
    for i = 0 to d - 1 do
      let x = Float.Array.unsafe_get a i -. Float.Array.unsafe_get b i in
      s := !s +. x *. x
    done;
    sqrt !s in

  (* log(Gamma(d/2 + 1)):
       d even -> integer argument k = d/2+1: log((d/2)!) = sum_{i=1}^{d/2} log(i)
       d odd  -> half-integer m+1/2, m=(d+1)/2:
                 log(Gamma(m+1/2)) = log(sqrt(pi)) + sum_{i=1}^{m} log(2i-1) - m*log(2) *)
  let log_gamma_d2_plus1 =
    if d mod 2 = 0 then begin
      let s = ref 0. in
      for i = 1 to d / 2 do
        s := !s +. log (float_of_int i)
      done;
      !s
    end else begin
      let m = (d + 1) / 2 in
      let s = ref (0.5 *. log Float.pi) in
      for i = 1 to m do
        s := !s +. log (float_of_int (2 * i - 1))
      done;
      s := !s -. float_of_int m *. log 2.;
      !s
    end in

  (* log(V(1, d)) where V(r, d) = pi^(d/2) / Gamma(d/2+1) * r^d
     k / V(d_k, d) = k * exp(-log_vol1 - d*log(d_k))               *)
  let log_vol1 =
    0.5 *. float_of_int d *. log Float.pi -. log_gamma_d2_plus1 in
  let fd = float_of_int d in

  (* Back-translate a hex k-mer hash to its DNA sequence and reverse complement,
     returned as "fwd/revcomp" in a single string.
     Encoding: A=0, C=1, G=2, T=3 in base 4, stored as a hex integer.
     Uses Int64 so hashes up to k=32 are handled correctly.

     Two formats are supported:
       Non-gapped  "3f7"  - single hex integer; length = kmer_len (from -k).
                            Returns "?/?" if -k was not given or hex is invalid.
       Gapped      "8_d"  - underscore-separated hex segments; bases are
                            distributed equally across segments (kmer_len / n_segs
                            each), so -k is required here too.  The RC reverses
                            segment order and takes the RC of each segment
                            individually.  Returns "?/?" if -k was not given or
                            any segment is not valid hex. *)
  let bases = [| 'A'; 'C'; 'G'; 'T' |] in
  let comp = [| 'T'; 'G'; 'C'; 'A' |] in
  (* Decode a single hex integer as n_bases base-4 digits; return (fwd, rc) *)
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

  (* Append a decoded DNA column after a hash name, or nothing if -k was not given *)
  let kmer_col name =
    if kmer_len = 0 then ""
    else "\t" ^ decode_hash name in

  (* Corresponding header fragment *)
  let kmer_hdr suffix =
    if kmer_len = 0 then ""
    else "\t" ^ suffix in

  (* Number of parts already printed, for inter-section spacing *)
  let parts_done = ref 0 in
  let section_sep () =
    if !parts_done > 0 then print_char '\n';
    incr parts_done in

  (* ----------------------------------------------------------------
     Part 1: sorted 1-NN distances and elbow
     ---------------------------------------------------------------- *)
  if do_part1 then begin
    if !Parameters.verbose then
      Printf.eprintf "(%s): Computing 1-NN distances for %d points in %d dimensions...\n%!"
        Sys.argv.(0) n d;

    let nn1 = Array.init n (fun i ->
      let best = ref infinity in
      for j = 0 to n - 1 do
        if j <> i then begin
          let dist = euclidean data.(i) data.(j) in
          if dist < !best then best := dist
        end
      done;
      names.(i), !best) in

    let sorted_nn1 = Array.copy nn1 in
    Array.sort (fun (_, d1) (_, d2) -> compare d1 d2) sorted_nn1;

    let elbow_rank = Clustering.kneedle_elbow (fun i -> snd sorted_nn1.(i)) n in

    section_sep ();
    Printf.printf
      "=== Part 1: sorted 1-NN distances (n=%d, D=%d) ===\n\
       # elbow at rank %d, distance %.15g\n\
       rank\tname%s\t1nn_distance\telbow\n"
      n d elbow_rank (snd sorted_nn1.(elbow_rank))
      (kmer_hdr "kmer");
    Array.iteri
      (fun rank (name, dist) ->
        Printf.printf "%d\t%s%s\t%.15g\t%s\n"
          rank name (kmer_col name) dist
          (if rank = elbow_rank then "*" else ""))
      sorted_nn1
  end;

  (* ----------------------------------------------------------------
     Part 2: running density estimate k / V(d_k, D) per point
     ---------------------------------------------------------------- *)
  if do_part2 then begin
    if !Parameters.verbose then
      Printf.eprintf
        "(%s): Computing density curves (%d points x %d neighbours = %d pairs)...\n%!"
        Sys.argv.(0) n (n - 1) (n * (n - 1));

    section_sep ();
    Printf.printf
      "=== Part 2: running density estimate k/V(d_k,D) (D=%d, log_vol(unit D-ball)=%.6g) ===\n\
       # density(k, d_k) = k * exp(-log_vol1 - D*log(d_k))\n\
       name%s\tk\tdistance\tneighbour%s\tdensity\tmax\n"
      d log_vol1
      (kmer_hdr "kmer") (kmer_hdr "neighbour_kmer");

    for i = 0 to n - 1 do
      let dists = Array.init (n - 1) (fun j ->
        let j' = if j < i then j else j + 1 in
        euclidean data.(i) data.(j'), j') in
      Array.sort (fun (d1, _) (d2, _) -> compare d1 d2) dists;
      let densities = Array.init (n - 1) (fun idx ->
        let dist, _ = dists.(idx) in
        if dist <= 0. then infinity
        else float_of_int (idx + 1) *. exp (-. log_vol1 -. fd *. log dist)) in
      let max_idx = ref 0 in
      Array.iteri
        (fun idx density ->
          if density > densities.(!max_idx) then max_idx := idx)
        densities;
      Array.iteri
        (fun idx density ->
          let dist, j' = dists.(idx) in
          Printf.printf "%s%s\t%d\t%.15g\t%s%s\t%.15g\t%s\n"
            names.(i) (kmer_col names.(i))
            (idx + 1) dist
            names.(j') (kmer_col names.(j'))
            density
            (if idx = !max_idx then "*" else ""))
        densities;
      if !Parameters.verbose && (i + 1) mod 100 = 0 then
        Printf.eprintf "%s\r(%s): density curves: done %d/%d points%!"
          String.TermIO.clear Sys.argv.(0) (i + 1) n
    done;
    if !Parameters.verbose then
      Printf.eprintf "%s\r(%s): density curves: done %d/%d points.\n%!"
        String.TermIO.clear Sys.argv.(0) n n
  end;

  (* ----------------------------------------------------------------
     Part 3: elbow in max running density, sampled
     ---------------------------------------------------------------- *)
  if do_part3 then begin
    (* Draw a random sample of min(n_sample, n) distinct starting points
       via a partial Fisher-Yates shuffle *)
    let eff_n_sample = min n_sample n in
    let idx = Array.init n Fun.id in
    for i = 0 to eff_n_sample - 1 do
      let j = i + Random.int (n - i) in
      let tmp = idx.(i) in idx.(i) <- idx.(j); idx.(j) <- tmp
    done;
    let sample = Array.sub idx 0 eff_n_sample in

    if !Parameters.verbose then
      Printf.eprintf
        "(%s): Computing max densities for %d sampled points (out of %d) in %d dimensions...\n%!"
        Sys.argv.(0) eff_n_sample n d;

    (* For each sampled point compute all n-1 distances, then find the
       maximum of k/V(d_k,D) over all ranks k, recording the rank k* where
       it is achieved.  k* allows inverting the density to a distance via
       epsilon = (k* / (max_density * V(1,D))) ^ (1/D). *)
    let max_densities = Array.init eff_n_sample (fun si ->
      let i = sample.(si) in
      let dists = Array.init (n - 1) (fun j ->
        let j' = if j < i then j else j + 1 in
        euclidean data.(i) data.(j')) in
      Array.sort compare dists;
      let best = ref 0. and best_k = ref 0 and best_dist = ref 0. in
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
      names.(i), !best, !best_k, !best_dist) in

    (* Sort by increasing dist_star; elbow estimated from these distances *)
    Array.sort (fun (_, _, _, r1) (_, _, _, r2) -> compare r1 r2) max_densities;

    let elbow_rank =
      Clustering.kneedle_elbow (fun i -> let _, _, _, r = max_densities.(i) in r) eff_n_sample in

    section_sep ();
    Printf.printf
      "=== Part 3: elbow in max running density, sampled \
       (n_sample=%d, n=%d, D=%d, log_vol(unit D-ball)=%.6g) ===\n\
       # elbow at rank %d, dist_star %.15g\n\
       rank\tname%s\tmax_density\tk_star\tdist_star\telbow\n"
      eff_n_sample n d log_vol1
      elbow_rank (let _, _, _, r = max_densities.(elbow_rank) in r)
      (kmer_hdr "kmer");
    Array.iteri
      (fun rank (name, max_d, k_star, dist_star) ->
        Printf.printf "%d\t%s%s\t%.15g\t%d\t%.15g\t%s\n"
          rank name (kmer_col name) max_d k_star dist_star
          (if rank = elbow_rank then "*" else ""))
      max_densities;
    if !Parameters.verbose then
      Printf.eprintf "(%s): done.\n%!" Sys.argv.(0)
  end
