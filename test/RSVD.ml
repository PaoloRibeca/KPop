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

(* Benchmark / accuracy comparison of CA.rsvd vs CA.twist.
   Loads a KPop count database, builds the chi-matrix as in KPopTwist (using the
   same defaults: no k-mer filtering, no subsampling, no threshold), then computes:
     - the exact SVD via CA.twist
     - a randomised truncated SVD via CA.rsvd
   and reports, per CA dimension:
     - the exact and approximate singular values
     - the per-dimension relative error in sv
     - the sign-aligned L2 relative error in column principal coordinates
   followed by a one-line summary of the maximum errors.

   NOTE: when --dimensions 0 (the default), rsvd uses all k_ca = min(m,n)-1
   dimensions.  With power iterations enabled, this requires
   dimensions + oversampling <= n_samples; for small databases use -q 0 or
   specify an explicit -d value. *)

open BiOCamLib
open KPop

module Defaults =
  struct
    (* 0 = all available CA dimensions *)
    let dimensions = 0
    let n_oversampling = 10
    let n_power_iter = 2
    let threads = Processes.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let input = ref ""
    let dimensions = ref Defaults.dimensions
    let n_oversampling = ref Defaults.n_oversampling
    let n_power_iter = ref Defaults.n_power_iter
    let threads = ref Defaults.threads
    let output = ref ""
    let verbose = ref Defaults.verbose
  end

let usage () =
  Printf.eprintf
    "Usage: %s <binary_db_prefix> [OPTIONS]\n\
     Options:\n\
     \  -d, --dimensions <n>     CA dimensions to request from rsvd (default: all)\n\
     \  -p, --oversampling <n>   sketch oversampling (default: %d)\n\
     \  -q, --power-iter <n>     subspace power iterations (default: %d)\n\
     \  -T, --threads <n>        computing threads (default: %d)\n\
     \  -o, --output <prefix>    write per-k-mer row inertia to <prefix>.KPopTwisted\n\
     \  -v, --verbose            verbose output\n\
     \  -h, --help               print this message and exit\n%!"
    Sys.argv.(0)
    Defaults.n_oversampling
    Defaults.n_power_iter
    Defaults.threads

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
    | "-d" | "--dimensions" ->
      Parameters.dimensions := int_of_string (get_next arg)
    | "-p" | "--oversampling" ->
      Parameters.n_oversampling := int_of_string (get_next arg)
    | "-q" | "--power-iter" ->
      Parameters.n_power_iter := int_of_string (get_next arg)
    | "-T" | "--threads" ->
      Parameters.threads := int_of_string (get_next arg)
    | "-o" | "--output" ->
      Parameters.output := get_next arg
    | "-v" | "--verbose" ->
      Parameters.verbose := true
    | "-h" | "--help" ->
      usage (); exit 0
    | s when s.[0] = '-' ->
      Printf.eprintf "Error: unknown option %s\n%!" s;
      usage (); exit 1
    | s ->
      if !Parameters.input = "" then
        Parameters.input := s
      else begin
        Printf.eprintf "Error: unexpected argument %s\n%!" s;
        usage (); exit 1
      end);
    incr i
  done;
  if !Parameters.input = "" then begin
    Printf.eprintf "Error: no input database specified\n%!";
    usage (); exit 1
  end

(* Load the k-mer database *)
let db = KMerDB.of_binary ~verbose:!Parameters.verbose !Parameters.input

(* Exact SVD via CA.twist *)
let twister_exact, twisted_exact, _ =
  CA.twist
    ~threads:!Parameters.threads
    ~verbose:!Parameters.verbose
    db

(* Determine how many dimensions to request *)
let k_ca =
  let tw = twister_exact.Twister.twister.Matrix.matrix in
  Array.length tw.Matrix.Base.row_names

let k =
  if !Parameters.dimensions = 0 then k_ca
  else begin
    if !Parameters.dimensions < 1 || !Parameters.dimensions > k_ca then begin
      Printf.eprintf
        "Error: --dimensions %d out of range [1, %d]\n%!"
        !Parameters.dimensions k_ca;
      exit 1
    end;
    !Parameters.dimensions
  end

(* Randomised SVD via CA.rsvd *)
let twister_rsvd, twisted_rsvd, _ =
  CA.rsvd
    ~n_oversampling:!Parameters.n_oversampling
    ~n_power_iter:!Parameters.n_power_iter
    ~threads:!Parameters.threads
    ~verbose:!Parameters.verbose
    db k

(* Helpers: singular value for dimension d from inertia array (sv = sqrt(inertia)) *)
let sv_of d inertia_arr =
  sqrt (Float.Array.get inertia_arr d)

(* Helpers: column principal coordinate G[j,d] from a Twisted.t *)
let g_jd j d (ts: Twisted.t) =
  Float.Array.get ts.Twisted.twisted.Matrix.matrix.Matrix.Base.data.(j) d

let () =
  let inertia_e = twister_exact.Twister.inertia.Matrix.matrix.Matrix.Base.data.(0) in
  let inertia_r = twister_rsvd.Twister.inertia.Matrix.matrix.Matrix.Base.data.(0) in

  let n_samples = Array.length twisted_exact.Twisted.twisted.Matrix.matrix.Matrix.Base.row_names in

  (* Per-dimension statistics *)
  (* Relative error in singular value *)
  let sv_rel = Array.make k 0. in
  (* Sign-aligned L2 relative error in column coords *)
  let coord_l2 = Array.make k 0. in

  for d = 0 to k - 1 do
    (* Singular value relative error *)
    let sv_e = sv_of d inertia_e in
    let sv_r = sv_of d inertia_r in
    sv_rel.(d) <- abs_float (sv_r -. sv_e) /. sv_e;

    (* Sign alignment: compare G_exact[:,d] dot G_rsvd[:,d] *)
    let dot = ref 0. in
    for j = 0 to n_samples - 1 do
      dot := !dot +. g_jd j d twisted_exact *. g_jd j d twisted_rsvd
    done;
    let sign = if !dot >= 0. then 1. else -1. in

    (* L2 norms *)
    let norm_sq_e = ref 0. in
    let diff_sq = ref 0. in
    for j = 0 to n_samples - 1 do
      let g_e = g_jd j d twisted_exact in
      let g_r = sign *. g_jd j d twisted_rsvd in
      norm_sq_e := !norm_sq_e +. g_e *. g_e;
      diff_sq := !diff_sq +. (g_r -. g_e) *. (g_r -. g_e)
    done;
    coord_l2.(d) <- sqrt (!diff_sq /. !norm_sq_e)
  done;

  (* Print header *)
  Printf.printf
    "%-6s  %14s  %14s  %12s  %12s\n"
    "dim" "sv_exact" "sv_rsvd" "sv_rel_err" "coord_L2_err";
  Printf.printf "%s\n%!" (String.make 64 '-');

  let max_sv_rel = ref 0. in
  let max_coord_l2 = ref 0. in
  for d = 0 to k - 1 do
    let sv_e = sv_of d inertia_e in
    let sv_r = sv_of d inertia_r in
    Printf.printf
      "%-6d  %14.6f  %14.6f  %12.2e  %12.2e\n"
      (d + 1) sv_e sv_r sv_rel.(d) coord_l2.(d);
    if sv_rel.(d) > !max_sv_rel then max_sv_rel := sv_rel.(d);
    if coord_l2.(d) > !max_coord_l2 then max_coord_l2 := coord_l2.(d)
  done;

  Printf.printf "%s\n" (String.make 64 '-');
  Printf.printf
    "Max sv relative error:             %12.2e (%.4f%%)\n\
     Max coord L2 relative error:       %12.2e (%.4f%%)\n%!"
    !max_sv_rel (!max_sv_rel *. 100.)
    !max_coord_l2 (!max_coord_l2 *. 100.);

  (* Per-k-mer row contribution to total inertia.
     Each k-mer i contributes CTR_i = r[i] * sum_d T[d,i]^2 * sv_d^2 to the
     total inertia, where r[i] is the row mass and T[d,i] is the twister entry.
     Row masses are not stored in the Twister output and must be recomputed from
     the raw counts using the same formula as CA.prepare_chi:
       r[i] = (1/n_samples) * sum_j X[i,j] / col_sums[j]
     where col_sums[j] covers exactly the k-mers that survived filtering. *)
  if !Parameters.output <> "" then begin
    let tw = twister_exact.Twister.twister.Matrix.matrix in
    let kmer_names = tw.Matrix.Base.col_names in
    let n_kmers = Array.length kmer_names in
    let n_samp = db.core.n_cols in
    let n_f = float_of_int n_samp in

    (* Map kmer names back to original db row indices *)
    let old_idx = Array.map
      (fun name -> Better.StringHashtbl.find db.row_names_to_idx name)
      kmer_names in

    (* Column sums over the filtered k-mer set *)
    let col_sums = Array.make n_samp 0. in
    Array.iter
      (fun ki ->
        for j = 0 to n_samp - 1 do
          col_sums.(j) <- col_sums.(j) +.
            KMerDB.CountBAVector.(db.core.data.(j).@(ki))
        done)
      old_idx;

    (* Row masses: r[i] = (1/n) * sum_j X[i,j] / col_sums[j] *)
    let row_masses = Array.init n_kmers (fun i ->
      let ki = old_idx.(i) in
      let s = ref 0. in
      for j = 0 to n_samp - 1 do
        if col_sums.(j) > 0. then
          s := !s +. KMerDB.CountBAVector.(db.core.data.(j).@(ki)) /. col_sums.(j)
      done;
      !s /. n_f) in

    (* Row contributions: CTR_i = r[i] * sum_d T[d,i]^2 * inertia[d].
       Use all k_ca dimensions from the exact SVD for the most accurate values. *)
    let k_ca_all = Array.length tw.Matrix.Base.row_names in
    let inertia_all = twister_exact.Twister.inertia.Matrix.matrix.Matrix.Base.data.(0) in
    let ctr = Array.init n_kmers (fun i ->
      let s = ref 0. in
      for d = 0 to k_ca_all - 1 do
        let t_di = Float.Array.get tw.Matrix.Base.data.(d) i in
        s := !s +. t_di *. t_di *. Float.Array.get inertia_all d
      done;
      row_masses.(i) *. !s) in

    (* Total inertia = sum_i CTR_i = sum_d sv_d^2 (sanity check) *)
    let total_inertia = ref 0. in
    for d = 0 to k_ca_all - 1 do
      total_inertia := !total_inertia +. Float.Array.get inertia_all d
    done;

    (* Write <output>.KPopTwisted:
       - inertia matrix: 1 row "inertia", 1 col "row_inertia", value = total inertia
       - twisted matrix: n_kmers rows, 1 col "row_inertia", values = CTR_i *)
    let contrib: Twisted.t = {
      inertia = {
        Matrix.which = Matrix.Type.Inertia;
        matrix = {
          Matrix.Base.col_names = [| "row_inertia" |];
          row_names = [| "inertia" |];
          data = [| Float.Array.make 1 !total_inertia |]
        }
      };
      twisted = {
        Matrix.which = Matrix.Type.Twisted;
        matrix = {
          Matrix.Base.col_names = [| "row_inertia" |];
          row_names = kmer_names;
          data = Array.init n_kmers (fun i -> Float.Array.make 1 ctr.(i))
        }
      }
    } in
    Twisted.to_binary ~verbose:!Parameters.verbose contrib !Parameters.output;
    Printf.printf
      "Row inertia written to '%s.KPopTwisted' (%d k-mers, total inertia = %.6f)\n%!"
      !Parameters.output n_kmers !total_inertia
  end
