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

(* Four-part test for correspondence_stubs.c / rsvd_stubs.c / Correspondence.ml:
     Part 1 - KPopCorrespondenceSVD C stub:
       Calls Correspondence.dgesdd directly on a 4x3 float64 matrix
       and verifies reconstruction accuracy and orthonormality.
     Part 2 - Correspondence.twist end-to-end:
       Builds a synthetic 4-kmer x 3-sample KMerDB, runs CA, and
       verifies that applying the returned twister to each sample's
       normalised spectrum recovers its column principal coordinate.
     Part 3 - KPopRSVD C stub:
       Calls Correspondence.rsvd_stub directly on an 8x6 float64 matrix,
       verifies orthonormality of U/VT, and checks that singular values
       are close to those from dgesdd.
     Part 4 - Correspondence.rsvd end-to-end:
       Runs rsvd on the same synthetic 4-kmer x 3-sample KMerDB as Part 2
       and verifies that the inertia values agree with those from twist. *)

open KPop

let make_ba n =
  Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout n

let fail fmt = Printf.ksprintf (fun s -> Printf.eprintf "FAIL: %s\n%!" s; exit 1) fmt

let check label max_err tol =
  Printf.printf "  %-50s %.2e  %s\n%!" label max_err
    (if max_err <= tol then "PASS" else (fail "%s: %.2e exceeds %.2e" label max_err tol; ""))

(* --------------------------------------------------------------------------
   Part 1: direct test of the KPopCorrespondenceSVD C stub
   Matrix A (4 rows = k-mers, 3 cols = samples, row-major):
     1  2  3
     4  5  6
     7  8  9    <- first three rows are rank-2 (row3 = 2*row2 - row1)
     2  0  4    <- breaks rank-2, so A has rank 3 and k=3 singular values
   -------------------------------------------------------------------------- *)
let () =
  Printf.printf "=== Part 1: dgesdd stub (4x3 matrix) ===\n%!";
  let m = 4 and n = 3 in
  let k = min m n in
  let a_orig = [| 1.; 2.; 3.;
                  4.; 5.; 6.;
                  7.; 8.; 9.;
                  2.; 0.; 4. |] in
  let a       = make_ba (m * n) in
  let u_flat  = make_ba (m * k) in
  let sv      = make_ba k in
  let vt_flat = make_ba (k * n) in
  Array.iteri (fun i x -> Bigarray.Array1.set a i x) a_orig;
  let info = Correspondence.dgesdd a m n u_flat sv vt_flat 1 in
  if info <> 0 then fail "dgesdd returned info = %d" info;
  Printf.printf "  info = 0:                                          PASS\n%!";

  (* Singular values must be non-negative and non-increasing *)
  for d = 0 to k - 1 do
    if Bigarray.Array1.get sv d < 0.
    then fail "sv[%d] = %.6f < 0" d (Bigarray.Array1.get sv d)
  done;
  for d = 0 to k - 2 do
    if Bigarray.Array1.get sv d < Bigarray.Array1.get sv (d + 1)
    then fail "sv[%d] = %.6f < sv[%d] = %.6f (not non-increasing)"
           d (Bigarray.Array1.get sv d) (d+1) (Bigarray.Array1.get sv (d+1))
  done;
  Printf.printf "  singular values non-negative and non-increasing:  PASS\n%!";

  (* Reconstruction: A_orig ≈ U * diag(sv) * VT *)
  let max_recon = ref 0. in
  for i = 0 to m - 1 do
    for j = 0 to n - 1 do
      let recon = ref 0. in
      for d = 0 to k - 1 do
        recon := !recon
          +. Bigarray.Array1.get u_flat  (i * k + d)
          *. Bigarray.Array1.get sv      d
          *. Bigarray.Array1.get vt_flat (d * n + j)
      done;
      let err = abs_float (!recon -. a_orig.(i * n + j)) in
      if err > !max_recon then max_recon := err
    done
  done;
  check "reconstruction max |A - U·diag(sv)·VT|" !max_recon 1e-9;

  (* U column-orthonormality: U^T U = I_k *)
  let max_u = ref 0. in
  for d1 = 0 to k - 1 do
    for d2 = 0 to k - 1 do
      let dot = ref 0. in
      for i = 0 to m - 1 do
        dot := !dot
          +. Bigarray.Array1.get u_flat (i * k + d1)
          *. Bigarray.Array1.get u_flat (i * k + d2)
      done;
      let err = abs_float (!dot -. (if d1 = d2 then 1. else 0.)) in
      if err > !max_u then max_u := err
    done
  done;
  check "U column-orthonormality max |U^T U - I|" !max_u 1e-9;

  (* VT row-orthonormality: VT VT^T = I_k *)
  let max_vt = ref 0. in
  for d1 = 0 to k - 1 do
    for d2 = 0 to k - 1 do
      let dot = ref 0. in
      for j = 0 to n - 1 do
        dot := !dot
          +. Bigarray.Array1.get vt_flat (d1 * n + j)
          *. Bigarray.Array1.get vt_flat (d2 * n + j)
      done;
      let err = abs_float (!dot -. (if d1 = d2 then 1. else 0.)) in
      if err > !max_vt then max_vt := err
    done
  done;
  check "VT row-orthonormality  max |VT VT^T - I|" !max_vt 1e-9

(* --------------------------------------------------------------------------
   Part 2: end-to-end test of Correspondence.twist on a synthetic KMerDB.
   Count table (4 k-mers x 3 samples):
          s1  s2  s3
     k1:   2   0   4
     k2:   0   3   1
     k3:   1   2   0
     k4:   3   1   2
   Column sums: 6, 6, 7.  Row sums: 6, 4, 3, 6 (all nonzero).
   k = min(4,3) = 3, k_ca = 2 non-trivial CA dimensions.
   Key property verified: applying the twister T[d,i] = U[i,d]/sqrt(r[i])
   to a sample's column-normalised spectrum recovers its column principal
   coordinate G[j,d] = VT[d,j]*sv[d]*sqrt(n_samples). *)
let () =
  Printf.printf "\n=== Part 2: Correspondence.twist (4 kmers x 3 samples) ===\n%!";
  let n_kmers   = 4 and n_samples = 3 in
  (* Count data: counts.(j).(i) = count of k-mer i in sample j *)
  let counts = [| [| 2.; 0.; 1.; 3. |];
                  [| 0.; 3.; 2.; 1. |];
                  [| 4.; 1.; 0.; 2. |] |] in
  (* Build synthetic KMerDB *)
  let kmer_names   = Array.init n_kmers   (fun i -> Printf.sprintf "k%d" (i + 1))
  and sample_names = Array.init n_samples (fun j -> Printf.sprintf "s%d" (j + 1)) in
  let data =
    Array.init n_samples (fun j ->
      let v = KMerDB.CountBAVector.make n_kmers KMerDB.CountBAVector.N.zero in
      for i = 0 to n_kmers - 1 do
        KMerDB.CountBAVector.(v.@(i) <- counts.(j).(i))
      done;
      v) in
  let db = KMerDB.of_core {
    KMerDB.n_cols             = n_samples;
    n_rows                    = n_kmers;
    n_meta                    = 0;
    idx_to_col_names          = sample_names;
    idx_to_row_names          = kmer_names;
    idx_to_meta_names         = [||];
    meta                      = [||];
    data
  } in
  let twister, twisted_samples, _twisted_kmers =
    Correspondence.twist ~threads:1 db in

  (* Shapes: k_ca = min(n_kmers, n_samples) - 1 = 2 *)
  let k_ca = min n_kmers n_samples - 1 in
  let tw = twister.Twister.twister.Matrix.matrix in
  if Array.length tw.Matrix.Base.row_names <> k_ca
  then fail "twister: expected %d rows (dims), got %d" k_ca (Array.length tw.Matrix.Base.row_names);
  if Array.length tw.Matrix.Base.col_names <> n_kmers
  then fail "twister: expected %d cols (k-mers), got %d" n_kmers (Array.length tw.Matrix.Base.col_names);
  let ts = twisted_samples.Twisted.twisted.Matrix.matrix in
  if Array.length ts.Matrix.Base.row_names <> n_samples
  then fail "twisted_samples: expected %d rows, got %d" n_samples (Array.length ts.Matrix.Base.row_names);
  Printf.printf "  output shapes correct:                             PASS\n%!";

  (* Verify: applying twister to each sample's normalised spectrum
     recovers its column principal coordinate.
     x_norm[i,j] = counts[j][i] / col_sum[j]
     G'[j,d]     = sum_i T[d,i] * x_norm[i,j]
     G[j,d]      = twisted_samples.data[j][d]  *)
  let col_sums = Array.init n_samples
    (fun j -> Array.fold_left ( +. ) 0. counts.(j)) in
  let max_apply = ref 0. in
  for j = 0 to n_samples - 1 do
    for d = 0 to k_ca - 1 do
      let g_prime = ref 0. in
      for i = 0 to n_kmers - 1 do
        let t_di     = Float.Array.get tw.Matrix.Base.data.(d) i in
        let x_norm_i = counts.(j).(i) /. col_sums.(j) in
        g_prime := !g_prime +. t_di *. x_norm_i
      done;
      let g = Float.Array.get ts.Matrix.Base.data.(j) d in
      let err = abs_float (!g_prime -. g) in
      if err > !max_apply then max_apply := err
    done
  done;
  check "twister application max |T·x_norm - G|" !max_apply 1e-9

(* --------------------------------------------------------------------------
   Part 3: direct test of the KPopRSVD C stub.
   Matrix A (8 rows, 6 cols, row-major) — digits of pi and e as entries,
   no special structure, full column rank.  We ask for k=2 singular triplets
   with n_oversampling=2 (sketch l=4) and n_power_iter=2, then compare with
   the exact SVD computed by dgesdd.
   -------------------------------------------------------------------------- *)
let () =
  Printf.printf "\n=== Part 3: rsvd_stub (8x6 matrix, k=2, p=2, q=2) ===\n%!";
  let m = 8 and n = 6 in
  let a_orig = [|
    3.; 1.; 4.; 1.; 5.; 9.;
    2.; 6.; 5.; 3.; 5.; 8.;
    9.; 7.; 9.; 3.; 2.; 3.;
    8.; 4.; 6.; 2.; 6.; 4.;
    3.; 3.; 8.; 3.; 2.; 7.;
    9.; 5.; 0.; 2.; 8.; 8.;
    4.; 1.; 9.; 7.; 1.; 6.;
    9.; 3.; 9.; 9.; 3.; 7.
  |] in

  (* Exact SVD via dgesdd *)
  let k_full = min m n in   (* = 6 *)
  let a_exact = make_ba (m * n) in
  Array.iteri (fun i x -> Bigarray.Array1.set a_exact i x) a_orig;
  let u_exact  = make_ba (m * k_full) in
  let sv_exact = make_ba k_full in
  let vt_exact = make_ba (k_full * n) in
  let info = Correspondence.dgesdd a_exact m n u_exact sv_exact vt_exact 1 in
  if info <> 0 then fail "dgesdd (exact) returned info = %d" info;

  (* Randomised SVD via rsvd_stub: k=2, oversampling p=2, power iter q=2 *)
  let k = 2 and p = 2 and q = 2 in
  let a_chi = make_ba (m * n) in
  Array.iteri (fun i x -> Bigarray.Array1.set a_chi i x) a_orig;
  let u_r  = make_ba (m * k) in
  let sv_r = make_ba k in
  let vt_r = make_ba (k * n) in
  let info = Correspondence.rsvd_stub a_chi m n k p q u_r sv_r vt_r 1 in
  if info <> 0 then fail "rsvd_stub returned info = %d" info;
  Printf.printf "  info = 0:                                          PASS\n%!";

  (* U column-orthonormality: U_r^T U_r = I_k *)
  let max_u = ref 0. in
  for d1 = 0 to k - 1 do
    for d2 = 0 to k - 1 do
      let dot = ref 0. in
      for i = 0 to m - 1 do
        dot := !dot
          +. Bigarray.Array1.get u_r (i * k + d1)
          *. Bigarray.Array1.get u_r (i * k + d2)
      done;
      let err = abs_float (!dot -. (if d1 = d2 then 1. else 0.)) in
      if err > !max_u then max_u := err
    done
  done;
  check "U column-orthonormality  max |U^T U - I|" !max_u 1e-9;

  (* VT row-orthonormality: VT_r VT_r^T = I_k *)
  let max_vt = ref 0. in
  for d1 = 0 to k - 1 do
    for d2 = 0 to k - 1 do
      let dot = ref 0. in
      for j = 0 to n - 1 do
        dot := !dot
          +. Bigarray.Array1.get vt_r (d1 * n + j)
          *. Bigarray.Array1.get vt_r (d2 * n + j)
      done;
      let err = abs_float (!dot -. (if d1 = d2 then 1. else 0.)) in
      if err > !max_vt then max_vt := err
    done
  done;
  check "VT row-orthonormality    max |VT VT^T - I|" !max_vt 1e-9;

  (* Singular value relative error vs exact: |sv_r[d] - sv_exact[d]| / sv_exact[d] *)
  let max_sv_rel = ref 0. in
  for d = 0 to k - 1 do
    let sv_e = Bigarray.Array1.get sv_exact d in
    let sv_a = Bigarray.Array1.get sv_r d in
    let rel = abs_float (sv_a -. sv_e) /. sv_e in
    if rel > !max_sv_rel then max_sv_rel := rel
  done;
  check "singular value relative error vs exact" !max_sv_rel 0.02

(* --------------------------------------------------------------------------
   Part 4: end-to-end rsvd on the same 4-kmer x 3-sample KMerDB as Part 2.
   We use dimensions=2, n_oversampling=1 (sketch l=3=n_samples), n_power_iter=1.
   We verify output shapes and check that the per-dimension inertia (sv^2) from
   rsvd matches the corresponding values from twist to within 10%.
   -------------------------------------------------------------------------- *)
let () =
  Printf.printf "\n=== Part 4: Correspondence.rsvd (4 kmers x 3 samples) ===\n%!";
  let n_kmers   = 4 and n_samples = 3 in
  let counts = [| [| 2.; 0.; 1.; 3. |];
                  [| 0.; 3.; 2.; 1. |];
                  [| 4.; 1.; 0.; 2. |] |] in
  let kmer_names   = Array.init n_kmers   (fun i -> Printf.sprintf "k%d" (i + 1))
  and sample_names = Array.init n_samples (fun j -> Printf.sprintf "s%d" (j + 1)) in
  let make_db () =
    let data =
      Array.init n_samples (fun j ->
        let v = KMerDB.CountBAVector.make n_kmers KMerDB.CountBAVector.N.zero in
        for i = 0 to n_kmers - 1 do
          KMerDB.CountBAVector.(v.@(i) <- counts.(j).(i))
        done;
        v) in
    KMerDB.of_core {
      KMerDB.n_cols             = n_samples;
      n_rows                    = n_kmers;
      n_meta                    = 0;
      idx_to_col_names          = sample_names;
      idx_to_row_names          = kmer_names;
      idx_to_meta_names         = [||];
      meta                      = [||];
      data
    } in
  let dimensions = 2 in
  let twister_r, _, _ =
    Correspondence.rsvd ~n_oversampling:1 ~n_power_iter:1 ~threads:1
      (make_db ()) dimensions in
  let twister_t, _, _ =
    Correspondence.twist ~threads:1 (make_db ()) in

  (* Shapes *)
  let tw_r = twister_r.Twister.twister.Matrix.matrix in
  if Array.length tw_r.Matrix.Base.row_names <> dimensions
  then fail "rsvd twister: expected %d rows (dims), got %d" dimensions
         (Array.length tw_r.Matrix.Base.row_names);
  if Array.length tw_r.Matrix.Base.col_names <> n_kmers
  then fail "rsvd twister: expected %d cols (k-mers), got %d" n_kmers
         (Array.length tw_r.Matrix.Base.col_names);
  Printf.printf "  output shapes correct:                             PASS\n%!";

  (* Compare per-dimension inertia (sv^2): |inertia_r[d] - inertia_t[d]| / inertia_t[d] *)
  let inertia_r = twister_r.Twister.inertia.Matrix.matrix.Matrix.Base.data.(0) in
  let inertia_t = twister_t.Twister.inertia.Matrix.matrix.Matrix.Base.data.(0) in
  let max_inertia_rel = ref 0. in
  for d = 0 to dimensions - 1 do
    let i_t = Float.Array.get inertia_t d in
    let i_r = Float.Array.get inertia_r d in
    let rel = abs_float (i_r -. i_t) /. i_t in
    if rel > !max_inertia_rel then max_inertia_rel := rel
  done;
  check "inertia relative error vs twist" !max_inertia_rel 0.10;

  Printf.printf "\nAll tests passed.\n%!"
