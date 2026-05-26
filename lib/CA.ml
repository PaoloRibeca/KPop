(*
    CA.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    CA.ml implements the correspondence-analysis pipelines that produce
    a Twister and Twisted database from a k-mer count database.  It
    provides the `twist` entry point (full LAPACK `dgesdd_` SVD via the
    `CA.c` stub) and the `rsvd` entry point (randomised truncated SVD
    via the `RSVD.c` stub; Halko, Martinsson & Tropp 2011), together
    with the shared `prepare_chi` step that handles k-mer keep-listing,
    random subsampling, count thresholding, condition-number filtering,
    and the chi-square matrix construction.

    This program was designed and developed by the author(s),
    with the assistance of the following AI tool(s):
      2026 Claude (Anthropic).
    The final logic and implementation were reviewed and verified in
    their entirety by the author(s).

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

(* We cannot open BiOCamLib here due to the ambiguity between BiOCamLib.Matrix and KPop.Matrix *)
open BiOCamLib.Better

include (
  struct
    (* K-mer filter selector, used for both --kmers-threshold (low-rowsum
       singleton removal) and --kmers-condition-number (low-CTR
       uniform-across-samples removal).  [Off] disables the filter
       (formerly indicated by passing 0.).  [Manual f] uses the legacy
       fractional cutoff (rowsum >= f * max_rowsum, or CTR >= max_CTR / f).
       [Auto] picks the cutoff at the Kneedle elbow of the sorted
       distribution, with no user-supplied magic number. *)
    module Filter =
      struct
        type t =
          | Off
          | Manual of float
          | Auto
        let of_string = function
          | "off" | "0" | "0." | "0.0" -> Off
          | "auto" -> Auto
          | s ->
            (try Manual (float_of_string s) with _ ->
              Exception.raise_unrecognized_initializer __FUNCTION__ "k-mer filter" s)
        let to_string = function
          | Off -> "off"
          | Manual f -> string_of_float f
          | Auto -> "auto"
      end
    (* C binding to LAPACK dgesdd_ for thin SVD via CA.c.
       The input matrix (a) is overwritten during computation.
       On successful return (info = 0):
        u[i * k + d] = left singular vector: k-mer i, dimension d (m x k row-major)
        s[d] = d-th singular value, non-increasing (k elements)
        vt[d * n_cols + j] = right singular vector: dimension d, sample j (k x n row-major)
       where k = min(m, n_cols). *)
    external _dgesdd_:
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int -> int ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int ->
      int
      = "OpenblasSVD_bytecode" "OpenblasSVD"
    (* C binding to the randomised truncated SVD via RSVD.c.
       The input matrix (a) is read-only.
       On successful return (info = 0):
        u[i * k + d] = left singular vector: k-mer i, dimension d (m x k row-major)
        s[d] = d-th singular value, non-increasing (k elements)
        vt[d * n_cols + j] = right singular vector: dimension d, sample j (k x n row-major) *)
    external _rsvd_:
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int -> int ->
      int -> int -> int ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int ->
      int
      = "OpenblasRSVD_bytecode" "OpenblasRSVD"
    let make_ba n =
      Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout n
    (* Read k-mer names to keep from a file (one per line, no header) *)
    let read_kmers_keep fname =
      let res = Hashtbl.create 4096 in
      let input = open_in fname in
      begin try
        while true do
          let line = input_line input |> String.trim in
          if line <> "" then
            Hashtbl.replace res line ()
        done
      with End_of_file -> ()
      end;
      close_in input;
      res
    (* Steps 1-9: filter k-mers and build the chi-matrix from a k-mer database.
       Returns (n_samples, m, kmer_indices, row_masses, chi) where:
       n_samples = number of samples (columns)
       m = number of k-mers retained after filtering
       kmer_indices = original indices of retained k-mers in db
       row_masses = r[i] for each retained k-mer
       chi = m x n_samples chi-matrix, row-major Bigarray *)
    let prepare_chi
        ?(kmers_keep = "")
        ?(kmers_sample = 1.)
        ?(threshold_kmers = Filter.Off)
        ?(condition_number = Filter.Off)
        ?(verbose = false)
        db =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let ( .@@!()<- ) = Bigarray.Array1.unsafe_set in
      let n_kmers_total = db.KMerDB.core.n_rows and n_samples = db.core.n_cols in
      let n_samples_f = float_of_int n_samples in
      (* Check on the number of samples *)
      if n_samples < 2 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf
             "Correspondence analysis requires at least 2 samples (found %d)" n_samples);
      (* Step 1: Apply optional keep-list filter *)
      let kmer_indices =
        if kmers_keep = "" then
          Array.init n_kmers_total (fun i -> i)
        else begin
          let keep_set = read_kmers_keep kmers_keep in
          let res = ref [] in
          for i = 0 to n_kmers_total - 1 do
            if Hashtbl.mem keep_set db.core.idx_to_row_names.(i) then
              List.accum res i
          done;
          let result = Array.of_rlist !res in
          if verbose then
            Printf.eprintf "%s K-mer keep-list: retained %d/%d k-mers.\n%!"
              prefix (Array.length result) n_kmers_total;
          result
        end in
      (* Step 2: Random subsampling without replacement (partial Fisher-Yates) *)
      let kmer_indices =
        let n = Array.length kmer_indices in
        let k = int_of_float (float_of_int n *. kmers_sample) in
        let k = max 0 (min k n) in
        if k >= n then
          kmer_indices
        else begin
          let a = Array.copy kmer_indices in
          for i = 0 to k - 1 do
            let j = i + Random.int (n - i) in
            let tmp = a.(i) in
            a.(i) <- a.(j);
            a.(j) <- tmp
          done;
          let result = Array.sub a 0 k in
          Array.sort compare result;
          if verbose then
            Printf.eprintf "%s K-mer sampling: retained %d/%d k-mers.\n%!"
              prefix k n;
          result
        end in
      (* Step 3: Compute row sums from original counts (used for threshold) *)
      let m0 = Array.length kmer_indices in
      let row_sums = Array.make m0 0. in
      Array.iteri
        (fun new_i old_i ->
          let s = ref 0. in
          for j = 0 to n_samples - 1 do
            s := !s +. KMerDB.CountBAVector.(db.core.data.(j).@(old_i))
          done;
          row_sums.(new_i) <- !s)
        kmer_indices;
      (* Step 4: Threshold - keep only k-mers whose row_sum is above a cutoff.
         For Manual f the cutoff is max_rowsum * f (legacy behaviour).
         For Auto the cutoff is the value at the Kneedle elbow of the
         sorted-ascending row_sum distribution: k-mers in the noise tail
         (singletons, rare k-mers) get dropped automatically without a
         user-supplied magic number.  For Off no filtering is performed. *)
      let kmer_indices =
        match threshold_kmers with
        | Filter.Off -> kmer_indices
        | Filter.Manual t ->
          let max_sum = Array.fold_left max 0. row_sums in
          let cutoff = max_sum *. t in
          let res = ref [] in
          Array.iteri
            (fun i old_i ->
              if row_sums.(i) >= cutoff then
                List.accum res old_i)
            kmer_indices;
          let result = Array.of_rlist !res in
          if verbose then
            Printf.eprintf "%s K-mer threshold (manual %g, cutoff %g): retained %d/%d k-mers.\n%!"
              prefix t cutoff (Array.length result) m0;
          result
        | Filter.Auto ->
          let sorted = Array.copy row_sums in
          Array.sort compare sorted;
          let elbow = Clustering.kneedle_elbow (fun i -> sorted.(i)) m0 in
          let cutoff = sorted.(elbow) in
          let res = ref [] in
          Array.iteri
            (fun i old_i ->
              if row_sums.(i) >= cutoff then
                List.accum res old_i)
            kmer_indices;
          let result = Array.of_rlist !res in
          if verbose then
            Printf.eprintf "%s K-mer threshold (auto, Kneedle elbow at row_sum %g): retained %d/%d k-mers.\n%!"
              prefix cutoff (Array.length result) m0;
          result in
      (* Step 5: Remove k-mers with zero total count (safety guard against
         division by zero in the row-mass computation) *)
      let kmer_indices =
        let res = ref [] and n_zero = ref 0 in
        Array.iter
          (fun old_i ->
            let s = ref 0. in
            for j = 0 to n_samples - 1 do
              s := !s +. KMerDB.CountBAVector.(db.core.data.(j).@(old_i))
            done;
            if !s > 0. then
              List.accum res old_i
            else
              incr n_zero)
          kmer_indices;
        if !n_zero > 0 && verbose then
          Printf.eprintf "%s Removed %d zero-count %s.\n%!"
            prefix !n_zero (String.pluralize_int "k-mer" !n_zero);
        Array.of_rlist !res in
      let m = Array.length kmer_indices in
      if m < 2 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf
             "Correspondence analysis requires at least 2 k-mers after \
              filtering (found %d)" m);
      (* Step 6: Compute column sums of the selected k-mers (for normalisation).
          col_sums[j] = sum_i X[kmer_indices[i], j] *)
      let col_sums = Array.make n_samples 0. in
      Array.iter
        (fun old_i ->
          for j = 0 to n_samples - 1 do
            col_sums.(j) <- col_sums.(j) +.
              KMerDB.CountBAVector.(db.core.data.(j).@(old_i))
          done)
        kmer_indices;
      Array.iteri
        (fun j cs ->
          if cs = 0. then
            Exception.raise __FUNCTION__ IO_Format
              (Printf.sprintf "Sample '%s' has zero total count for the selected k-mers"
                db.core.idx_to_col_names.(j)))
        col_sums;
      (* Step 7: Compute row masses r[i] from column-normalised counts.
          X_norm[i,j] = X[i,j] / col_sums[j] (each column sums to 1)
                    N = n_samples (grand total of X_norm)
                 r[i] = sum_j(X_norm[i,j]) / N (row mass) *)
      let row_masses = Array.make m 0. in
      Array.iteri
        (fun new_i old_i ->
          let s = ref 0. in
          for j = 0 to n_samples - 1 do
            s := !s +. KMerDB.CountBAVector.(db.core.data.(j).@(old_i)) /. col_sums.(j)
          done;
          row_masses.(new_i) <- !s /. n_samples_f)
        kmer_indices;
      (* Step 8: Condition-number filter (optional).
         Discard k-mers whose row contribution to total inertia,
          CTR_i = ||S[i,:]||^2
                = (sum_j x_norm[i,j]^2 - n * r[i]^2) / (n * r[i]),
         falls below max(CTR) / condition_number.
         Uses the current col_sums and row_masses to compute CTR_i;
         col_sums and row_masses are then recomputed for the surviving set
         so that the chi-matrix (Step 9) is correctly normalised. *)
      let kmer_indices, m, col_sums, row_masses =
        match condition_number with
        | Filter.Off ->
          kmer_indices, m, col_sums, row_masses
        | _ ->
          let row_norms = Array.init m (fun new_i ->
            let old_i = kmer_indices.(new_i) in
            let ri = row_masses.(new_i) in
            let sq = ref 0. in
            for j = 0 to n_samples - 1 do
              let xn = KMerDB.CountBAVector.(db.core.data.(j).@(old_i)) /. col_sums.(j) in
              sq := !sq +. xn *. xn
            done;
            (!sq -. n_samples_f *. ri *. ri) /. (n_samples_f *. ri)) in
          let max_norm = Array.fold_left max 0. row_norms in
          if max_norm = 0. then
            kmer_indices, m, col_sums, row_masses
          else begin
            (* For Manual the cutoff is max_norm / parameter (legacy semantics);
               for Auto it is the value at the Kneedle elbow of the
               sorted-ascending CTR distribution: nearly-uniform low-CTR k-mers
               in the noise tail get dropped automatically. *)
            let cutoff, log_tag =
              match condition_number with
              | Filter.Manual t -> max_norm /. t, Printf.sprintf "manual %g, cutoff %g" t (max_norm /. t)
              | Filter.Auto ->
                let sorted = Array.copy row_norms in
                Array.sort compare sorted;
                let elbow = Clustering.kneedle_elbow (fun i -> sorted.(i)) m in
                sorted.(elbow), Printf.sprintf "auto, Kneedle elbow at CTR %g" sorted.(elbow)
              | Filter.Off -> assert false (* handled above *) in
            let res = ref [] in
            Array.iteri
              (fun new_i old_i ->
                if row_norms.(new_i) >= cutoff then
                  List.accum res old_i)
              kmer_indices;
            let kmer_indices' = Array.of_rlist !res in
            let m' = Array.length kmer_indices' in
            if m' < 2 then
              Exception.raise __FUNCTION__ IO_Format
                (Printf.sprintf
                   "Condition-number filter too aggressive: only %d k-mer(s) survive \
                    (at least 2 required for CA)" m');
            if verbose then
              Printf.eprintf "%s K-mer condition-number filter (%s): retained %d/%d k-mers.\n%!"
                prefix log_tag m' m;
            (* Recompute column sums for the surviving set *)
            let col_sums' = Array.make n_samples 0. in
            Array.iter
              (fun old_i ->
                for j = 0 to n_samples - 1 do
                  col_sums'.(j) <- col_sums'.(j) +.
                    KMerDB.CountBAVector.(db.core.data.(j).@(old_i))
                done)
              kmer_indices';
            Array.iteri
              (fun j cs ->
                if cs = 0. then
                  Exception.raise __FUNCTION__ IO_Format
                    (Printf.sprintf
                       "Sample '%s' has zero total count after condition-number filter"
                       db.core.idx_to_col_names.(j)))
              col_sums';
            (* Recompute row masses for the surviving set *)
            let row_masses' = Array.make m' 0. in
            Array.iteri
              (fun new_i old_i ->
                let s = ref 0. in
                for j = 0 to n_samples - 1 do
                  s := !s +. KMerDB.CountBAVector.(db.core.data.(j).@(old_i)) /. col_sums'.(j)
                done;
                row_masses'.(new_i) <- !s /. n_samples_f)
              kmer_indices';
            kmer_indices', m', col_sums', row_masses'
          end in
      (* Step 9: Build the chi-matrix S (m x n_samples, row-major).
         With uniform column masses c[j] = 1/N = 1/n_samples:
          P[i,j] = X_norm[i,j] / n_samples
          S[i,j] = (P[i,j] - r[i] * c[j]) / sqrt(r[i] * c[j])
                 = (X_norm[i,j]/n - r[i]/n) / sqrt(r[i]/n) *)
      if verbose then
        Printf.eprintf "%s Building %d x %d chi-matrix...\n%!" prefix m n_samples;
      let chi = make_ba (m * n_samples) in
      Array.iteri
        (fun new_i old_i ->
          let ri = row_masses.(new_i) in
          let sqrt_ri_over_n = sqrt (ri /. n_samples_f) in
          for j = 0 to n_samples - 1 do
            let x_ij = KMerDB.CountBAVector.(db.core.data.(j).@(old_i)) in
            let x_norm_ij = x_ij /. col_sums.(j) in
            let p_ij = x_norm_ij /. n_samples_f in
            let s_ij = (p_ij -. ri /. n_samples_f) /. sqrt_ri_over_n in
            chi.@@!(new_i * n_samples + j) <- s_ij
          done)
        kmer_indices;
      n_samples, m, kmer_indices, row_masses, chi
    (* Assemble Twister.t and Twisted.t outputs from SVD results.
            k_ca = number of CA dimensions to output
        k_stride = number of columns in u_flat (= k_ca for rsvd, min(m,n) for twist)
       The two Twisted.t values share the same inertia matrix. *)
    let assemble_outputs
        k_ca k_stride kmer_indices row_masses
        u_flat sv vt_flat
        n_samples m db =
      let ( .@!() ) = Float.Array.( .@!() ) in
      (* Bigarray.Array1 unsafe accessors -- second naming tier ( .@@!() )
         to avoid colliding with the Float.Array ones above. *)
      let ( .@@!() ) = Bigarray.Array1.unsafe_get in
      let ( .@@!()<- ) = Bigarray.Array1.unsafe_set in
      let dim_names = Array.init k_ca (fun d -> Printf.sprintf "Dim%d" (d + 1))
      and kmer_names = Array.map (fun i -> db.KMerDB.core.idx_to_row_names.(i)) kmer_indices
      and sample_names = Array.sub db.core.idx_to_col_names 0 n_samples in
      let sqrt_n = float_of_int n_samples |> sqrt in
      (* Sign canonicalization: for each dimension d, if the element of U[:,d]
         with the largest absolute value is negative, negate both U[:,d] and VT[d,:]
         so that the dominant component of each left singular vector is positive.
         This makes the sign convention deterministic and consistent between
         the exact (dgesdd) and randomised SVD. *)
      for d = 0 to k_ca - 1 do
        let i_max = ref 0 and abs_max = ref 0. in
        for i = 0 to m - 1 do
          let v = Float.abs u_flat.@@!(i * k_stride + d) in
          if v > !abs_max then begin abs_max := v; i_max := i end
        done;
        if u_flat.@@!(!i_max * k_stride + d) < 0. then begin
          for i = 0 to m - 1 do
            let idx = i * k_stride + d in
            u_flat.@@!(idx) <- -. u_flat.@@!(idx)
          done;
          for j = 0 to n_samples - 1 do
            let idx = d * n_samples + j in
            vt_flat.@@!(idx) <- -. vt_flat.@@!(idx)
          done
        end
      done;
      (* Shared inertia matrix *)
      let inertia_matrix = {
        Matrix.which = Inertia;
        matrix = {
          col_names = dim_names;
          row_names = [| "inertia" |];
          data = [|
            Float.Array.init k_ca
              (fun d ->
                let sv_d = sv.@@!(d) in
                sv_d *. sv_d)
          |]
        }
      } in
      (* Pseudo-inverse of the vector of singular values *)
      let sqrt_eps = sqrt Float.epsilon in
      let sv_psinv =
        Float.Array.init k_ca
          (fun d ->
            let sv_d = sv.@@!(d) in
            if sv_d < sqrt_eps then
              0.
            else
              1. /. sv_d) in
      (* Column standard coordinates - one row per sample *)
      let twisted_samples = {
        Twisted.inertia = inertia_matrix;
        twisted = {
          Matrix.which = Twisted;
          matrix = {
            col_names = dim_names;
            row_names = sample_names;
            data =
              Array.init n_samples
                (fun j ->
                  Float.Array.init k_ca
                    (fun d ->
                      let vt_dj = vt_flat.@@!(d * n_samples + j) in
                      vt_dj *. sqrt_n))
          }
        }
      } in
      (* Row standard coordinates - one row per k-mer *)
      let twisted_kmers = {
        Twisted.inertia = inertia_matrix;
        twisted = {
          Matrix.which = Twisted;
          matrix = {
            col_names = dim_names;
            row_names = kmer_names;
            data =
              Array.init m
                (fun new_i ->
                  let sqrt_ri = sqrt row_masses.(new_i) in
                  Float.Array.init k_ca
                    (fun d ->
                      let u_id = u_flat.@@!(new_i * k_stride + d) in
                      u_id /. sqrt_ri))
          }
        }
      } in
      (* Twister: scaled standard row coordinates (Twister[d,i] = U[i,d] / (sqrt(r[i]) * sv[d])),
         stored with dimensions as rows and k-mers as columns *)
      let twister = {
        Twister.twister = {
          Matrix.which = Twister;
          matrix = {
            col_names = kmer_names;
            row_names = dim_names;
            data =
              Array.init k_ca
                (fun d ->
                  let sv_psinv_d = sv_psinv.@!(d) in
                  Float.Array.init m
                    (fun new_i ->
                      let sqrt_ri = sqrt row_masses.(new_i) in
                      let u_id = u_flat.@@!(new_i * k_stride + d) in
                      u_id /. sqrt_ri *. sv_psinv_d))
          }
        };
        inertia = inertia_matrix
      } in
      twister, twisted_samples, twisted_kmers
    (* Perform full Correspondence Analysis on a k-mer database using LAPACK dgesdd_.
       All min(m, n_samples) - 1 non-trivial CA dimensions are computed. *)
    let twist
        ?(kmers_keep = "")
        ?(kmers_sample = 1.)
        ?(threshold_kmers = Filter.Off)
        ?(condition_number = Filter.Off)
        ?(threads = 1)
        ?(verbose = false)
        db =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n_samples, m, kmer_indices, row_masses, chi =
        prepare_chi ~kmers_keep ~kmers_sample ~threshold_kmers ~condition_number ~verbose db in
      let k = min m n_samples in
      let k_ca = k - 1 in
      if k_ca < 1 then
        Exception.raise __FUNCTION__ IO_Format
          "Correspondence analysis yields no dimensions \
           (at least 2 k-mers and 2 samples are required)";
      let u_flat = make_ba (m * k) in
      let sv = make_ba k in
      let vt_flat = make_ba (k * n_samples) in
      if verbose then
        Printf.eprintf "%s Computing full SVD (%d x %d)...\n%!" prefix m n_samples;
      let info = _dgesdd_ chi m n_samples u_flat sv vt_flat threads in
      if info < 0 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf
             "SVD failed: dgesdd_ returned info = %d \
              (argument %d has an illegal value)" info (-info))
      else if info > 0 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "SVD failed to converge: dgesdd_ returned info = %d" info);
      if verbose then
        Printf.eprintf "%s SVD complete.\n%!" prefix;
      assemble_outputs k_ca k kmer_indices row_masses u_flat sv vt_flat n_samples m db
    (* Perform randomised truncated Correspondence Analysis via subspace iteration
       (Halko, Martinsson & Tropp 2011).  Only the top [dimensions] CA dimensions
       are computed, which is faster and uses less memory than [twist] when
       [dimensions] << min(m, n_samples).
       [n_oversampling]: extra sketch dimensions for accuracy (default 10).
       [n_power_iter]: subspace refinement iterations (default 2). *)
    let rsvd
        ?(kmers_keep = "")
        ?(kmers_sample = 1.)
        ?(threshold_kmers = Filter.Off)
        ?(condition_number = Filter.Off)
        ?(n_oversampling = 10)
        ?(n_power_iter = 2)
        ?(threads = 1)
        ?(verbose = false)
        db
        dimensions =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n_samples, m, kmer_indices, row_masses, chi =
        prepare_chi ~kmers_keep ~kmers_sample ~threshold_kmers ~condition_number ~verbose db in
      let k_ca_max = min m n_samples - 1 in
      if k_ca_max < 1 then
        Exception.raise __FUNCTION__ IO_Format
          "Correspondence analysis yields no dimensions \
           (at least 2 k-mers and 2 samples are required)";
      if dimensions < 1 || dimensions > k_ca_max then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf
             "dimensions must be in [1, %d] (got %d)" k_ca_max dimensions);
      let l = dimensions + n_oversampling in
      if n_power_iter > 0 && l > n_samples then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf
             "dimensions + oversampling (%d) exceeds n_samples (%d) with power \
              iterations enabled; reduce dimensions, reduce oversampling, set \
              n_power_iter=0, or use twist instead"
             l n_samples);
      let u_flat = make_ba (m * dimensions) in
      let sv = make_ba dimensions in
      let vt_flat = make_ba (dimensions * n_samples) in
      if verbose then
        Printf.eprintf
          "%s Computing randomised SVD (target rank %d, sketch %d, power iter %d)...\n%!"
          prefix dimensions l n_power_iter;
      let info =
        _rsvd_ chi m n_samples dimensions n_oversampling n_power_iter
          u_flat sv vt_flat threads in
      if info < 0 then
        Exception.raise __FUNCTION__ IO_Format
          "Randomised SVD failed: memory allocation error in C stub"
      else if info > 0 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Randomised SVD failed: LAPACK returned info = %d" info);
      if verbose then
        Printf.eprintf "%s Randomised SVD complete.\n%!" prefix;
      (* k_stride = dimensions: u_flat is m×dimensions row-major (no padding) *)
      assemble_outputs dimensions dimensions kmer_indices row_masses
        u_flat sv vt_flat n_samples m db
  end: sig
    (* K-mer filter selector (see implementation for full semantics).
       [Off] disables the filter; [Manual f] uses the legacy fractional
       cutoff; [Auto] picks the cutoff at the Kneedle elbow of the
       sorted-ascending distribution. *)
    module Filter:
      sig
        type t =
          | Off
          | Manual of float
          | Auto
        val of_string: string -> t
        val to_string: t -> string
      end
    (* Direct C binding to LAPACK dgesdd_ (thin SVD).  Exposed for unit
       testing of the C stub from test/CA.ml; production code should call
       [twist] instead, which handles all the bookkeeping. *)
    val _dgesdd_:
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int -> int ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int ->
      int
    (* Direct C binding to the randomised truncated SVD.  Exposed for unit
       testing of the C stub from test/CA.ml; production code should call
       [rsvd] instead, which handles all the bookkeeping. *)
    val _rsvd_:
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int -> int ->
      int -> int -> int ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      int ->
      int
    (* Perform full Correspondence Analysis on a k-mer database using LAPACK dgesdd_.
       Applies optional k-mer filtering (keep-list, random subsampling, row-sum
       threshold) and computes all min(m, n_samples) - 1 non-trivial CA dimensions.
       Returns:
         - a Twister.t (the transformation matrix and inertia)
         - a Twisted.t (column standard coordinates, one entry per sample)
         - a Twisted.t (row standard coordinates, one entry per k-mer).
       The two Twisted.t values share the same inertia matrix.
       Singular vector signs are canonicalised so that the dominant component
       of each left singular vector (U[:,d]) is positive, making results
       consistent between exact and randomised SVD. *)
    val twist:
      ?kmers_keep:string ->
      ?kmers_sample:float ->
      ?threshold_kmers:Filter.t ->
      ?condition_number:Filter.t ->
      ?threads:int ->
      ?verbose:bool ->
      KMerDB.t ->
      Twister.t * Twisted.t * Twisted.t
    (* Perform randomised truncated Correspondence Analysis (HMT 2011).
       Computes only the top [dimensions] CA dimensions, which is faster and
       uses less memory than [twist] when dimensions << min(m, n_samples).
       [n_oversampling] (default 10) and [n_power_iter] (default 2) control
       the accuracy/cost trade-off of the randomised SVD. *)
    val rsvd:
      ?kmers_keep:string ->
      ?kmers_sample:float ->
      ?threshold_kmers:Filter.t ->
      ?condition_number:Filter.t ->
      ?n_oversampling:int ->
      ?n_power_iter:int ->
      ?threads:int ->
      ?verbose:bool ->
      KMerDB.t ->
      int ->
      Twister.t * Twisted.t * Twisted.t
  end
)

