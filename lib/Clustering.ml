(*
    Clustering.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Clustering.ml implements clustering algorithms on top of CA-embedded
    coordinates.  Two algorithms are exposed:

     * greedy leader clustering for k-mer and sample spaces, with two
       epsilon-estimation methods (`firstNN`: Kneedle elbow in sorted
       1-NN distances; `density`: Kneedle elbow in sorted `dist_star`
       values) and three processing orders (`inertia`, `firstNN`,
       `density`).  FAISS identifies candidate nearest neighbours;
       exact generalised distances are recomputed for the epsilon
       comparison.

     * HDBSCAN* (Campello et al., 2013/2015), with sparse / dense /
       auto MST modes.  Consumed by the phylo-tree dispatcher in
       Twisted.ml (--phylo-method hdbscan) and by the flat per-point
       cluster assignment used by --clusters-method hdbscan.

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

(* We cannot open BiOCamLib here due to the ambiguity between BiOCamLib.Matrix
   and KPop.Matrix *)
module Numbers = BiOCamLib.Numbers
module Processes = BiOCamLib.Processes
module Tools = BiOCamLib.Tools
module Trees = BiOCamLib.Trees
open BiOCamLib.Better

include (
  struct
    (* Clustering algorithm selector *)
    module Algorithm =
      struct
        type t =
          | Greedy
          | Hdbscan
          | Montecarlo
        let of_string = function
          | "greedy" -> Greedy
          | "hdbscan" -> Hdbscan
          | "montecarlo" -> Montecarlo
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "clustering algorithm" s
        let to_string = function
          | Greedy -> "greedy"
          | Hdbscan -> "hdbscan"
          | Montecarlo -> "montecarlo"
      end
    (* Strategy used to estimate the epsilon threshold for greedy leader clustering *)
    module GreedyEpsilon =
      struct
        type t =
          (* Kneedle elbow in sorted FAISS 1-NN distances *)
          | FirstNN
          (* Kneedle elbow in sorted dist_star values *)
          | Density
          (* Per-point dist_star as the absorption threshold; the processing
             order is forced to ascending dist_star regardless of [order_] *)
          | Adaptive
        let of_string = function
          | "firstNN" -> FirstNN
          | "density" -> Density
          | "adaptive" -> Adaptive
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "greedy epsilon strategy" s
        let to_string = function
          | FirstNN -> "firstNN"
          | Density -> "density"
          | Adaptive -> "adaptive"
      end
    (* Order in which points are presented to the greedy leader clusterer *)
    module Order =
      struct
        type t =
          (* Decreasing row inertia proxy: most informative points first *)
          | Inertia
          (* Increasing FAISS 1-NN distance: densest regions first, O(n log n) *)
          | FirstNN
          (* Increasing dist_star: densest regions first, O(n^2) *)
          | Density
        let of_string = function
          | "inertia" -> Inertia
          | "firstNN" -> FirstNN
          | "density" -> Density
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "clustering order" s
        let to_string = function
          | Inertia -> "inertia"
          | FirstNN -> "firstNN"
          | Density -> "density"
      end
    (* Kneedle elbow: find the rank in [0, m - 1] that maximises the
       deviation between the normalised rank i/(m - 1) and the normalised
       value (get_val i - get_val 0) / (get_val (m - 1) - get_val 0).
       [get_val] is assumed non-decreasing on [0, m - 1].  Returns 0 when
       the value range is zero (flat distribution) *)
    let kneedle_elbow get_val m =
      let v0 = get_val 0 and vm = get_val (m - 1) in
      let range = vm -. v0 in
      if range = 0. then 0
      else begin
        let fm1 = float_of_int (m - 1) in
        let best = ref 0 and best_dev = ref neg_infinity in
        for i = 0 to m - 1 do
          let xi = float_of_int i /. fm1 in
          let yi = (get_val i -. v0) /. range in
          let dev = abs_float (xi -. yi) in
          if dev > !best_dev then begin best_dev := dev; best := i end
        done;
        !best
      end
    (* FAISS self-query: build an index over the embedded rows, query it with
       the same rows asking for the k nearest neighbours of each, delete the
       index, and return the (offsets, distances) result.  Used both for the
       1-NN elbow in [run] and for the sparse k-NN graph in [Hdbscan] *)
    let faiss_self_query ~index_type ~k data =
      let n = Array.length data in
      let d = if n = 0 then 0 else Float.Array.length data.(0) in
      let data_ba = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout n d in
      Matrix.Base.to_bigarray data data_ba;
      let index = Interfaiss.create ~index_type d in
      Interfaiss.train index data_ba;
      Interfaiss.add index data_ba;
      let offsets, distances = Interfaiss.query index data_ba k in
      Interfaiss.delete index;
      offsets, distances
    (* Greedy leader clustering via FAISS nearest-neighbour search.
       Epsilon is estimated automatically from a kneedle elbow in either
       sorted 1-NN distances (GreedyEpsilon.FirstNN) or sorted dist_star values
       (GreedyEpsilon.Density).  Estimation tables and the final cluster assignment
       are printed to stdout.
        [coords] : n × d array of standard CA coordinates
        [names] : name for each of the n points
        [inertia_vec] : d-dimensional inertia vector (sv_d^2 for each CA dimension).
       Returns [rep_orig] where [rep_orig.(i)] is the original index of the
       representative that point i was assigned to (rep_orig.(i) = i for
       representatives themselves) *)
    let run_greedy
        ?(verbose = false)
        ~what_label
        ~epsilon_
        ~order_
        ~density_sample_number
        ~index_type
        ~metric
        ~distance
        ~distance_normalize
        coords names inertia_vec =
      let n = Array.length coords in
      let d = Float.Array.length inertia_vec in
      let idx_str = Interfaiss.Type.to_string index_type in
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let ( .@!() ) = Float.Array.( .@!() ) in
      (* Metric weights, embeddings, generalized distance.
         Pre-apply sqrt(metric_weight[d]) scaling so that L2 in embedding
         space equals the weighted generalized distance for all L2-compatible
         distance functions (Euclidean, cosine, angle).  Normalize additionally
         for cosine/angle.  The flat unit metric is then used with
         Space.Distance.compute for exact distance comparisons *)
      let mw = Space.Distance.Metric.compute metric inertia_vec in
      let flat = Float.Array.make d 1. in
      let compute_embedding coord =
        let v = Float.Array.init d (fun k -> coord.@!(k) *. sqrt mw.@!(k)) in
        if distance_normalize then begin
          let norm = Space.Distance.compute_norm distance flat v in
          if norm > 0. then
            Float.Array.init d (fun k -> v.@!(k) /. norm)
          else v
        end else v in
      let embeds = Array.map compute_embedding coords in
      (* Exact generalized distance between two embedding vectors.
         FAISS is used only to identify candidate nearest neighbours;
         the exact distance is always recomputed here for epsilon comparison *)
      let embed_dist a b = Space.Distance.compute distance flat a b in
      (* log(V(1,D)) for dist_star.
         Euclidean ball volume formula, used as a heuristic normalisation factor
         regardless of the chosen distance metric *)
      let log_gamma_d2p1 =
        if d mod 2 = 0 then begin
          let s = ref 0. in
          for i = 1 to d / 2 do s := !s +. log (float_of_int i) done; !s
        end else begin
          let m_half = (d + 1) / 2 in
          let s = ref (0.5 *. log Float.pi) in
          for i = 1 to m_half do
            s := !s +. log (float_of_int (2 * i - 1))
          done;
          s := !s -. float_of_int m_half *. log 2.; !s
        end in
      let log_vol1 = 0.5 *. float_of_int d *. log Float.pi -. log_gamma_d2p1 in
      let fd = float_of_int d in
      (* Single-row float32 Bigarray reused for all FAISS add/query calls *)
      let buf = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout 1 d in
      let load_buf i =
        Float.Array.iteri (fun j x -> buf.{0, j} <- x) embeds.(i) in
      (* Row inertia proxy: sum_d lambda_d * coord[d]^2.
         Based on raw standard CA coords, independent of the distance metric,
         so Order.Inertia always reflects CA variance *)
      let inertia_proxy () =
        Array.init n (fun i ->
          let s = ref 0. in
          for dim = 0 to d - 1 do
            let t = coords.(i).@!(dim) in
            s := !s +. inertia_vec.@!(dim) *. t *. t
          done;
          !s) in
      (* dist_star for each point in orig_indices: the distance at which k/V(d_k,D)
         is maximised, using the generalized embedding distance throughout *)
      let compute_dist_star orig_indices =
        let m = Array.length orig_indices in
        Array.init m (fun si ->
          let i = orig_indices.(si) in
          let dists = Array.init (n - 1) (fun j ->
            let j' = if j < i then j else j + 1 in
            embed_dist embeds.(i) embeds.(j')) in
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
          if verbose && (si + 1) mod 100 = 0 then
            Printf.eprintf "%s\r%s dist_star for %s: done %d/%d%!"
              String.TermIO.clear prefix what_label (si + 1) m;
          names.(i), !best, !best_k, !best_dist) in
      (* Greedy leader: FAISS index grows as representatives are added.
         FAISS identifies the nearest existing representative (by L2 on embeddings);
         the exact generalized distance is used for the epsilon comparison.
         [epsilon_of] maps an original point index to its absorption threshold:
         a constant function for the global-epsilon strategies (FirstNN, Density),
         a per-point lookup for the Adaptive strategy *)
      let greedy_leader order epsilon_of =
        let rep_orig = Array.make n (-1) in
        let rep_of_faiss = Array.make n (-1) in
        let n_reps = ref 0 in
        let index = Interfaiss.create ~index_type d in
        Array.iter (fun i ->
          if !n_reps > 0 then begin
            load_buf i;
            let offsets, _ = Interfaiss.query index buf 1 in
            let faiss_slot = Int64.to_int offsets.{0, 0} in
            let rep_idx = rep_of_faiss.(faiss_slot) in
            let dist = embed_dist embeds.(i) embeds.(rep_idx) in
            if dist < epsilon_of i then
              rep_orig.(i) <- rep_idx
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
      (* Helper: run a FAISS batch 1-NN query over all n embedded points and
         return an array of exact generalized 1-NN distances indexed by
         original index.  Defensively queries [min n 4] columns and scans
         each row for the first non-self entry: approximate indices like
         HNSW are not guaranteed to return self at column 0 (it can land
         later, or be missed altogether), so a naive `col 0 if non-self
         else col 1' lookup is vulnerable -- in the worst case self lands
         at col 1 and we'd silently use a wrong candidate, or at both col
         0 and col 1 (duplicate-vector quirk) and we'd report distance 0
         as the 1-NN distance. *)
      let faiss_nn1_distances () =
        let k_safe = min n 4 in
        let offsets, _ = faiss_self_query ~index_type ~k:k_safe embeds in
        Array.init n (fun i ->
          let nn_i = ref (-1) in
          let k = ref 0 in
          while !nn_i < 0 && !k < k_safe do
            let j = Int64.to_int offsets.{i, !k} in
            if j >= 0 && j <> i then nn_i := j;
            incr k
          done;
          if !nn_i < 0 then 0.
          else embed_dist embeds.(i) embeds.(!nn_i)) in
      (* Stdout section separator *)
      let parts_done = ref 0 in
      let section_sep () =
        if !parts_done > 0 then print_char '\n'; incr parts_done in
      let metric_str = Space.Distance.Metric.to_string metric in
      let distance_str = Space.Distance.to_string distance in
      let order_str = Order.to_string order_ in
      (* Epsilon estimation.
         Returns (epsilon, ds_by_orig option, nn1_by_orig option).
         ds_by_orig: dist_star indexed by original point index (reused for Order.Density).
         nn1_by_orig: 1-NN distance indexed by original index  (reused for Order.FirstNN) *)
      let epsilon, ds_all_opt, nn1_from_eps =
        match epsilon_ with
        | GreedyEpsilon.Adaptive ->
          (* Per-point absorption threshold: compute dist_star for all n points;
             each point's own dist_star becomes its absorption epsilon, and
             the processing order is forced to ascending dist_star regardless
             of [order_].  No global elbow is computed *)
          if verbose then
            Printf.eprintf
              "%s Computing per-point dist_star for all %d %s (O(n^2))...\n%!"
              prefix n what_label;
          let raw = compute_dist_star (Array.init n Fun.id) in
          if verbose then
            Printf.eprintf "%s\r%s dist_star for %s: done %d/%d.\n%!"
              String.TermIO.clear prefix what_label n n;
          let ds_by_orig = Array.map (fun (_, _, _, ds) -> ds) raw in
          let sorted = Array.copy raw in
          Array.sort (fun (_, _, _, a) (_, _, _, b) -> compare a b) sorted;
          section_sep ();
          Printf.printf
            "=== Adaptive epsilon table for %s: per-point dist_star, all points \
             (n=%d, D=%d, metric=%s, distance=%s, \
             log_vol(unit D-ball)=%.6g) ===\n\
             rank\tname\tmax_density\tk_star\tdist_star\n"
            what_label n d metric_str distance_str log_vol1;
          Array.iteri
            (fun rank (name, max_dens, k_star, dist_star) ->
              Printf.printf "%d\t%s\t%.15g\t%d\t%.15g\n"
                rank name max_dens k_star dist_star)
            sorted;
          0., Some ds_by_orig, None
        | GreedyEpsilon.Density when order_ = Order.Density ->
          (* Compute dist_star for ALL n points; reuse for both elbow and ordering *)
          if verbose then
            Printf.eprintf
              "%s Computing dist_star for all %d %s (O(n^2))...\n%!"
              prefix n what_label;
          let raw = compute_dist_star (Array.init n Fun.id) in
          if verbose then
            Printf.eprintf "%s\r%s dist_star for %s: done %d/%d.\n%!"
              String.TermIO.clear prefix what_label n n;
          let ds_by_orig = Array.map (fun (_, _, _, ds) -> ds) raw in
          let sorted = Array.copy raw in
          Array.sort (fun (_, _, _, a) (_, _, _, b) -> compare a b) sorted;
          let elbow =
            kneedle_elbow (fun i -> let (_, _, _, r) = sorted.(i) in r) n in
          let eps = let (_, _, _, r) = sorted.(elbow) in r in
          section_sep ();
          Printf.printf
            "=== Epsilon estimation for %s: elbow in dist_star, all points \
             (n=%d, D=%d, metric=%s, distance=%s, \
             log_vol(unit D-ball)=%.6g) ===\n\
             # elbow at rank %d, dist_star %.15g\n\
             rank\tname\tmax_density\tk_star\tdist_star\telbow\n"
            what_label n d metric_str distance_str log_vol1 elbow eps;
          Array.iteri
            (fun rank (name, max_dens, k_star, dist_star) ->
              Printf.printf "%d\t%s\t%.15g\t%d\t%.15g\t%s\n"
                rank name max_dens k_star dist_star
                (if rank = elbow then "*" else ""))
            sorted;
          eps, Some ds_by_orig, None
        | GreedyEpsilon.Density ->
          (* Compute dist_star for a random sample; covers Order.Inertia and Order.FirstNN *)
          let eff = min density_sample_number n in
          let idx = Array.init n Fun.id in
          for i = 0 to eff - 1 do
            let j = i + Random.int (n - i) in
            let tmp = idx.(i) in idx.(i) <- idx.(j); idx.(j) <- tmp
          done;
          let sample = Array.sub idx 0 eff in
          if verbose then
            Printf.eprintf
              "%s Computing dist_star for %d sampled %s (out of %d)...\n%!"
              prefix eff what_label n;
          let raw = compute_dist_star sample in
          if verbose then
            Printf.eprintf "%s\r%s dist_star for %s: done %d/%d.\n%!"
              String.TermIO.clear prefix what_label eff eff;
          let sorted = Array.copy raw in
          Array.sort (fun (_, _, _, a) (_, _, _, b) -> compare a b) sorted;
          let elbow =
            kneedle_elbow (fun i -> let (_, _, _, r) = sorted.(i) in r) eff in
          let eps = let (_, _, _, r) = sorted.(elbow) in r in
          section_sep ();
          Printf.printf
            "=== Epsilon estimation for %s: elbow in dist_star, sampled \
             (n_sample=%d, n=%d, D=%d, metric=%s, distance=%s, \
             log_vol(unit D-ball)=%.6g) ===\n\
             # elbow at rank %d, dist_star %.15g\n\
             rank\tname\tmax_density\tk_star\tdist_star\telbow\n"
            what_label eff n d metric_str distance_str log_vol1 elbow eps;
          Array.iteri
            (fun rank (name, max_dens, k_star, dist_star) ->
              Printf.printf "%d\t%s\t%.15g\t%d\t%.15g\t%s\n"
                rank name max_dens k_star dist_star
                (if rank = elbow then "*" else ""))
            sorted;
          eps, None, None
        | GreedyEpsilon.FirstNN ->
          (* FAISS batch 1-NN; nn1_by_orig is also returned for Order.FirstNN *)
          if verbose then
            Printf.eprintf
              "%s Computing 1-NN distances for %d %s via FAISS %s...\n%!"
              prefix n what_label idx_str;
          let nn1_by_orig = faiss_nn1_distances () in
          if verbose then Printf.eprintf "%s done.\n%!" prefix;
          let sorted_nn1 = Array.init n (fun i -> names.(i), nn1_by_orig.(i)) in
          Array.sort (fun (_, a) (_, b) -> compare a b) sorted_nn1;
          let elbow = kneedle_elbow (fun i -> snd sorted_nn1.(i)) n in
          let eps = snd sorted_nn1.(elbow) in
          section_sep ();
          Printf.printf
            "=== Epsilon estimation for %s: sorted 1-NN distances \
             (n=%d, D=%d, metric=%s, distance=%s, index=%s) ===\n\
             # elbow at rank %d, distance %.15g\n\
             rank\tname\t1nn_distance\telbow\n"
            what_label n d metric_str distance_str idx_str elbow eps;
          Array.iteri
            (fun rank (name, dist) ->
              Printf.printf "%d\t%s\t%.15g\t%s\n"
                rank name dist (if rank = elbow then "*" else ""))
            sorted_nn1;
          (* For Order.Density: dist_star for all n still needed *)
          let ds_opt =
            match order_ with
            | Order.Density ->
              if verbose then
                Printf.eprintf
                  "%s Computing dist_star for all %d %s for ordering (O(n^2))...\n%!"
                  prefix n what_label;
              let raw = compute_dist_star (Array.init n Fun.id) in
              if verbose then
                Printf.eprintf "%s\r%s dist_star for %s: done %d/%d.\n%!"
                  String.TermIO.clear prefix what_label n n;
              Some (Array.map (fun (_, _, _, ds) -> ds) raw)
            | _ -> None in
          eps, ds_opt, Some nn1_by_orig in
      (* If Order.FirstNN but GreedyEpsilon.Density: FAISS 1-NN wasn't run yet, do it now.
         Adaptive doesn't need it because the order is forced to dist_star *)
      let nn1_opt =
        match epsilon_, order_, nn1_from_eps with
        | GreedyEpsilon.Adaptive, _, _ -> None
        | _, Order.FirstNN, None ->
          if verbose then
            Printf.eprintf
              "%s Computing 1-NN distances for %d %s for ordering \
               via FAISS %s...\n%!"
              prefix n what_label idx_str;
          let nn1 = faiss_nn1_distances () in
          if verbose then Printf.eprintf "%s done.\n%!" prefix;
          Some nn1
        | _, _, opt -> opt in
      (* Build processing order.  Adaptive forces ascending-dist_star order
         regardless of [order_] *)
      let order_arr = Array.init n Fun.id in
      (match epsilon_, order_, ds_all_opt, nn1_opt with
      | GreedyEpsilon.Adaptive, _, Some ds, _ ->
        Array.sort (fun i j -> compare ds.(i) ds.(j)) order_arr
      | _, Order.Inertia, _, _ ->
        let ip = inertia_proxy () in
        Array.sort (fun i j -> compare ip.(j) ip.(i)) order_arr
      | _, Order.FirstNN, _, Some nn1 ->
        Array.sort (fun i j -> compare nn1.(i) nn1.(j)) order_arr
      | _, Order.Density, Some ds, _ ->
        Array.sort (fun i j -> compare ds.(i) ds.(j)) order_arr
      | _ -> assert false);
      (* Threshold lookup: global constant for FirstNN/Density, per-point for Adaptive *)
      let epsilon_of = match epsilon_, ds_all_opt with
        | GreedyEpsilon.Adaptive, Some ds -> (fun i -> ds.(i))
        | _ -> (fun _ -> epsilon) in
      let effective_order_str = match epsilon_ with
        | GreedyEpsilon.Adaptive -> "density (forced by adaptive)"
        | _ -> order_str in
      let epsilon_log_str = match epsilon_ with
        | GreedyEpsilon.Adaptive -> "per-point dist_star"
        | _ -> Printf.sprintf "%.15g" epsilon in
      if verbose then
        Printf.eprintf
          "%s Running greedy leader for %s \
           (epsilon=%s, order=%s, metric=%s, distance=%s)...\n%!"
          prefix what_label epsilon_log_str effective_order_str metric_str distance_str;
      let rep_orig, n_reps = greedy_leader order_arr epsilon_of in
      let n_abs = n - n_reps in
      (* Print cluster assignment table *)
      let epsilon_str = match epsilon_ with
        | GreedyEpsilon.FirstNN -> Printf.sprintf "1-NN elbow epsilon=%.15g" epsilon
        | GreedyEpsilon.Density -> Printf.sprintf "dist_star elbow epsilon=%.15g" epsilon
        | GreedyEpsilon.Adaptive -> "adaptive (per-point dist_star)" in
      section_sep ();
      Printf.printf
        "=== Clustering of %s: greedy leader \
         (%s, metric=%s, distance=%s, index=%s, order=%s, D=%d) ===\n\
         # n=%d n_representatives=%d n_absorbed=%d compression=%.1f%%\n\
         name\trepresentative\tstatus\n"
        what_label epsilon_str metric_str distance_str idx_str effective_order_str d
        n n_reps n_abs
        (100. *. float_of_int n_abs /. float_of_int n);
      for i = 0 to n - 1 do
        let rep = rep_orig.(i) in
        Printf.printf "%s\t%s\t%s\n"
          names.(i) names.(rep) (if rep = i then "rep" else "abs")
      done;
      if verbose then
        Printf.eprintf
          "%s Clustering of %s done. \
           %d representatives, %d absorbed (%.1f%%).\n%!"
          prefix what_label n_reps n_abs
          (100. *. float_of_int n_abs /. float_of_int n);
      rep_orig
    (* MST construction strategy for HDBSCAN.
       Auto tries Sparse first and falls back to Dense via per-edge MST
       completion if the sparse k-NN graph does not cover all MST edges
       (the MST is over mreach distances rather than Euclidean ones, so
       a Euclidean k-NN graph can miss MST edges even with an exact index).
       Sparse fails loudly on disconnection; Dense always succeeds at O(n^2) *)
    module HdbscanMstMode =
      struct
        type t =
          | Auto
          | Sparse
          | Dense
        let of_string = function
          | "auto" -> Auto
          | "sparse" -> Sparse
          | "dense" -> Dense
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "mst mode" s
        let to_string = function
          | Auto -> "auto"
          | Sparse -> "sparse"
          | Dense -> "dense"
      end
    module Hdbscan =
      struct
        (* HDBSCAN* (Campello et al., 2013/2015).  See KPop-PhyloSplits.tex
           section 3.4.gamma for the design rationale.  Pipeline:
             1. Pairwise Euclidean distance matrix on the embedded rows (parallel).
             2. Per-point core distance = distance to min_samples-th NN.
             3. Mutual-reachability edges mreach(i,j) = max(core(i), core(j), d(i,j)).
             4. Kruskal MST under mreach via a small union-find.
             5. Build the binary merge tree from the MST.
             6. Top-down condensation by min_cluster_size, computing per-leaf
                lambda_death within its current condensed cluster.
             7. Per-condensed-cluster persistence
                  = sum_{p in C} (lambda_death(p) - lambda_birth(C)).
             8. Emit each condensed cluster as a Trees.Splits split with
                weight = persistence *)
        (* Small union-find with path compression and rank *)
        module UF =
          struct
            type t = { parent: int array; rank: int array }
            let make n =
              { parent = Array.init n Fun.id; rank = Array.make n 0 }
            let rec find t x =
              let p = t.parent.(x) in
              if p <> x then begin
                let r = find t p in
                t.parent.(x) <- r;
                r
              end else
                x
            let union t x y =
              let px = find t x and py = find t y in
              if px = py then
                None
              else begin
                let new_root =
                  if t.rank.(px) < t.rank.(py) then begin
                    t.parent.(px) <- py;
                    py
                  end else if t.rank.(px) > t.rank.(py) then begin
                    t.parent.(py) <- px;
                    px
                  end else begin
                    t.parent.(py) <- px;
                    t.rank.(px) <- t.rank.(px) + 1;
                    px
                  end in
                Some (px, py, new_root)
              end
          end
        (* Unsafe Float.Array indexing for the hot inner loops below: every
           embedded row has the same length d, guaranteed by the caller *)
        let ( .@!() ) = Float.Array.( .@!() )
        (* Squared Euclidean distance between two embedded rows *)
        let sq_dist a b =
          let n = Float.Array.length a in
          let s = ref 0. in
          for k = 0 to n - 1 do
            let d = a.@!(k) -. b.@!(k) in
            s := !s +. d *. d
          done;
          !s [@@inline]
        (* Euclidean distance between two embedded rows *)
        let eucl_dist a b = sqrt (sq_dist a b) [@@inline]
        (* Mutual-reachability distance:
             mreach(i,j) = max(d(i,j), core(i), core(j)).
           This is the canonical HDBSCAN mreach formula *)
        let mreach_of core d_ij i j =
          max d_ij (max core.(i) core.(j)) [@@inline]
        (* Compute all upper-triangle pairwise distances in parallel.
           Returns a flat array of (dist, i, j) tuples of length n*(n-1)/2 *)
        let pairwise_distances ?(threads = 1) ?(verbose = false) data =
          let n = Array.length data in
          let total_pairs = n * (n - 1) / 2 in
          if total_pairs = 0 then
            [||]
          else begin
            let open String.TermIO in
            let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
            let edges = Array.make total_pairs (0., 0, 0) in
            let next_row = ref 0 and done_pairs = ref 0 in
            let work_threads = max 1 (min threads (n - 1)) in
            Processes.Parallel.process_stream_chunkwise
              (fun () ->
                if !next_row < n - 1 then begin
                  let i = !next_row in
                  incr next_row;
                  i
                end else
                  raise End_of_file)
              (fun i ->
                let row = Array.init (n - i - 1) (fun k ->
                  let j = i + 1 + k in
                  eucl_dist data.(i) data.(j), i, j) in
                i, row)
              (fun (i, row) ->
                (* Base offset for row i in the flat array:
                   sum_{k=0}^{i-1} (n - k - 1) = i*(2*n - i - 1) / 2 *)
                let base = i * (2 * n - i - 1) / 2 in
                Array.blit row 0 edges base (Array.length row);
                done_pairs := !done_pairs + Array.length row;
                if verbose then
                  Printf.eprintf "%s\r%s Pairwise %d/%d edges%!"
                    clear prefix !done_pairs total_pairs)
              work_threads;
            if verbose then
              Printf.eprintf "%s\r%s Pairwise %d/%d edges.\n%!"
                clear prefix !done_pairs total_pairs;
            edges
          end
        (* Per-point core distance: distance to the k-th NN.
           k is clamped to [1, n-1] so it always indexes a real neighbour *)
        let core_distances n k edges =
          let kk = max 1 (min k (n - 1)) in
          (* For each point, gather all its incident edge distances *)
          let per_point = Array.make n [] in
          Array.iter (fun (d, i, j) ->
            per_point.(i) <- d :: per_point.(i);
            per_point.(j) <- d :: per_point.(j))
            edges;
          Array.init n (fun i ->
            let arr = Array.of_list per_point.(i) in
            Array.sort compare arr;
            if Array.length arr >= kk then
              arr.(kk - 1)
            else if Array.length arr > 0 then
              arr.(Array.length arr - 1)
            else
              0.)
        (* Convert (d, i, j) edges into mreach edges using the core distances,
           then sort ascending *)
        let mreach_edges core edges =
          let n = Array.length edges in
          let m = Array.init n (fun k ->
            let d, i, j = edges.(k) in
            mreach_of core d i j, i, j) in
          Array.sort (fun (a, _, _) (b, _, _) -> compare a b) m;
          m
        (* Kruskal MST: walk sorted edges, accept those joining different
           components.  Returns the MST array (size n - 1) together with the
           number of edges actually accepted; the caller checks whether the
           full MST was built and raises with a path-specific message on
           disconnection *)
        let kruskal_mst n sorted_edges =
          let uf = UF.make n in
          let mst = Array.make (max 0 (n - 1)) (0., 0, 0) in
          let count = ref 0 in
          let len = Array.length sorted_edges in
          let k = ref 0 in
          while !count < n - 1 && !k < len do
            let (_, i, j) as edge = sorted_edges.(!k) in
            if UF.union uf i j <> None then begin
              mst.(!count) <- edge;
              incr count
            end;
            incr k
          done;
          mst, !count
        (* Merge-tree node.  Leaves are nodes [0; n-1] (size 1, lambda = infinity,
           left = right = -1).  Internal nodes are [n; 2n - 2] in creation order
           (ascending d, descending lambda) *)
        type merge_node_t = { left: int; right: int; lambda: float; size: int }
        let build_merge_tree n mst =
          let total = max 1 (2 * n - 1) in
          let tree = Array.make total { left = -1; right = -1; lambda = infinity; size = 1 } in
          (* Singletons are leaves; initialised above *)
          let uf = UF.make n in
          let comp_node = Array.init n Fun.id in
          let next_id = ref n in
          Array.iter (fun (d, a, b) ->
            let lambda = if d <= 0. then infinity else 1. /. d in
            let ra = UF.find uf a and rb = UF.find uf b in
            if ra <> rb then begin
              let left = comp_node.(ra) and right = comp_node.(rb) in
              let size = tree.(left).size + tree.(right).size in
              let id = !next_id in
              incr next_id;
              tree.(id) <- { left; right; lambda; size };
              match UF.union uf a b with
              | Some (_, _, new_root) -> comp_node.(new_root) <- id
              | None -> ()
            end)
            mst;
          tree, !next_id - 1
        (* Collect leaves under each merge-tree node, memoised *)
        let make_leaves_cache tree total =
          let cache = Array.make total None in
          let rec leaves node =
            match cache.(node) with
            | Some xs -> xs
            | None ->
              let xs =
                if tree.(node).left = -1 then
                  [| node |]
                else
                  Array.append (leaves tree.(node).left) (leaves tree.(node).right) in
              cache.(node) <- Some xs;
              xs in
          leaves
        (* Top-down condensation.  Emits one entry per condensed cluster as
           (members_sorted_array, persistence).  Returns (clusters, lambda_death)
           so callers that want per-leaf branch lengths (mreach-style trees)
           can reuse the lambda_death information without re-walking the tree. *)
        let condense tree root min_cluster_size leaves_of n =
          let lambda_death = Array.make n 0. in
          let emitted = ref [] in
          (* descend node: walk the merge tree within the CURRENT condensed
             cluster (whose membership is everything below `node` that hasn't
             already been noise-stripped), updating lambda_death for each
             leaf to record the lambda at which it leaves THIS cluster *)
          let rec descend node =
            (* Singleton leaf -- stays in the current cluster *)
            if tree.(node).left = -1 then
              ()
            else begin
              let { left; right; lambda; _ } = tree.(node) in
              let left_small = tree.(left).size < min_cluster_size
              and right_small = tree.(right).size < min_cluster_size in
              match left_small, right_small with
              | true, true ->
                (* Both children too small: the current cluster terminates
                   here.  All leaves under `node` die at lambda *)
                Array.iter
                  (fun p -> lambda_death.(p) <- lambda)
                  (leaves_of node)
              | true, false ->
                (* Left noise-stripped; right is the same cluster, continuing.
                   Set lambda_death for ALL leaves of this subtree at this
                   lambda (a "if nothing deeper happens, you die here"
                   provisional value); the descent into right will overwrite
                   right's leaves with deeper events as appropriate *)
                Array.iter
                  (fun p -> lambda_death.(p) <- lambda)
                  (leaves_of node);
                descend right
              | false, true ->
                Array.iter
                  (fun p -> lambda_death.(p) <- lambda)
                  (leaves_of node);
                descend left
              | false, false ->
                (* Real cluster split: every leaf below `node` leaves the
                   current cluster at lambda (those not already stripped
                   by an earlier noise event) *)
                Array.iter
                  (fun p ->
                    if lambda_death.(p) = 0. then
                      lambda_death.(p) <- lambda)
                  (leaves_of node)
            end in
          (* condense_cluster node birth: enter the new condensed cluster
             rooted at `node` with cluster-birth lambda = birth; emit it (if
             non-trivial), then recurse into any real-split children *)
          let rec condense_cluster node birth =
            let members = leaves_of node in
            Array.iter (fun p -> lambda_death.(p) <- 0.) members;
            descend node;
            let persistence = ref 0. in
            Array.iter (fun p ->
              if lambda_death.(p) > birth then
                persistence := !persistence +. (lambda_death.(p) -. birth))
              members;
            (* Emit only proper, non-trivial clusters: at least min_cluster_size
               members AND strictly fewer than n (the all-elements "root"
               cluster canonicalises to the empty split, which downstream
               consumers reject) *)
            if !persistence > 0.
               && Array.length members >= min_cluster_size
               && Array.length members < n then begin
              let sorted = Array.copy members in
              Array.sort compare sorted;
              emitted := (sorted, !persistence) :: !emitted
            end;
            (* Find the descendant nodes that are real splits and recurse *)
            let rec find_splits node =
              if tree.(node).left <> -1 then begin
                let { left; right; lambda; _ } = tree.(node) in
                let left_small = tree.(left).size < min_cluster_size
                and right_small = tree.(right).size < min_cluster_size in
                match left_small, right_small with
                | true, true -> ()
                | true, false -> find_splits right
                | false, true -> find_splits left
                | false, false ->
                  condense_cluster left lambda;
                  condense_cluster right lambda
              end in
            find_splits node in
          condense_cluster root 0.;
          !emitted, lambda_death
        (* Sparse k-NN graph via FAISS.  Returns (offsets, distances) of
           shape (n, k_query) where k_query = num_neighbors + 1 (the +1
           accounts for the self-match at column 0) *)
        let compute_knn_graph ~index_type ~num_neighbors data =
          let n = Array.length data in
          let k_query = min n (num_neighbors + 1) in
          faiss_self_query ~index_type ~k:k_query data
        (* Sparse mreach edges from the FAISS k-NN graph.  For every (i, j)
           pair appearing in i's neighbour list (or symmetrically j's) we
           emit one canonical edge with i < j (deduplicated), recomputing
           the Euclidean distance in float64 from the source embedding so
           that mreach is bit-stable against the dense path *)
        let sparse_mreach_edges data core offsets =
          let n = Array.length data in
          let k_query = Bigarray.Array2.dim2 offsets in
          let acc = ref [] in
          (* FAISS is not guaranteed to return self at column 0 for approximate
             indices like HNSW; iterating from k = 0 with the [i < j] guard
             filters self (i < i is false) regardless of where it lands *)
          for i = 0 to n - 1 do
            for k = 0 to k_query - 1 do
              let j = Int64.to_int offsets.{i, k} in
              if i < j then begin
                let d_ij = eucl_dist data.(i) data.(j) in
                let mr = mreach_of core d_ij i j in
                acc := (mr, i, j) :: !acc
              end
            done
          done;
          let arr = Array.of_list !acc in
          Array.sort (fun (a, _, _) (b, _, _) -> compare a b) arr;
          arr
        (* Default neighbourhood size for the sparse path: enough above
           min_samples to capture the typical mreach MST edges, capped at
           n - 1 (the maximum meaningful neighbourhood) *)
        let default_num_neighbors ~min_samples n =
          max (min_samples + 1) (min (n - 1) 30)
        (* Build dense mreach edges (all n(n-1)/2 pairs).  Always covers
           every MST edge by construction *)
        let dense_mreach ~prefix ~threads ~verbose ~min_samples n data =
          if verbose then
            Printf.eprintf "%s Computing %d pairwise distances...\n%!"
              prefix (n * (n - 1) / 2);
          let edges = pairwise_distances ~threads ~verbose data in
          if verbose then
            Printf.eprintf "%s Computing core distances (min_samples=%d)...\n%!"
              prefix min_samples;
          let core = core_distances n min_samples edges in
          if verbose then
            Printf.eprintf "%s Building mutual-reachability graph and MST...\n%!" prefix;
          mreach_edges core edges
        (* Build sparse mreach edges from a FAISS k-NN graph.  May miss MST
           edges if k is too small or if the index is approximate; the caller
           detects this via kruskal_mst's edge count.  Returns both the core
           distances and the sorted mreach edges so the caller can reuse the
           cores when completing the MST via inter-component edges *)
        let sparse_mreach ~prefix ~verbose ~min_samples ~num_neighbors ~index_type n data =
          if num_neighbors < min_samples then
            Exception.raise __FUNCTION__ Initialize
              (Printf.sprintf
                 "--phylo-hdbscan-num-neighbors (%d) must be >= --phylo-hdbscan-min-samples (%d) \
                  so the core distance can be read from the k-NN result"
                 num_neighbors min_samples);
          if verbose then
            Printf.eprintf "%s Building FAISS %s k-NN graph (k=%d)...\n%!"
              prefix (Interfaiss.Type.to_string index_type) num_neighbors;
          let offsets, _ = compute_knn_graph ~index_type ~num_neighbors data in
          let k_query = Bigarray.Array2.dim2 offsets in
          (* Walk the FAISS result row, skipping any self-match (which is not
             guaranteed to be at column 0 with approximate indices like HNSW),
             and record the distance to the min_samples-th non-self neighbour.
             Because k_query = num_neighbors + 1 and num_neighbors >= min_samples,
             this loop always finds at least min_samples non-self entries *)
          let core =
            Array.init n (fun i ->
              let count = ref 0 and dist = ref 0. in
              let k = ref 0 in
              while !count < min_samples && !k < k_query do
                let j = Int64.to_int offsets.{i, !k} in
                if j <> i then begin
                  dist := eucl_dist data.(i) data.(j);
                  incr count
                end;
                incr k
              done;
              !dist) in
          if verbose then
            Printf.eprintf "%s Building sparse mutual-reachability graph and MST...\n%!"
              prefix;
          let edges = sparse_mreach_edges data core offsets in
          core, edges
        (* Per-edge MST completion: when the sparse Kruskal leaves the graph
           disconnected, the only candidate edges that can complete the MST
           are those crossing between the surviving components.  We scan
           exactly those inter-component pairs, sort by mreach, and continue
           the Kruskal walk -- avoiding both the full O(n^2) dense pairwise
           recomputation and the work already done in the sparse pass.
           Reuses the cores derived from the sparse path; if the FAISS
           index was approximate, these may differ marginally from the
           exact dense-path cores, so Auto's output may differ from Dense
           by a tiny epsilon when HNSW approximates a few k-NN queries *)
        let complete_mst_per_edge ~prefix ~verbose n sparse_mst sparse_count data core =
          let uf = UF.make n in
          for k = 0 to sparse_count - 1 do
            let _, i, j = sparse_mst.(k) in
            ignore (UF.union uf i j)
          done;
          let component_of = Array.init n (fun i -> UF.find uf i) in
          let stack = Tools.ArrayStack.empty () in
          for i = 0 to n - 1 do
            for j = i + 1 to n - 1 do
              if component_of.(i) <> component_of.(j) then begin
                let d_ij = eucl_dist data.(i) data.(j) in
                let mr = mreach_of core d_ij i j in
                Tools.ArrayStack.push stack (mr, i, j)
              end
            done
          done;
          let inter = Tools.ArrayStack.contents stack in
          if verbose then
            Printf.eprintf "%s Completing MST: scanning %d inter-component edges...\n%!"
              prefix (Array.length inter);
          Array.sort (fun (a, _, _) (b, _, _) -> compare a b) inter;
          let mst = Array.make (max 0 (n - 1)) (0., 0, 0) in
          Array.blit sparse_mst 0 mst 0 sparse_count;
          let count = ref sparse_count in
          let kk = ref 0 in
          let len = Array.length inter in
          while !count < n - 1 && !kk < len do
            let (_, i, j) as edge = inter.(!kk) in
            if UF.union uf i j <> None then begin
              mst.(!count) <- edge;
              incr count
            end;
            incr kk
          done;
          mst, !count
        (* Shared MST -> merge-tree -> condense pipeline.  Both [make_splits]
           and [make_clusters] call this and then differ only in how they
           consume the resulting (members, persistence) tuples *)
        let pipeline ~threads ~verbose ~mst_mode ~num_neighbors ~index_type
            ~min_cluster_size ~min_samples ~prefix n data =
          let k_sparse () =
            match num_neighbors with
            | Some k -> min (n - 1) (max 1 k)
            | None -> default_num_neighbors ~min_samples n in
          let try_sparse_mst () =
            let k = k_sparse () in
            let core, mreach =
              sparse_mreach ~prefix ~verbose ~min_samples
                ~num_neighbors:k ~index_type n data in
            let mst, count = kruskal_mst n mreach in
            mst, count, core in
          let dense_mst () =
            let mreach = dense_mreach ~prefix ~threads ~verbose ~min_samples n data in
            let mst, count = kruskal_mst n mreach in
            if count <> n - 1 then
              Exception.raise __FUNCTION__ Initialize
                (Printf.sprintf "Dense MST is disconnected (got %d edges, expected %d) \
                                 -- this indicates a data anomaly (degenerate or NaN \
                                 coordinates)"
                  count (n - 1));
            mst in
          let mst =
            match mst_mode with
            | HdbscanMstMode.Dense -> dense_mst ()
            | HdbscanMstMode.Sparse ->
              let mst, count, _ = try_sparse_mst () in
              if count <> n - 1 then
                Exception.raise __FUNCTION__ Initialize
                  (Printf.sprintf
                     "MST is disconnected (got %d edges, expected %d): the sparse k-NN \
                      graph with num_neighbors=%d does not cover all \
                      necessary mreach edges; try a larger value, or switch to \
                      mst_mode 'dense' (always works) or 'auto' (sparse \
                      with per-edge completion)"
                     count (n - 1) (k_sparse ()));
              mst
            | HdbscanMstMode.Auto ->
              let mst, count, core = try_sparse_mst () in
              if count = n - 1 then
                mst
              else begin
                if verbose then
                  Printf.eprintf
                    "%s Sparse MST incomplete (%d/%d edges); completing per-edge.\n%!"
                    prefix count (n - 1);
                let mst', count' =
                  complete_mst_per_edge ~prefix ~verbose n mst count data core in
                if count' <> n - 1 then
                  Exception.raise __FUNCTION__ Initialize
                    (Printf.sprintf "Per-edge MST completion failed (got %d edges, \
                                     expected %d) -- this is a bug, please report"
                      count' (n - 1));
                mst'
              end in
          let tree, root = build_merge_tree n mst in
          let leaves_of = make_leaves_cache tree (2 * n - 1) in
          if verbose then
            Printf.eprintf "%s Condensing tree (min_cluster_size=%d)...\n%!"
              prefix min_cluster_size;
          let clusters, lambda_death =
            condense tree root min_cluster_size leaves_of n in
          clusters, tree, root, lambda_death
        (* Branch-length convention selector for [make_tree].

           [Mreach]: build a tree whose topology matches the full binary
              HDBSCAN merge tree, with edge lengths set to the
              mutual-reachability distance interval over which each
              merge node existed.  Branch above an internal node u
              with parent p is 1/lambda(u) - 1/lambda(p) (the mreach
              distance grows from u's birth at 1/lambda(u) up to
              p's birth at 1/lambda(p)); for a leaf the branch is
              1/lambda(p).

           [Persistence]: build a tree whose topology corresponds to
              the condensed-cluster hierarchy (small noise subclusters
              absorbed into their parents as polytomies, matching
              HDBSCAN's standard cluster output), with each cluster's
              edge length set to its persistence score (sum over its
              leaves of lambda_death - lambda_birth).  This is the
              same per-edge weight that the legacy splits emission
              used; the produced tree is what Yggdrasill would
              compatibility-assemble from those splits. *)
        module LengthsMode =
          struct
            type t =
              | Persistence
              | Mreach
            let of_string = function
              | "persistence" -> Persistence
              | "mreach" -> Mreach
              | s ->
                Exception.raise_unrecognized_initializer __FUNCTION__
                  "HDBSCAN branch-length convention" s
            let to_string = function
              | Persistence -> "persistence"
              | Mreach -> "mreach"
          end
        let make_tree
            ?(threads = 1) ?(verbose = false)
            ?(mst_mode = HdbscanMstMode.Auto)
            ?num_neighbors
            ?(index_type = Interfaiss.Type.of_string "hnsw(32)")
            ?(lengths_mode = LengthsMode.Mreach)
            ~min_cluster_size ~min_samples m =
          let row_names = m.Matrix.Base.row_names in
          let n = Array.length row_names in
          if n < 2 then
            (* Degenerate: a single-leaf tree.  Make a 1-leaf Newick. *)
            Trees.Newick.leaf (if n = 1 then row_names.(0) else "")
          else begin
            let open String.TermIO in
            let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
            let clusters, tree, root, lambda_death =
              pipeline ~threads ~verbose ~mst_mode ~num_neighbors ~index_type
                ~min_cluster_size ~min_samples ~prefix n m.data in
            match lengths_mode with
            | LengthsMode.Persistence ->
              (* Emit one (cluster_members, persistence) pair per condensed
                 cluster and let Trees.of_splits assemble them into a
                 hierarchy.  The of_splits polytomy-aware reconstruction
                 handles noise-stripped clusters cleanly. *)
              let splits = Trees.Splits.create row_names in
              List.iter (fun (members, persistence) ->
                let split = Trees.Splits.Split.of_array members in
                Trees.Splits.add_split splits split persistence)
                clusters;
              if verbose then begin
                let num_clusters = List.length clusters in
                Printf.eprintf "%s Persistence mode: %d condensed %s.\n%!"
                  prefix num_clusters (String.pluralize_int "cluster" num_clusters)
              end;
              let _used, nwk, _residual = Trees.of_splits ~verbose splits in
              nwk
            | LengthsMode.Mreach ->
              ignore clusters;
              ignore lambda_death;
              (* Walk the full binary merge tree, emit Newick.  Branch
                 above each child of a merge node is 1/lambda(parent) -
                 1/lambda(child); leaves get 1/lambda(parent). *)
              let rec build_node u =
                if u < n then
                  Trees.Newick.leaf row_names.(u)
                else begin
                  let m_u = tree.(u) in
                  let parent_birth = 1. /. m_u.lambda in
                  let child_branch child =
                    if child < n then
                      max 0. parent_birth
                    else
                      max 0. (parent_birth -. 1. /. tree.(child).lambda) in
                  let left_branch = child_branch m_u.left in
                  let right_branch = child_branch m_u.right in
                  Trees.Newick.join
                    [| Trees.Newick.edge ~length:left_branch (), build_node m_u.left;
                       Trees.Newick.edge ~length:right_branch (), build_node m_u.right |]
                end in
              let nwk = build_node root in
              if verbose then
                Printf.eprintf "%s Mreach mode: built binary merge tree (%d leaves).\n%!"
                  prefix n;
              nwk
          end
        (* Flat per-point cluster assignment.  Returns an int array indexed by
           original point index: each entry is either a cluster id (>= 0,
           contiguous from 0) or -1 for noise (point not in any condensed
           cluster).  Each point is assigned to its smallest containing
           condensed cluster (canonical HDBSCAN flat output) *)
        let make_clusters
            ?(threads = 1) ?(verbose = false)
            ?(mst_mode = HdbscanMstMode.Auto)
            ?num_neighbors
            ?(index_type = Interfaiss.Type.of_string "hnsw(32)")
            ~min_cluster_size ~min_samples data =
          let n = Array.length data in
          let cluster_of = Array.make n (-1) in
          if n < 2 then
            cluster_of
          else begin
            let open String.TermIO in
            let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
            let clusters, _tree, _root, _lambda_death =
              pipeline ~threads ~verbose ~mst_mode ~num_neighbors ~index_type
                ~min_cluster_size ~min_samples ~prefix n data in
            (* condense emits clusters in reverse DFS-preorder (deepest first
               on top of the accumulator list).  Walking via List.rev gives
               root-to-leaf order; overwriting cluster_of as we go means
               deeper / smaller clusters end up winning the assignment *)
            List.iteri (fun cid (members, _) ->
              Array.iter (fun p -> cluster_of.(p) <- cid) members)
              (List.rev clusters);
            (* Relabel cluster ids to a contiguous 0..k-1 range (some of the
               original cid values may have been fully overwritten by deeper
               clusters and so no longer appear) *)
            let cid_remap = Hashtbl.create 16 in
            let next_cid = ref 0 in
            Array.iteri (fun i cid ->
              if cid >= 0 then begin
                let new_cid =
                  match Hashtbl.find_opt cid_remap cid with
                  | Some k -> k
                  | None ->
                    let k = !next_cid in
                    incr next_cid;
                    Hashtbl.add cid_remap cid k;
                    k in
                cluster_of.(i) <- new_cid
              end)
              cluster_of;
            if verbose then begin
              let n_assigned = n - Array.fold_left
                (fun acc c -> if c < 0 then acc + 1 else acc) 0 cluster_of in
              Printf.eprintf "%s Assigned %d/%d points to %d %s (%d noise).\n%!"
                prefix n_assigned n !next_cid
                (String.pluralize_int "cluster" !next_cid)
                (n - n_assigned)
            end;
            cluster_of
          end
      end
    (* HDBSCAN clustering on top of CA-embedded coordinates.  Applies the
       same metric / distance / normalisation pre-scaling as [run] so that
       --metric and --distance are honoured uniformly across the two
       clustering algorithms.  Prints a cluster-assignment table to stdout
       and returns an int array indexed by original point index: each entry
       is a cluster id (>= 0) for assigned points, -1 for noise *)
    let run_hdbscan
        ?(verbose = false)
        ?(threads = 1)
        ~what_label
        ~min_cluster_size
        ?min_samples
        ~mst_mode
        ?num_neighbors
        ~index_type
        ~metric
        ~distance
        ~distance_normalize
        coords names inertia_vec =
      let n = Array.length coords in
      let d = Float.Array.length inertia_vec in
      let idx_str = Interfaiss.Type.to_string index_type in
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let ( .@!() ) = Float.Array.( .@!() ) in
      (* Same pre-scaling as run_greedy: sqrt(metric_weight[d]) so that L2
         in embedding space equals the metric-weighted generalised distance,
         then per-row normalisation for cosine/angle *)
      let mw = Space.Distance.Metric.compute metric inertia_vec in
      let flat = Float.Array.make d 1. in
      let compute_embedding coord =
        let v = Float.Array.init d (fun k -> coord.@!(k) *. sqrt mw.@!(k)) in
        if distance_normalize then begin
          let norm = Space.Distance.compute_norm distance flat v in
          if norm > 0. then
            Float.Array.init d (fun k -> v.@!(k) /. norm)
          else v
        end else v in
      let embeds = Array.map compute_embedding coords in
      let effective_min_samples =
        match min_samples with
        | Some k -> k
        | None -> min_cluster_size in
      if verbose then
        Printf.eprintf
          "%s Running HDBSCAN for %s (min_cluster_size=%d, min_samples=%d, \
           mst_mode=%s, metric=%s, distance=%s, index=%s)...\n%!"
          prefix what_label min_cluster_size effective_min_samples
          (HdbscanMstMode.to_string mst_mode)
          (Space.Distance.Metric.to_string metric)
          (Space.Distance.to_string distance) idx_str;
      let cluster_of =
        Hdbscan.make_clusters ~threads ~verbose ~mst_mode ?num_neighbors
          ~index_type ~min_cluster_size ~min_samples:effective_min_samples embeds in
      (* Print cluster assignment table *)
      let metric_str = Space.Distance.Metric.to_string metric in
      let distance_str = Space.Distance.to_string distance in
      let n_noise =
        Array.fold_left (fun acc c -> if c < 0 then acc + 1 else acc) 0 cluster_of in
      let n_assigned = n - n_noise in
      let n_clusters =
        Array.fold_left max (-1) cluster_of + 1 in
      Printf.printf
        "=== Clustering of %s: HDBSCAN \
         (min_cluster_size=%d, min_samples=%d, mst_mode=%s, metric=%s, distance=%s, \
          index=%s, D=%d) ===\n\
         # n=%d n_clusters=%d n_noise=%d coverage=%.1f%%\n\
         name\tcluster\tstatus\n"
        what_label min_cluster_size effective_min_samples
        (HdbscanMstMode.to_string mst_mode) metric_str distance_str idx_str d
        n n_clusters n_noise
        (100. *. float_of_int n_assigned /. float_of_int n);
      for i = 0 to n - 1 do
        let cid = cluster_of.(i) in
        if cid >= 0 then
          Printf.printf "%s\tC@%d\tmember\n" names.(i) cid
        else
          Printf.printf "%s\tnoise\toutlier\n" names.(i)
      done;
      if verbose then
        Printf.eprintf
          "%s Clustering of %s done. \
           %d points in %d %s, %d noise (%.1f%% coverage).\n%!"
          prefix what_label n_assigned n_clusters
          (String.pluralize_int "cluster" n_clusters) n_noise
          (100. *. float_of_int n_assigned /. float_of_int n);
      cluster_of
    (* The radius the distances themselves imply.
       Where a corpus is a set of groups separated by a divergence cutoff, the distances
       between two randomly drawn points are bimodal -- a low mode of same-group pairs, a
       high mode of different-group ones, and an empty band between -- and any radius in
       that band recovers the groups, which is why the answer is insensitive to where in
       the band it falls.  The band is a property of the distribution and needs no labels
       to find.  Note that the 1-NN distances cannot show it: a point's nearest neighbour
       belongs to its own group, so they lie wholly inside the low mode, which is why
       [GreedyEpsilon.FirstNN] returns a radius an order of magnitude too small on a corpus
       of this shape.  Returns the midpoint of the band, or [None] when the distances are
       unimodal and no radius is implied *)
    let valley_radius ?(pairs = 100000) ?(bins = 200) ?(smooth = 5) ?(floor_ = 0.05) state
        embed_dist embeds =
      let n = Array.length embeds in
      if n < 50 then None
      else begin
        let sampled = Array.make pairs 0. and m = ref 0 in
        for _ = 1 to pairs do
          let a = Random.State.int state n and b = Random.State.int state n in
          if a <> b then begin
            sampled.(!m) <- embed_dist embeds.(a) embeds.(b); incr m
          end
        done;
        let sampled = Array.sub sampled 0 !m in
        Array.sort compare sampled;
        let hi = sampled.(int_of_float (float_of_int !m *. 0.99)) in
        if hi <= 0. then None
        else begin
          let w = hi /. float_of_int bins in
          let hist = Array.make bins 0. in
          Array.iter
            (fun x ->
              let b = int_of_float (x /. w) in
              if b >= 0 && b < bins then hist.(b) <- hist.(b) +. 1.)
            sampled;
          (* A moving average, so that a ragged histogram does not invent extrema *)
          let sm =
            Array.init bins
              (fun b ->
                let s = ref 0. and c = ref 0 in
                for j = b - smooth to b + smooth do
                  if j >= 0 && j < bins then begin s := !s +. hist.(j); incr c end
                done;
                !s /. float_of_int !c) in
          let top = ref 0 in
          Array.iteri (fun b v -> if v > sm.(!top) then top := b) sm;
          (* The tallest OTHER mode, accepted only when the lowest point between it and the
             tallest is under half its height -- otherwise it is a shoulder, not a mode *)
          let second = ref (-1) in
          for b = 0 to bins - 1 do
            if abs (b - !top) > 2 * smooth then begin
              let lo = min b !top and hi = max b !top in
              let v = ref sm.(lo) in
              for j = lo to hi do if sm.(j) < !v then v := sm.(j) done;
              if !v <= 0.5 *. sm.(b) && (!second < 0 || sm.(b) > sm.(!second)) then second := b
            end
          done;
          if !second < 0 then None
          else begin
            let lo_b = min !second !top and hi_b = max !second !top in
            let vb = ref lo_b in
            for j = lo_b to hi_b do if sm.(j) < sm.(!vb) then vb := j done;
            (* The midpoint of the band, not its lowest point.  Where the groups really are
               separated the valley is not a dip but a gap, so its minimum is arbitrary
               within it, and it sits nearer the low mode only because most random pairs are
               between-group and their tail dominates *)
            let lower = Float.min sm.(lo_b) sm.(hi_b) in
            let lo_e = ref !vb and hi_e = ref !vb in
            while !lo_e > lo_b && sm.(!lo_e - 1) <= floor_ *. lower do decr lo_e done;
            while !hi_e < hi_b && sm.(!hi_e + 1) <= floor_ *. lower do incr hi_e done;
            Some ((float_of_int (!lo_e + !hi_e) /. 2. +. 0.5) *. w)
          end
        end
      end
    (* Monte-Carlo search over partitions, scored by the simplified silhouette.
       Greedy leader clustering and HDBSCAN both apply a rule and accept whatever partition
       it yields; neither optimises anything, so neither can recover from a starting point
       that is too fine or too coarse, and on a corpus whose groups differ wildly in size
       both settle far from the answer.  This searches instead: merge and split moves are
       proposed, scored, and accepted or rejected by Metropolis, so the chain can cross a
       ridge that a rule cannot.
       THE SCORE IS THE SIMPLIFIED SILHOUETTE -- each sampled point against class CENTROIDS
       rather than against every other point, which costs O(sample * k * d) instead of
       O(sample * n * d) -- and it is the rare unsupervised criterion that punishes BOTH
       errors, since splitting a real cluster drives the nearest-other-class distance down
       towards the own-class one while merging two drives the own-class distance up.  A
       criterion rewarding only tightness, such as a separation ratio or the depth of the
       valley above, is maximised by shattering the corpus and cannot be optimised directly.
       THE STARTING PARTITION is one leader pass at [valley_radius], deliberately finer than
       the answer since the search merges far more often than it splits; where the distances
       imply no radius, every point starts on its own *)
    (* THE SPACE EVERY CLUSTERER HERE WORKS IN.  A twisted register stores standard
       coordinates, so the metric has to be applied before a distance is taken: scaling each
       axis by the square root of its weight makes a plain distance over the result the
       weighted generalised distance the user asked for *)
    let make_embeddings ~metric ~distance ~distance_normalize coords inertia_vec =
      let ( .@!() ) = Float.Array.( .@!() ) in
      let d = Float.Array.length inertia_vec in
      let mw = Space.Distance.Metric.compute metric inertia_vec in
      let flat = Float.Array.make d 1. in
      let compute_embedding coord =
        let v = Float.Array.init d (fun k -> coord.@!(k) *. sqrt mw.@!(k)) in
        if distance_normalize then begin
          let norm = Space.Distance.compute_norm distance flat v in
          if norm > 0. then Float.Array.init d (fun k -> v.@!(k) /. norm) else v
        end else v in
      Array.map compute_embedding coords,
      fun a b -> Space.Distance.compute distance flat a b
    (* THE DIMENSION THE DATA OCCUPIES, which is not the number of axes it is stored in, by the
       two-nearest-neighbour estimator of Facco, d'Errico, Rodriguez and Laio (2017).  For each
       point the ratio of the distances to its second and first neighbours is Pareto with the
       dimension as its parameter, and -- this being the whole reason only two are used -- the
       local density cancels, so one estimate holds over a set whose density varies from place
       to place.  Read it as the dimension of a NEIGHBOURHOOD: on a set that clusters, both of
       a point's neighbours are usually in its own group, so this measures the spread within a
       group and not the arrangement of the groups, which is a larger number *)
    let intrinsic_dimension ?(sample = 800) state embed_dist embeds =
      let n = Array.length embeds in
      let idx = Array.init (min sample n) (fun _ -> Random.State.int state n) in
      let tot = ref 0. and cnt = ref 0 in
      Array.iter
        (fun i ->
          let r1 = ref infinity and r2 = ref infinity in
          Array.iter
            (fun j ->
              if i <> j then begin
                let dd = embed_dist embeds.(i) embeds.(j) in
                if dd < !r1 then begin r2 := !r1; r1 := dd end
                else if dd < !r2 then r2 := dd
              end)
            idx;
          if !r1 > 0. && !r2 < infinity then begin
            tot := !tot +. log (!r2 /. !r1); incr cnt
          end)
        idx;
      if !cnt > 0 && !tot > 0. then float_of_int !cnt /. !tot else 0.
    (* A VALLEY IN THE DISTANCES BETWEEN PAIRS, as the calibrated detector reports one: [radius] is
       the centre of its bin, [share] the share of pairs below it, [z] its prominence in standard
       deviations of its own resampling distribution, [reappearance] the share of resamples in
       which it is still at least a quarter deep, and [depth] how deep it is, as a share of the
       lower of the two modes it separates *)
    type valley_t = { radius: float; share: float; z: float; reappearance: float; depth: float }
    module ValleyStats = Numbers.OnlineStats (Numbers.Float)
    (* The detector's constants.  The histogram has [valley_bins] bins up to the 99th percentile,
       which is read off [valley_fine_bins] finer ones rather than by sorting; a feature scoring
       under [valley_floor] is simplified away; a valley is kept when its [z] is at least
       [valley_z_min] and its share lies between [valley_share_min] and [valley_share_max]; and
       kept valleys whose shares are closer than [valley_chain] are one valley *)
    let valley_bins = 200
    and valley_fine_bins = 1_000_000
    and valley_floor = 0.25
    and valley_z_min = 3.
    and valley_share_min = 0.002
    and valley_share_max = 0.9
    and valley_chain = 0.015
    (* A moving average over one bin either side, the window shrinking at the ends *)
    let smooth_counts counts =
      let n = Float.Array.length counts in
      Float.Array.init n
        (fun b ->
          let lo = max 0 (b - 1) and hi = min (n - 1) (b + 1) and s = ref 0. in
          for j = lo to hi do s := !s +. Float.Array.get counts j done;
          !s /. float_of_int (hi - lo + 1))
    (* THE EXTREMA OF A SMOOTHED HISTOGRAM, as an alternating peak-trough-...-peak sequence of
       (is a peak, bin, height).  A plateau counts once, at its middle -- the left of the two
       middle bins when it is even -- and a run of like extrema keeps its most extreme member, so
       that a trough always has a peak on each side and neither end of the range is a trough *)
    let valley_extrema sm =
      let n = Float.Array.length sm and ext = ref [] and b = ref 0 in
      while !b < n do
        let e = ref !b and v = Float.Array.get sm !b in
        while !e + 1 < n && Float.Array.get sm (!e + 1) = v do incr e done;
        let left = if !b > 0 then Float.Array.get sm (!b - 1) else neg_infinity
        and right = if !e + 1 < n then Float.Array.get sm (!e + 1) else neg_infinity
        and mid = (!b + !e) / 2 in
        if v > left && v > right then ext := (true, mid, v) :: !ext
        else if v < left && v < right then ext := (false, mid, v) :: !ext;
        b := !e + 1
      done;
      let alt = ref [] in
      List.iter
        (fun ((is_peak, _, v) as x) ->
          match !alt with
          | [] -> if is_peak then alt := [ x ]
          | (last_is_peak, _, last_v) :: rest ->
            if last_is_peak = is_peak then begin
              if (is_peak && v > last_v) || (not is_peak && v < last_v) then alt := x :: rest
            end else alt := x :: !alt)
        (List.rev !ext);
      (match !alt with (false, _, _) :: rest -> rest | l -> l) |> List.rev |> Array.of_list
    (* PERSISTENCE SIMPLIFICATION BY RELATIVE SIZE, smallest score first, until every feature
       scores at least [valley_floor].  A trough T between peaks L and R scores
       (min(L,R) - T) / min(L,R) and goes with the lower of L and R, the right one on a tie.  An
       interior peak P between troughs A and B scores (P - max(A,B)) / min(L',R'), L' and R' being
       the peaks beyond A and B, and goes with the higher of A and B, the right one on a tie -- so
       a bump inside a gap is judged against the modes the gap separates and not against the floor
       beside it, and cannot split one valley into two.  On a tie between scores the first feature
       found goes, troughs before peaks and left to right *)
    let simplify_extrema arr =
      let a = ref arr and go = ref true in
      while !go do
        let arr = !a in
        let len = Array.length arr in
        let v k = let _, _, x = arr.(k) in x in
        let best = ref infinity and lo = ref (-1) and k = ref 1 in
        while !k < len - 1 do
          let p = Float.min (v (!k - 1)) (v (!k + 1)) in
          let s = if p > 0. then (p -. v !k) /. p else 0. in
          if s < !best then begin
            best := s;
            lo := if v (!k - 1) < v (!k + 1) then !k - 1 else !k
          end;
          k := !k + 2
        done;
        k := 2;
        while !k < len - 2 do
          let outer = Float.min (v (!k - 2)) (v (!k + 2)) in
          let s =
            if outer > 0. then (v !k -. Float.max (v (!k - 1)) (v (!k + 1))) /. outer else 0. in
          if s < !best then begin
            best := s;
            lo := if v (!k - 1) > v (!k + 1) then !k - 1 else !k
          end;
          k := !k + 2
        done;
        if !lo < 0 || !best >= valley_floor then go := false
        else a := Array.append (Array.sub arr 0 !lo) (Array.sub arr (!lo + 2) (len - !lo - 2))
      done;
      !a
    let troughs_of_counts counts =
      let arr = smooth_counts counts |> valley_extrema |> simplify_extrema in
      List.init ((Array.length arr - 1) / 2)
        (fun i ->
          let k = (2 * i) + 1 in
          let _, tb, t = arr.(k) and _, lb, l = arr.(k - 1) and _, rb, r = arr.(k + 1) in
          let p = Float.min l r in
          tb, lb, rb, p -. t, if p > 0. then (p -. t) /. p else 0.)
    (* THE PAIRS A HISTOGRAM IS BUILT FROM: every one of them when there are few enough, and
       otherwise a sample drawn with replacement, a point never being paired with itself *)
    type valley_pairs_t =
      | AllPairs of int
      | SampledPairs of int array * int array
    let valley_pairs ~pairs_max ~pairs_sample state n =
      if n * (n - 1) / 2 <= pairs_max then AllPairs n
      else begin
        let pi = Array.make pairs_sample 0 and pj = Array.make pairs_sample 0 in
        for p = 0 to pairs_sample - 1 do
          (* In sequence and not under [and], whose order of evaluation is unspecified and would
             make the draws depend on the compiler *)
          let i = Random.State.int state n in
          let j = ref (Random.State.int state n) in
          while !j = i do j := Random.State.int state n done;
          pi.(p) <- i;
          pj.(p) <- !j
        done;
        SampledPairs (pi, pj)
      end
    let valley_distances embed_dist embeds = function
      | AllPairs n ->
        let d = Float.Array.make (n * (n - 1) / 2) 0. and id = ref 0 in
        for i = 0 to n - 2 do
          for j = i + 1 to n - 1 do
            Float.Array.unsafe_set d !id (embed_dist embeds.(i) embeds.(j));
            incr id
          done
        done;
        d
      | SampledPairs (pi, pj) ->
        Float.Array.init (Array.length pi) (fun p -> embed_dist embeds.(pi.(p)) embeds.(pj.(p)))
    (* A resample of the points: each is drawn [n] times with replacement, and what is kept is how
       many times each was drawn *)
    let valley_multiplicities state resamples n =
      Array.init resamples
        (fun _ ->
          let m = Array.make n 0 in
          for _ = 1 to n do
            let x = Random.State.int state n in
            m.(x) <- m.(x) + 1
          done;
          m)
    (* EVERY CANDIDATE TROUGH OF ONE SET OF DISTANCES, each scored against the resamples.  A
       resample weights every pair by the product of the multiplicities of its two points, is
       binned on the bins the full data set, rescaled to the full number of pairs and smoothed
       alike, and has the prominence and depth of each trough read at that trough's own bins, so
       that what varies from one resample to the next is the data and never where the trough is
       looked for.  Returns each candidate with its depth in every resample *)
    let calibrate_troughs ~multiplicities pairs dists =
      let np = Float.Array.length dists in
      let maximum = Float.Array.fold_left Float.max 0. dists in
      if np = 0 || maximum <= 0. then []
      else begin
        (* The 99th percentile as the upper edge of the first fine bin at which the cumulative count
           passes 99% of the pairs *)
        let fw = maximum *. (1. +. 1e-9) /. float_of_int valley_fine_bins
        and fine = Array.make valley_fine_bins 0 in
        Float.Array.iter
          (fun x ->
            let b = min (valley_fine_bins - 1) (int_of_float (x /. fw)) in
            fine.(b) <- fine.(b) + 1)
          dists;
        let target = 0.99 *. float_of_int np and cum = ref 0 and fb = ref 0 in
        while !fb < valley_fine_bins - 1 && float_of_int (!cum + fine.(!fb)) <= target do
          cum := !cum + fine.(!fb);
          incr fb
        done;
        (* The last bin is the overflow: it counts towards all pairs but takes no part in the
           smoothing or in the search for extrema *)
        let w = float_of_int (!fb + 1) *. fw /. float_of_int valley_bins
        and idx = Bytes.create np and counts = Float.Array.make (valley_bins + 1) 0. in
        Float.Array.iteri
          (fun p x ->
            let b = min valley_bins (int_of_float (x /. w)) in
            Bytes.unsafe_set idx p (Char.unsafe_chr b);
            Float.Array.set counts b (Float.Array.get counts b +. 1.))
          dists;
        let cands = troughs_of_counts (Float.Array.sub counts 0 valley_bins) |> Array.of_list in
        let nc = Array.length cands and nr = Array.length multiplicities
        and total = float_of_int np in
        let proms = Array.make_matrix nc nr 0. and depths = Array.make_matrix nc nr 0. in
        Array.iteri
          (fun r m ->
            let acc = Float.Array.make (valley_bins + 1) 0. in
            let add p weight =
              let b = Char.code (Bytes.unsafe_get idx p) in
              Float.Array.unsafe_set acc b (Float.Array.unsafe_get acc b +. weight) in
            (match pairs with
             | AllPairs n ->
               let id = ref 0 in
               for i = 0 to n - 2 do
                 let mi = m.(i) in
                 if mi = 0 then id := !id + (n - 1 - i)
                 else
                   for j = i + 1 to n - 1 do
                     let mj = Array.unsafe_get m j in
                     if mj > 0 then add !id (float_of_int (mi * mj));
                     incr id
                   done
               done
             | SampledPairs (pi, pj) ->
               Array.iteri
                 (fun p i ->
                   let mm = m.(i) * m.(pj.(p)) in
                   if mm > 0 then add p (float_of_int mm))
                 pi);
            let mass = Float.Array.fold_left ( +. ) 0. acc in
            let scale = if mass > 0. then total /. mass else 0. in
            let smr =
              Float.Array.init valley_bins (fun b -> Float.Array.get acc b *. scale)
              |> smooth_counts in
            Array.iteri
              (fun k (tb, lb, rb, _, _) ->
                let p = Float.min (Float.Array.get smr lb) (Float.Array.get smr rb)
                and t = Float.Array.get smr tb in
                proms.(k).(r) <- p -. t;
                depths.(k).(r) <- if p > 0. then (p -. t) /. p else 0.)
              cands)
          multiplicities;
        let below = Float.Array.make (valley_bins + 1) 0. in
        for b = 1 to valley_bins do
          Float.Array.set below b (Float.Array.get below (b - 1) +. Float.Array.get counts (b - 1))
        done;
        Array.mapi
          (fun k (tb, _, _, prominence, depth) ->
            let stats = ValleyStats.make () in
            Array.iter (ValleyStats.add stats) proms.(k);
            let sd = ValleyStats.sample_standard_deviation stats in
            let z =
              if sd > 0. then prominence /. sd else if prominence > 0. then infinity else 0.
            and reappearing =
              Array.fold_left (fun acc x -> if x >= valley_floor then acc + 1 else acc) 0 depths.(k) in
            { radius = (float_of_int tb +. 0.5) *. w;
              share = (Float.Array.get below tb +. (0.5 *. Float.Array.get counts tb)) /. total;
              z; reappearance = float_of_int reappearing /. float_of_int nr; depth },
            depths.(k))
          cands
        |> Array.to_list
      end
    let keep_valleys valleys =
      List.filter
        (fun v -> v.z >= valley_z_min && v.share >= valley_share_min && v.share <= valley_share_max)
        valleys
      |> List.stable_sort (fun a b -> compare a.share b.share)
      (* A chain of shares each within [valley_chain] of the one before is one valley, and keeps
         its highest z -- the first reaching it, which is the lower share on a tie *)
      |> List.fold_left
           (fun acc v ->
             match acc with
             | (last, best) :: rest when v.share -. last.share < valley_chain ->
               (v, if v.z > best.z then v else best) :: rest
             | _ -> (v, v) :: acc)
           []
      |> List.rev_map snd
    (* The random state the detector draws from.  It is its own, so that choosing the calibrated
       detector does not move a single draw of any other stream seeded from the same number *)
    let valley_state seed = Random.State.make [| seed; 0x5CA1E |]
    let check_valley_arguments ~resamples ?multiplicities n =
      match multiplicities with
      | None ->
        if resamples < 2 then
          Exception.raise __FUNCTION__ Initialize
            (Printf.sprintf "at least 2 resamples are needed for a standard deviation, not %d"
               resamples)
      | Some m ->
        if Array.length m < 2 then
          Exception.raise __FUNCTION__ Initialize
            (Printf.sprintf "at least 2 resamples are needed for a standard deviation, not %d"
               (Array.length m));
        Array.iter
          (fun a ->
            if Array.length a <> n then
              Exception.raise __FUNCTION__ Initialize
                (Printf.sprintf "a resample has %d multiplicities for %d points" (Array.length a) n))
          m
    let calibrated_troughs ?(seed = 17) ?(pairs_max = 10_000_000) ?(pairs_sample = 2_000_000)
        ?(resamples = 50) ?multiplicities ~metric ~distance ~distance_normalize coords
        inertia_vec =
      let n = Array.length coords in
      check_valley_arguments ~resamples ?multiplicities n;
      if n < 50 then []
      else begin
        let embeds, embed_dist =
          make_embeddings ~metric ~distance ~distance_normalize coords inertia_vec in
        let state = valley_state seed in
        let pairs = valley_pairs ~pairs_max ~pairs_sample state n in
        let multiplicities =
          match multiplicities with
          | Some m -> m
          | None -> valley_multiplicities state resamples n in
        valley_distances embed_dist embeds pairs
        |> calibrate_troughs ~multiplicities pairs |> List.map fst
      end
    let find_valleys ?(verbose = false) ?seed ?pairs_max ?pairs_sample ?resamples ?multiplicities
        ~metric ~distance ~distance_normalize coords inertia_vec =
      let prefix = String.TermIO.grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let candidates =
        calibrated_troughs ?seed ?pairs_max ?pairs_sample ?resamples ?multiplicities ~metric
          ~distance ~distance_normalize coords inertia_vec in
      let kept = keep_valleys candidates in
      if verbose then
        Printf.eprintf "%s %d of %d candidate %s kept%s.\n%!" prefix (List.length kept)
          (List.length candidates)
          (String.pluralize_int "trough" (List.length candidates))
          (List.map
             (fun v ->
               Printf.sprintf "%.1f%% of pairs below (z %.1f, reappearing in %.0f%%)"
                 (100. *. v.share) v.z (100. *. v.reappearance))
             kept
           |> String.concat "; "
           |> fun s -> if s = "" then "" else ": " ^ s);
      kept
    (* HOW FAR APART TWO PARTITIONS OF THE SAME THINGS ARE, given as labellings aligned by
       index.  Three numbers rather than one, because the obvious single number cannot be read
       alone: the adjusted Rand index charges a partition BOTH for splitting a class of the
       other and for merging two of them, so a clean refinement and a real disagreement can
       score alike and the number of clusters ends up tangled into the verdict.  Homogeneity is
       1 when every cluster is pure, whatever their number, and so sees only mixing;
       completeness is 1 when every class is intact, and so falls with every split, cleanly made
       or not.  The pair separates what the index conflates (Rosenberg and Hirschberg 2007) --
       measured on norovirus VP1, one partition scores 0.31 against the genotypes and 0.67
       against the drift variants with homogeneity unmoved at 0.88, which is the whole story *)
    let compare_partitions a b =
      let n = min (Array.length a) (Array.length b) in
      let cell = Hashtbl.create 64 and ca = Hashtbl.create 64 and cb = Hashtbl.create 64 in
      let bump t k = Hashtbl.replace t k (1 + try Hashtbl.find t k with Not_found -> 0) in
      for i = 0 to n - 1 do
        bump cell (a.(i), b.(i)); bump ca a.(i); bump cb b.(i)
      done;
      let choose2 x = let x = float_of_int x in x *. (x -. 1.) /. 2. in
      let sum t f = Hashtbl.fold (fun _ v acc -> acc +. f v) t 0. in
      let index = sum cell choose2 and sa = sum ca choose2 and sb = sum cb choose2
      and total = choose2 n and fn = float_of_int n in
      let expected = if total > 0. then sa *. sb /. total else 0. in
      let maximum = (sa +. sb) /. 2. in
      let ari =
        if maximum -. expected <> 0. then (index -. expected) /. (maximum -. expected) else 1. in
      let log2 x = log x /. log 2. in
      let entropy t = sum t (fun v -> let p = float_of_int v /. fn in -.(p *. log2 p)) in
      let h_a = entropy ca and h_b = entropy cb in
      (* H(A|B) when [on_b], H(B|A) otherwise *)
      let conditional on_b =
        Hashtbl.fold
          (fun (ka, kb) v acc ->
            let p = float_of_int v /. fn
            and marginal =
              float_of_int (Hashtbl.find (if on_b then cb else ca) (if on_b then kb else ka)) in
            acc -. (p *. log2 (float_of_int v /. marginal)))
          cell 0. in
      (* WHICH ARGUMENT IS WHICH MATTERS, and the two are not symmetric: [a] is the partition
         being judged and [b] the reference it is judged against.  Homogeneity asks whether a
         cluster holds one class, so it is the reference that is conditioned on the partition;
         completeness asks whether a class stays in one cluster, so it is the other way about.
         Reversing the arguments exchanges the two *)
      let homogeneity = if h_b > 0. then 1. -. (conditional false /. h_b) else 1.
      and completeness = if h_a > 0. then 1. -. (conditional true /. h_a) else 1. in
      ari, homogeneity, completeness
    (* WHETHER THERE IS ANYTHING HERE TO PARTITION AT ALL.  A set that falls into groups has two
       modes in its pairwise-distance distribution -- one for pairs within a group, one for
       pairs across -- with a valley between them, and a set that does not has one mode and no
       valley.  The radius the search starts from is read off that valley, so its absence is
       already known and is worth SAYING rather than swallowing, because the search returns a
       partition either way: some partition always maximises the criterion, and over a set with
       no groups in it that partition is an artefact of the criterion rather than a finding *)
    let assess_structure ?(verbose = false) ?(seed = 17) ?sample ~what_label ~metric ~distance
        ~distance_normalize coords inertia_vec =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let embeds, embed_dist =
        make_embeddings ~metric ~distance ~distance_normalize coords inertia_vec in
      let state = Random.State.make [| seed |] in
      let valley = valley_radius state embed_dist embeds in
      let dim = intrinsic_dimension ?sample state embed_dist embeds in
      (* A VALLEY OUTSIDE THE DATA IS NOT A VALLEY.  The detector reads a minimum off a smoothed
         histogram and can put it below every distance there is, which is an empty left mode and
         so no mode at all; the giveaway is that a radius like that admits nothing, the leader
         clustering it seeds returning one cluster per sample.  Measuring the share of pairs
         that actually fall below it is what tells a boundary between two modes from a dip in
         the left tail, and only the first is evidence of groups *)
      let below =
        match valley with
        | None -> 0.
        | Some r ->
          let n = Array.length embeds in
          let hits = ref 0 and drawn = ref 0 in
          for _ = 1 to 4000 do
            let i = Random.State.int state n and j = Random.State.int state n in
            if i <> j then begin
              incr drawn;
              if embed_dist embeds.(i) embeds.(j) <= r then incr hits
            end
          done;
          if !drawn > 0 then float_of_int !hits /. float_of_int !drawn else 0. in
      let structured = valley <> None && below >= 0.002 in
      if verbose then begin
        Printf.eprintf
          "%s %s: a neighbourhood occupies about %.1f dimensions of the %d they are stored in.\n%!"
          prefix what_label dim (Float.Array.length inertia_vec);
        (* Said of THESE EMBEDDINGS and not of the data, the two being different claims: a
           projection that does not span the corpus can hide groups that are really there *)
        match valley with
        | Some r when structured ->
          Printf.eprintf
            "%s %s: their distances are two-moded, the valley lying at %.4g with %.1f%% of \
             pairs below it.\n%!"
            prefix what_label r (100. *. below)
        | Some r ->
          Printf.eprintf
            "%s %s: %s.\n%!" prefix what_label
            (red
               (Printf.sprintf
                  "a valley was placed at %.4g but only %.2f%% of pairs lie below it, so there \
                   is no left mode and this embedding shows no groups" r (100. *. below)))
        | None ->
          Printf.eprintf
            "%s %s: %s.\n%!" prefix what_label
            (red "their distances have a single mode, so this embedding shows no groups")
      end;
      structured, dim
    let run_montecarlo
        ?(verbose = false) ?(seed = 17) ?(threads = 1) ?(replicas = 1) ?(exchange_every = 250)
        ?(decades = 3.) ~what_label ~steps ~sample ~temperature ~cooling ~metric ~distance
        ~distance_normalize coords names inertia_vec =
      let n = Array.length coords and d = Float.Array.length inertia_vec in
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let ( .@!() ) = Float.Array.( .@!() ) in
      let embeds, embed_dist =
        make_embeddings ~metric ~distance ~distance_normalize coords inertia_vec in
      let state = Random.State.make [| seed |] in
      (* NO VALLEY MEANS NO GROUPS, and the search is about to return a partition regardless,
         so it is said here rather than left to be inferred from a starting point of one cluster
         per sample.  A radius of zero is the honest consequence and not a fallback: nothing is
         within it, so every sample leads its own cluster and the search begins from the
         partition that assumes nothing *)
      let radius =
        match valley_radius state embed_dist embeds with
        | Some r -> r
        | None ->
          if verbose then
            Printf.eprintf "%s %s: %s.\n%!" prefix what_label
              (String.TermIO.red
                 "the distances have no valley, so there are no groups here to find -- what \
                  follows maximises the criterion over a set that does not cluster");
          0. in
      let start = Array.make n (-1) and leaders = ref [] and n_lead = ref 0 in
      for i = 0 to n - 1 do
        let best = ref (-1) and best_d = ref infinity in
        List.iter
          (fun l ->
            let dd = embed_dist embeds.(i) embeds.(l) in
            if dd <= radius && dd < !best_d then begin best_d := dd; best := l end)
          !leaders;
        if !best >= 0 then start.(i) <- start.(!best)
        else begin leaders := i :: !leaders; start.(i) <- !n_lead; incr n_lead end
      done;
      let n_start = !n_lead and n_samp = min sample n in
      (* ONE CHAIN, run from [from_] at [temp0] for [nsteps], scoring on [samp].  Everything it
         needs is allocated here rather than shared, so that a replica may run in a process of
         its own: the embeddings are read-only and survive the fork, and only the assignment
         travels back.  Returns the best partition the chain saw and the fraction of worsening
         proposals it accepted, the score being left to the caller to recompute exactly -- the
         silhouette of a partition is a different number on a different sample, and the samples
         differ between epochs *)
      let run_chain chain_state samp from_ temp0 nsteps =
        let assign = Array.copy from_ in
        let cent = Array.init n (fun _ -> Float.Array.make d 0.) and size = Array.make n 0 in
        let live = Array.make n 0 and n_live = ref 0 in
        let recentre () =
          Array.fill size 0 n 0;
          for c = 0 to n - 1 do Float.Array.fill cent.(c) 0 d 0. done;
          for i = 0 to n - 1 do
            let c = assign.(i) in
            size.(c) <- size.(c) + 1;
            for k = 0 to d - 1 do
              Float.Array.unsafe_set cent.(c) k (cent.(c).@!(k) +. embeds.(i).@!(k))
            done
          done;
          n_live := 0;
          for c = 0 to n - 1 do
            if size.(c) > 0 then begin
              let f = float_of_int size.(c) in
              for k = 0 to d - 1 do
                Float.Array.unsafe_set cent.(c) k (cent.(c).@!(k) /. f)
              done;
              live.(!n_live) <- c; incr n_live
            end
          done in
        (* A POINT ALONE IN ITS CLUSTER SCORES ZERO, and is counted.  Its own centroid is itself,
           so the quantity the others are scored by is undefined for it, and leaving it out of
           the average instead is what makes shattering pay: every singleton removes a term, so
           a partition of nothing but singletons is judged on the handful of points that are
           not, and can be made to score almost anything.  Zero is the convention and it is the
           right one here, being the value of a point equally close to its own cluster and the
           next, which is what a cluster of one says about it *)
        let silhouette () =
          let tot = ref 0. and cnt = ref 0 in
          Array.iter
            (fun i ->
              let c = assign.(i) in
              if size.(c) >= 2 then begin
                let a = embed_dist embeds.(i) cent.(c) and b = ref infinity in
                for li = 0 to !n_live - 1 do
                  let cc = live.(li) in
                  if cc <> c then begin
                    let s = embed_dist embeds.(i) cent.(cc) in
                    if s < !b then b := s
                  end
                done;
                if !b < infinity then begin
                  tot := !tot +. (!b -. a) /. Float.max a !b; incr cnt
                end
              end else
                incr cnt)
            samp;
          if !cnt > 0 then !tot /. float_of_int !cnt else -1. in
        let merge () =
          let ca = live.(Random.State.int chain_state !n_live) in
          let cb = ref (-1) and best = ref infinity in
          for li = 0 to !n_live - 1 do
            let cc = live.(li) in
            if cc <> ca then begin
              let s = embed_dist cent.(ca) cent.(cc) in
              if s < !best then begin best := s; cb := cc end
            end
          done;
          if !cb >= 0 then
            for i = 0 to n - 1 do if assign.(i) = ca then assign.(i) <- !cb done in
        let split () =
          let ca = live.(Random.State.int chain_state !n_live) in
          if size.(ca) >= 6 then begin
            let members = Array.make size.(ca) 0 and m = ref 0 in
            for i = 0 to n - 1 do
              if assign.(i) = ca then begin members.(!m) <- i; incr m end
            done;
            let k1 = members.(Random.State.int chain_state !m)
            and k2 = members.(Random.State.int chain_state !m) in
            if k1 <> k2 then begin
              let fresh = ref 0 in
              while size.(!fresh) > 0 do incr fresh done;
              let m1 = Float.Array.make d 0. and m2 = Float.Array.make d 0. in
              Float.Array.blit embeds.(k1) 0 m1 0 d;
              Float.Array.blit embeds.(k2) 0 m2 0 d;
              for _ = 1 to 5 do
                let n1 = ref 0 and n2 = ref 0 in
                for j = 0 to !m - 1 do
                  let i = members.(j) in
                  if embed_dist embeds.(i) m1 <= embed_dist embeds.(i) m2 then begin
                    assign.(i) <- ca; incr n1
                  end else begin assign.(i) <- !fresh; incr n2 end
                done;
                if !n1 > 0 && !n2 > 0 then begin
                  Float.Array.fill m1 0 d 0.;
                  Float.Array.fill m2 0 d 0.;
                  for j = 0 to !m - 1 do
                    let i = members.(j) in
                    let t = if assign.(i) = ca then m1 else m2 in
                    for k = 0 to d - 1 do
                      Float.Array.unsafe_set t k (t.@!(k) +. embeds.(i).@!(k))
                    done
                  done;
                  let f1 = float_of_int !n1 and f2 = float_of_int !n2 in
                  for k = 0 to d - 1 do
                    Float.Array.unsafe_set m1 k (m1.@!(k) /. f1);
                    Float.Array.unsafe_set m2 k (m2.@!(k) /. f2)
                  done
                end
              done
            end
          end in
        recentre ();
        let start_score = silhouette () in
        let cur = ref start_score and temp = ref temp0 in
        let best = Array.copy assign and best_score = ref start_score and saved = Array.make n 0 in
        (* Proposals that would worsen the score, and how many of them were taken anyway.  This
           ratio is what the temperature actually controls, and the caller tunes the ladder on
           it *)
        let worse = ref 0 and worse_taken = ref 0 in
        for _ = 1 to nsteps do
          Array.blit assign 0 saved 0 n;
          if Random.State.float chain_state 1. < 0.5 && !n_live > 2 then merge () else split ();
          recentre ();
          let proposed = silhouette () in
          let take =
            (* A PROPOSAL THAT CHANGES NOTHING IS NOT A WORSENING ONE, and counting it as one
               corrupts the very rate the ladder is tuned on: exp(0/T) is 1 at every
               temperature, so a tie is always taken, and a partition offering many of them --
               which one made mostly of singletons does, the score being blind to where a
               singleton moves -- reports an acceptance near 100% however cold the chain is.
               The ladder then reads itself as far too hot and slides until it is not tuning
               anything *)
            if proposed >= !cur then true
            else begin
              incr worse;
              let ok = Random.State.float chain_state 1. < exp ((proposed -. !cur) /. !temp) in
              if ok then incr worse_taken;
              ok
            end in
          if take then begin
            cur := proposed;
            if proposed > !best_score then begin
              best_score := proposed; Array.blit assign 0 best 0 n
            end
          end else begin
            Array.blit saved 0 assign 0 n; recentre ()
          end;
          temp := !temp *. cooling
        done;
        best, if !worse > 0 then float_of_int !worse_taken /. float_of_int !worse else 0. in
      (* A LADDER OF CHAINS, RESET TO THE BEST OF THEM AT EVERY EPOCH.  A single chain has to
         be cold enough to refine and hot enough to escape, and cannot be both; a ladder of
         temperatures can, the hot replicas wandering while the cold ones settle.  What
         propagates a good partition between them is not the exchange of classical parallel
         tempering but plain selection: every epoch, each replica restarts from the best
         partition ANY of them has found, at its own temperature.  That abandons detailed
         balance and with it any claim to be sampling a distribution -- which we do not want,
         since the object here is the best partition and not the ensemble -- and in exchange
         it propagates a good state to every replica in one epoch rather than diffusing it
         along the ladder one neighbour at a time -- which matters because the run is a few
         tens of epochs long and a diffusing state may not cross the ladder within it.
         THE LADDER TUNES ITSELF ON THE ACCEPTANCE RATE, which is the fraction of WORSENING
         proposals a chain takes anyway.  That is the quantity the temperature controls: a chain
         taking none of them is greedy, one taking all of them is a random walk, and the useful
         range lies between.  It is also monotone in T, so a ladder that brackets the target
         band stops moving, and only a ladder entirely above or entirely below it slides.
         Which rung won the epoch will not serve instead, however the wins are weighted: a chain
         at T -> 0 takes only improving moves, so over a single epoch from a shared start it
         out-climbs a chain that spends the same steps wandering, whatever the landscape.  That
         bias is towards cold and it is unbounded -- a ladder following it descends until it is
         doing greedy descent, at which point the temperature it reports is an estimate of
         nothing.  Where this ladder stops is reported with the rates, so that the temperature
         need not be guessed before the run that would measure it.  Each epoch is one parallel
         map, and the only thing that travels back from a replica is an assignment.
         THE RUNGS ARE SPACED IN THE LOGARITHM of the temperature and not in the temperature.
         Metropolis compares a score difference to T as a ratio, so it is the ratio between two
         rungs that decides how differently they behave: an arithmetic ladder crowds every rung
         into one regime, a geometric one visits several.  Spanning [decades] decades over the
         replicas puts the rungs at 0.1, 0.01, 0.001 and so on *)
      let replicas = max 1 replicas in
      let epochs = if replicas = 1 then 1 else max 1 (steps / exchange_every) in
      let per_epoch = if replicas = 1 then steps else steps / epochs in
      let rung = if replicas > 1 then 10. ** (decades /. float_of_int (replicas - 1)) else 1. in
      let ladder = Array.init replicas (fun r -> temperature *. (rung ** float_of_int r)) in
      (* The band the ladder is asked to straddle, in fraction of worsening proposals accepted.
         It is wide because the ladder only has to bracket it, not to sit in it: a rung below
         the floor is doing greedy descent and one above the ceiling is walking at random *)
      let accept_min = 0.1 and accept_max = 0.5 in
      let seeds = Array.init replicas (fun r -> seed + 1013 * (r + 1)) in
      let wins = Array.make replicas 0 and rates = Array.make replicas 0. in
      (* The silhouette over every point rather than over the sample the chains score on.
         Selection between epochs is decided on this, because the sample moves from epoch to
         epoch and the differences being judged are thousandths, well inside its noise: a
         sampled comparison would discard the better partition about as often as it kept it.
         One exact evaluation costs what a single Metropolis step costs, against the hundreds
         of them each replica runs per epoch, so exactness here is free *)
      let full_silhouette assign =
        let cent = Array.init n (fun _ -> Float.Array.make d 0.) and size = Array.make n 0 in
        for i = 0 to n - 1 do
          let c = assign.(i) in
          size.(c) <- size.(c) + 1;
          for k = 0 to d - 1 do
            Float.Array.unsafe_set cent.(c) k (cent.(c).@!(k) +. embeds.(i).@!(k))
          done
        done;
        let live = Array.make n 0 and n_live = ref 0 in
        for c = 0 to n - 1 do
          if size.(c) > 0 then begin
            let f = float_of_int size.(c) in
            for k = 0 to d - 1 do Float.Array.unsafe_set cent.(c) k (cent.(c).@!(k) /. f) done;
            live.(!n_live) <- c; incr n_live
          end
        done;
        let tot = ref 0. and cnt = ref 0 in
        for i = 0 to n - 1 do
          let c = assign.(i) in
          if size.(c) >= 2 then begin
            let a = embed_dist embeds.(i) cent.(c) and b = ref infinity in
            for li = 0 to !n_live - 1 do
              let cc = live.(li) in
              if cc <> c then begin
                let s = embed_dist embeds.(i) cent.(cc) in
                if s < !b then b := s
              end
            done;
            if !b < infinity then begin
              tot := !tot +. (!b -. a) /. Float.max a !b; incr cnt
            end
          end else
            incr cnt
        done;
        if !cnt > 0 then !tot /. float_of_int !cnt else -1. in
      let best = ref (Array.copy start) and best_score = ref (full_silhouette start) in
      if verbose then
        Printf.eprintf "%s %s: starting from %d clusters, %d %s over %d %s of %d steps.\n%!"
          prefix what_label n_start replicas (String.pluralize_int "replica" replicas)
          epochs (String.pluralize_int "epoch" epochs) per_epoch;
      for epoch = 1 to epochs do
        (* Every replica starts the epoch from the best partition found so far, and differs
           from its neighbours only in temperature.  The scoring sample is redrawn here and
           shared by all of them: a chain scored on one fixed sample for thousands of steps
           optimises that sample rather than the corpus, and the gap it opens between the two is
           of the order of the improvements being chased -- four hundredths, where 250 points
           stand for 3027 -- whereas a sample that moves under the chain cannot be learnt at
           all.  It is redrawn per epoch rather than per step
           because a chain needs a fixed yardstick over the run of moves it is comparing, and
           shared across replicas so that they all face the same problem *)
        let from_ = !best in
        let samp = Array.init n_samp (fun _ -> Random.State.int state n) in
        let cand = Array.make replicas from_ in
        let take r chain_best rate = cand.(r) <- chain_best; rates.(r) <- rate in
        if replicas = 1 then begin
          let chain_best, rate =
            run_chain (Random.State.make [| seeds.(0) + epoch |]) samp from_ ladder.(0)
              per_epoch in
          take 0 chain_best rate
        end else begin
          let next = ref 0 in
          Processes.Parallel.process_stream_chunkwise
            (fun () -> if !next < replicas then begin let r = !next in incr next; r end
                       else raise End_of_file)
            (fun r ->
              let chain_best, rate =
                run_chain (Random.State.make [| seeds.(r) + epoch |]) samp from_ ladder.(r)
                  per_epoch in
              r, chain_best, rate)
            (fun (r, chain_best, rate) -> take r chain_best rate)
            (max 1 (min threads replicas))
        end;
        (* Selection, on the exact score, so that [best_score] never falls *)
        let exact = Array.map full_silhouette cand in
        let arg = ref 0 in
        Array.iteri (fun r s -> if s > exact.(!arg) then arg := r) exact;
        if exact.(!arg) > !best_score then begin
          best_score := exact.(!arg); best := cand.(!arg); wins.(!arg) <- wins.(!arg) + 1
        end;
        if verbose then
          Printf.eprintf "%s %s: epoch %d/%d, best silhouette %.4f, epoch won at T=%.3g.\n%!"
            prefix what_label epoch epochs !best_score ladder.(!arg);
        (* The ladder slides only when it lies wholly outside the band, which cannot happen once
           it brackets it, acceptance being monotone in the temperature.  Half a rung at a time,
           so that where it comes to rest is not coarser than the spacing it came from *)
        if replicas > 1 then begin
          let lo = Array.fold_left Float.min infinity rates
          and hi = Array.fold_left Float.max neg_infinity rates in
          let shift =
            if lo > accept_max then 1. /. sqrt rung
            else if hi < accept_min then sqrt rung else 1. in
          if shift <> 1. then
            for r = 0 to replicas - 1 do ladder.(r) <- ladder.(r) *. shift done
        end
      done;
      if verbose && replicas > 1 then begin
        let acc = Buffer.create 64 in
        Array.iteri
          (fun r w ->
            Buffer.add_string acc
              (Printf.sprintf " T=%.3g:%.0f%%/%d" ladder.(r) (100. *. rates.(r)) w))
          wins;
        Printf.eprintf
          "%s %s: ladder settled at, with the worsening moves each rung accepted and the epochs \
           it won:%s.\n%!"
          prefix what_label (Buffer.contents acc)
      end;
      let assign = !best in
      (* The centroids of the partition being returned, for the representatives below *)
      let cent = Array.init n (fun _ -> Float.Array.make d 0.) and size = Array.make n 0 in
      Array.fill size 0 n 0;
      for i = 0 to n - 1 do
        let c = assign.(i) in
        size.(c) <- size.(c) + 1;
        for k = 0 to d - 1 do
          Float.Array.unsafe_set cent.(c) k (cent.(c).@!(k) +. embeds.(i).@!(k))
        done
      done;
      let n_live = ref 0 in
      for c = 0 to n - 1 do
        if size.(c) > 0 then begin
          let f = float_of_int size.(c) in
          for k = 0 to d - 1 do
            Float.Array.unsafe_set cent.(c) k (cent.(c).@!(k) /. f)
          done;
          incr n_live
        end
      done;
      let cur = best_score in
      let rep = Array.make n (-1) and rep_d = Array.make n infinity in
      for i = 0 to n - 1 do
        let c = assign.(i) in
        let dd = embed_dist embeds.(i) cent.(c) in
        if dd < rep_d.(c) then begin rep_d.(c) <- dd; rep.(c) <- i end
      done;
      let rep_orig = Array.init n (fun i -> rep.(assign.(i))) in
      let metric_str = Space.Distance.Metric.to_string metric
      and distance_str = Space.Distance.to_string distance in
      Printf.printf
        "=== Clustering of %s: Monte-Carlo \
         (steps=%d, sample=%d, temperature=%.15g, cooling=%.15g, metric=%s, distance=%s, \
         D=%d) ===\n\
         # n=%d n_clusters=%d silhouette=%.15g started_from=%d\n\
         name\trepresentative\tstatus\n"
        what_label steps (min sample n) temperature cooling metric_str distance_str d
        n !n_live !cur n_start;
      for i = 0 to n - 1 do
        Printf.printf "%s\t%s\t%s\n"
          names.(i) names.(rep_orig.(i)) (if rep_orig.(i) = i then "rep" else "abs")
      done;
      if verbose then
        Printf.eprintf "%s Clustering of %s done. %d %s, silhouette %.4f.\n%!"
          prefix what_label !n_live (String.pluralize_int "cluster" !n_live) !cur;
      rep_orig
  end: sig
    module Algorithm:
      sig
        type t =
          | Greedy
          | Hdbscan
          | Montecarlo
        val of_string: string -> t
        val to_string: t -> string
      end
    module GreedyEpsilon:
      sig
        type t =
          | FirstNN
          | Density
          | Adaptive
        val of_string: string -> t
        val to_string: t -> string
      end
    module Order:
      sig
        type t =
          | Inertia
          | FirstNN
          | Density
        val of_string: string -> t
        val to_string: t -> string
      end
    val kneedle_elbow: (int -> float) -> int -> int
    val run_greedy:
      ?verbose:bool ->
      what_label:string ->
      epsilon_:GreedyEpsilon.t ->
      order_:Order.t ->
      density_sample_number:int ->
      index_type:Interfaiss.Type.t ->
      metric:Space.Distance.Metric.t ->
      distance:Space.Distance.t ->
      distance_normalize:bool ->
      Float.Array.t array ->
      string array ->
      Float.Array.t ->
      int array
    module HdbscanMstMode:
      sig
        type t =
          | Auto
          | Sparse
          | Dense
        val of_string: string -> t
        val to_string: t -> string
      end
    module Hdbscan:
      sig
        module LengthsMode:
          sig
            type t =
              | Persistence
              | Mreach
            val of_string: string -> t
            val to_string: t -> string
          end
        val make_tree:
          ?threads:int -> ?verbose:bool ->
          ?mst_mode:HdbscanMstMode.t ->
          ?num_neighbors:int ->
          ?index_type:Interfaiss.Type.t ->
          ?lengths_mode:LengthsMode.t ->
          min_cluster_size:int -> min_samples:int ->
          Matrix.Base.t -> Trees.Newick.t
        val make_clusters:
          ?threads:int -> ?verbose:bool ->
          ?mst_mode:HdbscanMstMode.t ->
          ?num_neighbors:int ->
          ?index_type:Interfaiss.Type.t ->
          min_cluster_size:int -> min_samples:int ->
          Float.Array.t array -> int array
      end
    val valley_radius:
      ?pairs:int -> ?bins:int -> ?smooth:int -> ?floor_:float ->
      Random.State.t ->
      (Float.Array.t -> Float.Array.t -> float) ->
      Float.Array.t array ->
      float option
    (* A valley in the pairwise distances, as the calibrated detector reports one *)
    type valley_t = { radius: float; share: float; z: float; reappearance: float; depth: float }
    (* The troughs of a histogram once smoothed and simplified, each as its bin, the bins of the
       peaks either side of it, its prominence and its relative depth *)
    val troughs_of_counts: Float.Array.t -> (int * int * int * float * float) list
    (* The candidates that pass the detector's thresholds, a chain of near-identical shares
       counting as one valley *)
    val keep_valleys: valley_t list -> valley_t list
    (* Every candidate trough of the pairwise distances, calibrated against resamples of the points
       and not yet thresholded.  [multiplicities], when given, replaces the resamples drawn *)
    val calibrated_troughs:
      ?seed:int ->
      ?pairs_max:int ->
      ?pairs_sample:int ->
      ?resamples:int ->
      ?multiplicities:int array array ->
      metric:Space.Distance.Metric.t ->
      distance:Space.Distance.t ->
      distance_normalize:bool ->
      Float.Array.t array ->
      Float.Array.t ->
      valley_t list
    val find_valleys:
      ?verbose:bool ->
      ?seed:int ->
      ?pairs_max:int ->
      ?pairs_sample:int ->
      ?resamples:int ->
      ?multiplicities:int array array ->
      metric:Space.Distance.Metric.t ->
      distance:Space.Distance.t ->
      distance_normalize:bool ->
      Float.Array.t array ->
      Float.Array.t ->
      valley_t list
    (* Adjusted Rand index, homogeneity and completeness of two labellings of the same items,
       given aligned by index.  Three numbers because one cannot be read alone.  The first array
       is the partition being judged, the second the reference: reversing them exchanges
       homogeneity and completeness *)
    val compare_partitions: string array -> string array -> float * float * float
    (* Whether the pairwise distances are two-moded, which is to say whether there are groups
       here at all, and the dimension a neighbourhood occupies *)
    val assess_structure:
      ?verbose:bool ->
      ?seed:int ->
      ?sample:int ->
      what_label:string ->
      metric:Space.Distance.Metric.t ->
      distance:Space.Distance.t ->
      distance_normalize:bool ->
      Float.Array.t array ->
      Float.Array.t ->
      bool * float
    val run_montecarlo:
      ?verbose:bool ->
      ?seed:int ->
      ?threads:int ->
      ?replicas:int ->
      ?exchange_every:int ->
      ?decades:float ->
      what_label:string ->
      steps:int ->
      sample:int ->
      temperature:float ->
      cooling:float ->
      metric:Space.Distance.Metric.t ->
      distance:Space.Distance.t ->
      distance_normalize:bool ->
      Float.Array.t array ->
      string array ->
      Float.Array.t ->
      int array
    val run_hdbscan:
      ?verbose:bool ->
      ?threads:int ->
      what_label:string ->
      min_cluster_size:int ->
      ?min_samples:int ->
      mst_mode:HdbscanMstMode.t ->
      ?num_neighbors:int ->
      index_type:Interfaiss.Type.t ->
      metric:Space.Distance.Metric.t ->
      distance:Space.Distance.t ->
      distance_normalize:bool ->
      Float.Array.t array ->
      string array ->
      Float.Array.t ->
      int array
  end
)

