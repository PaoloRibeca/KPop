(*
    SparseNJ.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    SparseNJ.ml implements a sparse, locally-regularised variant of
    neighbour-joining.  Two differences from classical NJ:
      (1) at each merge iteration, only pairs (i, j) where j appears
          among i's K nearest active neighbours (or vice versa,
          controlled by the symmetry flag) are considered as merge
          candidates.  At fixed K this lifts the per-iteration cost from
          the classical O(n_act^2) to O(n_act K).
      (2) the row sum used in the NJ Q-formula is replaced by a local
          K-NN estimator, r_hat(i) = (n_act - 1) / K * sum over i's K
          nearest active neighbours of d(i, j).  On non-additive distance
          matrices this regularises against noise contributions from
          distant leaves and empirically beats classical NJ's full row
          sum by a measurable margin on the protein-k = 5 evaluation
          dataset (see ../../DocsYard/KPop/PhyloSplits/
          KPop-PhyloSplits-Evaluation.tex \S 9.2).

    The primary entry point [compute] returns a Trees.Newick.t directly:
    each merge of two active clusters yields a Newick join node carrying
    the NJ-computed branch lengths above its two children.  The final
    three-way merge of the last three active clusters produces an
    unrooted-tree-style trifurcating root.

    Two implementation modes are available via [Mode.t]:
      - [Dense] (default, validated): a persistent symmetric distance
        cache stores the current Saitou-Nei distance between every pair
        of active clusters, updated symmetrically at each merge.  K-NN
        selection at every iteration brute-scans active candidates and
        ranks them by their current Saitou-Nei distance.  Time
        O(n^3); memory O(n^2).  Reproduces the test prototype
        (test/Trees/sparse_nj.py) byte-identical on prot_k5_kt0.035.
      - [Subquadratic] (experimental): a per-cluster K-NN list with
        a reverse index, FAISS-driven candidate expansion when a
        list shrinks below K, and an explicit Saitou-Nei recursion
        over the merge tree for distances not in the K-NN cache.
        Time O(n K^2 + n K log n); memory O(n K).  Reproduces
        the [Dense] result whenever the FAISS expansion captures the
        true Saitou-Nei K-NN of the merged cluster -- which requires
        the centroid embedding to approximate the Saitou-Nei
        neighbourhood structure well.  Validate empirically before
        switching to [Subquadratic] as the default.

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

module Trees = BiOCamLib.Trees
open BiOCamLib.Better

include (
  struct
    (* Row-sum estimator used in the NJ Q-formula.  [Knn] approximates
       r(i) as (n_act - 1) / K times the sum of i's K-NN distances;
       [Topk] is the same up to a constant scaling (uses the mean
       and rescales by n_act - 1); [Full] uses the exact global row
       sum, recomputed lazily from the Saitou-Nei distance cache
       (with centroid-Euclidean fallback) at each iteration. *)
    module RowSum =
      struct
        type t =
          | Knn
          | Topk
          | Full
        let of_string = function
          | "knn" -> Knn
          | "topk" -> Topk
          | "full" -> Full
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ row-sum estimator" s
        let to_string = function
          | Knn -> "knn"
          | Topk -> "topk"
          | Full -> "full"
      end
    (* K-NN symmetry policy.  [One] includes pair (i, j) in the
       candidate set if either j is in i's K-NN list or i is in j's;
       [Both] requires both.  One-sided is consistently better on
       real, non-additive distance matrices. *)
    module Symmetry =
      struct
        type t =
          | One
          | Both
        let of_string = function
          | "one" -> One
          | "both" -> Both
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ symmetry" s
        let to_string = function
          | One -> "one"
          | Both -> "both"
      end
    (* Distance model fed to the NJ Q-formula and row-sum estimator.
       [SaitouNei] is the exact NJ update d(u, x) = 1/2 (d(i, x) +
       d(j, x) - d(i, j)), computed by memoised recursion over the
       merge history -- the validated default, but whose recursion
       tiles the full O(n^2) distance-pair space (see
       PhyloSplits-Subquadratic.tex).  [Centroid] is an O(d) geometric
       proxy: the Euclidean distance between the size-weighted cluster
       centroids already maintained for K-NN retrieval, with no
       recursion (and so no n^2 cache fill).  [Hyperbolic] is also an
       O(d) recursion-free proxy, but tracks Saitou-Nei far better than
       the centroid: leaf positions are lifted onto the upper
       hyperboloid and each merge places the new cluster on the
       geodesic between its parents at the NJ branch length, so
       hyperbolic distances follow the exact update closely across the
       whole run (Phase 0: Pearson 0.97+ vs Saitou-Nei, where the
       centroid collapses to ~0.33 mid-run).  [HyperbolicFrechet] is
       the same hyperboloid distance, but places the merged cluster at
       the size-weighted hyperbolic barycenter (geodesic point at
       fraction size_j/(size_i+size_j) of d(i,j) from i) instead of at
       the NJ branch length b_i.  The barycenter is bounded in
       [0, d(i,j)] and depends only on cluster sizes, not on the noisy
       b_i, so it does not amplify branch-length noise on a steeply
       curved manifold --- the failure mode that makes plain
       [Hyperbolic] crash at large radial scale.  [Hybrid] uses the
       centroid proxy for the cheap bulk work (retrieval, K-NN lists,
       row sums, heap ordering) but confirms the winning merge by
       exact Saitou-Nei over the popped shortlist only, and stores the
       exact pairwise distance so the recursion stays consistent --- a
       bounded amount of exact recursion per merge, intended to recover
       the proxy's quality gap while keeping coverage sub-quadratic
       (whether it does is the empirical question the cache-coverage
       counter answers).  All non-default models experimental and only
       wired into the RP-forest mode. *)
    module Distance =
      struct
        type t =
          | SaitouNei
          | Centroid
          | Hyperbolic
          | HyperbolicFrechet
          | Hybrid
        let of_string = function
          | "saitou-nei" | "saitou_nei" | "saitounei" | "nj" -> SaitouNei
          | "centroid" -> Centroid
          | "hyperbolic" -> Hyperbolic
          | "hyperbolic-frechet" | "hyperbolic_frechet" | "hyp-frechet" ->
            HyperbolicFrechet
          | "hybrid" | "centroid-exact" | "hybrid-centroid" -> Hybrid
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ distance model" s
        let to_string = function
          | SaitouNei -> "saitou-nei"
          | Centroid -> "centroid"
          | Hyperbolic -> "hyperbolic"
          | HyperbolicFrechet -> "hyperbolic-frechet"
          | Hybrid -> "hybrid"
      end
    (* Implementation mode.  [Dense] is the validated reference; see
       compute_dense.  [Subquadratic] is the experimental
       O(n K^2 + n K log n) centroid-based path; see compute_subquadratic.
       [Hyperbolic] is the geodesic-tracked hyperbolic-embedding path
       motivated by the Phase 0 result showing hyperbolic geodesic
       placement tracks Saitou-Nei distances across the whole NJ run
       (whereas size-weighted centroids drift catastrophically in
       middle iterations); see compute_hyperbolic. *)
    module Mode =
      struct
        type t =
          | Dense
          | Subquadratic
          | Hyperbolic
          | CoverTree
          | PeriodicRebuild
          | RPForestRebuild
        let of_string = function
          | "dense" -> Dense
          | "subquadratic" -> Subquadratic
          | "hyperbolic" -> Hyperbolic
          | "cover-tree" | "cover_tree" | "covertree" -> CoverTree
          | "periodic-rebuild" | "periodic_rebuild" | "periodicrebuild" ->
            PeriodicRebuild
          | "rp-forest" | "rp_forest" | "rpforest"
          | "rp-forest-rebuild" | "rp_forest_rebuild" ->
            RPForestRebuild
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ mode" s
        let to_string = function
          | Dense -> "dense"
          | Subquadratic -> "subquadratic"
          | Hyperbolic -> "hyperbolic"
          | CoverTree -> "cover-tree"
          | PeriodicRebuild -> "periodic-rebuild"
          | RPForestRebuild -> "rp-forest"
      end
    (* Squared L2 distance between two Float.Array vectors, accumulated
       in float64. *)
    let sq_dist a b =
      let dim = Float.Array.length a in
      let acc = ref 0. in
      for k = 0 to dim - 1 do
        let delta =
          Float.Array.unsafe_get a k -. Float.Array.unsafe_get b k in
        acc := !acc +. delta *. delta
      done;
      !acc
    let eucl_dist a b = sqrt (sq_dist a b) [@@inline]
    let cache_key i j = if i < j then (i, j) else (j, i) [@@inline]
    (* Hyperboloid model, shared by [compute_hyperbolic] (where it is
       the retrieval metric) and by the [Hyperbolic] distance proxy of
       [compute_rp_forest_rebuild] (where it is the Q-formula distance).

       A point p in H^d is stored as a (d + 1)-vector with p.(0) the
       time-like coord (positive) and p.(1..d) the spatial coords,
       satisfying the Lorentz constraint -p.(0)^2 + sum p.(k)^2 = -1.
       The Lorentz inner product is <p, q>_L = -p.(0) q.(0) + sum
       p.(k) q.(k); hyperbolic distance is acosh(-<p, q>_L).

       Leaf positions are lifted from the Euclidean embedding via
           p_i -> (cosh(s |p_i|), sinh(s |p_i|) * p_i / |p_i|)
       with [scale] = s a tunable radial scale; per-merge the new
       cluster is placed at the geodesic point at hyperbolic distance
       b_i from p_i (toward p_j).  Under tree-metric additivity this
       yields d_H(p_u, p_x) = d_NJ(u, x) exactly; on real data it is
       approximate but tracks Saitou-Nei closely (Phase 0: Pearson
       0.97+; see PhyloSplits-Subquadratic.tex Path 3). *)
    let lorentz_inner p q =
      let acc = ref (~-. (Float.Array.unsafe_get p 0 *. Float.Array.unsafe_get q 0)) in
      let dim = Float.Array.length p in
      for k = 1 to dim - 1 do
        acc := !acc +. Float.Array.unsafe_get p k *. Float.Array.unsafe_get q k
      done;
      !acc
    let hyp_dist p q =
      let inner = ~-. (lorentz_inner p q) in
      let inner = max inner 1. in
      acosh inner
    let hyp_lift ~scale x =
      let d = Float.Array.length x in
      let norm =
        let acc = ref 0. in
        for k = 0 to d - 1 do
          let v = Float.Array.unsafe_get x k in
          acc := !acc +. v *. v
        done;
        sqrt !acc in
      let safe_norm = if norm < 1e-12 then 1e-12 else norm in
      let s = scale *. safe_norm in
      let cosh_s = cosh s and sinh_s = sinh s in
      let coeff = sinh_s /. safe_norm in
      let out = Float.Array.create (d + 1) in
      Float.Array.unsafe_set out 0 cosh_s;
      for k = 0 to d - 1 do
        Float.Array.unsafe_set out (k + 1) (coeff *. Float.Array.unsafe_get x k)
      done;
      out
    (* Geodesic point at hyperbolic distance [t] from [p] toward [q] on
       the hyperboloid.  Both p and q must satisfy the Lorentz
       constraint -<x, x>_L = 1.  Output also satisfies it. *)
    let hyp_geodesic p q t =
      let d = hyp_dist p q in
      let dim = Float.Array.length p in
      if d < 1e-12 then Float.Array.copy p
      else begin
        let cosh_d = cosh d and sinh_d = sinh d in
        let cosh_t = cosh t and sinh_t = sinh t in
        let inv_sinh_d = 1. /. sinh_d in
        let out = Float.Array.create dim in
        for k = 0 to dim - 1 do
          let pk = Float.Array.unsafe_get p k and qk = Float.Array.unsafe_get q k in
          let v_k = (qk -. cosh_d *. pk) *. inv_sinh_d in
          Float.Array.unsafe_set out k (cosh_t *. pk +. sinh_t *. v_k)
        done;
        out
      end
    (* Compute the candidate pair set: pairs (i, j) where j is among
       i's K-NN (One-sided) or both i and j are in each other's K-NN
       (Both).  Returns a list of (i, j) with i < j, deduplicated;
       entries pointing at inactive slots are filtered. *)
    let candidate_pairs sym nbrs active =
      let n = Array.length nbrs in
      let seen = Hashtbl.create (16 * n) in
      let res = ref [] in
      let add a b =
        if a <> b && active.(a) && active.(b) then begin
          let lo, hi = if a < b then a, b else b, a in
          if not (Hashtbl.mem seen (lo, hi)) then begin
            Hashtbl.add seen (lo, hi) ();
            res := (lo, hi) :: !res
          end
        end in
      (match sym with
       | Symmetry.One ->
         for i = 0 to n - 1 do
           if active.(i) then
             Array.iter (fun j -> add i j) nbrs.(i)
         done
       | Symmetry.Both ->
         let nbrs_sets = Array.map (fun arr ->
           let h = Hashtbl.create (Array.length arr) in
           Array.iter (fun j -> Hashtbl.replace h j ()) arr;
           h) nbrs in
         for i = 0 to n - 1 do
           if active.(i) then
             Array.iter
               (fun j ->
                 if Hashtbl.mem nbrs_sets.(j) i then
                   add i j)
               nbrs.(i)
         done);
      !res
    (* Private min-heap of (q', vi, vj, i, j) entries ordered by q'.
       Used by the heap-driven main loop to avoid the per-merge
       O(n_act * K) scan over candidate pairs.  Stored q' is an
       approximation (Q' = d(i,j) - r̂(i) - r̂(j), no n_act-dependent
       coefficient); freshness is validated on pop by comparing
       stored versions against the per-cluster version array. *)
    module QHeap =
      struct
        type entry = {
          q_prime: float;
          vi: int;
          vj: int;
          i: int;
          j: int
        }
        let dummy = { q_prime = 0.; vi = 0; vj = 0; i = 0; j = 0 }
        type t = {
          mutable arr: entry array;
          mutable n: int
        }
        let create ?(cap = 64) () = { arr = Array.make cap dummy; n = 0 }
        let length h = h.n
        let ensure h cap =
          if cap > Array.length h.arr then begin
            let new_cap = max cap (2 * Array.length h.arr) in
            let new_arr = Array.make new_cap dummy in
            Array.blit h.arr 0 new_arr 0 h.n;
            h.arr <- new_arr
          end
        let push h e =
          ensure h (h.n + 1);
          h.arr.(h.n) <- e;
          h.n <- h.n + 1;
          let i = ref (h.n - 1) in
          let going = ref true in
          while !going && !i > 0 do
            let p = (!i - 1) / 2 in
            if h.arr.(!i).q_prime < h.arr.(p).q_prime then begin
              let t = h.arr.(!i) in
              h.arr.(!i) <- h.arr.(p);
              h.arr.(p) <- t;
              i := p
            end else going := false
          done
        let pop h =
          if h.n = 0 then None
          else begin
            let top = h.arr.(0) in
            h.n <- h.n - 1;
            if h.n > 0 then begin
              h.arr.(0) <- h.arr.(h.n);
              let i = ref 0 in
              let going = ref true in
              while !going do
                let l = 2 * !i + 1 and r = 2 * !i + 2 in
                let smallest = ref !i in
                if l < h.n
                   && h.arr.(l).q_prime < h.arr.(!smallest).q_prime then
                  smallest := l;
                if r < h.n
                   && h.arr.(r).q_prime < h.arr.(!smallest).q_prime then
                  smallest := r;
                if !smallest = !i then going := false
                else begin
                  let t = h.arr.(!i) in
                  h.arr.(!i) <- h.arr.(!smallest);
                  h.arr.(!smallest) <- t;
                  i := !smallest
                end
              done
            end;
            Some top
          end
      end
    (* === Dense reference implementation === *)
    let compute_dense ?(verbose = false)
        ?(index_type = Interfaiss.Type.of_string "hnsw(32)")
        ?(k_nn = 10) ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n = Array.length names in
      if Array.length data <> n then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sample count mismatch: %d names vs %d data rows"
             n (Array.length data));
      if n < 3 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sparse-NJ requires at least 3 leaves (got %d)" n);
      let dim = Float.Array.length data.(0) in
      let trees = Array.init n (fun i -> Trees.Newick.leaf names.(i)) in
      let active = Array.make n true in
      let n_active = ref n in
      let size = Array.make n 1 in
      let embedding = Array.init n (fun i -> Float.Array.copy data.(i)) in
      let dist_cache : (int * int, float) Hashtbl.t = Hashtbl.create (n * 8) in
      let dist_of i j =
        if i = j then 0.
        else
          let key = cache_key i j in
          match Hashtbl.find_opt dist_cache key with
          | Some d -> d
          | None ->
            eucl_dist embedding.(i) embedding.(j) in
      let touch_dist i j =
        if i <> j then begin
          let key = cache_key i j in
          if not (Hashtbl.mem dist_cache key) then
            Hashtbl.add dist_cache key (eucl_dist embedding.(i) embedding.(j))
        end in
      let bootstrap_cache () =
        if n >= 2 then begin
          let ba = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout n dim in
          for i = 0 to n - 1 do
            let row = embedding.(i) in
            for kk = 0 to dim - 1 do
              ba.{i, kk} <- Float.Array.unsafe_get row kk
            done
          done;
          let index = Interfaiss.create ~index_type dim in
          Interfaiss.train index ba;
          Interfaiss.add index ba;
          let k_query = min n (k_nn + 1) in
          let offsets, _ = Interfaiss.query index ba k_query in
          Interfaiss.delete index;
          for i = 0 to n - 1 do
            for kk = 0 to k_query - 1 do
              let j = Int64.to_int offsets.{i, kk} in
              if j >= 0 && j < n && j <> i then touch_dist i j
            done
          done
        end in
      let nbrs = Array.make n [||] in
      let dists = Array.make n [||] in
      let refresh_knn () =
        for i = 0 to n - 1 do
          if active.(i) then begin
            let acc = ref [] in
            for j = 0 to n - 1 do
              if j <> i && active.(j) then
                acc := (dist_of i j, j) :: !acc
            done;
            let arr = Array.of_list !acc in
            Array.sort (fun (a, _) (b, _) -> compare a b) arr;
            let m = min k_nn (Array.length arr) in
            nbrs.(i) <- Array.init m (fun p -> snd arr.(p));
            dists.(i) <- Array.init m (fun p -> fst arr.(p))
          end else begin
            nbrs.(i) <- [||];
            dists.(i) <- [||]
          end
        done in
      let row_sums () =
        let s = Float.Array.create n in
        Float.Array.fill s 0 n 0.;
        let n_active_minus_1 = float_of_int (!n_active - 1) in
        (match row_sum with
         | RowSum.Full ->
           for i = 0 to n - 1 do
             if active.(i) then begin
               let acc = ref 0. in
               for j = 0 to n - 1 do
                 if j <> i && active.(j) then
                   acc := !acc +. dist_of i j
               done;
               Float.Array.unsafe_set s i !acc
             end
           done
         | RowSum.Knn | RowSum.Topk ->
           for i = 0 to n - 1 do
             if active.(i) then begin
               let arr_d = dists.(i) in
               let kk = Array.length arr_d in
               if kk > 0 then begin
                 let acc = ref 0. in
                 for p = 0 to kk - 1 do
                   acc := !acc +. arr_d.(p)
                 done;
                 Float.Array.unsafe_set s i
                   (!acc *. n_active_minus_1 /. float_of_int kk)
               end
             end
           done);
        s in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (dense): index=%s, K=%d, rowsum=%s, sym=%s.\n%!" prefix
          (Interfaiss.Type.to_string index_type) k_nn
          (RowSum.to_string row_sum) (Symmetry.to_string symmetry);
      if verbose then
        Printf.eprintf "%s Bootstrapping K-NN graph via FAISS over %d points...\n%!"
          prefix n;
      bootstrap_cache ();
      if verbose then
        Printf.eprintf "%s Starting %d merges with Saitou-Nei K-NN refresh per iteration.\n%!"
          prefix (n - 3);
      while !n_active > 3 do
        refresh_knn ();
        let s = row_sums () in
        let cand = candidate_pairs symmetry nbrs active in
        let cand =
          if cand = [] then begin
            let acc = ref [] in
            for i = 0 to n - 1 do
              if active.(i) then
                for j = i + 1 to n - 1 do
                  if active.(j) then
                    List.accum acc (i, j)
                done
            done;
            !acc
          end else
            cand in
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        let best_q = ref infinity and best_pair = ref (-1, -1) in
        List.iter
          (fun (i, j) ->
            let q =
              n_act_minus_2 *. dist_of i j
              -. Float.Array.unsafe_get s i -. Float.Array.unsafe_get s j in
            if q < !best_q then begin
              best_q := q;
              best_pair := (i, j)
            end)
          cand;
        let i, j = !best_pair in
        let d_ij = dist_of i j in
        let s_i = Float.Array.unsafe_get s i and s_j = Float.Array.unsafe_get s j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i = 0.5 *. d_ij +. (s_i -. s_j) /. (2. *. n_act_minus_2_max) in
        let b_j = d_ij -. b_i in
        let b_i = max 0. b_i and b_j = max 0. b_j in
        let merged =
          Trees.Newick.join
            [| Trees.Newick.edge ~length:b_i (), trees.(i);
               Trees.Newick.edge ~length:b_j (), trees.(j) |] in
        trees.(i) <- merged;
        for x = 0 to n - 1 do
          if active.(x) && x <> i && x <> j then begin
            let d_ix = dist_of i x and d_jx = dist_of j x in
            let new_d = max 0. (0.5 *. (d_ix +. d_jx -. d_ij)) in
            Hashtbl.replace dist_cache (cache_key i x) new_d;
            Hashtbl.remove dist_cache (cache_key j x)
          end
        done;
        Hashtbl.remove dist_cache (cache_key i j);
        let si = size.(i) and sj = size.(j) in
        let total = si + sj in
        let si_f = float_of_int si and sj_f = float_of_int sj
        and total_f = float_of_int total in
        let new_emb = Float.Array.init dim (fun k ->
          (si_f *. Float.Array.unsafe_get embedding.(i) k
           +. sj_f *. Float.Array.unsafe_get embedding.(j) k)
          /. total_f) in
        embedding.(i) <- new_emb;
        size.(i) <- total;
        active.(j) <- false;
        decr n_active
      done;
      let active_three =
        let acc = ref [] in
        for i = n - 1 downto 0 do
          if active.(i) then List.accum acc i
        done;
        Array.of_rlist !acc in
      assert (Array.length active_three = 3);
      let i1, i2, i3 = active_three.(0), active_three.(1), active_three.(2) in
      let d12 = dist_of i1 i2
      and d13 = dist_of i1 i3
      and d23 = dist_of i2 i3 in
      let b1 = max 0. (0.5 *. (d12 +. d13 -. d23))
      and b2 = max 0. (0.5 *. (d12 +. d23 -. d13))
      and b3 = max 0. (0.5 *. (d13 +. d23 -. d12)) in
      let root =
        Trees.Newick.join
          [| Trees.Newick.edge ~length:b1 (), trees.(i1);
             Trees.Newick.edge ~length:b2 (), trees.(i2);
             Trees.Newick.edge ~length:b3 (), trees.(i3) |] in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (dense): built unrooted tree on %d leaves.\n%!"
          prefix n;
      root
    (* === Subquadratic experimental implementation ===

       State: each cluster (leaf or merged) gets its own slot in
       [0, 2n - 3); leaves occupy [0, n), merges fill [n, 2n - 3).
       Slots are NOT reused -- so the merge history is uniquely
       indexable and the Saitou-Nei recursion is well-defined.

       Saitou-Nei distance lookup: cache-hit on first try; otherwise
       recurse via the merge formula.  Cached entries are pruned
       lazily (we don't actively evict).

       K-NN maintenance: each active cluster has [nbrs.(v) : (slot, dist)
       array].  A reverse index [rev_nbrs.(v)] tracks who has v in
       their K-NN list, so on merge we update only the affected
       clusters in O(K^2) worst case.  When a K-NN list shrinks
       below K (because both i and j were neighbours and got replaced
       by a single u entry), we re-fill from a FAISS centroid query
       over the active set.

       FAISS: rebuilt from scratch each merge for the prototype
       (O(n_active * dim) per merge).  A periodic-rebuild or
       incremental insert/remove variant is the obvious next
       optimisation but is not load-bearing for the algorithmic
       complexity claim (which is the work-per-merge in the K-NN
       maintenance layer). *)
    let compute_subquadratic ?(verbose = false)
        ?(index_type = Interfaiss.Type.of_string "hnsw(32)")
        ?(k_nn = 10) ?(k_query_factor = 3)
        ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n = Array.length names in
      if Array.length data <> n then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sample count mismatch: %d names vs %d data rows"
             n (Array.length data));
      if n < 3 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sparse-NJ requires at least 3 leaves (got %d)" n);
      let dim = Float.Array.length data.(0) in
      let max_slots = 2 * n - 3 in
      (* Per-slot state *)
      let trees = Array.make max_slots (Trees.Newick.leaf "") in
      let active = Array.make max_slots false in
      let size = Array.make max_slots 0 in
      let embedding = Array.make max_slots (Float.Array.make 0 0.) in
      (* Merge history (only meaningful for slots >= n) *)
      let merge_left = Array.make max_slots (-1) in
      let merge_right = Array.make max_slots (-1) in
      let merge_dist = Array.make max_slots 0. in
      (* K-NN per cluster: sorted ascending by distance, length <= k_nn.
         Stored as (slot, dist) pairs.  Empty for inactive slots. *)
      let nbrs : (int * float) array array = Array.make max_slots [||] in
      (* Reverse index: rev_nbrs.(v) lists slots whose nbrs contains v.
         Maintained on every nbrs update so that on merge of (i, j) we
         know exactly which other clusters need patching. *)
      let rev_nbrs : (int, unit) Hashtbl.t array =
        Array.init max_slots (fun _ -> Hashtbl.create 8) in
      let rev_add v u =
        if v >= 0 && v < max_slots then Hashtbl.replace rev_nbrs.(v) u () in
      let rev_remove v u =
        if v >= 0 && v < max_slots then Hashtbl.remove rev_nbrs.(v) u in
      (* Saitou-Nei distance cache.  Populated lazily; entries for pairs
         involving inactive slots are not actively evicted. *)
      let dist_cache : (int * int, float) Hashtbl.t = Hashtbl.create (n * k_nn * 8) in
      let rec dist_of i j =
        if i = j then 0.
        else
          let key = cache_key i j in
          match Hashtbl.find_opt dist_cache key with
          | Some d -> d
          | None ->
            let d =
              if i < n && j < n then
                eucl_dist data.(i) data.(j)
              else begin
                (* Expand the younger (higher-index) of the two via its
                   merge.  Younger is always a merged slot (index >= n)
                   since min(i, j) might be a leaf but max must be >= n
                   when we reach this branch. *)
                let younger = max i j and older = min i j in
                let yl = merge_left.(younger)
                and yr = merge_right.(younger)
                and dlr = merge_dist.(younger) in
                0.5 *. (dist_of older yl +. dist_of older yr -. dlr)
              end in
            Hashtbl.add dist_cache key d;
            d in
      (* Pack the current active centroid embeddings (excluding [exclude]
         if [exclude >= 0]) into a (n_act, dim) Bigarray suitable for a
         FAISS index.  Returns the Bigarray + a slot-index map. *)
      let pack_active_embeddings ~exclude =
        let slots = Array.make max_slots (-1) in
        let n_act = ref 0 in
        for s = 0 to max_slots - 1 do
          if active.(s) && s <> exclude then begin
            slots.(!n_act) <- s;
            incr n_act
          end
        done;
        let m = !n_act in
        let ba = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout (max m 1) dim in
        for p = 0 to m - 1 do
          let s = slots.(p) in
          let row = embedding.(s) in
          for k = 0 to dim - 1 do
            ba.{p, k} <- Float.Array.unsafe_get row k
          done
        done;
        ba, slots, m in
      (* FAISS K_QUERY-NN of v's current centroid embedding (excluding
         v itself).  Returns a list of slot indices. *)
      let k_query = k_nn * k_query_factor in
      let faiss_expand v =
        let ba_active, slots, m = pack_active_embeddings ~exclude:v in
        if m = 0 then []
        else begin
          let query = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout 1 dim in
          let row = embedding.(v) in
          for k = 0 to dim - 1 do
            query.{0, k} <- Float.Array.unsafe_get row k
          done;
          let index = Interfaiss.create ~index_type dim in
          Interfaiss.train index ba_active;
          Interfaiss.add index ba_active;
          let k_use = min m k_query in
          let offsets, _ = Interfaiss.query index query k_use in
          Interfaiss.delete index;
          let res = ref [] in
          for k = 0 to k_use - 1 do
            let p = Int64.to_int offsets.{0, k} in
            if p >= 0 && p < m then res := slots.(p) :: !res
          done;
          !res
        end in
      (* Recompute nbrs.(v) from the given candidate list, ranking by
         current Saitou-Nei distance and keeping the top K.  Also
         updates the rev_nbrs index for both the removed-old and
         added-new neighbours. *)
      let rebuild_nbrs v candidates =
        let seen = Hashtbl.create 16 in
        let pairs = ref [] in
        List.iter (fun c ->
          if c <> v && c >= 0 && c < max_slots && active.(c)
             && not (Hashtbl.mem seen c) then begin
            Hashtbl.add seen c ();
            pairs := (c, dist_of v c) :: !pairs
          end) candidates;
        let arr = Array.of_list !pairs in
        Array.sort (fun (_, a) (_, b) -> compare a b) arr;
        let m = min k_nn (Array.length arr) in
        let new_nbrs = Array.sub arr 0 m in
        Array.iter (fun (j, _) -> rev_remove j v) nbrs.(v);
        nbrs.(v) <- new_nbrs;
        Array.iter (fun (j, _) -> rev_add j v) new_nbrs in
      (* Initialise leaves *)
      for i = 0 to n - 1 do
        trees.(i) <- Trees.Newick.leaf names.(i);
        active.(i) <- true;
        size.(i) <- 1;
        embedding.(i) <- Float.Array.copy data.(i)
      done;
      let n_active = ref n in
      let next_slot = ref n in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (subquadratic): index=%s, K=%d, K_QUERY=%d, rowsum=%s, sym=%s.\n%!"
          prefix (Interfaiss.Type.to_string index_type) k_nn k_query
          (RowSum.to_string row_sum) (Symmetry.to_string symmetry);
      (* Bootstrap: build initial K-NN lists for the leaves via a single
         FAISS pass.  The candidate set is K_QUERY-nearest by centroid
         (= Euclidean at bootstrap, since no merges have happened yet);
         distances are exact. *)
      if verbose then
        Printf.eprintf "%s Bootstrapping K-NN graph via FAISS over %d points...\n%!"
          prefix n;
      begin
        let ba = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout n dim in
        for i = 0 to n - 1 do
          let row = embedding.(i) in
          for k = 0 to dim - 1 do
            ba.{i, k} <- Float.Array.unsafe_get row k
          done
        done;
        let index = Interfaiss.create ~index_type dim in
        Interfaiss.train index ba;
        Interfaiss.add index ba;
        let kq = min n (k_query + 1) in
        let offsets, _ = Interfaiss.query index ba kq in
        Interfaiss.delete index;
        for i = 0 to n - 1 do
          let cands = ref [] in
          for k = 0 to kq - 1 do
            let j = Int64.to_int offsets.{i, k} in
            if j >= 0 && j < n && j <> i then cands := j :: !cands
          done;
          rebuild_nbrs i !cands
        done
      end;
      (* Build a slot-indexed nbrs view for [candidate_pairs] reuse;
         candidate_pairs expects an [int array array] keyed by slot. *)
      let nbrs_idx () =
        Array.init max_slots (fun v ->
          let a = nbrs.(v) in
          Array.init (Array.length a) (fun k -> fst a.(k))) in
      let row_sums () =
        let s = Float.Array.create max_slots in
        Float.Array.fill s 0 max_slots 0.;
        let n_active_minus_1 = float_of_int (!n_active - 1) in
        (match row_sum with
         | RowSum.Full ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let acc = ref 0. in
               for j = 0 to max_slots - 1 do
                 if j <> i && active.(j) then
                   acc := !acc +. dist_of i j
               done;
               Float.Array.unsafe_set s i !acc
             end
           done
         | RowSum.Knn | RowSum.Topk ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let arr = nbrs.(i) in
               let kk = Array.length arr in
               if kk > 0 then begin
                 let acc = ref 0. in
                 for p = 0 to kk - 1 do
                   acc := !acc +. snd arr.(p)
                 done;
                 Float.Array.unsafe_set s i
                   (!acc *. n_active_minus_1 /. float_of_int kk)
               end
             end
           done);
        s in
      (* Main loop *)
      while !n_active > 3 do
        let s = row_sums () in
        let cand = candidate_pairs symmetry (nbrs_idx ()) active in
        let cand =
          if cand = [] then begin
            let acc = ref [] in
            for i = 0 to max_slots - 1 do
              if active.(i) then
                for j = i + 1 to max_slots - 1 do
                  if active.(j) then
                    List.accum acc (i, j)
                done
            done;
            !acc
          end else
            cand in
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        let best_q = ref infinity and best_pair = ref (-1, -1) in
        List.iter
          (fun (i, j) ->
            let q =
              n_act_minus_2 *. dist_of i j
              -. Float.Array.unsafe_get s i -. Float.Array.unsafe_get s j in
            if q < !best_q then begin
              best_q := q;
              best_pair := (i, j)
            end)
          cand;
        let i, j = !best_pair in
        let d_ij = dist_of i j in
        let s_i = Float.Array.unsafe_get s i and s_j = Float.Array.unsafe_get s j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i = 0.5 *. d_ij +. (s_i -. s_j) /. (2. *. n_act_minus_2_max) in
        let b_j = d_ij -. b_i in
        let b_i = max 0. b_i and b_j = max 0. b_j in
        let u = !next_slot in
        incr next_slot;
        trees.(u) <- Trees.Newick.join
          [| Trees.Newick.edge ~length:b_i (), trees.(i);
             Trees.Newick.edge ~length:b_j (), trees.(j) |];
        active.(u) <- true;
        size.(u) <- size.(i) + size.(j);
        merge_left.(u) <- i;
        merge_right.(u) <- j;
        merge_dist.(u) <- d_ij;
        let si_f = float_of_int size.(i) and sj_f = float_of_int size.(j)
        and tot_f = float_of_int size.(u) in
        embedding.(u) <- Float.Array.init dim (fun k ->
          (si_f *. Float.Array.unsafe_get embedding.(i) k
           +. sj_f *. Float.Array.unsafe_get embedding.(j) k)
          /. tot_f);
        (* Step A: Build K-NN(u) from parent inheritance + FAISS expansion.
           Explicitly exclude i and j -- both are still flagged active at
           this point (we deactivate below to keep their nbrs / rev_nbrs
           snapshots readable) but they must not enter u's K-NN. *)
        let faiss_cands = faiss_expand u in
        let cands_u = ref [] in
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(i);
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(j);
        List.iter (fun x ->
          if x <> i && x <> j then cands_u := x :: !cands_u) faiss_cands;
        rebuild_nbrs u !cands_u;
        (* Snapshot the set of clusters whose K-NN list could be invalidated
           by this merge.  Two sources:
             - rev_nbrs.(i) and rev_nbrs.(j): clusters that had i or j in
               their K-NN list (they need i / j removed and replaced by u
               with a Saitou-Nei distance).
             - FAISS-K-NN(u): clusters that may now want u as a new K-NN
               entry because their centroid is near u's centroid (the
               "reverse insertion" candidate set; without this step,
               unaffected v's never get to consider u and we plateau well
               below the dense quality). *)
        let affected = Hashtbl.create 16 in
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(i);
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(j);
        List.iter (fun v ->
          if v <> i && v <> j && v <> u && active.(v) then
            Hashtbl.replace affected v ()) faiss_cands;
        (* Step B: Deactivate i and j (before patching others, so the
           patches won't reintroduce them).  Their nbrs are cleared. *)
        Array.iter (fun (x, _) -> rev_remove x i) nbrs.(i);
        Array.iter (fun (x, _) -> rev_remove x j) nbrs.(j);
        nbrs.(i) <- [||];
        nbrs.(j) <- [||];
        active.(i) <- false;
        active.(j) <- false;
        decr n_active;
        (* Step C: Patch each affected v.  v's old K-NN contained i
           and/or j; replace with u (if not already there) and rebuild.
           If v's neighbour list would drop below K, expand via FAISS. *)
        Hashtbl.iter (fun v () ->
          let cands_v = ref [] in
          Array.iter (fun (x, _) ->
            if active.(x) && x <> v then cands_v := x :: !cands_v) nbrs.(v);
          if not (List.mem u !cands_v) then cands_v := u :: !cands_v;
          (* If we expect to drop below K, expand from FAISS to refill *)
          let already = List.length !cands_v in
          if already < k_nn + 1 then begin
            let fx = faiss_expand v in
            List.iter (fun x ->
              if x <> v && not (List.mem x !cands_v) then
                cands_v := x :: !cands_v) fx
          end;
          rebuild_nbrs v !cands_v) affected
      done;
      let active_three =
        let acc = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then List.accum acc s
        done;
        Array.of_rlist !acc in
      assert (Array.length active_three = 3);
      let i1, i2, i3 = active_three.(0), active_three.(1), active_three.(2) in
      let d12 = dist_of i1 i2
      and d13 = dist_of i1 i3
      and d23 = dist_of i2 i3 in
      let b1 = max 0. (0.5 *. (d12 +. d13 -. d23))
      and b2 = max 0. (0.5 *. (d12 +. d23 -. d13))
      and b3 = max 0. (0.5 *. (d13 +. d23 -. d12)) in
      let root =
        Trees.Newick.join
          [| Trees.Newick.edge ~length:b1 (), trees.(i1);
             Trees.Newick.edge ~length:b2 (), trees.(i2);
             Trees.Newick.edge ~length:b3 (), trees.(i3) |] in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (subquadratic): built unrooted tree on %d leaves.\n%!"
          prefix n;
      root
    (* ====================================================================
       Periodic-rebuild FAISS implementation.

       Same algorithmic skeleton as [compute_subquadratic] but the
       per-merge FAISS rebuild is replaced by a persistent index plus
       tombstoning of merged-away slots and a side-scan over the
       merged-in slots that have appeared since the last rebuild.
       Every [rebuild_interval] merges (default ceil(sqrt n)) the
       index is rebuilt over the current active set and the tombstone /
       new-slot bookkeeping resets.

       With FAISS-Flat the per-query cost stays O(n_act * d) regardless,
       so this saves only the per-merge rebuild constant (the K-NN
       query cost dominates).  With HNSW the per-query cost would drop
       to O(K log n) and the algorithm reaches O(n^{1.5} log n)
       overall, but soft-deletion in FAISS HNSW is well known to
       degenerate queries toward linear as dead nodes accumulate in
       the graph; the periodic rebuild every sqrt(n) merges caps the
       tombstone fraction in the index at O(1/sqrt(n)), bounding the
       degeneration. *)
    let compute_periodic_rebuild ?(verbose = false)
        ?(index_type = Interfaiss.Type.of_string "flat")
        ?(k_nn = 10) ?(k_query_factor = 3)
        ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n = Array.length names in
      if Array.length data <> n then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sample count mismatch: %d names vs %d data rows"
             n (Array.length data));
      if n < 3 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sparse-NJ requires at least 3 leaves (got %d)" n);
      let dim = Float.Array.length data.(0) in
      let max_slots = 2 * n - 3 in
      let trees = Array.make max_slots (Trees.Newick.leaf "") in
      let active = Array.make max_slots false in
      let size = Array.make max_slots 0 in
      let embedding = Array.make max_slots (Float.Array.make 0 0.) in
      let merge_left = Array.make max_slots (-1) in
      let merge_right = Array.make max_slots (-1) in
      let merge_dist = Array.make max_slots 0. in
      let nbrs : (int * float) array array = Array.make max_slots [||] in
      let rev_nbrs : (int, unit) Hashtbl.t array =
        Array.init max_slots (fun _ -> Hashtbl.create 8) in
      let rev_add v u =
        if v >= 0 && v < max_slots then Hashtbl.replace rev_nbrs.(v) u () in
      let rev_remove v u =
        if v >= 0 && v < max_slots then Hashtbl.remove rev_nbrs.(v) u in
      let dist_cache : (int * int, float) Hashtbl.t = Hashtbl.create (n * k_nn * 8) in
      let rec dist_of i j =
        if i = j then 0.
        else
          let key = cache_key i j in
          match Hashtbl.find_opt dist_cache key with
          | Some d -> d
          | None ->
            let d =
              if i < n && j < n then
                eucl_dist data.(i) data.(j)
              else
                let younger = max i j and older = min i j in
                let yl = merge_left.(younger)
                and yr = merge_right.(younger)
                and dlr = merge_dist.(younger) in
                0.5 *. (dist_of older yl +. dist_of older yr -. dlr) in
            Hashtbl.add dist_cache key d;
            d in
      (* Persistent FAISS state.  [cur_index] / [cur_index_slots] hold
         the index built at the last rebuild and the slot-id for each
         row in it; [tombstoned] marks rows whose slot has since been
         merged away; [new_slots] are slot-ids that did not exist at
         the last rebuild (so they are not in the index and must be
         side-scanned during each query); [n_tombstones] counts
         tombstoned rows to size the over-fetch.  Rebuild every
         [rebuild_interval] merges (default ceil(sqrt n)). *)
      let cur_index: Interfaiss.t option ref = ref None in
      let cur_index_slots: int array ref = ref [||] in
      let tombstoned: bool array ref = ref [||] in
      let n_tombstones = ref 0 in
      let new_slots: int list ref = ref [] in
      let merges_since_rebuild = ref 0 in
      let rebuild_interval =
        max 1 (int_of_float (ceil (sqrt (float_of_int n)))) in
      let n_rebuilds = ref 0 in
      let rebuild_index () =
        (match !cur_index with
         | Some idx -> Interfaiss.delete idx
         | None -> ());
        let act = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then act := s :: !act
        done;
        let arr = Array.of_list !act in
        let m = Array.length arr in
        if m > 0 then begin
          let ba = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout m dim in
          for p = 0 to m - 1 do
            let row = embedding.(arr.(p)) in
            for k = 0 to dim - 1 do
              ba.{p, k} <- Float.Array.unsafe_get row k
            done
          done;
          let idx = Interfaiss.create ~index_type dim in
          Interfaiss.train idx ba;
          Interfaiss.add idx ba;
          cur_index := Some idx;
          cur_index_slots := arr;
          tombstoned := Array.make m false
        end else begin
          cur_index := None;
          cur_index_slots := [||];
          tombstoned := [||]
        end;
        n_tombstones := 0;
        new_slots := [];
        merges_since_rebuild := 0;
        incr n_rebuilds in
      let tombstone_slot s =
        let arr = !cur_index_slots in
        let tomb = !tombstoned in
        let m = Array.length arr in
        let p = ref 0 in
        while !p < m && (arr.(!p) <> s || tomb.(!p)) do
          incr p
        done;
        if !p < m then begin
          tomb.(!p) <- true;
          incr n_tombstones
        end in
      let k_query = k_nn * k_query_factor in
      let faiss_expand v =
        let res = ref [] in
        (match !cur_index with
         | None -> ()
         | Some idx ->
           let m = Array.length !cur_index_slots in
           if m > 0 then begin
             let query =
               Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout 1 dim in
             let row = embedding.(v) in
             for k = 0 to dim - 1 do
               query.{0, k} <- Float.Array.unsafe_get row k
             done;
             (* Over-fetch to account for tombstones and self-hit. *)
             let k_use = min m (k_query + !n_tombstones + 1) in
             let offsets, _ = Interfaiss.query idx query k_use in
             let tomb = !tombstoned in
             let slots = !cur_index_slots in
             for k = 0 to k_use - 1 do
               let p = Int64.to_int offsets.{0, k} in
               if p >= 0 && p < m && not tomb.(p) then begin
                 let s = slots.(p) in
                 if s <> v && active.(s) then res := s :: !res
               end
             done
           end);
        (* Side-scan the new slots not yet in the index. *)
        List.iter (fun s ->
          if s <> v && active.(s) then res := s :: !res) !new_slots;
        !res in
      let rebuild_nbrs v candidates =
        let seen = Hashtbl.create 16 in
        let pairs = ref [] in
        List.iter (fun c ->
          if c <> v && c >= 0 && c < max_slots && active.(c)
             && not (Hashtbl.mem seen c) then begin
            Hashtbl.add seen c ();
            pairs := (c, dist_of v c) :: !pairs
          end) candidates;
        let arr = Array.of_list !pairs in
        Array.sort (fun (_, a) (_, b) -> compare a b) arr;
        let m = min k_nn (Array.length arr) in
        let new_nbrs = Array.sub arr 0 m in
        Array.iter (fun (j, _) -> rev_remove j v) nbrs.(v);
        nbrs.(v) <- new_nbrs;
        Array.iter (fun (j, _) -> rev_add j v) new_nbrs in
      for i = 0 to n - 1 do
        trees.(i) <- Trees.Newick.leaf names.(i);
        active.(i) <- true;
        size.(i) <- 1;
        embedding.(i) <- Float.Array.copy data.(i)
      done;
      let n_active = ref n in
      let next_slot = ref n in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (periodic-rebuild): index=%s, K=%d, K_QUERY=%d, rebuild_interval=%d, rowsum=%s, sym=%s.\n%!"
          prefix (Interfaiss.Type.to_string index_type) k_nn k_query
          rebuild_interval (RowSum.to_string row_sum) (Symmetry.to_string symmetry);
      rebuild_index ();
      (* Bootstrap: build initial K-NN lists for the leaves via the
         persistent index. *)
      (match !cur_index with
       | None -> ()
       | Some idx ->
         let ba = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout n dim in
         for i = 0 to n - 1 do
           let row = embedding.(i) in
           for k = 0 to dim - 1 do
             ba.{i, k} <- Float.Array.unsafe_get row k
           done
         done;
         let kq = min n (k_query + 1) in
         let offsets, _ = Interfaiss.query idx ba kq in
         let slots = !cur_index_slots in
         let m = Array.length slots in
         for i = 0 to n - 1 do
           let cands = ref [] in
           for k = 0 to kq - 1 do
             let p = Int64.to_int offsets.{i, k} in
             if p >= 0 && p < m then begin
               let s = slots.(p) in
               if s <> i && active.(s) then cands := s :: !cands
             end
           done;
           rebuild_nbrs i !cands
         done);
      let nbrs_idx () =
        Array.init max_slots (fun v ->
          let a = nbrs.(v) in
          Array.init (Array.length a) (fun k -> fst a.(k))) in
      let row_sums () =
        let s = Float.Array.create max_slots in
        Float.Array.fill s 0 max_slots 0.;
        let n_active_minus_1 = float_of_int (!n_active - 1) in
        (match row_sum with
         | RowSum.Full ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let acc = ref 0. in
               for j = 0 to max_slots - 1 do
                 if j <> i && active.(j) then
                   acc := !acc +. dist_of i j
               done;
               Float.Array.unsafe_set s i !acc
             end
           done
         | RowSum.Knn | RowSum.Topk ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let arr = nbrs.(i) in
               let kk = Array.length arr in
               if kk > 0 then begin
                 let acc = ref 0. in
                 for p = 0 to kk - 1 do
                   acc := !acc +. snd arr.(p)
                 done;
                 Float.Array.unsafe_set s i
                   (!acc *. n_active_minus_1 /. float_of_int kk)
               end
             end
           done);
        s in
      while !n_active > 3 do
        let s = row_sums () in
        let cand = candidate_pairs symmetry (nbrs_idx ()) active in
        let cand =
          if cand = [] then begin
            let acc = ref [] in
            for i = 0 to max_slots - 1 do
              if active.(i) then
                for j = i + 1 to max_slots - 1 do
                  if active.(j) then
                    List.accum acc (i, j)
                done
            done;
            !acc
          end else
            cand in
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        let best_q = ref infinity and best_pair = ref (-1, -1) in
        List.iter
          (fun (i, j) ->
            let q =
              n_act_minus_2 *. dist_of i j
              -. Float.Array.unsafe_get s i -. Float.Array.unsafe_get s j in
            if q < !best_q then begin
              best_q := q;
              best_pair := (i, j)
            end)
          cand;
        let i, j = !best_pair in
        let d_ij = dist_of i j in
        let s_i = Float.Array.unsafe_get s i and s_j = Float.Array.unsafe_get s j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i = 0.5 *. d_ij +. (s_i -. s_j) /. (2. *. n_act_minus_2_max) in
        let b_j = d_ij -. b_i in
        let b_i = max 0. b_i and b_j = max 0. b_j in
        let u = !next_slot in
        incr next_slot;
        trees.(u) <- Trees.Newick.join
          [| Trees.Newick.edge ~length:b_i (), trees.(i);
             Trees.Newick.edge ~length:b_j (), trees.(j) |];
        active.(u) <- true;
        size.(u) <- size.(i) + size.(j);
        merge_left.(u) <- i;
        merge_right.(u) <- j;
        merge_dist.(u) <- d_ij;
        let si_f = float_of_int size.(i) and sj_f = float_of_int size.(j)
        and tot_f = float_of_int size.(u) in
        embedding.(u) <- Float.Array.init dim (fun k ->
          (si_f *. Float.Array.unsafe_get embedding.(i) k
           +. sj_f *. Float.Array.unsafe_get embedding.(j) k)
          /. tot_f);
        let faiss_cands = faiss_expand u in
        let cands_u = ref [] in
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(i);
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(j);
        List.iter (fun x ->
          if x <> i && x <> j then cands_u := x :: !cands_u) faiss_cands;
        rebuild_nbrs u !cands_u;
        let affected = Hashtbl.create 16 in
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(i);
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(j);
        List.iter (fun v ->
          if v <> i && v <> j && v <> u && active.(v) then
            Hashtbl.replace affected v ()) faiss_cands;
        Array.iter (fun (x, _) -> rev_remove x i) nbrs.(i);
        Array.iter (fun (x, _) -> rev_remove x j) nbrs.(j);
        nbrs.(i) <- [||];
        nbrs.(j) <- [||];
        active.(i) <- false;
        active.(j) <- false;
        decr n_active;
        (* Update persistent-index bookkeeping for this merge: i and j
           tombstoned in the existing index; u added to new_slots
           (it will be in the index after the next rebuild). *)
        tombstone_slot i;
        tombstone_slot j;
        new_slots := u :: !new_slots;
        Hashtbl.iter (fun v () ->
          let cands_v = ref [] in
          Array.iter (fun (x, _) ->
            if active.(x) && x <> v then cands_v := x :: !cands_v) nbrs.(v);
          if not (List.mem u !cands_v) then cands_v := u :: !cands_v;
          let already = List.length !cands_v in
          if already < k_nn + 1 then begin
            let fx = faiss_expand v in
            List.iter (fun x ->
              if x <> v && not (List.mem x !cands_v) then
                cands_v := x :: !cands_v) fx
          end;
          rebuild_nbrs v !cands_v) affected;
        incr merges_since_rebuild;
        if !merges_since_rebuild >= rebuild_interval then
          rebuild_index ()
      done;
      (match !cur_index with
       | Some idx -> Interfaiss.delete idx; cur_index := None
       | None -> ());
      let active_three =
        let acc = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then List.accum acc s
        done;
        Array.of_rlist !acc in
      assert (Array.length active_three = 3);
      let i1, i2, i3 = active_three.(0), active_three.(1), active_three.(2) in
      let d12 = dist_of i1 i2
      and d13 = dist_of i1 i3
      and d23 = dist_of i2 i3 in
      let b1 = max 0. (0.5 *. (d12 +. d13 -. d23))
      and b2 = max 0. (0.5 *. (d12 +. d23 -. d13))
      and b3 = max 0. (0.5 *. (d13 +. d23 -. d12)) in
      let root =
        Trees.Newick.join
          [| Trees.Newick.edge ~length:b1 (), trees.(i1);
             Trees.Newick.edge ~length:b2 (), trees.(i2);
             Trees.Newick.edge ~length:b3 (), trees.(i3) |] in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (periodic-rebuild): built unrooted tree on %d leaves (%d rebuilds).\n%!"
          prefix n !n_rebuilds;
      root
    (* ====================================================================
       Random-projection-forest periodic-rebuild implementation.

       Same skeleton as [compute_periodic_rebuild] but the spatial
       index is a native OCaml RP-forest (lib/RPForest.ml) instead of
       FAISS.  Build cost O(M n log n) per rebuild; K-NN query cost
       O(M log n + M * leaf_size * d), sub-linear in n at fixed M.
       No FFI overhead, which is the regime where FAISS-HNSW failed
       to beat Flat in our scaling sweep.  Forest contains the
       active set at the last rebuild; merged-away slots are filtered
       by [active.(s)] at query time, and new merged slots are
       side-scanned.  *)
    let compute_rp_forest_rebuild ?(verbose = false)
        ?(k_nn = 10) ?(k_query_factor = 3)
        ?(n_trees = 25) ?(leaf_size = 16)
        ?(distance = Distance.SaitouNei) ?(hyp_scale = 1.0)
        ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n = Array.length names in
      if Array.length data <> n then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sample count mismatch: %d names vs %d data rows"
             n (Array.length data));
      if n < 3 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sparse-NJ requires at least 3 leaves (got %d)" n);
      let dim = Float.Array.length data.(0) in
      let max_slots = 2 * n - 3 in
      let trees = Array.make max_slots (Trees.Newick.leaf "") in
      let active = Array.make max_slots false in
      let size = Array.make max_slots 0 in
      let embedding = Array.make max_slots (Float.Array.make 0 0.) in
      let merge_left = Array.make max_slots (-1) in
      let merge_right = Array.make max_slots (-1) in
      let merge_dist = Array.make max_slots 0. in
      let use_hyp =
        distance = Distance.Hyperbolic
        || distance = Distance.HyperbolicFrechet in
      (* Hyperbolic positions on the upper hyperboloid (dim + 1 coords),
         maintained by geodesic placement; only populated under a
         hyperbolic distance model. *)
      let hyp_pos = Array.make max_slots (Float.Array.make 0 0.) in
      let nbrs : (int * float) array array = Array.make max_slots [||] in
      let rev_nbrs : (int, unit) Hashtbl.t array =
        Array.init max_slots (fun _ -> Hashtbl.create 8) in
      let rev_add v u =
        if v >= 0 && v < max_slots then Hashtbl.replace rev_nbrs.(v) u () in
      let rev_remove v u =
        if v >= 0 && v < max_slots then Hashtbl.remove rev_nbrs.(v) u in
      let dist_cache : (int * int, float) Hashtbl.t = Hashtbl.create (n * k_nn * 8) in
      let n_dist_of_internal = ref 0 in
      let n_unique_dist_pairs = ref 0 in
      (* [Centroid] proxy: O(d) Euclidean distance between the
         size-weighted cluster centroids in [embedding], no recursion
         and no memoisation (the embedding is the authoritative state).
         [Hyperbolic] proxy: O(d) hyperboloid distance between the
         geodesically-placed positions in [hyp_pos].  [SaitouNei]:
         exact NJ update by memoised recursion. *)
      let dist_of_centroid i j =
        if i = j then 0. else eucl_dist embedding.(i) embedding.(j) in
      let dist_of_hyp i j =
        if i = j then 0. else hyp_dist hyp_pos.(i) hyp_pos.(j) in
      let rec dist_of_saitou_nei i j =
        incr n_dist_of_internal;
        if i = j then 0.
        else
          let key = cache_key i j in
          match Hashtbl.find_opt dist_cache key with
          | Some d -> d
          | None ->
            incr n_unique_dist_pairs;
            let d =
              if i < n && j < n then
                eucl_dist data.(i) data.(j)
              else
                let younger = max i j and older = min i j in
                let yl = merge_left.(younger)
                and yr = merge_right.(younger)
                and dlr = merge_dist.(younger) in
                0.5 *. (dist_of_saitou_nei older yl
                        +. dist_of_saitou_nei older yr -. dlr) in
            Hashtbl.add dist_cache key d;
            d in
      (* Bulk distance used by retrieval / K-NN lists / row sums / heap
         ordering.  Hybrid uses the centroid proxy here (cheap). *)
      let dist_of =
        match distance with
        | Distance.SaitouNei -> dist_of_saitou_nei
        | Distance.Centroid -> dist_of_centroid
        | Distance.Hyperbolic | Distance.HyperbolicFrechet -> dist_of_hyp
        | Distance.Hybrid -> dist_of_centroid in
      (* Confirmation distance used only to rank the popped shortlist
         and to set the chosen merge's branch lengths / stored
         merge_dist.  Hybrid upgrades this to exact Saitou-Nei; the
         other models reuse their bulk distance. *)
      let use_exact_confirm = (distance = Distance.Hybrid) in
      let dist_confirm i j =
        if use_exact_confirm then dist_of_saitou_nei i j else dist_of i j in
      let module Pt = struct
        type t = int * Float.Array.t
        let equal (a, _) (b, _) = a = b
        let hash (a, _) = a
        let embedding (_, e) = e
      end in
      let module Forest = RPForest.Make (Pt) in
      let cur_forest: Forest.t option ref = ref None in
      let new_slots: int list ref = ref [] in
      let merges_since_rebuild = ref 0 in
      let rebuild_interval =
        max 1 (int_of_float (ceil (sqrt (float_of_int n)))) in
      let n_rebuilds = ref 0 in
      let n_rebuild_nbrs_calls = ref 0 in
      let n_forest_queries = ref 0 in
      let n_dist_of_calls = ref 0 in
      let n_dist_of_misses = ref 0 in
      let total_affected = ref 0 in
      let n_merges = ref 0 in
      (* Heap pop accounting: every pop_best invocation classifies each
         popped entry as dead (a.i or a.j inactive), stale (version
         mismatch -> re-push), or fresh (accepted into the budget). *)
      let n_pops_total = ref 0 in
      let n_pops_dead = ref 0 in
      let n_pops_stale = ref 0 in
      let n_pops_fresh = ref 0 in
      let max_heap_size = ref 0 in
      let t_dist_of = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.dist_of" in
      let t_active_scan = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.active_scan" in
      let t_forest_build = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.forest_build" in
      let t_forest_query = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.forest_query" in
      let t_rebuild_nbrs = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.rebuild_nbrs" in
      let t_pop_best = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.pop_best" in
      let t_push_pairs = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.push_pairs" in
      let t_step_c = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.step_c" in
      let t_bootstrap = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.bootstrap" in
      let t_main_loop = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.main_loop" in
      let t_total = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_rp_forest_rebuild.total" in
      BiOCamLib.Tools.Timer.start t_total;
      let dist_of_raw = dist_of in
      let dist_of i j =
        BiOCamLib.Tools.Timer.start t_dist_of;
        incr n_dist_of_calls;
        (* Cache-miss accounting is only meaningful for the SaitouNei
           model; the Centroid proxy keeps no cache (unique_pairs
           stays 0, which is the signal that coverage collapsed). *)
        if distance = Distance.SaitouNei && i <> j then begin
          let key = cache_key i j in
          if not (Hashtbl.mem dist_cache key) then incr n_dist_of_misses
        end;
        let d = dist_of_raw i j in
        BiOCamLib.Tools.Timer.stop t_dist_of;
        d in
      let rebuild_forest () =
        BiOCamLib.Tools.Timer.start t_active_scan;
        let act = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then act := (s, embedding.(s)) :: !act
        done;
        let arr = Array.of_list !act in
        BiOCamLib.Tools.Timer.stop t_active_scan;
        BiOCamLib.Tools.Timer.start t_forest_build;
        if Array.length arr > 0 then
          cur_forest :=
            Some
              (Forest.build ~n_trees ~leaf_size ~seed:(0xC0FFEE + !n_rebuilds)
                 arr)
        else
          cur_forest := None;
        BiOCamLib.Tools.Timer.stop t_forest_build;
        new_slots := [];
        merges_since_rebuild := 0;
        incr n_rebuilds in
      let k_query = k_nn * k_query_factor in
      let forest_expand v =
        BiOCamLib.Tools.Timer.start t_forest_query;
        incr n_forest_queries;
        let res = ref [] in
        (match !cur_forest with
         | None -> ()
         | Some forest ->
           let k_use = k_query + 2 * !merges_since_rebuild + 1 in
           let cands = Forest.knn forest (v, embedding.(v)) k_use in
           List.iter
             (fun (s, _) ->
               if s <> v && s >= 0 && s < max_slots && active.(s) then
                 res := s :: !res)
             cands);
        List.iter (fun s ->
          if s <> v && active.(s) then res := s :: !res) !new_slots;
        BiOCamLib.Tools.Timer.stop t_forest_query;
        !res in
      let rebuild_nbrs v candidates =
        BiOCamLib.Tools.Timer.start t_rebuild_nbrs;
        incr n_rebuild_nbrs_calls;
        let seen = Hashtbl.create 16 in
        let pairs = ref [] in
        List.iter (fun c ->
          if c <> v && c >= 0 && c < max_slots && active.(c)
             && not (Hashtbl.mem seen c) then begin
            Hashtbl.add seen c ();
            pairs := (c, dist_of v c) :: !pairs
          end) candidates;
        let arr = Array.of_list !pairs in
        Array.sort (fun (_, a) (_, b) -> compare a b) arr;
        let m = min k_nn (Array.length arr) in
        let new_nbrs = Array.sub arr 0 m in
        Array.iter (fun (j, _) -> rev_remove j v) nbrs.(v);
        nbrs.(v) <- new_nbrs;
        Array.iter (fun (j, _) -> rev_add j v) new_nbrs;
        BiOCamLib.Tools.Timer.stop t_rebuild_nbrs in
      for i = 0 to n - 1 do
        trees.(i) <- Trees.Newick.leaf names.(i);
        active.(i) <- true;
        size.(i) <- 1;
        embedding.(i) <- Float.Array.copy data.(i);
        if use_hyp then hyp_pos.(i) <- hyp_lift ~scale:hyp_scale data.(i)
      done;
      let n_active = ref n in
      let next_slot = ref n in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (rp-forest): n_trees=%d, leaf_size=%d, K=%d, K_QUERY=%d, rebuild_interval=%d, distance=%s, rowsum=%s, sym=%s.\n%!"
          prefix n_trees leaf_size k_nn k_query rebuild_interval
          (Distance.to_string distance)
          (RowSum.to_string row_sum) (Symmetry.to_string symmetry);
      (* Diagnostic: how faithfully does the leaf lift reproduce the
         base (= leaf-level Saitou-Nei = Euclidean) distances?  At the
         leaves, before any merge, Saitou-Nei equals the Euclidean CA
         distance, so the Pearson correlation between hyp_dist on the
         lifted leaves and eucl_dist on the raw leaves directly measures
         the lift's distance fidelity.  Sampled over random leaf pairs. *)
      if verbose && use_hyp then begin
        let n_pairs = min 4000 (n * (n - 1) / 2) in
        let seed = ref 0x9E3779B1 in
        let next_rand m =
          (* xorshift, deterministic; avoids Random-state nondeterminism *)
          let x = !seed in
          let x = x lxor (x lsl 13) in
          let x = x lxor (x lsr 17) in
          let x = x lxor (x lsl 5) in
          seed := x land 0x3FFFFFFF;
          !seed mod m in
        let sx = ref 0. and sy = ref 0. and sxx = ref 0. and syy = ref 0.
        and sxy = ref 0. and cnt = ref 0 in
        for _ = 1 to n_pairs do
          let i = next_rand n and j = next_rand n in
          if i <> j then begin
            let dh = hyp_dist hyp_pos.(i) hyp_pos.(j) in
            let de = eucl_dist data.(i) data.(j) in
            sx := !sx +. dh; sy := !sy +. de;
            sxx := !sxx +. dh *. dh; syy := !syy +. de *. de;
            sxy := !sxy +. dh *. de; incr cnt
          end
        done;
        let nf = float_of_int !cnt in
        let cov = !sxy /. nf -. (!sx /. nf) *. (!sy /. nf) in
        let vx = !sxx /. nf -. (!sx /. nf) ** 2.
        and vy = !syy /. nf -. (!sy /. nf) ** 2. in
        let pearson =
          if vx > 0. && vy > 0. then cov /. sqrt (vx *. vy) else nan in
        Printf.eprintf "%s   leaf-lift fidelity (hyp vs Euclidean, scale=%g): Pearson=%.4f over %d pairs\n%!"
          prefix hyp_scale pearson !cnt
      end;
      (* Per-cluster cached mean K-NN distance r̂(v).  Used to compute
         row sums on the fly (r(v) = (n_act - 1) * r̂(v)) and as part
         of the heap key Q' = d(i,j) - r̂(i) - r̂(j). *)
      let r_hat = Float.Array.create max_slots in
      Float.Array.fill r_hat 0 max_slots 0.;
      let update_r_hat v =
        let arr = nbrs.(v) in
        let kk = Array.length arr in
        if kk = 0 then Float.Array.unsafe_set r_hat v 0.
        else begin
          let acc = ref 0. in
          for p = 0 to kk - 1 do
            acc := !acc +. snd arr.(p)
          done;
          Float.Array.unsafe_set r_hat v
            (!acc /. float_of_int kk)
        end in
      (* Per-cluster version: bumped each time the cluster's nbrs list
         is rewritten.  Heap entries store the versions at insertion
         time and are treated as stale (and recomputed) on pop if
         versions no longer match. *)
      let version = Array.make max_slots 0 in
      let heap = QHeap.create ~cap:(8 * n * k_nn) () in
      let q_prime i j =
        dist_of i j
        -. Float.Array.unsafe_get r_hat i
        -. Float.Array.unsafe_get r_hat j in
      let push_pair_unordered a b =
        if a <> b && a >= 0 && b >= 0
           && a < max_slots && b < max_slots
           && active.(a) && active.(b) then begin
          let i, j = if a < b then a, b else b, a in
          let qp = q_prime i j in
          QHeap.push heap
            { QHeap.q_prime = qp;
              vi = version.(i); vj = version.(j); i; j }
        end in
      let push_pairs_of v =
        BiOCamLib.Tools.Timer.start t_push_pairs;
        Array.iter (fun (x, _) -> push_pair_unordered v x) nbrs.(v);
        BiOCamLib.Tools.Timer.stop t_push_pairs in
      BiOCamLib.Tools.Timer.start t_bootstrap;
      rebuild_forest ();
      (match !cur_forest with
       | None -> ()
       | Some forest ->
         for i = 0 to n - 1 do
           let cands = Forest.knn forest (i, embedding.(i)) (k_query + 1) in
           let acc = ref [] in
           List.iter
             (fun (s, _) ->
               if s <> i && s >= 0 && s < max_slots && active.(s) then
                 acc := s :: !acc)
             cands;
           rebuild_nbrs i !acc;
           update_r_hat i
         done);
      for i = 0 to n - 1 do
        push_pairs_of i
      done;
      BiOCamLib.Tools.Timer.stop t_bootstrap;
      (* Skip the RowSum CLI selector for this mode: K-NN mean is the
         only estimator the heap supports.  RowSum.Full / Topk would
         be possible extensions but are not the validated path. *)
      let _ = row_sum in
      let _ = symmetry in
      (* Pop top [budget] fresh entries by Q', then pick the exact-Q
         minimum among them.  Q' = d - r̂_i - r̂_j is monotone-ish with
         Q = (n-2)d - (n-1)(r̂_i + r̂_j) but the orderings diverge
         (especially late in the run when n_act shrinks), so we
         survey a small batch and decide by the true Q. *)
      let q_exact i j =
        let n_act_minus_1 = float_of_int (!n_active - 1) in
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        n_act_minus_2 *. dist_confirm i j
        -. n_act_minus_1
           *. (Float.Array.unsafe_get r_hat i
               +. Float.Array.unsafe_get r_hat j) in
      let pop_best () =
        BiOCamLib.Tools.Timer.start t_pop_best;
        let cur_size = QHeap.length heap in
        if cur_size > !max_heap_size then max_heap_size := cur_size;
        let budget = max 4 (4 * k_nn) in
        let fresh = ref [] in
        let n_fresh = ref 0 in
        let going = ref true in
        while !going do
          match QHeap.pop heap with
          | None -> going := false
          | Some e ->
            incr n_pops_total;
            if not active.(e.QHeap.i) || not active.(e.QHeap.j) then
              incr n_pops_dead
            else if version.(e.QHeap.i) <> e.QHeap.vi
                 || version.(e.QHeap.j) <> e.QHeap.vj then begin
              incr n_pops_stale;
              let qp = q_prime e.QHeap.i e.QHeap.j in
              QHeap.push heap
                { e with QHeap.q_prime = qp;
                         vi = version.(e.QHeap.i);
                         vj = version.(e.QHeap.j) }
            end else begin
              incr n_pops_fresh;
              fresh := e :: !fresh;
              incr n_fresh;
              if !n_fresh >= budget then going := false
            end
        done;
        let best_q = ref infinity and best_pair = ref (-1, -1) in
        List.iter
          (fun (e: QHeap.entry) ->
            let q = q_exact e.QHeap.i e.QHeap.j in
            if q < !best_q then begin
              best_q := q;
              best_pair := (e.QHeap.i, e.QHeap.j)
            end;
            QHeap.push heap e)
          !fresh;
        BiOCamLib.Tools.Timer.stop t_pop_best;
        !best_pair in
      BiOCamLib.Tools.Timer.start t_main_loop;
      while !n_active > 3 do
        let i, j = pop_best () in
        let i, j =
          if i = -1 then begin
            (* Heap drained.  Fallback: any two active clusters --
               should not happen at our K. *)
            let a = ref (-1) and b = ref (-1) in
            let s = ref 0 in
            while !s < max_slots && !a = -1 do
              if active.(!s) then a := !s;
              incr s
            done;
            while !s < max_slots && !b = -1 do
              if active.(!s) then b := !s;
              incr s
            done;
            !a, !b
          end else
            i, j in
        (* Under Hybrid this is the exact Saitou-Nei distance (a cache
           hit, since pop_best already confirmed this pair); it both
           sets the branch lengths and is stored in merge_dist so the
           recursion remains consistent for future confirmations. *)
        let d_ij = dist_confirm i j in
        let n_act_minus_1 = float_of_int (!n_active - 1) in
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        let s_i = n_act_minus_1 *. Float.Array.unsafe_get r_hat i in
        let s_j = n_act_minus_1 *. Float.Array.unsafe_get r_hat j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i = 0.5 *. d_ij +. (s_i -. s_j) /. (2. *. n_act_minus_2_max) in
        let b_j = d_ij -. b_i in
        let b_i = max 0. b_i and b_j = max 0. b_j in
        let u = !next_slot in
        incr next_slot;
        trees.(u) <- Trees.Newick.join
          [| Trees.Newick.edge ~length:b_i (), trees.(i);
             Trees.Newick.edge ~length:b_j (), trees.(j) |];
        active.(u) <- true;
        size.(u) <- size.(i) + size.(j);
        merge_left.(u) <- i;
        merge_right.(u) <- j;
        merge_dist.(u) <- d_ij;
        let si_f = float_of_int size.(i) and sj_f = float_of_int size.(j)
        and tot_f = float_of_int size.(u) in
        embedding.(u) <- Float.Array.init dim (fun k ->
          (si_f *. Float.Array.unsafe_get embedding.(i) k
           +. sj_f *. Float.Array.unsafe_get embedding.(j) k)
          /. tot_f);
        (* Hyperbolic Q-distance: place u on the geodesic between its
           parents.  [Hyperbolic] uses the NJ branch length b_i from i
           (faithful under additivity, but b_i is noisy and the
           curvature amplifies that noise); [HyperbolicFrechet] uses the
           size-weighted barycenter --- distance
           (size_j / (size_i + size_j)) * d(i,j) from i --- which is
           bounded in [0, d(i,j)] and noise-free.  The Euclidean
           centroid above is still kept for RP-forest retrieval; only
           the Q-formula / row-sum distances are hyperbolic. *)
        if use_hyp then begin
          let t =
            match distance with
            | Distance.HyperbolicFrechet -> (sj_f /. tot_f) *. d_ij
            | _ -> b_i in
          hyp_pos.(u) <- hyp_geodesic hyp_pos.(i) hyp_pos.(j) t
        end;
        let forest_cands = forest_expand u in
        let cands_u = ref [] in
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(i);
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(j);
        List.iter (fun x ->
          if x <> i && x <> j then cands_u := x :: !cands_u) forest_cands;
        rebuild_nbrs u !cands_u;
        update_r_hat u;
        let affected = Hashtbl.create 16 in
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(i);
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(j);
        List.iter (fun v ->
          if v <> i && v <> j && v <> u && active.(v) then
            Hashtbl.replace affected v ()) forest_cands;
        total_affected := !total_affected + Hashtbl.length affected;
        incr n_merges;
        Array.iter (fun (x, _) -> rev_remove x i) nbrs.(i);
        Array.iter (fun (x, _) -> rev_remove x j) nbrs.(j);
        nbrs.(i) <- [||];
        nbrs.(j) <- [||];
        active.(i) <- false;
        active.(j) <- false;
        decr n_active;
        new_slots := u :: !new_slots;
        BiOCamLib.Tools.Timer.start t_step_c;
        Hashtbl.iter (fun v () ->
          let cands_v = ref [] in
          Array.iter (fun (x, _) ->
            if active.(x) && x <> v then cands_v := x :: !cands_v) nbrs.(v);
          if not (List.mem u !cands_v) then cands_v := u :: !cands_v;
          let already = List.length !cands_v in
          if already < k_nn + 1 then begin
            let fx = forest_expand v in
            List.iter (fun x ->
              if x <> v && not (List.mem x !cands_v) then
                cands_v := x :: !cands_v) fx
          end;
          rebuild_nbrs v !cands_v;
          update_r_hat v;
          version.(v) <- version.(v) + 1) affected;
        BiOCamLib.Tools.Timer.stop t_step_c;
        push_pairs_of u;
        Hashtbl.iter (fun v () -> push_pairs_of v) affected;
        incr merges_since_rebuild;
        if !merges_since_rebuild >= rebuild_interval then
          rebuild_forest ()
      done;
      BiOCamLib.Tools.Timer.stop t_main_loop;
      cur_forest := None;
      let active_three =
        let acc = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then List.accum acc s
        done;
        Array.of_rlist !acc in
      assert (Array.length active_three = 3);
      let i1, i2, i3 = active_three.(0), active_three.(1), active_three.(2) in
      let d12 = dist_of i1 i2
      and d13 = dist_of i1 i3
      and d23 = dist_of i2 i3 in
      let b1 = max 0. (0.5 *. (d12 +. d13 -. d23))
      and b2 = max 0. (0.5 *. (d12 +. d23 -. d13))
      and b3 = max 0. (0.5 *. (d13 +. d23 -. d12)) in
      let root =
        Trees.Newick.join
          [| Trees.Newick.edge ~length:b1 (), trees.(i1);
             Trees.Newick.edge ~length:b2 (), trees.(i2);
             Trees.Newick.edge ~length:b3 (), trees.(i3) |] in
      BiOCamLib.Tools.Timer.stop t_total;
      if verbose then begin
        Printf.eprintf "%s Sparse-NJ (rp-forest): built unrooted tree on %d leaves (%d rebuilds).\n%!"
          prefix n !n_rebuilds;
        let mean_affected =
          if !n_merges = 0 then 0.
          else float_of_int !total_affected /. float_of_int !n_merges in
        let fmerges = max 1 !n_merges in
        Printf.eprintf
          "%s  counts: n=%d, merges=%d, rebuild_nbrs=%d, forest_queries=%d, dist_of_top=%d, dist_of_internal=%d, unique_pairs=%d, mean_affected/merge=%.2f\n%!"
          prefix n !n_merges !n_rebuild_nbrs_calls !n_forest_queries
          !n_dist_of_calls !n_dist_of_internal !n_unique_dist_pairs
          mean_affected;
        Printf.eprintf
          "%s  pop_best:  total=%d (per_merge=%.1f), dead=%d, stale=%d, fresh=%d, ratio_stale=%.2f, ratio_dead=%.2f\n%!"
          prefix !n_pops_total
          (float_of_int !n_pops_total /. float_of_int fmerges)
          !n_pops_dead !n_pops_stale !n_pops_fresh
          (float_of_int !n_pops_stale
           /. float_of_int (max 1 !n_pops_fresh))
          (float_of_int !n_pops_dead
           /. float_of_int (max 1 !n_pops_fresh));
        Printf.eprintf
          "%s  heap:      max_size=%d, final_size=%d\n%!"
          prefix !max_heap_size (QHeap.length heap);
        Printf.eprintf
          "%s  dist_of:   recursion_factor=%.2f (internal/top), cache_density=%.3f (unique_pairs/n^2 with n=%d)\n%!"
          prefix
          (float_of_int !n_dist_of_internal
           /. float_of_int (max 1 !n_dist_of_calls))
          (float_of_int !n_unique_dist_pairs
           /. float_of_int (n * n))
          n;
        let pfx_t = "SparseNJ.compute_rp_forest_rebuild." in
        let pn = String.length pfx_t in
        Printf.eprintf "%s --- Timer breakdown ---\n%!" prefix;
        List.iter (fun (name, secs) ->
          let nn = String.length name in
          if nn >= pn && String.sub name 0 pn = pfx_t then
            Printf.eprintf "%s   %-30s %8.3f s\n%!" prefix
              (String.sub name pn (nn - pn)) secs)
          (BiOCamLib.Tools.Timer.snapshot ())
      end;
      root
    (* ====================================================================
       Cover-tree-backed subquadratic implementation.

       Same algorithmic skeleton as [compute_subquadratic] (rev_nbrs +
       FAISS-K_QUERY-expansion with Saitou-Nei reranking) but with the
       per-merge FAISS rebuild replaced by a persistent cover tree
       over the active cluster centroid embeddings.  Cover tree
       supports true insert / remove / K-NN in O(c^{O(1)} log n)
       amortised under bounded doubling dimension, so the K-NN
       queries that drove the empirical O(n^2) wall-clock of the FAISS-
       Flat path drop to O(K log n).  No FAISS index is built at all
       in this mode; the cover tree is the sole spatial index. *)
    let compute_cover_tree ?(verbose = false)
        ?(k_nn = 10) ?(k_query_factor = 3)
        ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n = Array.length names in
      if Array.length data <> n then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sample count mismatch: %d names vs %d data rows"
             n (Array.length data));
      if n < 3 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sparse-NJ requires at least 3 leaves (got %d)" n);
      let dim = Float.Array.length data.(0) in
      let max_slots = 2 * n - 3 in
      let trees = Array.make max_slots (Trees.Newick.leaf "") in
      let active = Array.make max_slots false in
      let size = Array.make max_slots 0 in
      let embedding = Array.make max_slots (Float.Array.make 0 0.) in
      let merge_left = Array.make max_slots (-1) in
      let merge_right = Array.make max_slots (-1) in
      let merge_dist = Array.make max_slots 0. in
      let nbrs : (int * float) array array = Array.make max_slots [||] in
      let rev_nbrs : (int, unit) Hashtbl.t array =
        Array.init max_slots (fun _ -> Hashtbl.create 8) in
      let rev_add v u =
        if v >= 0 && v < max_slots then Hashtbl.replace rev_nbrs.(v) u () in
      let rev_remove v u =
        if v >= 0 && v < max_slots then Hashtbl.remove rev_nbrs.(v) u in
      let dist_cache : (int * int, float) Hashtbl.t = Hashtbl.create (n * k_nn * 8) in
      let rec dist_of i j =
        if i = j then 0.
        else
          let key = cache_key i j in
          match Hashtbl.find_opt dist_cache key with
          | Some d -> d
          | None ->
            let d =
              if i < n && j < n then
                eucl_dist data.(i) data.(j)
              else
                let younger = max i j and older = min i j in
                let yl = merge_left.(younger)
                and yr = merge_right.(younger)
                and dlr = merge_dist.(younger) in
                0.5 *. (dist_of older yl +. dist_of older yr -. dlr) in
            Hashtbl.add dist_cache key d;
            d in
      (* Per-operation timers.  All accumulate to global Tools.Timer
         buckets that compute_cover_tree dumps via BiOCamLib.Tools.Timer.snapshot
         at exit.  Timers are non-nesting: each guards exactly its own
         function's body; we don't try to attribute nested costs (e.g.,
         dist_of inside rebuild_nbrs charges to rebuild_nbrs, not
         dist_of, because the latter's timer is started AFTER the
         former).  BiOCamLib.Tools.Timer.start is a no-op when the timer is
         already running, so nested starts don't double-count -- they
         just lose the inner attribution. *)
      let t_dist_of = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.dist_of" in
      let t_rebuild_nbrs = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.rebuild_nbrs" in
      let t_ct_insert = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.ct_insert" in
      let t_ct_remove = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.ct_remove" in
      let t_ct_knn = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.ct_knn" in
      let t_row_sums = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.row_sums" in
      let t_candidate_pairs = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.candidate_pairs" in
      let t_q_search = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.q_search" in
      let t_step_c = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.step_c" in
      let t_bootstrap = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.bootstrap" in
      let t_main_loop = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.main_loop" in
      let t_total = BiOCamLib.Tools.Timer.of_string "SparseNJ.compute_cover_tree.total" in
      BiOCamLib.Tools.Timer.start t_total;
      (* Wrap dist_of so its compute time is attributed (cache hits are
         essentially free; we want to know how much the recursive
         Saitou-Nei chain costs). *)
      let dist_of_raw = dist_of in
      let dist_of i j =
        BiOCamLib.Tools.Timer.start t_dist_of;
        let d = dist_of_raw i j in
        BiOCamLib.Tools.Timer.stop t_dist_of;
        d in
      let rebuild_nbrs v candidates =
        BiOCamLib.Tools.Timer.start t_rebuild_nbrs;
        let seen = Hashtbl.create 16 in
        let pairs = ref [] in
        List.iter (fun c ->
          if c <> v && c >= 0 && c < max_slots && active.(c)
             && not (Hashtbl.mem seen c) then begin
            Hashtbl.add seen c ();
            pairs := (c, dist_of v c) :: !pairs
          end) candidates;
        let arr = Array.of_list !pairs in
        Array.sort (fun (_, a) (_, b) -> compare a b) arr;
        let m = min k_nn (Array.length arr) in
        let new_nbrs = Array.sub arr 0 m in
        Array.iter (fun (j, _) -> rev_remove j v) nbrs.(v);
        nbrs.(v) <- new_nbrs;
        Array.iter (fun (j, _) -> rev_add j v) new_nbrs;
        BiOCamLib.Tools.Timer.stop t_rebuild_nbrs in
      for i = 0 to n - 1 do
        trees.(i) <- Trees.Newick.leaf names.(i);
        active.(i) <- true;
        size.(i) <- 1;
        embedding.(i) <- Float.Array.copy data.(i)
      done;
      let n_active = ref n in
      let next_slot = ref n in
      let module M = struct
        type t = int * Float.Array.t
        let equal (i, _) (j, _) = i = j
        let dist (_, a) (_, b) =
          let dm = Float.Array.length a in
          let acc = ref 0. in
          for k = 0 to dm - 1 do
            let delta = Float.Array.unsafe_get a k -. Float.Array.unsafe_get b k in
            acc := !acc +. delta *. delta
          done;
          sqrt !acc
      end in
      let module CT = CoverTree.Make (M) in
      let ct = CT.create () in
      BiOCamLib.Tools.Timer.start t_ct_insert;
      for i = 0 to n - 1 do
        CT.insert ct (i, embedding.(i))
      done;
      BiOCamLib.Tools.Timer.stop t_ct_insert;
      let k_query = k_nn * k_query_factor in
      let ct_expand v =
        if !n_active <= 1 then []
        else begin
          let q_point = (v, embedding.(v)) in
          BiOCamLib.Tools.Timer.start t_ct_knn;
          let knn = CT.knn ct q_point (k_query + 1) in
          BiOCamLib.Tools.Timer.stop t_ct_knn;
          List.filter_map (fun (s, _) ->
            if s = v || not active.(s) then None else Some s
          ) knn
        end in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (cover-tree): K=%d, K_QUERY=%d, rowsum=%s, sym=%s.\n%!"
          prefix k_nn k_query
          (RowSum.to_string row_sum) (Symmetry.to_string symmetry);
      (* Bootstrap: K-NN per leaf via the cover tree *)
      BiOCamLib.Tools.Timer.start t_bootstrap;
      for i = 0 to n - 1 do
        let cands = ct_expand i in
        rebuild_nbrs i cands
      done;
      BiOCamLib.Tools.Timer.stop t_bootstrap;
      let nbrs_idx () =
        Array.init max_slots (fun v ->
          let a = nbrs.(v) in
          Array.init (Array.length a) (fun k -> fst a.(k))) in
      let row_sums () =
        BiOCamLib.Tools.Timer.start t_row_sums;
        let s = Float.Array.create max_slots in
        Float.Array.fill s 0 max_slots 0.;
        let n_active_minus_1 = float_of_int (!n_active - 1) in
        (match row_sum with
         | RowSum.Full ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let acc = ref 0. in
               for j = 0 to max_slots - 1 do
                 if j <> i && active.(j) then
                   acc := !acc +. dist_of i j
               done;
               Float.Array.unsafe_set s i !acc
             end
           done
         | RowSum.Knn | RowSum.Topk ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let arr = nbrs.(i) in
               let kk = Array.length arr in
               if kk > 0 then begin
                 let acc = ref 0. in
                 for p = 0 to kk - 1 do
                   acc := !acc +. snd arr.(p)
                 done;
                 Float.Array.unsafe_set s i
                   (!acc *. n_active_minus_1 /. float_of_int kk)
               end
             end
           done);
        BiOCamLib.Tools.Timer.stop t_row_sums;
        s in
      BiOCamLib.Tools.Timer.start t_main_loop;
      while !n_active > 3 do
        let s = row_sums () in
        BiOCamLib.Tools.Timer.start t_candidate_pairs;
        let cand = candidate_pairs symmetry (nbrs_idx ()) active in
        let cand =
          if cand = [] then begin
            let acc = ref [] in
            for i = 0 to max_slots - 1 do
              if active.(i) then
                for j = i + 1 to max_slots - 1 do
                  if active.(j) then List.accum acc (i, j)
                done
            done;
            !acc
          end else cand in
        BiOCamLib.Tools.Timer.stop t_candidate_pairs;
        BiOCamLib.Tools.Timer.start t_q_search;
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        let best_q = ref infinity and best_pair = ref (-1, -1) in
        List.iter
          (fun (i, j) ->
            let q =
              n_act_minus_2 *. dist_of i j
              -. Float.Array.unsafe_get s i -. Float.Array.unsafe_get s j in
            if q < !best_q then begin
              best_q := q;
              best_pair := (i, j)
            end)
          cand;
        BiOCamLib.Tools.Timer.stop t_q_search;
        let i, j = !best_pair in
        let d_ij = dist_of i j in
        let s_i = Float.Array.unsafe_get s i and s_j = Float.Array.unsafe_get s j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i = 0.5 *. d_ij +. (s_i -. s_j) /. (2. *. n_act_minus_2_max) in
        let b_j = d_ij -. b_i in
        let b_i = max 0. b_i and b_j = max 0. b_j in
        let u = !next_slot in
        incr next_slot;
        trees.(u) <- Trees.Newick.join
          [| Trees.Newick.edge ~length:b_i (), trees.(i);
             Trees.Newick.edge ~length:b_j (), trees.(j) |];
        active.(u) <- true;
        size.(u) <- size.(i) + size.(j);
        merge_left.(u) <- i;
        merge_right.(u) <- j;
        merge_dist.(u) <- d_ij;
        let si_f = float_of_int size.(i) and sj_f = float_of_int size.(j)
        and tot_f = float_of_int size.(u) in
        embedding.(u) <- Float.Array.init dim (fun k ->
          (si_f *. Float.Array.unsafe_get embedding.(i) k
           +. sj_f *. Float.Array.unsafe_get embedding.(j) k)
          /. tot_f);
        let ct_cands = ct_expand u in
        let cands_u = ref [] in
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(i);
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) nbrs.(j);
        List.iter (fun x ->
          if x <> i && x <> j then cands_u := x :: !cands_u) ct_cands;
        rebuild_nbrs u !cands_u;
        let affected = Hashtbl.create 16 in
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(i);
        Hashtbl.iter (fun v () ->
          if v <> i && v <> j && active.(v) then
            Hashtbl.replace affected v ()) rev_nbrs.(j);
        List.iter (fun v ->
          if v <> i && v <> j && v <> u && active.(v) then
            Hashtbl.replace affected v ()) ct_cands;
        Array.iter (fun (x, _) -> rev_remove x i) nbrs.(i);
        Array.iter (fun (x, _) -> rev_remove x j) nbrs.(j);
        nbrs.(i) <- [||];
        nbrs.(j) <- [||];
        active.(i) <- false;
        active.(j) <- false;
        decr n_active;
        (* Update the cover tree: add u, remove i and j.  The add of
           u must precede the affected-v patching loop (so Step C
           queries see u); the removes of i / j must precede it too
           (so dead slots don't keep coming back as candidates). *)
        BiOCamLib.Tools.Timer.start t_ct_insert;
        CT.insert ct (u, embedding.(u));
        BiOCamLib.Tools.Timer.stop t_ct_insert;
        BiOCamLib.Tools.Timer.start t_ct_remove;
        let _ = CT.remove ct (i, embedding.(i)) in
        let _ = CT.remove ct (j, embedding.(j)) in
        BiOCamLib.Tools.Timer.stop t_ct_remove;
        BiOCamLib.Tools.Timer.start t_step_c;
        Hashtbl.iter (fun v () ->
          let cands_v = ref [] in
          Array.iter (fun (x, _) ->
            if active.(x) && x <> v then cands_v := x :: !cands_v) nbrs.(v);
          if not (List.mem u !cands_v) then cands_v := u :: !cands_v;
          let already = List.length !cands_v in
          if already < k_nn + 1 then begin
            let fx = ct_expand v in
            List.iter (fun x ->
              if x <> v && not (List.mem x !cands_v) then
                cands_v := x :: !cands_v) fx
          end;
          rebuild_nbrs v !cands_v) affected;
        BiOCamLib.Tools.Timer.stop t_step_c
      done;
      BiOCamLib.Tools.Timer.stop t_main_loop;
      let active_three =
        let acc = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then List.accum acc s
        done;
        Array.of_rlist !acc in
      assert (Array.length active_three = 3);
      let i1, i2, i3 = active_three.(0), active_three.(1), active_three.(2) in
      let d12 = dist_of i1 i2
      and d13 = dist_of i1 i3
      and d23 = dist_of i2 i3 in
      let b1 = max 0. (0.5 *. (d12 +. d13 -. d23))
      and b2 = max 0. (0.5 *. (d12 +. d23 -. d13))
      and b3 = max 0. (0.5 *. (d13 +. d23 -. d12)) in
      let root =
        Trees.Newick.join
          [| Trees.Newick.edge ~length:b1 (), trees.(i1);
             Trees.Newick.edge ~length:b2 (), trees.(i2);
             Trees.Newick.edge ~length:b3 (), trees.(i3) |] in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (cover-tree): built unrooted tree on %d leaves.\n%!"
          prefix n;
      BiOCamLib.Tools.Timer.stop t_total;
      if verbose then begin
        Printf.eprintf "%s --- Timer breakdown (raw snapshot) ---\n%!" prefix;
        let snap = BiOCamLib.Tools.Timer.snapshot () in
        Printf.eprintf "%s   snapshot has %d entries\n%!" prefix (List.length snap);
        List.iter (fun (name, secs) ->
          Printf.eprintf "%s   %-50s %8.3f s\n%!" prefix name secs) snap
      end;
      root
    (* ====================================================================
       Hyperbolic-embedding implementation.  Hyperboloid positions are
       maintained by geodesic placement and used directly as the K-NN
       retrieval metric (brute-force scan).  The hyperboloid helpers
       [lorentz_inner], [hyp_dist], [hyp_lift], [hyp_geodesic] are
       shared module-level functions defined near the top. *)
    let compute_hyperbolic ?(verbose = false)
        ?(k_nn = 10) ?(k_query_factor = 1) ?(hyp_scale = 1.0)
        ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n = Array.length names in
      if Array.length data <> n then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sample count mismatch: %d names vs %d data rows"
             n (Array.length data));
      if n < 3 then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Sparse-NJ requires at least 3 leaves (got %d)" n);
      let max_slots = 2 * n - 3 in
      (* Per-slot state.  Slots [0, n) are leaves; [n, max_slots) are
         merged clusters created in order. *)
      let trees = Array.make max_slots (Trees.Newick.leaf "") in
      let active = Array.make max_slots false in
      let size = Array.make max_slots 0 in
      let hyp_pos = Array.make max_slots (Float.Array.make 0 0.) in
      let merge_left = Array.make max_slots (-1) in
      let merge_right = Array.make max_slots (-1) in
      let merge_dist = Array.make max_slots 0. in
      let nbrs : (int * float) array array = Array.make max_slots [||] in
      let dist_cache : (int * int, float) Hashtbl.t = Hashtbl.create (n * k_nn * 8) in
      let rec dist_of i j =
        if i = j then 0.
        else
          let key = cache_key i j in
          match Hashtbl.find_opt dist_cache key with
          | Some d -> d
          | None ->
            let d =
              if i < n && j < n then
                eucl_dist data.(i) data.(j)
              else
                let younger = max i j and older = min i j in
                let yl = merge_left.(younger)
                and yr = merge_right.(younger)
                and dlr = merge_dist.(younger) in
                0.5 *. (dist_of older yl +. dist_of older yr -. dlr) in
            Hashtbl.add dist_cache key d;
            d in
      (* Brute-force [k_nn * k_query_factor] nearest by hyperbolic
         distance over the active set, then rerank by Saitou-Nei
         (via [dist_of]) and keep the top [k_nn].  The k_query_factor
         tunes how many hyperbolic candidates we consider; factor = 1
         uses hyperbolic ranking directly (the cheapest variant), but
         empirically the hyperbolic Spearman against Saitou-Nei isn't
         quite tight enough at top-K for sparse-NJ to land on the
         best Q-pair, so a small constant factor (2-3) is needed to
         widen the candidate pool.  O(n_active * k_query_factor) per
         query; Phase 2 will replace this with a cover-tree query for
         O(k_nn * k_query_factor * log n). *)
      let k_query = k_nn * k_query_factor in
      let hyp_knn v =
        let p_v = hyp_pos.(v) in
        let acc = ref [] in
        for s = 0 to max_slots - 1 do
          if s <> v && active.(s) then
            acc := (hyp_dist p_v hyp_pos.(s), s) :: !acc
        done;
        let arr = Array.of_list !acc in
        Array.sort (fun (a, _) (b, _) -> compare a b) arr;
        let m = min k_query (Array.length arr) in
        if k_query_factor <= 1 then
          Array.sub arr 0 (min k_nn m)
        else begin
          (* Rerank the top k_query hyperbolic candidates by Saitou-Nei *)
          let cand =
            Array.init m (fun p ->
              let (_d_hyp, x) = arr.(p) in
              (dist_of v x, x)) in
          Array.sort (fun (a, _) (b, _) -> compare a b) cand;
          let kk = min k_nn (Array.length cand) in
          Array.sub cand 0 kk
        end in
      (* Lift the n leaf principal-coord vectors into the hyperboloid. *)
      for i = 0 to n - 1 do
        trees.(i) <- Trees.Newick.leaf names.(i);
        active.(i) <- true;
        size.(i) <- 1;
        hyp_pos.(i) <- hyp_lift ~scale:hyp_scale data.(i)
      done;
      let n_active = ref n in
      let next_slot = ref n in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (hyperbolic): K=%d, scale=%g, rowsum=%s, sym=%s.\n%!"
          prefix k_nn hyp_scale (RowSum.to_string row_sum)
          (Symmetry.to_string symmetry);
      (* Bootstrap K-NN lists for the leaves *)
      for i = 0 to n - 1 do
        let knn = hyp_knn i in
        nbrs.(i) <- Array.map (fun (d, j) ->
          (* Replace hyp distance with the true Saitou-Nei (= Euclidean
             at leaves) so dists.(i) is consistent with what the
             Q-formula expects.  Hyperbolic was only used for ranking. *)
          let _ = d in
          (j, dist_of i j)) knn
      done;
      let nbrs_idx () =
        Array.init max_slots (fun v ->
          let a = nbrs.(v) in
          Array.init (Array.length a) (fun k -> fst a.(k))) in
      let row_sums () =
        let s = Float.Array.create max_slots in
        Float.Array.fill s 0 max_slots 0.;
        let n_active_minus_1 = float_of_int (!n_active - 1) in
        (match row_sum with
         | RowSum.Full ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let acc = ref 0. in
               for j = 0 to max_slots - 1 do
                 if j <> i && active.(j) then
                   acc := !acc +. dist_of i j
               done;
               Float.Array.unsafe_set s i !acc
             end
           done
         | RowSum.Knn | RowSum.Topk ->
           for i = 0 to max_slots - 1 do
             if active.(i) then begin
               let arr = nbrs.(i) in
               let kk = Array.length arr in
               if kk > 0 then begin
                 let acc = ref 0. in
                 for p = 0 to kk - 1 do
                   acc := !acc +. snd arr.(p)
                 done;
                 Float.Array.unsafe_set s i
                   (!acc *. n_active_minus_1 /. float_of_int kk)
               end
             end
           done);
        s in
      while !n_active > 3 do
        let s = row_sums () in
        let cand = candidate_pairs symmetry (nbrs_idx ()) active in
        let cand =
          if cand = [] then begin
            let acc = ref [] in
            for i = 0 to max_slots - 1 do
              if active.(i) then
                for j = i + 1 to max_slots - 1 do
                  if active.(j) then
                    List.accum acc (i, j)
                done
            done;
            !acc
          end else
            cand in
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        let best_q = ref infinity and best_pair = ref (-1, -1) in
        List.iter
          (fun (i, j) ->
            let q =
              n_act_minus_2 *. dist_of i j
              -. Float.Array.unsafe_get s i -. Float.Array.unsafe_get s j in
            if q < !best_q then begin
              best_q := q;
              best_pair := (i, j)
            end)
          cand;
        let i, j = !best_pair in
        let d_ij = dist_of i j in
        let s_i = Float.Array.unsafe_get s i and s_j = Float.Array.unsafe_get s j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i = 0.5 *. d_ij +. (s_i -. s_j) /. (2. *. n_act_minus_2_max) in
        let b_j = d_ij -. b_i in
        let b_i = max 0. b_i and b_j = max 0. b_j in
        let u = !next_slot in
        incr next_slot;
        trees.(u) <- Trees.Newick.join
          [| Trees.Newick.edge ~length:b_i (), trees.(i);
             Trees.Newick.edge ~length:b_j (), trees.(j) |];
        active.(u) <- true;
        size.(u) <- size.(i) + size.(j);
        merge_left.(u) <- i;
        merge_right.(u) <- j;
        merge_dist.(u) <- d_ij;
        (* Hyperbolic geodesic placement: p_u at distance b_i from p_i
           toward p_j on the geodesic.  Under additivity this gives
           d_H(p_u, p_x) = d_NJ(u, x) for every other active x. *)
        hyp_pos.(u) <- hyp_geodesic hyp_pos.(i) hyp_pos.(j) b_i;
        (* Deactivate i and j BEFORE the K-NN recomputation so the
           brute-force scan filters them out automatically. *)
        nbrs.(i) <- [||];
        nbrs.(j) <- [||];
        active.(i) <- false;
        active.(j) <- false;
        decr n_active;
        (* Rebuild K-NN of u and of every active cluster whose K-NN
           might have changed.  For the brute-force Phase 1 we
           recompute the K-NN of every active cluster on every merge
           via [hyp_knn]: O(n_active^2) per iteration, O(n^3) total --
           same as the dense reference but tests the algorithmic
           idea cleanly before the cover-tree work.  Phase 2 swaps
           [hyp_knn] for a cover-tree query (O(log n)) and refreshes
           only the clusters whose K-NN actually changed. *)
        for v = 0 to max_slots - 1 do
          if active.(v) then begin
            let knn = hyp_knn v in
            nbrs.(v) <- Array.map (fun (_d_hyp, x) ->
              (x, dist_of v x)) knn
          end
        done
      done;
      let active_three =
        let acc = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then List.accum acc s
        done;
        Array.of_rlist !acc in
      assert (Array.length active_three = 3);
      let i1, i2, i3 = active_three.(0), active_three.(1), active_three.(2) in
      let d12 = dist_of i1 i2
      and d13 = dist_of i1 i3
      and d23 = dist_of i2 i3 in
      let b1 = max 0. (0.5 *. (d12 +. d13 -. d23))
      and b2 = max 0. (0.5 *. (d12 +. d23 -. d13))
      and b3 = max 0. (0.5 *. (d13 +. d23 -. d12)) in
      let root =
        Trees.Newick.join
          [| Trees.Newick.edge ~length:b1 (), trees.(i1);
             Trees.Newick.edge ~length:b2 (), trees.(i2);
             Trees.Newick.edge ~length:b3 (), trees.(i3) |] in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (hyperbolic): built unrooted tree on %d leaves.\n%!"
          prefix n;
      root
    (* Dispatcher *)
    let compute ?(verbose = false)
        ?(mode = Mode.Dense)
        ?(index_type = Interfaiss.Type.of_string "hnsw(32)")
        ?(k_nn = 10) ?(k_query_factor = 3) ?(hyp_scale = 1.0)
        ?(distance = Distance.SaitouNei)
        ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      match mode with
      | Mode.Dense ->
        compute_dense ~verbose ~index_type ~k_nn ~row_sum ~symmetry names data
      | Mode.Subquadratic ->
        compute_subquadratic ~verbose ~index_type ~k_nn ~k_query_factor
          ~row_sum ~symmetry names data
      | Mode.Hyperbolic ->
        compute_hyperbolic ~verbose ~k_nn ~k_query_factor ~hyp_scale
          ~row_sum ~symmetry names data
      | Mode.CoverTree ->
        compute_cover_tree ~verbose ~k_nn ~k_query_factor
          ~row_sum ~symmetry names data
      | Mode.PeriodicRebuild ->
        compute_periodic_rebuild ~verbose ~index_type ~k_nn ~k_query_factor
          ~row_sum ~symmetry names data
      | Mode.RPForestRebuild ->
        compute_rp_forest_rebuild ~verbose ~k_nn ~k_query_factor
          ~distance ~hyp_scale ~row_sum ~symmetry names data
  end: sig
    module RowSum:
      sig
        type t =
          | Knn
          | Topk
          | Full
        val of_string: string -> t
        val to_string: t -> string
      end
    module Symmetry:
      sig
        type t =
          | One
          | Both
        val of_string: string -> t
        val to_string: t -> string
      end
    module Distance:
      sig
        type t =
          | SaitouNei
          | Centroid
          | Hyperbolic
          | HyperbolicFrechet
          | Hybrid
        val of_string: string -> t
        val to_string: t -> string
      end
    module Mode:
      sig
        type t =
          | Dense
          | Subquadratic
          | Hyperbolic
          | CoverTree
          | PeriodicRebuild
          | RPForestRebuild
        val of_string: string -> t
        val to_string: t -> string
      end
    (* Build an unrooted Newick tree by sparse-NJ on the row-stored
       embedding [data].  [names] is the array of leaf names matching
       the rows of [data].  Each internal merge contributes a join node
       carrying the NJ-computed branch lengths above its children; the
       final three-way merge produces a trifurcating unrooted root.
       [mode] selects between the validated O(n^3)/O(n^2) dense
       reference, the experimental centroid-based subquadratic
       prototype, and the hyperbolic-embedding mode (Phase 1: with
       brute-force K-NN scan, same O(n^3) asymptotic as the dense
       reference; Phase 2 will swap in a cover-tree query for
       O(n K log n)).  [index_type] picks the FAISS index;
       [k_query_factor] sets the FAISS expansion factor in
       Subquadratic mode; [hyp_scale] sets the radial scale s for the
       initial hyperboloid lift in Hyperbolic mode (empirically
       s in [0.5, 1.0] works; default 1.0). *)
    val compute:
      ?verbose:bool ->
      ?mode:Mode.t ->
      ?index_type:Interfaiss.Type.t ->
      ?k_nn:int ->
      ?k_query_factor:int ->
      ?hyp_scale:float ->
      ?distance:Distance.t ->
      ?row_sum:RowSum.t ->
      ?symmetry:Symmetry.t ->
      string array ->
      Float.Array.t array ->
      Trees.Newick.t
  end
)
