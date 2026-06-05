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

    Three implementation modes are available via [Mode.t]:
      - [Dense] (default, validated): a persistent symmetric distance
        cache stores the current Saitou-Nei distance between every pair
        of active clusters, updated symmetrically at each merge.  K-NN
        selection at every iteration brute-scans active candidates and
        ranks them by their current Saitou-Nei distance.  Time
        O(n^3); memory O(n^2).  Reproduces the test prototype
        (test/Trees/sparse_nj.py) byte-identical on prot_k5_kt0.035.
      - [PeriodicRebuild]: K-NN-restricted scan-NJ with exact Saitou-Nei
        distances served from a persistent FAISS index rebuilt every
        sqrt(n) merges.  Theta(n^2) time and memory; the recommended
        exact engine at scale.
      - [RPForestRebuild]: K-NN retrieval from a native RP-forest with a
        heap-driven merge loop; under the [Hyperbolic] geometric-proxy
        Q-distance it drops the Theta(n^2) Saitou-Nei recursion for
        sub-quadratic memory, at a quality cost -- the fallback for very
        large n.

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
       [Full] uses the exact global row sum, recomputed lazily from the
       Saitou-Nei distance cache (with centroid-Euclidean fallback) at
       each iteration. *)
    module RowSum =
      struct
        type t =
          | Knn
          | Full
        let of_string = function
          | "knn" -> Knn
          | "full" -> Full
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ row-sum estimator" s
        let to_string = function
          | Knn -> "knn"
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
       KPopPhylo-Subquadratic.tex).  [Centroid] is an O(d) geometric
       proxy: the Euclidean distance between the size-weighted cluster
       centroids already maintained for K-NN retrieval, with no
       recursion (and so no n^2 cache fill).  [Hyperbolic] is also an
       O(d) recursion-free proxy, but tracks Saitou-Nei far better than
       the centroid: leaf positions are lifted onto the upper
       hyperboloid and each merge places the new cluster on the
       geodesic between its parents at the NJ branch length, so
       hyperbolic distances follow the exact update closely across the
       whole run (Phase 0: Pearson 0.97+ vs Saitou-Nei, where the
       centroid collapses to ~0.33 mid-run).  The proxy models
       ([Centroid], [Hyperbolic]) are only honoured by the RP-forest
       mode; [SaitouNei] is the default everywhere. *)
    module Distance =
      struct
        type t =
          | SaitouNei
          | Centroid
          | Hyperbolic
        let of_string = function
          | "saitou-nei" | "saitou_nei" | "saitounei" | "nj" -> SaitouNei
          | "centroid" -> Centroid
          | "hyperbolic" -> Hyperbolic
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ distance model" s
        let to_string = function
          | SaitouNei -> "saitou-nei"
          | Centroid -> "centroid"
          | Hyperbolic -> "hyperbolic"
      end
    (* Implementation mode.  [Dense] is the validated O(n^3) reference
       (compute_dense); [PeriodicRebuild] is the recommended exact
       Theta(n^2) engine (persistent-FAISS scan-NJ,
       compute_periodic_rebuild); [RPForestRebuild] is the RP-forest
       heap-driven path whose [Hyperbolic] proxy distance gives the
       sub-quadratic-memory fallback (compute_rp_forest_rebuild). *)
    module Mode =
      struct
        type t =
          | Dense
          | PeriodicRebuild
          | RPForestRebuild
        let of_string = function
          | "dense" -> Dense
          | "periodic-rebuild" | "periodic_rebuild" | "periodicrebuild" ->
            PeriodicRebuild
          | "rp-forest" | "rp_forest" | "rpforest"
          | "rp-forest-rebuild" | "rp_forest_rebuild" ->
            RPForestRebuild
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ mode" s
        let to_string = function
          | Dense -> "dense"
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
    (* Hyperboloid model, used by the [Hyperbolic] distance proxy of
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
       0.97+; see KPopPhylo-Subquadratic.tex Path 3). *)
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
       stored versions against the per-cluster version array.  The
       backing store is BiOCamLib's growable [Tools.ArrayStack]; [push]
       and [pop] sift up / down over its random-access slots. *)
    module QHeap =
      struct
        type entry = {
          q_prime: float;
          vi: int;
          vj: int;
          i: int;
          j: int
        }
        module Store = BiOCamLib.Tools.ArrayStack
        type t = entry Store.t
        let create () = Store.empty ()
        let length = Store.length
        (* The sift arithmetic only ever indexes slots known to be in
           [0, length), so the random-access reads/writes use the
           no-bounds-check [.@!()] accessors. *)
        let ( .@!() ) = Store.unsafe_get
        let ( .@!()<- ) = Store.unsafe_set
        let push h e =
          Store.push h e;
          let i = ref (Store.length h - 1) in
          let going = ref true in
          while !going && !i > 0 do
            let p = (!i - 1) / 2 in
            if h.@!(!i).q_prime < h.@!(p).q_prime then begin
              let t = h.@!(!i) in
              h.@!(!i) <- h.@!(p);
              h.@!(p) <- t;
              i := p
            end else going := false
          done
        let pop h =
          if Store.is_empty h then None
          else begin
            let top = h.@!(0) in
            let last = Store.pop h in
            if not (Store.is_empty h) then begin
              h.@!(0) <- last;
              let n = Store.length h in
              let i = ref 0 in
              let going = ref true in
              while !going do
                let l = 2 * !i + 1 and r = 2 * !i + 2 in
                let smallest = ref !i in
                if l < n
                   && h.@!(l).q_prime < h.@!(!smallest).q_prime then
                  smallest := l;
                if r < n
                   && h.@!(r).q_prime < h.@!(!smallest).q_prime then
                  smallest := r;
                if !smallest = !i then going := false
                else begin
                  let t = h.@!(!i) in
                  h.@!(!i) <- h.@!(!smallest);
                  h.@!(!smallest) <- t;
                  i := !smallest
                end
              done
            end;
            Some top
          end
      end
    (* Validate the inputs shared by every compute_* entry point and
       return the leaf count; [func] is the caller's __FUNCTION__ so the
       raised errors name the actual mode. *)
    let check_inputs func names data =
      let n = Array.length names in
      if Array.length data <> n then
        Exception.raise func IO_Format
          (Printf.sprintf "Sample count mismatch: %d names vs %d data rows"
             n (Array.length data));
      if n < 3 then
        Exception.raise func IO_Format
          (Printf.sprintf "Sparse-NJ requires at least 3 leaves (got %d)" n);
      n
    (* Saitou-Nei split of an edge length [d_ij] into its two child
       branch lengths, clamped non-negative. *)
    let nj_branch_lengths ~d_ij ~s_i ~s_j ~n_act_minus_2_max =
      let b_i = 0.5 *. d_ij +. (s_i -. s_j) /. (2. *. n_act_minus_2_max) in
      let b_j = d_ij -. b_i in
      max 0. b_i, max 0. b_j
    (* Scan the candidate pairs for the minimum Q = (n_act - 2) d(i, j)
       - s_i - s_j, returning the argmin pair (first minimum wins). *)
    let best_q_pair ~n_act_minus_2 ~s ~dist_of cand =
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
      !best_pair
    (* The closing three-leaf join shared by every mode: collect the
       three still-active slots (descending), compute the trifurcating
       branch lengths, and build the unrooted root.  [bound] is the
       active array's length (n for dense, 2n - 3 for the slotted
       modes); [dist_of] is the mode's own distance closure. *)
    let finalize_root ~bound ~active ~trees ~dist_of =
      let active_three =
        let acc = ref [] in
        for s = bound - 1 downto 0 do
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
      Trees.Newick.join
        [| Trees.Newick.edge ~length:b1 (), trees.(i1);
           Trees.Newick.edge ~length:b2 (), trees.(i2);
           Trees.Newick.edge ~length:b3 (), trees.(i3) |]
    (* The memoised Saitou-Nei distance recursion shared by the slotted
       modes: a cache hit, else the leaf Euclidean distance or the NJ
       merge update expanded over the younger cluster's history (the
       younger is always a merged slot, index >= n). *)
    let make_dist_of ~n ~data ~merge_left ~merge_right ~merge_dist ~dist_cache =
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
      dist_of
    (* Deactivate the merged pair (i, j): drop them from the reverse
       index, clear their K-NN lists, mark them inactive and decrement
       the active count.  Shared by the slotted modes. *)
    let deactivate_pair ~i ~j ~nbrs ~rev_remove ~active ~n_active =
      Array.iter (fun (x, _) -> rev_remove x i) nbrs.(i);
      Array.iter (fun (x, _) -> rev_remove x j) nbrs.(j);
      nbrs.(i) <- [||];
      nbrs.(j) <- [||];
      active.(i) <- false;
      active.(j) <- false;
      decr n_active
    (* The active clusters whose K-NN list a merge of (i, j) into u may
       invalidate: everyone who had i or j as a neighbour, plus the
       index's K-NN of u ([idx_cands]) as reverse-insertion candidates. *)
    let collect_affected ~i ~j ~u ~active ~rev_nbrs idx_cands =
      let affected = Hashtbl.create 16 in
      Hashtbl.iter (fun v () ->
        if v <> i && v <> j && active.(v) then
          Hashtbl.replace affected v ()) rev_nbrs.(i);
      Hashtbl.iter (fun v () ->
        if v <> i && v <> j && active.(v) then
          Hashtbl.replace affected v ()) rev_nbrs.(j);
      List.iter (fun v ->
        if v <> i && v <> j && v <> u && active.(v) then
          Hashtbl.replace affected v ()) idx_cands;
      affected
    (* Size-weighted centroid of two cluster embeddings. *)
    let weighted_centroid ~dim ~si_f ~sj_f ~tot_f emb_i emb_j =
      Float.Array.init dim (fun k ->
        (si_f *. Float.Array.unsafe_get emb_i k
         +. sj_f *. Float.Array.unsafe_get emb_j k)
        /. tot_f)
    (* The shared per-cluster K-NN view over the [nbrs] slot arrays,
       projected to plain slot indices for [candidate_pairs]. *)
    let make_nbrs_idx ~max_slots ~nbrs =
      fun () ->
        Array.init max_slots (fun v ->
          let a = nbrs.(v) in
          Array.init (Array.length a) (fun k -> fst a.(k)))
    (* The shared K-NN rebuild for the slotted modes: rank the deduped
       candidate slots by [dist_of] from v, keep the top [k_nn], and
       keep the reverse index consistent. *)
    let make_rebuild_nbrs ~k_nn ~max_slots ~active ~nbrs ~rev_add ~rev_remove ~dist_of =
      fun v candidates ->
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
        Array.iter (fun (j, _) -> rev_add j v) new_nbrs
    (* The shared row-sum estimator for the slotted scan modes: [Full]
       is the exact active row sum, [Knn] the K-NN mean rescaled by
       n_act - 1.  [bound] is the active array's length. *)
    let make_row_sums ~bound ~active ~nbrs ~dist_of ~row_sum ~n_active =
      fun () ->
        let s = Float.Array.create bound in
        Float.Array.fill s 0 bound 0.;
        let n_active_minus_1 = float_of_int (!n_active - 1) in
        (match row_sum with
         | RowSum.Full ->
           for i = 0 to bound - 1 do
             if active.(i) then begin
               let acc = ref 0. in
               for j = 0 to bound - 1 do
                 if j <> i && active.(j) then
                   acc := !acc +. dist_of i j
               done;
               Float.Array.unsafe_set s i !acc
             end
           done
         | RowSum.Knn ->
           for i = 0 to bound - 1 do
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
        s
    (* The per-slot state shared by the slotted scan modes.  Each mode
       allocates these arrays/refs locally and bundles them here by field
       punning, so [run_core] can drive the merge loop over them. *)
    type slot_state = {
      trees: Trees.Newick.t array;
      active: bool array;
      size: int array;
      embedding: Float.Array.t array;
      merge_left: int array;
      merge_right: int array;
      merge_dist: float array;
      nbrs: (int * float) array array;
      rev_nbrs: (int, unit) Hashtbl.t array;
      n_active: int ref;
      next_slot: int ref
    }
    (* The spatial-index strategy that distinguishes the scan modes:
       [expand v] is the index K-NN of slot v, [on_merge] keeps the index
       consistent after a merge (i, j) -> u, [rebuild] is the periodic
       index rebuild (a no-op where the index is incremental). *)
    type strategy = {
      expand: int -> int list;
      on_merge: i:int -> j:int -> u:int -> unit;
      rebuild: unit -> unit
    }
    (* Create the merged slot u from (i, j): a join node carrying the two
       NJ branch lengths, its merge-history entry, and its size-weighted
       centroid embedding.  Returns u. *)
    let do_merge ~dim ~st ~i ~j ~d_ij ~b_i ~b_j =
      let u = !(st.next_slot) in
      incr st.next_slot;
      st.trees.(u) <- Trees.Newick.join
        [| Trees.Newick.edge ~length:b_i (), st.trees.(i);
           Trees.Newick.edge ~length:b_j (), st.trees.(j) |];
      st.active.(u) <- true;
      st.size.(u) <- st.size.(i) + st.size.(j);
      st.merge_left.(u) <- i;
      st.merge_right.(u) <- j;
      st.merge_dist.(u) <- d_ij;
      let si_f = float_of_int st.size.(i) and sj_f = float_of_int st.size.(j)
      and tot_f = float_of_int st.size.(u) in
      st.embedding.(u) <-
        weighted_centroid ~dim ~si_f ~sj_f ~tot_f st.embedding.(i) st.embedding.(j);
      u
    (* All (i, j), i < j, both active, in the canonical nested-loop order
       -- the fallback candidate set when the K-NN graph yields no pairs.
       [bound] is the active array's length. *)
    let all_active_pairs ~bound ~active =
      let acc = ref [] in
      for i = 0 to bound - 1 do
        if active.(i) then
          for j = i + 1 to bound - 1 do
            if active.(j) then List.accum acc (i, j)
          done
      done;
      !acc
    (* Step C: patch every affected cluster v whose K-NN list a merge into
       u may have invalidated -- keep its surviving neighbours, splice in
       u, refill from [expand] if it would drop below k_nn, and rebuild. *)
    let patch_affected ~u ~k_nn ~active ~nbrs ~rebuild_nbrs ~expand affected =
      Hashtbl.iter (fun v () ->
        let cands_v = ref [] in
        Array.iter (fun (x, _) ->
          if active.(x) && x <> v then cands_v := x :: !cands_v) nbrs.(v);
        if not (List.mem u !cands_v) then cands_v := u :: !cands_v;
        let already = List.length !cands_v in
        if already < k_nn + 1 then begin
          let fx = expand v in
          List.iter (fun x ->
            if x <> v && not (List.mem x !cands_v) then
              cands_v := x :: !cands_v) fx
        end;
        rebuild_nbrs v !cands_v) affected
    (* The shared scan-merge driver for the slotted modes.  While more
       than three clusters remain: estimate row sums, pick the minimum-Q
       candidate pair, merge it, rebuild u's K-NN from its parents plus the
       index, then patch every affected cluster -- the [strategy] supplies
       the index-specific expand / on_merge / rebuild.  Returns the
       trifurcating unrooted root. *)
    let run_core ~max_slots ~dim ~k_nn ~symmetry ~st ~dist_of ~rebuild_nbrs
        ~row_sums ~nbrs_idx ~rev_remove ~strategy =
      while !(st.n_active) > 3 do
        let s = row_sums () in
        let cand = candidate_pairs symmetry (nbrs_idx ()) st.active in
        let cand =
          if cand = [] then all_active_pairs ~bound:max_slots ~active:st.active
          else cand in
        let n_act_minus_2 = float_of_int (!(st.n_active) - 2) in
        let i, j = best_q_pair ~n_act_minus_2 ~s ~dist_of cand in
        let d_ij = dist_of i j in
        let s_i = Float.Array.unsafe_get s i and s_j = Float.Array.unsafe_get s j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i, b_j = nj_branch_lengths ~d_ij ~s_i ~s_j ~n_act_minus_2_max in
        let u = do_merge ~dim ~st ~i ~j ~d_ij ~b_i ~b_j in
        let idx_cands = strategy.expand u in
        let cands_u = ref [] in
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) st.nbrs.(i);
        Array.iter (fun (x, _) ->
          if x <> i && x <> j then cands_u := x :: !cands_u) st.nbrs.(j);
        List.iter (fun x ->
          if x <> i && x <> j then cands_u := x :: !cands_u) idx_cands;
        rebuild_nbrs u !cands_u;
        let affected =
          collect_affected ~i ~j ~u ~active:st.active ~rev_nbrs:st.rev_nbrs idx_cands in
        deactivate_pair ~i ~j ~nbrs:st.nbrs ~rev_remove ~active:st.active
          ~n_active:st.n_active;
        strategy.on_merge ~i ~j ~u;
        patch_affected ~u ~k_nn ~active:st.active ~nbrs:st.nbrs
          ~rebuild_nbrs ~expand:strategy.expand affected;
        strategy.rebuild ()
      done;
      finalize_root ~bound:max_slots ~active:st.active ~trees:st.trees ~dist_of
    (* Dense reference implementation *)
    let compute_dense ?(verbose = false)
        ?(index_type = Interfaiss.Type.of_string "hnsw(32)")
        ?(k_nn = 10) ?(row_sum = RowSum.Knn) ?(symmetry = Symmetry.One)
        names data =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let n = check_inputs __FUNCTION__ names data in
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
         | RowSum.Knn ->
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
        let i, j = best_q_pair ~n_act_minus_2 ~s ~dist_of cand in
        let d_ij = dist_of i j in
        let s_i = Float.Array.unsafe_get s i and s_j = Float.Array.unsafe_get s j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i, b_j = nj_branch_lengths ~d_ij ~s_i ~s_j ~n_act_minus_2_max in
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
        let new_emb = weighted_centroid ~dim ~si_f ~sj_f ~tot_f:total_f embedding.(i) embedding.(j) in
        embedding.(i) <- new_emb;
        size.(i) <- total;
        active.(j) <- false;
        decr n_active
      done;
      let root = finalize_root ~bound:n ~active ~trees ~dist_of in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (dense): built unrooted tree on %d leaves.\n%!"
          prefix n;
      root
    (* Periodic-rebuild FAISS implementation.
       K-NN-restricted scan-NJ (the shared [run_core] driver) backed by
       a persistent FAISS index plus
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
      let n = check_inputs __FUNCTION__ names data in
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
      let dist_of = make_dist_of ~n ~data ~merge_left ~merge_right ~merge_dist ~dist_cache in
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
      let rebuild_nbrs = make_rebuild_nbrs ~k_nn ~max_slots ~active ~nbrs ~rev_add ~rev_remove ~dist_of in
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
      let nbrs_idx = make_nbrs_idx ~max_slots ~nbrs in
      let row_sums = make_row_sums ~bound:max_slots ~active ~nbrs ~dist_of ~row_sum ~n_active in
      let st =
        { trees; active; size; embedding; merge_left; merge_right; merge_dist;
          nbrs; rev_nbrs; n_active; next_slot } in
      let strategy =
        { expand = faiss_expand;
          on_merge =
            (fun ~i ~j ~u ->
              tombstone_slot i;
              tombstone_slot j;
              new_slots := u :: !new_slots);
          rebuild =
            (fun () ->
              incr merges_since_rebuild;
              if !merges_since_rebuild >= rebuild_interval then
                rebuild_index ()) } in
      let root =
        run_core ~max_slots ~dim ~k_nn ~symmetry ~st ~dist_of ~rebuild_nbrs
          ~row_sums ~nbrs_idx ~rev_remove ~strategy in
      (match !cur_index with
       | Some idx -> Interfaiss.delete idx; cur_index := None
       | None -> ());
      if verbose then
        Printf.eprintf "%s Sparse-NJ (periodic-rebuild): built unrooted tree on %d leaves (%d rebuilds).\n%!"
          prefix n !n_rebuilds;
      root
    (* Random-projection-forest periodic-rebuild implementation.
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
      let n = check_inputs __FUNCTION__ names data in
      let dim = Float.Array.length data.(0) in
      let max_slots = 2 * n - 3 in
      let trees = Array.make max_slots (Trees.Newick.leaf "") in
      let active = Array.make max_slots false in
      let size = Array.make max_slots 0 in
      let embedding = Array.make max_slots (Float.Array.make 0 0.) in
      let merge_left = Array.make max_slots (-1) in
      let merge_right = Array.make max_slots (-1) in
      let merge_dist = Array.make max_slots 0. in
      let use_hyp = distance = Distance.Hyperbolic in
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
      let dist_of_saitou_nei = make_dist_of ~n ~data ~merge_left ~merge_right ~merge_dist ~dist_cache in
      (* Bulk distance used by retrieval / K-NN lists / row sums / heap
         ordering. *)
      let dist_of =
        match distance with
        | Distance.SaitouNei -> dist_of_saitou_nei
        | Distance.Centroid -> dist_of_centroid
        | Distance.Hyperbolic -> dist_of_hyp in
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
      let rebuild_forest () =
        let act = ref [] in
        for s = max_slots - 1 downto 0 do
          if active.(s) then act := (s, embedding.(s)) :: !act
        done;
        let arr = Array.of_list !act in
        if Array.length arr > 0 then
          cur_forest :=
            Some
              (Forest.build ~n_trees ~leaf_size ~seed:(0xC0FFEE + !n_rebuilds)
                 arr)
        else
          cur_forest := None;
        new_slots := [];
        merges_since_rebuild := 0;
        incr n_rebuilds in
      let k_query = k_nn * k_query_factor in
      let forest_expand v =
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
        !res in
      let rebuild_nbrs = make_rebuild_nbrs ~k_nn ~max_slots ~active ~nbrs ~rev_add ~rev_remove ~dist_of in
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
      let heap = QHeap.create () in
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
        Array.iter (fun (x, _) -> push_pair_unordered v x) nbrs.(v) in
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
      (* Skip the RowSum CLI selector for this mode: K-NN mean is the
         only estimator the heap supports.  RowSum.Full would be a
         possible extension but is not the validated path. *)
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
        n_act_minus_2 *. dist_of i j
        -. n_act_minus_1
           *. (Float.Array.unsafe_get r_hat i
               +. Float.Array.unsafe_get r_hat j) in
      let pop_best () =
        let budget = max 4 (4 * k_nn) in
        let fresh = ref [] in
        let n_fresh = ref 0 in
        let going = ref true in
        while !going do
          match QHeap.pop heap with
          | None -> going := false
          | Some e ->
            if not active.(e.QHeap.i) || not active.(e.QHeap.j) then
              ()
            else if version.(e.QHeap.i) <> e.QHeap.vi
                 || version.(e.QHeap.j) <> e.QHeap.vj then begin
              let qp = q_prime e.QHeap.i e.QHeap.j in
              QHeap.push heap
                { e with QHeap.q_prime = qp;
                         vi = version.(e.QHeap.i);
                         vj = version.(e.QHeap.j) }
            end else begin
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
        !best_pair in
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
        let d_ij = dist_of i j in
        let n_act_minus_1 = float_of_int (!n_active - 1) in
        let n_act_minus_2 = float_of_int (!n_active - 2) in
        let s_i = n_act_minus_1 *. Float.Array.unsafe_get r_hat i in
        let s_j = n_act_minus_1 *. Float.Array.unsafe_get r_hat j in
        let n_act_minus_2_max = max n_act_minus_2 1. in
        let b_i, b_j = nj_branch_lengths ~d_ij ~s_i ~s_j ~n_act_minus_2_max in
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
        embedding.(u) <- weighted_centroid ~dim ~si_f ~sj_f ~tot_f embedding.(i) embedding.(j);
        (* Hyperbolic Q-distance: place u on the geodesic between its
           parents at the NJ branch length b_i from i.  The Euclidean
           centroid above is still kept for RP-forest retrieval; only
           the Q-formula / row-sum distances are hyperbolic. *)
        if use_hyp then
          hyp_pos.(u) <- hyp_geodesic hyp_pos.(i) hyp_pos.(j) b_i;
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
        let affected = collect_affected ~i ~j ~u ~active ~rev_nbrs forest_cands in
        deactivate_pair ~i ~j ~nbrs ~rev_remove ~active ~n_active;
        new_slots := u :: !new_slots;
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
        push_pairs_of u;
        Hashtbl.iter (fun v () -> push_pairs_of v) affected;
        incr merges_since_rebuild;
        if !merges_since_rebuild >= rebuild_interval then
          rebuild_forest ()
      done;
      cur_forest := None;
      let root = finalize_root ~bound:max_slots ~active ~trees ~dist_of in
      if verbose then
        Printf.eprintf "%s Sparse-NJ (rp-forest): built unrooted tree on %d leaves (%d rebuilds).\n%!"
          prefix n !n_rebuilds;
      root
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
        val of_string: string -> t
        val to_string: t -> string
      end
    module Mode:
      sig
        type t =
          | Dense
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
       [mode] selects the engine: [Dense] is the validated O(n^3)-time /
       O(n^2)-memory reference (best quality, the default);
       [PeriodicRebuild] is the recommended exact Theta(n^2) engine
       (K-NN-restricted scan-NJ over a persistent FAISS index with a
       sqrt(n)-periodic rebuild); [RPForestRebuild] uses a native
       RP-forest with a heap-driven merge loop and, under the
       [Hyperbolic] geometric-proxy [distance], gives the
       sub-quadratic-memory fallback.  [index_type] picks the FAISS index
       (PeriodicRebuild); [k_query_factor] sets the K-NN over-fetch
       factor; [distance] selects the rp-forest Q-distance model;
       [hyp_scale] sets the radial scale s for the hyperboloid lift
       (empirically s in [0.5, 1.0]; default 1.0); [row_sum] and
       [symmetry] tune the NJ row-sum estimator and the K-NN candidate
       symmetry. *)
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
