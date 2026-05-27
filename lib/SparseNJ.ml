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
        let of_string = function
          | "dense" -> Dense
          | "subquadratic" -> Subquadratic
          | "hyperbolic" -> Hyperbolic
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "sparse-NJ mode" s
        let to_string = function
          | Dense -> "dense"
          | Subquadratic -> "subquadratic"
          | Hyperbolic -> "hyperbolic"
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
       Hyperbolic-embedding implementation.

       Hyperboloid model conventions:  a point p in H^d is stored as a
       (d + 1)-vector with p.(0) the time-like coord (positive) and
       p.(1..d) the spatial coords, satisfying the Lorentz constraint
       -p.(0)^2 + p.(1)^2 + ... + p.(d)^2 = -1.  The Lorentz inner
       product is then  <p, q>_L = -p.(0) q.(0) + p.(1) q.(1) + ... ;
       hyperbolic distance is  acosh(-<p, q>_L).

       Initial leaf positions are lifted from the principal-coord
       Euclidean embedding via
           p_i -> (cosh(s |p_i|), sinh(s |p_i|) * p_i / |p_i|)
       with [scale] = s a tunable radial scale.  Empirically s in
       [0.5, 1.0] gives Spearman >= 0.9 against true Saitou-Nei
       distances across the whole NJ run (Phase 0 diagnostic on
       prot_k5_kt0_035; see PhyloSplits-Subquadratic.tex Path 3).

       Per-merge update places the new cluster at the geodesic point
       at hyperbolic distance b_i from p_i (toward p_j) using the
       closed-form expression
           p_u = cosh(b_i) p_i + sinh(b_i) v_{ij}
           v_{ij} = (p_j - cosh(d_H(p_i, p_j)) p_i) / sinh(d_H(p_i, p_j))
       which keeps p_u on the hyperboloid by construction.  Under
       perfect tree-metric additivity p_u then satisfies
       d_H(p_u, p_x) = d_NJ(u, x) for every other active x; under
       near-additivity (real data) the equality is approximate, and
       Phase 0 confirms the approximation tracks Saitou-Nei
       distances closely enough for K-NN selection. *)
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
    module Mode:
      sig
        type t =
          | Dense
          | Subquadratic
          | Hyperbolic
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
      ?row_sum:RowSum.t ->
      ?symmetry:Symmetry.t ->
      string array ->
      Float.Array.t array ->
      Trees.Newick.t
  end
)
