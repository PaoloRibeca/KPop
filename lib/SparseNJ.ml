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
    unrooted-tree-style trifurcating root.  [make_splits] is a thin
    wrapper that converts the result into a Trees.Splits.t (one entry
    per non-trivial bipartition, weighted by the corresponding edge
    length), provided for ensemble or compatibility-filtered downstream
    consumers.

    This prototype uses the full pairwise distance matrix internally
    (O(n^2) memory) for simplicity.  A FAISS-driven version that
    maintains an incrementally-updated K-NN graph and never materialises
    the full distance matrix is the natural next step; the algorithm
    structure here is already set up for that swap.

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
       sum, falling back to classical NJ behaviour on a restricted
       candidate set. *)
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
    (* Compute the full pairwise Euclidean distance matrix from the
       row-stored embedding data.  Each row of [data] is a sample's
       coordinate vector in the (c - 1)-D twisted space. *)
    let pairwise_distances data =
      let n = Array.length data in
      let d_mat = Array.make_matrix n n 0. in
      for i = 0 to n - 1 do
        for j = i + 1 to n - 1 do
          let row_i = data.(i) and row_j = data.(j) in
          let dim = Float.Array.length row_i in
          let acc = ref 0. in
          for k = 0 to dim - 1 do
            let delta =
              Float.Array.unsafe_get row_i k -. Float.Array.unsafe_get row_j k in
            acc := !acc +. delta *. delta
          done;
          let v = sqrt !acc in
          d_mat.(i).(j) <- v;
          d_mat.(j).(i) <- v
        done
      done;
      d_mat
    (* For each row of the distance matrix, find the indices of its
       [k] nearest active neighbours.  [active] is the boolean mask
       of which row/column indices are still in play.  Returns an
       int array per row with at most [k] entries, sorted by ascending
       distance.  Self-index is always excluded. *)
    let knn_active d_mat active k_nn =
      let n = Array.length d_mat in
      Array.init n
        (fun i ->
          if not active.(i) then
            [||]
          else begin
            let buf = ref [] and count = ref 0 in
            for j = 0 to n - 1 do
              if j <> i && active.(j) then begin
                buf := (d_mat.(i).(j), j) :: !buf;
                incr count
              end
            done;
            let arr = Array.of_list !buf in
            Array.sort (fun (a, _) (b, _) -> compare a b) arr;
            let k = min k_nn !count in
            Array.init k (fun pos -> snd arr.(pos))
          end)
    (* Local row sum estimator for the surviving active set.  All sums
       use the current distance matrix d_mat and the active mask. *)
    let row_sums kind d_mat active nbrs n_active =
      let n = Array.length d_mat in
      let s = Float.Array.create n in
      Float.Array.fill s 0 n 0.;
      let n_active_minus_1 = float_of_int (n_active - 1) in
      (match kind with
       | RowSum.Full ->
         for i = 0 to n - 1 do
           if active.(i) then begin
             let acc = ref 0. in
             for j = 0 to n - 1 do
               if j <> i && active.(j) then
                 acc := !acc +. d_mat.(i).(j)
             done;
             Float.Array.unsafe_set s i !acc
           end
         done
       | RowSum.Knn ->
         (* sum over i's K-NN scaled by (n_active - 1) / K *)
         for i = 0 to n - 1 do
           if active.(i) then begin
             let acc = ref 0. in
             let arr = nbrs.(i) in
             let kk = Array.length arr in
             if kk > 0 then begin
               for p = 0 to kk - 1 do
                 acc := !acc +. d_mat.(i).(arr.(p))
               done;
               Float.Array.unsafe_set s i (!acc *. n_active_minus_1 /. float_of_int kk)
             end
           end
         done
       | RowSum.Topk ->
         (* mean over i's K-NN scaled by (n_active - 1) *)
         for i = 0 to n - 1 do
           if active.(i) then begin
             let acc = ref 0. in
             let arr = nbrs.(i) in
             let kk = Array.length arr in
             if kk > 0 then begin
               for p = 0 to kk - 1 do
                 acc := !acc +. d_mat.(i).(arr.(p))
               done;
               Float.Array.unsafe_set s i (!acc /. float_of_int kk *. n_active_minus_1)
             end
           end
         done);
      s
    (* Compute the candidate pair set: pairs (i, j) where j is among
       i's K-NN (One-sided) or both i and j are in each other's K-NN
       (Both).  Returns a list of (i, j) with i < j and no duplicates. *)
    let candidate_pairs sym nbrs active =
      let n = Array.length nbrs in
      let seen = Hashtbl.create (16 * n) in
      let res = ref [] in
      let add a b =
        let lo, hi = if a < b then a, b else b, a in
        if not (Hashtbl.mem seen (lo, hi)) then begin
          Hashtbl.add seen (lo, hi) ();
          res := (lo, hi) :: !res
        end in
      (match sym with
       | Symmetry.One ->
         for i = 0 to n - 1 do
           if active.(i) then
             Array.iter (fun j -> add i j) nbrs.(i)
         done
       | Symmetry.Both ->
         (* Materialise membership tests as sets per row *)
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
    (* The merge-time NJ algorithm: build a Newick tree directly.
       Each active cluster carries the current Newick subtree of
       leaves it spans; at each merge of (i, j) into a new node we
       join trees.(i) and trees.(j) under the new edge lengths
       b_i and b_j computed from the NJ Q-formula.  The loop runs
       while at least 4 active clusters remain; the final
       three-way merge of the last 3 active clusters produces the
       unrooted-tree-style trifurcating root. *)
    let compute ?(verbose = false)
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
      if verbose then
        Printf.eprintf "%s Computing %d x %d distance matrix...\n%!" prefix n n;
      let d_mat = pairwise_distances data in
      (* Per-active-cluster state: the Newick subtree spanning its
         leaves, plus an active mask.  Slots are reused as merges
         deactivate one of the two children. *)
      let trees = Array.init n (fun i -> Trees.Newick.leaf names.(i)) in
      let active = Array.make n true in
      let n_active = ref n in
      if verbose then
        Printf.eprintf "%s Sparse-NJ: K=%d, rowsum=%s, sym=%s.\n%!" prefix
          k_nn (RowSum.to_string row_sum) (Symmetry.to_string symmetry);
      (* NJ merge loop: continue while at least 4 active clusters
         remain so that we can finish with a clean 3-way root. *)
      while !n_active > 3 do
        let nbrs = knn_active d_mat active k_nn in
        let s = row_sums row_sum d_mat active nbrs !n_active in
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
              n_act_minus_2 *. d_mat.(i).(j)
              -. Float.Array.unsafe_get s i -. Float.Array.unsafe_get s j in
            if q < !best_q then begin
              best_q := q;
              best_pair := (i, j)
            end)
          cand;
        let i, j = !best_pair in
        let d_ij = d_mat.(i).(j) in
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
        (* Update the distance matrix row/col for the new node:
              d(u, x) = (d(i, x) + d(j, x) - d(i, j)) / 2.
           Then deactivate j. *)
        for x = 0 to n - 1 do
          if x <> i && x <> j && active.(x) then begin
            let new_d = 0.5 *. (d_mat.(i).(x) +. d_mat.(j).(x) -. d_ij) in
            d_mat.(i).(x) <- new_d;
            d_mat.(x).(i) <- new_d
          end
        done;
        active.(j) <- false;
        decr n_active
      done;
      (* Final three-way merge.  Identify the 3 remaining active
         indices and compute their pendant branch lengths from the
         pairwise distances among them, then join under a single
         trifurcating root.  This is the standard unrooted-NJ
         termination. *)
      let active_three =
        let acc = ref [] in
        for i = n - 1 downto 0 do
          if active.(i) then List.accum acc i
        done;
        Array.of_rlist !acc in
      assert (Array.length active_three = 3);
      let i1, i2, i3 = active_three.(0), active_three.(1), active_three.(2) in
      let d12 = d_mat.(i1).(i2) and d13 = d_mat.(i1).(i3) and d23 = d_mat.(i2).(i3) in
      let b1 = max 0. (0.5 *. (d12 +. d13 -. d23))
      and b2 = max 0. (0.5 *. (d12 +. d23 -. d13))
      and b3 = max 0. (0.5 *. (d13 +. d23 -. d12)) in
      let root =
        Trees.Newick.join
          [| Trees.Newick.edge ~length:b1 (), trees.(i1);
             Trees.Newick.edge ~length:b2 (), trees.(i2);
             Trees.Newick.edge ~length:b3 (), trees.(i3) |] in
      if verbose then
        Printf.eprintf "%s Sparse-NJ: built unrooted tree on %d leaves.\n%!"
          prefix n;
      root
    (* Thin wrapper: build the Newick tree first, then read off its
       internal bipartitions as a Trees.Splits.t.  The weights are
       the branch lengths of the corresponding edges, matching the
       Yggdrasill-side convention. *)
    let make_splits ?verbose ?k_nn ?row_sum ?symmetry names data =
      let tree = compute ?verbose ?k_nn ?row_sum ?symmetry names data in
      Trees.Splits.of_newick tree
  end: sig
    (* Row-sum estimator selector. *)
    module RowSum:
      sig
        type t =
          | Knn
          | Topk
          | Full
        val of_string: string -> t
        val to_string: t -> string
      end
    (* K-NN symmetry policy selector. *)
    module Symmetry:
      sig
        type t =
          | One
          | Both
        val of_string: string -> t
        val to_string: t -> string
      end
    (* Build an unrooted Newick tree by sparse-NJ on the row-stored
       embedding [data].  [names] is the array of leaf names matching
       the rows of [data].  Each internal merge contributes a join node
       carrying the NJ-computed branch lengths above its children; the
       final three-way merge produces a trifurcating unrooted root. *)
    val compute:
      ?verbose:bool ->
      ?k_nn:int ->
      ?row_sum:RowSum.t ->
      ?symmetry:Symmetry.t ->
      string array ->
      Float.Array.t array ->
      Trees.Newick.t
    (* Convenience wrapper: build the Newick tree via [compute] then
       extract its internal bipartitions as a Trees.Splits.t, weighted
       by edge length.  Provided for callers that want to feed the
       sparse-NJ output into the existing PhyloSplits / Yggdrasill
       ensembling pipeline.  Direct Newick output is the preferred
       path (see [compute]) and avoids the splits-to-tree round-trip. *)
    val make_splits:
      ?verbose:bool ->
      ?k_nn:int ->
      ?row_sum:RowSum.t ->
      ?symmetry:Symmetry.t ->
      string array ->
      Float.Array.t array ->
      Trees.Splits.t
  end
)
