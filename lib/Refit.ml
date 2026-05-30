(*
    Refit.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Refit.ml re-estimates the branch lengths of a fixed tree topology by
    sparse ordinary least squares against a set of target pairwise
    distances, holding the topology fixed.  This is the optional fourth
    stage of the sparse-NJ pipeline: sparse-NJ produces the topology
    (and native twisted-L2 branch lengths) sub-quadratically; this
    module replaces the branch lengths with ones fitted to a better-
    calibrated distance (e.g. mash-like Jaccard, see Jaccard.ml) while
    keeping the topology.

    The branch lengths b solve  X b ~= d  in least squares, where each
    row of X is a leaf pair on a chosen O(n K) subset, each column is a
    tree edge, X[pair, edge] = 1 iff the edge lies on the path between
    the pair, and d is the target distance for that pair.  The
    over-determined sparse system (2n - 3 unknowns, O(n K) equations) is
    solved by LSQR (Paige-Saunders 1982), which needs only sparse
    matrix-vector products and so stays O(n K) per iteration.

    The target distances are supplied through the [dist] callback, so
    this module is independent of how distances are computed.

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
    (* Pair-selection strategy for the OLS system.  [Knn]: each leaf's
       K topologically-nearest leaves (reuses the sparse-NJ regime,
       O(n K) equations).  [PerBranch]: one pair per edge, the minimum
       set that covers every branch, O(n) equations.  [All]: every
       pair, O(n^2) (gold standard, small n only). *)
    module Strategy =
      struct
        type t =
          | Knn
          | PerBranch
          | All
        let of_string = function
          | "knn" | "reuse-knn" -> Knn
          | "per-branch" | "perbranch" -> PerBranch
          | "all" -> All
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "refit pair-selection strategy" s
        let to_string = function
          | Knn -> "knn"
          | PerBranch -> "per-branch"
          | All -> "all"
      end
    (* Sparse over-determined least-squares solve  min ||A x - b||  by
       LSQR (Paige & Saunders, ACM TOMS 1982).  [matvec] computes A v,
       [rmatvec] computes A^T u; [m] rows, [n] cols.  Damping-free,
       fixed iteration cap with a relative-residual stopping test. *)
    let lsqr ~m ~n ~matvec ~rmatvec b ~iters ~tol =
      let x = Float.Array.make n 0. in
      let u = Float.Array.copy b in
      let beta = ref (sqrt (let s = ref 0. in
        Float.Array.iter (fun v -> s := !s +. v *. v) u; !s)) in
      if !beta > 0. then
        Float.Array.iteri (fun i v -> Float.Array.unsafe_set u i (v /. !beta)) u;
      let v = rmatvec u in
      let alpha = ref (sqrt (let s = ref 0. in
        Float.Array.iter (fun e -> s := !s +. e *. e) v; !s)) in
      if !alpha > 0. then
        Float.Array.iteri (fun i e -> Float.Array.unsafe_set v i (e /. !alpha)) v;
      let w = Float.Array.copy v in
      let phi_bar = ref !beta and rho_bar = ref !alpha in
      let b_norm = !beta in
      (try
        for _ = 1 to iters do
          (* Continue the bidiagonalisation: u = A v - alpha u *)
          let au = matvec v in
          Float.Array.iteri
            (fun i a -> Float.Array.unsafe_set u i (a -. !alpha *. Float.Array.unsafe_get u i)) au;
          beta := sqrt (let s = ref 0. in
            Float.Array.iter (fun e -> s := !s +. e *. e) u; !s);
          if !beta > 0. then
            Float.Array.iteri (fun i e -> Float.Array.unsafe_set u i (e /. !beta)) u;
          (* v = A^T u - beta v *)
          let av = rmatvec u in
          Float.Array.iteri
            (fun i a -> Float.Array.unsafe_set v i (a -. !beta *. Float.Array.unsafe_get v i)) av;
          alpha := sqrt (let s = ref 0. in
            Float.Array.iter (fun e -> s := !s +. e *. e) v; !s);
          if !alpha > 0. then
            Float.Array.iteri (fun i e -> Float.Array.unsafe_set v i (e /. !alpha)) v;
          (* Givens rotation on the bidiagonal *)
          let rho = sqrt (!rho_bar *. !rho_bar +. !beta *. !beta) in
          let c = !rho_bar /. rho and s = !beta /. rho in
          let theta = s *. !alpha in
          rho_bar := -. c *. !alpha;
          let phi = c *. !phi_bar in
          phi_bar := s *. !phi_bar;
          (* Update x and w *)
          let phi_over_rho = phi /. rho and theta_over_rho = theta /. rho in
          Float.Array.iteri
            (fun i wi ->
              Float.Array.unsafe_set x i
                (Float.Array.unsafe_get x i +. phi_over_rho *. wi);
              Float.Array.unsafe_set w i
                (Float.Array.unsafe_get v i -. theta_over_rho *. wi))
            w;
          (* |phi_bar| is the current residual norm *)
          if !phi_bar <= tol *. b_norm then raise Exit
        done
      with Exit -> ());
      x
    (* Build, from a Newick tree, the data the refit needs:
         leaf_names : array of leaf names in flat order;
         leaf_node  : array of leaf node indices into the flat array;
         root_path  : per leaf, the list of edge ids (= the child node
                      index of each edge) from the leaf up to the root;
         n_edges    : number of edges (non-root nodes). *)
    let analyse tree =
      let fl = Trees.Newick.dfs_flatten tree in
      let n_nodes = Array.length fl in
      (* edge id for a node = its own index (its parent edge); the root
         (parent_idx = -1) owns no edge.  Renumber edges 0.. densely. *)
      let parent = Array.make n_nodes (-1) in
      let is_leaf = Array.make n_nodes false in
      let leaf_names = ref [] and leaf_nodes = ref [] in
      Array.iteri
        (fun i (_, par, node, descs) ->
          parent.(i) <- par;
          if Array.length descs = 0 then begin
            is_leaf.(i) <- true;
            leaf_names := Trees.Newick.get_node_name node :: !leaf_names;
            leaf_nodes := i :: !leaf_nodes
          end)
        fl;
      let leaf_nodes = Array.of_rlist !leaf_nodes
      and leaf_names = Array.of_rlist !leaf_names in
      (* dense edge ids: every non-root node gets one *)
      let edge_id = Array.make n_nodes (-1) in
      let n_edges = ref 0 in
      for i = 0 to n_nodes - 1 do
        if parent.(i) <> -1 then begin
          edge_id.(i) <- !n_edges;
          incr n_edges
        end
      done;
      let root_path lf =
        let acc = ref [] and cur = ref lf in
        while parent.(!cur) <> -1 do
          acc := edge_id.(!cur) :: !acc;
          cur := parent.(!cur)
        done;
        !acc in
      let paths = Array.map root_path leaf_nodes in
      fl, parent, edge_id, leaf_names, leaf_nodes, paths, !n_edges
    (* Topological K-NN among leaves: distances with every edge length
       set to 1, via dijkstra from each leaf.  Returns, for each leaf
       index, the indices of its K nearest other leaves. *)
    let topological_knn tree leaf_nodes k =
      let unit_tree =
        Trees.Newick.dfs_map
          (fun node _ -> node)
          (fun _ edge -> Trees.Newick.set_edge_length edge 1.)
          tree in
      let fl = Trees.Newick.dfs_flatten unit_tree in
      let n_leaves = Array.length leaf_nodes in
      Array.init n_leaves
        (fun i ->
          let d = Trees.Newick.dijkstra fl leaf_nodes.(i) in
          let cand =
            Array.init n_leaves
              (fun j -> Float.Array.get d leaf_nodes.(j), j) in
          Array.sort (fun (a, _) (b, _) -> compare a b) cand;
          (* skip self (distance 0 at position 0) *)
          let kk = min k (n_leaves - 1) in
          Array.init kk (fun p -> snd cand.(p + 1)))
    (* Refit branch lengths of [tree] by sparse OLS against [dist].
       [dist a b] returns the target distance between leaf names a, b. *)
    let refit ?(strategy = Strategy.Knn) ?(k = 10) ?(iters = 200)
        ?(tol = 1e-8) ?(verbose = false) ~dist tree =
      let open String.TermIO in
      let prefix = grey (Printf.sprintf "(%s):" __FUNCTION__) in
      let fl, _parent, edge_id, leaf_names, leaf_nodes, paths, n_edges =
        analyse tree in
      let n_leaves = Array.length leaf_names in
      let name_of_leaf = leaf_names in
      (* Select the pair set as index pairs into leaf_names. *)
      let pairs =
        match strategy with
        | Strategy.All ->
          let acc = ref [] in
          for i = 0 to n_leaves - 1 do
            for j = i + 1 to n_leaves - 1 do
              acc := (i, j) :: !acc
            done
          done;
          Array.of_rlist !acc
        | Strategy.Knn ->
          let nbrs = topological_knn tree leaf_nodes k in
          let seen = Hashtbl.create (n_leaves * k) in
          let acc = ref [] in
          Array.iteri
            (fun i row ->
              Array.iter
                (fun j ->
                  let lo, hi = if i < j then i, j else j, i in
                  if lo <> hi && not (Hashtbl.mem seen (lo, hi)) then begin
                    Hashtbl.add seen (lo, hi) ();
                    acc := (lo, hi) :: !acc
                  end)
                row)
            nbrs;
          Array.of_rlist !acc
        | Strategy.PerBranch ->
          (* nearest cross-edge pair per edge, approximated by the
             topological 1-NN of each leaf (covers every pendant edge);
             internal edges are covered by the union of K=1 NN pairs. *)
          let nbrs = topological_knn tree leaf_nodes 1 in
          let seen = Hashtbl.create n_leaves in
          let acc = ref [] in
          Array.iteri
            (fun i row ->
              if Array.length row > 0 then begin
                let j = row.(0) in
                let lo, hi = if i < j then i, j else j, i in
                if not (Hashtbl.mem seen (lo, hi)) then begin
                  Hashtbl.add seen (lo, hi) ();
                  acc := (lo, hi) :: !acc
                end
              end)
            nbrs;
          Array.of_rlist !acc in
      let n_eq = Array.length pairs in
      (* Sparse incidence rows: for pair p, the path edges = symmetric
         difference of the two leaves' root-paths. *)
      let row_edges =
        Array.map
          (fun (i, j) ->
            let si = List.fold_left (fun s e -> IntSet.add e s) IntSet.empty paths.(i) in
            let sj = List.fold_left (fun s e -> IntSet.add e s) IntSet.empty paths.(j) in
            let diff = IntSet.union (IntSet.diff si sj) (IntSet.diff sj si) in
            IntSet.elements diff |> Array.of_list)
          pairs in
      let target =
        Float.Array.init n_eq
          (fun p ->
            let i, j = pairs.(p) in
            dist name_of_leaf.(i) name_of_leaf.(j)) in
      if verbose then
        Printf.eprintf "%s OLS refit: %d edges, %d pair-equations (%s, K=%d).\n%!"
          prefix n_edges n_eq (Strategy.to_string strategy) k;
      (* A v : for each pair-row, sum of v over the row's edges. *)
      let matvec v =
        Float.Array.init n_eq
          (fun p ->
            let s = ref 0. in
            Array.iter (fun e -> s := !s +. Float.Array.unsafe_get v e) row_edges.(p);
            !s) in
      (* A^T u : for each edge, sum of u over the rows touching it. *)
      let rmatvec u =
        let out = Float.Array.make n_edges 0. in
        Array.iteri
          (fun p edges ->
            let up = Float.Array.unsafe_get u p in
            Array.iter
              (fun e ->
                Float.Array.unsafe_set out e (Float.Array.unsafe_get out e +. up))
              edges)
          row_edges;
        out in
      let b = lsqr ~m:n_eq ~n:n_edges ~matvec ~rmatvec target ~iters ~tol in
      (* Clip negatives to zero. *)
      Float.Array.iteri
        (fun i v -> if v < 0. then Float.Array.unsafe_set b i 0.) b;
      (* Write the fitted lengths back by reconstructing the tree
         bottom-up from the flat representation, so the result is
         independent of any traversal-order assumption.  Children have
         higher pre-order index than their parent, so processing in
         reverse index order builds each subtree before its parent.
         The edge above node ci carries length b.(edge_id.(ci)). *)
      let n_nodes = Array.length fl in
      let subtree = Array.make n_nodes (Trees.Newick.leaf "") in
      for i = n_nodes - 1 downto 0 do
        let _, _, node, children = fl.(i) in
        let name = Trees.Newick.get_node_name node in
        subtree.(i) <-
          if Array.length children = 0 then
            Trees.Newick.leaf name
          else
            Trees.Newick.join ~name
              (Array.map
                 (fun (_, ci) ->
                   let len = Float.Array.get b edge_id.(ci) in
                   Trees.Newick.edge ~length:len (), subtree.(ci))
                 children)
      done;
      Trees.Newick.set_is_root subtree.(0) (Trees.Newick.get_is_root tree)
  end: sig
    module Strategy:
      sig
        type t =
          | Knn
          | PerBranch
          | All
        val of_string: string -> t
        val to_string: t -> string
      end
    val refit:
      ?strategy:Strategy.t -> ?k:int -> ?iters:int -> ?tol:float ->
      ?verbose:bool -> dist:(string -> string -> float) ->
      BiOCamLib.Trees.Newick.t -> BiOCamLib.Trees.Newick.t
  end
)
