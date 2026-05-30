(*
    CoverTree.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    CoverTree.ml implements a Beygelzimer-Langford-Ravikumar 2006 cover
    tree as a functor over an arbitrary metric.  Cover trees support
    insert / delete / K-NN in O(c^{O(1)} log n) amortised under bounded
    doubling dimension; for the centroid-Euclidean metric KPop's
    sparse-NJ uses, c is small and the bounds are effectively
    O(d log n) per operation.

    The structure maintains, at every integer level L, a 2^L-net of the
    current point set: a maximal set of points pairwise > 2^L apart
    such that every other point lies within 2^L of some net point.
    Higher levels (large 2^L) are sparser; lower levels approach
    the full set as L -> -infinity.

    Representation: explicit cover tree.  Each node carries a point
    and a hash table of children indexed by their level; a child at
    level L means d(parent, child) <= 2^parent_level and child is the
    representative at level L of its subtree.  The level-indexed table
    plus a cached min_child_level lets the K-NN descent fetch
    next-level children in O(|children at next level|) instead of
    O(|all children|), and lets us decide whether to descend further
    in O(1).

    The delete operation here uses the straightforward
    "collect-descendants and reinsert" approach: simpler to verify
    correct than the paper's promote-children algorithm, same big-O
    in expectation.

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

open BiOCamLib.Better

include (
  struct
    module type Metric =
      sig
        type t
        val equal: t -> t -> bool
        val dist: t -> t -> float
      end
    module Make (M: Metric) =
      struct
        type point = M.t
        type node = {
          point: point;
          (* Children indexed by their attachment level.  Entries are
             never empty: when the last child at a level is removed
             the key is dropped from the table. *)
          children: (int, node list) Hashtbl.t;
          (* Minimum level among any child; max_int if no children.
             Lets the K-NN descent check "any child below cur_level"
             in O(1) without iterating the table. *)
          mutable min_child_level: int
        }
        type t = {
          mutable root: (node * int) option (* (root_node, root_level) *)
        }
        let create () = { root = None }
        let is_empty t = t.root = None
        let make_node p = {
          point = p;
          children = Hashtbl.create 4;
          min_child_level = max_int
        }
        (* 2^level as an IEEE float; ldexp is much cheaper than (**). *)
        let cover_radius level = Stdlib.ldexp 1.0 level
        let add_child parent child child_level =
          let existing =
            match Hashtbl.find_opt parent.children child_level with
            | Some l -> l
            | None -> [] in
          Hashtbl.replace parent.children child_level (child :: existing);
          if child_level < parent.min_child_level then
            parent.min_child_level <- child_level
        let recompute_min_child_level n =
          let m = ref max_int in
          Hashtbl.iter (fun l _ -> if l < !m then m := l) n.children;
          n.min_child_level <- !m
        let rec iter_subtree f n =
          f n.point;
          Hashtbl.iter
            (fun _ cs -> List.iter (iter_subtree f) cs) n.children
        let size t =
          match t.root with
          | None -> 0
          | Some (r, _) ->
            let c = ref 0 in
            iter_subtree (fun _ -> incr c) r;
            !c
        (* Descend the cover.  q is the current cover set (nodes at
           cur_level that all satisfy d(p, n) <= 2^cur_level).  At
           each step we expand q to its children at the next level
           (and q itself implicitly), filter to those within
           2^(cur_level - 1) of p, and recurse.  When the filtered
           set is empty p attaches as a new child of any node in q
           at level cur_level - 1. *)
        let rec do_insert ~q ~cur_level p =
          let next_level = cur_level - 1 in
          let next_radius = cover_radius next_level in
          let next_cands = ref [] in
          List.iter
            (fun parent ->
              if M.dist parent.point p <= next_radius then
                next_cands := parent :: !next_cands;
              match Hashtbl.find_opt parent.children next_level with
              | None -> ()
              | Some cs ->
                List.iter
                  (fun c ->
                    if M.dist c.point p <= next_radius then
                      next_cands := c :: !next_cands)
                  cs)
            q;
          match !next_cands with
          | [] ->
            let parent =
              match
                List.find_opt
                  (fun n -> M.dist n.point p <= cover_radius cur_level) q
              with
              | Some n -> n
              | None -> List.hd q in
            add_child parent (make_node p) next_level
          | next_q ->
            do_insert ~q:next_q ~cur_level:next_level p
        let insert t p =
          match t.root with
          | None ->
            t.root <- Some (make_node p, 0)
          | Some (root, root_level) ->
            let d = M.dist root.point p in
            if d = 0. && M.equal root.point p then
              ()
            else begin
              (* Choose a top level large enough for the root to cover p:
                 2^top_level >= d. *)
              let top_level =
                let rec find l =
                  if cover_radius l >= d then l else find (l + 1) in
                find root_level in
              if top_level > root_level then
                t.root <- Some (root, top_level);
              do_insert ~q:[root] ~cur_level:top_level p
            end
        (* K-NN query.  Maintains best as parallel arrays of size k
           sorted ascending by distance (O(k) insert by shift, O(1)
           kth-distance lookup, no allocations).  The cover set q
           descends by level, carrying each node's cached distance
           to query so pruning is recomputation-free. *)
        let knn t query k =
          if k <= 0 then []
          else
            match t.root with
            | None -> []
            | Some (root, root_level) ->
              let best_d = Array.make k infinity in
              let best_p = Array.make k root.point in
              let n_best = ref 0 in
              let kth_dist () =
                if !n_best < k then infinity else best_d.(k - 1) in
              let add_candidate d p =
                if d < kth_dist () then begin
                  let pos = ref !n_best in
                  while !pos > 0 && best_d.(!pos - 1) > d do
                    if !pos < k then begin
                      best_d.(!pos) <- best_d.(!pos - 1);
                      best_p.(!pos) <- best_p.(!pos - 1)
                    end;
                    decr pos
                  done;
                  if !pos < k then begin
                    best_d.(!pos) <- d;
                    best_p.(!pos) <- p
                  end;
                  if !n_best < k then incr n_best
                end in
              let root_d = M.dist root.point query in
              add_candidate root_d root.point;
              (* q carries (node, d(node, query)) to avoid recomputing
                 distance during pruning. *)
              let q = ref [(root, root_d)] in
              let cur_level = ref root_level in
              let going = ref true in
              while !going do
                let next_level = !cur_level - 1 in
                let next_radius = cover_radius next_level in
                (* Expand q to its self-edges plus explicit children at
                   next_level. *)
                let expanded = ref [] in
                List.iter
                  (fun (n, d_n) ->
                    expanded := (n, d_n) :: !expanded;
                    match Hashtbl.find_opt n.children next_level with
                    | None -> ()
                    | Some cs ->
                      List.iter
                        (fun c ->
                          let d = M.dist c.point query in
                          add_candidate d c.point;
                          expanded := (c, d) :: !expanded)
                        cs)
                  !q;
                let prune_bound = kth_dist () +. next_radius in
                let pruned =
                  List.filter
                    (fun (_, d) -> d <= prune_bound) !expanded in
                if pruned = [] then
                  going := false
                else begin
                  q := pruned;
                  cur_level := next_level;
                  let any_below =
                    List.exists
                      (fun (n, _) -> n.min_child_level < next_level) pruned in
                  if not any_below then going := false
                end
              done;
              let result = ref [] in
              for i = !n_best - 1 downto 0 do
                result := best_p.(i) :: !result
              done;
              !result
        (* Remove [p] from the tree: find the unique node carrying p,
           splice it from its parent, collect all descendants, then
           reinsert each one.  O(|subtree(p)|) reinsertions in the
           worst case; O(log n) amortised under bounded doubling
           dimension. *)
        let collect_subtree_points n =
          let acc = ref [] in
          Hashtbl.iter
            (fun _ cs ->
              List.iter
                (fun c -> iter_subtree (fun pt -> acc := pt :: !acc) c) cs)
            n.children;
          !acc
        let remove t p =
          match t.root with
          | None -> false
          | Some (root, _) when M.equal root.point p ->
            let descendants = collect_subtree_points root in
            t.root <- None;
            List.iter (fun pt -> insert t pt) descendants;
            true
          | Some (root, _) ->
            let found = ref false in
            let descendants = ref [] in
            let rec walk parent =
              if !found then ()
              else begin
                let levels =
                  Hashtbl.fold (fun lvl _ acc -> lvl :: acc) parent.children [] in
                List.iter
                  (fun lvl ->
                    if !found then ()
                    else begin
                      let cs = Hashtbl.find parent.children lvl in
                      let kept =
                        List.filter
                          (fun c ->
                            if (not !found) && M.equal c.point p then begin
                              found := true;
                              List.iter
                                (fun pt -> descendants := pt :: !descendants)
                                (collect_subtree_points c);
                              false
                            end else begin
                              walk c;
                              true
                            end)
                          cs in
                      if kept = [] then
                        Hashtbl.remove parent.children lvl
                      else if List.length kept <> List.length cs then
                        Hashtbl.replace parent.children lvl kept
                    end)
                  levels;
                if !found then recompute_min_child_level parent
              end in
            walk root;
            if !found then
              List.iter (fun pt -> insert t pt) !descendants;
            !found
      end
  end: sig
    module type Metric =
      sig
        type t
        val equal: t -> t -> bool
        val dist: t -> t -> float
      end
    module Make (M: Metric):
      sig
        type point = M.t
        type t
        val create: unit -> t
        val is_empty: t -> bool
        val size: t -> int
        val insert: t -> point -> unit
        val remove: t -> point -> bool
        val knn: t -> point -> int -> point list
      end
  end
)
