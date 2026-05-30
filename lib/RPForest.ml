(*
    RPForest.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    RPForest.ml implements a random projection forest as a functor over
    a point type carrying a Float.Array embedding.  Each tree partitions
    the point set by a sequence of random unit-vector projections,
    splitting at the median of the projected coordinate; the forest is
    an ensemble of M independent such trees.  K-NN queries descend each
    tree to a leaf, union the candidate sets, then re-rank by exact
    embedding distance and return the top-K.

    The structure is intentionally static: build-once, query-many.  In
    the sparse-NJ periodic-rebuild loop a fresh forest is built every
    sqrt(n) merges over the current active centroid set, so dynamic
    insert/delete are not needed.  Build cost is O(M n log n);
    K-NN query cost is O(M log n + M * leaf_size * d) which is
    sub-linear in n at fixed (M, leaf_size, d).

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
    module type Point =
      sig
        type t
        val equal: t -> t -> bool
        val hash: t -> int
        val embedding: t -> Float.Array.t
      end
    module Make (P: Point) =
      struct
        module PSet = Hashtbl.Make (struct
          type t = P.t
          let equal = P.equal
          let hash = P.hash
        end)
        type point = P.t
        type node =
          | Leaf of point array
          | Split of {
              direction: Float.Array.t;
              threshold: float;
              left: node;
              right: node
            }
        type t = {
          trees: node array;
          n_points: int;
          dim: int
        }
        let dot a b =
          let n = Float.Array.length a in
          let acc = ref 0. in
          for k = 0 to n - 1 do
            acc :=
              !acc
              +. Float.Array.unsafe_get a k *. Float.Array.unsafe_get b k
          done;
          !acc
        let sq_dist a b =
          let n = Float.Array.length a in
          let acc = ref 0. in
          for k = 0 to n - 1 do
            let delta =
              Float.Array.unsafe_get a k -. Float.Array.unsafe_get b k in
            acc := !acc +. delta *. delta
          done;
          !acc
        let eucl_dist a b = sqrt (sq_dist a b)
        (* Box-Muller normal sample.  Stateless; call twice for two
           independent samples.  Guards against log(0). *)
        let normal rng =
          let u1 = max 1e-15 (Random.State.float rng 1.0) in
          let u2 = Random.State.float rng 1.0 in
          sqrt (-2. *. log u1) *. cos (2. *. Float.pi *. u2)
        let random_unit_vector rng dim =
          let v = Float.Array.init dim (fun _ -> normal rng) in
          let norm = sqrt (dot v v) in
          if norm > 0. then
            for k = 0 to dim - 1 do
              Float.Array.unsafe_set v k
                (Float.Array.unsafe_get v k /. norm)
            done;
          v
        let rec build_node ~leaf_size ~dim rng points =
          if Array.length points <= leaf_size then Leaf points
          else begin
            let direction = random_unit_vector rng dim in
            let projs =
              Array.map (fun p -> dot direction (P.embedding p)) points in
            let sorted = Array.copy projs in
            Array.sort compare sorted;
            let med = sorted.(Array.length sorted / 2) in
            let l = ref [] and r = ref [] in
            Array.iteri
              (fun i p ->
                if projs.(i) < med then l := p :: !l
                else r := p :: !r)
              points;
            let la = Array.of_list !l and ra = Array.of_list !r in
            (* Guard against degenerate splits (all projections equal). *)
            if Array.length la = 0 || Array.length ra = 0 then Leaf points
            else
              Split {
                direction;
                threshold = med;
                left = build_node ~leaf_size ~dim rng la;
                right = build_node ~leaf_size ~dim rng ra
              }
          end
        let build ?(leaf_size = 16) ?(n_trees = 10) ?(seed = 0xC0FFEE)
            points =
          let n_points = Array.length points in
          if n_points = 0 then
            { trees = [||]; n_points = 0; dim = 0 }
          else begin
            let dim = Float.Array.length (P.embedding points.(0)) in
            let base_rng = Random.State.make [| seed |] in
            let trees =
              Array.init n_trees
                (fun _ ->
                  let s = Random.State.bits base_rng in
                  let tree_rng = Random.State.make [| s |] in
                  build_node ~leaf_size ~dim tree_rng points) in
            { trees; n_points; dim }
          end
        let size t = t.n_points
        let rec descend node q_emb =
          match node with
          | Leaf pts -> pts
          | Split { direction; threshold; left; right } ->
            let p = dot direction q_emb in
            if p < threshold then descend left q_emb
            else descend right q_emb
        let knn t query k =
          if k <= 0 || t.n_points = 0 then []
          else begin
            let q_emb = P.embedding query in
            let cand = PSet.create 32 in
            Array.iter
              (fun tree ->
                let leaf = descend tree q_emb in
                Array.iter
                  (fun p ->
                    if not (PSet.mem cand p) then PSet.add cand p ())
                  leaf)
              t.trees;
            let pairs = ref [] in
            PSet.iter
              (fun p () ->
                let d = eucl_dist (P.embedding p) q_emb in
                pairs := (d, p) :: !pairs)
              cand;
            let arr = Array.of_list !pairs in
            Array.sort (fun (a, _) (b, _) -> compare a b) arr;
            let m = min k (Array.length arr) in
            Array.init m (fun i -> snd arr.(i)) |> Array.to_list
          end
      end
  end: sig
    module type Point =
      sig
        type t
        val equal: t -> t -> bool
        val hash: t -> int
        val embedding: t -> Float.Array.t
      end
    module Make (P: Point):
      sig
        type point = P.t
        type t
        val build:
          ?leaf_size:int -> ?n_trees:int -> ?seed:int ->
          point array -> t
        val size: t -> int
        val knn: t -> point -> int -> point list
      end
  end
)
