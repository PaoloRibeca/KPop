(*
    Space.ml -- (c) 2022-2025 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Space.ml provides the distance and metric primitives used throughout
    the KPop suite.  It defines `Space.Distance.t` (`euclidean`,
    `cosine`, `angle`, `minkowski(p)`), `Space.Distance.Metric.t` for
    the inertia-weighting transformation (`flat`, `powers(...)`), and
    the pairwise `compute` / `compute_norm` helpers operating on
    `Float.Array.t` row vectors.

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

open BiOCamLib
open Better

(* A number of distance functions.
   They all have signature: float array -> float array -> float array -> float *)
module Distance:
  sig
    module Metric:
      (* Functions to deduce a metric from a vector *)
      sig
        (* We make the type private to implement constraints *)
        type t = private
          | Flat
          (* Parameters are: internal power, threshold (on accumulated transformed inertia), external power.
            Powers must be >= 0., threshold between 0. and 1. *)
          | Powers of float * float * float
        val compute: t -> Float.Array.t -> Float.Array.t
        val of_string: string -> t
        val to_string: t -> string
      end
    (* We make the type private to implement constraints *)
    type t = private
      | Euclidean
      | Cosine (* Same as Euclidean^2 / 2, or 1 - cos theta (in [0, 2]) *)
      | Angle (* Same as theta in radians, i.e. arccos (1 - Euclidean^2 / 2) *)
      | Minkowski of float (* Theoretically speaking, the parameter should be an integer *)
    (* Arguments are distance function, metric, v *)
    val compute_norm: t -> Float.Array.t -> Float.Array.t -> float
    (* Arguments are distance function, metric, a, b.
       We also parameterise the computation with two arbitrary adaptor functions for a and b,
        which allow, for instance, to implement normalised differences *)
    val compute: ?adaptor_a:(float -> float) -> ?adaptor_b:(float -> float) ->
                 t -> Float.Array.t -> Float.Array.t -> Float.Array.t -> float
    val of_string: string -> t
    val to_string: t -> string
    module Iterator:
      sig
        type distance_t = t
        type t (* Mutable *)
        val range: t -> float * float
        val make: ?max_distance_component:float -> distance_t -> float -> (int -> float) -> int -> t
        (* Returns indices and 1D-component of the distance *)
        val get_opt: t -> (int * int * float) option
        val incr: ?max_distance_component:float -> t -> unit (* The operator is mutable *)
        val output_summary: t -> unit
      end
  end
= struct
    (* Metric *)
    module Metric =
      struct
        type t =
          | Flat
          | Powers of float * float * float
        (* Internal type used to cluster identical consecutive elements *)
        module RFreqs = Numbers.RFloatFreqsVector
        let compute = function
          | Flat ->
            (fun m ->
              (* The flat metric leaves distances unchanged *)
              let l = Float.Array.length m in
              if l = 0 then
                m
              else
                Float.Array.make l 1.)
          | Powers (power_int, threshold, power_ext) ->
            (fun m ->
              (* We assume that the elements of the metric are not necessarily normalised
                  but always non-negative, as inertia (which is the square of SVs) would be.
                 In order to guarantee that the order of elements does not change after transformation,
                  (1) elements must also be sorted (in decreasing order, as SVs would be), and
                  (2) powers must be non-negative.
                 Requirement (2) is checked at construction time as the type is private,
                  so we just assert it; same for the threshold, which must be
                  between 0. and 1. .
                 Requirement (1) is immaterial, as of_floatarray() internally re-sorts elements.
                 This leaves only the positivity of elements to be checked, which is done
                  by of_floatarray() as well.
                 Note that under these assumptions the "neutral" transformation
                  leaving the flat metric invariant is "powers(1,1,1)".
                 Finally, we multiply the result by the number of dimensions,
                  for consistency with the case of flat metric *)
              assert (power_int >= 0. && threshold >= 0. && threshold <= 1. && power_ext >= 0.);
              let l = Float.Array.length m |> float_of_int in
              RFreqs.(of_floatarray ~non_negative:true m
                       |> normalize_abs |> pow_abs power_int
                       |> threshold_accum_abs threshold |> pow_abs power_ext
                       |> normalize_abs |> to_floatarray)
                |> Float.Array.map (( *. ) l))
        let of_string_re = Str.regexp "[(,)]"
        let of_string s =
          let raise kind = Exception.raise __FUNCTION__ IO_Format (Printf.sprintf "%s metric '%s'" kind s) in
          match Str.full_split of_string_re s with
          | [ Text "flat" ] ->
            Flat
          | [ Text "powers"; Delim "(";
              Text power_int; Delim ","; Text threshold; Delim ","; Text power_ext;
              Delim ")" ] ->
            (* Powers must be non-negative, and the threshold between 0. and 1. *)
            let power_int, threshold, power_ext =
              try
                float_of_string power_int, float_of_string threshold, float_of_string power_ext
              with _ ->
                raise "Invalid" in
            if power_int < 0. || (threshold < 0. || threshold > 1.) || power_ext < 0. then
              raise "Invalid";
            Powers (power_int, threshold, power_ext)
          | _ ->
            raise "Unknown"
        let to_string = function
          | Flat ->
            "flat"
          | Powers (power_int, threshold, power_ext) ->
            Printf.sprintf "powers(%.15g,%.15g,%.15g)" power_int threshold power_ext
      end
    (* Distance *)
    type t =
      | Euclidean
      | Cosine
      | Angle
      | Minkowski of float
    let compute_norm_unscaled_component f metr diff =
      match f with
      | Euclidean | Cosine | Angle ->
        diff *. diff *. metr
      | Minkowski power ->
        ((abs_float diff) ** power) *. metr
    let ( .@() ) = Float.Array.( .@() )
    let compute_norm f m v =
      let acc = ref 0. in
      match f with
      | Euclidean | Cosine | Angle ->
        Float.Array.iteri
          (fun i el ->
            acc := !acc +. (el *. el *. m.@(i)))
          v;
        sqrt !acc
      | Minkowski power ->
        Float.Array.iteri
          (fun i el ->
            acc := !acc +. (((abs_float el) ** power) *. m.@(i)))
          v;
        !acc ** (1. /. power)
    let compute ?(adaptor_a = Fun.id) ?(adaptor_b = Fun.id) f m a b =
      let length_a = Float.Array.length a and length_m = Float.Array.length m and length_b = Float.Array.length b in
      if length_a <> length_m || length_m <> length_b then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Vectors have incompatible lengths (found [%d, %d, %d] for left, metric, right)"
            length_a length_m length_b);
      (* We could define everything in terms of compute_unscaled_component,
          but we rewrite things in order to optimise the switch out of the cycle *)
      let acc = ref 0. in
      match f with
      | Minkowski power ->
        Float.Array.iter2i
          (fun i el_a el_b ->
            let diff = abs_float (adaptor_a el_a -. adaptor_b el_b) in
            acc := !acc +. ((diff ** power) *. m.@(i)))
          a b;
        !acc ** (1. /. power)
      | _ ->
        Float.Array.iter2i
          (fun i el_a el_b ->
            let diff = adaptor_a el_a -. adaptor_b el_b in
            acc := !acc +. (diff *. diff *. m.@(i)))
          a b;
        begin match f with
        | Euclidean -> sqrt !acc
        | Cosine -> !acc /. 2.
        | Angle -> acos (1. -. !acc /. 2.)
        | _ -> assert false
        end
    let of_string_re = Str.regexp "[()]"
    let of_string s =
      let fail kind =
        Exception.raise __FUNCTION__ Initialize (Printf.sprintf "%s distance '%s'" kind s) in
      match Str.full_split of_string_re s with
      | [ Text "euclidean" ] ->
        Euclidean
      | [ Text "cosine" ] ->
        Cosine
      | [ Text "angle" ] ->
        Angle
      | [ Text "minkowski"; Delim "("; Text power; Delim ")" ] ->
          let power =
            try
              float_of_string power
            with _ ->
              fail "Invalid" in
          if power < 0. then
            fail "Invalid";
          Minkowski power
        | _ ->
          fail "Unknown"
    let to_string = function
      | Euclidean -> "euclidean"
      | Cosine -> "cosine"
      | Angle -> "angle"
      | Minkowski power -> Printf.sprintf "minkowski(%.15g)" power
    module Iterator =
      struct
        module FloatIntMultimap = Tools.Multimap (ComparableFloat) (ComparableInt)
        type distance_t = t
        type t = {
          (* Function of the _difference_ between components of the same vector *)
          compute_distance_component: float -> float;
          n: int;
          sorted: FloatIntMultimap.t;
          (* There is, in the worst case, one minimum per stride *)
          mutable state: state_t IntMap.t;
          (* The invariant is that the minimum for the highest stride is always global.
             If the level is unbroken, one can avoid computing the following one *)
        } and state_t = {
          lo: float * int; (* Value in the multimap *)
          hi: float * int  (* Value in the multimap *)
        }
        let range it =
          if it.n = 0 then
            0., 0.
          else begin
            let coord_lo, _ = FloatIntMultimap.KeyMap.min_binding it.sorted
            and coord_hi, _ = FloatIntMultimap.KeyMap.max_binding it.sorted in
            coord_lo, coord_hi
          end
        let output_summary it =
          Printf.printf "Distance.Iterator( n=%d state={" it.n;
          IntMap.iter
            (fun i {lo = (lo_coord, lo_idx); hi = (hi_coord, hi_idx)} ->
              Printf.printf " %d->[d=%.14g|%d->%.14g|%d->%.14g]"
                i (hi_coord -. lo_coord) lo_idx lo_coord hi_idx hi_coord)
            it.state;
          Printf.printf " } )\n%!"
        (* Auxiliary function *)
        let sorted_to_state lo hi =
          let lo_coord, lo_set = lo and hi_coord, hi_set = hi in
          { lo = lo_coord, FloatIntMultimap.ValSet.min_elt lo_set;
            hi = hi_coord, FloatIntMultimap.ValSet.min_elt hi_set }
        (* Computes the minimum state, larger than a given bound for the difference, across the whole stride.
           Can be None if the minimum difference induces a component above max_distance_component.
           In the special case stride = 0, we return something only if there are coinciding points.
           If stride > 0 and this is None, then all the strides above will be as well,
            as we are considering the minimum over all possible intervals *)
        let get_minimum_opt ?(max_distance_component = infinity) dist sorted stride diff_bound =
          try
            let lo = FloatIntMultimap.KeyMap.min_binding sorted |> ref
            and max_coord, _ = FloatIntMultimap.KeyMap.max_binding sorted and min_state = ref None in
            if stride = 0 then begin
              let process_lo () =
                let lo_coord, lo_set = !lo in
                if FloatIntMultimap.ValSet.cardinal lo_set > 1 then begin
                  let lo_first = FloatIntMultimap.ValSet.min_elt lo_set in
                  let lo_second = FloatIntMultimap.ValSet.find_next lo_first lo_set in
                  min_state := Some { lo = lo_coord, lo_first; hi = lo_coord, lo_second }
                end in
              process_lo ();
              while fst !lo <> max_coord && !min_state = None do
                lo := FloatIntMultimap.KeyMap.find_next (fst !lo) sorted;
                process_lo ()
              done;
              !min_state
            end else begin
              (* First, we generate the first interval *)
              let hi = ref !lo and i = ref stride in
              while !i > 0 do
                hi := FloatIntMultimap.KeyMap.find_next (fst !hi) sorted;
                decr i
              done;
              let min_diff = ref infinity and min_state = ref None in
              let process_interval () =
                let diff = fst !hi -. fst !lo in
                (* We also have to satisfy the lower bound *)
                if diff > diff_bound && diff < !min_diff then begin
                  min_diff := diff;
                  min_state := Some (sorted_to_state !lo !hi)
                end in
              process_interval ();
              (* Second, we iterate over all possible intervals *)
              while fst !hi <> max_coord do
                lo := FloatIntMultimap.KeyMap.find_next (fst !lo) sorted;
                hi := FloatIntMultimap.KeyMap.find_next (fst !hi) sorted;
                process_interval ()
              done;
              if !min_state <> None && dist !min_diff <= max_distance_component then
                !min_state
              else
                None
            end
          with _ ->
            None
        (* Computes the next valid interval within a stride. The current interval is assumed to be valid.
           The interval can be after the current one, but then the distance must be the same.
           It can also be before the current one, but then the distance must be greater *)
        let get_next_opt ?(max_distance_component = infinity) dist sorted stride state =
          let lo_coord, lo_idx = state.lo and hi_coord, hi_idx = state.hi in
          let diff = hi_coord -. lo_coord
          and lo_set = FloatIntMultimap.KeyMap.find lo_coord sorted in
          let max_lo_set = FloatIntMultimap.ValSet.max_elt lo_set
          and max_coord, _ = FloatIntMultimap.KeyMap.max_binding sorted in
          try
            (* First, we try and see if there are other intervals with the same distance.
               They will come _after_ the current one *)
            if stride = 0 then begin
              (* Is there one more pair in this group?
                 The last of the group is defined by lo = next_to_last, hi = last
                  (we exclude pairs made of identical elements).
                 If that has been reached, one has to move to the next group *)
              let lo_idx = ref lo_idx and hi_idx = ref hi_idx in
              while begin
                if !hi_idx = max_lo_set then begin
                  lo_idx := FloatIntMultimap.ValSet.find_next !lo_idx lo_set;
                  if !lo_idx <> max_lo_set then
                    hi_idx := FloatIntMultimap.ValSet.find_next !lo_idx lo_set
                end else
                  hi_idx := FloatIntMultimap.ValSet.find_next !hi_idx lo_set;
                !lo_idx <> max_lo_set && !lo_idx = !hi_idx
              end do
                ()
              done;
              if !lo_idx <> max_lo_set then
                Some { lo = lo_coord, !lo_idx; hi = lo_coord, !hi_idx }
              else begin
                (* Is there another group with more than one element? *)
                let lo = ref (lo_coord, lo_set) in
                while begin
                  lo := FloatIntMultimap.KeyMap.find_next (fst !lo) sorted;
                  fst !lo <> max_coord && snd !lo |> FloatIntMultimap.ValSet.cardinal = 1
                end do
                  ()
                done;
                let lo_coord, lo_set = !lo in
                if lo_coord <> max_coord then begin
                  (* Here cardinal > 1 *)
                  let lo_first = FloatIntMultimap.ValSet.min_elt lo_set in
                  let lo_second = FloatIntMultimap.ValSet.find_next lo_first lo_set in
                  Some { lo = lo_coord, lo_first; hi = lo_coord, lo_second }
                end else
                  None
              end
            end else begin
              let hi_set = FloatIntMultimap.KeyMap.find hi_coord sorted in
              let max_hi_set = FloatIntMultimap.ValSet.max_elt hi_set in
              if lo_idx = max_lo_set && hi_idx = max_hi_set && hi_coord = max_coord then
                raise_notrace Exit; (* This one is OK as it will be caught *)
              if hi_idx <> max_hi_set then
                Some { state with hi = hi_coord, FloatIntMultimap.ValSet.find_next hi_idx hi_set }
              else if lo_idx <> max_lo_set then
                Some
                  { lo = lo_coord, FloatIntMultimap.ValSet.find_next lo_idx lo_set;
                    hi = hi_coord, FloatIntMultimap.ValSet.min_elt hi_set }
              else begin
                (* Here hi_coord <> max_coord *)
                let lo = ref (lo_coord, lo_set) and hi = ref (hi_coord, hi_set) and hi_coord = ref hi_coord in
                while begin
                  lo := FloatIntMultimap.KeyMap.find_next (fst !lo) sorted;
                  hi := FloatIntMultimap.KeyMap.find_next (fst !hi) sorted;
                  hi_coord := fst !hi;
                  !hi_coord <> max_coord && !hi_coord -. fst !lo <> diff
                end do
                  ()
                done;
                (* Remember that for this to be valid, the distance must be exactly the same *)
                if !hi_coord = max_coord then
                  (* We have to restart from the beginning *)
                  raise_notrace Exit (* This one is OK as it will be caught *)
                else begin
                  (* Here !hi_coord -. !lo_coord = diff *)
                  assert (!hi_coord -. fst !lo = diff);
                  Some (sorted_to_state !lo !hi)
                end
              end
            end
          with Exit ->
            (* If that has not worked, we proceed to the next valid interval within the same stride, if any *)
            get_minimum_opt ~max_distance_component dist sorted stride diff
        let make ?(max_distance_component = infinity) dist metr init n =
          let sorted = ref FloatIntMultimap.empty in
          for i = 0 to n - 1 do
            sorted := FloatIntMultimap.add (init i) i !sorted
          done;
          let res = {
            compute_distance_component = compute_norm_unscaled_component dist metr;
            n;
            sorted = !sorted;
            state = IntMap.empty
          } in
          begin match
            get_minimum_opt ~max_distance_component res.compute_distance_component res.sorted 0 neg_infinity
          with
          | None ->
            (* Means there are no coinciding points, we try with 1 *)
            begin match
              get_minimum_opt ~max_distance_component res.compute_distance_component res.sorted 1 neg_infinity
            with
            | None ->
              (* This might only occur when there are no points. Nothing to do, as the iterator is currently invalid *)
              ()
            | Some w ->
              res.state <- IntMap.singleton 1 w (* At the moment, level 1 is complete *)
            end
          | Some w ->
            res.state <- IntMap.singleton 0 w (* At the moment, level 0 is complete *)
          end;
          res
        (* Auxiliary function. It assumes there is a valid interval *)
        let find_minimum_interval it =
          let min_stride = ref it.n and min_diff = ref infinity in
          (* First we find what is the minimum across all active layers *)
          IntMap.iter
            (fun i state ->
              let coord_lo, _ = state.lo and coord_hi, _ = state.hi in
              let diff = coord_hi -. coord_lo in
              if diff < !min_diff then begin
                min_stride := i;
                min_diff := diff
              end)
            it.state;
          let min_stride = !min_stride in
          (* If there are no valid intervals, the iterator must have been previously invalidated *)
          assert (it.state <> IntMap.empty && min_stride <> it.n);
          (* There must be at least one valid interval *)
          min_stride, IntMap.find min_stride it.state
        let get_opt it =
          if it.state = IntMap.empty then
            None
          else begin
            let _, min_state = find_minimum_interval it in
            let coord_lo, idx_lo = min_state.lo and coord_hi, idx_hi = min_state.hi in
            Some (min idx_lo idx_hi, max idx_lo idx_hi, coord_hi -. coord_lo |> it.compute_distance_component)
          end
        let incr ?(max_distance_component = infinity) it =
          (* We find and update the minimum interval.
             If the update happens in the topmost stride, we compute the next one, if any is still available *)
          if it.state <> IntMap.empty then begin
            let min_stride, min_state = find_minimum_interval it in
            begin match begin
              get_next_opt ~max_distance_component it.compute_distance_component it.sorted min_stride min_state
            end with
            | None -> it.state <- IntMap.remove min_stride it.state
            | Some w -> it.state <- IntMap.add min_stride w it.state
            end;
            (* At this point the iterator can legitimately be invalid *)
            if it.state <> IntMap.empty then begin
              let aug_min_stride = min_stride + 1 and stride_hi, _ = IntMap.max_binding it.state in
              if min_stride = stride_hi && aug_min_stride <> it.n then begin
                (* We can use the difference at this level as a bound for the next one, I think :-) *)
                let coord_lo, _ = min_state.lo and coord_hi, _ = min_state.hi in
                let diff = coord_hi -. coord_lo in
                match begin
                  get_minimum_opt ~max_distance_component it.compute_distance_component it.sorted aug_min_stride diff
                end with
                | None -> ()
                | Some w -> it.state <- IntMap.add aug_min_stride w it.state
              end
            end
          end
    end
  end

