(*
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

(* A number of distance functions.
   They all have signature: float array -> float array -> float array -> float *)
module Distance:
  sig
    module Metric:
      (* Functions to deduce a metric from a vector *)
      sig
        type t (* We hide the type to implement constraints *)
        exception Negative_metric_element of float
        exception Unsorted_metric_vector of Float.Array.t
        val compute: t -> Float.Array.t -> Float.Array.t
        exception Unknown_metric of string
        exception Negative_power of float
        exception Negative_tightness of float
        val of_string: string -> t
        val to_string: t -> string
      end
    type t (* We hide the type to implement constraints *)
    exception Invalid_dimension of int
    val compute_component: t -> Float.Array.t -> int -> float -> float
    exception Incompatible_lengths of int * int * int
    (* What happens when the vectors have incompatible lengths *)
    type mode_t =
      | Fail
      | Infinity
    val set_mode: mode_t -> unit
    val compute: t -> Float.Array.t -> Float.Array.t -> Float.Array.t -> float
    exception Unknown_distance of string
    exception Negative_power of float
    val of_string: string -> t
    val to_string: t -> string
    module Iterator:
      sig
        type distance_t = t
        type t (* Mutable *)
        val make: ?max_distance_component:float -> distance_t -> Float.Array.t -> int -> (int -> float) -> int -> t
        (* Returns indices and 1D-component of the distance *)
        val get_opt: t -> (int * int * float) option
        val incr: ?max_distance_component:float -> t -> unit (* The operator is mutable *)
      end
  end
= struct
    module Metric =
      struct
        type t =
          | Flat
          | Power of float
          (* Parameters are: threshold (1/n/t), tightness left, tightness right *)
          | Sigmoid of float * float * float * float
        exception Negative_metric_element of float
        exception Unsorted_metric_vector of Float.Array.t
        let compute = function
          | Flat ->
            Float.Array.map
             (fun el ->
               if el < 0. then
                 Negative_metric_element el |> raise;
               1.)
          | Power power ->
            (fun v ->
              (* We normalise with respect to the transformed max *)
              let m = ref 0. in
              Float.Array.iter
                (fun el ->
                  if el < 0. then
                   Negative_metric_element el |> raise;
                  m := el ** power |> max !m)
                v;
              let m = !m in
              if m > 0. then
                Float.Array.map
                  (fun el ->
                    (el ** power) /. m)
                  v
              else
                (* All zeros *)
                v)
          | Sigmoid (power, threshold, l_tightness, r_tightness) ->
            (fun v ->
              (* We normalise to one *)
              let acc = ref 0. in
              Float.Array.iteri
                (fun i el ->
                  if el < 0. then
                    Negative_metric_element el |> raise;
                  if i > 0 && el > Float.Array.get v (i - 1) then
                    Unsorted_metric_vector v |> raise;
                  acc := !acc +. el ** power)
                v;
              let acc = !acc in
              if acc > 0. then begin
                (* Determine inflection point *)
                let f_n = float_of_int (Float.Array.length v) in
                let threshold = 1. /. threshold /. f_n and t = ref (-1) in
                Float.Array.iteri
                  (fun i el ->
                    let el = (el ** power) /. acc in
                    if el < threshold && !t = (-1) then
                      t := i)
                  v;
                if !t <> (-1) then begin
                  let t = (float_of_int !t +. 0.5) /. f_n in
                  Float.Array.mapi
(*
f<-function(x,t=0.5,kl=10,kr=100){a<-ifelse(x<t,x/t,(x-t)/(1-t)); y<-ifelse(x<t,(-(2*a-3)*a*a-1)+1/2/(1-t)*((a-1)*a*a),1/t/2*((a-2)*a*a+a)-(2*a-3)*a*a); (ifelse(y>0,-y*((1+kr)/kr)/((1+y*kr)/kr),-y*((1+kl)/kl)/((1-y*kl)/kl))+1)/2}
*)
                    (fun i _ ->
                      let x = (float_of_int i +. 0.5) /. f_n in
                      let y =
                        if x < t then begin
                          let a = x /. t in
                          let a_a = a *. a in
                          ((-2. *. a +. 3.) *. a_a -. 1.) +. 0.5 /. (1. -. t) *. ((a -. 1.) *. a_a)
                        end else begin
                          let a = (x -. t) /. (1. -. t) in
                          let a_a = a *. a in
                          0.5 /. t *. ((a -. 2.) *. a_a +. a) -. (2. *. a -. 3.) *. a_a
                        end in
                      let res =
                        if y > 0. then
                          -. y *. ((1. +. r_tightness) /. r_tightness) /. ((1. +. y *. r_tightness) /. r_tightness)
                        else
                          -. y *. ((1. +. l_tightness) /. l_tightness) /. ((1. -. y *. l_tightness) /. l_tightness) in
                      (1. +. res) /. 2. |> max 0. |> min 1.)
                    v
                end else
                  v
              end else
                (* All zeros *)
                v)
        exception Unknown_metric of string
        exception Negative_power of float
        exception Negative_tightness of float
        let of_string_re = Str.regexp "[(,)]"
        let of_string name =
          match name with
          | "flat" ->
            Flat
          | s ->
            match Str.full_split of_string_re s with
            | [ Text "power"; Delim "("; Text power; Delim ")" ] ->
              let power =
                try
                  float_of_string power
                with _ ->
                  Unknown_metric s |> raise in
              if power < 0. then
                Negative_power power |> raise;
              Power power
            | [ Text "sigmoid"; Delim "(";
                Text power; Delim ","; Text thresh; Delim ","; Text l_tight; Delim ","; Text r_tight;
                Delim ")" ] ->
              let power, thresh, l_tight, r_tight =
                try
                  float_of_string power, float_of_string thresh, float_of_string l_tight, float_of_string r_tight
                with _ ->
                  Unknown_metric s |> raise in
              if power < 0. then
                Negative_power power |> raise;
              Sigmoid (power, thresh, l_tight, r_tight)
            | _ ->
              Unknown_metric s |> raise
        let to_string = function
          | Flat ->
            "flat"
          | Power power ->
            Printf.sprintf "power(%.15g)" power
          | Sigmoid (power, thresh, l_tight, r_tight) ->
            Printf.sprintf "sigmoid(%.15g,%.15g,%.15g,%.15g)" power thresh l_tight r_tight
      end
    type t =
      | Euclidean
      | Minkowski of float (* Theoretically speaking, the parameter should be an integer *)
    exception Incompatible_lengths of int * int * int
    type mode_t =
      | Fail
      | Infinity
    let mode = ref Fail
    let set_mode new_mode = mode := new_mode
    let compute f m a b =
      let length_a = Float.Array.length a and length_m = Float.Array.length m and length_b = Float.Array.length b in
      if length_a <> length_m || length_m <> length_b then begin
        match !mode with
        | Fail -> Incompatible_lengths (length_a, length_m, length_b) |> raise
        | Infinity -> infinity
      end else
        match f with
        | Euclidean ->
          let acc = ref 0. in
          Float.Array.iteri
            (fun i el_a ->
              let diff = el_a -. Float.Array.get b i in
              acc := !acc +. (diff *. diff *. Float.Array.get m i))
            a;
          sqrt !acc
        | Minkowski power ->
          let acc = ref 0. in
          Float.Array.iteri
            (fun i el_a ->
              let diff = el_a -. Float.Array.get b i |> abs_float in
              acc := !acc +. ((diff ** power) *. Float.Array.get m i))
            a;
          !acc ** (1. /. power)
    exception Invalid_dimension of int
    let compute_component f m dim diff =
      let length_m = Float.Array.length m in
      if dim >= length_m then
        Invalid_dimension dim |> raise;
      match f with
      | Euclidean ->
        diff *. diff *. Float.Array.get m dim
      | Minkowski power ->
        (diff ** power) *. Float.Array.get m dim
    exception Unknown_distance of string
    exception Negative_power of float
    let of_string_re = Str.regexp "[()]"
    let of_string = function
      | "euclidean" ->
        Euclidean
      | s ->
        match Str.full_split of_string_re s with
        | [ Text "minkowski"; Delim "("; Text power; Delim ")" ] ->
          let power =
            try
              float_of_string power
            with _ ->
              Unknown_distance s |> raise in
          if power < 0. then
            Negative_power power |> raise;
          Minkowski power
        | _ ->
          Unknown_distance s |> raise
    let to_string = function
      | Euclidean -> "euclidean"
      | Minkowski power -> Printf.sprintf "minkowski(%.15g)" power
    module Iterator =
      struct
        module FloatIntMultimap = Tools.Multimap (Tools.ComparableFloat) (Tools.ComparableInt)
        type distance_t = t
        type t = {
          (* Function of the _difference_ between components *)
          compute_distance_component: float -> float;
          n: int;
          sorted: FloatIntMultimap.t;
          (* There is, in the worst case, one minimum per stride *)
          state: state_t option array;
          mutable stride_lo: int;
          (* The invariant is that the minimum for the highest stride is always global.
             If the level is unbroken, one can avoid computing the following one *)
          mutable stride_hi: int
        } and state_t = {
          lo: float * int; (* Value in the multimap *)
          hi: float * int  (* Value in the multimap *)
        }
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
            and max_coord, _ = FloatIntMultimap.max_binding sorted and min_state = ref None in
            if stride = 0 then begin
              let process_lo () =
                if snd !lo |> FloatIntMultimap.ValSet.cardinal > 1 then
                  min_state := Some (sorted_to_state !lo !lo) in
              process_lo ();
              while fst !lo <> max_coord && !min_state = None do
                lo := FloatIntMultimap.find_next (fst !lo) sorted;
                process_lo ()
              done;
              !min_state
            end else begin
              (* First, we generate the first interval *)
              let hi = ref !lo and i = ref stride in
              while !i > 0 do
                hi := FloatIntMultimap.find_next (fst !hi) sorted;
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
                lo := FloatIntMultimap.find_next (fst !lo) sorted;
                hi := FloatIntMultimap.find_next (fst !hi) sorted;
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
          and lo_set = FloatIntMultimap.KeyMap.find lo_coord sorted
          and hi_set = FloatIntMultimap.KeyMap.find hi_coord sorted in
          let max_lo_set = FloatIntMultimap.ValSet.max_elt lo_set
          and max_hi_set = FloatIntMultimap.ValSet.max_elt hi_set
          and max_coord, _ = FloatIntMultimap.max_binding sorted in
          try
            (* First, we try and see if there are other intervals with the same distance.
               They will come _after_ the current one *)
            if lo_idx = max_lo_set && hi_idx = max_hi_set && hi_coord = max_coord then
              raise_notrace Exit;
            if hi_idx <> max_hi_set then
              Some
                { state with
                  hi =
                    hi_coord,
                    FloatIntMultimap.ValSet.find_first (fun k -> FloatIntMultimap.ValOrd.compare k hi_idx > 0) hi_set }
            else if lo_idx <> max_lo_set then
              Some
                { state with
                  lo =
                    lo_coord,
                    FloatIntMultimap.ValSet.find_first (fun k -> FloatIntMultimap.ValOrd.compare k lo_idx > 0) lo_set }
            else begin
              (* Here hi_coord <> max_coord *)
              let lo = ref (lo_coord, lo_set) and hi = ref (hi_coord, hi_set) and hi_coord = ref hi_coord in
              while begin
                lo := FloatIntMultimap.find_next (fst !lo) sorted;
                hi := FloatIntMultimap.find_next (fst !hi) sorted;
                hi_coord := fst !hi;
                !hi_coord <> max_coord && !hi_coord -. fst !lo <> diff
              end do
                ()
              done;
              (* Remember that for this to be valid, the distance must be exactly the same *)
              if !hi_coord = max_coord then
                (* We have to restart from the beginning *)
                raise_notrace Exit
              else begin
                (* Here !hi_coord -. !lo_coord = diff *)
                assert (!hi_coord -. fst !lo = diff);
                Some (sorted_to_state !lo !hi)
              end
            end
          with Exit ->
            (* If that has not worked, we proceed to the next valid interval within the same stride, if any *)
            get_minimum_opt ~max_distance_component dist sorted stride diff
        let make ?(max_distance_component = infinity) dist metr d init n =
          let red_n = n - 1 and sorted = ref FloatIntMultimap.empty in
          for i = 0 to red_n do
            sorted := FloatIntMultimap.add (init i) i !sorted
          done;
          let res = {
            compute_distance_component = compute_component dist metr d;
            n;
            sorted = !sorted;
            state = Array.make n None;
            stride_lo = 0;
            stride_hi = 0
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
              (* This might only occur when there are no points. We invalidate the iterator *)
              res.stride_lo <- n;
              res.stride_hi <- n - 1
            | w ->
              res.state.(1) <- w;
              res.stride_lo <- 1;
              res.stride_hi <- 1 (* At the moment, level 1 is complete *)
            end
          | w ->
            res.state.(0) <- w;
            res.stride_lo <- 0;
            res.stride_hi <- 0 (* At the moment, level 0 is complete *)
          end;
          res
        (* Auxiliary function. It assumes there is a valid interval *)
        let find_minimum_interval it =
          let min_stride = ref it.n and min_diff = ref infinity in
          (* First we find what is the minimum across all active layers *)
          for i = it.stride_lo to it.stride_hi do
            match it.state.(i) with
            | None -> ()
            | Some state ->
              let coord_lo, _ = state.lo and coord_hi, _ = state.hi in
              let diff = coord_hi -. coord_lo in
              if diff < !min_diff then begin
                min_stride := i;
                min_diff := diff
              end
          done;
          let min_stride = !min_stride in
          (* If there are no valid intervals, the iterator must have been previously invalidated *)
          assert (min_stride <> it.n);
          (* There must be at least one valid interval *)
          match it.state.(min_stride) with
          | None ->
            assert false
          | Some min_state ->
            min_stride, min_state
        let get_opt it =
          if it.stride_lo = it.n then
            None
          else begin
            let _, min_state = find_minimum_interval it in
            let coord_lo, idx_lo = min_state.lo and coord_hi, idx_hi = min_state.hi in
            Some (min idx_lo idx_hi, max idx_lo idx_hi, coord_hi -. coord_lo |> it.compute_distance_component)
          end
        let incr ?(max_distance_component = infinity) it =
          (* We find and update the minimum interval.
             If the update happens in the topmost stride, we compute the next one, if any is still available *)
          if it.stride_lo <> it.n then begin
            let min_stride, min_state = find_minimum_interval it in
            it.state.(min_stride) <-
              get_next_opt ~max_distance_component it.compute_distance_component it.sorted min_stride min_state;
            let aug_min_stride = min_stride + 1 in
            if min_stride = it.stride_hi && aug_min_stride <> it.n then begin
              (* We can use the difference at this level as a bound for the next one, I hope :-) *)
              let coord_lo, _ = min_state.lo and coord_hi, _ = min_state.hi in
              let diff = coord_hi -. coord_lo in
              it.state.(aug_min_stride) <-
                get_minimum_opt ~max_distance_component it.compute_distance_component it.sorted aug_min_stride diff;
              it.stride_hi <- aug_min_stride
            end;
            (* Should we invalidate the iterator? *)
            let stride = ref 0 in
            while !stride < it.n && it.state.(!stride) = None do
              incr stride
            done;
            if !stride = it.n then begin
              it.stride_lo <- it.n;
              it.stride_hi <- it.n - 1
            end else begin
              (* Case it.state.(!stride) <> None *)
              it.stride_lo <- !stride;
              while !stride < it.n do
                if it.state.(!stride) <> None then
                  it.stride_hi <- !stride;
                incr stride
              done
            end
          end
    end
  end

