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
  end

