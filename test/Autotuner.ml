(*
    Autotuner.ml -- (c) 2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Autotuner.ml exercises what KPop-autotuner builds its choice of space
    and level from, starting with the calibrated valley detector: the
    simplification of a histogram's extrema against known series, the
    thresholds that keep a valley, and the detector end to end on point
    sets whose valleys are known in advance.

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

open BiOCamLib
open Better
open KPop

(* Helpers. *)

(* Under the flat metric the weights are ones whatever the inertia, so the detector measures the
   coordinates exactly as given *)
let flat = Space.Distance.Metric.of_string "flat"
and euclidean = Space.Distance.of_string "euclidean"

let valleys ?seed ?pairs_max ?pairs_sample ?resamples ?multiplicities coords =
  let d = if Array.length coords > 0 then Float.Array.length coords.(0) else 1 in
  Clustering.find_valleys ?seed ?pairs_max ?pairs_sample ?resamples ?multiplicities
    ~metric:flat ~distance:euclidean ~distance_normalize:false coords (Float.Array.make d 1.)

let check_near name ~tolerance ~expected got =
  Testing.verify name (fun () ->
    Float.abs (got -. expected) <= tolerance,
    Printf.sprintf "expected %g within %g, got %.17g" expected tolerance got)

(* A standard normal deviate by Box and Muller, the two uniform draws taken in sequence *)
let gaussian state =
  let u = 1. -. Random.State.float state 1. in
  let v = Random.State.float state 1. in
  sqrt (-2. *. log u) *. cos (2. *. Float.pi *. v)

let blob state ~centre ~sd n =
  Array.init n
    (fun _ ->
      Float.Array.init (Array.length centre) (fun k -> centre.(k) +. (sd *. gaussian state)))

let choose2 n = float_of_int (n * (n - 1) / 2)

(* A series of runs of equal counts, long enough that one bin of smoothing either side leaves every
   run's interior exactly as it was *)
let runs spec = List.concat_map (fun (len, v) -> List.init len (Fun.const v)) spec |> Float.Array.of_list

(* The troughs of a series, checked bin for bin and, within the precision the reference was printed
   with, in prominence and depth *)
let check_troughs name counts expected =
  let got = Clustering.troughs_of_counts counts in
  Testing.check_int (name ^ ": number of troughs") ~expected:(List.length expected) (List.length got);
  if List.length got = List.length expected then
    List.iter2
      (fun (tb, lb, rb, prominence, depth) (etb, elb, erb, eprominence, edepth) ->
        let what = Printf.sprintf "%s: the trough expected at bin %d" name etb in
        Testing.check_equal (what ^ " and its peaks")
          ~to_string:(fun (t, l, r) -> Printf.sprintf "%d between %d and %d" t l r)
          ~expected:(etb, elb, erb) (tb, lb, rb);
        check_near (what ^ ", prominence") ~tolerance:0.05 ~expected:eprominence prominence;
        check_near (what ^ ", relative depth") ~tolerance:0.005 ~expected:edepth depth)
      got expected

(* Troughs. *)

let test_troughs () =
  Testing.section "Troughs of a histogram" (fun () ->
    (* The synthetic cases of the scale-rule prototype, test/ScaleRule/detect3.ml --synthetic, with
       the values it printed *)
    let g c h s b = h *. exp (-.((float_of_int b -. c) ** 2.) /. (2. *. s *. s)) in
    let base b = 20. +. g 40. 1000. 10. b +. g 140. 2000. 15. b in
    let bump h b = base b +. if b >= 78 && b <= 82 then h else 0. in
    let series f = Float.Array.init 200 f in
    check_troughs "two modes" (series base) [ 80, 40, 140, 995.6, 0.98 ];
    List.iter
      (fun h ->
        check_troughs (Printf.sprintf "a bump of +%g inside the gap" h) (series (bump h))
          [ 76, 40, 140, 994.9, 0.98 ])
      [ 5.; 8.; 10.; 30.; 100. ];
    check_troughs "a third mode of 100 inside the gap"
      (series (fun b -> base b +. g 85. 100. 5. b)) [ 72, 40, 140, 986.8, 0.97 ];
    check_troughs "a third mode of 400 inside the gap"
      (series (fun b -> base b +. g 85. 400. 5. b))
      [ 70, 40, 85, 380.8, 0.91; 98, 85, 140, 342.5, 0.82 ];
    check_troughs "a ripple on every bin"
      (series (fun b -> base b +. (3. *. sin (float_of_int b *. 2.)))) [ 81, 40, 140, 995.5, 0.98 ];
    check_troughs "a dip at the top of a mode"
      (series (fun b -> base b -. if b >= 39 && b <= 41 then 100. else 0.))
      [ 80, 43, 140, 952.1, 0.98 ];
    check_troughs "a flat floor" (series (fun b -> if b >= 70 && b <= 100 then 25. else base b))
      [ 85, 40, 140, 991.7, 0.98 ];
    check_troughs "a minimum at the edge of the range"
      (series (fun b -> 20. +. g 140. 2000. 15. b)) [];
    (* Two equal peaks around the shallowest trough: the right one goes with it, so the trough that
       survives is bounded on its left by the first peak and not the second *)
    check_troughs "equal peaks either side of the shallowest trough"
      (runs [ 3, 0.; 10, 100.; 10, 80.; 10, 100.; 10, 0.; 10, 100.; 3, 0. ])
      [ 37, 7, 47, 100., 1. ];
    (* A floor four bins wide after smoothing: the trough is the left of its two middle bins *)
    check_troughs "an even floor" (runs [ 3, 0.; 10, 100.; 6, 0.; 10, 100.; 3, 0. ])
      [ 15, 7, 23, 100., 1. ])

(* Thresholds. *)

let test_keep () =
  Testing.section "Which candidates are valleys" (fun () ->
    let v ?(z = 10.) share = { Clustering.radius = 1.; share; z; reappearance = 1.; depth = 1. } in
    let shares l = Clustering.keep_valleys l |> List.map (fun v -> v.Clustering.share) in
    let to_string l = List.map (Printf.sprintf "%g") l |> String.concat ", " in
    Testing.check_equal "the share window is inclusive at both ends" ~to_string
      ~expected:[ 0.002; 0.9 ] (shares [ v 0.0019; v 0.002; v 0.9; v 0.9001 ]);
    Testing.check_equal "z below 3 is dropped and z of exactly 3 kept" ~to_string
      ~expected:[ 0.5 ] (shares [ v ~z:2.999 0.1; v ~z:3. 0.5 ]);
    Testing.check_equal "a chain of close shares is one valley, the lower share on a tie" ~to_string
      ~expected:[ 0.1 ] (shares [ v 0.12; v 0.1; v 0.11 ]);
    Testing.check_equal "a chain keeps its highest z" ~to_string ~expected:[ 0.12 ]
      (shares [ v 0.1; v 0.11; v ~z:20. 0.12 ]);
    Testing.check_equal "shares 1.5 points apart are two valleys" ~to_string
      ~expected:[ 0.5; 0.515 ] (shares [ v 0.515; v 0.5 ]))

(* The detector end to end. *)

let test_detector () =
  Testing.section "The calibrated detector" (fun () ->
    let state = Random.State.make [| 42 |] in
    let two_blobs =
      Array.append
        (blob state ~centre:[| 0.; 0.; 0.; 0.; 0. |] ~sd:0.05 40)
        (blob state ~centre:[| 10.; 0.; 0.; 0.; 0. |] ~sd:0.05 30) in
    (match valleys two_blobs with
     | [ found ] ->
       check_near "two separated groups: the share below is the share of pairs within a group"
         ~tolerance:1e-3 ~expected:((choose2 40 +. choose2 30) /. choose2 70) found.Clustering.share
     | l -> Testing.check_int "two separated groups give exactly one valley" ~expected:1 (List.length l));
    Testing.check_int "one Gaussian cloud gives no valley" ~expected:0
      (valleys (blob state ~centre:[| 0.; 0.; 0.; 0.; 0. |] ~sd:1. 400) |> List.length);
    Testing.check_int "fewer than 50 points give no valley" ~expected:0
      (valleys (Array.sub two_blobs 0 49) |> List.length);
    Testing.check_int "identical points give no valley" ~expected:0
      (valleys (Array.init 60 (fun _ -> Float.Array.make 5 1.)) |> List.length);
    Testing.check "the same seed gives the same valleys" (fun () ->
      valleys ~seed:7 two_blobs = valleys ~seed:7 two_blobs);
    Testing.check_raises "a single resample is refused, having no standard deviation" (fun () ->
      valleys ~resamples:1 two_blobs);
    (* Resamples that weight every point alike reproduce the full data exactly once rescaled, so
       the prominence does not vary at all *)
    (match
       valleys ~multiplicities:[| Array.make 70 1; Array.make 70 2 |] two_blobs
     with
     | [ found ] ->
       Testing.check_bool "resamples weighting every point alike give an infinite z" ~expected:true
         (found.Clustering.z = infinity);
       Testing.check_float "and a valley present in every one of them" ~expected:1.
         found.Clustering.reappearance
     | l -> Testing.check_int "identical resamples give exactly one valley" ~expected:1 (List.length l));
    (* Above the pair cap the histogram is built from a sample of pairs, which must find the same
       structure *)
    let side = 10. in
    let three_blobs =
      Array.concat
        [ blob state ~centre:[| 0.; 0.; 0. |] ~sd:0.05 100;
          blob state ~centre:[| side; 0.; 0. |] ~sd:0.05 100;
          blob state ~centre:[| side /. 2.; side *. sqrt 3. /. 2.; 0. |] ~sd:0.05 100 ] in
    let all_pairs = valleys three_blobs
    and sampled = valleys ~pairs_max:40_000 ~pairs_sample:20_000 three_blobs in
    Testing.check_int "a sample of pairs finds as many valleys as all of them"
      ~expected:(List.length all_pairs) (List.length sampled);
    List.iter2
      (fun a s ->
        check_near "and at nearly the same share" ~tolerance:0.01 ~expected:a.Clustering.share
          s.Clustering.share)
      all_pairs
      (if List.length sampled = List.length all_pairs then sampled else all_pairs))

(* The axis ladder. *)

let test_ladder () =
  Testing.section "The rungs of the ladder and the one it picks" (fun () ->
    let rungs a = Array.to_list a |> List.map string_of_int |> String.concat " " in
    Testing.check_equal "rungs double from one and end at the most axes there are"
      ~to_string:rungs ~expected:[| 1; 2; 4; 8; 16; 20 |] (Clustering.ladder_rungs 20);
    Testing.check_equal "and end once where that is itself a power of two" ~to_string:rungs
      ~expected:[| 1; 2; 4; 8 |] (Clustering.ladder_rungs 8);
    Testing.check_equal "one axis is a ladder of one rung" ~to_string:rungs ~expected:[| 1 |]
      (Clustering.ladder_rungs 1);
    Testing.check_raises "no axes at all is refused" (fun () -> Clustering.ladder_rungs 0);
    let picked = function Some i -> Printf.sprintf "rung %d" i | None -> "none" in
    (* The right end of the FIRST plateau at the highest score, never a later rung as high *)
    Testing.check_equal "a later, unconnected rung of the same height does not win"
      ~to_string:picked ~expected:(Some 1) (Clustering.pick_rung [| 1; 3; 1; 3 |]);
    Testing.check_equal "a plateau is taken at its right end" ~to_string:picked ~expected:(Some 2)
      (Clustering.pick_rung [| 1; 3; 3; 1; 3 |]);
    Testing.check_equal "1 1 2 3 3 2 2 3 2 picks the second of the first two threes"
      ~to_string:picked ~expected:(Some 4) (Clustering.pick_rung [| 1; 1; 2; 3; 3; 2; 2; 3; 2 |]);
    Testing.check_equal "2 2 3 1 1 1 3 3 2 picks the first three" ~to_string:picked
      ~expected:(Some 2) (Clustering.pick_rung [| 2; 2; 3; 1; 1; 1; 3; 3; 2 |]);
    Testing.check_equal "no rung scoring is no pick" ~to_string:picked ~expected:None
      (Clustering.pick_rung [| 0; 0; 0 |]));
  Testing.section "Which valley is the level" (fun () ->
    let module L = Clustering.Level in
    Testing.check_equal "a level prints as it parses" ~to_string:(String.concat ", ")
      ~expected:[ "finest"; "coarsest"; "share(0.25)" ]
      (List.map (fun s -> L.to_string (L.of_string s)) [ "finest"; "coarsest"; "share(0.25)" ]);
    Testing.check_raises "a share above 1 is refused" (fun () -> L.of_string "share(1.5)");
    Testing.check_raises "and one that is not a number" (fun () -> L.of_string "share(x)");
    Testing.check_raises "and a level of no known kind" (fun () -> L.of_string "median");
    let v share = { Clustering.radius = 1.; share; z = 10.; reappearance = 1.; depth = 1. } in
    let level l valleys = Option.map (fun v -> v.Clustering.share) (L.select l valleys)
    and shown = function Some s -> Printf.sprintf "%g" s | None -> "none" in
    let three = [ v 0.75; v 0.5; v 0.25 ] in
    Testing.check_equal "the finest is the lowest share" ~to_string:shown ~expected:(Some 0.25)
      (level L.Finest three);
    Testing.check_equal "the coarsest stops at half of the pairs" ~to_string:shown
      ~expected:(Some 0.5) (level L.Coarsest three);
    Testing.check_equal "a share midway between two valleys takes the lower" ~to_string:shown
      ~expected:(Some 0.25) (level (L.Share 0.375) three);
    Testing.check_equal "and otherwise the nearer" ~to_string:shown ~expected:(Some 0.5)
      (level (L.Share 0.45) three);
    Testing.check_equal "no valley at or below half of the pairs is no level" ~to_string:shown
      ~expected:None (level L.Finest [ v 0.75 ]));
  Testing.section "The ladder end to end" (fun () ->
    let state = Random.State.make [| 43 |] and dmax = 8 in
    (* Three groups apart on the first three axes, with an inertia that falls away as a real
       analysis's does *)
    let inertia = Float.Array.init dmax (fun k -> 2. ** float_of_int (-k))
    and centre k = Array.init dmax (fun j -> if j = k then 3. else 0.) in
    let coords =
      Array.concat
        [ blob state ~centre:(centre 0) ~sd:0.3 60; blob state ~centre:(centre 1) ~sd:0.3 60;
          blob state ~centre:(centre 2) ~sd:0.3 60 ] in
    (* Every rung must find exactly what the detector finds on that many axes, whichever way the
       ladder takes its distances *)
    let against_detector what ~metric ~distance_normalize =
      let l =
        Clustering.ladder ~dmax ~metric ~distance:euclidean ~distance_normalize coords inertia in
      Array.iteri
        (fun r d ->
          let direct =
            Clustering.calibrated_troughs ~metric ~distance:euclidean ~distance_normalize coords
              (Float.Array.sub inertia 0 d)
            |> Clustering.keep_valleys
          and got = l.Clustering.valleys.(r) in
          let at = Printf.sprintf "%s, %d axes" what d in
          Testing.check_int (at ^ ": the valleys the detector finds there")
            ~expected:(List.length direct) (List.length got);
          if List.length got = List.length direct then
            List.iter2
              (fun a b ->
                check_near (at ^ ": the share") ~tolerance:1e-3 ~expected:a.Clustering.share
                  b.Clustering.share;
                check_near (at ^ ": the radius") ~tolerance:(1e-9 *. a.Clustering.radius)
                  ~expected:a.Clustering.radius b.Clustering.radius)
              direct got)
        l.Clustering.rungs;
      l in
    let flat_ladder = against_detector "running sums, flat" ~metric:flat ~distance_normalize:false in
    ignore
      (against_detector "running sums, powers(1,1,1)"
         ~metric:(Space.Distance.Metric.of_string "powers(1,1,1)") ~distance_normalize:false);
    ignore (against_detector "embeddings afresh" ~metric:flat ~distance_normalize:true);
    (* The level of three groups is the share of pairs inside one, and the rung picked must hold
       it.  Not every rung does: where the groups lie far enough apart that only a few stray pairs
       fall between the two modes, the gap breaks into troughs none of which is significant, and
       the valley is kept only at the fewer axes where the modes are closer *)
    let expected = 3. *. choose2 60 /. choose2 180 in
    (match flat_ladder.Clustering.pick with
     | None -> Testing.check_bool "three groups give the ladder a pick" ~expected:true false
     | Some i ->
       Clustering.Level.select Clustering.Level.Finest flat_ladder.Clustering.valleys.(i)
       |> Option.fold ~none:nan ~some:(fun v -> v.Clustering.share)
       |> check_near "the rung the ladder picks holds the level of the groups" ~tolerance:0.02
            ~expected);
    Testing.check_int "and every resample's choice is recorded" ~expected:50
      (Array.length flat_ladder.Clustering.resample_picks);
    Testing.check_raises "a ladder taller than the inertia is refused" (fun () ->
      Clustering.ladder ~dmax:9 ~metric:flat ~distance:euclidean ~distance_normalize:false coords
        inertia))

let () =
  test_troughs ();
  test_keep ();
  test_detector ();
  test_ladder ();
  Testing.summary ()
