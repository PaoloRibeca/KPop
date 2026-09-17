(* Distances between the rows of an embedding, as KPop's Space.Distance computes them, built rung by
   rung.  A metric weight w multiplies the squared component difference under Euclidean and angle, and
   the absolute one under Manhattan, Minkowski(1); angle is the angle between the two vectors after
   each is normalised under the same weights, acos (1 - |a/|a| - b/|b||^2 / 2).  The flat metric
   weighs every axis 1, powers(1,1,1) weighs an axis by its inertia -- a common factor, which the
   library applies and this does not, changes neither a silhouette nor which pairs fall below a
   valley. *)

type kind = Euclidean | Manhattan | Angle

let kind_of_string = function
  | "euclidean" -> Euclidean
  | "manhattan" -> Manhattan
  | "angle" -> Angle
  | s -> failwith ("unknown distance " ^ s)

let kind_to_string = function Euclidean -> "euclidean" | Manhattan -> "manhattan" | Angle -> "angle"

(* Coordinates scaled so that plain differences carry the metric: by the square root of the weight for
   the squared differences of Euclidean and angle, by the weight itself for Manhattan *)
let weigh kind ~flat iv x a =
  if flat then x else match kind with Euclidean | Angle -> x *. sqrt iv.(a) | Manhattan -> x *. iv.(a)

(* Adds axes [from_, to_) to the pair accumulations -- squared or absolute differences, or dot
   products -- and, for angle, to the rows' squared norms *)
let accumulate kind (pc: float array array) acc norms from_ to_ =
  let n = Array.length pc in
  let col = Array.make n 0. in
  for a = from_ to to_ - 1 do
    for i = 0 to n - 1 do col.(i) <- pc.(i).(a) done;
    if kind = Angle then for i = 0 to n - 1 do norms.(i) <- norms.(i) +. (col.(i) *. col.(i)) done;
    let p = ref 0 in
    for i = 0 to n - 2 do
      let x = col.(i) in
      match kind with
      | Euclidean ->
        for j = i + 1 to n - 1 do
          let y = x -. col.(j) in
          Float.Array.unsafe_set acc !p (Float.Array.unsafe_get acc !p +. (y *. y));
          incr p
        done
      | Manhattan ->
        for j = i + 1 to n - 1 do
          Float.Array.unsafe_set acc !p (Float.Array.unsafe_get acc !p +. abs_float (x -. col.(j)));
          incr p
        done
      | Angle ->
        for j = i + 1 to n - 1 do
          Float.Array.unsafe_set acc !p (Float.Array.unsafe_get acc !p +. (x *. col.(j)));
          incr p
        done
    done
  done

let angle_of dot ni nj =
  let den = sqrt (ni *. nj) in
  if den > 0. then acos (Float.max (-1.) (Float.min 1. (dot /. den))) else 0.

(* The distances the accumulations stand for, pair by pair in (i, j > i) order *)
let to_dist kind n acc norms dist =
  let p = ref 0 in
  for i = 0 to n - 2 do
    match kind with
    | Euclidean ->
      for _ = i + 1 to n - 1 do
        Float.Array.unsafe_set dist !p (sqrt (Float.Array.unsafe_get acc !p));
        incr p
      done
    | Manhattan ->
      for _ = i + 1 to n - 1 do
        Float.Array.unsafe_set dist !p (Float.Array.unsafe_get acc !p);
        incr p
      done
    | Angle ->
      for j = i + 1 to n - 1 do
        Float.Array.unsafe_set dist !p (angle_of (Float.Array.unsafe_get acc !p) norms.(i) norms.(j));
        incr p
      done
  done

(* The distance between two rows in their first d axes, stopping early once it exceeds [bound]
   where the accumulation only grows; [norms] are the rows' squared norms in those d axes *)
let between kind (pc: float array array) norms d i l bound =
  let a = pc.(i) and b = pc.(l) in
  match kind with
  | Euclidean ->
    let b2 = bound *. bound and s = ref 0. and k = ref 0 in
    while !k < d && !s <= b2 do
      let y = a.(!k) -. b.(!k) in
      s := !s +. (y *. y);
      incr k
    done;
    if !s <= b2 then sqrt !s else infinity
  | Manhattan ->
    let s = ref 0. and k = ref 0 in
    while !k < d && !s <= bound do
      s := !s +. abs_float (a.(!k) -. b.(!k));
      incr k
    done;
    if !s <= bound then !s else infinity
  | Angle ->
    let s = ref 0. in
    for k = 0 to d - 1 do s := !s +. (a.(k) *. b.(k)) done;
    let t = angle_of !s norms.(i) norms.(l) in
    if t <= bound then t else infinity

(* Leader clustering at a radius in the first d axes, in index order *)
let leader kind pc norms d radius =
  let n = Array.length pc in
  let assign = Array.make n (-1) and leaders = ref [] and nl = ref 0 in
  for i = 0 to n - 1 do
    let best = ref (-1) and best_d = ref infinity in
    List.iter
      (fun l ->
        let t = between kind pc norms d i l radius in
        if t < !best_d then begin
          best_d := t;
          best := l
        end)
      !leaders;
    if !best >= 0 then assign.(i) <- assign.(!best)
    else begin
      leaders := i :: !leaders;
      assign.(i) <- !nl;
      incr nl
    end
  done;
  assign, !nl

(* A radius no pair of rows exceeds in the first d axes *)
let top_radius kind (pc: float array array) d =
  match kind with
  | Angle -> Float.pi *. 1.000001
  | Euclidean | Manhattan ->
    let s = ref 0. in
    for a = 0 to d - 1 do
      let lo = ref infinity and hi = ref neg_infinity in
      Array.iter (fun row -> let v = row.(a) in if v < !lo then lo := v; if v > !hi then hi := v) pc;
      let r = !hi -. !lo in
      s := !s +. (if kind = Euclidean then r *. r else r)
    done;
    (if kind = Euclidean then sqrt !s else !s) *. 1.000001

(* The classical silhouette of a partition over the distances given, a singleton scoring 0 *)
let silhouette n dist assign k =
  let s = Array.make (n * k) 0. and cnt = Array.make k 0 in
  Array.iter (fun c -> cnt.(c) <- cnt.(c) + 1) assign;
  let p = ref 0 in
  for i = 0 to n - 2 do
    let ci = assign.(i) and base_i = i * k in
    for j = i + 1 to n - 1 do
      let dd = Float.Array.unsafe_get dist !p in
      incr p;
      let cj = assign.(j) in
      s.(base_i + cj) <- s.(base_i + cj) +. dd;
      s.((j * k) + ci) <- s.((j * k) + ci) +. dd
    done
  done;
  let tot = ref 0. in
  for i = 0 to n - 1 do
    let ci = assign.(i) in
    if cnt.(ci) > 1 then begin
      let a = s.((i * k) + ci) /. float_of_int (cnt.(ci) - 1) and b = ref infinity in
      for c = 0 to k - 1 do
        if c <> ci && cnt.(c) > 0 then begin
          let m = s.((i * k) + c) /. float_of_int cnt.(c) in
          if m < !b then b := m
        end
      done;
      let den = Float.max a !b in
      if !b < infinity && den > 0. then tot := !tot +. ((!b -. a) /. den)
    end
  done;
  !tot /. float_of_int n
