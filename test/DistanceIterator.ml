open BiOCamLib
open KPop.Space

let () =
  let distance = Distance.of_string "minkowski(1)"
  and init = [| 0.1; 0.1; 0.2; 0.2; 0.2; 0.7; 0.5; 0.99; 0.999; 0.05; 0.4; 0.05 |] in
  let len = Array.length init in
  let it = Distance.Iterator.make ~max_distance_component:0.3 distance 1. (fun i -> init.(i)) len in
  let res = Distance.Iterator.get_opt it |> ref in
  while !res <> None do
    let idx_lo, idx_hi, comp = Tools.Option.unbox !res in
    Printf.printf "(%d, %d): %.15g\n" idx_lo idx_hi comp;
    Distance.Iterator.output_summary it;
    Distance.Iterator.incr ~max_distance_component:0.3 it;
    res := Distance.Iterator.get_opt it
  done
