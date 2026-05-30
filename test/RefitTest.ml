(* Validation for lib/Refit.ml: load a KPopDMatrix and a Newick tree,
   refit branch lengths by sparse OLS against the matrix, and report
   the fit quality (Pearson of refit-tree patristic vs target on the
   fitted pairs).  Compared externally against
   test/Trees/ols_branch_refit_sparse.py on the same inputs. *)

open BiOCamLib.Better
module Trees = BiOCamLib.Trees
open KPop

let load_dmatrix path =
  let ic = open_in path in
  let header = input_line ic in
  let names = Array.of_list (List.tl (String.split_on_char '\t' header)) in
  let n = Array.length names in
  let d = Array.make_matrix n n 0. in
  for i = 0 to n - 1 do
    let line = input_line ic in
    let fields = Array.of_list (String.split_on_char '\t' line) in
    for j = 0 to n - 1 do
      d.(i).(j) <- float_of_string fields.(j + 1)
    done
  done;
  close_in ic;
  names, d

let () =
  let dmx = Sys.argv.(1) and tree_f = Sys.argv.(2) and out_f = Sys.argv.(3) in
  let strategy = if Array.length Sys.argv > 4 then Sys.argv.(4) else "knn" in
  let k = if Array.length Sys.argv > 5 then int_of_string Sys.argv.(5) else 10 in
  let names, d = load_dmatrix dmx in
  let idx = Hashtbl.create (Array.length names) in
  Array.iteri (fun i nm -> Hashtbl.replace idx nm i) names;
  let dist a b = d.(Hashtbl.find idx a).(Hashtbl.find idx b) in
  let tree =
    Trees.Newick.of_file
      ~negative_branches:Trees.Newick.NegativeBranchesPolicy.Zero tree_f in
  let refit =
    Refit.refit ~verbose:true
      ~strategy:(Refit.Strategy.of_string strategy) ~k ~dist tree in
  Trees.Newick.to_file ~rich_format:false refit out_f;
  Printf.printf "Refit written to %s (strategy=%s, K=%d, n=%d)\n%!"
    out_f strategy k (Array.length names)
