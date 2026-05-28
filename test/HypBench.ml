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

(* Sparse-NJ benchmark driver.  Loads a flat-tabular embedding file
   (one line per leaf: <name>\t<coord_1>\t...\t<coord_d>), runs
   SparseNJ.compute with the given parameters, writes the resulting
   Newick to stdout.  Used by the larger-n hyperbolic validation
   experiment in test/Trees/hyp_scaling.py: it bypasses the
   .KPopTwister / .KPopTwisted file format entirely, since the
   experiment uses synthetic embeddings that don't go through CA.

   Usage:
     HypBench.exe <emb_file> <mode> <K> <K_QUERY_factor> <hyp_scale>
                  [<rowsum>] [<symmetry>]

   where <mode> is one of "dense" | "subquadratic" | "hyperbolic". *)

open BiOCamLib.Better
module Trees = BiOCamLib.Trees
open KPop

let load_embeddings path =
  let names = ref [] and rows = ref [] in
  let ic = open_in path in
  (try
    while true do
      let line = input_line ic in
      let parts = String.split_on_char '\t' line in
      match parts with
      | [] -> ()
      | name :: coords ->
        names := name :: !names;
        let arr = Float.Array.of_list (List.map float_of_string coords) in
        rows := arr :: !rows
    done
  with End_of_file -> ());
  close_in ic;
  Array.of_rlist !names, Array.of_rlist !rows

let () =
  if Array.length Sys.argv < 6 then begin
    Printf.eprintf "Usage: %s <emb_file> <mode> <K> <K_QUERY_factor> <hyp_scale> [<rowsum>] [<symmetry>]\n%!"
      Sys.argv.(0);
    exit 2
  end;
  let emb_path = Sys.argv.(1) in
  let mode = SparseNJ.Mode.of_string Sys.argv.(2) in
  let k_nn = int_of_string Sys.argv.(3) in
  let k_query_factor = int_of_string Sys.argv.(4) in
  let hyp_scale = float_of_string Sys.argv.(5) in
  let row_sum =
    if Array.length Sys.argv >= 7 then SparseNJ.RowSum.of_string Sys.argv.(6)
    else SparseNJ.RowSum.Knn in
  let symmetry =
    if Array.length Sys.argv >= 8 then SparseNJ.Symmetry.of_string Sys.argv.(7)
    else SparseNJ.Symmetry.One in
  let names, data = load_embeddings emb_path in
  Printf.eprintf "Loaded %d leaves in %d-D embedding.\n%!"
    (Array.length names)
    (if Array.length data = 0 then 0 else Float.Array.length data.(0));
  let tree =
    SparseNJ.compute
      ~verbose:false ~mode
      ~index_type:(Interfaiss.Type.of_string "flat")
      ~k_nn ~k_query_factor ~hyp_scale ~row_sum ~symmetry
      names data in
  print_string (Trees.Newick.to_string tree);
  print_newline ()
