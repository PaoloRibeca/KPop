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

open BiOCamLib
open Better

let ( .@() ) = Float.Array.( .@() )
let ( .@()<- ) = Float.Array.( .@()<- )

(* Extends BiOCamLib matrix class with distance machinery.
   Encapsulation checks are not performed at this level *)
module Base:
  sig
    include module type of Matrix
    (* Compute row normalisations *)
    val get_normalizations: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                            Space.Distance.t -> Float.Array.t -> t -> Float.Array.t
    (* Get embeddings (principal coordinates) from standard coordinates *)
    val get_embeddings: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                        Space.Distance.t -> Float.Array.t -> t -> t
    val get_distance_matrix: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                             Space.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    val get_distance_rowwise: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                              Space.Distance.t -> Float.Array.t -> t -> t -> t
  end
= struct
    include Matrix
    (* Compute normalisations for rows *)
    let get_normalizations ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) distance metric m =
      let n_rows = Array.length m.row_names and n_cols = Array.length m.col_names in
      let res = Float.Array.create n_rows
      and rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0 in
      (* Generate points to be computed by the parallel processs *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_rows then
            let to_do = min rows_per_step (n_rows - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          let res = ref [] in
          (* We iterate backwards so as to avoid to have to reverse the list in the end *)
          for i = hi_row downto lo_row do
            Space.Distance.compute_norm distance metric m.data.(i) |> List.accum res
          done;
          lo_row, !res)
        (fun (lo_row, norms) ->
          List.iteri
            (fun offs_i norm_i ->
              res.@(lo_row + offs_i) <- if norm_i = 0. then 1. else norm_i;
              if verbose && !processed_rows mod elements_per_step = 0 then
                Printf.eprintf "%s\r(%s): Done %d/%d rows%!"
                  String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
              incr processed_rows)
            norms)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d rows.\n%!" String.TermIO.clear __FUNCTION__ !processed_rows n_rows;
      res
    (* Compute embeddings (principal coordinates) from standard coordinates *)
    let get_embeddings ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      let d = Float.Array.length metric in
      if Array.length m.col_names <> d then
        Incompatible_geometries (Array.make d "", m.col_names) |> raise;
      let inv_power =
        match distance with
        | Space.Distance.Euclidean | Cosine -> 0.5
        | Minkowski p -> 1. /. p in
      let normalized_metric = Float.Array.map (fun x -> x ** inv_power) metric
      and rows_per_step = max 1 (elements_per_step / d) and processed_rows = ref 0
      and n = Array.length m.row_names in
      let data = Array.make n (Float.Array.create 0) in
      (* Generate points to be computed by the parallel process *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n then
            let to_do = min rows_per_step (n - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          let res = ref [] in
          (* We iterate backwards so as to avoid to have to reverse the list in the end *)
          for i = hi_row downto lo_row do
            let data_row = m.data.(i) in
            let v = Float.Array.init d (fun col -> data_row.@(col) *. normalized_metric.@(col)) in
            if normalize then begin
              let norm = Space.Distance.compute_norm distance metric v in
              if norm <> 0. then
                Float.Array.iteri (fun i x -> v.@(i) <- (x /. norm)) v
            end;
            List.accum res v
          done;
          lo_row, !res)
        (fun (lo_row, rows) ->
          List.iteri
            (fun offs_i row_i ->
              data.(lo_row + offs_i) <- row_i;
              if verbose && !processed_rows mod rows_per_step = 0 then
                Printf.eprintf "%s\r(%s): Done %d/%d rows%!"
                  String.TermIO.clear __FUNCTION__ !processed_rows n;
              incr processed_rows)
            rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d rows.\n%!" String.TermIO.clear __FUNCTION__ !processed_rows n;
      { m with data = data }
    (* Compute rowwise distance *)
    let get_distance_matrix ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      let d = Array.length m.row_names in
      (* We compute normalisations *)
      let norms =
        if normalize then
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m
        else
          Float.Array.make d 1. in
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let data = Array.init d (fun _ -> Float.Array.create d) in
      (* Generate points to be computed by the parallel processs *)
      let total = (d * (d + 1)) / 2 and i = ref 0 and j = ref 0 and elts_done = ref 0 and end_reached = ref (d = 0) in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file;
          (* We only compute 1/2 of the matrix, and symmetrise it at the end of the computation *)
          let res = ref [] in
          begin try
            let cntr = ref 0 in
            while !cntr < elements_per_step do
              List.accum res (!i, !j);
              incr j;
              if !j > !i then begin
                incr i;
                if !i = d then begin
                  end_reached := true;
                  raise Exit
                end;
                j := 0
              end;
              incr cntr
            done
          with Exit -> ()
          end;
          List.rev !res)
        (List.map
          (* We decorate each matrix element coordinate with the respective distance *)
          (fun (i, j) ->
            i, j, begin
              Space.Distance.compute
                ~adaptor_a:(fun a -> a /. norms.@(i)) ~adaptor_b:(fun b -> b /. norms.@(j))
                distance metric m.data.(i) m.data.(j)
            end))
        (List.iter
          (fun (i, j, dist) ->
            (* Only here do we actually fill out the memory for the result *)
            data.(i).@(j) <- dist;
            (* We symmetrise the matrix *)
            data.(j).@(i) <- dist;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "%s\r(%s): Done %d/%d elements%!"
                String.TermIO.clear __FUNCTION__ !elts_done total;
            incr elts_done))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements.\n%!" String.TermIO.clear __FUNCTION__ !elts_done total;
      { col_names = m.row_names;
        row_names = m.row_names;
        data = data }
    let get_distance_rowwise ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m1 m2 =
      if m1.col_names <> m2.col_names then
        Incompatible_geometries (m1.col_names, m2.col_names) |> raise;
      let r1 = Array.length m1.row_names and r2 = Array.length m2.row_names in
      (* We compute normalisations *)
      let n1, n2 =
        if normalize then
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m1,
          get_normalizations ~threads ~elements_per_step ~verbose distance metric m2
        else
          Float.Array.make r1 1., Float.Array.make r2 1. in
      (*
      to_file ~verbose (transpose ~verbose {
        col_names = m1.row_names;
        row_names = [| "Normalizations" |];
        data = [| n1 |]
      }) "N1.txt";
      to_file ~verbose (transpose ~verbose {
        col_names = m2.row_names;
        row_names = [| "Normalizations" |];
        data = [| n2 |]
      }) "N2.txt";
      *)
      (* We immediately allocate all the needed memory, as we already know how much we will need *)
      let data = Array.init r2 (fun _ -> Float.Array.create r1) in
      (* Generate points to be computed by the parallel processs *)
      let prod = r1 * r2 in
      let i = ref 0 and j = ref 0 and elts_done = ref 0 and end_reached = ref (prod = 0) in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !end_reached then
            raise End_of_file;
          let res = ref [] in
          begin try
            let cntr = ref 0 in
            while !cntr < elements_per_step do
              List.accum res (!i, !j);
              incr j;
              if !j = r2 then begin
                incr i;
                if !i = r1 then begin
                  end_reached := true;
                  raise Exit
                end;
                j := 0
              end;
              incr cntr
            done
          with Exit -> ()
          end;
          List.rev !res)
        (List.map
          (* We decorate each matrix element coordinate with the respective distance *)
          (fun (i, j) ->
            i, j, begin
              Space.Distance.compute
                ~adaptor_a:(fun a -> a /. n1.@(i)) ~adaptor_b:(fun b -> b /. n2.@(j))
                distance metric m1.data.(i) m2.data.(j)
            end))
        (List.iter
          (fun (i, j, dist) ->
            data.(j).@(i) <- dist;
            if verbose && !elts_done mod elements_per_step = 0 then
              Printf.eprintf "%s\r(%s): Done %d/%d elements=%.3g%%%!"
                String.TermIO.clear __FUNCTION__
                !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod);
            incr elts_done))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Done %d/%d elements=%.3g%%.\n%!"
          String.TermIO.clear __FUNCTION__
          !elts_done prod (100. *. float_of_int !elts_done /. float_of_int prod);
      { col_names = m1.row_names;
        row_names = m2.row_names;
        data = data }
  end

(* KPop-specialised matrices.
   We include in order not to have a repeated module prefix *)
include (
  struct
    module Type =
      struct
        type t =
          | Distill
          | Twister
          | Inertia
          | Metrics
          | Twisted
          | Vectors
          | DMatrix
        let to_string = function
          | Distill -> "KPopDistill"
          | Twister -> "KPopTwister"
          | Inertia -> "KPopInertia"
          | Metrics -> "KPopMetrics"
          | Twisted -> "KPopTwisted"
          | Vectors -> "KPopVectors"
          | DMatrix -> "KPopDMatrix"
        exception Unknown_type of string
        let of_string = function
            "KPopDistill" -> Distill
          | "KPopTwister" -> Twister
          | "KPopInertia" -> Inertia
          | "KPopMetrics" -> Metrics
          | "KPopTwisted" -> Twisted
          | "KPopVectors" -> Vectors
          | "KPopDMatrix" -> DMatrix
          | w -> Unknown_type w |> raise
      end
    type t = {
      which: Type.t;
      matrix: Base.t
    }
    let empty which =
      { which; matrix = Base.empty }
    (* The three following functions implement automatic file naming *)
    let make_filename_table which = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ "." ^ Type.to_string which ^ ".txt"
    let make_filename_binary which name =
      match which, name with
      | _, _ when String.length name >= 5 && String.sub name 0 5 = "/dev/" -> name
      | Type.Distill, prefix | Metrics, prefix | Twisted, prefix | Vectors, prefix | DMatrix, prefix ->
        prefix ^ "." ^ Type.to_string which
      | Twister, _ | Inertia, _ -> assert false (* Should always be done through KPopTwister *)
    let make_filename_summary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopSummary.txt"
    (* We redefine the implementation for Matrix in order to set the correct KPop types *)
    let of_file ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) which prefix =
      { which; matrix = make_filename_table which prefix |> Base.of_file ~threads ~bytes_per_step ~verbose }
    let to_file ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) m prefix =
      make_filename_table m.which prefix |> Base.to_file ~precision ~threads ~elements_per_step ~verbose m.matrix
    let transpose_single_threaded ?(verbose = false) m =
      { m with matrix = Base.transpose_single_threaded ~verbose m.matrix }
    let transpose ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m =
      { m with matrix = Base.transpose ~threads ~elements_per_step ~verbose m.matrix }
    exception Incompatible_matrices of Type.t * Type.t
    let merge_rowwise ?(verbose = false) m1 m2 =
      if m1.which <> m2.which then
        Incompatible_matrices (m1.which, m2.which) |> raise;
      { which = m1.which; matrix = Base.merge_rowwise ~verbose m1.matrix m2.matrix }
    let multiply_matrix_vector_single_threaded ?(verbose = false) m =
      Base.multiply_matrix_vector_single_threaded ~verbose m.matrix
    let multiply_matrix_sparse_vector_single_threaded ?(verbose = false) m =
      Base.multiply_matrix_sparse_vector_single_threaded ~verbose m.matrix
    let multiply_matrix_vector ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) m v =
      Base.multiply_matrix_vector ~threads ~elements_per_step ~verbose m.matrix v
    let multiply_matrix_matrix ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) which m1 m2 =
      { which; matrix = Base.multiply_matrix_matrix ~threads ~elements_per_step ~verbose m1.matrix m2.matrix }
    exception Unexpected_type of Type.t * Type.t
    let get_embeddings ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      if m.which <> Twisted then
        Unexpected_type (m.which, Twisted) |> raise;
      { which = Vectors;
        matrix = Base.get_embeddings ~normalize ~threads ~elements_per_step ~verbose distance metric m.matrix }
    module SplitsAlgorithm =
      struct
        type t =
          | Gaps
          | Centroids
        let to_string = function
          | Gaps -> "gaps"
          | Centroids -> "centroids"
        exception Unknown_algorithm of string
        let of_string = function
          | "gaps" -> Gaps
          | "centroids" -> Centroids
          | w -> Unknown_algorithm w |> raise
        (* Implementation module *)
        module Bipartition =
          struct
            exception Invalid_acceptance_probability_at_zero of float
            exception Invalid_difference_magnification_factor of float
            exception Insufficient_elements of int
            exception Invalid_element of int * int
            let make
                ?(acceptance_probability_at_zero = 0.2) ?(difference_magnification_factor = 10.) ?(verbose = false)
                m init_set =
              if acceptance_probability_at_zero <= 0. || acceptance_probability_at_zero > 1. then
                Invalid_acceptance_probability_at_zero acceptance_probability_at_zero |> raise;
              if difference_magnification_factor <= 0. then
                Invalid_difference_magnification_factor difference_magnification_factor |> raise;
              if m.which <> Vectors then
                Unexpected_type (m.which, Vectors) |> raise;
              let inverse_acceptance = (1. -. acceptance_probability_at_zero) /. acceptance_probability_at_zero
              and negative_scale = -. difference_magnification_factor
              and elements = IntSet.elements_array init_set in
              let num_elements = Array.length elements in
              if num_elements < 2 then
                Insufficient_elements num_elements |> raise;
              let n = Array.length m.matrix.row_names
              and d = Array.length m.matrix.col_names
              and one = ref IntSet.empty and cardinal_one = ref 0
              and two = ref IntSet.empty and cardinal_two = ref 0 in
              let centroid_one = Float.Array.make d 0. |> ref
              and centroid_two = Float.Array.make d 0. |> ref in
              IntSet.iter
                (fun i ->
                  if i >= n then
                    Invalid_element (i, n) |> raise;
                  let v = m.matrix.data.(i) in
                  (* We randomly assign the element to either set *)
                  if Random.bool () then begin
                    two := IntSet.add i !two;
                    incr cardinal_two;
                    Float.Array.iteri
                      (fun j x -> !centroid_two.@(j) <- !centroid_two.@(j) +. x)
                      v
                  end else begin
                    one := IntSet.add i !one;
                    incr cardinal_one;
                    Float.Array.iteri
                      (fun j x -> !centroid_one.@(j) <- !centroid_one.@(j) +. x)
                      v
                  end)
                init_set;
              (* Temporary space *)
              let old_centroid_one = Float.Array.make d 0. |> ref
              and old_centroid_two = Float.Array.make d 0. |> ref
              and compute_objective () =
                let normalize sum card =
                  if card > 1. then
                    sum /. card
                  else
                    sum in
                let cardinal_one = float_of_int !cardinal_one
                and cardinal_two = float_of_int !cardinal_two and res = ref 0. in
                if cardinal_one > 0. && cardinal_two > 0. then
                  Float.Array.iter2
                    (fun sum_one sum_two ->
                      let min, max = min_max (normalize sum_one cardinal_one) (normalize sum_two cardinal_two) in
                      res := !res +. (max -. min))
                    !centroid_one !centroid_two;
                !res /. sqrt (1. +. Float.abs (cardinal_one -. cardinal_two)) in
              let objective = compute_objective () |> ref in

              if verbose then begin
                Printf.eprintf "(%s): Begin (objective=%.3g, one=[" __FUNCTION__ !objective;
                IntSet.iter (Printf.eprintf " %d") !one;
                Printf.eprintf " ], two=[";
                IntSet.iter (Printf.eprintf " %d") !two;
                Printf.eprintf " ])\n%!"
              end;

              (* All-time bests *)
              let max_objective = ref !objective and max_one = ref !one and max_two = ref !two
              and terminator = max num_elements 40 and rejected = ref 0 and steps = ref 0 in
              (* We stop if no improvement happens for as many moves as the number of elements *)
              while !rejected < terminator do

                if !steps mod 1000 = 0 then
                  Printf.eprintf " Step #%d: objective=%.3g, max_objective=%.3g\n%!" !steps !objective !max_objective;
                incr steps;

                (* We save the old state *)
                let old_one = !one and old_cardinal_one = !cardinal_one
                and old_two = !two and old_cardinal_two = !cardinal_two
                and old_objective = !objective in
                Float.Array.blit !centroid_one 0 !old_centroid_one 0 d;
                Float.Array.blit !centroid_two 0 !old_centroid_two 0 d;
                let selected = elements.(Random.int num_elements) in
                let v = m.matrix.data.(selected) in
                if IntSet.mem selected !one then begin
                  (* Move element from partition one to partition two *)
                  one := IntSet.remove selected !one;
                  decr cardinal_one;
                  two := IntSet.add selected !two;
                  incr cardinal_two;
                  Float.Array.iter2i
                    (fun i sum_one sum_two ->
                      let coord = v.@(i) in
                      !centroid_one.@(i) <- sum_one -. coord;
                      !centroid_two.@(i) <- sum_two +. coord)
                    !old_centroid_one !old_centroid_two
                end else begin
                  (* Move element from partition two to partition one *)
                  two := IntSet.remove selected !two;
                  decr cardinal_two;
                  one := IntSet.add selected !one;
                  incr cardinal_one;
                  Float.Array.iter2i
                    (fun i sum_one sum_two ->
                      let coord = v.@(i) in
                      !centroid_two.@(i) <- sum_two -. coord;
                      !centroid_one.@(i) <- sum_one +. coord)
                    !old_centroid_one !old_centroid_two
                end;
                objective := compute_objective ();
                (* Should we accept the move? *)
                let delta = !objective -. old_objective in
                let score = 1. /. (1. +. inverse_acceptance *. exp (negative_scale *. delta)) in
                if Random.float 1. <= score then begin
                  (* Accept *)
                  rejected := 0;
                  if !objective > !max_objective then begin
                    (* Update all-time minimum *)
                    max_objective := !objective;
                    max_one := !one;
                    max_two := !two
                  end
                end else begin
                  (* Reject *)
                  incr rejected;
                  (* We have to restore the previous state *)
                  one := old_one;
                  cardinal_one := old_cardinal_one;
                  two := old_two;
                  cardinal_two := old_cardinal_two;
                  let tmp = !centroid_one in
                  centroid_one := !old_centroid_one;
                  old_centroid_one := tmp;
                  let tmp = !centroid_two in
                  centroid_two := !old_centroid_two;
                  old_centroid_two := tmp;
                  objective := old_objective
                end
              done;

              if verbose then begin
                Printf.eprintf "(%s): End (objective=%.3g, max=%.3g, one=[" __FUNCTION__ !objective !max_objective;
                IntSet.iter (fun i -> m.matrix.row_names.(i) |> Printf.eprintf " '%s'") !max_one;
                Printf.eprintf " ], two=[";
                IntSet.iter (fun i -> m.matrix.row_names.(i) |> Printf.eprintf " '%s'") !max_two;
                Printf.eprintf " ], steps=%d)\n%!" !steps
              end;

              !max_one, !max_two, !max_objective, !steps
          end
      end
    let get_splits ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false) algorithm_type max_splits m =
      if m.which <> Vectors then
        Unexpected_type (m.which, Vectors) |> raise;
      match algorithm_type with
      | SplitsAlgorithm.Gaps ->
        (* Embeddings are stored rowwise.
          We begin by sorting coordinates along each dimension (i.e., by sorting columns) *)
        let n = Array.length m.matrix.row_names in
        let cols_per_step = max 1 (elements_per_step / n) and processed_cols = ref 0
        and d = Array.length m.matrix.col_names in
        let row_permutations = Array.make d [||] and gaps = Tools.StackArray.create () in
        (* Generate points to be computed by the parallel process *)
        Processes.Parallel.process_stream_chunkwise
          (fun () ->
            if !processed_cols < d then
              let to_do = min cols_per_step (d - !processed_cols) in
              let new_processed_cols = !processed_cols + to_do in
              let res = !processed_cols, new_processed_cols - 1 in
              processed_cols := new_processed_cols;
              res
            else
              raise End_of_file)
          (fun (lo_col, hi_col) ->
            let res = ref [] in
            (* We iterate backwards so as to avoid to have to reverse the list in the end *)
            for i = hi_col downto lo_col do
              (* We annotate the transposed value with its row index *)
              let coords__idxs = Array.init n (fun row -> m.matrix.data.(row).@(i), row) in
              (* We sort the vector *)
              Array.sort (fun (coord_1, _) (coord_2, _) -> compare coord_1 coord_2) coords__idxs;
              (* We compute gaps, i.e. differences between consecutive coordinates.
                Gaps are annotated with their indices *)
              let gaps__idxs = Array.init (n - 1) (fun j -> fst coords__idxs.(j + 1) -. fst coords__idxs.(j), i, j) in
              (* We sort the vector *)
              Array.sort (fun (gap_1, _, _) (gap_2, _, _) -> compare gap_1 gap_2) gaps__idxs;
              (* We return the permutation of row indices and the gap vector *)
              List.accum res (Array.init n (fun row -> snd coords__idxs.(row)), gaps__idxs)
            done;
            lo_col, !res)
          (fun (lo_col, cols) ->
            List.iteri
              (fun offs_i (perm_i, gaps_i) ->
                row_permutations.(lo_col + offs_i) <- perm_i;
                Tools.StackArray.push_array gaps gaps_i;
                if verbose && !processed_cols mod cols_per_step = 0 then
                  Printf.eprintf "%s\r(%s): Done %d/%d cols%!"
                    String.TermIO.clear __FUNCTION__ !processed_cols n;
                incr processed_cols)
              cols)
          threads;
        if verbose then
          Printf.eprintf "%s\r(%s): Done %d/%d cols.\n%!"
            String.TermIO.clear __FUNCTION__ !processed_cols n;
        (* We sort the gaps *)
        let gaps = Tools.StackArray.contents gaps in
        Array.sort
          (fun (gap_1, dim_1, idx_1) (gap_2, dim_2, idx_2) ->
            (* We sort splits by decreasing gap size first, then by increasing dimension (and row index) *)
            let rgap = compare gap_2 gap_1 in
            if rgap <> 0 then
              rgap
            else begin
              let dim = compare dim_1 dim_2 in
              if dim <> 0 then
                dim
              else
                compare idx_1 idx_2
            end)
          gaps;
        (* We generate splits from the selected number of gaps *)
        let res = Trees.Splits.create m.matrix.row_names in
        for i = 0 to (Array.length gaps |> min max_splits) - 1 do
          let gap, dim, idx = gaps.(i) in
          let split = Array.sub row_permutations.(dim) 0 (idx + 1) |> Trees.Splits.Split.of_array in
          Trees.Splits.add_split res split gap
        done;
        res
      | Centroids ->
        let res = Trees.Splits.create m.matrix.row_names in
        let rec refine_by_bipartition set =
          if IntSet.cardinal set > 1 then begin
            (* Bipartition.evolve () should work fine provided that there are at least 2 elements *)
            let one, two, objective, _ = SplitsAlgorithm.Bipartition.make ~verbose m set in
            Trees.Splits.add_split res (IntSet.elements_array one |> Trees.Splits.Split.of_array) objective;
            refine_by_bipartition one;
            refine_by_bipartition two
          end else
            Trees.Splits.add_split res (IntSet.elements_array set |> Trees.Splits.Split.of_array) 0. in
        Seq.init (Array.length m.matrix.row_names) (fun i -> i) |> IntSet.of_seq |> refine_by_bipartition;
        res
    let get_distance_matrix ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m =
      if m.which <> Twisted then
        Unexpected_type (m.which, Twisted) |> raise;
      { which = DMatrix;
        matrix = Base.get_distance_matrix ~normalize ~threads ~elements_per_step ~verbose distance metric m.matrix }
    (* Compute distances between the rows of two matrices - more general version of the previous one *)
    let get_distance_rowwise ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m1 m2 =
      if m1.which <> Twisted then
        Unexpected_type (m1.which, Twisted) |> raise;
      if m2.which <> Twisted then
        Unexpected_type (m2.which, Twisted) |> raise;
      { which = DMatrix;
        matrix =
          Base.get_distance_rowwise ~normalize ~threads ~elements_per_step ~verbose
            distance metric m1.matrix m2.matrix }
    module FloatIntMultimap = Tools.Multimap (ComparableFloat) (ComparableInt)
    let summarize_distance_matrix_row req_len row_name row col_names buf =
      let n_cols = Float.Array.length row in
      let f_n_cols = float_of_int n_cols and distr = ref FloatIntMultimap.empty in
      (* We find the median and filter the result *)
      Float.Array.iteri
        (fun col_idx dist ->
          distr := FloatIntMultimap.add dist col_idx !distr)
        row;
      let eff_len = ref 0 and median_pos = n_cols / 2 |> ref and median = ref 0. and acc = ref 0. in
      FloatIntMultimap.iter_set
        (fun dist set ->
          let set_len = FloatIntMultimap.ValSet.cardinal set in
          acc := !acc +. (float_of_int set_len *. dist);
          if !median_pos >= 0 && !median_pos - set_len < 0 then
            median := dist;
          median_pos := !median_pos - set_len;
          if !eff_len < req_len then
            eff_len := !eff_len + set_len)
        !distr;
      let eff_len = !eff_len and median = !median and mean =
        if n_cols > 0 then
          !acc /. f_n_cols
        else
          0. in
      (* We compute standard deviation e MAD *)
      acc := 0.;
      let ddistr = ref FloatMap.empty in
      Float.Array.iteri
        (fun _ dist ->
          let d = dist -. mean in
          acc := !acc +. (d *. d);
          let d = (dist -. median) |> abs_float in
          ddistr :=
            match FloatMap.find_opt d !ddistr with
            | None ->
              FloatMap.add d 1 !ddistr
            | Some n ->
              FloatMap.add d (n + 1) !ddistr)
        row;
      median_pos := n_cols / 2;
      let mad = ref 0. in
      FloatMap.iter
        (fun d occs ->
          if !median_pos >= 0 && !median_pos - occs < 0 then
            mad := d;
          median_pos := !median_pos - occs)
        !ddistr;
      let mad = !mad and stddev =
        if n_cols > 1 then
          !acc /. (f_n_cols -. 1.) |> sqrt
        else
          0. in
      Printf.bprintf buf "%s\t%.15g\t%.15g\t%.15g\t%.15g" row_name mean stddev median mad;
      FloatIntMultimap.iteri
        (fun i dist col_idx ->
          if i < eff_len then
            Printf.bprintf buf "\t%s\t%.15g\t%.15g" col_names.(col_idx) dist ((dist -. mean) /. stddev))
        !distr;
      Printf.bprintf buf "\n"
    let summarize_rowwise
        ?(normalize = true) ?(keep_at_most = Some 2) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric m1 m2 prefix =
      if m1.which <> Twisted then
        Unexpected_type (m1.which, Twisted) |> raise;
      if m2.which <> Twisted then
        Unexpected_type (m2.which, Twisted) |> raise;
      if m1.matrix.col_names <> m2.matrix.col_names then
        Base.Incompatible_geometries (m1.matrix.col_names, m2.matrix.col_names) |> raise;
      let r1 = Array.length m1.matrix.row_names and r2 = Array.length m2.matrix.row_names in
      (* We compute normalisations *)
      let n1, n2 =
        if normalize then
          Base.get_normalizations ~threads ~elements_per_step ~verbose distance metric m1.matrix,
          Base.get_normalizations ~threads ~elements_per_step ~verbose distance metric m2.matrix
        else
          Float.Array.make r1 1., Float.Array.make r2 1. in
      (*
      Base.to_file ~verbose (Base.transpose ~verbose {
        col_names = m1.matrix.row_names;
        row_names = [| "Normalizations" |];
        data = [| n1 |]
      }) "N1.txt";
      Base.to_file ~verbose (Base.transpose ~verbose {
        col_names = m2.matrix.row_names;
        row_names = [| "Normalizations" |];
        data = [| n2 |]
      }) "N2.txt";
      *)
      let fname = make_filename_summary prefix in
      let output = open_out fname
      and n_cols = Array.length m1.matrix.col_names in
      let req_len =
        match keep_at_most with
        | None -> r1
        | Some at_most -> at_most in
      let rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0
      and buf = Buffer.create 1048576 in
      (* Parallel section *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < r2 then
            let to_do = min rows_per_step (r2 - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          Buffer.clear buf;
          for j = lo_row to hi_row do
            (* For each row number of m2, we compute the respective distances from the rows of m1... *)
            let distances =
              Float.Array.init r1
                (fun i ->
                  Space.Distance.compute
                    ~adaptor_a:(fun a -> a /. n1.@(i)) ~adaptor_b:(fun b -> b /. n2.@(j))
                    distance metric m1.matrix.data.(i) m2.matrix.data.(j)) in
            (* ...and summarise them *)
            summarize_distance_matrix_row
              req_len m2.matrix.row_names.(j) distances m1.matrix.row_names buf
          done;
          hi_row - lo_row + 1, Buffer.contents buf)
        (fun (n_processed, block) ->
          Printf.fprintf output "%s" block;
          let new_processed_rows = !processed_rows + n_processed in
          if verbose && new_processed_rows / rows_per_step > !processed_rows / rows_per_step then
            Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows%!"
              String.TermIO.clear __FUNCTION__ fname new_processed_rows r2;
          processed_rows := new_processed_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows.\n%!"
          String.TermIO.clear __FUNCTION__ fname r2 r2;
      close_out output
    let summarize_distance ?(keep_at_most = Some 2) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        m prefix =
      if m.which <> DMatrix then
        Unexpected_type (m.which, DMatrix) |> raise;
      let fname = make_filename_summary prefix in
      let output = open_out fname
      and n_cols = Array.length m.matrix.col_names in
      let req_len =
        match keep_at_most with
        | None -> n_cols
        | Some at_most -> at_most
      and n_rows = Array.length m.matrix.row_names in
      (* Parallel section *)
      let rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0
      and buf = Buffer.create 1048576 in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_rows then
            let to_do = min rows_per_step (n_rows - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          Buffer.clear buf;
          for i = lo_row to hi_row do
            summarize_distance_matrix_row
              req_len m.matrix.row_names.(i) m.matrix.data.(i) m.matrix.col_names buf
          done;
          hi_row - lo_row + 1, Buffer.contents buf)
        (fun (n_processed, block) ->
          Printf.fprintf output "%s" block;
          let new_processed_rows = !processed_rows + n_processed in
          if verbose && new_processed_rows / rows_per_step > !processed_rows / rows_per_step then
            Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows%!"
              String.TermIO.clear __FUNCTION__ fname new_processed_rows n_rows;
          processed_rows := new_processed_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows.\n%!"
          String.TermIO.clear __FUNCTION__ fname n_rows n_rows;
      close_out output
    (* *)
    let archive_version = "2022-04-03"
    (* *)
    let to_channel output m =
      Type.to_string m.which |> output_value output;
      archive_version |> output_value output;
      output_value output m.matrix
    let to_binary ?(verbose = false) m prefix =
      let fname = make_filename_binary m.which prefix in
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      to_channel output m;
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    exception Incompatible_archive_version of string * string * string
    let of_channel input =
      let which = (input_value input: string) in
      let version = (input_value input: string) in
      if version <> archive_version then
        Incompatible_archive_version (which, version, archive_version) |> raise;
      { which = Type.of_string which; matrix = (input_value input: Base.t) }
    let of_binary ?(verbose = false) which prefix =
      let fname = make_filename_binary which prefix in
      let input = open_in fname in
      if verbose then
        Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let res = of_channel input in
      close_in input;
      if which <> res.which then
        Unexpected_type (which, res.which) |> raise;
      if verbose then
        Printf.eprintf " done.\n%!";
      res
  end: sig
    module Type:
      sig
        type t =
          | Distill
          | Twister
          | Inertia
          | Metrics
          | Twisted
          | Vectors
          | DMatrix
        val to_string: t -> string
        exception Unknown_type of string
        val of_string: string -> t
      end
    type t = {
      which: Type.t;
      matrix: Base.t
    }
    val empty: Type.t -> t
    (* All file name arguments are in fact _prefixes_ *)
    val of_file: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> Type.t -> string -> t
    (* This one discards type information - use at your own risk *)
    val to_file: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val transpose_single_threaded: ?verbose:bool -> t -> t
    val transpose: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> t
    exception Incompatible_matrices of Type.t * Type.t
    (* Merge two matrices - the type of the two inputs must be the same *)
    val merge_rowwise: ?verbose:bool -> t -> t -> t
    (* TODO: No type checks are performed (yet) when multiplying matrices *)
    val multiply_matrix_vector_single_threaded: ?verbose:bool -> t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_sparse_vector_single_threaded: ?verbose:bool -> t -> Base.sparse_vector_t -> Float.Array.t
    val multiply_matrix_vector: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                                t -> Float.Array.t -> Float.Array.t
    val multiply_matrix_matrix: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> Type.t -> t -> t -> t
    exception Unexpected_type of Type.t * Type.t
    (* Compute embeddings (principal coordinates) from standard coordinates - input type must be Twisted *)
    val get_embeddings: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                        Space.Distance.t -> Float.Array.t -> t -> t
    module SplitsAlgorithm:
      sig
        type t =
          | Gaps
          | Centroids
        val to_string: t -> string
        exception Unknown_algorithm of string
        val of_string: string -> t
      end
    (* Compute splits from embeddings - input type must be Vectors *)
    val get_splits: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                    SplitsAlgorithm.t -> int -> t -> Trees.Splits.t
    (* Compute distances between the rows of a matrix - input type must be Twisted *)
    val get_distance_matrix: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                             Space.Distance.t -> Float.Array.t -> t -> t
    (* Compute distances between the rows of two matrices - input types must be Twisted.
       More general version of the previous one *)
    val get_distance_rowwise: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                              Space.Distance.t -> Float.Array.t -> t -> t -> t
    (* Compute distances between the rows of two matrices and summarise them - input types must be Twisted *)
    val summarize_rowwise: ?normalize:bool -> ?keep_at_most:int option -> ?threads:int -> ?elements_per_step:int ->
                           ?verbose:bool ->
                           Space.Distance.t -> Float.Array.t -> t -> t -> string -> unit
    (* Summarise distances - input type must be DMatrix *)
    val summarize_distance: ?keep_at_most:int option -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                            t -> string -> unit
    (* Binary marshalling of the matrix *)
    val to_channel: out_channel -> t -> unit
    exception Incompatible_archive_version of string * string * string
    val of_channel: in_channel -> t
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> Type.t -> string -> t
  end
)

