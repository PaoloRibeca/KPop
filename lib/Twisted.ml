(*
    Twisted.ml -- (c) 2022-2026 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    Twisted.ml implements the Twisted database (CA-embedded sample
    coordinates) together with the operations that consume it: FAISS-
    backed nearest-neighbour search, distance summarisation, embedding
    output, the three phylogenetic-splits algorithms (gaps, centroids
    with multi-seed bootstrap, and HDBSCAN with dense / sparse / auto
    MST modes), and the binary `.KPopTwisted` file format.

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

(* We cannot open BiOCamLib here due to the ambiguity between BiOCamLib.Matrix and KPop.Matrix *)
module Numbers = BiOCamLib.Numbers
module Processes = BiOCamLib.Processes
module Tools = BiOCamLib.Tools
module Trees = BiOCamLib.Trees
open BiOCamLib.Better

(* Policy mandating how many additional neighbours should be visited to compute statistics
    in addition to the requested ones *)
module NeighborsPolicy:
  sig
    (* We make the type private to implement constraints *)
    type t = private
      | Proportional of float (* Factor must be >= 1. *)
      | Additional of int (* Number must be >= 0 *)
    val to_string: t -> string
    val of_string: string -> t
    val get_to_be_visited: t -> int -> int option -> int
  end
= struct
    type t =
      | Proportional of float (* Factor must be >= 1. *)
      | Additional of int (* Number must be >= 0 *)
    let to_string = function
      | Proportional f -> "times(" ^ string_of_float f ^ ")"
      | Additional n -> "plus(" ^ string_of_int n ^ ")"
    let of_string_re = Str.regexp "[()]"
    let of_string s =
      let fail kind =
        Exception.raise __FUNCTION__ Initialize (Printf.sprintf "%s policy '%s'" kind s) in
      match Str.full_split of_string_re s with
      | [ Text "times"; Delim "("; Text mult; Delim ")" ] ->
        let mult =
          try
            float_of_string mult
          with _ ->
            fail "Invalid" in
        if mult < 1. then
          fail "Invalid";
        Proportional mult
      | [ Text "plus"; Delim "("; Text add; Delim ")" ] ->
        let add =
          try
            int_of_string add
          with _ ->
            fail "Unknown" in
        if add < 0 then
          fail "Invalid";
        Additional add
      | _ ->
        fail "Unknown"
    let get_to_be_visited p n k =
      match k, p with
      | None, _ ->
        n
      | Some k, Proportional f ->
        f *. float_of_int k +. 0.5 |> int_of_float |> min n
      | Some k, Additional i ->
        min (k + i) n
  end

module SplitsAlgorithm:
  sig
    type t =
      | Gaps
      | Centroids
      | Hdbscan
    val to_string: t -> string
    val of_string: string -> t
  end
= struct
    type t =
      | Gaps
      | Centroids
      | Hdbscan
    let of_string = function
      | "gaps" -> Gaps
      | "centroids" -> Centroids
      | "hdbscan" -> Hdbscan
      | s ->
        Exception.raise_unrecognized_initializer __FUNCTION__ "algorithm" s
    let to_string = function
      | Gaps -> "gaps"
      | Centroids -> "centroids"
      | Hdbscan -> "hdbscan"
  end

module BalancePenalty:
  sig
    (* The cardinality-imbalance penalty used as the denominator of the
       centroids objective.  Given the two side cardinalities c1, c2,
       the denominator is computed as (1 + |c1 - c2|^alpha)^beta.
       alpha = 1 and beta = 0.5 reproduce the original sqrt(1 + |c1 - c2|)
       soft penalty exactly. *)
    type t = { alpha: float; beta: float }
    val default: t
    val to_string: t -> string
    val of_string: string -> t
  end
= struct
    type t = { alpha: float; beta: float }
    let default = { alpha = 1.; beta = 0.5 }
    let of_string_re = Str.regexp "[(,)]"
    let of_string s =
      let raise kind =
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "%s balance penalty '%s'" kind s) in
      match Str.full_split of_string_re s with
      | [ Text "penalty"; Delim "(";
          Text alpha; Delim ","; Text beta;
          Delim ")" ] ->
        let alpha, beta =
          try float_of_string alpha, float_of_string beta
          with _ -> raise "Invalid" in
        if alpha < 0. || beta < 0. then raise "Invalid";
        { alpha; beta }
      | _ ->
        raise "Unknown"
    let to_string { alpha; beta } =
      Printf.sprintf "penalty(%.15g,%.15g)" alpha beta
  end

(* Twisted vectors (in standard coordinates) are a combination of KPop matrices.
   We include in order not to have a repeated module prefix *)
include (
  struct
    let ( .@() ) = Float.Array.( .@() )
    let ( .@()<- ) = Float.Array.( .@()<- )
    type t = {
      inertia: Matrix.t;
      twisted: Matrix.t
    }
    let empty = {
      inertia = Matrix.empty Inertia;
      twisted = Matrix.empty Twisted
    }
    let raise_incompatible_inertias __FUNCTION__ =
      Exception.raise_incompatible_arrays __FUNCTION__ "matrices" "inertias" Float.Array.iter string_of_float
    let get_inertias t = t.inertia.matrix.data.(0)
    let merge_rowwise t1 t2 =
      (* Perform some compatibility checks *)
      let inertias1 = get_inertias t1 and inertias2 = get_inertias t2 in
      if inertias1 <> inertias2 then
        raise_incompatible_inertias __FUNCTION__ inertias1 inertias2;
      { inertia = t1.inertia; twisted = Matrix.merge_rowwise t1.twisted t2.twisted }
    (* *)
    let to_embeddings ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric t =
      let metrics = get_inertias t |> Space.Distance.Metric.compute metric in
      if verbose then begin
        Printf.eprintf "(%s): Metrics=(" __FUNCTION__;
        Float.Array.iter (Printf.eprintf " %.6g") metrics;
        Printf.eprintf " )\n%!"
      end;
      { Matrix.which = Vectors;
        matrix =
          Matrix.Base.get_embeddings ~normalize ~threads ~elements_per_step ~verbose
            distance metrics t.twisted.matrix }
    (* *)
    module FloatIntMultimap = Tools.Multimap (ComparableFloat) (ComparableInt)
    let summarize_distance_matrix_row ?(col_names_mapping = fun idx -> idx) req_len row_name row col_names buf =
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
      (* We compute standard deviation and MAD *)
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
            Printf.bprintf buf "\t%s\t%.15g\t%.15g"
              col_names.(col_names_mapping(col_idx)) dist ((dist -. mean) /. stddev))
        !distr;
      Printf.bprintf buf "\n"
    (* *)
    let make_filename_summary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopSummary.txt"
    (* Summarise distances between the rows of two matrices, and output the distance matrix if needed *)
    let summarize_distances_rowwise
        ?(normalize = true) ?(keep_at_most = Some 2) ?(output_distance_matrix = false) ?(precision = 15)
        ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        distance metric t1 t2 prefix =
      (* Perform some compatibility checks *)
      let inertias1 = get_inertias t1 and inertias2 = get_inertias t2 in
      if inertias1 <> inertias2 then
        raise_incompatible_inertias __FUNCTION__ inertias1 inertias2;
      let col_names_1 = t1.twisted.matrix.col_names and col_names_2 = t2.twisted.matrix.col_names in
      if col_names_1 <> col_names_2 then
        Matrix.Exception.raise_incompatible_geometries __FUNCTION__ col_names_1 col_names_2;
      let metrics = Space.Distance.Metric.compute metric inertias1 in
      let row_names_1 = t1.twisted.matrix.row_names and row_names_2 = t2.twisted.matrix.row_names in
      let r1 = Array.length row_names_1 and r2 = Array.length row_names_2 in
      (* We compute normalisations *)
      let n1, n2 =
        if normalize then
          Matrix.Base.get_normalizations ~threads ~elements_per_step ~verbose distance metrics t1.twisted.matrix,
          Matrix.Base.get_normalizations ~threads ~elements_per_step ~verbose distance metrics t2.twisted.matrix
        else
          Float.Array.make r1 1., Float.Array.make r2 1. in
      (*
      Base.to_file ~verbose (Base.transpose ~verbose {
        col_names = t1.matrix.row_names;
        row_names = [| "Normalizations" |];
        data = [| n1 |]
      }) "N1.txt";
      Base.to_file ~verbose (Base.transpose ~verbose {
        col_names = t2.matrix.row_names;
        row_names = [| "Normalizations" |];
        data = [| n2 |]
      }) "N2.txt";
      *)
      let path = make_filename_summary prefix in
      let output = open_out path
      and n_cols = Array.length t1.twisted.matrix.col_names in
      let req_len =
        match keep_at_most with
        | None -> r1
        | Some at_most -> at_most in
      let rows_per_step = max 1 (elements_per_step / n_cols) and processed_rows = ref 0
      and buf = Buffer.create 1048576 and empty = Float.Array.create 0 in
      let data =
        if output_distance_matrix then
          Fun.const empty |> Array.init r2
        else
          [||] in
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
          let rows = ref [] in
          for j = lo_row to hi_row do
            (* For each row index of t2, we compute the respective distances from the rows of t1... *)
            let distances =
              Float.Array.init r1
                (fun i ->
                  Space.Distance.compute
                    ~adaptor_a:(fun a -> a /. n1.@(i)) ~adaptor_b:(fun b -> b /. n2.@(j))
                    distance metrics t1.twisted.matrix.data.(i) t2.twisted.matrix.data.(j)) in
            (* ...and summarise them *)
            summarize_distance_matrix_row
              req_len t2.twisted.matrix.row_names.(j) distances t1.twisted.matrix.row_names buf;
            if output_distance_matrix then
              List.accum rows (j, distances)
          done;
          hi_row - lo_row + 1, Buffer.contents buf, !rows)
        (fun (n_processed, block, rows) ->
          Printf.fprintf output "%s" block;
          List.iter (fun (i, row) -> data.(i) <- row) rows;
          let new_processed_rows = !processed_rows + n_processed in
          if verbose && new_processed_rows / rows_per_step > !processed_rows / rows_per_step then
            Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows%!"
              String.TermIO.clear __FUNCTION__ path new_processed_rows r2;
          processed_rows := new_processed_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows.\n%!"
          String.TermIO.clear __FUNCTION__ path r2 r2;
      close_out output;
      if output_distance_matrix then
        Matrix.to_file ~precision ~threads ~elements_per_step ~verbose {
          which = DMatrix;
          matrix = {
            row_names = row_names_2;
            col_names = row_names_1;
            data
          }
        } prefix
    let summarize_neighbors
        ?(normalize = true) ?(how_many = Some 2) ?(policy = NeighborsPolicy.of_string "times(2)")
        ?(index_type = Interfaiss.Type.of_string "flat") ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        metric db qs prefix =
      (* Perform some compatibility checks *)
      let inertias_qs = get_inertias qs and inertias_db = get_inertias db in
      if inertias_qs <> inertias_db then
        raise_incompatible_inertias __FUNCTION__ inertias_qs inertias_db;
      let col_names_qs = qs.twisted.matrix.col_names and col_names_db = db.twisted.matrix.col_names in
      if col_names_qs <> col_names_db then
        Matrix.Exception.raise_incompatible_geometries __FUNCTION__ col_names_qs col_names_db;
      let metrics = Space.Distance.Metric.compute metric inertias_qs in
      (* Here distance is by definition euclidean *)
      let distance = Space.Distance.of_string "euclidean" and d = Array.length col_names_qs in
      let n_db = Array.length db.twisted.matrix.row_names and n_qs = Array.length qs.twisted.matrix.row_names in
      (* We compute embeddings... *)
      let embeddings_db = to_embeddings ~normalize ~elements_per_step ~threads ~verbose distance metric db
      and embeddings_qs = to_embeddings ~normalize ~elements_per_step ~threads ~verbose distance metric qs
      (* ...and copy them to bigarrays *)
      and data_db = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout n_db d
      and data_qs = Bigarray.Array2.create Bigarray.Float32 Bigarray.C_layout n_qs d in
      Matrix.Base.to_bigarray embeddings_db.matrix.data data_db;
      Matrix.Base.to_bigarray embeddings_qs.matrix.data data_qs;
      (* We generate and train the Faiss index *)
      if verbose then
        Printf.eprintf "(%s): Generating and training vector index '%s'...%!"
          __FUNCTION__ (Interfaiss.Type.to_string index_type);
      let index = Interfaiss.create ~index_type d in
      Interfaiss.train index data_db;
      Interfaiss.add index data_db;
      if verbose then
        Printf.eprintf " done.\n%!";
      (* This is a well-behaved number - cannot be more than n_db *)
      let eff_nn_number = NeighborsPolicy.get_to_be_visited policy n_db how_many in
      if verbose then
        Printf.eprintf "(%s): Querying and deleting vector index...%!" __FUNCTION__;
      let res_idxs, _ = Interfaiss.query index data_qs eff_nn_number in (* We'll recompute distances *)
      (* We delete the index *)
      Interfaiss.delete index;
      if verbose then
        Printf.eprintf " done.\n%!";
      (* We compute normalisations *)
      let norm_db, norm_qs =
        if normalize then
          Matrix.Base.get_normalizations ~threads ~elements_per_step ~verbose distance metrics db.twisted.matrix,
          Matrix.Base.get_normalizations ~threads ~elements_per_step ~verbose distance metrics qs.twisted.matrix
        else
          Float.Array.make n_db 1., Float.Array.make n_qs 1. in
      (* We compute statistics and output results *)
      let path = make_filename_summary prefix in
      let output = open_out path
      and rows_per_step = max 1 (elements_per_step / eff_nn_number) and processed_rows = ref 0
      and buf = Buffer.create 1048576 in
      (* Parallel section *)
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_qs then
            let to_do = min rows_per_step (n_qs - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          Buffer.clear buf;
          for j = lo_row to hi_row do
            let idxs = Bigarray.Array2.slice_left res_idxs j in
            (* For each query, we (re)compute the respective distances from the matches.
               As many of Faiss indices are not exhaustive, as a first thing
                we have to compute the actual number of neighbours found *)
            let eff_dim_idxs =
              let dim_idxs = Bigarray.Array1.dim idxs and res = ref 0 in
              begin try
                for i = 0 to dim_idxs - 1 do
                  if idxs.{i} = -1L then begin
                    res := i;
                    raise_notrace Exit (* This one is OK as it will be caught *)
                  end
                done;
                dim_idxs
              with Exit ->
                !res
              end in
            let distances =
              Float.Array.init eff_dim_idxs
                (fun i ->
                  let idx = Int64.to_int idxs.{i} in
                  Space.Distance.compute
                    ~adaptor_a:(fun a -> a /. norm_db.@(idx)) ~adaptor_b:(fun b -> b /. norm_qs.@(j))
                    distance metrics db.twisted.matrix.data.(idx) qs.twisted.matrix.data.(j)) in
            (* ...and summarise them *)
            summarize_distance_matrix_row ~col_names_mapping:(fun i -> Int64.to_int idxs.{i})
              eff_nn_number qs.twisted.matrix.row_names.(j) distances db.twisted.matrix.row_names buf
          done;
          hi_row - lo_row + 1, Buffer.contents buf)
        (fun (n_processed, block) ->
          Printf.fprintf output "%s" block;
          let new_processed_rows = !processed_rows + n_processed in
          if verbose && new_processed_rows / rows_per_step > !processed_rows / rows_per_step then
            Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows%!"
              String.TermIO.clear __FUNCTION__ path new_processed_rows n_qs;
          processed_rows := new_processed_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Writing distance digest to file '%s': done %d/%d rows.\n%!"
          String.TermIO.clear __FUNCTION__ path n_qs n_qs;
      close_out output
    (* Implementation module *)
    module Bipartition =
      struct
        let make
            ?(acceptance_probability_at_zero = 0.2) ?(difference_magnification_factor = 10.)
            ?(balance_penalty = BalancePenalty.default)
            ?rng_state ?(verbose = false)
            m init_set =
          (* If no RNG state is supplied, derive one deterministically from
             the input matrix's row names.  This makes runs reproducible by
             default (same input -> same partition) and data-dependent
             (different inputs use different initialisations). *)
          let rng_state =
            match rng_state with
            | Some s -> s
            | None -> Random.State.make [| Hashtbl.hash m.Matrix.Base.row_names |] in
          if acceptance_probability_at_zero <= 0. || acceptance_probability_at_zero > 1. then
            Exception.raise __FUNCTION__ Initialize
              (Printf.sprintf "Invalid acceptance probability at zero (expected float between 0. and 1., found %.16g)"
                acceptance_probability_at_zero);
          if difference_magnification_factor <= 0. then
            Exception.raise __FUNCTION__ Initialize
              (Printf.sprintf "Difference magnification factor cannot be negative (found %.16g)"
                difference_magnification_factor);
          let inverse_acceptance = (1. -. acceptance_probability_at_zero) /. acceptance_probability_at_zero
          and negative_scale = -. difference_magnification_factor
          and elements = IntSet.elements_array init_set in
          let num_elements = Array.length elements in
          if num_elements < 2 then
            Exception.raise __FUNCTION__ Initialize
              (Printf.sprintf "There must be at least two elements (found %d)" num_elements);
          let n = Array.length m.Matrix.Base.row_names
          and d = Array.length m.Matrix.Base.col_names
          and one = ref IntSet.empty and cardinal_one = ref 0
          and two = ref IntSet.empty and cardinal_two = ref 0 in
          let centroid_one = Float.Array.make d 0. |> ref
          and centroid_two = Float.Array.make d 0. |> ref in
          IntSet.iter
            (fun i ->
              if i >= n then
                Exception.raise_index_out_of_range __FUNCTION__ i "set" n;
              let v = m.data.(i) in
              (* We randomly assign the element to either set *)
              if Random.State.bool rng_state then begin
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
            (* Cardinality-imbalance penalty.  Denominator is
               (1 + |c1 - c2|^alpha)^beta; alpha = 1, beta = 0.5 reproduces
               the original sqrt(1 + |c1 - c2|) exactly. *)
            let imbalance = Float.abs (cardinal_one -. cardinal_two) in
            let denom =
              (1. +. (imbalance ** balance_penalty.BalancePenalty.alpha))
              ** balance_penalty.BalancePenalty.beta in
            !res /. denom in
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
            if verbose && !steps mod 1000 = 0 then (* Maybe remove in the future *)
              Printf.eprintf " Step #%d: objective=%.3g, max_objective=%.3g\n%!" !steps !objective !max_objective;
            incr steps;
            (* We save the old state *)
            let old_one = !one and old_cardinal_one = !cardinal_one
            and old_two = !two and old_cardinal_two = !cardinal_two
            and old_objective = !objective in
            Float.Array.blit !centroid_one 0 !old_centroid_one 0 d;
            Float.Array.blit !centroid_two 0 !old_centroid_two 0 d;
            let selected = elements.(Random.State.int rng_state num_elements) in
            let v = m.data.(selected) in
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
            if Random.State.float rng_state 1. <= score then begin
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
            IntSet.iter (fun i -> m.row_names.(i) |> Printf.eprintf " '%s'") !max_one;
            Printf.eprintf " ], two=[";
            IntSet.iter (fun i -> m.row_names.(i) |> Printf.eprintf " '%s'") !max_two;
            Printf.eprintf " ], steps=%d)\n%!" !steps
          end;
          !max_one, !max_two, !max_objective, !steps
      end
    (* Dense HDBSCAN* against the embedded coordinates *)
    let get_splits
        ?(normalize = true) ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        ?(balance_penalty = BalancePenalty.default)
        ?(gaps_prefilter_kneedle = false)
        ?(num_seeds = 1)
        ?seed
        ?(hdbscan_min_cluster_size = 5)
        ?hdbscan_min_samples
        ?(hdbscan_mst_mode = Clustering.HdbscanMstMode.Auto)
        ?hdbscan_num_neighbors
        ?hdbscan_index_type
        distance metric algorithm_type max_splits t =
      (* We compute embeddings *)
      let m = to_embeddings ~normalize ~elements_per_step ~threads ~verbose distance metric t in
      match algorithm_type with
      | SplitsAlgorithm.Gaps ->
        (* Embeddings are stored rowwise.
          We begin by sorting coordinates along each dimension (i.e., by sorting columns) *)
        let n = Array.length m.matrix.row_names in
        let cols_per_step = max 1 (elements_per_step / n) and processed_cols = ref 0
        and d = Array.length m.matrix.col_names in
        (* Per-dimension inertia weights: sigma_d = sqrt(I_d).  We multiply each
           gap width by sigma_d before the global sort so that the cross-
           dimension ranking is inertia-aware: a gap of the same physical width
           on a high-inertia dimension carries more signal than on a low-
           inertia one, and the ranking the compatibility filter sees should
           reflect that.  Within a single dimension the multiplier is a
           positive constant, so the ordering of the within-dimension gaps is
           unchanged. *)
        let inertias = get_inertias t in
        let row_permutations = Array.make d [||] and gaps = Tools.ArrayStack.empty () in
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
              (* We compute gaps, i.e. differences between consecutive coordinates,
                 each multiplied by sigma_d = sqrt(I_d) so that high-inertia
                 dimensions' gaps rank above low-inertia ones in the global
                 cross-dimension sort.  Gaps are annotated with their indices. *)
              let sigma_d = sqrt inertias.@(i) in
              let gaps__idxs = Array.init (n - 1) (fun j -> (fst coords__idxs.(j + 1) -. fst coords__idxs.(j)) *. sigma_d, i, j) in
              (* We sort the vector *)
              Array.sort (fun (gap_1, _, _) (gap_2, _, _) -> compare gap_1 gap_2) gaps__idxs;
              (* Optional per-dimension Kneedle filter on the sorted gap array:
                 keep only gaps at or above the elbow rank.  When the gap
                 distribution is flat we drop every gap (treating them all as
                 noise) -- this is the only place where the local semantics
                 diverge from Clustering.kneedle_elbow's rank-0 fallback. *)
              let gaps__idxs =
                if gaps_prefilter_kneedle && Array.length gaps__idxs >= 2 then begin
                  let get_gap idx = let (g, _, _) = gaps__idxs.(idx) in g in
                  let mm = Array.length gaps__idxs in
                  let g0 = get_gap 0 and gm = get_gap (mm - 1) in
                  if gm -. g0 = 0. then [||]
                  else begin
                    let elbow = Clustering.kneedle_elbow get_gap mm in
                    Array.sub gaps__idxs elbow (mm - elbow)
                  end
                end else
                  gaps__idxs in
              (* We return the permutation of row indices and the gap vector *)
              List.accum res (Array.init n (fun row -> snd coords__idxs.(row)), gaps__idxs)
            done;
            lo_col, !res)
          (fun (lo_col, cols) ->
            List.iteri
              (fun offs_i (perm_i, gaps_i) ->
                row_permutations.(lo_col + offs_i) <- perm_i;
                Tools.ArrayStack.push_array gaps gaps_i;
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
        let gaps = Tools.ArrayStack.contents gaps in
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
        (* Recursive divisive bipartition of the sample set, run for
           [num_seeds] independent RNG initialisations.  At each level
           we call Bipartition.make, emit the found split as a
           candidate (weighted by the objective), and recurse on each
           side.  The set of recursive splits emitted by a single seed
           is mutually compatible by construction (each is a
           refinement of its parent), so within a seed the
           compatibility filter never rejects centroid splits against
           one another.  Across seeds, coincident bipartitions get
           their weights summed by Trees.Splits.add_split, giving
           bootstrap-support semantics: a split agreed upon by all
           seeds carries roughly [num_seeds] times the weight of a
           seed-specific one.
           The base seed is either the explicit user [?seed] or a
           hash of the input row names; per-iteration RNG states are
           derived via [Random.State.make [| base_seed; i |]], so the
           recursion is deterministic per input/seed/iteration.
           The seeds are independent and are run in parallel: each
           worker builds a thread-local Trees.Splits.t, and the
           serial reducer merges those into the master [res] via
           [Trees.Splits.add_splits].  Float addition is associative
           and commutative on the canonical-split keys, so the merge
           result is independent of worker completion order --
           output is byte-identical across thread counts.  Per-call
           Bipartition.make verbose trace is suppressed (it would
           interleave unreadably); progress is reported by the
           reducer as one "Done j/K seeds" line. *)
        let res = Trees.Splits.create m.matrix.row_names in
        let base_seed =
          match seed with
          | Some n -> n
          | None -> Hashtbl.hash m.matrix.row_names in
        let full_set =
          Seq.init (Array.length m.matrix.row_names) Fun.id |> IntSet.of_seq in
        let next_seed = ref 0 and done_seeds = ref 0 in
        let centroid_threads = max 1 (min threads num_seeds) in
        Processes.Parallel.process_stream_chunkwise
          (fun () ->
            if !next_seed < num_seeds then begin
              let i = !next_seed in
              incr next_seed;
              i
            end else
              raise End_of_file)
          (fun i ->
            let rng_state = Random.State.make [| base_seed; i |] in
            let local = Trees.Splits.create m.matrix.row_names in
            let rec emit set =
              if IntSet.cardinal set > 1 then begin
                let one, two, objective, _steps =
                  Bipartition.make ~balance_penalty ~rng_state ~verbose:false m.matrix set in
                let split_arr =
                  IntSet.elements_array one |> Trees.Splits.Split.of_array in
                Trees.Splits.add_split local split_arr objective;
                emit one;
                emit two
              end in
            emit full_set;
            local)
          (fun local ->
            Trees.Splits.add_splits res local;
            incr done_seeds;
            if verbose then
              Printf.eprintf "%s\r(%s): Done %d/%d seeds%!"
                String.TermIO.clear __FUNCTION__ !done_seeds num_seeds)
          centroid_threads;
        if verbose then
          Printf.eprintf "%s\r(%s): Done %d/%d seeds.\n%!"
            String.TermIO.clear __FUNCTION__ !done_seeds num_seeds;
        res
      | Hdbscan ->
        let min_samples =
          match hdbscan_min_samples with
          | Some k -> k
          | None -> hdbscan_min_cluster_size in
        Clustering.Hdbscan.make_splits ~threads ~verbose ~mst_mode:hdbscan_mst_mode
          ?num_neighbors:hdbscan_num_neighbors
          ?index_type:hdbscan_index_type
          ~min_cluster_size:hdbscan_min_cluster_size ~min_samples m.matrix
    (* *)
    let to_files ?(precision = 15) ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) v prefix =
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose v.inertia prefix;
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose v.twisted prefix
    let of_files ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) prefix =
      let inertia = Matrix.of_file ~threads ~bytes_per_step ~verbose Inertia prefix
      and twisted = Matrix.of_file ~threads ~bytes_per_step ~verbose Twisted prefix in
      (* Let's run at least some checks *)
      if inertia.matrix.row_names <> [| "inertia" |] then
        Matrix.Exception.raise_unexpected_columns_in_inertia_file __FUNCTION__ inertia.matrix.row_names;
      if inertia.matrix.col_names <> twisted.matrix.col_names then
        Matrix.Exception.raise_incompatible_geometries __FUNCTION__ inertia.matrix.col_names twisted.matrix.col_names;
      { inertia; twisted }
    (* *)
    let archive_version = "2025-10-08"
    (* *)
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopTwisted"
    let to_binary ?(verbose = false) t prefix =
      let path = make_filename_binary prefix in
      let output = open_out path in
      if verbose then
        Printf.eprintf "(%s): Outputting vectors to file '%s'...%!" __FUNCTION__ path;
      output_value output "KPopTwisted";
      output_value output archive_version;
      Matrix.to_channel output t.inertia;
      Matrix.to_channel output t.twisted;
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) prefix =
      let path = make_filename_binary prefix in
      let input = open_in path in
      if verbose then
        Printf.eprintf "(%s): Reading vectors from file '%s'...%!" __FUNCTION__ path;
      let which = (input_value input: string) in
      let version = (input_value input: string) in
      if which <> "KPopTwisted" || version <> archive_version then
        Matrix.Exception.raise_incompatible_archive_version __FUNCTION__ which version;
      let inertia = Matrix.of_channel input in
      let twisted = Matrix.of_channel input in
      close_in input;
      if Matrix.Type.Inertia <> inertia.which then
        Matrix.Exception.raise_unexpected_type __FUNCTION__ Inertia inertia.which;
      if Matrix.Type.Twisted <> twisted.which then
        Matrix.Exception.raise_unexpected_type __FUNCTION__ Twisted twisted.which;
      if verbose then
        Printf.eprintf " done.\n%!";
      { inertia; twisted }
  end: sig
    type t = {
      (* The original inertia as computed by CA.
         This is a matrix of type Inertia *)
      inertia: Matrix.t;
      (* The twisted spectra.
         They are a matrix of type Twisted *)
      twisted: Matrix.t
    }
    val empty: t
    (* It fails if matrices have incompatible intertias *)
    val merge_rowwise: t -> t -> t
    (* Recompute vectors from the specified distance and metric functions *)
    val to_embeddings: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                       Space.Distance.t -> Space.Distance.Metric.t -> t -> Matrix.t
    (* Compute distances between the rows of two matrices and summarise them.
       It fails when matrices are incompatible *)
    val summarize_distances_rowwise: ?normalize:bool -> ?keep_at_most:int option ->
                                     ?output_distance_matrix:bool -> ?precision:int ->
                                     ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                                     Space.Distance.t -> Space.Distance.Metric.t -> t -> t -> string -> unit
    (* Find nearest neighbours of each row of a matrix among the rows of another matrix
        using Euclidean distance and the specified metric function, and summarise such neighbours
        while considering more elements according to the specified policy.
       It fails when matrices are incompatible *)
    val summarize_neighbors: ?normalize:bool -> ?how_many:int option ->
                             ?policy:NeighborsPolicy.t -> ?index_type:Interfaiss.Type.t ->
                             ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                             Space.Distance.Metric.t -> t -> t -> string -> unit
    (* Output splits for the vectors computed with the specified distance and metric functions *)
    val get_splits: ?normalize:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                    ?balance_penalty:BalancePenalty.t ->
                    ?gaps_prefilter_kneedle:bool ->
                    ?num_seeds:int ->
                    ?seed:int ->
                    ?hdbscan_min_cluster_size:int ->
                    ?hdbscan_min_samples:int ->
                    ?hdbscan_mst_mode:Clustering.HdbscanMstMode.t ->
                    ?hdbscan_num_neighbors:int ->
                    ?hdbscan_index_type:Interfaiss.Type.t ->
                    Space.Distance.t -> Space.Distance.Metric.t ->
                    SplitsAlgorithm.t -> int -> t -> Trees.Splits.t
    (* Input/Output *)
    val to_files: ?precision:int -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> unit
    val of_files: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t
  end
)

