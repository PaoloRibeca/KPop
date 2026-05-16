(*
    KMerDB.ml -- (c) 2022-2025 Paolo Ribeca, <paolo.ribeca@gmail.com>

    This file is part of KPop, a scalable method for comparative analysis
    of microbial genomes and environmental samples based on full k-mer
    spectra and correspondence analysis (CA).

    KMerDB.ml is the top-level wrapper around `KMerDB_Base`.  It exposes
    file I/O for the binary `.KPopSpectra` format, the
    metadata-augmentation operations consumed by `KPopCountDB`, register
    management, and the higher-level queries used by the KPop binaries.

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
(*module Numbers = BiOCamLib.Numbers
module Processes = BiOCamLib.Processes
module Tools = BiOCamLib.Tools*)
open BiOCamLib.Better

module Transformation:
  sig
    module Statistics:
      sig
        type t = {
          non_zero: int;
          min: int;
          max: int;
          sum: float;
          sum_log: float
        }
        type table_t = {
          col_stats: t array;
          row_stats: t array
        }
      end
    (* Transformation function *)
    type t =
      (* Thresholds in the interval (0.,1.) are taken as relative ones *)
      | Threshold of float
      | Power of float
      | CLR
      | Pseudo of float * bool (* The boolean is whether to quantise *)
    val of_string: string -> t
    val to_string: t -> string
    val compute: which:t -> col_num:int -> col_stats:Statistics.t -> row_num:int -> row_stats:Statistics.t ->
                 float -> float
  end
= struct
    module Statistics =
      struct
        type t = {
          non_zero: int;
          min: int;
          max: int;
          sum: float;
          sum_log: float
        }
        type table_t = {
          col_stats: t array;
          row_stats: t array
        }
      end
    type t =
      | Threshold of float
      | Power of float
      | CLR
      | Pseudo of float * bool
    let raise __FUNCTION__ kind s =
      Exception.raise __FUNCTION__ IO_Format (Printf.sprintf "%s transformation '%s'" kind s)
    let of_string_re = Str.regexp "[(,)]"
    let of_string s =
      let raise kind = raise __FUNCTION__ kind s in
      try
        match Str.full_split of_string_re s with
        | [ Text "threshold"; Delim "("; Text threshold; Delim ")" ] ->
          Threshold (float_of_string threshold)
        | [ Text "pow"; Delim "("; Text power; Delim ")" ]
        | [ Text "power"; Delim "("; Text power; Delim ")" ] ->
          Power (float_of_string power)
        | [ Text "binary" ] ->
          Power 0.
        | [ Text "clr" ] | [ Text "CLR" ] ->
          CLR
        | [ Text "pseudo"; Delim "("; Text power; Delim ","; Text quantize; Delim ")" ]
        | [ Text "pseudocounts"; Delim "("; Text power; Delim ","; Text quantize; Delim ")" ] ->
          let power = float_of_string power in
          let res = Pseudo (power, bool_of_string quantize) in
          if power < 0. then
            raise "Invalid";
          res
        | _ ->
          raise "Unknown"
      with _ ->
        raise "Unknown"
    let to_string = function
      | Threshold threshold ->
        Printf.sprintf "threshold(%g)" threshold
      | Power power ->
        Printf.sprintf "power(%g)" power
      | CLR ->
        "clr"
      | Pseudo (power, quantize) ->
        Printf.sprintf "pseudocounts(%g,%b)" power quantize
    let [@warning "-27"] compute ~which ~col_num ~col_stats ~row_num ~row_stats counts =
      match which with
      | Threshold threshold ->
        let threshold =
          if threshold < 1. then
            threshold *. col_stats.Statistics.sum
          else
            threshold in
        if counts >= threshold then
          counts
        else
          0.
      | Power 1. -> (* Optimisation *)
        counts
      | Power power ->
        counts ** power
      | CLR ->
        (float_of_int col_stats.max +. 1.) *. log (counts +. 1.)
      | Pseudo (power, quantize) ->
        if power < 0. then
          to_string which |> raise __FUNCTION__ "Invalid";
        let v =
          if power = 0. then
            (float_of_int col_stats.max +. 1.) *. log (counts +. 1.)
          else begin
            if power < 1. then
              (counts ** power) *. ((float_of_int col_stats.max) ** (1. -. power)) /. power
            else
              (counts ** power)
          end in
        if quantize then
          floor v
        else
          v
  end

module CombinationCriterion:
  sig
    type t = RescaledMean | RescaledMedian
    val of_string: string -> t
    val to_string: t -> string
  end
= struct
    type t = RescaledMean | RescaledMedian
    let of_string = function
      | "mean" -> RescaledMean
      | "median" -> RescaledMedian
      | s ->
        Exception.raise_unrecognized_initializer __FUNCTION__ "Combination criterion" s
    let to_string = function
      | RescaledMean -> "mean"
      | RescaledMedian -> "median"
  end

include (
  struct
    include KMerDB_Base
    let ( .@() ) = CountBAVector.( .@() )
    let ( .@()<- ) = CountBAVector.( .@()<- )
    (* Helper module *)
    module ColOrRow =
      struct
        type t =
        | Col
        | Row
        let to_string = function
        | Col -> "column"
        | Row -> "row"
      end
    (* Implementation function *)
    let stats_table_of_core_db ?(threads = 1) ?(verbose = false) core =
      let compute_one what n =
        let red_len =
          match what with
          | ColOrRow.Col ->
            core.n_rows - 1
          | Row ->
            core.n_cols - 1 in
        let non_zero = ref 0 and min = ref 0 and max = ref 0 and sum = ref 0. and sum_log = ref 0. in
        for i = 0 to red_len do
          let v =
            match what with
            | Col ->
              CountBAVector.N.to_int core.data.(n).@(i)
            | Row ->
              CountBAVector.N.to_int core.data.(i).@(n) in
          let f_v = float_of_int v in
          if f_v >= 0. then
            incr non_zero;
          min := Stdlib.min !min v;
          max := Stdlib.max !max v;
          sum := !sum +. f_v;
          sum_log := !sum_log +. log f_v
        done;
        { Transformation.Statistics.non_zero = !non_zero;
          min = !min;
          max = !max;
          sum = !sum;
          sum_log = !sum_log } in
      let compute_all what =
        let n =
          match what with
          | ColOrRow.Col ->
            core.n_cols
          | Row ->
            core.n_rows in
        let step = n / threads / 5 |> max 1 and processed = ref 0 and res = ref [] in
        Processes.Parallel.process_stream_chunkwise
          (fun () ->
            if verbose then
              Printf.eprintf "%s\r(%s): Computing %s statistics [%d/%d]%!"
                String.TermIO.clear __FUNCTION__ (ColOrRow.to_string what) !processed n;
            let to_do = n - !processed in
            if to_do > 0 then begin
              let to_do = min to_do step in
              let res = !processed, to_do in
              processed := !processed + to_do;
              res
            end else
              raise End_of_file)
          (fun (processed, to_do) ->
            let res = ref [] in
            for i = 0 to to_do - 1 do
              processed + i |> compute_one what |> List.accum res
            done;
            Array.of_rlist !res)
          (List.accum res)
          threads;
        let rec binary_merge_arrays processed to_do =
          match processed, to_do with
          | [], [] ->
            [||]
          | [], [a] ->
            a
          | [res], [] ->
            res
          | [res], [a] ->
            Array.append res a
          | _, [] ->
            List.rev processed |> binary_merge_arrays []
          | _, [a] ->
            a :: processed |> List.rev |> binary_merge_arrays []
          | _, a1 :: a2 :: tl ->
            binary_merge_arrays ((Array.append a1 a2) :: processed) tl in
        let res = List.rev !res |> binary_merge_arrays [] in
        if verbose then
          Printf.eprintf "%s\r(%s): Computing %s statistics [%d/%d]\n%!"
            String.TermIO.clear __FUNCTION__ (ColOrRow.to_string what) n n;
        res in
      { Transformation.Statistics.col_stats = compute_all Col;
        row_stats = compute_all Row }
    let transform ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) transformation db =
      let transform = Transformation.compute ~which:transformation
      and stats = stats_table_of_core_db ~threads ~verbose db.core
      and n_rows = db.core.n_rows and n_cols = db.core.n_cols
      and new_data = Array.copy db.core.data and processed_cols = ref 0 in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_cols < n_cols then
            let to_do = max 1 (elements_per_step / n_rows) |> min (n_cols - !processed_cols) in
            let new_processed_cols = !processed_cols + to_do in
            let res = !processed_cols, new_processed_cols - 1 in
            processed_cols := new_processed_cols;
            res
          else
            raise End_of_file)
        (fun (lo_col, hi_col) ->
          let res = ref [] in
          for col_idx = lo_col to hi_col do
            List.accum res begin
              col_idx,
              CountBAVector.mapi
                (fun row_idx n ->
                  transform
                    ~col_num:n_cols ~col_stats:stats.col_stats.(col_idx)
                    ~row_num:n_rows ~row_stats:stats.row_stats.(row_idx) n)
                db.core.data.(col_idx)
            end
          done;
          !res)
        (List.iter
          (fun (col_idx, col) ->
            new_data.(col_idx) <- col;
            incr processed_cols;
            if verbose then
              Printf.eprintf "%s\r(%s): Transforming database: done %d/%d %s%!"
                String.TermIO.clear __FUNCTION__ !processed_cols n_cols
                (String.pluralize_int "column" !processed_cols)))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Transforming database: done %d/%d %s.\n%!"
          String.TermIO.clear __FUNCTION__ n_cols n_cols
          (String.pluralize_int "column" n_cols);
      { db with core = { db.core with data = new_data } }
    (* It should be OK to have the same label on both LHS and RHS,
        as a temporary vector is used to generate the combination *)
    let add_combined_selected ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        db new_label selection criterion =
      let stats = stats_table_of_core_db ~threads ~verbose db.core and db = ref db in
      if verbose then
        Printf.eprintf "(%s): Adding/replacing spectrum '%s': [%!" __FUNCTION__ new_label;
      (* We allocate the result *)
      let new_col_idx = add_empty_column_if_needed db new_label in
      let new_col = !db.core.data.(new_col_idx) in
      (* Computing valid labels and maximum normalisation across columns *)
      let found_cols = ref [] and max_norm = ref 0. in
      StringSet.iter
        (fun label ->
          if verbose then
            Printf.eprintf " '%s'%!" label;
          (* Some labels might be invalid *)
          match StringHashtbl.find_opt !db.col_names_to_idx label with
          | Some col_idx ->
            List.accum found_cols col_idx;
            max_norm := max !max_norm stats.col_stats.(col_idx).sum
          | None ->
            if verbose then
              Printf.eprintf "(NOT FOUND)%!")
        selection;
      let found_cols = Array.of_list !found_cols in (* We don't really care about the order here *)
      let num_found_cols = Array.length found_cols and max_norm = !max_norm and norm = ref 0. in
      if verbose then
        Printf.eprintf " ] n_found=%d max_norm=%.16g.\n%!" num_found_cols max_norm;
      let n_rows = !db.core.n_rows and processed_rows = ref 0 in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_rows then
            let to_do = max 1 (elements_per_step / num_found_cols) |> min (n_rows - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          let module Freqs = Numbers.FloatFreqsVector in
          (* We need one histogram per row to be able to compute statistics such as the median *)
          let row_combinators = Array.init (hi_row - lo_row + 1) (fun _ -> Freqs.make ~non_negative:true ()) in
          for i = lo_row to hi_row do
            (* We iterate over valid columns *)
            Array.iter
              (fun col_idx ->
                (* We normalise columns *separately* before combining them *)
                let col = !db.core.data.(col_idx) and norm = stats.col_stats.(col_idx).sum in
                (* All counts are non-negative *)
                if norm > 0. then
                  (* We add the renormalised sum to the suitable row histogram *)
                  col.@(i) *. max_norm /. norm |> Freqs.add row_combinators.(i - lo_row))
              found_cols
          done;
          (* For each row histogram in the input range, we now generate a combination and pass it along *)
          lo_row,
          Array.map
            (fun combinator ->
              match criterion with
              | CombinationCriterion.RescaledMean ->
                Freqs.sum combinator
              | RescaledMedian ->
                Freqs.median combinator *. float_of_int num_found_cols)
            row_combinators)
        (fun (lo_row, block) ->
          let n_processed = Array.length block in
          for i = lo_row to lo_row + n_processed - 1 do
            let res_i = block.(i - lo_row) in
            Numbers.Float.(norm ++ res_i);
            (* Actual copy to data *)
            new_col.@(i) <- res_i
          done;
          let old_processed_rows = !processed_rows in
          processed_rows := !processed_rows + n_processed;
          if verbose && !processed_rows / 10000 > old_processed_rows / 10000 then
            Printf.eprintf "%s\r(%s): Combining spectra: done %d/%d lines%!"
              String.TermIO.clear __FUNCTION__ !processed_rows n_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Combining spectra: done %d/%d lines. Norm=%.16g\n%!"
          String.TermIO.clear __FUNCTION__ !processed_rows n_rows !norm;
      (* If metadata is present in the database, we generate some for the new column too *)
      if !db.core.n_meta > 0 then begin
        (* For each metadata field, we compute the intersection of the values across all selected columns *)
        let res = Array.make !db.core.n_meta StringSet.empty in
        StringSet.iter
          (fun label ->
            match StringHashtbl.find_opt !db.col_names_to_idx label with
            | Some col_idx ->
              let col = !db.core.meta.(col_idx) in
              for i = 0 to !db.core.n_meta - 1 do
                res.(i) <- StringSet.add col.(i) res.(i)
              done
            | None -> ())
          selection;
        !db.core.meta.(new_col_idx) <-
          Array.map
            (fun set ->
              if StringSet.cardinal set = 1 then
                StringSet.min_elt set
              else
                "")
            res
      end;
      !db
    let get_indicator_vector db classes_label =
      match StringHashtbl.find_opt db.meta_names_to_idx classes_label with
      | None ->
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Could not find specified metadata field '%s'" classes_label)
      | Some classes_label_idx ->
        (* We derive the class indicator vector *)
        let n_samples = db.core.n_cols in
        let num_classes = ref 0 and class_to_ind = ref StringMap.empty and ind_to_class = ref IntMap.empty
        and res = Array.make n_samples (-1) in
        (* We iterate explicitly as there might be trailing space at the end of vectors *)
        for i = 0 to n_samples - 1 do
          res.(i) <- begin
            let class_label = db.core.meta.(i).(classes_label_idx) in
            match StringMap.find_opt class_label !class_to_ind with
            | None ->
              (* We insert the new label *)
              let res = !num_classes in
              class_to_ind := StringMap.add class_label res !class_to_ind;
              ind_to_class := IntMap.add res class_label !ind_to_class;
              incr num_classes;
              res
            | Some ind ->
              ind
          end
        done;
        !num_classes, !ind_to_class, res
    let split_spectra ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        db classes_label criterion =
      let _, ind_to_class, ind_classes = get_indicator_vector db classes_label in
      if verbose then begin
        Printf.eprintf "(%s): Classes=[" __FUNCTION__;
        Array.iter (Printf.eprintf " %d") ind_classes;
        Printf.eprintf " ]\n%!"
      end;
      (* We split spectra names according to their class *)
      let module I2SMM = Tools.Multimap (ComparableInt) (ComparableString) in
      let split_names = ref I2SMM.empty in
      Array.iteri
        (fun i ind ->
          split_names := I2SMM.add ind db.core.idx_to_col_names.(i) !split_names)
        ind_classes;
      let res = ref db in
      I2SMM.iter_set
        (fun ind spectra_names ->
          let class_name = IntMap.find ind ind_to_class in
          if StringHashtbl.mem db.col_names_to_idx class_name then
            Exception.raise __FUNCTION__ IO_Format
              (Printf.sprintf "Specified metadata field '%s' is also the name of a spectrum in the database"
                class_name);
          res := add_combined_selected ~threads ~elements_per_step ~verbose !res class_name spectra_names criterion)
        !split_names;
      remove_selected !res (Array.to_seq db.core.idx_to_col_names |> StringSet.of_seq)
    (* *)
    module Stats = Numbers.OnlineStats(Numbers.Float) (*(CountBAVector.N)*)
    module Freqs = Numbers.FloatFreqsVector
    module FAVector = Numbers.FAVector
    module LF_FAVector = Numbers.LinearFit(FAVector)
    let distill_kmers ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        db classes_label summary_prefix =
      let stats_table = stats_table_of_core_db ~threads ~verbose db.core in
      let n_classes, _, ind_classes = get_indicator_vector db classes_label and n_samples = db.core.n_cols in
      if n_classes = 1 || n_classes = n_samples then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf
            "Invalid number of equivalence classes identified by metadata field '%s' (found %d, cannot be 1 or %d)"
            classes_label n_classes n_samples);
      if verbose then begin
        Printf.eprintf "(%s): Classes=[" __FUNCTION__;
        Array.iter (Printf.eprintf " %d") ind_classes;
        Printf.eprintf " ]\n%!"
      end;
      let n_kmers = db.core.n_rows and processed_kmers = ref 0 in
      let stats_classes = Array.init n_classes (fun _ -> Array.init n_classes (fun _ -> Stats.make ()))
      and avgs_on_mean = FAVector.(make n_kmers N.zero) and avgs_off_mean = FAVector.(make n_kmers N.zero)
      and avgs_on_median = FAVector.(make n_kmers N.zero) and avgs_off_median = FAVector.(make n_kmers N.zero)
      and vars_on_mean = FAVector.(make n_kmers N.zero) and vars_off_mean = FAVector.(make n_kmers N.zero)
      and vars_on_median = FAVector.(make n_kmers N.zero) and vars_off_median = FAVector.(make n_kmers N.zero)
      and covs_on_mean = FAVector.(make n_kmers N.zero) and covs_off_mean = FAVector.(make n_kmers N.zero)
      and covs_on_median = FAVector.(make n_kmers N.zero) and covs_off_median = FAVector.(make n_kmers N.zero) in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_kmers < n_kmers then
            let to_do =
              max 1 (elements_per_step / (n_samples * (n_samples - 1) / 2)) |> min (n_kmers - !processed_kmers) in
            let new_processed_kmers = !processed_kmers + to_do in
            let res = !processed_kmers, new_processed_kmers - 1 in
            processed_kmers := new_processed_kmers;
            res
          else
            raise End_of_file)
        (fun (lo_kmer, hi_kmer) ->
          let res = ref [] in
          for kmer = lo_kmer to hi_kmer do
            CountBAVector.(
              (* We iterate over all the couples of samples *)
              for i = 0 to n_samples - 1 do
                let ind_class_i = ind_classes.(i)
                (* Counts must be normalised *)
                and counts_i = N.to_float db.core.data.(i).@(kmer) /. stats_table.col_stats.(i).sum in
                for j = i + 1 to n_samples - 1 do
                  let ind_class_j = ind_classes.(j)
                  (* Counts must be normalised *)
                  and counts_j = N.to_float db.core.data.(j).@(kmer) /. stats_table.col_stats.(j).sum in
                  let diff = Float.abs (counts_i -. counts_j) in
                  if ind_class_i = ind_class_j then
                    Stats.add stats_classes.(ind_class_i).(ind_class_j) diff
                  else begin
                    if ind_class_i < ind_class_j then
                      Stats.add stats_classes.(ind_class_i).(ind_class_j) diff
                    else
                      Stats.add stats_classes.(ind_class_j).(ind_class_i) diff
                  end
                done
              done;
              let avgs_on = Freqs.make () and avgs_off = Freqs.make ()
              and vars_on = Freqs.make () and vars_off = Freqs.make ()
              and covs_on = Freqs.make () and covs_off = Freqs.make () in
              (* We iterate over all the couples of classes *)
              for i = 0 to n_classes - 1 do
                let stats = stats_classes.(i).(i) in
                Stats.mean stats |> Freqs.add avgs_on;
                Stats.sample_variance stats |> Freqs.add vars_on;
                Stats.sample_coefficient_of_variation stats |> Freqs.add covs_on;
                Stats.clear stats;
                for j = i + 1 to n_classes - 1 do
                  let stats = stats_classes.(i).(j) in
                  Stats.mean stats |> Freqs.add avgs_off;
                  Stats.sample_variance stats |> Freqs.add vars_off;
                  Stats.sample_coefficient_of_variation stats |> Freqs.add covs_off;
                  Stats.clear stats
                done
              done;
              begin
                kmer, Freqs.mean avgs_on, Freqs.mean avgs_off, Freqs.median avgs_on, Freqs.median avgs_off,
                      Freqs.mean vars_on, Freqs.mean vars_off, Freqs.median vars_on, Freqs.median vars_off,
                      Freqs.mean covs_on, Freqs.mean covs_off, Freqs.median covs_on, Freqs.median covs_off
              end |> List.accum res
            );
          done;
          !res)
        (fun block ->
          let old_processed_kmers = !processed_kmers in
          List.iter
            (fun
              (kmer,
               avgs_on_mean_kmer, avgs_off_mean_kmer, avgs_on_median_kmer, avgs_off_median_kmer,
               vars_on_mean_kmer, vars_off_mean_kmer, vars_on_median_kmer, vars_off_median_kmer,
               covs_on_mean_kmer, covs_off_mean_kmer, covs_on_median_kmer, covs_off_median_kmer) ->
              FAVector.(
                avgs_on_mean.@(kmer) <- avgs_on_mean_kmer;
                avgs_off_mean.@(kmer) <- avgs_off_mean_kmer;
                avgs_on_median.@(kmer) <- avgs_on_median_kmer;
                avgs_off_median.@(kmer) <- avgs_off_median_kmer;
                vars_on_mean.@(kmer) <- vars_on_mean_kmer;
                vars_off_mean.@(kmer) <- vars_off_mean_kmer;
                vars_on_median.@(kmer) <- vars_on_median_kmer;
                vars_off_median.@(kmer) <- vars_off_median_kmer;
                covs_on_mean.@(kmer) <- covs_on_mean_kmer;
                covs_off_mean.@(kmer) <- covs_off_mean_kmer;
                covs_on_median.@(kmer) <- covs_on_median_kmer;
                covs_off_median.@(kmer) <- covs_off_median_kmer
              );
              incr processed_kmers
          ) block;
          if verbose && !processed_kmers / 100 > old_processed_kmers / 100 then
            Printf.eprintf "%s\r(%s): Distilled %d/%d kmers%!"
              String.TermIO.clear __FUNCTION__ !processed_kmers n_kmers)
        threads;
      let fit_avgs_mean, _, residuals_avgs_mean = LF_FAVector.make avgs_on_mean avgs_off_mean
      and fit_avgs_median, _, residuals_avgs_median = LF_FAVector.make avgs_on_median avgs_off_median
      and fit_vars_mean, _, residuals_vars_mean = LF_FAVector.make vars_on_mean vars_off_mean
      and fit_vars_median, _, residuals_vars_median = LF_FAVector.make vars_on_median vars_off_median
      and fit_covs_mean, _, residuals_covs_mean = LF_FAVector.make covs_on_mean covs_off_mean
      and fit_covs_median, _, residuals_covs_median = LF_FAVector.make covs_on_median covs_off_median in
      if verbose then begin
        Printf.eprintf "%s\r(%s): Distilled %d/%d kmers.\n"
          String.TermIO.clear __FUNCTION__ !processed_kmers n_kmers;
        Printf.eprintf "(%s): Fit for avgs mean is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_avgs_mean) (LF_FAVector.get_slope fit_avgs_mean);
        Printf.eprintf "(%s): Fit for avgs median is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_avgs_median) (LF_FAVector.get_slope fit_avgs_median);
        Printf.eprintf "(%s): Fit for vars mean is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_vars_mean) (LF_FAVector.get_slope fit_vars_mean);
        Printf.eprintf "(%s): Fit for vars median is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_vars_median) (LF_FAVector.get_slope fit_vars_median);
        Printf.eprintf "(%s): Fit for covs mean is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_covs_mean) (LF_FAVector.get_slope fit_covs_mean);
        Printf.eprintf "(%s): Fit for covs median is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_covs_median) (LF_FAVector.get_slope fit_covs_median)
      end;
      (* We output the summary *)
      let summary = {
        Matrix.which = Distill;
        matrix = {
          col_names = Array.sub db.core.idx_to_row_names 0 db.core.n_rows;
          row_names =
            [| "InnerAvgMean"; "OuterAvgMean"; "ResidualAvgMean";
               "InnerAvgMedian"; "OuterAvgMedian"; "ResidualAvgMedian";
               "InnerVarMean"; "OuterVarMean"; "ResidualVarMean";
               "InnerVarMedian"; "OuterVarMedian"; "ResidualVarMedian";
               "InnerCOVMean"; "OuterCOVMean"; "ResidualCOVMean";
               "InnerCOVMedian"; "OuterCOVMedian"; "ResidualCOVMedian" |];
          data = FAVector.(
            [| to_floatarray avgs_on_mean; to_floatarray avgs_off_mean; to_floatarray residuals_avgs_mean;
               to_floatarray avgs_on_median; to_floatarray avgs_off_median; to_floatarray residuals_avgs_median;
               to_floatarray vars_on_mean; to_floatarray vars_off_mean; to_floatarray residuals_vars_mean;
               to_floatarray vars_on_median; to_floatarray vars_off_median; to_floatarray residuals_vars_median;
               to_floatarray covs_on_mean; to_floatarray covs_off_mean; to_floatarray residuals_covs_mean;
               to_floatarray covs_on_median; to_floatarray covs_off_median; to_floatarray residuals_covs_median |]
          )
        }
      } in
      Matrix.to_file ~threads (Matrix.transpose ~threads summary) summary_prefix
    let to_distances
        ?(precision = 15) ?(normalise = true) ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false)
        distance db selection_1 selection_2 prefix =
      let stats = stats_table_of_core_db ~threads ~verbose db.core
      and n_r = db.core.n_rows in (* Number of k-mers. It stays the same even if we select a subset of spectra *)
      let make_submatrix selection =
        (* We select spectra, i.e. columns *)
        let idxs = ref IntSet.empty in
        Array.iteri
          (fun col_idx col_name ->
            if StringSet.mem col_name selection then
              idxs := IntSet.add col_idx !idxs)
          db.core.idx_to_col_names;
        let idxs = IntSet.elements_array !idxs in
        let n_c = Array.length idxs in (* We now know the number of columns as well *)
        (* The distance matrix is computed rowwise, and k-mers are physically stored as rows in db:
            we need to transpose *)
        { Matrix.Base.col_names = Array.sub db.core.idx_to_row_names 0 n_r;
          row_names = Array.init n_c (fun i -> db.core.idx_to_col_names.(idxs.(i)));
          data =
            (* Here we just need to convert the counts to floats, and possibly normalise *)
            Array.init n_c
              (fun i ->
                let idx = idxs.(i) in
                let norm = stats.col_stats.(idx).sum in
                let norm =
                  if not normalise || norm = 0. then
                    1.
                  else
                    norm in
                (*
                Printf.eprintf "Norm@%d=%g\n%!" idx norm;
                *)
                Float.Array.init n_r
                  (fun j -> CountBAVector.N.to_float db.core.data.(idx).@(j) /. norm)) } in
      let metric = Float.Array.make n_r 1.
      and matrix_1 = make_submatrix selection_1 and matrix_2 = make_submatrix selection_2 in
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose {
        which = DMatrix;
        matrix =
          Matrix.Base.get_distance_rowwise ~threads ~elements_per_step ~verbose distance metric matrix_1 matrix_2
      } prefix
    (* The following functions implement tabular I/O.
       Due to historical reasons (mostly the fact that the CA implementation in R
        wants the table that way) these are laid out as the _transpose_ of the internal
        data, which complicates things a bit when reading things.
       This might disappear in the future, in particular if we switch to a native
        implementation of CA.
       Note that while the internal order of metadata labels, row labels and column names
        is determined by the order in which they get added to the database, we output them
        sorted to make things more comparable. For the time being, this reordering also
        gets automatically transmitted to CA, as tabular output is what R uses as input *)
    let make_filename_meta_table = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopMMatrix.txt"
    let make_filename_kmer_table = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopKMatrix.txt"
    let to_files ?(precision = 15) ?(output_zero_kmers = false)
                 ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) db prefix =
      (* Needed to understand whether to print k-mers or not *)
      let stats = stats_table_of_core_db ~threads ~verbose db.core
      and meta = ref [] and rows = ref [] and cols = ref [] in
      (* We determine which rows and colunms should be output after all filters have been applied *)
      (*  Rows: metadata and k-mers *)
      Array.iteri
        (fun i meta_name ->
          (* There might be additional storage *)
          if i < db.core.n_meta then
            List.accum meta (meta_name, i))
        db.core.idx_to_meta_names;
      Array.iteri
        (fun i row_name ->
          (* There might be additional storage.
             Also, we only print the row if it has non-zero elements or if we are explicitly requested to do so *)
          if i < db.core.n_rows && (stats.row_stats.(i).sum > 0. || output_zero_kmers) then
            List.accum rows (row_name, i))
        db.core.idx_to_row_names;
      (*  Columns *)
      Array.iteri
        (fun i col_name ->
          (* There might be additional storage *)
          if i < db.core.n_cols then
            List.accum cols (col_name, i))
        db.core.idx_to_col_names;
      (* We want to re-sort labels based on lexicographic order *)
      let sort a =
        let l = Array.length a and m = ref StringMap.empty in
        Array.iter
          (fun (name, src_i) ->
            m := StringMap.add name src_i !m)
          a;
        let sorted = Array.make l "" and permutation = Array.make l (-1) in
        StringMap.iteri
          (fun tgt_i label src_i ->
            sorted.(tgt_i) <- label;
            permutation.(tgt_i) <- src_i)
          !m;
        sorted, permutation in
      let meta, meta_perm = Array.of_rlist !meta |> sort
      and rows, rows_perm = Array.of_rlist !rows |> sort
      and cols, cols_perm = Array.of_rlist !cols |> sort in
      let n_rows = Array.length rows and n_meta = Array.length meta and n_cols = Array.length cols in
      (* There must be at least one metadata field label and one column label to print *)
      if n_meta > 0 && n_cols > 0 then begin
        let module MetaIO = Matrix.Base.IO.Make (
          struct
            module Element =
              struct
                type t = string
                let of_string _ = assert false
                let [@warning "-27"] to_string ?(precision = 15) s = s [@@inline]
              end
            type t = string array
            let create _ = assert false
            let set _ _ _ = assert false
            let get_col_name _ i = cols.(i) [@@inline]
            let get_row_name _ i = meta.(i) [@@inline]
            let get_datum _ i j =
              (* Remember that we are writing out the _transposed_ matrix *)
              db.core.meta.(cols_perm.(j)).(meta_perm.(i))
          end
        ) in
        let output = make_filename_meta_table prefix |> open_out in
        MetaIO.to_channel ~precision ~threads ~elements_per_step ~verbose (MetaIO.Virtual (n_meta, n_cols)) output;
        close_out output
      end;
      (* There must be at least one row label and one column label to print *)
      if n_rows > 0 && n_cols > 0 then begin
        let module KmerIO = Matrix.Base.IO.Make (
          struct
            module Element =
              struct
                type t = float
                let of_string _ = assert false
                let to_string ?(precision = 15) = Printf.sprintf "%.*g" precision
              end
            type t = CountBAVector.t
            let create _ = assert false
            let set _ _ _ = assert false
            let get_col_name _ i = cols.(i) [@@inline]
            let get_row_name _ i = rows.(i) [@@inline]
            let get_datum _ i j =
              (* Remember that we are writing out the _transposed_ matrix *)
              db.core.data.(cols_perm.(j)).@(rows_perm.(i))
          end
        ) in
        let output = make_filename_kmer_table prefix |> open_out in
        KmerIO.to_channel ~precision ~threads ~elements_per_step ~verbose (KmerIO.Virtual (n_rows, n_cols)) output;
        close_out output
      end
    let of_files ?(threads = 1) ?(bytes_per_step = 4194304) ?(elements_per_step = 10000) ?(verbose = false) prefix =
      let module MetaIO = Matrix.Base.IO.Make (
        struct
          module Element =
            struct
              type t = string
              let of_string s = s [@@inline]
              let [@warning "-27"] to_string ?(precision = 15) _ = assert false
            end
          type t = string array
          let create n = Array.make n ""
          let set = Array.set
          let get_col_name _ _ = assert false
          let get_row_name _ _ = assert false
          let get_datum m i j = m.(i).(j) [@@inline] (* Needed by transpose() *)
        end
      ) in
      let input = make_filename_meta_table prefix |> open_in in
      let meta = MetaIO.of_channel ~threads ~bytes_per_step ~verbose input in
      close_in input;
      let module KmerIO = Matrix.Base.IO.Make (
        struct
          module Element =
            struct
              type t = float
              let of_string = float_of_string
              let [@warning "-27"] to_string ?(precision = 15) _ = assert false
            end
          type t = CountBAVector.t
          let create n = CountBAVector.make n CountBAVector.N.zero
          let set = CountBAVector.set
          let get_col_name _ _ = assert false
          let get_row_name _ _ = assert false
          let get_datum m i j = m.(i).@(j) [@@inline] (* Needed by transpose() *)
        end
      ) in
      let input = make_filename_kmer_table prefix |> open_in in
      let kmer = KmerIO.of_channel ~threads ~bytes_per_step ~verbose input in
      close_in input;
      (* Remember that the matrices we have just read are the _transposed_ of what
          we need to store *)
      let meta = MetaIO.transpose ~threads ~elements_per_step ~verbose meta
      and kmer = KmerIO.transpose ~threads ~elements_per_step ~verbose kmer in
      (* Perform some compatibility check. Row names must be the same, or absent *)
      if meta.row_names <> kmer.row_names && meta.row_names <> [||] && kmer.row_names <> [||] then
        Matrix.Exception.raise_incompatible_geometries __FUNCTION__
          meta.row_names kmer.row_names;
      of_core {
        n_cols = Array.length meta.row_names;
        n_rows = Array.length kmer.col_names;
        n_meta = Array.length meta.col_names;
        idx_to_col_names = meta.row_names;
        idx_to_row_names = kmer.col_names;
        idx_to_meta_names = meta.col_names;
        meta = meta.data;
        data = kmer.data
      }
  end: sig
    include module type of KMerDB_Base
    (* Transformations *)
    val transform: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
      Transformation.t -> t -> t
    (* Combinations.
       Generate according to the specified criterion a combination of the spectra having the given labels,
        name the combination as directed, and add it to the database *)
    val add_combined_selected: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                               t -> string -> StringSet.t -> CombinationCriterion.t -> t
    (* Combine spectra into class representatives according to the specified class labels.
       It can fail if the specified metadata field label does not exist,
        or if a spectrum with the same name as the label is already present in the database *)
    val split_spectra: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                       t -> string -> CombinationCriterion.t -> t
    (* Distill k-mers, i.e., sort them by decreasing discriminative power according to the specified class labels.
       It can fail if the specified metadata field label does not exist,
        or if the number of equivalence classes identified by the labels is invalid
        (that is, 1 or the same as the number of k-mers) *)
    val distill_kmers: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> string -> unit
    (* Spectral distance matrix *)
    val to_distances: ?precision:int -> ?normalise:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                      Space.Distance.t -> t -> StringSet.t -> StringSet.t -> string -> unit
    (* Input/output *)
    val to_files: ?precision:int -> ?output_zero_kmers:bool ->
                  ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                  t -> string -> unit
    val of_files: ?threads:int -> ?bytes_per_step:int -> ?elements_per_step:int -> ?verbose:bool ->
                  string -> t
  end
)

