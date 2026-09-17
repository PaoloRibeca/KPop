let strip s =
  let s = String.trim s in
  let l = String.length s in
  if l >= 2 && s.[0] = '"' && s.[l - 1] = '"' then String.sub s 1 (l - 2) else s

let read_inertia fn =
  let ic = open_in fn in
  ignore (input_line ic);
  let l = input_line ic in
  close_in ic;
  match String.split_on_char '\t' (String.trim l) with
  | _ :: vs -> Array.of_list (List.map float_of_string vs)
  | [] -> failwith ("empty inertia line in " ^ fn)

let read_twisted fn =
  let ic = open_in fn in
  ignore (input_line ic);
  let names = ref [] and rows = ref [] in
  (try
     while true do
       let l = input_line ic in
       if String.trim l <> "" then
         match String.split_on_char '\t' (String.trim l) with
         | nm :: vs ->
           names := strip nm :: !names;
           rows := Array.of_list (List.map float_of_string vs) :: !rows
         | [] -> ()
     done
   with End_of_file -> ());
  close_in ic;
  Array.of_list (List.rev !names), Array.of_list (List.rev !rows)

let read_labels fn =
  let h = Hashtbl.create 4096 in
  let ic = open_in fn in
  (try
     while true do
       let l = String.trim (input_line ic) in
       if l <> "" then
         match String.split_on_char '\t' l with
         | [ nm; lab ] -> Hashtbl.replace h (strip nm) (strip lab)
         | _ -> failwith ("bad label line: " ^ l)
     done
   with End_of_file -> ());
  close_in ic;
  h
