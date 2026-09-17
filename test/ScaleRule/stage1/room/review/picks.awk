# Densest valley-partition rung per tag, per density variant and leader order, with no floor and with floor 8
BEGIN { FS = "\t"; split("dens_id dens_id_nn dens_id_ge5 dens_id_n32 dens_id_n32_nn", V, " ") }
$1 ~ /^#/ || $1 == "tag" { if ($1 == "tag") for (i = 1; i <= NF; i++) col[$i] = i; next }
$3 == "labels" { lwt[$1, $2] = $12; ll[$1, $2] = $14; tags[$1] = 1; next }
$3 == "valley" {
  for (v = 1; v <= 5; v++) {
    x = $(col[V[v]]); k = $1 SUBSEP $4 SUBSEP V[v]
    if (!(k in b0) || x > b0[k]) { b0[k] = x; p0[k] = $2 }
    if ($2 >= 8 && (!(k in b8) || x > b8[k])) { b8[k] = x; p8[k] = $2 }
  }
  keys[$1 SUBSEP $4] = 1
}
END {
  for (t in tags) for (o in keys) { split(o, a, SUBSEP); if (a[1] != t) continue
    line = t "\t" a[2]
    for (v = 1; v <= 5; v++) { k = t SUBSEP a[2] SUBSEP V[v]; line = line "\t" V[v] ": " p0[k] "/" p8[k] }
    print line }
}
