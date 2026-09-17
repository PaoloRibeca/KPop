# Input: ANSI-stripped stderr of one auto run; -v corpus= -v seed=
# Output: rungs to RUNGS file, rounds to ROUNDS file
function flush(r,   i) {
  for (i = 1; i <= nb; i++) print corpus, seed, r, buf[i] > RUNGS
  print corpus, seed, r, pick, level, cls, raxes > ROUNDS
  nb = 0; pick = ""; level = ""
}
/\(KPop__Clustering\.ladder\): [0-9]+ ax(is|es): / {
  line = $0
  sub(/^.*\(KPop__Clustering\.ladder\): /, "", line)
  d = line; sub(/ .*/, "", d)
  s = line; sub(/^[0-9]+ ax(is|es): /, "", s); score = s; sub(/ .*/, "", score)
  k = s; sub(/^.* of /, "", k); kept = k; sub(/ kept.*/, "", kept)
  lst = ""; nle = 0; n = 0
  rest = s; sub(/^[^:]*/, "", rest); sub(/^: /, "", rest); sub(/\.$/, "", rest)
  if (rest != "") {
    m = split(rest, a, /\), /)
    for (j = 1; j <= m; j++) {
      e = a[j]; sub(/\)$/, "", e)
      sh = e; sub(/%.*/, "", sh)
      z = e; sub(/^.*\(z /, "", z)
      lst = lst (lst == "" ? "" : ",") sh ":" z
      n++; if (sh + 0 <= 50) nle++
    }
  }
  buf[++nb] = d OFS score OFS kept OFS n OFS nle OFS lst
  next
}
/\(KPop__Clustering\.ladder\): The ladder picks / {
  pick = $0; sub(/^.*picks /, "", pick); sub(/ .*/, "", pick)
  level = $0; sub(/^.*level of /, "", level); sub(/%.*/, "", level)
  next
}
/\(Dune__exe__KPop_autotuner\): Round [0-9]+: / {
  r = $0; sub(/^.*Round /, "", r); sub(/:.*/, "", r)
  cls = $0; sub(/^.*Round [0-9]+: /, "", cls); sub(/ .*/, "", cls)
  raxes = ""; if ($0 ~ / in [0-9]+ ax/) { raxes = $0; sub(/^.* in /, "", raxes); sub(/ .*/, "", raxes) }
  flush(r)
}
