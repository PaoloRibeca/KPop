# reads de-coloured, \r-split stderr; prints one line per round:
# tag round pick | counts per rung | maxlow(<=4) maxhigh(>=16) | picked rung's counted valleys (share/z) | tie-break flag | resample line
function flushr(   i, d, maxl, maxh, maxall, cnt, first, out, tie) {
  if (nr == 0) return
  maxl = -1; maxh = -1; maxall = -1; out = ""
  for (i = 1; i <= nr; i++) {
    d = rd[i]; out = out d ":" rc[i] " "
    if (d + 0 <= 4 && rc[i] > maxl) maxl = rc[i]
    if (d + 0 >= 16 && rc[i] > maxh) maxh = rc[i]
    if (rc[i] > maxall) maxall = rc[i]
  }
  cnt = 0; for (i = 1; i <= nr; i++) if (rc[i] == maxall) cnt++
  tie = (cnt > 1) ? "TIE(" cnt " rungs at max " maxall ")" : "unique max " maxall
  printf "%s\tr%d\tpick=%s\t%s\tmaxLow=%d maxHigh=%d\t%s\tpicked:%s\t%s\n", tag, round, pick, out, maxl, maxh, tie, pv[pick], res
  # high-rung single-valley stats
  for (i = 1; i <= nr; i++) if (rd[i] + 0 >= 16) {
    printf "HIGH\t%s\tr%d\t%s\t%d\t%s\n", tag, round, rd[i], rc[i], rv[i] > "/dev/stderr"
  }
  nr = 0; delete rd; delete rc; delete rv; delete pv; res = ""
}
/ladder\): [0-9]+ ax(is|es): / {
  line = $0; sub(/.*ladder\): /, "", line)
  split(line, a, " "); d = a[1]
  match(line, /: ([0-9]+) valleys? with at most half/, m); c = m[1]
  vals = line; sub(/.* kept: /, "", vals); gsub(/\.$/, "", vals)
  # keep only counted valleys (the first c)
  n = split(vals, vv, /, /); keep = ""
  for (j = 1; j <= c && j <= n; j++) keep = keep (j > 1 ? "," : "") vv[j]
  nr++; rd[nr] = d; rc[nr] = c; rv[nr] = keep; pv[d] = keep
  next
}
/The ladder picks/ { match($0, /picks ([0-9]+) ax/, m); pick = m[1]; next }
/The resamples pick/ { res = $0; sub(/.*The resamples /, "resamples ", res); next }
/autotuner\): Round [0-9]+:/ { match($0, /Round ([0-9]+):/, m); round = m[1]; flushr(); next }
