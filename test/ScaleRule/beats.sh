#!/usr/bin/env bash
# beats.sh <tsv...>: every partition that scores above CDC under the same silhouette, guard and s.
# TSV columns: tag d cat name k f_p noise4 noise8 noise16 silhouette guard s score
set -u
cd /home/ribeca/.claude/jobs/29a54aa8/tmp/guard2
cat "$@" > all.tsv
echo "#### attacks above CDC (degen, grad, peel, noise); crit/proj/pert are references, listed after"
awk -F'\t' '
NR == FNR { if ($3 == "ref") c[$1, $10, $11, $12] = $13; next }
$3 == "degen" || $3 == "grad" || $3 == "peel" || $3 == "noise" {
  v = c[$1, $10, $11, $12]
  if ($13 > v + 1e-9)
    printf "%-10s %-10s %-3s s=%-3s %-5s %-58s %7.4f  CDC %7.4f  +%.4f\n", $1, $10, $11, $12, $3, substr($4, 1, 58), $13, v, $13 - v
}' all.tsv all.tsv
echo
echo "#### references above CDC (crit, proj, pert)"
awk -F'\t' '
NR == FNR { if ($3 == "ref") c[$1, $10, $11, $12] = $13; next }
$3 == "crit" || $3 == "proj" || $3 == "pert" {
  v = c[$1, $10, $11, $12]
  if ($13 > v + 1e-9)
    printf "%-10s %-10s %-3s s=%-3s %-5s %-58s %7.4f  CDC %7.4f  +%.4f\n", $1, $10, $11, $12, $3, substr($4, 1, 58), $13, v, $13 - v
}' all.tsv all.tsv
echo
echo "#### how many attacks beat CDC, and the worst of them, per embedding / silhouette / guard / s"
awk -F'\t' '
NR == FNR { if ($3 == "ref") c[$1, $10, $11, $12] = $13; next }
$3 == "degen" || $3 == "grad" || $3 == "peel" || $3 == "noise" {
  k = $1 SUBSEP $10 SUBSEP $11 SUBSEP $12
  v = c[$1, $10, $11, $12]
  if ($13 > v + 1e-9) { n[k]++; if ($13 - v > m[k]) { m[k] = $13 - v; w[k] = $4 } }
  seen[k] = 1
}
END {
  for (k in seen) {
    split(k, f, SUBSEP)
    printf "%-10s %-10s %-3s s=%-3s CDC %7.4f  beaten by %2d, worst +%.4f  %s\n",
      f[1], f[2], f[3], f[4], c[f[1], f[2], f[3], f[4]], n[k] + 0, m[k] + 0, (n[k] ? w[k] : "-")
  }
}' all.tsv all.tsv | sort -k1,1 -k2,2 -k4,4 -k3,3
