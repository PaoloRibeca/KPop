#!/usr/bin/env bash
# mstar.sh <tsv...>: how large a clump can still fool the objective, measured with NO guard.
#
# WHY THIS AND NOT A THRESHOLD SWEEP.  Sweeping the guard's threshold mostly measures the corpus:
# an attack whose clusters are smaller than the threshold is charged -1 per member and scores about
# 0, so a guarded estimate lands one defence rung above the largest group that already wins with no
# guard.  The quantity underneath is simpler: take a single group of m points against everything
# else, score it with no guard at all, and ask for which m it beats the reference.  The largest
# such m is m*.  It is NOT a safe threshold: under G3b a group of exactly s members is eligible, and
# groups larger than m* can beat the guarded reference, whose own score moves with s.
#
# Under G0 there is no charge and no eligibility.  G3b = G2 - noise/n, so a guarded reference is
# charged for its own small classes -- more heavily as the data is thinned -- which moves every
# guarded comparison for reasons that have nothing to do with the attacks.
#
# The probes are the families that put one group against the rest, whose name gives the group's
# size: "best seed+N NN group" and "farthest point + N NN as one group vs rest" (N+1 points), and
# "all but farthest N, as one group" (the N points farthest from the global centroid).  Every one
# is a two-cluster partition of periphery against bulk, coarser than any labelling, and every one
# is seeded at the periphery -- so m* is a lower bound, interval-censored by the size ladder.
set -u

gawk -F'\t' '
$3 == "end" { next }
$11 == "G0" {
  t = $1; kd = $10
  N[t] = $16; FR[t] = $14; SD[t] = $15
  if ($3 == "ref") { ref[t, kd] = $13 + 0; KR[t, kd] = $5 + 0; cells[t SUBSEP kd] = 1; next }
  m = 0
  if (match($4, /^best seed\+[0-9]+ NN group/)) { split($4, a, " "); sub(/^seed\+/, "", a[2]); m = a[2] + 1 }
  else if (match($4, /^farthest point \+ [0-9]+ NN/)) { split($4, a, " "); m = a[4] + 1 }
  else if (match($4, /^all but farthest [0-9]+, as one group$/)) { split($4, a, " "); m = a[4] + 0 }
  else next
  if (!((t, kd, m) in best) || $13 + 0 > best[t, kd, m]) best[t, kd, m] = $13 + 0
  if (!((t, kd, m) in seen)) { seen[t, kd, m] = 1; ms[t, kd] = ms[t, kd] " " m }
}
END {
  # k is printed because the granularity arm moves exactly that, and an m* with no k beside it says
  # nothing about whether m* follows the labelling
  printf "%-30s %-10s %6s %5s %3s %4s %5s  %s\n",
    "run", "silhouette", "n", "frac", "sd", "k", "m*", "advantage over the labels by group size"
  for (c in cells) {
    split(c, f, SUBSEP); t = f[1]; kd = f[2]
    k = split(ms[t, kd], v, " ")
    for (i = 1; i <= k; i++) for (j = i + 1; j <= k; j++) if (v[j] + 0 < v[i] + 0) { x = v[i]; v[i] = v[j]; v[j] = x }
    line = ""; star = 0
    for (i = 1; i <= k; i++) {
      m = v[i] + 0
      d = best[t, kd, m] - ref[t, kd]
      if (d > 0) star = m
      line = line sprintf(" %d:%+.4f%s", m, d, (d > 0 ? "*" : ""))
    }
    # A group that still wins at the largest size probed says only that m* is at least that large
    top = v[k] + 0
    printf "%-30s %-10s %6d %5.3f %3d %4d %5s  %s\n",
      substr(t, 1, 30), kd, N[t], FR[t], SD[t], KR[t, kd], (star ? (star == top ? ">=" star : star "") : "none"), line
  }
}
' "$@"
