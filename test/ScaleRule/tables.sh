#!/usr/bin/env bash
# tables.sh <s> <tsv...>: compact table at guard size <s>, the selected partitions only,
# with G0 G1 G2 G3a G3b for both silhouettes.
set -u
S="$1"
shift
awk -F'\t' -v S="$S" '
$11 == "G0" || $12 == S {
  key = $1 SUBSEP $4 SUBSEP $10
  if ($11 == "G0") {
    g0[key] = $13
    tn = $1 SUBSEP $4
    if (!(tn in seen)) { seen[tn] = 1; ord[++nn] = tn; K[tn] = $5; F[tn] = $6 }
  } else
    v[key SUBSEP $11] = $13
}
END {
  pat = "^CDC$|^CDC largest|^CDC two nearest|^CDC genogroups|^CDC polished|^CDC 1% moved|^CDC, its|^G3a projection|^all but farthest (1, as singletons|4, as one group|8, as one group|16, as one group)|^best seed\\+(3|7|15) NN group.*by si|^peel [2345] farthest from all|^peel 8 farthest from all|^peel farthest member \\+ (3|7|15) within-class NN, all|^greedy NN pairs \\(atomised|^greedy pairs of pairs|^greedy atomised, sizes >= (8|16)|^greedy NN pairs \\+ tightest (4|8|16)-group|^two tightest disjoint (4|8|16)-groups|^tightest (4|8|16)-group \\+ singletons|^closest pair|^all singletons|^one cluster"
  printf "#### guard size s=%s;  per silhouette: G0  G1  G2  G3a  G3b\n", S
  printf "%-9s %-56s %5s %6s | %-39s | %-39s\n", "tag", "partition", "k", "f_p", "classical", "simplified"
  for (i = 1; i <= nn; i++) {
    split(ord[i], f, SUBSEP); t = f[1]; nm = f[2]
    if (nm ~ pat) {
      printf "%-9s %-56s %5d %6.4f", t, substr(nm, 1, 56), K[ord[i]], F[ord[i]]
      for (j = 1; j <= 2; j++) {
        s = (j == 1 ? "classical" : "simplified"); k2 = t SUBSEP nm SUBSEP s
        printf " | %7.4f %7.4f %7.4f %7.4f %7.4f", g0[k2], v[k2 SUBSEP "G1"], v[k2 SUBSEP "G2"], v[k2 SUBSEP "G3a"], v[k2 SUBSEP "G3b"]
      }
      printf "\n"
    }
  }
}' "$@"
