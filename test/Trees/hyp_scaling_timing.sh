#!/usr/bin/env bash
# Timing-focused scaling measurement of the current (rebuild-per-merge)
# subquadratic mode vs dense reference.  No quality sweep -- just the
# wall-clock numbers and the J(vs dense) at the K_QUERY/K factor where
# subquadratic empirically lands ~ 1.000 on this fixture.
#
# Reuses synthetic Yule-tree + JC69 sequences through KPop's pipeline.
set -eu
cd "$(dirname "$0")/../.."
KPOP=$(pwd)
TMP=${HYP_TMP:-/tmp/hyp_scaling}
mkdir -p "$TMP"
SIM="$KPOP/test/Trees/simulate_phylo.py"
TWISTDB="$KPOP/build/KPopTwistDB"
COUNT="$KPOP/build/KPopCount"
TWIST="$KPOP/build/KPopTwist"
RF="$KPOP/test/Trees/tree_rf.py"

K=10
SEQLEN=800
BRANCH_SCALE=0.3
MU=0.05
KMER_K=10
NS=${NS:-"50 100 200 400 800"}
FACTOR=${FACTOR:-3}

printf "%-5s %-14s %-7s %-7s %-12s %-12s %-10s\n" "n" "mode" "factor" "index" "J(vs dense)" "J(vs truth)" "time_s"
for n in $NS; do
  pfx="$TMP/sim_n${n}"
  if [[ ! -s "${pfx}.fasta" ]]; then
    python3 "$SIM" "$pfx" "$n" "$SEQLEN" "$BRANCH_SCALE" 42 "$MU" 2>/dev/null
  fi
  if [[ ! -s "${pfx}_tw.KPopTwisted" ]]; then
    "$COUNT" -k $KMER_K -f "" "${pfx}.fasta" -o "$pfx" 2>/dev/null
    "$TWIST" -i "$pfx" -o "${pfx}_tw" 2>/dev/null
  fi
  # Dense reference
  t0=$(date +%s.%N)
  "$TWISTDB" -i T "${pfx}_tw" -i t "${pfx}_tw" -m "powers(1,1,1)" \
    --phylo-method sparse-nj --phylo-snj-mode dense \
    --phylo-snj-num-neighbors $K --phylo-snj-rowsum knn --phylo-snj-symmetry one \
    --phylo-snj-index-type flat -P "${pfx}_dense" 2>/dev/null
  t1=$(date +%s.%N)
  dt=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
  jt=$(python3 "$RF" "${pfx}.tree.nwk" "${pfx}_dense.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
  printf "%-5s %-14s %-7s %-7s %-12s %-12s %-10s\n" "$n" "dense" "--" "flat" "1.000" "$jt" "$dt"
  # Subquadratic at the chosen factor, for both flat and hnsw(32)
  for idx in flat "hnsw(32)"; do
    out="${pfx}_sq_${idx//[()]/_}"
    t0=$(date +%s.%N)
    "$TWISTDB" -i T "${pfx}_tw" -i t "${pfx}_tw" -m "powers(1,1,1)" \
      --phylo-method sparse-nj --phylo-snj-mode subquadratic \
      --phylo-snj-num-neighbors $K --phylo-snj-k-query-factor $FACTOR \
      --phylo-snj-index-type "$idx" \
      --phylo-snj-rowsum knn --phylo-snj-symmetry one \
      -P "$out" 2>/dev/null
    t1=$(date +%s.%N)
    dt=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
    jd=$(python3 "$RF" "${pfx}_dense.nwk" "${out}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
    jt=$(python3 "$RF" "${pfx}.tree.nwk" "${out}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
    printf "%-5s %-14s %-7s %-7s %-12s %-12s %-10s\n" "$n" "subquadratic" "$FACTOR" "$idx" "$jd" "$jt" "$dt"
  done
  echo
done
