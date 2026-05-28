#!/usr/bin/env bash
# Focused scaling sweep: dense + hyperbolic at multiple K_QUERY, with
# automatically picked scale (smaller-as-n-grows).  Subquadratic at
# matching factors for comparison.  Designed for n in the 100--800
# range with the 8 GB memory budget.
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
NS=${NS:-"50 100 200 400"}

# Scale tuning rule based on observed optima at n=50..200:
#   n=50 best ~0.7;  n=100 best ~0.5;  n=200 best ~0.3
# Roughly scale ~ 1 / sqrt(n / 50).  Use that as the auto-pick.
auto_scale() {
  local n=$1
  awk -v n="$n" 'BEGIN{printf "%.2f", 0.7 / sqrt(n/50.0)}'
}

printf "%-5s %-14s %-7s %-7s %-12s %-12s %-8s\n" "n" "mode" "factor" "scale" "J(vs dense)" "J(vs truth)" "time_s"
for n in $NS; do
  pfx="$TMP/sim_n${n}"
  if [[ ! -s "${pfx}.fasta" ]]; then
    python3 "$SIM" "$pfx" "$n" "$SEQLEN" "$BRANCH_SCALE" 42 "$MU" 2>/dev/null
  fi
  if [[ ! -s "${pfx}_tw.KPopTwisted" ]]; then
    "$COUNT" -k $KMER_K -f "" "${pfx}.fasta" -o "$pfx" 2>/dev/null
    "$TWIST" -i "$pfx" -o "${pfx}_tw" 2>/dev/null
  fi
  sc=$(auto_scale "$n")
  # Dense reference
  t0=$(date +%s.%N)
  "$TWISTDB" -i T "${pfx}_tw" -i t "${pfx}_tw" -m "powers(1,1,1)" \
    --phylo-method sparse-nj --phylo-snj-mode dense \
    --phylo-snj-num-neighbors $K --phylo-snj-rowsum knn --phylo-snj-symmetry one \
    --phylo-snj-index-type flat -P "${pfx}_dense" 2>/dev/null
  t1=$(date +%s.%N)
  dt=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
  jt=$(python3 "$RF" "${pfx}.tree.nwk" "${pfx}_dense.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
  printf "%-5s %-14s %-7s %-7s %-12s %-12s %-8s\n" "$n" "dense" "--" "--" "1.000" "$jt" "$dt"
  # One hyperbolic K_QUERY factor (the smallest that achieves J=1 at this n).
  # Try factors 1..5 and report each.
  for f in 1 2 3 5; do
    t0=$(date +%s.%N)
    "$TWISTDB" -i T "${pfx}_tw" -i t "${pfx}_tw" -m "powers(1,1,1)" \
      --phylo-method sparse-nj --phylo-snj-mode hyperbolic \
      --phylo-snj-num-neighbors $K --phylo-snj-k-query-factor $f \
      --phylo-snj-hyp-scale $sc --phylo-snj-rowsum knn --phylo-snj-symmetry one \
      -P "${pfx}_hyp_f${f}" 2>/dev/null
    t1=$(date +%s.%N)
    dt=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
    jd=$(python3 "$RF" "${pfx}_dense.nwk" "${pfx}_hyp_f${f}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
    jt=$(python3 "$RF" "${pfx}.tree.nwk" "${pfx}_hyp_f${f}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
    printf "%-5s %-14s %-7s %-7s %-12s %-12s %-8s\n" "$n" "hyperbolic" "$f" "$sc" "$jd" "$jt" "$dt"
  done
  for f in 1 2 3 5; do
    t0=$(date +%s.%N)
    "$TWISTDB" -i T "${pfx}_tw" -i t "${pfx}_tw" -m "powers(1,1,1)" \
      --phylo-method sparse-nj --phylo-snj-mode subquadratic \
      --phylo-snj-num-neighbors $K --phylo-snj-k-query-factor $f \
      --phylo-snj-index-type flat \
      --phylo-snj-rowsum knn --phylo-snj-symmetry one \
      -P "${pfx}_sq_f${f}" 2>/dev/null
    t1=$(date +%s.%N)
    dt=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
    jd=$(python3 "$RF" "${pfx}_dense.nwk" "${pfx}_sq_f${f}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
    jt=$(python3 "$RF" "${pfx}.tree.nwk" "${pfx}_sq_f${f}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
    printf "%-5s %-14s %-7s %-7s %-12s %-12s %-8s\n" "$n" "subquadratic" "$f" "--" "$jd" "$jt" "$dt"
  done
  echo
done
