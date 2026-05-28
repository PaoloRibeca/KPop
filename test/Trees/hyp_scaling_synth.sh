#!/usr/bin/env bash
# Scaling study: simulate phylogenetic sequences, run KPop pipeline,
# sweep sparse-NJ modes at increasing n.  Asks whether the K_QUERY/K
# factor needed for hyperbolic mode to match the dense reference
# stays bounded as n grows.
#
# Memory: simulator + KPop pipeline is sub-GB for n up to ~500 at
# seqlen=800.  Tested on a machine with 8 GB RAM.
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

# Sample n values bounded by 8 GB RAM
NS=${NS:-"50 100 200 400"}
FACTORS="1 2 3 5"

printf "%-5s %-14s %-7s %-12s %-12s %-8s\n" "n" "mode" "factor" "J(vs dense)" "J(vs truth)" "time_s"
for n in $NS; do
  pfx="$TMP/sim_n${n}"
  # Simulate
  if [[ ! -s "${pfx}.fasta" ]]; then
    python3 "$SIM" "$pfx" "$n" "$SEQLEN" "$BRANCH_SCALE" 42 "$MU" 2>/dev/null
  fi
  # KPop pipeline (skip KPopCountDB: KPopTwist accepts raw spectra
  # and treats each sample as its own class)
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
  printf "%-5s %-14s %-7s %-12s %-12s %-8s\n" "$n" "dense" "--" "1.000" "$jt" "$dt"
  # Modes
  for mode in subquadratic hyperbolic; do
    for f in $FACTORS; do
      out="${pfx}_${mode}_f${f}"
      t0=$(date +%s.%N)
      extra=""
      if [[ "$mode" == "hyperbolic" ]]; then
        extra="--phylo-snj-hyp-scale 0.7"
      else
        extra="--phylo-snj-index-type flat"
      fi
      "$TWISTDB" -i T "${pfx}_tw" -i t "${pfx}_tw" -m "powers(1,1,1)" \
        --phylo-method sparse-nj --phylo-snj-mode $mode \
        --phylo-snj-num-neighbors $K --phylo-snj-k-query-factor $f \
        $extra --phylo-snj-rowsum knn --phylo-snj-symmetry one \
        -P "$out" 2>/dev/null
      t1=$(date +%s.%N)
      dt=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
      jd=$(python3 "$RF" "${pfx}_dense.nwk" "${out}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
      jt=$(python3 "$RF" "${pfx}.tree.nwk" "${out}.nwk" 2>/dev/null | grep -oP '\d\.\d{3}' | tail -1)
      printf "%-5s %-14s %-7s %-12s %-12s %-8s\n" "$n" "$mode" "$f" "$jd" "$jt" "$dt"
    done
  done
  echo
done
