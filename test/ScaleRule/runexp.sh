#!/usr/bin/env bash
# runexp.sh <arm> [...]: the threshold-scale experiment, one guard3 run at a time.
#
# THE ARMS EXIST BECAUSE ONE OF THEM CANNOT ANSWER THE QUESTION ALONE.  Under proportional thinning
# n, n/k and every percentile of the class-size distribution fall by the same factor, so "the
# threshold is a fraction of n" and "the threshold tracks class size" predict the identical
# measurement.  Arm B moves n with the class sizes held fixed; arm C moves the class sizes with n,
# the points and the geometry held bit-identical.  Only across the three is the pair separable.
#
#   anchor          frac 0.95, one draw          -- the reference point, drawn rather than whole so
#                                                   that it carries a replicate band like the rest
#   A <frac> <reps> proportional thinning        -- known-confounded, run to measure the confound
#   B <frac> <reps> whole classes dropped        -- n falls, class sizes fixed
#   C <frac>        reference coarsened/refined  -- n fixed, class sizes move
#
# Runs are serial: each is single-threaded by construction and the machine is shared.
#
# WORK holds the guard3 binary and receives every output; TWISTED holds the text exports
# {vp1,rdrp,vp2}.T_rand.all.KPopTwisted.txt with their .KPopInertia.txt; LABELS holds the curated
# labels.  Completeness is judged by the sentinel row, so an interrupted arm can be rerun as it is.
set -u

G="${WORK:?set WORK to the directory holding the guard3 binary}"
R="${TWISTED:?set TWISTED to the directory holding the three .all.KPopTwisted.txt exports}"
L="${LABELS:-$HOME/Development/Git/NailIt/test/data/denovo}"

run () { # <twisted prefix> <labels> <axes> <base tag> <frac> <seed> <mode> <refmode>
  local tag="$4.f$5.r$6.$7.$8"
  # A killed run leaves a well-formed prefix, so completeness is judged on the sentinel row the
  # rig writes last, never on the file merely existing
  if [[ -s "$G/exp2.$tag.tsv" ]] && gawk -F'\t' '$3 == "end" { f = 1 } END { exit !f }' "$G/exp2.$tag.tsv"; then
    echo "=== $tag already complete, skipped"
    return
  fi
  echo "=== $tag started $(date +%T)"
  nice -n 19 "$G/guard3" "$1.all.KPopTwisted.txt" "$1.all.KPopInertia.txt" "$2" "$3" "$tag" \
    "$G/exp2.$tag.tsv" "$5" "$6" "$7" "$8" > "$G/exp2.$tag.out" 2> "$G/exp2.$tag.err"
  echo "=== $tag exit $? $(date +%T) n=$(gawk -F'\t' 'NR == 1 { print $16 }' "$G/exp2.$tag.tsv")"
}

each () { # <frac> <seed> <mode> <refmode>
  run $R/vp1.T_rand  $L/vp1.cdc  20 vp1.d20  "$1" "$2" "$3" "$4"
  run $R/rdrp.T_rand $L/rdrp.cdc  5 rdrp.d5  "$1" "$2" "$3" "$4"
  run $R/vp2.T_rand  $L/vp2.cdc  20 vp2.d20  "$1" "$2" "$3" "$4"
}

case "${1:-anchor}" in
  anchor) each 0.95 1 uniform none ;;
  A) for r in $(seq 1 "$3"); do each "$2" "$r" uniform none; done ;;
  B) for r in $(seq 1 "$3"); do each "$2" "$r" classes none; done ;;
  C) each "$2" 1 uniform merge5; each "$2" 1 uniform split ;;
  *) echo "unknown arm ${1:-}" >&2; exit 2 ;;
esac
echo "=== arm ${1:-anchor} done $(date +%T)"
