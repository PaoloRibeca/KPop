#!/usr/bin/env bash
# The valley cuts against KPopPhylo's trees, one case at a time on one core, then the JSON the
# page reads.  Each cut is at a valley of the pair-distance histogram (linear, +-5 bins), the
# same cuts the page lists as "At a pair-distance valley".  Those radii come from the
# calibrated detector's own run (detect3page.cand: tag, rung, radius, share, ..., status),
# and NOT from a scrape of hist.jsonl as they once did -- that file is written by vdump,
# so reading it here made this script depend on the order two programs happened to run in.
set -euo pipefail
cd "$(dirname "$0")"
R=/home/ribeca/.claude/jobs/29a54aa8/tmp/collapse/real-data
V=/home/ribeca/.claude/jobs/29a54aa8/tmp/simpson/vp2
L=/home/ribeca/Development/Git/NailIt/test/data/denovo
C=../detector/detect3union.cand
radii() { # <tag> <rung>: the kept valleys at or below half the pairs, in distance units
  gawk -F'\t' -v t="$1" -v d="$2" '$2 == t && $3 == d && $9 == "kept" && $5 <= 0.5 {
      out = out (out == "" ? "" : ",") $4 } END { print (out == "" ? "-" : out) }' "$C"
}
case_() { # <tag> <dir of the text exports> <d> <level names> <labels...>
  local tag=$1 dir=$2 d=$3 lv=$4
  shift 4
  if [ ! -s "$tag.d$d.nj.nwk" ]; then echo "missing $tag.d$d.nj.nwk" >&2; return 0; fi
  nice -n 19 ./splitcmp "$dir/$tag.all.KPopTwisted.txt" "$dir/$tag.all.KPopInertia.txt" "$d" "$tag" \
    "$(radii "$tag" "$d")" "$tag.d$d.nj.nwk" "$lv" "$@"
}
{
  case_ vp1.T_strat $R 20 genotype,variant $L/vp1.cdc $L/vp1.refined
  case_ vp1.T_strat $R 243 genotype,variant $L/vp1.cdc $L/vp1.refined
  case_ vp1.T_rand $R 10 genotype,variant $L/vp1.cdc $L/vp1.refined
  case_ vp1.T_core $R 20 genotype,variant $L/vp1.cdc $L/vp1.refined
  case_ vp2.T_rand $V 40 genotype $L/vp2.cdc
  case_ rdrp.T_rand $R 5 P-type $L/rdrp.cdc
  case_ rdrp.T_strat $R 10 P-type $L/rdrp.cdc
} > compare.out 2> compare.err
grep '^JSON ' compare.out | sed 's/^JSON //' \
  | gawk 'BEGIN { printf "[" } NR > 1 { printf "," } { printf "%s", $0 } END { print "]" }' > ../valleys/viz/phylo.json
