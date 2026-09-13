#!/usr/bin/env bash
# Regenerates both dumps on the union of the page's grid and the doubling ladder, cutting at the
# radii the calibrated detector found.  Written to a file rather than typed inline: the awk program
# inside radii() needs single quotes, and nesting those inside a quoted command line is what broke
# the last attempt.
set -u
G=1,2,4,5,8,10,16,20,32,40,64,80,120,128
R=/home/ribeca/.claude/jobs/29a54aa8/tmp/collapse/real-data
V=/home/ribeca/.claude/jobs/29a54aa8/tmp/simpson/vp2
L=/home/ribeca/Development/Git/NailIt/test/data/denovo
C=/home/ribeca/.claude/jobs/29a54aa8/tmp/simpson/detector/detect3union.cand

radii() { # <tag> <rung>
  gawk -F'\t' -v t="$1" -v d="$2" '$2 == t && $3 == d && $9 == "kept" && $5 <= 0.5 {
      out = out (out == "" ? "" : ",") $4 } END { print (out == "" ? "-" : out) }' "$C"
}
spec() { local tag=$1; shift; local s=""; local d; for d in "$@"; do s="$s${s:+;}$d=$(radii "$tag" "$d")"; done; echo "$s"; }

cd /home/ribeca/.claude/jobs/29a54aa8/tmp/simpson/valleys
rm -f hist.jsonl vdump.progress vdump.err
vd() { local tag=$1 dir=$2 ds=$3 lv=$4; shift 4
  nice -n 19 ./vdump "$dir/$tag.all.KPopTwisted.txt" "$dir/$tag.all.KPopInertia.txt" "$tag" "$ds" "$lv" \
    "$(spec "$tag" ${ds//,/ })" "$@" >> hist.jsonl 2>> vdump.err
  echo "$(date +%T) $tag exit=$?" >> vdump.progress
}
vd vp1.T_rand  $R $G,243 genotype,variant $L/vp1.cdc $L/vp1.refined
vd vp1.T_strat $R $G,243 genotype,variant $L/vp1.cdc $L/vp1.refined
vd vp1.T_core  $R $G,243 genotype,variant $L/vp1.cdc $L/vp1.refined
vd rdrp.T_rand  $R $G,199 P-type $L/rdrp.cdc
vd rdrp.T_strat $R $G,199 P-type $L/rdrp.cdc
vd vp2.T_rand   $V $G,219 genotype $L/vp2.cdc
echo "vdump lines=$(wc -l < hist.jsonl)"; cat vdump.progress; head -c 300 vdump.err

cd ../phylo
rm -f hist-tree.jsonl vtree.progress vtree.err
vt() { local tag=$1 dir=$2 ds=$3 lv=$4; shift 4
  nice -n 19 ./vtree "$dir/$tag.all.KPopTwisted.txt" "$dir/$tag.all.KPopInertia.txt" "$tag" "$ds" "$lv" \
    "$(spec "$tag" ${ds//,/ })" "$@" >> hist-tree.jsonl 2>> vtree.err
  echo "$(date +%T) $tag exit=$?" >> vtree.progress
}
vt vp1.T_strat $R $G,243 genotype,variant $L/vp1.cdc $L/vp1.refined
vt vp1.T_rand  $R $G,243 genotype,variant $L/vp1.cdc $L/vp1.refined
vt vp1.T_core  $R $G,243 genotype,variant $L/vp1.cdc $L/vp1.refined
vt rdrp.T_rand  $R $G,199 P-type $L/rdrp.cdc
vt rdrp.T_strat $R $G,199 P-type $L/rdrp.cdc
vt vp2.T_rand   $V $G,219 genotype $L/vp2.cdc
echo "vtree lines=$(wc -l < hist-tree.jsonl)"; cat vtree.progress; head -c 300 vtree.err
