#!/usr/bin/env bash
# VP2, the third system: the same measurement as the real-data strand's T_rand, reusing its
# compiled analysis, one core at nice 19.  220 spectra = the autotuner's automatic sample for 3027
set -euo pipefail
cd "$(dirname "$0")"
B=/home/ribeca/Development/Git/KPop.Claude/.build
DB=/home/ribeca/.claude/jobs/29a54aa8/tmp/vp2.spectra
L=/home/ribeca/Development/Git/NailIt/test/data/denovo
S=/home/ribeca/.claude/jobs/29a54aa8/tmp/collapse/real-data/src
cut -f1 $L/vp2.denovo.labels | shuf --random-source=<(yes 17) -n 220 > vp2.T_rand.sel
SEL="$(paste -sd, vp2.T_rand.sel)"
nice -n 19 $B/KPopCountDB -T 1 -i $DB -L "$SEL" -N -D -o vp2.T_rand.subset 2> countdb.err
nice -n 19 $B/KPopTwist -T 1 --kmers-threshold off -v -i vp2.T_rand.subset -o vp2.T_rand 2> twist.err
nice -n 19 $B/KPopTwistDB -T 1 -i T vp2.T_rand -t $DB -O t vp2.T_rand.all 2> twistdb.err
rm -f vp2.T_rand.subset.KPopSpectra
nice -n 19 $S/analyze vp2.T_rand.all.KPopTwisted.txt vp2.T_rand.all.KPopInertia.txt $L/vp2.cdc vp2.T_rand > vp2.T_rand.analyze.log 2>&1
echo done
