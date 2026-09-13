#!/usr/bin/env bash
# KPopPhylo's trees on the embeddings the valley cuts were measured in, one job at a time on one
# core.  Each projection becomes a binary twisted register of every spectrum; each tree keeps the
# first d axes through the metric's threshold, which retains axis k exactly when the inertia
# share of the axes before it is below the threshold -- so a threshold halfway between the shares
# at d-1 and d keeps d of them.  VP2 is restricted to the 3027-sequence corpus of the earlier runs
set -euo pipefail
cd "$(dirname "$0")"
B=/home/ribeca/Development/Git/KPop.Claude/.build
R=/home/ribeca/.claude/jobs/29a54aa8/tmp/collapse/real-data
V=/home/ribeca/.claude/jobs/29a54aa8/tmp/simpson/vp2
L=/home/ribeca/Development/Git/NailIt/test/data/denovo
REG=/home/ribeca/.claude/jobs/29a54aa8/tmp/regions
stamp() { printf '%s %s\n' "$(date +%H:%M:%S)" "$*" >> progress.log; }

if [ ! -f vp2.corpus.KPopSpectra ]; then
  nice -n 19 $B/KPopCountDB -T 1 -i /home/ribeca/.claude/jobs/29a54aa8/tmp/vp2.spectra \
    -L "$(cut -f1 $L/vp2.denovo.labels | paste -sd,)" -N -D -o vp2.corpus 2> vp2.corpus.err
  stamp "vp2 corpus database written"
fi

project() { # <tag> <twister dir> <spectra prefix>
  if [ ! -f "$1.all.KPopTwisted" ]; then
    nice -n 19 $B/KPopTwistDB -T 1 -i T "$2/$1" -t "$3" -o t "$1.all" 2> "$1.project.err"
    stamp "projected $1"
  fi
}
project vp1.T_strat $R $REG/vp1
project vp1.T_rand $R $REG/vp1
project vp1.T_core $R $REG/vp1
project rdrp.T_rand $R $REG/rdrp
project rdrp.T_strat $R $REG/rdrp
project vp2.T_rand $V vp2.corpus

threshold() { # <inertia.txt> <d>
  gawk -v d="$2" 'NR == 2 {
      n = NF - 1; tot = 0; for (i = 2; i <= NF; i++) tot += $i
      if (d >= n) { print 1; exit }
      a = 0; for (i = 2; i <= d; i++) a += $i
      printf "%.15g\n", (a + (a + $(d + 1))) / 2 / tot; exit }' "$1"
}

tree() { # <tag> <d> <inertia.txt>
  local out="$1.d$2.nj" t
  t="$(threshold "$3" "$2")"
  if [ ! -s "$out.nwk" ]; then
    stamp "tree $out starting (threshold $t)"
    nice -n 19 $B/KPopPhylo -T 1 -v -i t "$1.all" -m "powers(1,$t,1)" --method sparse-nj \
      --snj-mode periodic-rebuild -o "$out" > "$out.log" 2>&1
    stamp "tree $out done"
  fi
}
tree vp1.T_strat 20 $R/vp1.T_strat.all.KPopInertia.txt
tree rdrp.T_rand 5 $R/rdrp.T_rand.all.KPopInertia.txt
tree rdrp.T_strat 10 $R/rdrp.T_strat.all.KPopInertia.txt
tree vp1.T_rand 10 $R/vp1.T_rand.all.KPopInertia.txt
tree vp2.T_rand 40 $V/vp2.T_rand.all.KPopInertia.txt
tree vp1.T_core 20 $R/vp1.T_core.all.KPopInertia.txt
tree vp1.T_strat 243 $R/vp1.T_strat.all.KPopInertia.txt
stamp "all trees done"
