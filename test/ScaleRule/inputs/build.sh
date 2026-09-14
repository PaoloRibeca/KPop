#!/usr/bin/env bash
# build.sh <genome> <tag>: CA of the spectra named in <genome>.<tag>.sel, then all spectra projected
set -euo pipefail
G="$1"; T="$2"
B=/home/ribeca/Development/Git/KPop.Claude/.build
R=/home/ribeca/.claude/jobs/29a54aa8/tmp/regions
cd /home/ribeca/.claude/jobs/29a54aa8/tmp/collapse/real-data
SEL="$(paste -sd, "$G.$T.sel")"
nice -n 19 $B/KPopCountDB -T 1 -i "$R/$G" -L "$SEL" -N -D -o "$G.$T.subset" 2> "$G.$T.countdb.err"
nice -n 19 $B/KPopTwist -T 1 --kmers-threshold off -v -i "$G.$T.subset" -o "$G.$T" 2> "$G.$T.twist.err"
nice -n 19 $B/KPopTwistDB -T 1 -i T "$G.$T" -t "$R/$G" -O t "$G.$T.all" 2> "$G.$T.twistdb.err"
rm -f "$G.$T.subset.KPopSpectra"
echo done "$G" "$T"
