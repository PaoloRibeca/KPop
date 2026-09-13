#!/usr/bin/env bash
# run2.sh: every embedding in turn with the fixed binary, one program at a time
set -u
G=/home/ribeca/.claude/jobs/29a54aa8/tmp/guard2
R=/home/ribeca/.claude/jobs/29a54aa8/tmp/collapse/real-data
V=/home/ribeca/.claude/jobs/29a54aa8/tmp/simpson/vp2
L=/home/ribeca/Development/Git/NailIt/test/data/denovo
run () { # <twisted prefix> <labels> <axes> <tag>
  echo "=== $4 started $(date +%T)"
  nice -n 19 "$G/guard2" "$1.all.KPopTwisted.txt" "$1.all.KPopInertia.txt" "$2" "$3" "$4" \
    "$G/guard2.$4.tsv" > "$G/guard2.$4.out" 2> "$G/guard2.$4.err"
  echo "=== $4 exit $? $(date +%T)"
}
run $R/vp1.T_rand  $L/vp1.cdc   20 vp1.d20
run $R/vp1.T_rand  $L/vp1.cdc  243 vp1.d243
run $R/rdrp.T_rand $L/rdrp.cdc   5 rdrp.d5
run $R/rdrp.T_rand $L/rdrp.cdc 199 rdrp.d199
run $V/vp2.T_rand  $L/vp2.cdc   20 vp2.d20
run $V/vp2.T_rand  $L/vp2.cdc  219 vp2.d219
echo "=== all done $(date +%T)"
