#!/usr/bin/env bash

# Idempotent builder of the binary fixtures needed by the integration tests:
#   test/Primer/Train-5.KPopSpectra    -- consumed by test/RSVD.exe
#   test/Primer/Classes-5.KPopTwister  -- consumed by test/Splits.exe
#   test/Primer/Classes-5.KPopTwisted     and test/integration_splits.sh
# Run from the project root.  Skips work when the artefacts already exist.

set -eu

PRIMER="test/Primer"
COUNT="build/KPopCount"
COUNTDB="build/KPopCountDB"
TWIST="build/KPopTwist"
K=5

if [[ -s "$PRIMER/Train-$K.KPopSpectra" \
     && -s "$PRIMER/Classes-$K.KPopTwister" \
     && -s "$PRIMER/Classes-$K.KPopTwisted" ]]; then
  exit 0
fi

for bin in "$COUNT" "$COUNTDB" "$TWIST"; do
  if [[ ! -x "$bin" ]]; then
    echo "FAIL: $bin not built; run 'bash BUILD release-static' first" >&2
    exit 1
  fi
done

# Train-5.KPopSpectra is kept on disk because RSVD.exe reads it directly;
# Classes-5.{KPopTwister,KPopTwisted} pipe straight through from it
"$COUNT" -k $K -f "" "$PRIMER/Train/Train.fasta" -o "$PRIMER/Train-$K"
"$COUNTDB" -i "$PRIMER/Train-$K" -m "$PRIMER/Train/CLASS.txt" -c CLASS -o /dev/stdout \
  | "$TWIST" -i /dev/stdin -o "$PRIMER/Classes-$K"
