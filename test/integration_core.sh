#!/usr/bin/env bash

# Core-pipeline regression test: replays the Quick Start tutorial from the
# README, both without and with k-mer selection.  Catches breakage in the
# KPopCount -> KPopCountDB -> KPopTwist -> KPopTwistDB chain that the
# narrower splits-subsystem tests wouldn't surface.
#
# Run from the project root.  Assumes:
#   - .build/{KPopCount,KPopCountDB,KPopTwist,KPopTwistDB} exist
#     (from `bash BUILD release-static`)
#
# Exits 0 on success, 1 on any failure.

set -u

COUNT=".build/KPopCount"
COUNTDB=".build/KPopCountDB"
TWIST=".build/KPopTwist"
TWISTDB=".build/KPopTwistDB"
PRIMER="test/Primer"

for bin in "$COUNT" "$COUNTDB" "$TWIST" "$TWISTDB"; do
  if [[ ! -x "$bin" ]]; then
    echo "FAIL: $bin not built; run 'bash BUILD release-static' first" >&2
    exit 1
  fi
done
if [[ ! -d "$PRIMER/Train" || ! -d "$PRIMER/Test" ]]; then
  echo "FAIL: tutorial fixture missing under $PRIMER" >&2
  exit 1
fi

TMP="$(mktemp -d -t kpop-integration-core-XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

K=5
failed=0
pass() { printf "  %-60s PASS\n" "$1"; }
fail() { printf "  %-60s FAIL: %s\n" "$1" "$2"; failed=1; }

# Count misclassifications by parsing the KPopSummary: column 1 is the
# query name (e.g. "S2-C1"), column 6 is the predicted top-class label.
count_misclassified() {
  local summary="$1"
  awk -F '\t' 'BEGIN{c=0} {split($1,s,"-"); if (s[2]!=$6) ++c} END{print c}' "$summary"
}

# ----------------------------------------------------------------------------
# Pass 1: tutorial pipeline without k-mer selection
# ----------------------------------------------------------------------------
echo "=== Pass 1: tutorial pipeline, no k-mer selection (K=$K) ==="

$COUNT -k $K -f "" "$PRIMER/Train/Train.fasta" -o "$TMP/Train" >"$TMP/log1a.out" 2>"$TMP/log1a" \
  || { fail "KPopCount on Train.fasta" "exit $?; see $TMP/log1a"; }
[[ -s "$TMP/Train.KPopSpectra" ]] \
  && pass "KPopCount produced Train.KPopSpectra ($(stat -c%s "$TMP/Train.KPopSpectra") bytes)" \
  || fail "Train.KPopSpectra" "missing or empty"

$COUNTDB -i "$TMP/Train" -m "$PRIMER/Train/CLASS.txt" -c CLASS -o /dev/stdout 2>"$TMP/log1b" \
  | $TWIST -i /dev/stdin -o "$TMP/Classes" 2>"$TMP/log1c" \
  || { fail "KPopCountDB | KPopTwist" "see $TMP/log1{b,c}"; }
[[ -s "$TMP/Classes.KPopTwister" && -s "$TMP/Classes.KPopTwisted" ]] \
  && pass "KPopTwist produced Classes.KPopTwister + KPopTwisted" \
  || fail "Classes.KPopTwister/KPopTwisted" "one or both missing"

$COUNT -k $K -f "" "$PRIMER/Test/Test.fasta" -o /dev/stdout 2>"$TMP/log1d" \
  | $TWISTDB -i T "$TMP/Classes" -t /dev/stdin -d "$TMP/Classes" "$TMP/Test" 2>"$TMP/log1e" \
  || { fail "KPopCount | KPopTwistDB project+distances" "see $TMP/log1{d,e}"; }
[[ -s "$TMP/Test.KPopSummary.txt" ]] \
  && pass "KPopTwistDB produced Test.KPopSummary.txt" \
  || fail "Test.KPopSummary.txt" "missing or empty"

if [[ -s "$TMP/Test.KPopSummary.txt" ]]; then
  miscls=$(count_misclassified "$TMP/Test.KPopSummary.txt")
  if [[ "$miscls" == "0" ]]; then
    pass "Pass 1: 0 misclassified test sequences"
  else
    fail "Pass 1 misclassification count" "got $miscls, expected 0"
  fi
fi

# ----------------------------------------------------------------------------
# Pass 2: re-train with k-mer selection derived from the Pass 1 twister
# ----------------------------------------------------------------------------
echo
echo "=== Pass 2: tutorial pipeline, with k-mer selection (--keep) ==="

# Use the Pass 1 twister to derive a representative k-mer set via
# leader-clustering of k-mer standard coordinates (the documented
# --clusters T workflow in lib/Clustering.ml).
$TWISTDB -i T "$TMP/Classes" --clusters T "$TMP/Classes.kept_kmers" \
  >"$TMP/log2a.out" 2>"$TMP/log2a" \
  || { fail "KPopTwistDB --clusters T" "see $TMP/log2a"; }
[[ -s "$TMP/Classes.kept_kmers" ]] \
  && pass "k-mer selection: $(wc -l < "$TMP/Classes.kept_kmers") representative k-mers extracted" \
  || fail "Classes.kept_kmers" "missing or empty"

# Re-twist using only the selected k-mers
$COUNTDB -i "$TMP/Train" -m "$PRIMER/Train/CLASS.txt" -c CLASS -o /dev/stdout 2>"$TMP/log2b" \
  | $TWIST -i /dev/stdin --keep "$TMP/Classes.kept_kmers" -o "$TMP/ClassesK" 2>"$TMP/log2c" \
  || { fail "KPopTwist --keep ..." "see $TMP/log2{b,c}"; }
[[ -s "$TMP/ClassesK.KPopTwister" && -s "$TMP/ClassesK.KPopTwisted" ]] \
  && pass "KPopTwist --keep produced ClassesK.KPopTwister + KPopTwisted" \
  || fail "ClassesK.KPopTwister/KPopTwisted" "one or both missing"

# Project test data against the reduced twister
$COUNT -k $K -f "" "$PRIMER/Test/Test.fasta" -o /dev/stdout 2>"$TMP/log2d" \
  | $TWISTDB -i T "$TMP/ClassesK" -t /dev/stdin -d "$TMP/ClassesK" "$TMP/TestK" 2>"$TMP/log2e" \
  || { fail "KPopTwistDB project+distances (selected k-mers)" "see $TMP/log2{d,e}"; }
[[ -s "$TMP/TestK.KPopSummary.txt" ]] \
  && pass "KPopTwistDB produced TestK.KPopSummary.txt" \
  || fail "TestK.KPopSummary.txt" "missing or empty"

if [[ -s "$TMP/TestK.KPopSummary.txt" ]]; then
  misclsK=$(count_misclassified "$TMP/TestK.KPopSummary.txt")
  if [[ "$misclsK" == "0" ]]; then
    pass "Pass 2: 0 misclassified test sequences"
  else
    # A small handful of misclassifications is plausible with selection;
    # flag anything above a generous threshold as a regression.
    threshold=10
    if (( misclsK <= threshold )); then
      pass "Pass 2: $misclsK misclassified (<= $threshold tolerated under k-mer selection)"
    else
      fail "Pass 2 misclassification count" "got $misclsK, expected <= $threshold"
    fi
  fi
fi

# ----------------------------------------------------------------------------
# Pass 3: sample clustering with --clusters-greedy-epsilon adaptive
# ----------------------------------------------------------------------------
echo
echo "=== Pass 3: --clusters t with adaptive epsilon ==="

$TWISTDB -i t "$TMP/Classes" --clusters-greedy-epsilon adaptive \
  --clusters t "$TMP/Classes.adaptive.classes" >"$TMP/log3.out" 2>"$TMP/log3" \
  || { fail "KPopTwistDB --clusters t adaptive" "see $TMP/log3"; }
if [[ -s "$TMP/Classes.adaptive.classes" ]] \
     && [[ "$(wc -l < "$TMP/Classes.adaptive.classes")" == "2" ]] \
     && head -1 "$TMP/Classes.adaptive.classes" | grep -q $'\t' \
     && tail -1 "$TMP/Classes.adaptive.classes" | grep -q '^CLASS'; then
  n_reps=$(tail -1 "$TMP/Classes.adaptive.classes" | tr '\t' '\n' | grep '^C@' | sort -u | wc -l)
  pass "adaptive sample clustering produced 2-line class file ($n_reps representatives)"
else
  fail "adaptive class file" "missing or malformed"
fi

# ----------------------------------------------------------------------------
# Pass 4: sample clustering with --clusters-method hdbscan
# Covers (a) HDBSCAN produces a non-trivial partition and
# (b) metric pre-scaling matters (flat vs powers gives different assignments).
# ----------------------------------------------------------------------------
echo
echo "=== Pass 4: --clusters t with --clusters-method hdbscan ==="

$TWISTDB -i t "$TMP/Classes" --clusters-method hdbscan \
  --clusters-hdbscan-min-cluster-size 2 \
  --clusters t "$TMP/Classes.hdb.classes" >"$TMP/log4.out" 2>"$TMP/log4" \
  || { fail "KPopTwistDB --clusters t hdbscan" "see $TMP/log4"; }
if [[ -s "$TMP/Classes.hdb.classes" ]] \
     && [[ "$(wc -l < "$TMP/Classes.hdb.classes")" == "2" ]] \
     && tail -1 "$TMP/Classes.hdb.classes" | grep -qE 'C@[0-9]+'; then
  n_clusters=$(tail -1 "$TMP/Classes.hdb.classes" | tr '\t' '\n' | grep -oE '^C@[0-9]+$' | sort -u | wc -l)
  pass "HDBSCAN clustering produced 2-line class file ($n_clusters clusters)"
else
  fail "HDBSCAN class file" "missing or malformed"
fi

$TWISTDB -i t "$TMP/Classes" -m flat --clusters-method hdbscan \
  --clusters-hdbscan-min-cluster-size 2 \
  --clusters t "$TMP/Classes.hdb.flat.classes" >"$TMP/log4f.out" 2>"$TMP/log4f"
$TWISTDB -i t "$TMP/Classes" -m 'powers(1,1,1)' --clusters-method hdbscan \
  --clusters-hdbscan-min-cluster-size 2 \
  --clusters t "$TMP/Classes.hdb.powers.classes" >"$TMP/log4p.out" 2>"$TMP/log4p"
if ! cmp -s "$TMP/Classes.hdb.flat.classes" "$TMP/Classes.hdb.powers.classes"; then
  pass "metric pre-scaling matters for HDBSCAN (flat vs powers differ)"
else
  fail "metric pre-scaling for HDBSCAN" \
    "flat and powers(1,1,1) produced identical HDBSCAN clusters (unexpected)"
fi

echo
if [[ $failed -eq 0 ]]; then
  echo "All core-pipeline tests passed."
  exit 0
else
  echo "Some core-pipeline tests FAILED."
  exit 1
fi
