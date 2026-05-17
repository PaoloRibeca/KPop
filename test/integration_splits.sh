#!/usr/bin/env bash

# Integration tests for the splits subsystem: full pipeline from a
# KPopTwisted binary through KPopTwistDB --splits-algorithm and on into
# Yggdrasill -t for tree assembly.
#
# Run from the project root.  Assumes:
#   - build/KPopTwistDB exists (from `bash BUILD release-static`)
#   - _build/default/BiOCamLib/bin/Yggdrasill.exe exists
#     (from `dune build --profile=release-static BiOCamLib/bin/Yggdrasill.exe`)
#
# Exits 0 on full success, 1 on any failure.

set -u

BIN="build/KPopTwistDB"
YGG="_build/default/BiOCamLib/bin/Yggdrasill.exe"
DATA="test/Primer/Classes-5"

if [[ ! -x "$BIN" ]]; then
  echo "FAIL: $BIN not built; run 'bash BUILD release-static' first" >&2
  exit 1
fi
if [[ ! -x "$YGG" ]]; then
  echo "FAIL: $YGG not built; run 'dune build --profile=release-static BiOCamLib/bin/Yggdrasill.exe' first" >&2
  exit 1
fi

TMP="$(mktemp -d -t kpop-integration-splits-XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

failed=0
pass() { printf "  %-60s PASS\n" "$1"; }
fail() { printf "  %-60s FAIL: %s\n" "$1" "$2"; failed=1; }

# ----------------------------------------------------------------------------
# Part 1: every splits-algorithm runs cleanly on the fixture
# ----------------------------------------------------------------------------
echo "=== Part 1: each splits-algorithm produces an output file ==="
for algo_args in \
    "gaps" \
    "centroids --centroids-num-seeds 5 --centroids-seed 42" \
    "hdbscan --hdbscan-min-cluster-size 2"; do
  algo="${algo_args%% *}"
  if $BIN -i t $DATA --splits-algorithm $algo_args -S "$TMP/algo_$algo" >/dev/null 2>&1 \
       && [[ -s "$TMP/algo_$algo.PhyloSplits.txt" ]]; then
    pass "splits-algorithm $algo emits non-empty file"
  else
    fail "splits-algorithm $algo" "non-zero exit or empty output"
  fi
done

# ----------------------------------------------------------------------------
# Part 2: HDBSCAN output has 0 dropped weight in Yggdrasill
# ----------------------------------------------------------------------------
echo "=== Part 2: HDBSCAN splits are jointly compatible (0 dropped weight) ==="
for k in 1 2 3 4; do
  $BIN -i t $DATA --splits-algorithm hdbscan --hdbscan-min-cluster-size $k \
       -S "$TMP/hdb_K$k" >/dev/null 2>&1
  out="$($YGG -I "$TMP/hdb_K$k" -t "$TMP/hdb_K$k.tree" 2>&1)"
  dropped="$(printf '%s\n' "$out" | grep -oE '[0-9]+ incompatible splits dropped' | grep -oE '^[0-9]+')"
  if [[ "$dropped" == "0" ]]; then
    pass "HDBSCAN K=$k: 0 dropped splits in Yggdrasill"
  else
    fail "HDBSCAN K=$k" "dropped=$dropped (expected 0)"
  fi
done

# ----------------------------------------------------------------------------
# Part 3: empty-splits round-trip through Yggdrasill
# HDBSCAN with K=5 on the 10-sample fixture cannot emit any cluster; the
# splits file is "names;" (no colon), and Yggdrasill must still parse it
# and produce a star tree.
# ----------------------------------------------------------------------------
echo "=== Part 3: empty splits file parses cleanly ==="
$BIN -i t $DATA --splits-algorithm hdbscan --hdbscan-min-cluster-size 5 \
     -S "$TMP/empty" >/dev/null 2>&1
if grep -qE 'C1.*C2.*C10;$' "$TMP/empty.PhyloSplits.txt"; then
  pass "empty splits file emits the 'names;' shape"
else
  fail "empty splits file" "expected 'names;' format, got: $(cat "$TMP/empty.PhyloSplits.txt")"
fi
if $YGG -I "$TMP/empty" -t "$TMP/empty.tree" >/dev/null 2>&1 \
     && [[ -s "$TMP/empty.tree.nwk" ]]; then
  pass "Yggdrasill accepts the empty-splits file"
else
  fail "Yggdrasill on empty splits" "parse error or empty tree file"
fi

# ----------------------------------------------------------------------------
# Part 4: cross-thread reproducibility for deterministic paths
# ----------------------------------------------------------------------------
echo "=== Part 4: cross-thread reproducibility (deterministic modes only) ==="
$BIN -i t $DATA -T 1 --splits-algorithm centroids --centroids-num-seeds 10 \
     --centroids-seed 42 -S "$TMP/cen_T1" >/dev/null 2>&1
$BIN -i t $DATA -T 4 --splits-algorithm centroids --centroids-num-seeds 10 \
     --centroids-seed 42 -S "$TMP/cen_T4" >/dev/null 2>&1
if cmp -s "$TMP/cen_T1.PhyloSplits.txt" "$TMP/cen_T4.PhyloSplits.txt"; then
  pass "centroids: -T 1 == -T 4 (seed=42)"
else
  fail "centroids -T 1 vs -T 4" "outputs differ"
fi

$BIN -i t $DATA -T 1 --splits-algorithm hdbscan --hdbscan-mst-mode dense \
     -S "$TMP/dense_T1" >/dev/null 2>&1
$BIN -i t $DATA -T 4 --splits-algorithm hdbscan --hdbscan-mst-mode dense \
     -S "$TMP/dense_T4" >/dev/null 2>&1
if cmp -s "$TMP/dense_T1.PhyloSplits.txt" "$TMP/dense_T4.PhyloSplits.txt"; then
  pass "HDBSCAN dense: -T 1 == -T 4"
else
  fail "HDBSCAN dense -T 1 vs -T 4" "outputs differ"
fi

$BIN -i t $DATA -T 1 --splits-algorithm hdbscan --hdbscan-index-type flat \
     -S "$TMP/aflat_T1" >/dev/null 2>&1
$BIN -i t $DATA -T 4 --splits-algorithm hdbscan --hdbscan-index-type flat \
     -S "$TMP/aflat_T4" >/dev/null 2>&1
if cmp -s "$TMP/aflat_T1.PhyloSplits.txt" "$TMP/aflat_T4.PhyloSplits.txt"; then
  pass "HDBSCAN auto(flat): -T 1 == -T 4"
else
  fail "HDBSCAN auto(flat) -T 1 vs -T 4" "outputs differ"
fi

# ----------------------------------------------------------------------------
# Part 5: dense and auto-with-flat are byte-identical on small data
# ----------------------------------------------------------------------------
echo "=== Part 5: HDBSCAN Auto(flat) == Dense byte-for-byte (small n) ==="
for k in 1 2 3 4; do
  $BIN -i t $DATA --splits-algorithm hdbscan --hdbscan-min-cluster-size $k \
       --hdbscan-mst-mode dense -S "$TMP/d_$k" >/dev/null 2>&1
  $BIN -i t $DATA --splits-algorithm hdbscan --hdbscan-min-cluster-size $k \
       --hdbscan-mst-mode auto --hdbscan-index-type flat \
       -S "$TMP/a_$k" >/dev/null 2>&1
  if cmp -s "$TMP/d_$k.PhyloSplits.txt" "$TMP/a_$k.PhyloSplits.txt"; then
    pass "HDBSCAN K=$k: auto(flat) == dense"
  else
    fail "HDBSCAN K=$k: auto(flat) vs dense" "outputs differ"
  fi
done

# ----------------------------------------------------------------------------
# Part 6: sparse mode rejects under-sized num_neighbors with a clear message
# ----------------------------------------------------------------------------
echo "=== Part 6: sparse mode validates --hdbscan-num-neighbors ==="
out="$($BIN -i t $DATA --splits-algorithm hdbscan --hdbscan-num-neighbors 1 \
       --hdbscan-mst-mode sparse --hdbscan-min-samples 5 \
       -S "$TMP/bad" 2>&1 || true)"
if printf '%s' "$out" | grep -q 'must be >= --hdbscan-min-samples'; then
  pass "sparse mode rejects num_neighbors < min_samples with helpful message"
else
  fail "sparse-mode validation" "did not raise expected error message"
fi

echo
if [[ $failed -eq 0 ]]; then
  echo "All splits-integration tests passed."
  exit 0
else
  echo "Some splits-integration tests FAILED."
  exit 1
fi
