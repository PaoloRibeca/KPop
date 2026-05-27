#!/usr/bin/env bash

# Integration tests for the phylogenetic-tree subsystem: full pipeline
# from a KPopTwisted binary through KPopTwistDB --phylo-method and on
# into a Newick (.nwk) tree.
#
# Run from the project root.  Assumes:
#   - build/KPopTwistDB exists (from `bash BUILD release-static`)
#
# Exits 0 on full success, 1 on any failure.

set -u

BIN="build/KPopTwistDB"
DATA="test/Primer/Classes-5"

if [[ ! -x "$BIN" ]]; then
  echo "FAIL: $BIN not built; run 'bash BUILD release-static' first" >&2
  exit 1
fi

# Ensure the Classes-5 fixture exists (idempotent regeneration from raw FASTAs)
bash test/integration_build.sh

TMP="$(mktemp -d -t kpop-integration-phylo-XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

failed=0
pass() { printf "  %-60s PASS\n" "$1"; }
fail() { printf "  %-60s FAIL: %s\n" "$1" "$2"; failed=1; }

# ----------------------------------------------------------------------------
# Part 1: every phylo-algorithm runs cleanly on the fixture
# ----------------------------------------------------------------------------
echo "=== Part 1: each phylo-algorithm produces a Newick tree ==="
for algo_args in \
    "gaps" \
    "centroids --phylo-centroids-num-seeds 5 --phylo-centroids-seed 42" \
    "hdbscan --phylo-hdbscan-min-cluster-size 2" \
    "sparse-nj"; do
  algo="${algo_args%% *}"
  if $BIN -i t $DATA --phylo-method $algo_args -P "$TMP/algo_$algo" >/dev/null 2>&1 \
       && [[ -s "$TMP/algo_$algo.nwk" ]]; then
    pass "phylo-algorithm $algo emits non-empty .nwk file"
  else
    fail "phylo-algorithm $algo" "non-zero exit or empty output"
  fi
done

# ----------------------------------------------------------------------------
# Part 2: HDBSCAN trees parse as Newick
# ----------------------------------------------------------------------------
echo "=== Part 2: HDBSCAN .nwk output is well-formed ==="
for k in 1 2 3 4; do
  $BIN -i t $DATA --phylo-method hdbscan --phylo-hdbscan-min-cluster-size $k \
       -P "$TMP/hdb_K$k" >/dev/null 2>&1
  if [[ -s "$TMP/hdb_K$k.nwk" ]] && grep -qE '\);[[:space:]]*$' "$TMP/hdb_K$k.nwk"; then
    pass "HDBSCAN K=$k: Newick file ends with ');'"
  else
    fail "HDBSCAN K=$k" "missing or malformed .nwk output"
  fi
done

# ----------------------------------------------------------------------------
# Part 3: HDBSCAN with K too large for fixture still emits a (degenerate) tree
# ----------------------------------------------------------------------------
echo "=== Part 3: HDBSCAN with K=5 on 10-sample fixture emits a tree ==="
$BIN -i t $DATA --phylo-method hdbscan --phylo-hdbscan-min-cluster-size 5 \
     -P "$TMP/empty" >/dev/null 2>&1
if [[ -s "$TMP/empty.nwk" ]] && grep -qE '\);[[:space:]]*$' "$TMP/empty.nwk"; then
  pass "K=5 HDBSCAN emits a valid Newick file"
else
  fail "K=5 HDBSCAN" "no valid .nwk produced"
fi

# ----------------------------------------------------------------------------
# Part 4: cross-thread reproducibility for deterministic paths
# ----------------------------------------------------------------------------
echo "=== Part 4: cross-thread reproducibility (deterministic modes only) ==="
$BIN -i t $DATA -T 1 --phylo-method centroids --phylo-centroids-num-seeds 10 \
     --phylo-centroids-seed 42 -P "$TMP/cen_T1" >/dev/null 2>&1
$BIN -i t $DATA -T 4 --phylo-method centroids --phylo-centroids-num-seeds 10 \
     --phylo-centroids-seed 42 -P "$TMP/cen_T4" >/dev/null 2>&1
if cmp -s "$TMP/cen_T1.nwk" "$TMP/cen_T4.nwk"; then
  pass "centroids: -T 1 == -T 4 (seed=42)"
else
  fail "centroids -T 1 vs -T 4" "outputs differ"
fi

$BIN -i t $DATA -T 1 --phylo-method hdbscan --phylo-hdbscan-mst-mode dense \
     -P "$TMP/dense_T1" >/dev/null 2>&1
$BIN -i t $DATA -T 4 --phylo-method hdbscan --phylo-hdbscan-mst-mode dense \
     -P "$TMP/dense_T4" >/dev/null 2>&1
if cmp -s "$TMP/dense_T1.nwk" "$TMP/dense_T4.nwk"; then
  pass "HDBSCAN dense: -T 1 == -T 4"
else
  fail "HDBSCAN dense -T 1 vs -T 4" "outputs differ"
fi

$BIN -i t $DATA -T 1 --phylo-method hdbscan --phylo-hdbscan-index-type flat \
     -P "$TMP/aflat_T1" >/dev/null 2>&1
$BIN -i t $DATA -T 4 --phylo-method hdbscan --phylo-hdbscan-index-type flat \
     -P "$TMP/aflat_T4" >/dev/null 2>&1
if cmp -s "$TMP/aflat_T1.nwk" "$TMP/aflat_T4.nwk"; then
  pass "HDBSCAN auto(flat): -T 1 == -T 4"
else
  fail "HDBSCAN auto(flat) -T 1 vs -T 4" "outputs differ"
fi

# ----------------------------------------------------------------------------
# Part 5: dense and auto-with-flat are byte-identical on small data
# (Persistence-mode lengths; Mreach walks the raw merge tree and is
#  sensitive to MST tie-breaking when leaves merge at coincident lambda.)
# ----------------------------------------------------------------------------
echo "=== Part 5: HDBSCAN Auto(flat) == Dense byte-for-byte (small n) ==="
for k in 1 2 3 4; do
  $BIN -i t $DATA --phylo-method hdbscan --phylo-hdbscan-min-cluster-size $k \
       --phylo-hdbscan-lengths persistence \
       --phylo-hdbscan-mst-mode dense -P "$TMP/d_$k" >/dev/null 2>&1
  $BIN -i t $DATA --phylo-method hdbscan --phylo-hdbscan-min-cluster-size $k \
       --phylo-hdbscan-lengths persistence \
       --phylo-hdbscan-mst-mode auto --phylo-hdbscan-index-type flat \
       -P "$TMP/a_$k" >/dev/null 2>&1
  if cmp -s "$TMP/d_$k.nwk" "$TMP/a_$k.nwk"; then
    pass "HDBSCAN K=$k: auto(flat) == dense"
  else
    fail "HDBSCAN K=$k: auto(flat) vs dense" "outputs differ"
  fi
done

# ----------------------------------------------------------------------------
# Part 6: HDBSCAN persistence and mreach branch-length modes both work
# ----------------------------------------------------------------------------
echo "=== Part 6: HDBSCAN persistence/mreach branch-length modes ==="
for mode in persistence mreach; do
  $BIN -i t $DATA --phylo-method hdbscan --phylo-hdbscan-min-cluster-size 2 \
       --phylo-hdbscan-lengths $mode -P "$TMP/lm_$mode" >/dev/null 2>&1
  if [[ -s "$TMP/lm_$mode.nwk" ]] && grep -qE '\);[[:space:]]*$' "$TMP/lm_$mode.nwk"; then
    pass "HDBSCAN lengths=$mode emits valid .nwk"
  else
    fail "HDBSCAN lengths=$mode" "missing or malformed .nwk"
  fi
done

# ----------------------------------------------------------------------------
# Part 7: sparse mode rejects under-sized num_neighbors with a clear message
# ----------------------------------------------------------------------------
echo "=== Part 7: sparse mode validates --phylo-hdbscan-num-neighbors ==="
out="$($BIN -i t $DATA --phylo-method hdbscan --phylo-hdbscan-num-neighbors 1 \
       --phylo-hdbscan-mst-mode sparse --phylo-hdbscan-min-samples 5 \
       -P "$TMP/bad" 2>&1 || true)"
if printf '%s' "$out" | grep -q 'must be >= --phylo-hdbscan-min-samples'; then
  pass "sparse mode rejects num_neighbors < min_samples with helpful message"
else
  fail "sparse-mode validation" "did not raise expected error message"
fi

echo
if [[ $failed -eq 0 ]]; then
  echo "All phylo-integration tests passed."
  exit 0
else
  echo "Some phylo-integration tests FAILED."
  exit 1
fi
