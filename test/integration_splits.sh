#!/usr/bin/env bash

# Integration tests for the phylogenetic-tree subsystem, through the CLI:
# a twisted register in, a Newick (.nwk) tree out.
#
# Run from the project root.  Assumes:
#   - .build/KPopPhylo exists (from `bash BUILD release-static`)
#
# Exits 0 on full success, 1 on any failure.
#
# WHAT THIS SUITE IS FOR, now that there are two of them.  test/Phylo.ml calls
# Twisted.get_phylo_tree in process and can pass every knob the library has;
# this one drives the binary and can only pass what the CLI exposes.  So the
# division is not the same checks twice: here is what only the CLI reaches --
# that a file is written where -o says, that it is syntactically terminated,
# that argv errors are caught -- while the algorithmic invariants needing knobs
# KPopPhylo does not expose live in test/Phylo.ml Parts 2, 5 and 6.
#
# THREE CHECKS WERE DROPPED RATHER THAN PORTED, and this is where that is said
# instead of being left as a gap somebody rediscovers.  The tree surface moved
# out of KPopTwistDB into KPopPhylo, and the HDBSCAN mst-mode, lengths,
# index-type, num-neighbors and min-samples knobs, and the centroids seeding,
# have no CLI entry at all today.  The old Part 5 (auto(flat) == dense) and
# Part 6 (persistence vs mreach lengths) distinguished their two sides with
# exactly those flags: dropping the flags leaves `cmp` comparing a command with
# itself, which passes forever and tests nothing.  Both invariants are asserted
# in test/Phylo.ml Parts 2 and 5.  The old Part 7 asserted the message from the
# sparse-mreach guard in lib/Clustering.ml; that guard needs three of the
# missing knobs, and its text still names the old --phylo-hdbscan-* options, so
# a grep for it would pass while documenting a bug.  test/Phylo.ml Part 6
# provokes the same guard in process.

set -u

BIN=".build/KPopPhylo"
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

# EVERY RUN WRITES TO A PATH OF ITS OWN, EMPTIED FIRST.  Most checks below judge
# the file rather than the exit status, so a run that died would otherwise
# inherit the previous one's output and pass on it.
run() { # <output prefix> <args...>
  local out="$1"; shift
  rm -f "$TMP/$out.nwk"
  $BIN -i t "$DATA" "$@" -o "$TMP/$out" > /dev/null 2>&1
}

# Exists, is not empty, and ends as Newick does.  The leading [&U] comes from
# Trees.Newick's rich format, so the anchor is on the tail and not the head.
well_formed() { # <output prefix>
  [[ -s "$TMP/$1.nwk" ]] && grep -qE '\);[[:space:]]*$' "$TMP/$1.nwk"
}

# ----------------------------------------------------------------------------
# Part 1: every method runs cleanly and writes the file -o names
# ----------------------------------------------------------------------------
echo "=== Part 1: each method produces a Newick tree where -o says ==="
for method in gaps centroids hdbscan sparse-nj; do
  if run "algo_$method" --method "$method" && well_formed "algo_$method"; then
    pass "method $method: exit 0, well-formed tree at the -o prefix"
  else
    fail "method $method" "non-zero exit, or no well-formed .nwk at the -o prefix"
  fi
done

# ----------------------------------------------------------------------------
# Part 2: HDBSCAN across the cluster-size sweep
# ----------------------------------------------------------------------------
# --method hdbscan is named on every run: --hdbscan-min-cluster-size is honoured
# only under that method, and the default is sparse-nj, so omitting it would
# turn this sweep into four identical runs of a different algorithm.
echo "=== Part 2: HDBSCAN .nwk output is well-formed across K ==="
for k in 1 2 3 4; do
  run "hdb_K$k" --method hdbscan --hdbscan-min-cluster-size "$k"
  if well_formed "hdb_K$k"; then
    pass "HDBSCAN K=$k: Newick file ends with ');'"
  else
    fail "HDBSCAN K=$k" "missing or malformed .nwk output"
  fi
done

# ----------------------------------------------------------------------------
# Part 3: K=5 on the 10-leaf fixture -- half the leaves in one cluster
# ----------------------------------------------------------------------------
echo "=== Part 3: HDBSCAN with K=5 on the 10-leaf fixture emits a tree ==="
run empty --method hdbscan --hdbscan-min-cluster-size 5
if well_formed empty; then
  pass "K=5 HDBSCAN emits a valid Newick file"
else
  fail "K=5 HDBSCAN" "no valid .nwk produced"
fi

# ----------------------------------------------------------------------------
# Part 4: the same tree whatever the thread count
# ----------------------------------------------------------------------------
# Only pairs measured byte-identical belong under `cmp`.  Two runs differing in
# a genuine parameter can agree on splits and still order children differently,
# which `cmp` would call a failure; that comparison belongs in test/Phylo.ml,
# which compares split sets rather than bytes.
echo "=== Part 4: cross-thread reproducibility ==="
for method in centroids hdbscan; do
  run "t1_$method" --method "$method" -T 1
  run "t4_$method" --method "$method" -T 4
  if cmp -s "$TMP/t1_$method.nwk" "$TMP/t4_$method.nwk"; then
    pass "$method: -T 1 and -T 4 give byte-identical trees"
  else
    fail "$method -T 1 vs -T 4" "outputs differ"
  fi
done

# ----------------------------------------------------------------------------
# Part 5: the FAISS index is an implementation detail of sparse-NJ
# ----------------------------------------------------------------------------
echo "=== Part 5: sparse-NJ flat == hnsw(32) on the small fixture ==="
run snj_flat --method sparse-nj --snj-index-type flat
run snj_hnsw --method sparse-nj --snj-index-type "hnsw(32)"
if cmp -s "$TMP/snj_flat.nwk" "$TMP/snj_hnsw.nwk"; then
  pass "sparse-NJ: flat == hnsw(32) byte-for-byte (small n)"
else
  fail "sparse-NJ flat vs hnsw(32)" "outputs differ"
fi

# ----------------------------------------------------------------------------
# Part 6: argv errors, which only the CLI can get wrong
# ----------------------------------------------------------------------------
# Diagnostics go to stderr and carry ANSI colour even when piped, so these match
# the message body rather than anchoring on a prefix.
echo "=== Part 6: the command line is validated ==="
out="$($BIN -o "$TMP/no_input" 2>&1 || true)"
if printf '%s' "$out" | grep -q "mandatory"; then
  pass "a missing -i is refused as mandatory"
else
  fail "missing -i" "did not report a mandatory option"
fi

out="$($BIN -o "$TMP/wrong_order" -i t "$DATA" 2>&1 || true)"
if printf '%s' "$out" | grep -q "needs at least one of"; then
  pass "-o before -i is refused, actions running in the order given"
else
  fail "-o before -i" "did not report the missing register"
fi

out="$($BIN -i t "$TMP/does_not_exist" -o "$TMP/missing" 2>&1 || true)"
if printf '%s' "$out" | grep -q "Input file not found"; then
  pass "a missing register is reported by name"
else
  fail "missing register" "did not report the absent input file"
fi

echo
if [[ $failed -eq 0 ]]; then
  echo "All phylo-integration tests passed."
  exit 0
else
  echo "Some phylo-integration tests FAILED."
  exit 1
fi
