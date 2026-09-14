#!/usr/bin/env bash

# Integration tests for KPop-autotuner, through the CLI.
#
# Run from the project root.  Assumes:
#   - .build/KPop-autotuner exists (from `bash BUILD release-static`)
#
# Exits 0 on full success, 1 on any failure.
#   bash test/integration_autotuner.sh --record   rewrites the golden from the binary at hand
#
# WHAT THE GOLDEN IS FOR.  Changes to the search are meant to be inert under default flags, and
# "inert" is only a claim until some file says what the default output was.  test/Autotuner/ holds
# that for one small run: its stdout, which carries the search's per-round report, and the partition
# it writes.  Recording it is a deliberate act, done from a binary whose default path is the one
# being protected, never to make a failing comparison pass.
#
# EXACT ONLY ON THE MACHINE THAT RECORDED IT.  OpenBLAS is built with DYNAMIC_ARCH, so the last bits
# of a decomposition can depend on the CPU, and macOS CI runs this suite too.  The host the golden
# was recorded on is written beside it.  Anywhere else the byte comparison is reported as skipped,
# and the run is checked for what floating point cannot move: that it succeeds, how many rounds it
# reports, and which samples the partition names.

set -u

BIN=".build/KPop-autotuner"
TWISTDB=".build/KPopTwistDB"
DATA="test/Primer/Train-5"
GOLD="test/Autotuner"
# Small enough to run in seconds, and single-threaded so that nothing but the code is compared
FLAGS=(-T 1 --iterations 2 --montecarlo-steps 200 --montecarlo-replicas 2)

if [[ ! -x "$BIN" ]]; then
  echo "FAIL: $BIN not built; run 'bash BUILD release-static' first" >&2
  exit 1
fi

# Ensure the Train-5 fixture exists (idempotent regeneration from raw FASTAs)
bash test/integration_build.sh

TMP="$(mktemp -d -t kpop-integration-autotuner-XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

failed=0
pass() { printf "  %-60s PASS\n" "$1"; }
fail() { printf "  %-60s FAIL: %s\n" "$1" "$2"; failed=1; }
skip() { printf "  %-60s SKIP: %s\n" "$1" "$2"; }

host_fingerprint() {
  local cpu
  if [[ -r /proc/cpuinfo ]]; then
    cpu="$(grep -m 1 'model name' /proc/cpuinfo | sed 's/^[^:]*: *//')"
  else
    cpu="$(sysctl -n machdep.cpu.brand_string 2> /dev/null || echo unknown)"
  fi
  printf '%s | %s\n' "$(uname -sm)" "$cpu"
}

# The number of rounds a run reported, read off the header the search prints for each
rounds() { grep -c '^# n=' "$1"; }

if [[ "${1:-}" == "--record" ]]; then
  mkdir -p "$GOLD"
  "$BIN" -i "$DATA" -o "$TMP/t" "${FLAGS[@]}" > "$GOLD/Train-5.stdout" || {
    echo "FAIL: the run to be recorded did not succeed" >&2
    exit 1
  }
  cp "$TMP/t.KPopClasses.txt" "$GOLD/Train-5.KPopClasses.txt"
  host_fingerprint > "$GOLD/Train-5.host"
  echo "Recorded $GOLD/Train-5.{stdout,KPopClasses.txt,host} ($(rounds "$GOLD/Train-5.stdout") rounds)."
  exit 0
fi

# ----------------------------------------------------------------------------
# Part 1: default flags reproduce the recorded run
# ----------------------------------------------------------------------------
echo "=== Part 1: default flags reproduce the golden ==="
"$BIN" -i "$DATA" -o "$TMP/t" "${FLAGS[@]}" > "$TMP/stdout" 2> "$TMP/stderr"
rc=$?
if [[ $rc -eq 0 ]]; then
  pass "the run succeeds"
else
  fail "the run succeeds" "exit status $rc; $(tail -n 1 "$TMP/stderr")"
fi
if [[ "$(host_fingerprint)" == "$(cat "$GOLD/Train-5.host")" ]]; then
  if cmp -s "$TMP/stdout" "$GOLD/Train-5.stdout"; then
    pass "stdout is byte-identical to the golden"
  else
    fail "stdout is byte-identical to the golden" "differs"
  fi
  if cmp -s "$TMP/t.KPopClasses.txt" "$GOLD/Train-5.KPopClasses.txt"; then
    pass "the partition is byte-identical to the golden"
  else
    fail "the partition is byte-identical to the golden" "differs"
  fi
else
  skip "byte-identical stdout and partition" "recorded on $(cat "$GOLD/Train-5.host")"
  if [[ "$(rounds "$TMP/stdout")" == "$(rounds "$GOLD/Train-5.stdout")" ]]; then
    pass "the same number of rounds as the golden"
  else
    fail "the same number of rounds as the golden" \
      "$(rounds "$TMP/stdout") against $(rounds "$GOLD/Train-5.stdout")"
  fi
  if [[ "$(head -n 1 "$TMP/t.KPopClasses.txt")" == "$(head -n 1 "$GOLD/Train-5.KPopClasses.txt")" ]]; then
    pass "the partition names the same samples as the golden"
  else
    fail "the partition names the same samples as the golden" "sample names differ"
  fi
fi

# ----------------------------------------------------------------------------
# Part 2: the calibrated detector, and what its flags refuse
# ----------------------------------------------------------------------------
echo
echo "=== Part 2: the calibrated detector ==="
"$BIN" -i "$DATA" -o "$TMP/c" "${FLAGS[@]}" --valleys-method calibrated \
  > "$TMP/c.stdout" 2> "$TMP/c.stderr"
rc=$?
if [[ $rc -eq 0 ]]; then
  pass "--valleys-method calibrated succeeds"
else
  fail "--valleys-method calibrated succeeds" "exit status $rc; $(tail -n 1 "$TMP/c.stderr")"
fi
if [[ "$(rounds "$TMP/c.stdout")" -gt 0 \
      && "$(grep -c '^# n=.* level=finest fv=[^ ]* fp=' "$TMP/c.stdout")" == "$(rounds "$TMP/c.stdout")" ]]; then
  pass "every round's header gives the level and both shares"
else
  fail "every round's header gives the level and both shares" "$(grep -m 1 '^# n=' "$TMP/c.stdout")"
fi
if ! grep -q ' level=' "$TMP/stdout"; then
  pass "and the default header does not"
else
  fail "and the default header does not" "a level appears under default flags"
fi

# A refusal has to name what it refuses, or it is some other failure passing for one
refused() {
  local what="$1" pattern="$2"
  shift 2
  if "$BIN" -i "$DATA" -o "$TMP/r" "${FLAGS[@]}" "$@" > /dev/null 2> "$TMP/r.stderr"; then
    fail "$what" "accepted"
  elif grep -q -F -- "$pattern" "$TMP/r.stderr"; then
    pass "$what"
  else
    fail "$what" "refused, but not for that: $(tail -n 1 "$TMP/r.stderr")"
  fi
}
refused "--valleys-resamples 1 is refused" "needs at least 2 resamples" --valleys-resamples 1
refused "a level share above 1 is refused" "Invalid level 'share(1.5)'" --montecarlo-level 'share(1.5)'
refused "a valleys method of no known kind is refused" "Unrecognized valleys method 'bogus'" \
  --valleys-method bogus

# ----------------------------------------------------------------------------
# Part 3: the number of axes chosen by the valleys
# ----------------------------------------------------------------------------
echo
echo "=== Part 3: --dimensions auto ==="
AUTO=(--dimensions auto --valleys-method calibrated --partition-sample 50)
"$BIN" -i "$DATA" -o "$TMP/a" "${FLAGS[@]}" "${AUTO[@]}" --dimensions-inertia 0.5 \
  > "$TMP/a.stdout" 2> "$TMP/a.stderr"
rc=$?
if [[ $rc -eq 0 ]]; then
  pass "--dimensions auto succeeds"
else
  fail "--dimensions auto succeeds" "exit status $rc; $(tail -n 1 "$TMP/a.stderr")"
fi
if grep -q -F "Option '--dimensions-inertia' is ignored under '--dimensions auto'" "$TMP/a.stderr"; then
  pass "and says it ignores --dimensions-inertia"
else
  fail "and says it ignores --dimensions-inertia" "no warning on stderr"
fi
# The register written is in the axes the last round searched in, which its header gives as D
last_d="$(grep '^=== Clustering' "$TMP/a.stdout" | tail -n 1 | sed 's/.*D=\([0-9]*\).*/\1/')"
if [[ $rc -eq 0 ]] && "$TWISTDB" -i t "$TMP/a" -O t "$TMP/a.text" > /dev/null 2>&1 \
     && [[ "$(awk -F'\t' 'NR == 1 { print NF - 1; exit }' "$TMP/a.text.KPopTwisted.txt")" == "$last_d" ]]; then
  pass "the register written has the axes the last round searched in"
else
  fail "the register written has the axes the last round searched in" "expected ${last_d:-none}"
fi
refused "--dimensions auto without the calibrated detector is refused" \
  "requires '--valleys-method calibrated'" --dimensions auto --partition-sample 50
refused "--dimensions auto with --partition-sample 0 is refused" \
  "requires '--partition-sample' above 0" --dimensions auto --valleys-method calibrated
refused "--dimensions auto with --partition-sample 1 is refused" \
  "cannot work with '--partition-sample 1'" --dimensions auto --valleys-method calibrated \
  --partition-sample 1
refused "--dimensions auto with --projection-sample 1 is refused" \
  "cannot work with '--projection-sample 1'" "${AUTO[@]}" --projection-sample 1
refused "--dimensions auto with --report-anyway is refused" "refuses '--report-anyway'" \
  "${AUTO[@]}" --report-anyway
refused "a number of dimensions of no known kind is refused" \
  "Unrecognized number of dimensions 'sideways'" --dimensions sideways

echo
if [[ $failed -eq 0 ]]; then
  echo "All autotuner integration tests passed."
  exit 0
else
  echo "Some autotuner integration tests FAILED."
  exit 1
fi
