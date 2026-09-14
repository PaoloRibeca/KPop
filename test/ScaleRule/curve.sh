#!/usr/bin/env bash
# curve.sh <tsv...>: the preregistered estimator for the minimum-cluster threshold.
#
# THIS IS THE SECONDARY MEASURE, AND IT ADDS NOTHING TO mstar.sh.  An attack whose clusters are
# smaller than s is charged -1 per member and scores about 0, so once every member of A is below s
# the margin is minus the reference's own G3b, and s* lands on the first defence rung above the
# largest group that wins with no guard.  The G0 filter below also hides groups of exactly s that
# lose unguarded yet beat the guarded reference, whose own score moves with s.
#
# THE ATTACK SET IS NOT EVERY PARTITION THAT OUTSCORES THE REFERENCE UNDER THE GUARD.  It is
#   A = { P : family in (degen, grad, peel) } AND { G0(P) > G0(CDC) under the same silhouette }
# -- a partition the unguarded objective already rejects is not something a threshold is needed to
# suppress.  Without the second clause the binding attacker on half the cells is a copy of the
# reference with a few points peeled off, which loses by about 0.25 unguarded and which the guard
# restores to within 0.002 of the reference, so that the guard is credited for undoing damage no
# search would propose.  Membership is independent of s by construction (G0 carries no threshold),
# so the max is over a constant set at every rung.
#
# M(s) = max_{P in A} G3b_s(P) - G3b_s(CDC);  s* = min{ s : M(s) <= -delta },  delta = 0.002.
# Reported separately against the ABSOLUTE-size and RELATIVE-size sub-corpora (column 20 is a mask:
# bit 0 absolute, bit 1 relative).  If those two s* differ by two rungs or more, the measurement is
# a property of the corpus rather than of the guard, and no hypothesis may be read off it.
#
# Outcomes are not all numbers: VOID (|A| < 5), LEFT (already protected at the bottom rung) and
# RIGHT (protected at no rung) are outcomes in their own right and are never collapsed to a number.
set -u

DELTA="${DELTA:-0.002}"
MINA="${MINA:-5}"

gawk -F'\t' -v DELTA="$DELTA" -v MINA="$MINA" -v NIN=$# '
function rung(s,   i) { for (i = 1; i <= nl; i++) if (lad[i] + 0 == s + 0) return i; return 0 }

# --- pass 1: completeness, the unguarded scores, and the ladder.  The file list is given twice, so
# the first pass is every argument up to NIN -- FNR == NR would only ever catch the first file.
ARGIND <= NIN {
  if ($3 == "end") { done[$1] = 1; next }
  if ($11 == "G0") { g0[$1, $10, $4] = $13 + 0; if ($3 == "ref") ref0[$1, $10] = $13 + 0 }
  else if (!seenlad[$12 + 0]) { seenlad[$12 + 0] = 1; lad[++nl] = $12 + 0 }
  next
}

# --- pass 2: the guarded scores
$3 == "end" { next }
$11 == "G3b" {
  t = $1; kd = $10; s = $12 + 0; v = $13 + 0; fam = $20 + 0
  N[t] = $16; FR[t] = $14; SD[t] = $15; MD[t] = $19
  if ($3 == "ref") { r[t, kd, s] = v; RN[t, kd, s] = $17 + 0; KR[t, kd] = $5 + 0; cells[t SUBSEP kd] = 1; next }
  if ($3 != "degen" && $3 != "grad" && $3 != "peel") next
  if (g0[t, kd, $4] <= ref0[t, kd]) next                 # the objective already rejects it
  if (!((t, kd, $4) in mem)) { mem[t, kd, $4] = 1; A[t, kd]++ }
  if (v > r[t, kd, s] + 1e-9) beat[t, kd, s]++
  if (!((t, kd, s) in ball) || v > ball[t, kd, s]) { ball[t, kd, s] = v; wall[t, kd, s] = $4 }
  if (and(fam, 1) && (!((t, kd, s) in babs) || v > babs[t, kd, s])) babs[t, kd, s] = v
  if (and(fam, 2) && (!((t, kd, s) in brel) || v > brel[t, kd, s])) brel[t, kd, s] = v
}
$11 == "G1" && $3 == "ref" { r1[$1, $10, $12 + 0] = $13 + 0 }
$11 == "G2" && $3 == "ref" { r2[$1, $10, $12 + 0] = $13 + 0 }
$11 == "G0" && $3 == "ref" { r0[$1, $10] = $13 + 0 }

END {
  n = asort(lad, lad, "@val_num_asc")
  nl = n
  printf "#### estimator: s* = min{ s : max_A G3b_s(P) - G3b_s(CDC) <= -%s };  |A| < %s is VOID\n", DELTA, MINA
  printf "%-30s %-10s %5s %4s %5s | %-9s %-9s %-9s | %s\n",
    "run", "silhouette", "n", "|A|", "k_ref", "s*(all)", "s*(abs)", "s*(rel)", "gate"
  for (c in cells) {
    split(c, f, SUBSEP); t = f[1]; kd = f[2]
    na = A[t, kd] + 0
    sa = verdict(t, kd, "ball"); sb = verdict(t, kd, "babs"); sc = verdict(t, kd, "brel")
    # The gate also fires when one sub-corpus gives a threshold and the other gives none at all,
    # which is a larger disagreement than two rungs rather than no disagreement
    g = "-"
    numb = (sb ~ /^[0-9]/); numc = (sc ~ /^[0-9]/)
    if (na < MINA) g = "VOID"
    else if (numb != numc) g = "GATE FIRES: one sub-corpus protected, the other never"
    else if (numb && numc && (rung(sb) - rung(sc) >= 2 || rung(sc) - rung(sb) >= 2))
      g = "GATE FIRES: corpus, not guard"
    printf "%-30s %-10s %5d %4d %5d | %-9s %-9s %-9s | %s\n",
      substr(t, 1, 30), kd, N[t], na, KR[t, kd], sa, sb, sc, g
  }
  printf "\n#### the curve, with the reference decomposed (handicap = G2-G0, charge = noise/n)\n"
  printf "%-24s %-10s %4s %8s %8s %8s %8s %8s %4s  %s\n",
    "run", "silhouette", "s", "ref G3b", "maxA", "M(s)", "hcap", "charge", "#>", "binding attacker"
  for (c in cells) {
    split(c, f, SUBSEP); t = f[1]; kd = f[2]
    if (A[t, kd] + 0 < MINA) continue
    for (i = 1; i <= nl; i++) {
      s = lad[i]
      if (!((t, kd, s) in r)) continue
      m = ((t, kd, s) in ball) ? ball[t, kd, s] - r[t, kd, s] : -99
      printf "%-24s %-10s %4d %8.4f %8s %8s %8.4f %8.4f %4d  %s\n",
        substr(t, 1, 24), kd, s, r[t, kd, s],
        ((t, kd, s) in ball) ? sprintf("%8.4f", ball[t, kd, s]) : "     -",
        (m == -99) ? "     -" : sprintf("%+8.4f", m),
        r2[t, kd, s] - r0[t, kd], RN[t, kd, s] / N[t], beat[t, kd, s] + 0,
        substr(wall[t, kd, s], 1, 40)
    }
  }
}

function verdict(t, kd, which,   i, s, m, got) {
  if (A[t, kd] + 0 < MINA) return "VOID"
  for (i = 1; i <= nl; i++) {
    s = lad[i]
    if (!((t, kd, s) in r)) continue
    if (which == "ball") { if (!((t, kd, s) in ball)) continue; m = ball[t, kd, s] - r[t, kd, s] }
    else if (which == "babs") { if (!((t, kd, s) in babs)) continue; m = babs[t, kd, s] - r[t, kd, s] }
    else { if (!((t, kd, s) in brel)) continue; m = brel[t, kd, s] - r[t, kd, s] }
    got = 1
    if (m <= -DELTA) return (i == 1) ? "2 LEFT" : sprintf("%d", s)
  }
  return got ? "RIGHT" : "none"
}
' "$@" "$@"
