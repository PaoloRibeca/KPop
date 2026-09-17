# replay.awk -- alternative axis-ladder rung rules, replayed on the stage-1 auto rounds
#
#   gawk -F'\t' -f replay.awk rungs.tsv valleys.tsv picks.tsv fixed_outcomes.tsv
#
# Every rule sees, per round, the rungs in ladder order and each rung's kept valleys (share, z,
# usable).  Output: one line per rule and round (the pick), then per rule and corpus the summary
# metrics.  Rounds 2-8 are path-dependent; round 1 is the clean counterfactual.

function first_plateau(k, sc,   i, best, n) {
  n = nr[k]; best = 0
  for (i = 1; i <= n; i++) if (sc[i] > best) best = sc[i]
  if (best <= 0) return 0
  i = 1; while (sc[i] != best) i++
  while (i < n && sc[i + 1] == best) i++
  return rung[k, i]
}
function last_max(k, sc,   i, best, n, p) {
  n = nr[k]; best = 0; p = 0
  for (i = 1; i <= n; i++) if (sc[i] > best) best = sc[i]
  if (best <= 0) return 0
  for (i = 1; i <= n; i++) if (sc[i] == best) p = rung[k, i]
  return p
}
# Rule scores for rung index i of round k, filling sc[]
# A rule is a base scoring (current, z10, z20, sumz, maxz, persist) optionally joined with '+' to
# modifiers: floorN and ceilN, which give no score to rungs outside [N, M], and ties-right, which
# takes the last rung at the maximum instead of the right end of the first plateau
function scores(k, rule,   i, j, v, s, m, r, near, u, w, base, lo, top, g) {
  delete sc
  base = rule; sub(/\+.*/, "", base)
  if (base ~ /^(floor|ceil)/ || base == "ties-right") base = "current"
  lo = match(rule, /floor([0-9]+)/, g) ? g[1] + 0 : 0
  top = match(rule, /ceil([0-9]+)/, g) ? g[1] + 0 : 1e9
  for (i = 1; i <= nr[k]; i++) {
    r = rung[k, i]; s = 0; m = 0
    for (j = 1; j <= nv[k, i]; j++) {
      if (!vu[k, i, j]) continue
      if (base == "current") s++
      # relzNN: a valley counts when its z is at least NN% of the strongest usable valley anywhere
      # on this round's ladder, so the scale of z comes from the round and not from a constant
      else if (base ~ /^relz[0-9]+$/) { if (vz[k, i, j] >= substr(base, 5) / 100 * ladder_maxz(k)) s++ }
      else if (base == "z10") { if (vz[k, i, j] >= 10) s++ }
      else if (base == "z20") { if (vz[k, i, j] >= 20) s++ }
      else if (base == "sumz") s += vz[k, i, j]
      else if (base == "maxz") { if (vz[k, i, j] > m) m = vz[k, i, j]; s = m }
      else if (base == "persist") {
        # A usable valley counts when an adjacent rung has a usable valley within 0.03 of its share
        near = 0
        for (w = i - 1; w <= i + 1; w += 2) {
          if (w < 1 || w > nr[k]) continue
          for (u = 1; u <= nv[k, w]; u++)
            if (vu[k, w, u] && vs[k, w, u] - vs[k, i, j] <= 0.03 && vs[k, i, j] - vs[k, w, u] <= 0.03) near = 1
        }
        s += near
      }
    }
    if (r < lo || r > top) s = 0
    sc[i] = s
  }
}
function ladder_maxz(k,   i, j, m) {
  m = 0
  for (i = 1; i <= nr[k]; i++)
    for (j = 1; j <= nv[k, i]; j++) if (vu[k, i, j] && vz[k, i, j] > m) m = vz[k, i, j]
  return m
}
# The finest usable valley of rung i whose z is at least frac of the ladder's strongest; -1 if none
function finest(k, i, frac,   j, s, m) {
  s = -1; m = frac * ladder_maxz(k)
  for (j = 1; j <= nv[k, i]; j++)
    if (vu[k, i, j] && vz[k, i, j] >= m && (s < 0 || vs[k, i, j] < s)) s = vs[k, i, j]
  return s
}
# stabNN: the first rung from which the finest valley holding at least half the ladder's strongest
# z stays put when the axes double -- its share at the next rung within NN% of its own -- so the
# ladder stops where adding room no longer moves the level.  No such rung: the current rule
function stable(k, tol,   i, a, b) {
  for (i = 1; i < nr[k]; i++) {
    a = finest(k, i, 0.5); b = finest(k, i + 1, 0.5)
    if (a > 0 && b > 0 && (a - b <= tol * a) && (b - a <= tol * a)) return rung[k, i]
  }
  scores(k, "current")
  return first_plateau(k, sc)
}
function pick(k, rule) {
  if (rule ~ /^stab[0-9]+$/) return stable(k, substr(rule, 5) / 100)
  scores(k, rule)
  return rule ~ /ties-right/ ? last_max(k, sc) : first_plateau(k, sc)
}
function median(a, n,   b, i) {
  if (n == 0) return "-"
  for (i = 1; i <= n; i++) b[i] = a[i]
  asort(b)
  return n % 2 ? b[(n + 1) / 2] : (b[n / 2] + b[n / 2 + 1]) / 2
}

FNR == 1 { file++; next }
file == 1 {
  k = $1 SUBSEP $2 SUBSEP $3
  if (!(k in nr)) { order[++nk] = k }
  nr[k]++; rung[k, nr[k]] = $4; idx[k, $4] = nr[k]
  next
}
file == 2 {
  k = $1 SUBSEP $2 SUBSEP $3; i = idx[k, $4]
  nv[k, i]++; j = nv[k, i]; vs[k, i, j] = $5; vz[k, i, j] = ($6 == "inf" ? 1e9 : $6); vu[k, i, j] = $7
  next
}
file == 3 {
  k = $1 SUBSEP $2 SUBSEP $3
  logged[k] = $4
  # outcome of the round that actually searched at this rung on this corpus
  on[$1, $4]++; oh[$1, $4, on[$1, $4]] = $8; oc[$1, $4, on[$1, $4]] = $9
  next
}
file == 4 { next }

END {
  if (RULES == "") RULES = "current z10 z20 sumz maxz ties-right floor8 ceil64 persist"
  nrules = split(RULES, rules, " ")
  for (q = 1; q <= nrules; q++) {
    rule = rules[q]
    delete cnt; delete le4; delete hi; delete chg; delete good; delete r1; delete mv; delete nmv; delete picks; delete np
    for (t = 1; t <= nk; t++) {
      k = order[t]; split(k, f, SUBSEP); c = f[1]; seed = f[2]; rd = f[3]
      p = pick(k, rule)
      printf "PICK\t%s\t%s\t%s\t%s\t%s\n", rule, c, seed, rd, p
      cnt[c]++; np[c]++; picks[c, np[c]] = p
      if (p > 0 && p <= 4) le4[c]++
      if (p >= 128) hi[c]++
      if (p != logged[k]) chg[c]++
      if (p >= 16 && p <= 32) good[c]++
      if (rd == 1) r1[c] = r1[c] " s" seed ":" p
      prev = pk[rule, c, seed, rd - 1]
      pk[rule, c, seed, rd] = p
      if (rd > 1 && p > 0 && prev > 0) { nmv[c]++; d = log(p / prev) / log(2); mv[c, nmv[c]] = d < 0 ? -d : d }
      # estimated outcome: auto rounds on this corpus that searched at this rung
      for (e = 1; e <= on[c, p]; e++) { nh[rule, c]++; eh[rule, c, nh[rule, c]] = oh[c, p, e]; ec[rule, c, nh[rule, c]] = oc[c, p, e] }
      if (on[c, p] == 0) unk[rule, c]++
    }
    for (c in cnt) {
      delete pp; delete mm; delete hh; delete cc
      for (i = 1; i <= np[c]; i++) pp[i] = picks[c, i]
      for (i = 1; i <= nmv[c]; i++) mm[i] = mv[c, i]
      for (i = 1; i <= nh[rule, c]; i++) { hh[i] = eh[rule, c, i]; cc[i] = ec[rule, c, i] }
      asort(pp)
      printf "SUMMARY\t%s\t%s\tle4=%d/%d\t16to32=%d\tmedian=%s\trange=%s-%s\thigh128=%d\tchanged=%d\tmovement=%.2f\tround1=%s\test_hom=%s\test_comp=%s\tunmatched=%d\n",
        rule, c, le4[c] + 0, cnt[c], good[c] + 0, median(pp, np[c]), pp[1], pp[np[c]], hi[c] + 0, chg[c] + 0,
        median(mm, nmv[c]), r1[c], median(hh, nh[rule, c]), median(cc, nh[rule, c]), unk[rule, c] + 0
    }
  }
}
