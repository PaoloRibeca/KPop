# usage: gawk -v tag=TAG -f pairpr.awk LABELS.cdc RUN.stdout
# pair precision = same-cluster pairs that are same-class / same-cluster pairs (1 for a refinement of the classes)
# pair recall    = same-class pairs that are same-cluster / same-class pairs
# hbound = H(K)/H(C), an upper bound on homogeneity given the cluster sizes
# big = largest class; bigTop = share of it in its largest cluster; bigK = clusters holding >= 5% of it
function flush(   s, c, k, key, N, Hc, Hk, pk, pc, pb, big, bigN, top, nb, p) {
  if (round == 0) return
  N = 0; delete cc; delete kc; delete ck
  for (s in rep) if (s in lab) { N++; cc[lab[s]]++; kc[rep[s]]++; ck[lab[s] SUBSEP rep[s]]++ }
  Hc = 0; for (c in cc) Hc -= cc[c]/N * log(cc[c]/N)
  Hk = 0; for (k in kc) Hk -= kc[k]/N * log(kc[k]/N)
  pk = 0; for (k in kc) pk += kc[k]*(kc[k]-1)/2
  pc = 0; bigN = 0; for (c in cc) { pc += cc[c]*(cc[c]-1)/2; if (cc[c] > bigN) { bigN = cc[c]; big = c } }
  pb = 0; top = 0; nb = 0
  for (key in ck) {
    pb += ck[key]*(ck[key]-1)/2
    split(key, p, SUBSEP)
    if (p[1] == big) { if (ck[key] > top) top = ck[key]; if (ck[key] >= 0.05*bigN) nb++ }
  }
  printf "%s\tr%d\tD=%s\tk=%d\tpairPrec=%.3f\tpairRec=%.3f\tHc=%.3f\tHk=%.3f\thbound=%.3f\tbig=%s(%d)\tbigTop=%.3f\tbigK=%d\n", tag, round, D, length(kc), pb/pk, pb/pc, Hc, Hk, Hk/Hc, big, bigN, top/bigN, nb
}
FNR == NR { lab[$1] = $2; next }
/^=== Clustering/ { flush(); round++; delete rep; match($0, /D=([0-9]+)/, m); D = m[1]; next }
/^# n=/ || /^name\trepresentative/ { next }
{ rep[$1] = $2 }
END { flush() }
