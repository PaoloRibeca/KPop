# usage: gawk -v tag=TAG -f score.awk LABELS.cdc RUN.stdout
function flush(   c, k, key, Hc, Hk, Hck, Hkc, N, h, cpl, pl, pa, pc, na, nk, nkl) {
  if (round == 0) return
  N = 0; delete cc; delete kc; delete ck; delete ka
  for (s in rep) {
    ka[rep[s]]++
    if (s in lab) { N++; cc[lab[s]]++; kc[rep[s]]++; ck[lab[s] SUBSEP rep[s]]++ }
  }
  Hc = 0; for (c in cc) Hc -= cc[c]/N * log(cc[c]/N)
  Hk = 0; for (k in kc) Hk -= kc[k]/N * log(kc[k]/N)
  Hck = 0; Hkc = 0
  for (key in ck) {
    split(key, p, SUBSEP)
    Hck -= ck[key]/N * log(ck[key]/kc[p[2]])
    Hkc -= ck[key]/N * log(ck[key]/cc[p[1]])
  }
  h = (Hc == 0) ? 1 : 1 - Hck/Hc
  cpl = (Hk == 0) ? 1 : 1 - Hkc/Hk
  # within-cluster pair share over labelled, and over all sequences
  pl = 0; for (k in kc) pl += kc[k]*(kc[k]-1)
  pl /= N*(N-1)
  na = 0; pa = 0; nk = 0; for (k in ka) { na += ka[k]; pa += ka[k]*(ka[k]-1); nk++ }
  pa /= na*(na-1)
  pc = 0; for (c in cc) pc += cc[c]*(cc[c]-1)
  pc /= N*(N-1)
  nkl = length(kc)
  printf "%s\t%d\tD=%s\tn=%d\tlab=%d\tk_all=%d\tk_lab=%d\thom=%.4f\tcom=%.4f\tpairK_lab=%.4f\tpairK_all=%.4f\tpairC_lab=%.4f\ttypes=%d\tfp=%s\tfv=%s\n", tag, round, D, na, N, nk, nkl, h, cpl, pl, pa, pc, length(cc), fp, fv
}
FNR == NR { lab[$1] = $2; next }
/^=== Clustering/ { flush(); round++; delete rep; match($0, /D=([0-9]+)/, m); D = m[1]; fp = "-"; fv = "-"; next }
/^# n=/ { if (match($0, / fp=([0-9.e-]+)/, m)) fp = m[1]; if (match($0, / fv=([0-9.e-]+)/, m)) fv = m[1]; next }
/^name\trepresentative/ { next }
{ rep[$1] = $2 }
END { flush() }
