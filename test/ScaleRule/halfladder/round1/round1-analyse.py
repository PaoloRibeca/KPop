#!/usr/bin/env python3
"""round1-analyse -- seed agreement and label scores of the first round's search, rung by rung

    round1-analyse <round1_directory> <controls_directory|-> <labels.cdc> <output.json> [<corpus>]

For every rung: the partition each seed's first-round search returned (from its stdout), the adjusted
Rand index between every pair of seeds, and each seed's homogeneity, completeness and V-measure
against the labels, which judge and take no part.  Then a golden-section search on the median seed
agreement over the rungs, as a label-blind procedure would run it, recording which rungs it probed
and where it stopped.  Also checks that the first round at the rungs the controls ran matches them.
"""
import json, lzma, math, os, sys
from collections import Counter

R1, CTRL, LABELS, OUT = sys.argv[1:5]
CORPUS = sys.argv[5] if len(sys.argv) > 5 else 'rdrp'
RUNGS = [1, 2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64]
SEEDS = [17, 18, 19]


def first_round(text):
    part, inside = {}, False
    for line in text.splitlines():
        if line.startswith('=== Clustering'):
            if inside:
                break
            inside = True
            continue
        if not inside or line.startswith('#') or line.startswith('name\t'):
            continue
        p = line.split('\t')
        if len(p) == 3:
            part[p[0]] = p[1]
    return part


def contingency(a, b, names):
    return Counter((a[n], b[n]) for n in names), Counter(a[n] for n in names), Counter(b[n] for n in names)


def ari(a, b):
    names = [n for n in a if n in b]
    cell, ca, cb = contingency(a, b, names)
    c2 = lambda x: x * (x - 1) / 2
    index = sum(c2(v) for v in cell.values())
    sa, sb, tot = sum(c2(v) for v in ca.values()), sum(c2(v) for v in cb.values()), c2(len(names))
    exp = sa * sb / tot if tot else 0
    mx = (sa + sb) / 2
    return (index - exp) / (mx - exp) if mx != exp else 1.0


def scores(part, lab):
    names = [n for n in part if n in lab]
    N = len(names)
    cell, cp, cl = contingency(part, lab, names)
    H = lambda c: -sum(v / N * math.log(v / N) for v in c.values())
    hl, hp = H(cl), H(cp)
    h_l_given_p = -sum(v / N * math.log(v / cp[k[0]]) for k, v in cell.items())
    h_p_given_l = -sum(v / N * math.log(v / cl[k[1]]) for k, v in cell.items())
    hom = 1 - h_l_given_p / hl if hl else 1.0
    comp = 1 - h_p_given_l / hp if hp else 1.0
    v = 2 * hom * comp / (hom + comp) if hom + comp else 0.0
    return hom, comp, v, ari(part, {n: lab[n] for n in names})


lab = {}
for line in open(LABELS):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 2:
        lab[p[0]] = p[1]
status = {}
if os.path.exists(os.path.join(R1, 'runs.log')):
    for line in open(os.path.join(R1, 'runs.log')):
        p = line.rstrip('\n').split('\t')
        status[p[0]] = p[1]

rungs = []
for d in RUNGS:
    parts, row = {}, {'d': d, 'seeds': {}}
    for s in SEEDS:
        tag = f'{CORPUS}.r1d{d}.s{s}'
        path = os.path.join(R1, tag + '.stdout')
        st = status.get(tag)
        part = first_round(open(path).read()) if os.path.exists(path) else {}
        if st == 'exit=0' and part:
            parts[s] = part
            hom, comp, v, la = scores(part, lab)
            row['seeds'][s] = {'k': len(set(part.values())), 'hom': hom, 'comp': comp, 'v': v, 'ariLabels': la}
        else:
            row['seeds'][s] = {'status': st or 'not run'}
    pairs = [(a, b) for i, a in enumerate(SEEDS) for b in SEEDS[i + 1:] if a in parts and b in parts]
    row['agreement'] = [ari(parts[a], parts[b]) for a, b in pairs]
    med = lambda v: None if not v else sorted(v)[len(v) // 2] if len(v) % 2 else sum(sorted(v)[len(v) // 2 - 1:len(v) // 2 + 1]) / 2
    row['agreementMedian'] = med(row['agreement'])
    got = [x for x in row['seeds'].values() if 'v' in x]
    for key in ('hom', 'comp', 'v', 'ariLabels', 'k'):
        row[key + 'Median'] = med([x[key] for x in got])
    row['finished'] = len(got)
    rungs.append(row)

# Golden-section search over the rung index, maximising median seed agreement; a rung with no
# agreement (fewer than two partitions) scores -inf
f = {i: (r['agreementMedian'] if r['agreementMedian'] is not None else -math.inf) for i, r in enumerate(rungs)}
lo, hi, probes = 0, len(rungs) - 1, []
def probe(i):
    if i not in probes:
        probes.append(i)
    return f[i]
phi = (math.sqrt(5) - 1) / 2
c = round(hi - phi * (hi - lo)); dd = round(lo + phi * (hi - lo))
while hi - lo > 2:
    if c == dd:
        dd = c + 1
    if probe(c) >= probe(dd):
        hi = dd
    else:
        lo = c
    c = round(hi - phi * (hi - lo)); dd = round(lo + phi * (hi - lo))
best_i = max(range(lo, hi + 1), key=probe)
golden = {'probes': [rungs[i]['d'] for i in probes], 'pick': rungs[best_i]['d']}
judged = [r for r in rungs if r['vMedian'] is not None]
best_v = max(judged, key=lambda r: r['vMedian'])['d'] if judged else None
best_agree = max((r for r in rungs if r['agreementMedian'] is not None), key=lambda r: r['agreementMedian'], default=None)

# The controls ran 8, 16, 32 and 64 axes with the same settings: their first rounds must match
checks = []
for d in (8, 16, 32, 64):
    for s in SEEDS:
        cp = os.path.join(CTRL, f'{CORPUS}.fixed{d}.s{s}.stdout.xz')
        rp = os.path.join(R1, f'{CORPUS}.r1d{d}.s{s}.stdout')
        if os.path.exists(cp) and os.path.exists(rp):
            a, b = first_round(lzma.open(cp, 'rt').read()), first_round(open(rp).read())
            if a and b:
                checks.append({'d': d, 'seed': s, 'identical': a == b})

json.dump({'rungs': rungs, 'golden': golden, 'bestByLabels': best_v,
           'bestByAgreement': None if best_agree is None else best_agree['d'], 'controlChecks': checks},
          open(OUT, 'w'), indent=1)
fmt = lambda v: '   —  ' if v is None else f'{v:6.3f}'
print(' axes  done  agreement (pairs)            k    homog  compl  V      ARI-labels')
for r in rungs:
    ag = ' '.join(f'{x:.3f}' for x in r['agreement'])
    print(f"{r['d']:5d}  {r['finished']}/3  {fmt(r['agreementMedian'])} ({ag:20s})  {fmt(r['kMedian'])} {fmt(r['homMedian'])} {fmt(r['compMedian'])} {fmt(r['vMedian'])} {fmt(r['ariLabelsMedian'])}")
print('golden-section probes', golden['probes'], '-> picks', golden['pick'], '| best by agreement', best_agree['d'] if best_agree else None, '| best by labels (V)', best_v)
print('first rounds identical to the controls:', sum(1 for c in checks if c['identical']), 'of', len(checks))
