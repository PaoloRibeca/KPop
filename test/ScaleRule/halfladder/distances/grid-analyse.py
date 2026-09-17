#!/usr/bin/env python3
"""grid-analyse -- every distance and metric setting, reduced to what the page draws

    grid-analyse <grid_directory> <output.json>

For each setting directory <metric>.<distance> and each embedding in it: the curves per rung (the
labels' silhouette, each rung's own partition, the shuffled baseline and the excess over it, the
valley strength, the misordered share), the half-powers ladder's pick, the upper bound -- the
largest downward step in the log of the strongest usable valley's z, permutation-tested -- and
where the labels' and the excess curves peak.  Then a scoreboard per setting.  A scratch analysis;
the page reads the JSON.
"""
import json, math, os, random, sys

GRID, OUT = sys.argv[1], sys.argv[2]
TAGS = ['vp1.T_rand', 'vp1.T_strat', 'vp1.T_core', 'rdrp.T_rand', 'rdrp.T_strat', 'vp2.T_rand', 'vp1s300', 'rdrps300']
SMALL = {'vp1s300', 'rdrps300'}
SETTINGS = ['powers.euclidean', 'flat.euclidean', 'powers.angle', 'flat.angle', 'powers.manhattan', 'flat.manhattan']
rng = random.Random(20260916)


def tstat(x, i):
    a, b = x[:i], x[i:]
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    ss = sum((v - ma) ** 2 for v in a) + sum((v - mb) ** 2 for v in b)
    dof = len(x) - 2
    sp = math.sqrt(ss / dof) if dof > 0 and ss > 0 else 1e-9
    return (ma - mb) / (sp * math.sqrt(1 / len(a) + 1 / len(b)))


def best_split(x):
    best = (-1e9, None)
    for i in range(2, len(x) - 1):
        t = tstat(x, i)
        if t > best[0]:
            best = (t, i)
    return best


def spearman(a, b):
    if len(a) < 4:
        return None
    def ranks(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
        for k, i in enumerate(o): r[i] = k
        return r
    ra, rb = ranks(a), ranks(b); n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
    return None if den == 0 else num / den


def read_embedding(d, tag):
    tsv, out, cand = (os.path.join(d, tag + ext) for ext in ('.tsv', '.out', '.cand'))
    if not (os.path.exists(tsv) and os.path.getsize(tsv) > 0):
        return None
    f = lambda x: None if x == '-' else float(x)
    rows = {}
    for line in open(tsv):
        if line.startswith('#') or line.startswith('tag'):
            continue
        p = line.rstrip('\n').split('\t')
        rows[int(p[1])] = {'d': int(p[1]), 'lab': f(p[2]), 'ownK': None if p[3] == '-' else int(p[3]), 'own': f(p[4]),
                           'null': f(p[5]), 'nullSd': f(p[6]), 'excess': f(p[8]), 'z': None, 'mis': None, 'score': None}
    pick = None
    for line in open(out):
        if line.startswith('  d='):
            dd = int(line.split('d=')[1].split()[0])
            if dd in rows:
                rows[dd]['mis'] = float(line.split('misordered ')[1].split('%')[0])
                rows[dd]['score'] = int(line.split('score ')[1].split()[0])
        if 'RULE half-powers ladder, most valleys, first plateau' in line and 'pick none' not in line:
            q = line.split('pick')[1].split('(misordered')
            pick = {'d': int(q[0]), 'mis': float(q[1].split('%')[0])}
    for line in open(cand):
        p = line.rstrip('\n').split('\t')
        if len(p) < 9 or p[0] != 'CAND' or p[8] != 'kept' or float(p[4]) > 0.5:
            continue
        dd = int(p[2]); zv = 999.0 if p[5] in ('inf', 'infinity') else float(p[5])
        if dd in rows:
            rows[dd]['z'] = max(rows[dd]['z'] or 0.0, zv)
    rows = [rows[k] for k in sorted(rows)]
    ds = [r['d'] for r in rows]
    # The upper bound, from the first rung with a usable valley on
    first = next((i for i, r in enumerate(rows) if r['z'] is not None), None)
    bound, p, t_obs = None, None, None
    if first is not None and len(rows) - first >= 5:
        lz = [math.log(max(r['z'], 3.0)) if r['z'] is not None else math.log(3.0) for r in rows[first:]]
        t_obs, i = best_split(lz)
        null = [best_split(rng.sample(lz, len(lz)))[0] for _ in range(2000)]
        p = (1 + sum(1 for v in null if v >= t_obs)) / (1 + len(null))
        if i is not None and p < 0.05:
            bound = ds[first + i - 1]
    lab_peak = max(rows, key=lambda r: r['lab'])['d']
    upto = [r for r in rows if bound is None or r['d'] <= bound]
    ex = [r for r in upto if r['excess'] is not None]
    ex_peak = max(ex, key=lambda r: r['excess'])['d'] if ex else None
    own = [r for r in upto if r['own'] is not None]
    mis_rows = [r for r in rows if r['mis'] is not None]
    mis_best = min(mis_rows, key=lambda r: r['mis']) if mis_rows else None
    idx = {dd: k for k, dd in enumerate(ds)}
    return {
        'tag': tag, 'small': tag in SMALL, 'rows': rows, 'pick': pick, 'bound': bound, 'p': p,
        'labPeak': lab_peak, 'excessPeak': ex_peak,
        'excessNearLabels': None if ex_peak is None else abs(idx[ex_peak] - idx[lab_peak]) <= 1,
        'labInside': None if bound is None else lab_peak <= bound,
        'rhoExcess': spearman([r['excess'] for r in ex], [r['lab'] for r in ex]),
        'rhoOwn': spearman([r['own'] for r in own], [r['lab'] for r in own]),
        'misBest': None if mis_best is None else {'d': mis_best['d'], 'mis': mis_best['mis']},
    }


def mean(v):
    v = [x for x in v if x is not None]
    return None if not v else sum(v) / len(v)


result = {'settings': [], 'order': TAGS}
for s in SETTINGS:
    d = os.path.join(GRID, s)
    embs = [e for e in (read_embedding(d, t) for t in TAGS) if e is not None] if os.path.isdir(d) else []
    full = [e for e in embs if not e['small']]
    small = [e for e in embs if e['small']]
    bounded = [e for e in embs if e['bound'] is not None]
    result['settings'].append({
        'name': s, 'metric': s.split('.')[0], 'distance': s.split('.')[1], 'n': len(embs),
        'score': {
            'bounded': len(bounded),
            'labInside': sum(1 for e in bounded if e['labInside']),
            'misPickFull': mean([e['pick']['mis'] for e in full if e['pick']]),
            'misPickSmall': mean([e['pick']['mis'] for e in small if e['pick']]),
            'excessNear': sum(1 for e in embs if e['excessNearLabels']),
            'rhoExcess': mean([e['rhoExcess'] for e in embs]),
            'rhoOwn': mean([e['rhoOwn'] for e in embs]),
        },
        'embeddings': embs,
    })
json.dump(result, open(OUT, 'w'), separators=(',', ':'))
for s in result['settings']:
    sc = s['score']
    fmt = lambda v, k=2: '—' if v is None else f'{v:.{k}f}'
    print(f"{s['name']:18s} n={s['n']} bounded={sc['bounded']} labInside={sc['labInside']} "
          f"misPick full={fmt(sc['misPickFull'])}% small={fmt(sc['misPickSmall'])}% "
          f"excessNear={sc['excessNear']} rhoExcess={fmt(sc['rhoExcess'])} rhoOwn={fmt(sc['rhoOwn'])}")
    for e in s['embeddings']:
        print(f"    {e['tag']:12s} bound {str(e['bound']):>4s} labPeak {e['labPeak']:>3d} excessPeak {str(e['excessPeak']):>4s} "
              f"pick {e['pick']['d'] if e['pick'] else '—':>3} ({fmt(e['pick']['mis']) if e['pick'] else '—'}%) "
              f"rhoEx {fmt(e['rhoExcess'])} rhoOwn {fmt(e['rhoOwn'])}")
