#!/usr/bin/env python3
"""hold-analyse -- the loop held at the searched number of axes, against stage 1's arms

    hold-analyse <hold_directory> <stage1_runs_directory> <controls_directory> <labels_directory> <output.json> [<sample_directory>]

For every run of every arm -- today (half-height, --dimensions-inertia 0.5), auto (the doubling
ladder), and the loop held at a fixed number of axes -- round by round: the number of axes, the
clusters, the level, homogeneity and completeness against the labels, and the adjusted Rand index
between the round's partition and the previous round's.  Then per arm the stage-1 summary -- each
run's median over its last four rounds, the median over seeds and the spread -- and whether it
differs from today by more than twice the larger spread.  A scratch analysis; the page reads the
JSON.
"""
import json, lzma, math, os, re, sys
from collections import Counter

HOLD, RUNS, CTRL, LABELS, OUT = sys.argv[1:6]
SAMPLE = sys.argv[6] if len(sys.argv) > 6 else None
ARMS = {
    'rdrp': [('today', RUNS), ('auto', RUNS), ('fixed16', CTRL)],
    'vp2': [('today', RUNS), ('auto', RUNS), ('fixed11', HOLD)],
    'vp1': [('today', RUNS), ('auto', RUNS), ('fixed8', HOLD), ('fixed11', HOLD)],
}
if SAMPLE:
    for corpus in ARMS:
        ARMS[corpus] += [('bysize', SAMPLE), ('random', SAMPLE), ('bysizerandom', SAMPLE)]
LOGGED = {HOLD, SAMPLE} - {None}
SEEDS = [17, 18, 19]


def read(path):
    if os.path.exists(path):
        return open(path).read()
    if os.path.exists(path + '.xz'):
        return lzma.open(path + '.xz', 'rt').read()
    return None


def rounds(text):
    out, cur = [], None
    for line in text.splitlines():
        if line.startswith('=== Clustering'):
            m = re.search(r'D=(\d+)', line)
            cur = {'d': int(m.group(1)) if m else None, 'part': {}}
            out.append(cur)
            continue
        if cur is None:
            continue
        if line.startswith('# n='):
            m = re.search(r'n_clusters=(\d+)', line); cur['k'] = int(m.group(1)) if m else None
            m = re.search(r' fv=([0-9.eE-]+)', line); cur['fv'] = float(m.group(1)) if m else None
            continue
        if line.startswith('name\t'):
            continue
        p = line.split('\t')
        if len(p) == 3:
            cur['part'][p[0]] = p[1]
    return out


def ari(a, b):
    names = [n for n in a if n in b]
    cell, ca, cb = Counter((a[n], b[n]) for n in names), Counter(a[n] for n in names), Counter(b[n] for n in names)
    c2 = lambda x: x * (x - 1) / 2
    idx = sum(c2(v) for v in cell.values())
    sa, sb, tot = sum(c2(v) for v in ca.values()), sum(c2(v) for v in cb.values()), c2(len(names))
    exp = sa * sb / tot if tot else 0
    mx = (sa + sb) / 2
    return (idx - exp) / (mx - exp) if mx != exp else 1.0


def hom_comp(part, lab):
    names = [n for n in part if n in lab]
    N = len(names)
    cell, cp, cl = Counter((part[n], lab[n]) for n in names), Counter(part[n] for n in names), Counter(lab[n] for n in names)
    H = lambda c: -sum(v / N * math.log(v / N) for v in c.values())
    hl, hp = H(cl), H(cp)
    hlp = -sum(v / N * math.log(v / cp[k[0]]) for k, v in cell.items())
    hpl = -sum(v / N * math.log(v / cl[k[1]]) for k, v in cell.items())
    return (1 - hlp / hl if hl else 1.0), (1 - hpl / hp if hp else 1.0)


def median(v):
    v = sorted(x for x in v if x is not None)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


result = {'corpora': {}}
for corpus, arms in ARMS.items():
    lab = {}
    for line in open(os.path.join(LABELS, corpus + '.cdc')):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 2:
            lab[p[0]] = p[1]
    carms = []
    for arm, where in arms:
        runs = []
        done = set()
        log = os.path.join(where, 'runs.log')
        if where in LOGGED and os.path.exists(log):
            done = {l.split('\t')[0] for l in open(log) if '\texit=0\t' in l}
        for s in SEEDS:
            if where in LOGGED and f'{corpus}.{arm}.s{s}' not in done:
                continue
            text = read(os.path.join(where, f'{corpus}.{arm}.s{s}.stdout'))
            if text is None:
                continue
            rs = rounds(text)
            prev, rows = None, []
            for i, r in enumerate(rs):
                h, c = hom_comp(r['part'], lab)
                rows.append({'r': i + 1, 'd': r['d'], 'k': r.get('k'), 'fv': r.get('fv'), 'hom': h, 'comp': c,
                             'ariPrev': None if prev is None else ari(r['part'], prev)})
                prev = r['part']
            last = rows[-4:] if len(rows) >= 4 else rows
            runs.append({'seed': s, 'rounds': rows, 'hom': median([x['hom'] for x in last]),
                         'comp': median([x['comp'] for x in last]),
                         'stability': median([x['ariPrev'] for x in rows[1:]])})
        if not runs:
            continue
        spread = lambda key: max(r[key] for r in runs) - min(r[key] for r in runs)
        carms.append({'arm': arm, 'runs': runs, 'hom': median([r['hom'] for r in runs]), 'comp': median([r['comp'] for r in runs]),
                      'homSpread': spread('hom'), 'compSpread': spread('comp'),
                      'stability': median([r['stability'] for r in runs]), 'settled': sum(1 for r in runs if len(r['rounds']) < 8)})
    today = next((a for a in carms if a['arm'] == 'today'), None)
    for a in carms:
        if today and a is not today:
            for key in ('hom', 'comp'):
                bar = 2 * max(a[key + 'Spread'], today[key + 'Spread'])
                a[key + 'Delta'] = a[key] - today[key]
                a[key + 'Effect'] = abs(a[key] - today[key]) > bar
                a[key + 'Bar'] = bar
    result['corpora'][corpus] = carms
json.dump(result, open(OUT, 'w'), separators=(',', ':'))
for corpus, carms in result['corpora'].items():
    print(f'=== {corpus}')
    for a in carms:
        eff = '' if a['arm'] == 'today' else (
            f"  vs today hom {a['homDelta']:+.3f} ({'effect' if a['homEffect'] else 'no effect'}, bar {a['homBar']:.3f})"
            f"  comp {a['compDelta']:+.3f} ({'effect' if a['compEffect'] else 'no effect'}, bar {a['compBar']:.3f})")
        ks = ' | '.join(','.join(str(x['k']) for x in r['rounds']) for r in a['runs'])
        print(f"  {a['arm']:8s} runs {len(a['runs'])} hom {a['hom']:.3f} (±{a['homSpread']:.3f}) comp {a['comp']:.3f} (±{a['compSpread']:.3f})"
              f" round-to-round ARI {a['stability']:.3f} settled early {a['settled']}{eff}")
        print(f"           clusters per round: {ks}")
