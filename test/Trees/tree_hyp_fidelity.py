#!/usr/bin/env python3
"""Demonstration: hyperbolically embed a TREE and show its fidelity
holds at high curvature (the Test-C property), in direct contrast to a
raw (non-tree) distance matrix, which degrades as curvature rises.

Usage:
  tree_hyp_fidelity.py <newick> [<KPopDMatrix.txt for contrast>]

The Newick is parsed to a patristic (path-length) leaf-distance matrix
-- an exact tree metric -- then embedded into the hyperboloid by
Lorentzian classical MDS at a sweep of curvature scales c and target
dimensions d.  Fidelity is the Pearson correlation between the
recovered hyperbolic distances and the target tree distances.
"""

import sys
import numpy as np
from lorentzian_mds import lorentzian_mds, hyp_dist_matrix, fit_stats


# ---------- minimal Newick -> patristic ----------

def parse_newick(s):
    s = s.strip()
    for tag in ("[&U]", "[&R]"):
        if s.startswith(tag):
            s = s[len(tag):]
    s = s.strip().rstrip(";")
    pos = [0]; cnt = [0]; edges = []; leaves = {}

    def nid():
        i = cnt[0]; cnt[0] += 1; return i

    def name():
        st = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ":,()":
            pos[0] += 1
        return s[st:pos[0]]

    def length():
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1; st = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",()":
                pos[0] += 1
            try:
                return float(s[st:pos[0]])
            except ValueError:
                return 0.0
        return 0.0

    def clade():
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = []
            while True:
                kids.append(clade())
                if s[pos[0]] == ",":
                    pos[0] += 1
                elif s[pos[0]] == ")":
                    pos[0] += 1; break
            me = nid()
            for (cid, clen) in kids:
                edges.append((me, cid, clen))
            name()
            return (me, length())
        else:
            nm = name(); me = nid(); leaves[me] = nm
            return (me, length())

    clade()
    return edges, leaves, cnt[0]


def patristic_matrix(newick_path):
    s = open(newick_path).read()
    edges, leaves, N = parse_newick(s)
    adj = [[] for _ in range(N)]
    for u, v, w in edges:
        adj[u].append((v, w)); adj[v].append((u, w))
    leaf_ids = sorted(leaves)
    idx = {lid: i for i, lid in enumerate(leaf_ids)}
    n = len(leaf_ids)
    D = np.zeros((n, n))
    for lid in leaf_ids:                       # BFS over the tree from each leaf
        dist = {lid: 0.0}; stack = [lid]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + w; stack.append(v)
        for lid2 in leaf_ids:
            D[idx[lid], idx[lid2]] = dist[lid2]
    return D, [leaves[l] for l in leaf_ids]


def load_dmatrix(path):
    rows = []
    with open(path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append([float(v) for v in p[1:]])
    return np.array(rows)


def sweep(label, D):
    D = 0.5 * (D + D.T)
    D = D / D[np.triu_indices_from(D, 1)].mean()
    print(f"--- {label} (n={D.shape[0]}) ---")
    print(f"  {'c':>5} {'d':>4} {'neg':>5} {'offsheet':>9} {'Pearson':>9}")
    for c in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        for d in (5, 10, 20):
            Xr, w, diag = lorentzian_mds(c * D, d)
            Dr = hyp_dist_matrix(Xr)
            pear, _ = fit_stats(c * D, Dr)
            print(f"  {c:>5.1f} {d:>4d} {diag['n_neg_eig']:>5d} {diag['n_offsheet']:>9d} {pear:>9.4f}")
    print()


if __name__ == "__main__":
    Dtree, names = patristic_matrix(sys.argv[1])
    print(f"=== Hyperbolic embedding of a TREE (patristic of {sys.argv[1]}) ===")
    print("    fidelity should stay HIGH as curvature c rises (the tree property)\n")
    sweep("TREE patristic", Dtree)
    if len(sys.argv) > 2:
        Dca = load_dmatrix(sys.argv[2])
        print("=== Contrast: the raw (non-tree) distance matrix ===")
        print("    fidelity should PEAK low and FALL as c rises\n")
        sweep("raw distance matrix", Dca)
