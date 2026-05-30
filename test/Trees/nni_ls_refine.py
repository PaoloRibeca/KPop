#!/usr/bin/env python3
"""NNI local search scored by least-squares additive (tree-metric) fit.

This is the operational form of "project the distance matrix onto the
closest tree": starting from a given tree, try Nearest-Neighbour
Interchange moves and accept those that lower the residual
||D - patristic_OLS(T)||^2 (the closest-tree objective).  At each step
we also report the topology Jaccard against a reference tree, to see
whether moving toward the closest tree metric moves toward or away from
the true topology.

Usage:
  nni_ls_refine.py <KPopDMatrix.txt> <start.nwk> <ref.nwk> [max_sweeps]
"""

import sys
import numpy as np


# ---------- IO ----------

def load_dmatrix(path):
    names, rows = [], []
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")[1:]
        for line in f:
            p = line.rstrip("\n").split("\t")
            names.append(p[0]); rows.append([float(v) for v in p[1:]])
    return header, np.array(rows)


def parse_newick(path):
    s = open(path).read().strip()
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
        return None
    def clade():
        if s[pos[0]] == "(":
            pos[0] += 1; kids = []
            while True:
                kids.append(clade())
                if s[pos[0]] == ",": pos[0] += 1
                elif s[pos[0]] == ")": pos[0] += 1; break
            me = nid()
            for k in kids: edges.append((me, k))
            name(); length()
            return me
        else:
            nm = name(); length(); me = nid(); leaves[me] = nm
            return me
    clade()
    return edges, leaves


# ---------- tree as undirected adjacency ----------

def build_adj(edges, leaves):
    adj = {}
    for u, v in edges:
        adj.setdefault(u, set()).add(v); adj.setdefault(v, set()).add(u)
    # contract degree-2 nodes (e.g. an artificial rooted-binary root)
    changed = True
    while changed:
        changed = False
        for x in list(adj):
            if x not in leaves and len(adj[x]) == 2:
                a, b = tuple(adj[x])
                adj[a].discard(x); adj[b].discard(x)
                adj[a].add(b); adj[b].add(a)
                del adj[x]; changed = True; break
    return adj


def leaf_splits(adj, leaves):
    """Set of internal-edge bipartitions, each as frozenset of the
    lexicographically-smaller side's leaf names."""
    all_leaves = frozenset(leaves.values())
    splits = set()
    for u in adj:
        for v in adj[u]:
            if u < v and u not in leaves and v not in leaves:
                side = leaves_on_side(adj, leaves, u, v)
                names = frozenset(leaves[l] for l in side)
                other = all_leaves - names
                if 1 < len(names) < len(all_leaves) - 1:
                    splits.add(min(names, other, key=lambda s: (len(s), sorted(s))))
    return splits


def leaves_on_side(adj, leaves, u, v):
    """Leaves reachable from u without crossing edge u-v."""
    seen = {v, u}; stack = [u]; out = []
    if u in leaves: out.append(u)
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                if y in leaves: out.append(y)
                stack.append(y)
    return out


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a or b) else 1.0


# ---------- least-squares additive fit ----------

_PAIR_SUBSAMPLE = None   # set to an int to score on a fixed random pair subset

def ls_rss(adj, leaves, name_idx, D):
    """Residual sum of squares of the OLS additive fit of D to this
    topology (optimal branch lengths)."""
    edge_list = []
    for u in adj:
        for v in adj[u]:
            if u < v:
                edge_list.append((u, v))
    edge_id = {e: k for k, e in enumerate(edge_list)}
    E = len(edge_list)
    lf = sorted(leaves)                      # leaf node ids
    # per-leaf edge path to a fixed root leaf, as a set of edge ids
    root = lf[0]
    paths = {}
    for lid in lf:
        # BFS from root, record parent edges, then walk back from lid
        pass
    # build per-node parent edge via BFS from root
    parent_edge = {root: None}; order = [root]; seen = {root}
    while order:
        x = order.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                e = (min(x, y), max(x, y))
                parent_edge[y] = (e, x)
                order.append(y)
    def edges_to_root(lid):
        s = set(); x = lid
        while parent_edge[x] is not None:
            e, p = parent_edge[x]; s.add(edge_id[e]); x = p
        return s
    root_path = {lid: edges_to_root(lid) for lid in lf}
    pairs = [(i, j) for i in range(len(lf)) for j in range(i + 1, len(lf))]
    if _PAIR_SUBSAMPLE is not None and len(pairs) > _PAIR_SUBSAMPLE:
        rng = np.random.default_rng(12345)        # fixed -> comparable across candidates
        sel = rng.choice(len(pairs), size=_PAIR_SUBSAMPLE, replace=False)
        pairs = [pairs[k] for k in sel]
    M = np.zeros((len(pairs), E))
    d = np.zeros(len(pairs))
    for r, (i, j) in enumerate(pairs):
        li, lj = lf[i], lf[j]
        for eidx in root_path[li] ^ root_path[lj]:
            M[r, eidx] = 1.0
        d[r] = D[name_idx[leaves[li]], name_idx[leaves[lj]]]
    b, res, rank, sv = np.linalg.lstsq(M, d, rcond=None)
    fit = M @ b
    return float(np.sum((fit - d) ** 2))


# ---------- NNI ----------

def internal_edges(adj, leaves):
    out = []
    for u in adj:
        for v in adj[u]:
            if u < v and u not in leaves and v not in leaves:
                out.append((u, v))
    return out


def apply_nni(adj, u, v, which):
    """In place: swap a u-side subtree with a v-side subtree."""
    a = [x for x in adj[u] if x != v]
    c = [x for x in adj[v] if x != u]
    su = a[0]; sv = c[0] if which == 0 else c[1]
    adj[u].discard(su); adj[su].discard(u)
    adj[v].discard(sv); adj[sv].discard(v)
    adj[u].add(sv); adj[sv].add(u)
    adj[v].add(su); adj[su].add(v)


def copy_adj(adj):
    return {k: set(v) for k, v in adj.items()}


# ---------- main ----------

def main():
    global _PAIR_SUBSAMPLE
    dmx, start, ref = sys.argv[1], sys.argv[2], sys.argv[3]
    max_sweeps = int(sys.argv[4]) if len(sys.argv) > 4 else 30
    if len(sys.argv) > 5:
        _PAIR_SUBSAMPLE = int(sys.argv[5])
    names, D = load_dmatrix(dmx)
    name_idx = {nm: i for i, nm in enumerate(names)}
    e0, lv0 = parse_newick(start); adj = build_adj(e0, lv0)
    er, lvr = parse_newick(ref); radj = build_adj(er, lvr)
    ref_splits = leaf_splits(radj, lvr)

    def report(tag, adj, leaves):
        rss = ls_rss(adj, leaves, name_idx, D)
        jac = jaccard(leaf_splits(adj, leaves), ref_splits)
        print(f"  {tag:<22} RSS={rss:12.4f}  J(vs ref)={jac:.4f}")
        return rss, jac

    leaves = lv0
    print(f"=== NNI-LS refinement: {start}  (n={len(leaves)}) ===")
    rss, jac0 = report("start (sparse-NJ)", adj, leaves)
    for sweep in range(1, max_sweeps + 1):
        best = None
        for (u, v) in internal_edges(adj, leaves):
            for which in (0, 1):
                cand = copy_adj(adj)
                apply_nni(cand, u, v, which)
                r = ls_rss(cand, leaves, name_idx, D)
                if r < rss - 1e-9 and (best is None or r < best[0]):
                    best = (r, cand)
        if best is None:
            print(f"  (converged after {sweep-1} sweeps)")
            break
        rss, adj = best[0], best[1]
        jac = jaccard(leaf_splits(adj, leaves), ref_splits)
        print(f"  sweep {sweep:<2d}            RSS={rss:12.4f}  J(vs ref)={jac:.4f}", flush=True)
    rssf, jacf = report("final (LS-closest)", adj, leaves)
    print(f"  --> Jaccard vs ref: {jac0:.4f} (start) -> {jacf:.4f} (LS-refined); "
          f"{'UP' if jacf>jac0 else 'DOWN' if jacf<jac0 else 'SAME'}")


if __name__ == "__main__":
    main()
