#!/usr/bin/env python3
"""Hyperbolic sparse-NJ scaling study.

Two complementary sweeps:

  (a) Subsampled real data.  Take varying-n subsets of the
      prot_k5_kt0_035 principal-coord embedding (n in {15, 25, 40, 57})
      and ask: what's the smallest K_QUERY/K factor at which each mode
      matches the dense reference on that subset?

  (b) Synthetic BM-along-tree.  Generate random trees with n leaves
      and Brownian-motion-along-branches embeddings, varying n in
      {50, 100, 200, 400}.  Add Gaussian per-leaf noise to make the
      embedding non-additive (matches real-data character).

Topology agreement is measured via the existing tree_rf.py script
(split-set Jaccard against a reference Newick).

Headline question: does the K_QUERY/K factor required for the
hyperbolic mode to match the dense reference stay bounded as n grows?
"""
import os
import random
import re
import subprocess
import sys
import tempfile
import time

import numpy as np


KPOP_DIR = "/home/ribeca/Development/Git/KPop.Claude"
HYPBENCH = f"{KPOP_DIR}/_build/default/test/HypBench.exe"
TREE_RF = f"{KPOP_DIR}/test/Trees/tree_rf.py"
PROT_K5_EMB = "/tmp/phase0_emb.txt"  # set up by an earlier step


# ---------- Random tree generator (with noisy embedding) ----------

class Node:
    __slots__ = ("name", "children")
    def __init__(self, name=None):
        self.name = name
        self.children = []


def gen_random_tree(n, rng):
    if n < 2:
        raise ValueError("Need at least 2 leaves")
    leaves = [Node(name=f"L{0}"), Node(name=f"L{1}")]
    leaves[0].children = [(leaves[1], rng.expovariate(1.0))]
    edges = [(leaves[0], 0)]
    next_leaf_id = 2
    while next_leaf_id < n:
        parent, ci = rng.choice(edges)
        child, blen = parent.children[ci]
        new_leaf = Node(name=f"L{next_leaf_id}")
        next_leaf_id += 1
        new_bl = rng.expovariate(1.0)
        u = rng.random()
        blen1, blen2 = blen * u, blen * (1 - u)
        M = Node()
        M.children = [(child, blen2), (new_leaf, new_bl)]
        parent.children[ci] = (M, blen1)
        edges.remove((parent, ci))
        edges.append((parent, ci))
        edges.append((M, 0))
        edges.append((M, 1))
    return leaves[0]


def collect_leaves_and_paths(root):
    leaves = {}
    counter = [0]
    def dfs(node, path):
        if not node.children:
            leaves[node.name] = list(path)
        else:
            for child, blen in node.children:
                bid = counter[0]
                counter[0] += 1
                path.append((blen, bid))
                dfs(child, path)
                path.pop()
    dfs(root, [])
    return leaves, counter[0]


def to_newick_string(root):
    """Convert our Node tree to Newick with a top-level unrooted shape:
    we treat the root's single subtree as the entire tree.  Add a trailing
    semicolon and rooted-balanced shape so tree_rf.py can read it as a
    standard rooted Newick."""
    def rec(node):
        if not node.children:
            return node.name
        inner = ",".join(f"{rec(c)}:{bl:.6f}" for c, bl in node.children)
        return f"({inner})"
    # Top-level: our root has one or more children.  If exactly one,
    # treat it as just "subtree;" without enclosing parens.
    if len(root.children) == 1:
        c, _ = root.children[0]
        return rec(c) + ";"
    return rec(root) + ";"


def patristic_matrix(leaf_paths, names):
    n = len(names)
    D = np.zeros((n, n))
    paths = [leaf_paths[nm] for nm in names]
    for i in range(n):
        bi = {bid: bl for bl, bid in paths[i]}
        for j in range(i + 1, n):
            bj = {bid: bl for bl, bid in paths[j]}
            d = 0.0
            for bid, bl in bi.items():
                if bid not in bj:
                    d += bl
            for bid, bl in bj.items():
                if bid not in bi:
                    d += bl
            D[i, j] = d
            D[j, i] = d
    return D


def mds_embed_noisy(leaf_paths, names, dim, noise_sigma, rng_np):
    """Classical MDS embedding of the tree's patristic matrix into R^dim,
    with additive Gaussian noise per leaf coord to inject non-additivity.

    This gives ||p_i - p_j||_2 approx patristic(i, j), which is what NJ
    expects -- so dense NJ can actually recover the tree (modulo MDS
    distortion + the added noise).  Empirically dim ~ n / 2 gives
    Jaccard ~ 0.9+ for dense vs truth on the unnoised version."""
    D = patristic_matrix(leaf_paths, names)
    n = len(names)
    D2 = D ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D2 @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    d_use = min(dim, n - 1)
    pos_eigvals = np.maximum(eigvals[:d_use], 0)
    P = eigvecs[:, :d_use] * np.sqrt(pos_eigvals)
    # Pad with zeros if dim > n - 1
    if dim > P.shape[1]:
        P = np.hstack([P, np.zeros((n, dim - P.shape[1]))])
    if noise_sigma > 0:
        P = P + rng_np.normal(0, noise_sigma, P.shape)
    return P


# ---------- Subsampled real data ----------

def load_real_embeddings(path):
    """Load KPopVectors.txt format: header line, then name + coords per row."""
    with open(path) as f:
        f.readline()  # header
        names = []
        rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            names.append(parts[0])
            rows.append([float(x) for x in parts[1:]])
    return names, np.array(rows)


def subsample(names, X, n, rng):
    """Pick n random indices."""
    idx = sorted(rng.sample(range(len(names)), n))
    return [names[i] for i in idx], X[idx]


# ---------- Calls to HypBench + tree_rf ----------

def write_embeddings(path, names, P):
    with open(path, "w") as f:
        for i, nm in enumerate(names):
            row = "\t".join(f"{v:.10g}" for v in P[i])
            f.write(f"{nm}\t{row}\n")


def run_hypbench(emb_path, mode, k_nn, k_query_factor, hyp_scale, distance="saitou-nei"):
    out = subprocess.run(
        [HYPBENCH, emb_path, mode, str(k_nn), str(k_query_factor),
         str(hyp_scale), distance],
        check=True, capture_output=True, text=True)
    return out.stdout.strip()


# Surviving sparse-NJ candidate engines, benchmarked against the dense
# reference.  (label, mode, distance, hyp_scale):
#   rp-forest + saitou-nei  = exact, quadratic-time best
#   rp-forest + hyperbolic  = sub-quadratic-memory proxy (radial scale ~0.08)
#   periodic-rebuild        = recommended exact engine (persistent FAISS)
CANDIDATES = [
    ("rpf-saitou", "rp-forest", "saitou-nei", 1.0),
    ("rpf-hyper", "rp-forest", "hyperbolic", 0.08),
    ("periodic", "periodic-rebuild", "saitou-nei", 1.0),
]


_jacc_re = re.compile(r"^\s*\S+\s+(\d+)\s+(\d+|—)\s+(\d+|—)\s+([\d.]+|—)\s+([\d.]+|—)\s*$",
                      re.MULTILINE)


def jaccard_via_tree_rf(ref_nwk_path, test_nwk_path):
    out = subprocess.run(
        ["python3", TREE_RF, ref_nwk_path, test_nwk_path],
        check=True, capture_output=True, text=True)
    # Parse last numeric field (Jacc) from the non-reference row
    for line in out.stdout.splitlines():
        m = re.match(r"\s*\S+\s+\d+\s+\d+\s+\d+\s+[\d.]+\s+([\d.]+)\s*$", line)
        if m:
            return float(m.group(1))
    return None


def write_nwk(path, nwk):
    with open(path, "w") as f:
        f.write(nwk + "\n")


# ---------- Sweeps ----------

def sweep_real_subsamples():
    print("=" * 72)
    print("(a) Subsampled real data: prot_k5_kt0_035 principal-coord embedding")
    print("=" * 72)
    rng = random.Random(42)
    names_all, X_all = load_real_embeddings(PROT_K5_EMB)
    K = 10
    factors = [1, 2, 3, 5]
    print(f"{'n':>4} {'mode':<14} {'factor':>7} "
          f"{'J(vs dense)':>12} {'time_s':>8}")
    for n in [15, 25, 40, 57]:
        names, X = subsample(names_all, X_all, n, rng)
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False, mode="w") as tf:
            emb_path = tf.name
        write_embeddings(emb_path, names, X)
        # Dense reference
        t0 = time.time()
        dense_nwk = run_hypbench(emb_path, "dense", K, 1, 1.0)
        dt = time.time() - t0
        dense_nwk_path = emb_path + ".dense.nwk"
        write_nwk(dense_nwk_path, dense_nwk)
        print(f"{n:>4} {'dense':<14} {'--':>7} {1.000:>12.3f} {dt:>8.2f}")
        for label, mode, dist, hs in CANDIDATES:
            for f in factors:
                t0 = time.time()
                cand_nwk = run_hypbench(emb_path, mode, K, f, hs, dist)
                ct = time.time() - t0
                cand_path = emb_path + f".{label}.f{f}.nwk"
                write_nwk(cand_path, cand_nwk)
                j = jaccard_via_tree_rf(dense_nwk_path, cand_path)
                print(f"{n:>4} {label:<14} {f:>7} {j:>12.3f} {ct:>8.2f}")
                os.unlink(cand_path)
        os.unlink(dense_nwk_path)
        os.unlink(emb_path)
        print()


def sweep_synthetic():
    print("=" * 72)
    print("(b) Synthetic BM-along-tree with per-leaf Gaussian noise")
    print("=" * 72)
    rng = random.Random(43)
    rng_np = np.random.default_rng(43)
    K = 10
    factors = [1, 2, 3, 5]
    dim = 50
    noise_sigma = 0.5  # injects non-additivity
    print(f"{'n':>5} {'mode':<14} {'factor':>7} "
          f"{'J(vs dense)':>12} {'J(vs truth)':>12} {'time_s':>8}")
    for n in [50, 100, 200, 400]:
        tree = gen_random_tree(n, rng)
        leaf_paths, num_branches = collect_leaves_and_paths(tree)
        leaf_names = sorted(leaf_paths.keys(), key=lambda s: int(s[1:]))
        P = mds_embed_noisy(leaf_paths, leaf_names, dim, noise_sigma, rng_np)
        true_nwk = to_newick_string(tree)
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False, mode="w") as tf:
            emb_path = tf.name
        write_embeddings(emb_path, leaf_names, P)
        true_path = emb_path + ".true.nwk"
        write_nwk(true_path, true_nwk)
        t0 = time.time()
        dense_nwk = run_hypbench(emb_path, "dense", K, 1, 1.0)
        dt = time.time() - t0
        dense_path = emb_path + ".dense.nwk"
        write_nwk(dense_path, dense_nwk)
        j_truth = jaccard_via_tree_rf(true_path, dense_path)
        print(f"{n:>5} {'dense':<14} {'--':>7} {1.000:>12.3f} {j_truth:>12.3f} {dt:>8.2f}")
        for label, mode, dist, hs in CANDIDATES:
            for f in factors:
                t0 = time.time()
                cand_nwk = run_hypbench(emb_path, mode, K, f, hs, dist)
                ct = time.time() - t0
                cand_path = emb_path + f".{label}.f{f}.nwk"
                write_nwk(cand_path, cand_nwk)
                j_dense = jaccard_via_tree_rf(dense_path, cand_path)
                j_truth_cand = jaccard_via_tree_rf(true_path, cand_path)
                print(f"{n:>5} {label:<14} {f:>7} {j_dense:>12.3f} {j_truth_cand:>12.3f} {ct:>8.2f}")
                os.unlink(cand_path)
        os.unlink(true_path)
        os.unlink(dense_path)
        os.unlink(emb_path)
        print()


def main():
    if not os.path.exists(HYPBENCH):
        print(f"FAIL: {HYPBENCH} not found", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(PROT_K5_EMB):
        # Re-create it from the OCaml binary
        kbin = f"{KPOP_DIR}/build/KPopTwistDB"
        prot = f"{KPOP_DIR}/test/Trees/Tree.prot_k5_kt0_035"
        subprocess.run(
            [kbin, "-i", "T", prot, "-i", "t", prot, "-m", "powers(1,1,1)",
             "-e", PROT_K5_EMB.replace(".KPopVectors.txt", "").replace(".txt", "")],
            check=True, capture_output=True)
        if not os.path.exists(PROT_K5_EMB):
            # KPopVectors writes to <prefix>.KPopVectors.txt by default
            target = PROT_K5_EMB.replace(".txt", ".KPopVectors.txt")
            if os.path.exists(target):
                os.rename(target, PROT_K5_EMB)
            else:
                print(f"FAIL: couldn't materialise {PROT_K5_EMB}", file=sys.stderr)
                sys.exit(1)
    sweep_real_subsamples()
    sweep_synthetic()


if __name__ == "__main__":
    main()
