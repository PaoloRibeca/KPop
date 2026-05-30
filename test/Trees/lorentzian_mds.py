#!/usr/bin/env python3
"""Standalone prototype of Lorentzian (hyperbolic) classical MDS.

Hyperboloid model: points x in R^{1,d} with Lorentzian inner product
    <x,y>_L = -x0*y0 + sum_{k>=1} xk*yk
on the upper sheet <x,x>_L = -1, x0 > 0.  Hyperbolic distance is
    d(x,y) = acosh(-<x,y>_L).

Lorentzian MDS (the hyperbolic analogue of classical/Torgerson MDS):
given a target distance matrix D, the Lorentzian Gram matrix is
    M_ij = <x_i,x_j>_L = -cosh(D_ij).
For genuine hyperboloid points M has signature (1, d): one negative
eigenvalue (time-like) and d positive (space-like).  So eigendecompose
M, take the negative eigenpair for the time coordinate and the top-d
positive eigenpairs for the spatial coordinates.

Two tests:
  A. round-trip on TRUE hyperboloid points  -> should be ~machine eps
  B. embed EUCLIDEAN distances at various curvature scales c and target
     dimensions d  -> measures the distortion (and whether curvature and
     fidelity can be high simultaneously, the "decoupling" question).
"""

import numpy as np


# ---------- hyperboloid helpers ----------

def lorentz_gram(X):
    """M_ij = <x_i, x_j>_L for rows x_i of X (n, dim+1)."""
    J = np.ones(X.shape[1]); J[0] = -1.0
    return (X * J) @ X.T

def lorentz_sqnorm(X):
    J = np.ones(X.shape[1]); J[0] = -1.0
    return np.sum((X * J) * X, axis=1)

def hyp_dist_matrix(X):
    G = -lorentz_gram(X)
    np.fill_diagonal(G, 1.0)
    G = np.maximum(G, 1.0)               # acosh domain guard
    return np.arccosh(G)


# ---------- the method ----------

def lorentzian_mds(D, d):
    """Embed distance matrix D (n,n) into the d-dim hyperboloid.
    Returns (coords (n,d+1), eigenvalues ascending, diagnostics dict)."""
    n = D.shape[0]
    M = -np.cosh(D)                      # Lorentzian Gram target
    w, U = np.linalg.eigh(M)             # ascending eigenvalues, M symmetric
    lam_time = w[0]                      # most negative -> time-like
    u_time = U[:, 0]
    pos = np.argsort(w)[::-1][:d]        # d largest -> space-like
    coords = np.zeros((n, d + 1))
    coords[:, 0] = np.sqrt(max(-lam_time, 0.0)) * u_time
    for k, idx in enumerate(pos):
        coords[:, k + 1] = np.sqrt(max(w[idx], 0.0)) * U[:, idx]
    # the time eigenvector's global sign is arbitrary; fix it so the bulk
    # of the points land on the upper sheet (x0 > 0)
    if np.sum(coords[:, 0]) < 0:
        coords[:, 0] *= -1.0
    # project onto the sheet <x,x>_L = -1 where possible
    sq = lorentz_sqnorm(coords)
    n_offsheet = int(np.sum(sq >= -1e-9))   # points that fail to be time-like
    scale = np.where(sq < -1e-12, 1.0 / np.sqrt(-np.minimum(sq, -1e-12)), 1.0)
    coords = coords * scale[:, None]
    # any point left on the lower sheet -> flip onto the upper one
    coords[coords[:, 0] < 0] *= -1.0
    diag = {
        "lam_time": lam_time,
        "n_neg_eig": int(np.sum(w < -1e-9)),
        "n_pos_eig": int(np.sum(w > 1e-9)),
        "n_offsheet": n_offsheet,
        "top_pos": np.sort(w)[::-1][:min(d + 2, n)],
    }
    return coords, w, diag


def fit_stats(D_target, D_hat):
    iu = np.triu_indices_from(D_target, k=1)
    a = D_target[iu]; b = D_hat[iu]
    pear = np.corrcoef(a, b)[0, 1]
    # relative stress (scale-free): min over a global scale alpha of
    # ||alpha*b - a|| / ||a||
    alpha = np.dot(a, b) / np.dot(b, b) if np.dot(b, b) > 0 else 0.0
    stress = np.linalg.norm(alpha * b - a) / np.linalg.norm(a)
    return pear, stress


# ---------- Test A: round trip on true hyperboloid points ----------

def test_A(n=60, d=5, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal((n, d)) * 0.7      # spatial parts
    x0 = np.sqrt(1.0 + np.sum(v * v, axis=1))
    X = np.column_stack([x0, v])               # true hyperboloid points
    D = hyp_dist_matrix(X)
    Xr, w, diag = lorentzian_mds(D, d)
    Dr = hyp_dist_matrix(Xr)
    pear, stress = fit_stats(D, Dr)
    print(f"=== Test A: round-trip on TRUE hyperboloid points (n={n}, d={d}) ===")
    print(f"  signature: {diag['n_neg_eig']} negative, {diag['n_pos_eig']} positive eigenvalues "
          f"(expected 1 / {d})")
    print(f"  off-sheet points: {diag['n_offsheet']}")
    print(f"  recovered-vs-true distances: Pearson={pear:.6f}, rel.stress={stress:.2e}")
    print(f"  max abs distance error: {np.max(np.abs(D - Dr)):.3e}")
    print()


# ---------- Test B: embed Euclidean distances at various curvature/dim ----------

def test_B(n=80, p=10, seed=1):
    rng = np.random.default_rng(seed)
    Y = rng.standard_normal((n, p))
    De = np.sqrt(((Y[:, None, :] - Y[None, :, :]) ** 2).sum(-1))
    De /= De[np.triu_indices_from(De, 1)].mean()   # normalise mean pair dist to 1
    print(f"=== Test B: embed EUCLIDEAN distances into hyperbolic (n={n}, source dim={p}) ===")
    print("  target D = c * D_euclid ; measure fidelity of recovered hyp distances vs target")
    print(f"  {'c':>5} {'d':>4} {'neg':>4} {'offsheet':>9} {'Pearson':>9} {'rel.stress':>11}")
    for c in (0.1, 0.5, 1.0, 2.0):
        for d in (2, 5, 10):
            D = c * De
            Xr, w, diag = lorentzian_mds(D, d)
            Dr = hyp_dist_matrix(Xr)
            pear, stress = fit_stats(D, Dr)
            print(f"  {c:>5.1f} {d:>4d} {diag['n_neg_eig']:>4d} {diag['n_offsheet']:>9d} "
                  f"{pear:>9.4f} {stress:>11.4f}")
    print()


# ---------- Test C: embed a TREE (patristic) metric ----------

def random_tree_patristic(n, seed):
    """Patristic (path-length) leaf distances of a random binary tree."""
    rng = np.random.default_rng(seed)
    nodes = list(range(n)); nid = n; edges = []
    while len(nodes) > 1:
        i, j = rng.choice(len(nodes), size=2, replace=False)
        a, b = nodes[i], nodes[j]; p = nid; nid += 1
        edges.append((p, a, rng.uniform(0.1, 1.0)))
        edges.append((p, b, rng.uniform(0.1, 1.0)))
        nodes = [nodes[k] for k in range(len(nodes)) if k not in (i, j)] + [p]
    N = nid; INF = 1e9
    G = np.full((N, N), INF); np.fill_diagonal(G, 0.0)
    for u, v, w in edges:
        G[u, v] = w; G[v, u] = w
    for k in range(N):                    # Floyd-Warshall (tree -> unique paths)
        G = np.minimum(G, G[:, k:k + 1] + G[k:k + 1, :])
    return G[:n, :n]

def test_C(n=80, seed=2):
    D = random_tree_patristic(n, seed)
    D /= D[np.triu_indices_from(D, 1)].mean()
    print(f"=== Test C: embed a TREE (patristic) metric into hyperbolic (n={n}) ===")
    print("  the case hyperbolic space is *supposed* to hold well (Sarkar)")
    print(f"  {'c':>5} {'d':>4} {'neg':>4} {'offsheet':>9} {'Pearson':>9} {'rel.stress':>11}")
    for c in (0.5, 1.0, 2.0, 4.0):
        for d in (2, 5, 10):
            Dc = c * D
            Xr, w, diag = lorentzian_mds(Dc, d)
            Dr = hyp_dist_matrix(Xr)
            pear, stress = fit_stats(Dc, Dr)
            print(f"  {c:>5.1f} {d:>4d} {diag['n_neg_eig']:>4d} {diag['n_offsheet']:>9d} "
                  f"{pear:>9.4f} {stress:>11.4f}")
    print()


# ---------- Test D: embed a REAL KPop distance matrix ----------

def load_kpop_dmatrix(path):
    rows = []
    with open(path) as f:
        next(f)                                  # header of names
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows.append([float(v) for v in parts[1:]])   # drop row name
    return np.array(rows)

def test_D(path):
    D = load_kpop_dmatrix(path)
    n = D.shape[0]
    D = 0.5 * (D + D.T)                           # symmetrise (paranoia)
    mean = D[np.triu_indices_from(D, 1)].mean()
    D /= mean
    print(f"=== Test D: embed a REAL KPop distance matrix ({path}, n={n}) ===")
    print("  the actual near-tree distances sparse-NJ uses (CA-Euclidean,")
    print("  ~0.996 Pearson vs patristic) -- Test-B (random Euclidean) or Test-C (tree) regime?")
    print(f"  {'c':>5} {'d':>4} {'neg':>4} {'offsheet':>9} {'Pearson':>9} {'rel.stress':>11}")
    for c in (0.5, 1.0, 2.0, 4.0):
        for d in (2, 5, 10, 20):
            Dc = c * D
            Xr, w, diag = lorentzian_mds(Dc, d)
            Dr = hyp_dist_matrix(Xr)
            pear, stress = fit_stats(Dc, Dr)
            print(f"  {c:>5.1f} {d:>4d} {diag['n_neg_eig']:>4d} {diag['n_offsheet']:>9d} "
                  f"{pear:>9.4f} {stress:>11.4f}")
    print()


if __name__ == "__main__":
    import sys
    test_A()
    test_B()
    test_C()
    if len(sys.argv) > 1:
        test_D(sys.argv[1])
