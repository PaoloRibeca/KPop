#!/usr/bin/env python3
"""Simple phylogenetic sequence simulator.

Generates a random binary tree with n leaves and exponentially-distributed
branch lengths, then simulates DNA sequences along the tree under the
Jukes-Cantor model.  Outputs:
  - <prefix>.fasta : multi-FASTA, one record per leaf
  - <prefix>.tree.nwk : the true tree in Newick
  - <prefix>.class.txt : per-sample class table (each leaf gets its own class,
                        in the KPopCountDB tabular format expected by -m)

Memory-conscious: holds n * seqlen bytes of sequence state at once, no
intermediate copies.  At n=1000, L=500: 500 KB.  At n=1000, L=2000: 2 MB.

Usage:
  simulate_phylo.py <prefix> <n_leaves> <seq_length> <branch_scale> [<seed>]
"""
import random
import sys


BASES = "ACGT"


def gen_random_tree(n, rng):
    """Random binary tree with n leaves via a Yule (pure-birth) process.
    Each active leaf splits at constant rate; branch lengths are the time
    each branch survives between speciations.  Produces trees with a more
    uniform branch-length distribution than the edge-splitting recipe
    (which gives lots of near-zero branches near the root).

    Returns the (internal) root Node; leaves named S0..S{n-1}."""
    class Node:
        __slots__ = ("name", "children", "_pending_blen")
        def __init__(self, name=None):
            self.name = name
            self.children = []
            self._pending_blen = 0.0

    # Start with root + two child leaves at time 0
    root = Node()
    leaves_pool = [Node(name=f"S{i}") for i in range(n)]
    next_id = 2
    s0, s1 = leaves_pool[0], leaves_pool[1]
    s0._pending_blen = 0.0
    s1._pending_blen = 0.0
    root.children = [(s0, 0.0), (s1, 0.0)]
    active_leaves = [s0, s1]
    while next_id < n:
        # Time until next split (Yule: rate = #active leaves * lambda; lambda=1)
        rate = float(len(active_leaves))
        dt = rng.expovariate(rate)
        # All active leaves age by dt
        for leaf in active_leaves:
            leaf._pending_blen += dt
        # Pick one to split
        splitter = rng.choice(active_leaves)
        # Finalize its pending branch length
        # Find the (parent, ci) for splitter
        new_internal = Node()
        new_internal._pending_blen = 0.0
        # Splitter becomes an internal node: replace it everywhere
        # First, find its parent and replace it
        def find_parent(node):
            for k, (c, _) in enumerate(node.children):
                if c is splitter:
                    return node, k
                if c.children:
                    r = find_parent(c)
                    if r is not None:
                        return r
            return None
        parent, ci = find_parent(root)
        parent.children[ci] = (new_internal, splitter._pending_blen)
        # Now new_internal gets splitter (as a fresh leaf) plus a new leaf
        splitter._pending_blen = 0.0
        new_leaf = leaves_pool[next_id]
        new_leaf._pending_blen = 0.0
        next_id += 1
        new_internal.children = [(splitter, 0.0), (new_leaf, 0.0)]
        # Update active set: splitter stays, new_leaf joins
        active_leaves.append(new_leaf)
    # Finalize: any remaining pending branch lengths get applied as
    # the splitter's branch to its parent.  Since active leaves still
    # have pending time, walk one more "epoch" to extend each branch
    # by a small amount before finalising.  Simpler: just write the
    # pending values as branch lengths from each leaf's parent.
    def finalize(node):
        for k, (c, bl) in enumerate(node.children):
            if not c.children:  # leaf
                node.children[k] = (c, bl + c._pending_blen)
            else:
                finalize(c)
                node.children[k] = (c, bl + c._pending_blen)
    finalize(root)
    return root


def to_newick(root):
    def rec(node):
        if not node.children:
            return node.name
        inner = ",".join(f"{rec(c)}:{bl:.6f}" for c, bl in node.children)
        return f"({inner})"
    return rec(root) + ";"


def jc_mutate(seq, branch_len, mu, rng):
    """Apply Jukes-Cantor substitutions: each site has probability
    p_change = 0.75 * (1 - exp(-4 * mu * branch_len / 3)) of being
    substituted to one of the 3 other bases uniformly."""
    if branch_len <= 0:
        return seq
    p_change = 0.75 * (1.0 - 2.71828182845905 ** (-4 * mu * branch_len / 3.0))
    new_seq = list(seq)
    for i in range(len(new_seq)):
        if rng.random() < p_change:
            old = new_seq[i]
            others = [b for b in BASES if b != old]
            new_seq[i] = rng.choice(others)
    return "".join(new_seq)


def simulate(root, seq_length, mu, rng):
    """Walk the tree starting from a random root sequence, accumulating
    substitutions along each branch.  Returns dict leaf_name -> sequence."""
    root_seq = "".join(rng.choice(BASES) for _ in range(seq_length))
    leaves = {}
    def dfs(node, seq):
        if not node.children:
            leaves[node.name] = seq
        else:
            for child, blen in node.children:
                child_seq = jc_mutate(seq, blen, mu, rng)
                dfs(child, child_seq)
    dfs(root, root_seq)
    return leaves


def main():
    if len(sys.argv) < 5:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    prefix = sys.argv[1]
    n = int(sys.argv[2])
    seqlen = int(sys.argv[3])
    branch_scale = float(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) >= 6 else 42
    mu = float(sys.argv[6]) if len(sys.argv) >= 7 else 0.05

    rng = random.Random(seed)
    print(f"Generating tree (n={n}, seed={seed})...", file=sys.stderr)
    root = gen_random_tree(n, rng)

    # Scale all branch lengths
    def rescale(node):
        for k, (c, bl) in enumerate(node.children):
            node.children[k] = (c, bl * branch_scale)
            rescale(c)
    rescale(root)

    print(f"Simulating JC69 sequences (length={seqlen}, mu={mu})...",
          file=sys.stderr)
    seqs = simulate(root, seqlen, mu, rng)

    # Write FASTA
    fasta_path = f"{prefix}.fasta"
    with open(fasta_path, "w") as f:
        for name in sorted(seqs.keys(), key=lambda s: int(s[1:])):
            f.write(f">{name}\n{seqs[name]}\n")
    print(f"Wrote {fasta_path} ({len(seqs)} seqs, {seqlen} bp each)",
          file=sys.stderr)

    # Write CLASS table (each sample its own class so the CA twister gets
    # full per-sample resolution).  Class labels are "CLS_<name>" so they
    # are guaranteed not to collide with any sample name in the database.
    class_path = f"{prefix}.class.txt"
    sorted_names = sorted(seqs.keys(), key=lambda s: int(s[1:]))
    with open(class_path, "w") as f:
        f.write("\t" + "\t".join(sorted_names) + "\n")
        f.write("CLASS\t" + "\t".join(f"CLS_{nm}" for nm in sorted_names) + "\n")
    print(f"Wrote {class_path}", file=sys.stderr)

    # Write true tree
    tree_path = f"{prefix}.tree.nwk"
    with open(tree_path, "w") as f:
        f.write(to_newick(root) + "\n")
    print(f"Wrote {tree_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
