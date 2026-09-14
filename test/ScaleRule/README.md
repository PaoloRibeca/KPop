# The scale-rule harness

The programs behind `DocsYard/docs/KPop/Clustering/KPop-autotuner-scale-rule.html`, which proposes
how `KPop-autotuner` should choose the number of axes it searches in, which level of structure it
holds the search to, and how the search should resist both atomising and over-merging. They are kept
here because the document cites their measurements and they were written in a scratch directory that
does not survive; nothing here is part of the build or of any test suite.

**None of this is implemented in `lib/`.** These are prototypes that measure what a proposal would
do. Where one duplicates something the libraries already have, that is a fault of the prototype and
is noted below.

## Building and running

Each program is standalone OCaml over `kio.ml`, a reader for the text exports:

    ocamlfind ocamlopt -O3 -unsafe -package unix -linkpkg kio.ml <program>.ml -o <program>

The inputs are text exports of a twisted register and its inertia, as
`KPopTwistDB -i t <prefix> -O t <prefix>.all` writes them, plus tab-separated label files
(`<name>\t<class>`). The norovirus runs used `NailIt/test/data/denovo/{vp1,rdrp,vp2}.cdc` and
`vp1.refined`.

## What each program does

| File | What it is |
|---|---|
| `kio.ml` | Reads a twisted register, an inertia vector and a label file. **Duplicates the library**: `Matrix`/`Twisted` already parse these. |
| `detect3.ml` | The calibrated valley detector of §1, and the axis rules of §2. Persistence simplification judges a bump against the modes the gap separates; each surviving trough is calibrated by resampling the sequences. Prints kept valleys per rung with radius, share below, z and reappearance, plus the pick under each candidate rule. |
| `vdump.ml` | The pair-distance histogram for the plot page: raw counts, counts of same-class pairs, two smoothed curves, and — at radii `detect3` found — the share of pairs below and of each label level caught. Detects nothing itself. |
| `vtree.ml` | Average linkage (UPGMA) by nearest-neighbour chain with Lance-Williams updates; the join-height histogram; cuts at given radii and over a 120-height grid, each scored against the labels. |
| `treecut.ml` | The first version of the same, superseded by `vtree.ml` and kept only because the document's earliest tree numbers came from it. |
| `splitcmp.ml` | Compares a valley cut with a KPopPhylo sparse-NJ tree: splits, Robinson-Foulds, per-clade Jaccard. **Duplicates the library**: `Trees.Newick` parses *and writes* Newick and `Trees.Splits` compares splits. |
| `guard2.ml` | The minimum-cluster guard of §4.3 as the document measured it: four variants scored against some sixty partitions on six embeddings, at thresholds 4, 8 and 16, on both silhouettes. Kept because §4.3's numbers came from it; `guard3.ml` supersedes it. |
| `cs.ml` | The couples silhouette, the push-back scale, and split seeding. |
| `ap.ml` | Whether one pooled variance factor can calibrate every trough. It cannot: 0.34 to 3.7 times. |
| `phylo.sh` | Builds the projections and the sparse-NJ trees (~1 h; not repeated when only the detector changes). |
| `compare.sh` | Runs `splitcmp` over the seven cases and writes the page's `phylo.json`. |
| `regen.sh` | Regenerates both histogram dumps on the union of the page's grid and the doubling ladder. |
| `run2.sh`, `beats.sh`, `tables.sh` | Drive `guard2` and tabulate it. |
| `guard3.ml` | The same guards on a corpus built once and apart from the thresholds (2 to 32): attacks at sizes both absolute and relative to n, each partition scored once however many constructions reach it, and no partition built by projecting through a guard. Optionally subsamples the labelled rows (`uniform`, `strat`, `classes`) and coarsens or refines the reference in place (`merge<j>`, `split`). Writes a 20-column TSV that a complete run ends with a row of family `end`. |
| `mstar.sh` | m\*, the largest single group against the rest that beats the labels with no guard at all. A lower bound: every probe is a periphery-against-bulk split seeded at the periphery. |
| `curve.sh` | The smallest threshold at which no unguarded-winning attack beats the labels under G3b. It lands one defence rung above m\* and adds nothing to it; its header says why. |
| `runexp.sh` | Drives `guard3` over the three embeddings: the anchor, proportional thinning (arm A), whole classes dropped (arm B, not run) and the reference coarsened and refined in place (arm C). |
| `inputs/` | How the three embeddings were built: `build.sh` (VP1, RdRp) and `build-vp2.sh`, verbatim with the scratch paths they ran against, and the spectra each twister was computed from (`*.T_rand.sel`). The spectra databases are not kept. |
| `results/` | `g0.tsv.xz`, the unguarded rows of all 54 complete runs, which is enough to recompute every m\* and every comparison of references; `mstar.txt` and `curve-anchors.txt`, the two scripts' output. The guarded rows (81 MB) are not kept; `runexp.sh` regenerates them. |
| `valley-page-template.html` | The plot page, with `__DATA__`, `__TREE__` and `__PHYLO__` standing in for the three JSON payloads. |

## What is new here

The work these programs exist for, and the part worth keeping whatever happens to the rest:

- **A valley detector that can be argued with.** Candidate troughs are simplified in order of
  relative size, and — this is the point — an interior peak is judged against the modes beyond the
  troughs it sits between rather than against the floor beside it, so a bump inside a gap is
  measured against the modes the gap separates and no longer splits it in two. Each surviving
  trough is then calibrated on its own, by resampling the *sequences* and reweighting every pair by
  the product of the two multiplicities, which gives a prominence in standard deviations rather
  than a threshold someone chose. Replacing the half-height rule with it changed 29 of 105 rung
  valley lists and turned RdRp's 5-axis pair at z 5.5 and z 13.3 into one valley at z 27, on the
  labelled P-type level. `detect3.ml`.
- **Choosing the number of axes from the structure the embedding resolves.** Doubling from one
  axis, counting the valleys at or below half the pairs, and taking the right end of the first
  plateau. `detect3.ml` evaluates that against the alternatives on seven embeddings.
- **Why a minimum-cluster guard is the wrong defence.** Of four variants, barring small clusters
  from being anyone's nearest other cluster is *worse* than no guard (atomised pairs then score
  0.81–0.96 against the labels' 0.26–0.62), and charging each member of a sub-threshold cluster −1
  is the best of them (`guard2.ml`). Measured properly, none is worth having: see *The
  minimum-cluster guard, measured* below (`guard3.ml`, `mstar.sh`).
- **The bridge between a valley and a tree level**, which is what the whole question was about:
  cutting the average-linkage tree at each valley and scoring the cut against the labels, set
  against the best cut anywhere on the tree. `vtree.ml`.

## The minimum-cluster guard, measured

The question was what threshold s a minimum-cluster guard should use. Measured on the three
norovirus embeddings at the axes the rule picks (VP1 d = 20, RdRp d = 5, VP2 d = 20; `guard3.ml`,
`results/`), the answer is that no threshold is safe and no guard should be added.

**What beats the labels is not a small cluster.** With no guard at all, the partitions that outscore
the curated labels are two-cluster splits of periphery against bulk: the k points farthest from the
centroid against the rest (f_p 0.98–0.999), coarser than the labels. One outlier against the rest
already wins by +0.048/+0.071 on VP1 and +0.051/+0.079 on VP2 (classical/simplified), and loses only
on RdRp at d = 5. The largest such group that still wins, m\* (`mstar.sh`, 95% draw):

| | classical | simplified |
|---|---|---|
| VP1 | 2 | 4 |
| RdRp | none | none |
| VP2 | ≥ 32 | 16 |

Every m\* is a lower bound: all probes start at the periphery, only eight seeds are tried, and the
size ladder has gaps.

**No constant threshold can work.** Under G3b a group of exactly s members is eligible, and on VP2
classical the s farthest points beat the guarded labels at s = 16, 22 and 32. A threshold above 32
would charge 27 of 39 VP1 classes, 45 of 53 RdRp classes and 27 of 36 VP2 classes as noise; s = 8
already charges 15, 32 and 14 at full n.

**The labels are not the optimum in either direction.** Splitting the largest class raises the
unguarded score (VP1 +0.10, VP2 +0.09, classical), but five nearest-pair merges raise it too, in all
six cells. Refining the reference in place (`refmode split`) makes every m\* vanish, but only
arithmetically: the probes score identically under every reference mode, and only the reference's
own score moves.

**Scale.** Where m\* is measurable it grows with n — VP2 simplified per draw: 3–4 at n = 377, 5–7 at
755, 8–11 at 1510, 16 at 2869 — roughly in proportion. On VP1 and RdRp most cells are empty, and no
functional form can be read off. The draws are nested within each seed and reuse the full register's
embedding.

**Geometry dominates.** At full dimension (d = 243, 199, 219, the older corpus of `guard2.ml`) one
outlier against the rest beats the labels by +0.10 to +0.20 in all six cells, RdRp included.

**What follows.** The degenerate to defend against is a periphery-against-bulk bipartition. The
push-back as drafted, `silhouette − λ·max(0, ln(f_v/f_p) − ln 2)`, penalises only partitions finer
than the valley, so these pass it untouched. Whether the search under the level rule ever reaches
them, and whether a two-sided push-back would stop it, needs autotuner runs: this rig scores
partitions and cannot say.

## What is standard, and only absent from these libraries

Worth contributing, but nobody should read them as novel: average linkage by nearest-neighbour
chain with Lance-Williams updates is textbook, and so are the adjusted Rand index, homogeneity and
completeness — the last three are in `Clustering.compare_partitions` already.

## What the libraries already provide

Written before this was checked, which is the lesson rather than an aside:

- `Clustering.compare_partitions` returns the (ARI, homogeneity, completeness) triple that
  `treecut.ml` and `vtree.ml` each compute by hand — judged partition first, reference second.
- `Trees.Newick` and `Trees.Splits` do what `splitcmp.ml` reimplements.
- A `Cophenetic` program already exists for join-height work.
- `Numbers` has mean, variance, sample variance and median.
- `test/repl.sh` gives a toplevel with `BiOCamLib` and `KPop` linked, which is where a probe
  belongs.

## Reproducing the document's numbers

The detector over seven embeddings, fifteen rungs, fifty resamples, single-threaded (~3 min):

    ./detect3 <prefix>.all.KPopTwisted.txt <prefix>.all.KPopInertia.txt <labels> <tag> \
      1,2,4,5,8,10,16,20,32,40,64,80,120,128,999 50 1 > detect3.out 2> detect3.cand

`detect3.cand` carries one line per candidate trough: tag, rung, radius, share below, z,
reappearance, relative depth, and whether it was kept. `regen.sh` turns those radii into the
histogram dumps; `compare.sh` cuts the NJ trees at them.

The minimum-cluster experiment, single-threaded (a full-size run is about three minutes). Build in
the work directory, since compiling here would leave the objects in the repository:

    cp kio.ml guard3.ml "$WORK" && (cd "$WORK" && ocamlfind ocamlopt -O3 -unsafe -package unix \
      -linkpkg kio.ml guard3.ml -o guard3)
    export WORK TWISTED=<directory with the three .all.KPopTwisted.txt exports>
    ./runexp.sh anchor
    for f in 0.5 0.25 0.125; do ./runexp.sh A $f 5; done
    ./runexp.sh C 0.95
    ./mstar.sh "$WORK"/exp2.*.tsv

The anchors in `results/` were run by hand with the same arguments, so their tags lack the mode
that `runexp.sh` puts in; the numbers are the same.
