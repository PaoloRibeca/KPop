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
| `guard2.ml` | The minimum-cluster guard of §4.3: four variants scored against some sixty partitions on six embeddings, at three thresholds, on both silhouettes. |
| `cs.ml` | The couples silhouette, the push-back scale, and split seeding. |
| `ap.ml` | Whether one pooled variance factor can calibrate every trough. It cannot: 0.34 to 3.7 times. |
| `phylo.sh` | Builds the projections and the sparse-NJ trees (~1 h; not repeated when only the detector changes). |
| `compare.sh` | Runs `splitcmp` over the seven cases and writes the page's `phylo.json`. |
| `regen.sh` | Regenerates both histogram dumps on the union of the page's grid and the doubling ladder. |
| `run2.sh`, `beats.sh`, `tables.sh` | Drive `guard2` and tabulate it. |
| `valley-page-template.html` | The plot page, with `__DATA__`, `__TREE__` and `__PHYLO__` standing in for the three JSON payloads. |

## What the libraries already provide

Written before this was checked, which is the lesson rather than an aside:

- `Clustering.compare_partitions` returns the (ARI, homogeneity, completeness) triple that
  `treecut.ml` and `vtree.ml` each compute by hand — judged partition first, reference second.
- `Trees.Newick` and `Trees.Splits` do what `splitcmp.ml` reimplements.
- A `Cophenetic` program already exists for join-height work.
- `Numbers` has mean, variance, sample variance and median.
- `test/repl.sh` gives a toplevel with `BiOCamLib` and `KPop` linked, which is where a probe
  belongs.

Genuinely absent, and worth contributing when the proposal is settled: average linkage with
cophenetic heights; cutting a tree at a height and scoring the cut; the calibrated detector; the
minimum-cluster guard.

## Reproducing the document's numbers

The detector over seven embeddings, fifteen rungs, fifty resamples, single-threaded (~3 min):

    ./detect3 <prefix>.all.KPopTwisted.txt <prefix>.all.KPopInertia.txt <labels> <tag> \
      1,2,4,5,8,10,16,20,32,40,64,80,120,128,999 50 1 > detect3.out 2> detect3.cand

`detect3.cand` carries one line per candidate trough: tag, rung, radius, share below, z,
reappearance, relative depth, and whether it was kept. `regen.sh` turns those radii into the
histogram dumps; `compare.sh` cuts the NJ trees at them.
