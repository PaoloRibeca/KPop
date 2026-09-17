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
| `detect3.ml` | The calibrated valley detector of §1, and the axis rules of §2. Persistence simplification judges a bump against the modes the gap separates; each surviving trough is calibrated by resampling the sequences. Prints kept valleys per rung with radius, share below, z and reappearance, plus the pick under each candidate rule. The rules run on three ladders read off one set of detections: the 5, 10, 20, 40 … one, the doubling one, and the half powers of two the autotuner now takes — 1, 2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64, 91, 128, 181 and the most axes there are. |
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
| `guard3.ml` | The same guards on a corpus built once and apart from the thresholds (2 to 32): attacks at sizes both absolute and relative to n, duplicates among the shared constructions scored once (the best-of-eight-seeds probes write their own rows, so they can repeat another construction), and no partition built by projecting through a guard. Optionally subsamples the labelled rows (`uniform`, `strat`, `classes`) and coarsens or refines the reference in place (`merge<j>`, `split`). Writes a 20-column TSV that a complete run ends with a row of family `end`. |
| `mstar.sh` | m\*, the largest single group against the rest that beats the labels with no guard at all. A lower bound: every probe is a periphery-against-bulk split seeded at the periphery. |
| `curve.sh` | The smallest threshold at which no unguarded-winning attack beats the labels under G3b. It lands one defence rung above m\* and adds nothing to it; its header says why. |
| `runexp.sh` | Drives `guard3` over the three embeddings: the anchor, proportional thinning (arm A), whole classes dropped (arm B, not run) and the reference coarsened and refined in place (arm C). |
| `inputs/` | How the three embeddings were built: `build.sh` (VP1, RdRp) and `build-vp2.sh`, verbatim with the scratch paths they ran against, and the spectra each twister was computed from (`*.T_rand.sel`). The spectra databases are not kept. |
| `results/` | `g0.tsv.xz`, the unguarded rows of all 54 complete runs, which is enough to recompute every m\* and every comparison of references; `mstar.txt` and `curve-anchors.txt`, the two scripts' output. The guarded rows (81 MB) are not kept; `runexp.sh` regenerates them. |
| `valley-page-template.html` | The plot page, with `__DATA__`, `__TREE__` and `__PHYLO__` standing in for the three JSON payloads. |
| `stage1/` | Stage 1 run with the implemented autotuner (see *Stage 1, run* below). `stage1-run` drives the 18 runs two at a time from a frozen copy of the binary; `control-run` the six runs of the fixed-rung control on RdRp; `stage1-score` scores every round against the CDC labels and takes each run's median over its last four rounds; `stage1-summary` gives per corpus and arm the median over seeds, its spread, and whether each arm differs from `today` by more than twice the larger spread. `runs/` and `control/` hold the logs (xz), final partitions, scores and provenance; `verify/` the independent rechecks of stage 1. The spectra databases are not kept. |

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
norovirus embeddings next to the axes the rule picks (VP1 d = 20, RdRp d = 5, VP2 d = 20; `guard3.ml`,
`results/`), the answer is that no threshold works across the three and no guard should be added.

**What beats the labels is the periphery set apart from the bulk.** With no guard at all, every
degenerate probe that outscores the curated labels sets some of the outermost points apart from a bulk
holding almost every spectrum, as singletons or as one group (f_p 0.968–0.999 over all 54 runs),
coarser than the labels. As singletons they are the small clusters a guard is meant for; as one group
they are an eligible cluster once they number s. One outlier against the rest already wins by
+0.048/+0.070 on VP1 and +0.051/+0.079 on VP2 (classical/simplified), and loses only on RdRp at d = 5.
Near-copies of the labels win too (below), and no small cluster inside a class was ever probed. The largest such group that still wins, m\* (`mstar.sh`, 95% draw):

| | classical | simplified |
|---|---|---|
| VP1 | 2 | 4 |
| RdRp | none | none |
| VP2 | ≥ 32 | 16 |

Every m\* is a lower bound: all probes start at the periphery, only eight seeds are tried, and the
size ladder has gaps.

**No constant threshold can work.** Under G3b a group of exactly s members is eligible. On VP2 one
beats the guarded labels at every threshold tried, on one silhouette or the other: on the simplified
silhouette at every s up to 22 (+0.037 at s = 8), on the classical one at s = 2–6, 16, 22 and 32,
resisting only at s = 8 and 11, by 0.001. On VP1 no degenerate probe beats the guarded labels above
s = 4, and on RdRp none ever does, so the size that works depends on the data and finding it takes
labels. A threshold large enough for VP2 would charge 27 of 39 VP1 classes, 45 of 53 RdRp classes and
27 of 36 VP2 classes as noise at s = 32; s = 8 already charges 15, 32 and 14 at full n.

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

**What follows.** The degenerate partitions to defend against set the periphery apart from the bulk.
The push-back as drafted, `silhouette − λ·max(0, ln(f_v/f_p) − ln 2)`, penalises only partitions
finer than the valley, so these pass it untouched. Whether the search under the level rule ever reaches
them, and whether a two-sided push-back would stop it, needs autotuner runs: this rig scores
partitions and cannot say.

## Stage 1, run

The document's stage 1 as the autotuner now implements it (KPop.Claude `dd13d32`..`bc19be9`): two arms,
three corpora, three seeds (17, 18, 19), eight rounds each, one thread at nice 19. `today` is the
half-height detector with `--dimensions-inertia 0.5`; `auto` is the calibrated detector with
`--dimensions auto`. Both take projection and partition samples of 244 (VP1), 200 (RdRp) and 220 (VP2)
spectra. All 18 runs finished.

Median over seeds (spread across seeds), homogeneity / completeness against the CDC types, each run
contributing its median over rounds 5–8:

| | today | auto | auto − today (bar) |
|---|---|---|---|
| VP1 | 0.941 (0.040) / 0.659 (0.028) | 0.907 (0.168) / 0.658 (0.058) | −0.034 (0.337) / −0.001 (0.116) |
| RdRp | 0.918 (0.014) / 0.820 (0.024) | 0.698 (0.015) / 0.957 (0.006) | −0.220 (0.029) / +0.137 (0.049) |
| VP2 | 0.926 (0.008) / 0.781 (0.060) | 0.962 (0.063) / 0.719 (0.259) | +0.036 (0.126) / −0.063 (0.517) |

The bar is twice the larger spread. Every figure reproduces through `KPopTwistDB --clusters-compare`.
`rdrp.cdc` labels 9 sequences `-`, which the scorer counts as a class in both arms.

**RdRp: one change, coarser partitions.** Under `auto` the ladder picks 4 axes or fewer in 18 of 24
rounds, and those rounds return 3–11 clusters for 53 P-types. The clusters are near-exact unions of
P-types (pair recall 0.95–0.996), merging even across genogroups, so completeness rises and homogeneity
falls to what the cluster count allows. In 13 of the 18 a rung of 1–4 axes counts more valleys at or
below half the pairs than any rung of 16 axes or more; in the other 5 a higher rung ties and the
first-plateau tie-break, which favours fewer axes, decides. The winning low-rung valleys are mostly weak
(median z 9, often a few points of share apart), where rungs of 8–32 axes usually show one strong valley
(z 20–28) at 8–15% of pairs, and rungs of 64–199 axes one or two at 6–9%. The logs list only kept
valleys, so whether the low rungs' extra valleys are the gap fragmentation seen on synthetic groups
cannot be settled from them.

**VP1 and VP2: undetermined, not equal.** The bars there are set by `auto`'s own spread over three seeds,
each driven by one seed, and are too wide for a difference smaller than a third to a half of the scale to
count.

**The rung moves between rounds, not within one.** Pooled over seeds, `auto` searches VP1 in 8–64 axes,
RdRp in 1–64 and VP2 in 2–219 (2–128 after round 1), against `today`'s 21–26, 17–25 and 24–34; cluster
counts follow the axes. Within a round the 50 resamples agree with the pick in 71 of 72 rounds, so the
change comes from the projection being rebuilt from the previous partition. It does not explain the seed
spread: RdRp `auto` moves as much from round to round and has the smallest spread of all.

**Cost.** VP2 `auto` took 1.4–3.1 h against 20 min for `today`, and RdRp `auto` 1.6–2.9 times `today`'s
time while mostly searching in 4 axes or fewer, so the per-round ladder costs time of its own. The logs
carry no per-round timings, and the runs went two at a time on a shared machine.

**The fixed-rung control: the loss on RdRp is the pick rule's.** The calibrated detector with
`--dimensions` fixed at 16 or at 32 — rungs the ladder turned down — on the same seeds and settings
(`control-run`, `control/`):

| RdRp | homogeneity / completeness | against `today` (bar) |
|---|---|---|
| calibrated, 16 axes | 0.915 (0.017) / 0.821 (0.041) | −0.003 (0.034) / +0.001 (0.081) |
| calibrated, 32 axes | 0.937 (0.010) / 0.805 (0.022) | +0.019 (0.028) / −0.015 (0.049) |

Neither differs from `today`, and every round of every seed stays at 18–50 clusters (16 axes) or 31–77
(32 axes), homogeneity 0.86–0.97. The detector and the data are therefore not what costs `auto` on
RdRp: it is the choice of rung — counting valleys regardless of their z, and breaking ties towards fewer
axes. A fixed `--dimensions` asks the randomised decomposition for that many axes directly, where `auto`
decomposes into all of them and truncates; both keep the leading axes of the same analysis.

## The ladder's spacing, measured

The axis ladder steps by half powers of two — 1, 2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64, 91, 128,
181 and the most axes there are, each rung 2^(k/2) rounded — rather than doubling. `detect3` scores
both ladders off one set of detections, on the six embeddings the document measures, with the
misordered share of (same-type pair, different-type pair) comparisons at the pick as the judge.
`halfladder/` holds the run: `half-ladder-run`, `half.out` and the candidate lines, `half.cand.xz`.

| Embedding | Doubling: pick, misordered | Half powers: pick, misordered |
|---|---|---|
| VP1 uniform | 16, 1.70% | 23, 1.15% |
| VP1 stratified | 16, 1.87% | 23, 1.11% |
| VP1 cores | 16, 0.53% | 16, 0.53% |
| RdRp uniform | 4, 0.83% | 4, 0.83% |
| RdRp stratified | 4, 1.19% | 4, 1.19% |
| VP2 uniform | 16, 0.91% | 16, 0.91% |
| **Mean** | **1.17%** | **0.95%** |

Three things follow. The tie rule is exercised for the first time: on both VP1 samples 16 and 23
axes tie, and taking the right end of that plateau is what gains the 0.22 points, so the leftmost
end would give them back. The detector costs about 1.8 times as much a round, sixteen rungs against
nine at 243 axes, at 1–4 s a rung. And RdRp is untouched: it picks 4 axes on both samples because
4 holds three valleys and no other rung holds more, where the fixed-rung control above puts its
best real-search range at 16–32, so the finer ladder does not address that.

### On 300 sequences, where the axes matter most

Both ladders under `--dimensions auto` on the 300-sequence VP1 and RdRp subsets of
`stage1/subsample/`, three seeds each, otherwise the stage-1 settings, the frozen doubling binary
against the half-powers one (`stage1/subsample/auto/`, driven by `subsample-auto`). On 68 projected
spectra the doubling ladder offers 1, 2, 4, 8, 16, 32, 64 and 67 rungs, the half-powers one 1, 2, 3,
4, 6, 8, 11, 16, 23, 32, 45, 64 and 67.

| Corpus | Arm | Homogeneity (spread) | Completeness (spread) |
|---|---|---|---|
| RdRp | doubling | 0.870 (0.103) | 0.919 (0.062) |
| RdRp | half powers | 0.797 (0.171) | 0.941 (0.114) |
| VP1 | doubling | 0.791 (0.133) | 0.718 (0.106) |
| VP1 | half powers | 0.751 (0.009) | 0.737 (0.190) |

No difference counts: every difference is far inside twice the larger spread. What the runs do show
is that `auto` removes the stopping problem these subsets had under fixed rungs — all 12 runs
finished their 8 rounds, where 27 of the 30 fixed-rung runs stopped at "no groups at a usable
level" — and that this is `auto` itself, not the spacing, since the doubling arm finished too.

Judged the way the six big embeddings were judged, though, the finer rungs do earn their place on a
small corpus. Each subset was embedded as round 1 embeds one — a correspondence analysis of a
68-spectrum uniform sample with all 300 projected through it, leaving 67 axes — and `detect3` ran
both ladders on it (`halfladder/small/`, driven by `small-embed`):

| Embedding | Doubling: pick, misordered | Half powers: pick, misordered | Best rung by labels |
|---|---|---|---|
| VP1, 300 | 16, 3.09% | 16, 3.09% | 4, 1.53% |
| RdRp, 300 | 2, 3.81% | 3, 2.32% | 11, 0.23% |

RdRp's pick moves to a rung the doubling ladder does not have, gaining 1.49 points and becoming
unanimous over the 50 resamples where rung 2 held 94% of them. Both rules stay well above the best
rung either ladder offers, which is the price of not looking at labels, and that gap is much wider
here than on the full corpora.

The rung wanders between rounds on both ladders, which is the behaviour to fix for small sets: VP1
under doubling runs 67, 32, 4, 32, 64, 4, 32, 4 in one seed. The finer ladder keeps VP1 lower and
less extreme (mean rung 15–16 against 18–30, highest 45 against 67) and makes its homogeneity
consistent across seeds, 0.009 of spread against 0.133; on RdRp it does the opposite, one seed
sitting at 1–4 axes throughout and scoring 0.628/0.987.

## Searching for the number of axes

What `KPop-autotuner --dimensions search` does rests on five measurements, each kept under
`halfladder/`. The page drawing all of them is the artifact *Axis Ladder Evidence*.

**The silhouette falls as axes are added, and the valleys' strength says where** (`silhouette/`:
`silh.ml`, `silh-run`, one TSV per embedding, `bounds.json`). On the six full-corpus embeddings and
the two 300-sequence ones, the classical silhouette of the CDC labels rises and then falls with the
number of axes on every one of them, peaking at 3–6 axes on the full corpora, 4 and 11 on the
small ones. A partition carried unchanged from the ladder's pick peaks at its own rung or the next,
so it cannot bound anything by itself. What can is the strength of each rung's strongest usable
valley: the largest downward step in its log z, from the first rung with a usable valley on, is
significant under 4,000 shuffles of the rungs' order on seven of the eight embeddings, and the
labels' silhouette peaks at or before that step on all seven. That is `Clustering.ladder_bound`.

**No other distance or metric does better** (`distances/`: `dist.ml` mirrors `Space.Distance`,
`detect4.ml` and `silhgrid.ml` are `detect3.ml` and `silhnull.ml` with the metric and the distance as
arguments and reproduce them byte for byte under `powers(1,1,1)` and Euclidean; `grid-run`,
`grid-analyse.py`, one directory per setting). Across `flat` and `powers(1,1,1)` crossed with
Euclidean, Manhattan and angle, `flat` gives a significant step on all eight embeddings and never
before the labels' peak, but its ladder picks average 3.3–10.9% misordered on the full corpora
against 0.95% for the default. The excess of each rung's own silhouette over what the same number of
clusters gets on per-axis shuffled coordinates peaks within one rung of the labels' on at best four
of eight, so no static quantity found the number of axes inside the bound.

**Inside the bound, first rounds with different seeds agree best where the labels score best**
(`round1/`, `round1-vp2/`, `round1-vp1/`: `round1-scan`, `round1-analyse.py`, logs and one JSON per
corpus). The first round at every rung up to 64 axes, seeds 17–19. A golden-section search on the
median agreement, restricted to the bound, probes four rungs — twelve first rounds — and lands on
16 axes on RdRp (the labels' best), 11 on VP2 (best by ARI, on a plateau running 8–23) and 8 on VP1
(one below the labels' best, 11, and tied with it within the spread between seeds). Agreement is
not single-peaked: two axes and three clusters agree almost perfectly on VP1 and VP2.

**Rounds after the first make the partition worse** (`hold/`: `hold-run`, `hold-analyse.py`,
`hold.json`, logs). Held at the searched axes for eight rounds, completeness against the CDC types
falls between round 1 and round 2 on every corpus — RdRp 0.93 to 0.83, VP2 0.93 to 0.81, VP1 0.94 to
0.77 — and stays down, as it does under today's settings too.

**However the later rounds' sample is drawn** (`sample/`: `sample-run`, logs, and
`sample-options.patch` with `KPop_autotuner.ml.with-sample-options`, the two options the runs used —
`--partition-sample-allocation every-class|by-size` and `--partition-sample-members
extremes|random` — which were then removed, having helped nowhere). Shared by size, with members at
random, or both, the drop at round 2 is the same and never recovers; random members also make VP1's
rounds unsteady. Both together are close to a uniform sample, so what makes the first round good is
its sample spanning the corpus, not the way a later one uses the partition.

**What the implementation does with all that.** `--dimensions search`, run on each corpus with seed
17 and stage 1's settings, bounds itself at 32 axes on RdRp (p 0.001) and 23 on VP2 (p 0.044), finds
no significant step on VP1, probes five or six rungs, and takes the finest rung the seeds cannot
tell from the best:

| Corpus | Axes taken | Clusters | Homogeneity | Completeness | V | ARI |
|---|---|---|---|---|---|---|
| RdRp | 16 | 21 | 0.911 | 0.936 | 0.923 | 0.965 |
| VP2 | 11 | 19 | 0.917 | 0.946 | 0.931 | 0.985 |
| VP1 | 11 | 17 | 0.884 | 0.948 | 0.915 | 0.979 |

Stage 1's `today` arm scores 0.866, 0.848 and 0.777 by V-measure over its last four rounds, so the
gain is 0.06 to 0.14, nearly all of it completeness: 0.94 against 0.66–0.82. Each run costs about
fifteen first rounds — 25 minutes on RdRp, 25 on VP2 and 66 on VP1, single-threaded.

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
