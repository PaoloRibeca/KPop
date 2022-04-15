
# `KPop`: Unleash the full power of your *k*-mers!



`KPop` is implemented for the most part in [OCaml](https://ocaml.org), an industry-strength programming language that offers a number of advantages &mdash; amazing concision and symbolic power, static typing, incredible robustness and superior compiled speed. Due mostly to historical, prototyping reasons, a small part of `KPoP` is still in R, although we hope to evntually migrate everything to OCaml. All programs, both OCaml and R, are parallelised and will automatically use as many CPUs as are available on your machine, in order to speed up the wallclock execution time of your tasks as much as possible.

Depending on the problem at hand, `KPop` analysis can require a large amount of computational resources, memory in particular. That is a feature &mdash; i.e., a conscious design choice &mdash; rather than a bug. We recommend the use of a relatively large HPC node (at least 16 CPU cores and 256 GB of RAM) as a starting point for exploration.

## Table of contents

1. [Installation](#installation)
2. [Overview](#overview)
3. [Command line syntax](#command-line-syntax)
   - [`KPopCount`](#kpopcount)
   - [`KPopCountDB`](#kpopcountdb)
   - [`KPopTwist`](#kpoptwist)
   - [`KPopTwistDB`](#kpoptwistdb)
   - [`Parallel`](#parallel)
4. [Examples](#examples)
   - [Sequence classification](#sequence-classification)
     - [Classifier for simulated COVID-19 sequencing reads](#classifier-for-simulated-covid-19-sequencing-reads)
     - [Classifier for COVID-19 sequences (Hyena)](#classifier-for-covid-19-sequences-(Hyena))
   - [Pseudo-phylogenetic trees](#pseudo-phylogenetic-trees)

## Installation

There are several ways of installing the software on your machine:

> :warning: Note that the only operating system we support is Linux. :warning:

### Conda channel

### Pre-compiled binaries 

### Manual compilation

Alternatively, you can install `KPop` manually by compiling its sources. You will need an up-to-date distribution of the OCaml compiler and the [Dune package manager](https://github.com/ocaml/dune) for that. Both can be installed through [OPAM](https://opam.ocaml.org/), the official OCaml distribution system. Once you have a working OPAM distribution you will also have a working OCaml compiler, and Dune can be installed with the command
```
$ opam install dune
```
if it is not already present. Make sure that you install OCaml version 4.12 or later.

Then go to the directory into which you have downloaded the latest `KPop` sources, and type
```
$ . BUILD
```

That should generate all the executables you'll need (as of this writing, `Parallel`, `KPopCount`, `KPopCounterDB`, `kPopTwist`, `KPopTwist`, `KPopTwistDB`). Copy them to some favourite location in your PATH, for instance `~/.local/bin`.

## Overview

`KPop` comes as a number of different programs, each one

## Command line syntax

### `KPopCount`

This is the list of command line options available for the program `KPopCount`. You can visualise the list by typing
```bash
$ KPopCount -h
```
in your terminal. You will see a header containing information about the version:
```
This is the KPopCount program (version 0.3)
 (c) 2017-2022 Paolo Ribeca, <paolo.ribeca@gmail.com>
```
followed by detailed information. The general form(s) the command can be used is:
```
KPopCount -l|--label <output_vector_label> [OPTIONS]
```

Algorithmic parameters:
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-k`<br>`-K`<br>`--k-mer-size`<br>`--k-mer-length` | _&lt;k\_mer\_length&gt;_ |  k\-mer length \(must be positive and &lt;= 30\) | <ins>default=<mark>_12_</mark></ins> |
| `-m`<br>`-M`<br>`--max-results-size` | _&lt;positive\_integer&gt;_ |  maximum number of k\-mer signatures to be kept in memory at any given time\.<br>If more are present, the ones corresponding to the lowest cardinality will be removed from memory and printed out, and there will be repeated signatures in the output | <ins>default=<mark>_16777216_</mark></ins> |

Input/Output:
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-f`<br>`-F`<br>`--fasta` | _&lt;fasta\_file\_name&gt;_ |  FASTA input file containing sequences |  |
| `-s`<br>`-S`<br>`--single-end` | _&lt;fastq\_file\_name&gt;_ |  FASTQ input file containing single\-end sequencing reads |  |
| `-p`<br>`-P`<br>`--paired-end` | _&lt;fastq\_file\_name1&gt; &lt;fastq\_file\_name2&gt;_ |  FASTQ input files containing paired\-end sequencing reads |  |
| `-l`<br>`--label` | _&lt;output\_vector\_label&gt;_ |  label of the k\-mer vector in the output file | *(mandatory)* |
| `-o`<br>`--output` | _&lt;output\_file\_name&gt;_ |  name of generated output file | <ins>default=<mark>_&lt;stdout&gt;_</mark></ins> |

Miscellaneous:
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-v`<br>`--verbose` |  |  set verbose execution | <ins>default=<mark>_false_</mark></ins> |
| `-h`<br>`--help` |  |  print syntax and exit |  |

### `KPopCountDB`

This is the list of command line options available for the program `KPopCountDB`. You can visualise the list by typing
```bash
$ KPopCountDB -h
```
in your terminal. You will see a header containing information about the version:
```
This is the KPopCountDB program (version 0.27)
 (c) 2020-2022 Paolo Ribeca, <paolo.ribeca@gmail.com>
```
followed by detailed information. The general form(s) the command can be used is:
```
KPopCountDB [ACTIONS]
```

Actions \(executed delayed and in order of specification\):
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-e`<br>`-E`<br>`--empty` |  |  put an empty database into the register |  |
| `-i`<br>`-I`<br>`--input` | _&lt;binary\_file\_prefix&gt;_ |  load to the register the database present in the specified file  \(which must have extension \.KPopCounter\) |  |
| `-m`<br>`-M`<br>`--metadata`<br>`--add-metadata` | _&lt;metadata\_table\_file\_name&gt;_ |  add to the database present in the register metadata from the specified file |  |
| `-f`<br>`-F`<br>`--files`<br>`--add-files` | _&lt;k\-mer\_table\_file\_name&gt;\[','\.\.\.','&lt;k\-mer\_table\_file\_name&gt;\]_ |  add to the database present in the register k\-mers from the specified files |  |
| `--summary` |  |  print a summary of the database present in the register |  |
| `-l`<br>`-L`<br>`--labels`<br>`--selection-from-labels` | _&lt;vector\_label&gt;\[','\.\.\.','&lt;vector\_label&gt;\]_ |  put into the selection register the specified labels |  |
| `-r`<br>`-R`<br>`--regexps`<br>`--selection-from-regexps` | _&lt;metadata\_field&gt;'\~'&lt;regexp&gt;\[','\.\.\.','&lt;metadata\_field&gt;'\~'&lt;regexp&gt;\]_ |  put into the selection register the labels of the vectors whose metadata fields match the specified regexps\.<br>An empty metadata field matches the labels |  |
| `-a`<br>`-A`<br>`--add-sum-selection`<br>`--selection-add-sum` | _&lt;new\_vector\_label&gt;_ |  add to the database present in the register \(or replace if the new label exists\) a linear combination of the vectors whose labels are in the selection register |  |
| `-d`<br>`-D`<br>`--delete`<br>`--selection-remove` |  |  remove from the table the vectors whose labels are in the selection register |  |
| `-n`<br>`-N`<br>`--selection-negate` |  |  negate the labels that are present in the selection register |  |
| `-p`<br>`-P`<br>`--selection-print` |  |  print the labels that are present in the selection register |  |
| `-c`<br>`-C`<br>`--selection-clear` |  |  purges the selection register |  |
| `-s`<br>`-S`<br>`--selection-to-table-filter` |  |  filters out vectors whose labels are present in the selection register when writing the database as a tab\-separated file |  |
| `--table-emit-row-names` | _'true'&#124;'false'_ |  whether to emit row names for the database present in the register when writing it as a tab\-separated file | <ins>default=<mark>_true_</mark></ins> |
| `--table-emit-col-names` | _'true'&#124;'false'_ |  whether to emit column names for the database present in the register when writing it as a tab\-separated file | <ins>default=<mark>_true_</mark></ins> |
| `--table-emit-metadata` | _'true'&#124;'false'_ |  whether to emit metadata for the database present in the register when writing it as a tab\-separated file | <ins>default=<mark>_false_</mark></ins> |
| `--table-transpose` | _'true'&#124;'false'_ |  whether to transpose the database present in the register before writing it as a tab\-separated file  \(if 'true' : rows are vector names, columns are \(metadata and\) k\-mer names;   if 'false': rows are \(metadata and\) k\-mer names, columns are vector names\) | <ins>default=<mark>_false_</mark></ins> |
| `--table-threshold` | _&lt;non\_negative\_integer&gt;_ |  set to zero all counts that are less than this threshold before transforming and outputting them | <ins>default=<mark>_1_</mark></ins> |
| `--table-power` | _&lt;non\_negative\_float&gt;_ |  raise counts to this power before transforming and outputting them\.<br>A power of 0 when the 'pseudocount' method is used performs a logarithmic transformation | <ins>default=<mark>_1\._</mark></ins> |
| `--table-transform`<br>`--table-transformation` | _'none'&#124;'normalize'&#124;'pseudocount'&#124;'clr'_ |  transformation to apply to table elements before outputting them | <ins>default=<mark>_normalize_</mark></ins> |
| `--table-emit-zero-rows` | _'true'&#124;'false'_ |  whether to emit rows whose elements are all zero when writing the database as a tab\-separated file | <ins>default=<mark>_false_</mark></ins> |
| `--table-set-precision` | _&lt;positive\_integer&gt;_ |  set the number of precision digits to be used when outputting counts | <ins>default=<mark>_15_</mark></ins> |
| `-t`<br>`--table` | _&lt;file\_prefix&gt;_ |  write the database present in the register as a tab\-separated file  \(rows are k\-mer names, columns are vector names;   the file will be given extension \.KPopCounter\.txt\) |  |
| `-o`<br>`-O`<br>`--output` | _&lt;binary\_file\_prefix&gt;_ |  dump the database present in the register to the specified file  \(which will be given extension \.KPopCounter\) |  |

Miscellaneous \(executed immediately\):
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-T`<br>`--threads` | _&lt;computing\_threads&gt;_ |  number of concurrent computing threads to be spawned  \(default automatically detected from your configuration\) | <ins>default=<mark>_4_</mark></ins> |
| `-v`<br>`--verbose` |  |  set verbose execution | <ins>default=<mark>_false_</mark></ins> |
| `-h`<br>`--help` |  |  print syntax and exit |  |

### `KPopTwist`

This is the list of command line options available for the program `KPopTwist`. You can visualise the list by typing
```bash
$ KPopTwist -h
```
in your terminal. You will see a header containing information about the version:
```
This is the KPopTwist program (version 0.11)
 (c) 2022 Paolo Ribeca, <paolo.ribeca@gmail.com>
```
followed by detailed information. The general form(s) the command can be used is:
```
KPopTwist -i|--input <input_table_prefix> [OPTIONS]
```

Algorithmic parameters:
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-f`<br>`-F`<br>`-s`<br>`-S`<br>`--fraction`<br>`--sampling`<br>`--sampling-fraction` | _&lt;non\_negative\_float&gt;_ |  fraction of the rows to be considered and resampled before twisting | <ins>default=<mark>_1\._</mark></ins> |

Input/Output:
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-i`<br>`--input` | _&lt;input\_table\_prefix&gt;_ |  load the specified k\-mer database in the register and twist it\.<br>File extension is automatically determined  \(will be \.KPopCounter\)\.<br>The prefix is then re\-used for output  \(and the output file will be given prefix \.KPopTwisted\) | *(mandatory)* |
| `-T`<br>`--threads` | _&lt;computing\_threads&gt;_ |  number of concurrent computing threads to be spawned  \(default automatically detected from your configuration\) | <ins>default=<mark>_32_</mark></ins> |
| `-v`<br>`--verbose` |  |  set verbose execution | <ins>default=<mark>_false_</mark></ins> |
| `-h`<br>`--help` |  |  print syntax and exit |  |

### `KPopTwistDB`

This is the list of command line options available for the program `KPopTwistDB`. You can visualise the list by typing
```bash
$ KPopTwistDB -h
```
in your terminal. You will see a header containing information about the version:
```
This is the KPopTwistDB program (version 0.9)
 (c) 2022 Paolo Ribeca, <paolo.ribeca@gmail.com>
```
followed by detailed information. The general form(s) the command can be used is:
```
KPopTwistDB [ACTIONS]
```

Actions \(executed delayed and in order of specification\):
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-e`<br>`--empty` | _T&#124;t&#124;d_ |  load an empty twisted database into the specified register  \(T=twister; t=twisted; d=distance\) |  |
| `-i`<br>`--input` | _T&#124;t&#124;d &lt;binary\_file\_prefix&gt;_ |  load the specified binary database into the specified register  \(T=twister; t=twisted; d=distance\)\.<br>File extension is automatically determined depending on database type  \(will be: \.KPopTwister; \.KPopTwisted; or \.KPopDMatrix, respectively\) |  |
| `-I`<br>`--Input` | _T&#124;t&#124;d &lt;table\_file\_prefix&gt;_ |  load the specified tabular database\(s\) into the specified register  \(T=twister; t=twisted; d=distance\)\.<br>File extension is automatically determined depending on database type  \(will be: \.KPopTwister\.txt and \.KPopInertia\.txt; \.KPopTwisted\.txt;   or \.KPopDMatrix, respectively\) |  |
| `-f`<br>`-F`<br>`-k`<br>`-K`<br>`-a`<br>`-A`<br>`--add`<br>`--files`<br>`--kmers`<br>`--add-files`<br>`--add-kmers` | _&lt;k\-mer\_table\_file\_name&gt;\[','\.\.\.','&lt;k\-mer\_table\_file\_name&gt;\]_ |  twist k\-mers from the specified files through the transformation present in the twister register, and add the results to the database present in the twisted register |  |
| `--distance`<br>`--distance-function`<br>`--set-distance`<br>`--set-distance-function` | _'euclidean'&#124;'minkowski'_ |  set the function to be used when computing distances | <ins>default=<mark>_euclidean_</mark></ins> |
| `--distance-power`<br>`--set-distance-power` | _&lt;non\_negative\_float&gt;_ |  set the power to be used when computing \(Minkowski\) distances | <ins>default=<mark>_2\._</mark></ins> |
| `--metric`<br>`--metric-function`<br>`--set-metric`<br>`--set-metric-function` | _'flat'_ |  set the metric function to be used when computing distances | <ins>default=<mark>_flat_</mark></ins> |
| `-d`<br>`--distances`<br>`--compute-distances`<br>`--compute-distances-twisted` | _&lt;twisted\_binary\_file\_prefix&gt;_ |  compute distances between all the vectors present in the twisted register and all the vectors present in the specified twisted binary file  \(which must have extension \.KPopTwisted\) |  |
| `-o`<br>`--output` | _T&#124;t&#124;d &lt;binary\_file\_prefix&gt;_ |  dump the database present in the specified register  \(T=twister; t=twisted; d=distance\) to the specified binary file\.<br>File extension is automatically determined depending on database type  \(will be: \.KPopTwister; \.KPopTwisted; or \.KPopDMatrix, respectively\) |  |
| `--precision`<br>`--set-precision`<br>`--set-table-precision` | _&lt;positive\_integer&gt;_ |  set the number of precision digits to be used when outputting numbers | <ins>default=<mark>_15_</mark></ins> |
| `-O`<br>`--Output` | _T&#124;t&#124;d &lt;table\_file\_prefix&gt;_ |  dump the database present in the specified register  \(T=twister; t=twisted; d=distance\) to the specified tabular file\(s\)\.<br>File extension is automatically determined depending on database type  \(will be: \.KPopTwister\.txt and \.KPopInertia\.txt; \.KPopTwisted\.txt;   or \.KPopDMatrix, respectively\) |  |

Miscellaneous \(executed immediately\):
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-T`<br>`--threads` | _&lt;computing\_threads&gt;_ |  number of concurrent computing threads to be spawned  \(default automatically detected from your configuration\) | <ins>default=<mark>_4_</mark></ins> |
| `-v`<br>`--verbose` |  |  set verbose execution | <ins>default=<mark>_false_</mark></ins> |
| `-h`<br>`--help` |  |  print syntax and exit |  |

### `Parallel`

This is the list of command line options available for the program `KPopTwistDB`. You can visualise the list by typing
```bash
$ Parallel -h
```
in your terminal. You will see a header containing information about the version:
```
This is the Parallel program (version 0.3)
 (c) 2019-2022 Paolo Ribeca, <paolo.ribeca@gmail.com>
```
followed by detailed information. The general form(s) the command can be used is:
```
Parallel [OPTIONS] -- [COMMAND TO PARALLELIZE AND ITS OPTIONS]
```

Command to parallelize
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `--` |  |  consider all the subsequent parameters as the command to be executed in parallel\.<br>At least one command must be specified | *(mandatory)* |

Input/Output
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-l`<br>`--lines-per-block` | _&lt;positive\_integer&gt;_ |  number of lines to be processed per block | <ins>default=<mark>_10000_</mark></ins> |
| `-i`<br>`--input` | _&lt;input\_file&gt;_ |  name of input file | <ins>default=<mark>_stdin_</mark></ins> |
| `-o`<br>`--output` | _&lt;output\_file&gt;_ |  name of output file | <ins>default=<mark>_stdout_</mark></ins> |

Miscellaneous
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `-t`<br>`--threads` | _&lt;positive\_integer&gt;_ |  number of concurrent computing threads to be spawned  \(default automatically detected from your configuration\) | <ins>default=<mark>_32_</mark></ins> |
| `-v`<br>`--verbose` |  |  set verbose execution | <ins>default=<mark>_false_</mark></ins> |
| `-d`<br>`--debug` |  |  output debugging information | <ins>default=<mark>_false_</mark></ins> |
| `-h`<br>`--help` |  |  print syntax and exit |  |

## Examples

### Sequence classification

#### Classifier for simulated COVID-19 sequencing reads

As described in our [bioRxiv preprint](https://bioRxiv.org), we simulated 

In the following, we assume that the input files derived from the simulation have been organised into directories relative to your current location, and their placement reflects the set (training/test) and sequence cluster each sample belongs to.

So, we'll have two directories,
```
./Train
./Test
```
and each directory will contain subdirectories, one for each sequence cluster:
```
./Train/001
./Train/002
...
./Train/100
```
So for instance, directory `./Train/58` will contain input files for all the samples used as training data for sequence class 58, and command
```
$ ls ./Train/058
```
will return
```
05674_1.fastq  05686_2.fastq  05700_1.fastq  05712_2.fastq  05726_1.fastq  05738_2.fastq  05752_1.fastq  05764_2.fastq
05674_2.fastq  05688_1.fastq  05700_2.fastq  05714_1.fastq  05726_2.fastq  05740_1.fastq  05752_2.fastq  05766_1.fastq
05676_1.fastq  05688_2.fastq  05702_1.fastq  05714_2.fastq  05728_1.fastq  05740_2.fastq  05754_1.fastq  05766_2.fastq
05676_2.fastq  05690_1.fastq  05702_2.fastq  05716_1.fastq  05728_2.fastq  05742_1.fastq  05754_2.fastq  05768_1.fastq
05678_1.fastq  05690_2.fastq  05704_1.fastq  05716_2.fastq  05730_1.fastq  05742_2.fastq  05756_1.fastq  05768_2.fastq
05678_2.fastq  05692_1.fastq  05704_2.fastq  05718_1.fastq  05730_2.fastq  05744_1.fastq  05756_2.fastq  05770_1.fastq
05680_1.fastq  05692_2.fastq  05706_1.fastq  05718_2.fastq  05732_1.fastq  05744_2.fastq  05758_1.fastq  05770_2.fastq
05680_2.fastq  05694_1.fastq  05706_2.fastq  05720_1.fastq  05732_2.fastq  05746_1.fastq  05758_2.fastq  05772_1.fastq
05682_1.fastq  05694_2.fastq  05708_1.fastq  05720_2.fastq  05734_1.fastq  05746_2.fastq  05760_1.fastq  05772_2.fastq
05682_2.fastq  05696_1.fastq  05708_2.fastq  05722_1.fastq  05734_2.fastq  05748_1.fastq  05760_2.fastq  05774_1.fastq
05684_1.fastq  05696_2.fastq  05710_1.fastq  05722_2.fastq  05736_1.fastq  05748_2.fastq  05762_1.fastq  05774_2.fastq
05684_2.fastq  05698_1.fastq  05710_2.fastq  05724_1.fastq  05736_2.fastq  05750_1.fastq  05762_2.fastq  05776_1.fastq
05686_1.fastq  05698_2.fastq  05712_1.fastq  05724_2.fastq  05738_1.fastq  05750_2.fastq  05764_1.fastq  05776_2.fastq

```
A similar structure will have been put in place for test data under the `./Test` directory, with files equally split under `./Train` and `./Test` for cross-validation purposes (but different choices would be possible).

In order to analyse sequences in parallel fashion and decrease waiting times, we first prepare a short `bash` script named `process_one_class`, as follows:
```bash
#!/usr/bin/env bash

while read DIR; do
  CLASS=$(echo "$DIR" | awk '{l=split(gensub("[/]+$","",1),s,"/"); print s[l]}')
  echo "Processing class '${CLASS}'..." > /dev/stderr
  ls "$DIR"/*_1.fastq |
    Parallel --lines-per-block 1 -- awk '{l=split($0,s,"/"); system("KPopCount -k 10 -l "gensub("_1.fastq$","",1,s[l])" -p "$0" "gensub("_1.fastq$","_2.fastq",1))}' |
    KPopCountDB -f /dev/stdin -r "~." -a "$CLASS" -p -l "$CLASS" -n -p -d --summary --table-transform none -t /dev/stdout 2> /dev/null
done
```
The script takes as input a list of directories, each one on a single line
So, for instance,
```bash
$ echo Train/058 | ./process_one_class
```
would process all files present in directory `./Train/058`, generate the 10-mer spectra for each of them, combine them into an in-memory KPopCount database, generate their linear combination, discard the database, and output the combination to standard output. 

Note that the script is implicitly parallelised, in that each of the programs used will check for the number of available processors, and start an adequate number of computing threads to take full advantage of them.

At this point, in order to complete the "training" phase, we need to issue the two commands
```bash
$ ls -d Train/*/ | Parallel --lines-per-block 1 -- ./process_one_class | KPopCountDB -f /dev/stdin -o Classes
$ KPopTwist -i Classes
```

That will


#### Classifier for COVID-19 sequences (Hyena)

This is a rather more complex example, that showcases 

```
CHUNK_BLOCKS=125; BLOCK_SIZE=$(( 1024 * 1024 )); echo ../GISAID/2022_04_11/sequences.fasta | awk -v CHUNK_BLOCKS="$CHUNK_BLOCKS" -v BLOCK_SIZE="$BLOCK_SIZE" 'END{get_size="ls -l \047"$0"\047 | awk \047{print $5}\047" |& getline size; close(get_size); blocks=int((size+BLOCK_SIZE)/BLOCK_SIZE); chunks=int((blocks+CHUNK_BLOCKS)/CHUNK_BLOCKS); for (i=0;i<chunks;++i) print $0"\t"i*CHUNK_BLOCKS}' | Parallel -l 1 -t $(( $(nproc) * 3 )) -- awk -F '\t' -v BLOCK_SIZE="$BLOCK_SIZE" '{get_chunk="dd bs="BLOCK_SIZE" if=\047"$1"\047 skip="$2" 2> /dev/null"; overhang=""; line=""; while (get_chunk |& getline line) {if (line==""||line~"^>") break; overhang=overhang line"\n"} print $1"\t"($2*BLOCK_SIZE)+length(overhang); close(get_chunk)}' | awk -F '\t' '{if (NR>1) print $1"\t"old"\t"($2-old); old=$2}' | Parallel -l 1 -- awk -F '\t' -v BLOCK_SIZE="$BLOCK_SIZE" 'function remove_spaces(name,     s){split(name,s,"/"); return gensub("[ _]","","g",s[1])"/"s[2]"/"s[3]} BEGIN{nr=0; while (getline < "lineages.csv") {++nr; if (nr>1) {split($0,s,","); t[remove_spaces(gensub("^BHR/","Bahrain/",1,s[1]))]=s[2]}}} function output_sequence(){++counts[class]; print ">"name"\n"seq > (counts[class]%2==1?"Train":"Test")"/"offset"_"class".fasta"; return} {offset=$2; get_chunk="dd bs="BLOCK_SIZE" if=\047"$1"\047 iflag=\"skip_bytes,count_bytes\" skip="offset" count="$3" 2> /dev/null"; first=0; while (get_chunk |& getline) {if ($0~"^>") {if (first&&active) output_sequence(); first=1; active=0; split(substr($0,10),s,"[|]"); res=remove_spaces(s[1]); if (res in t) {active=1; name=res; class=t[res]; seq=""}} else {if (active) seq=seq $0}} if (active) output_sequence()} END{for (i in counts) print i"\t"counts[i]}' | awk -F '\t' '{counts[$1]+=$2} END{for (i in counts) print i"\t"counts[i]}' > STATS.txt

```

### Pseudo-phylogenetic trees


