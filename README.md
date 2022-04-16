
# `KPop`: Unleash the full power of your *k*-mers!



`KPop` is implemented for the most part in [OCaml](https://ocaml.org), an industry-strength programming language that offers a number of advantages &mdash; amazing concision and symbolic power, static typing, incredible robustness and superior compiled speed. Due mostly to historical, prototyping reasons, a small part of `KPoP` is still in R, although we hope to evntually migrate everything to OCaml. All programs, both OCaml and R, are parallelised and will automatically use as many CPUs as are available on your machine, in order to speed up the wallclock execution time of your tasks as much as possible.

Depending on the problem at hand, `KPop` analysis can require a large amount of computational resources, memory in particular. That is a feature &mdash; i.e., a conscious design choice &mdash; rather than a bug. We recommend the use of a relatively large HPC node (at least 16 CPU cores and 256 GB of RAM) as a starting point for exploration.

## Table of contents

[1. Installation](#1-installation)<br>
[2. Overview](#2-overview)<br>
[3. Command line syntax](#3-command-line-syntax)<br>
&emsp; [3.1. `KPopCount`](#31-kpopcount)<br>
&emsp; [3.2. `KPopCountDB`](#32-kpopcountdb)<br>
&emsp; [3.3. `KPopTwist`](#33-kpoptwist)<br>
&emsp; [3.4. `KPopTwistDB`](#34-kpoptwistdb)<br>
&emsp; [3.5. `Parallel`](#35-parallel)<br>
[4. Examples](#4-examples)<br>
&emsp; [4.1. Sequence classification](#41-sequence-classification)<br>
&emsp; &emsp; [4.1.1. Classifier for simulated COVID-19 sequencing reads](#411-classifier-for-simulated-covid-19-sequencing-reads)<br>
&emsp; &emsp; [4.1.2. Classifier for COVID-19 sequences (Hyena)](#412-classifier-for-covid-19-sequences-(Hyena))<br>
&emsp; [4.2. Pseudo-phylogenetic trees](#42-pseudo-phylogenetic-trees)<br>

## 1. Installation

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

## 2. Overview

`KPop` comes as a number of different programs, each one

## 3. Command line syntax

### 3.1. `KPopCount`

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

### 3.2. `KPopCountDB`

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

### 3.3. `KPopTwist`

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

### 3.4. `KPopTwistDB`

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

### 3.5. `Parallel`

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

## 4. Examples

### 4.1. Sequence classification

#### 4.1.1. Classifier for simulated COVID-19 sequencing reads

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

At this point, in order to perform the "training" phase, we need to issue the two commands
```bash
$ ls -d Train/*/ | Parallel --lines-per-block 1 -- ./process_one_class | KPopCountDB -f /dev/stdin -o Classes
$ KPopTwist -i Classes
```
The first command will generate 

The second command will twist



At this point, the command
```bash
ls Test/*/*_1.fastq | Parallel --lines-per-block 1 -- awk '{l=split($0,s,"/"); system("KPopCount -k 10 -l "gensub("_1.fastq$","",1,s[l])" -p "$0" "gensub("_1.fastq$","_2.fastq",1))}' | KPopTwistDB -i T Classes -f /dev/stdin -o t Test -v
```
will generate separate _k_-mer spectra for each file in the test set, and twist them according to the "classifying" transformation stored in file `Classes.KPopTwister`. The results will be stored in an additional file, `Test.KPopTwisted`.

All the files generated so far are binary &mdash; they 

Should you be curious to see what the content looks like in human-readable format, you can always use `KPopTwistDB` to convert the file (or any other file generated by `KPopTwist` or `KPopTwistDB`) to a tab-separate textual table as those you can import into, or export from, R. For instance, the command
```bash
$ KPopTwistDB -i t Test -O t Test
```
will load a binary file (hence option `-i`) of the "twisted" type (hence option `-i t`) with prefix `Test` (hence option `-i t Test`), i.e., according to the automatic naming conventions enforced by `KPop`, it will load file `Test.KPopTwisted` into the "twisted" register of `KPopTwistDB`. After that, the command will output the content of the "twisted" register in tabular form (hence option `-O t`) to a file with prefix `Test`, i.e., according to `KPop`'s automatic naming conventions, to file `Test.KPopTwisted.txt`. You can see how such naming conventions free you from thinking about files extensions, and avoid name clashes between files having the same prefix but a different content.

The file `Test.KPopTwisted.txt` will contain a header
```
""    "Dim1"    "Dim2"    "Dim3"    "Dim4"    "Dim5"    "Dim6"    "Dim7"    "Dim8"    "Dim9"    "Dim10"    "Dim11"    "Dim12"    "Dim13"    "Dim14"    "Dim15"    "Dim16"    "Dim17"    "Dim18"    "Dim19"    "Dim20"    "Dim21"    "Dim22"    "Dim23"    "Dim24"    "Dim25"    "Dim26"    "Dim27"    "Dim28"    "Dim29"    "Dim30"    "Dim31"    "Dim32"    "Dim33"    "Dim34"    "Dim35"    "Dim36"    "Dim37"    "Dim38"    "Dim39"    "Dim40"    "Dim41"    "Dim42"    "Dim43"    "Dim44"    "Dim45"    "Dim46"    "Dim47"    "Dim48"    "Dim49"    "Dim50"    "Dim51"    "Dim52"    "Dim53"    "Dim54"    "Dim55"    "Dim56"    "Dim57"    "Dim58"    "Dim59"    "Dim60"    "Dim61"    "Dim62"    "Dim63"    "Dim64"    "Dim65"    "Dim66"    "Dim67"    "Dim68"    "Dim69"    "Dim70"    "Dim71"    "Dim72"    "Dim73"    "Dim74"    "Dim75"    "Dim76"    "Dim77"    "Dim78"    "Dim79"    "Dim80"    "Dim81"    "Dim82"    "Dim83"    "Dim84"    "Dim85"    "Dim86"    "Dim87"    "Dim88"    "Dim89"    "Dim90"    "Dim91"    "Dim92"    "Dim93"    "Dim94"    "Dim95"    "Dim96"    "Dim97"    "Dim98"    "Dim99"
```
followed by rows such as
```
"00002"    0.208303753652877    -0.204365699709935    0.474116247298613    -0.329722324765537    -0.385408051653039    0.0787792685599588    -0.851349902900387    0.30390279629666    -0.570446452179573    0.0852710206894159    -0.14790774401408    0.166607494778671    -0.0546181793656805    0.0476882959310556    0.0770788371437616    -0.0444686466010732    0.444721244127617    -0.673498184325522    0.362464609963262    -0.23454642728241    -1.34485439984727    0.658186285078834    -2.71800078600714    -1.27972655963794    2.27870884130667    1.24896517540089    0.561025490167653    0.450790128719747    -1.26132001288061    -0.250979150371771    -0.835539884518506    -0.443631637558216    -0.423180456184031    0.430020246411259    0.46866195496528    0.387402309339635    0.612014444008186    -0.834351696527018    -0.307770640074338    -0.727900404893945    0.142020888265305    0.0135700796514422    0.00172479916704904    -0.0386401755418014    0.273868140313799    0.118584665052629    -0.323070417267769    -0.315861011471072    -0.28873815035643    0.485502657444407    0.0980163084187126    -0.241195808461705    -0.0902438356160271    -0.348076849358953    0.149833659280451    0.101820897520809    0.436106979627955    -0.278760531142856    0.07642997624802    0.0734403482595658    -0.0976549111054914    -0.0662557825632169    0.130576096666569    0.19461712397724    0.131184003676229    -0.0274033695221185    0.0189332375778051    -0.0434489000554688    0.194493500351839    -0.0330209932696197    -0.0582006902983816    -0.015835571654948    0.0599251621030898    -0.263070425755192    0.163512702555629    0.0841590618810513    0.0296529536401996    0.0069817958564591    -0.0415246013366838    -0.0429221341747157    0.0686729171428239    -0.0526129530232641    0.12576164944361    0.0894051484230737    0.129241266260151    0.247597048821563    -0.336945850645509    -0.0319494824195588    0.0551447733877536    -0.0345965874113023    0.0448018576557607    -0.0463083593768541    -0.0246429930707945    -0.010635515038952    0.0414757895743002    0.037329550406034    0.130848249442206    0.0476197328655712    0.0871447416286987
```
expressing the coordinates of twisted sequences or samples (in this case, sample `00002`) in twisted space.

At this point we have several files in our directory, namely
```
```

What is left to do in order to classify the sequences in our test set is to compute the distance in twisted space between each sequence (rows of file `Test.KPopTwisted`) and the representative of each equivalence class of the training set (rows of file `Classes.KPopTwisted`); the classification for each given sequence will be the closest equivalence class (provided that the closest and second closest match are separated by some reasonable margin). Computing all pairwise distances in twisted space between sequences and classes is accomplished by the command
```bash
$ KPopTwistDB -i t Test -d Classes -O d Test-vs-Classes
``` 

#### 4.1.2. Classifier for COVID-19 sequences (Hyena)

This is a rather more complex example, that showcases 

```
CHUNK_BLOCKS=125; BLOCK_SIZE=$(( 1024 * 1024 )); echo ../GISAID/2022_04_11/sequences.fasta | awk -v CHUNK_BLOCKS="$CHUNK_BLOCKS" -v BLOCK_SIZE="$BLOCK_SIZE" 'END{get_size="ls -l \047"$0"\047 | awk \047{print $5}\047" |& getline size; close(get_size); blocks=int((size+BLOCK_SIZE)/BLOCK_SIZE); chunks=int((blocks+CHUNK_BLOCKS)/CHUNK_BLOCKS); for (i=0;i<chunks;++i) print $0"\t"i*CHUNK_BLOCKS}' | Parallel -l 1 -t $(( $(nproc) * 3 )) -- awk -F '\t' -v BLOCK_SIZE="$BLOCK_SIZE" '{get_chunk="dd bs="BLOCK_SIZE" if=\047"$1"\047 skip="$2" 2> /dev/null"; overhang=""; line=""; while (get_chunk |& getline line) {if (line==""||line~"^>") break; overhang=overhang line"\n"} print $1"\t"($2*BLOCK_SIZE)+length(overhang); close(get_chunk)}' | awk -F '\t' '{if (NR>1) print $1"\t"old"\t"($2-old); old=$2}' | Parallel -l 1 -- awk -F '\t' -v BLOCK_SIZE="$BLOCK_SIZE" 'function remove_spaces(name,     s){split(name,s,"/"); return gensub("[ _]","","g",s[1])"/"s[2]"/"s[3]} BEGIN{nr=0; while (getline < "lineages.csv") {++nr; if (nr>1) {split($0,s,","); t[remove_spaces(gensub("^BHR/","Bahrain/",1,s[1]))]=s[2]}}} function output_sequence(){++counts[class]; print ">"name"\n"seq > (counts[class]%2==1?"Train":"Test")"/"offset"_"class".fasta"; return} {offset=$2; get_chunk="dd bs="BLOCK_SIZE" if=\047"$1"\047 iflag=\"skip_bytes,count_bytes\" skip="offset" count="$3" 2> /dev/null"; first=0; while (get_chunk |& getline) {if ($0~"^>") {if (first&&active) output_sequence(); first=1; active=0; split(substr($0,10),s,"[|]"); res=remove_spaces(s[1]); if (res in t) {active=1; name=res; class=t[res]; seq=""}} else {if (active) seq=seq $0}} if (active) output_sequence()} END{for (i in counts) print i"\t"counts[i]}' | awk -F '\t' '{counts[$1]+=$2} END{for (i in counts) print i"\t"counts[i]}' > STATS.txt

```


### 4.2. Pseudo-phylogenetic trees


