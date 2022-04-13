
# KPop: Unleash the full power of your *k*-mers!

## Overview


## Command line syntax

### `KPopCount`

This is the list of command line options available for the program `KPopCount`. You can visualise the list by typing
```
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
```
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

### `KPopTwistDB`

This is the list of command line options available for the program `KPopTwistDB`. You can visualise the list by typing
```
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

## Examples

### Sequence classification

#### Classifier for simulated COVID reads



#### COVID-19 classifier (Hyena)

This is a rather more complex example, that showcases 

```


```
