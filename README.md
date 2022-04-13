
# KPop: Unleash the full power of your *k*-mers!

## Overview


## Command line syntax

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
Actions (executed delayed and in order of specification):
| Option | Argument(s) | Default | Effect |
|-|-|-|-|
|`-e`<br>`-E`<br>`--empty`|||put an empty database into the register|
|`-i`<br>`-I`<br>`--input`|_<binary_file_prefix>_<br>&nbsp;<br>&nbsp;||load to the register the database present in the specified file<br> (which must have extension .KPopCounter)<br>&nbsp;|

## Examples

### Sequence classification

#### Classifier for simulated COVID reads



#### COVID-19 classifier (Hyena)

This is a rather more complex example, that showcases 

```


```
