
# KPop: Unleash the full power of your *k*-mers!

## Overview


## Command line syntax

### `KPopCountDB`

This is the list of commandline options available for `KPopCountDB`. You can visualise it by typing
```
KPopCountDB -h
```
in your terminal.
```
This is the KPopCountDB program (version 0.27)
 (c) 2020-2022 Paolo Ribeca, <paolo.ribeca@gmail.com>
Usage:
  KPopCountDB [ACTIONS]
 Actions (executed delayed and in order of specification):
  -e|-E|--empty
   | put an empty database into the register
  -i|-I|--input
    <binary_file_prefix>
   | load to the register the database present in the specified file
   |  (which must have extension .KPopCounter)
  -m|-M|--metadata|--add-metadata
    <metadata_table_file_name>
   | add to the database present in the register metadata from the specified file
```


## Examples

### Sequence classification

#### Classifier for simulated COVID reads



#### COVID-19 classifier (Hyena)

This is a rather more complex example, that showcases 

```


```
