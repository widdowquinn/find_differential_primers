# README.md - tests

This directory contains test code and data for the `diagnostic_primers` module.

Tests are divided conceptually into two kinds:

1. tests of `pdp.py` function from the command-line to confirm that all is working
2. tests for code correctness in the module

## 1. Tests of `pdp.py` command-line function

These tests comprise a series of 'known good' commands that can be run at the terminal. All the commands below should run without error. Due to the nature of the primer prediction tools used, it cannot be guaranteed that the predicted diagnostic primers will be the same on each run.

### 1a: `pdp.py process`

The `process` subcommand checks validity of the input file for:

* presence of the input sequence
* whether the input sequence is in multiple parts
* whether the input sequence contains ambiguity codons other than `N` (this are not tolerated by the primer prediction package `ePrimer3`)

If the `--validate` option is provided, `pdp.py` does no more than check for the above, and report to the terminal.

```
pdp.py process -v --validate test_data/testin.conf test_data/testprocess.conf
```

If the `--validate` option is not provided, then sequences that require stitching or ambiguity symbol replacement are modified, and the output written to file.

```
pdp.py process -v test_data/testin.conf test_data/testprocess.conf
```

The above command should create a new file `test_data/testprocess.conf`, and also create new files in the `test_input` subdirectory.

### 1b: `pdp.py prodigal`

The `prodigal` subcommand runs bacterial CDS prediction on the input sequences, generating one job per sequence and passing them to the requested scheduler. As the `multiprocessing` scheduler is expected always to be present, the examples below will use only this scheduler. To test the same command with an SGE-based system, use the `--scheduler SGE` option.

This subcommand should fail if the output directory already exists,

```
mkdir -p test_input/prodigal
pdp.py prodigal -v test_data/testprocess.conf test_data/testprodigal.conf
```

Using the `-f` or `--force` option should cause the Prodigal runs to go to completion, and the new config file to be written.

```
pdp.py prodigal -v test_data/testprocess.conf test_data/testprodigal.conf -f
```

### 1c: `pdp.py eprimer3`

The `eprimer3` subcommand runs EMBOSS ePrimer3 PCR primer pair prediction on the input sequences, generating one job per sequence and passing them to the requested scheduler. The commands below use the `multiprocessing` scheduler, but the `--scheduler SGE` option can be used to employ an SGE_like scheduler.

This subcommand should fail if the output directory already exists,

```
mkdir -p test_input/ePrimer3
pdp.py eprimer3 -v test_data/testprodigal.conf test_data/testeprimer3.conf
```

Using the `-f` or `--force` option should cause the ePrimer3 runs to go to completion, and a new config file to be written.

```
pdp.py eprimer3 -v test_data/testprodigal.conf test_data/testeprimer3.conf -f
```

### 1d: `pdp.py blastscreen`

The `blastscreen` subcommand screens the primer files associated with each input sequence against a BLAST-formatted nucleotide database. If a match is found to either the forward or reverse primer in a pair, the primer set is 'marked' as suspect, and removed from the primerset file (the internal oligo is ignored). The primer sets that pass the screen are written to a new file, and a new configuration file produced accordingly.

```
pdp.py blastscreen -v test_data/testeprimer3.conf test_data/testscreened.conf
```