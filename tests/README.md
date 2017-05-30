# README.md - tests

This directory contains test code and data for the `diagnostic_primers` module.

Tests are divided conceptually into two kinds:

1. tests of `pdp.py` function from the command-line to confirm that all is working
2. tests for code correctness in the module

## 1. Tests of `pdp.py` command-line function

These tests comprise a series of 'known good' commands that can be run at the terminal. All the commands below should run without error. Due to the nature of the primer prediction tools used, it cannot be guaranteed that the predicted diagnostic primers will be the same on each run.

### 1a: `pdp.py config`

The `config` subcommand provides for conversion and validation of configuration files.

```bash
pdp.py config <infile> --validate
pdp.py config <infile> --to_json <outfile>
pdp.py config <infile> --fix_sequences <outfile>
```

Validation will check either `.tab` or `.json` config input for

* presence of the input sequence
* whether the input sequence is in multiple parts
* whether the input sequence contains ambiguity codons other than `N` (this are not tolerated by the primer prediction package `ePrimer3`)

If the `--validate` option is provided, `pdp.py` does no more than check for the above, and report to the terminal. No output is generated

```
pdp.py config tests/test_input/config/testin.conf -v --validate  # will fail, too many columns
pdp.py config tests/test_input/config/testconf.tab -v --validate
pdp.py config tests/test_input/config/testconf.json -v --validate
```

If the `--to_json` option is provided, `pdp.py` converts input from `.tab` to `.json` format, and writes to the specified file, but **does not fix sequence problems**.

```
pdp.py config tests/test_input/config/testconf.tab -v --to_json tests/test_output/config/tab_converted_conf.json
```

If the `--fix_sequences` option is provided, `pdp.py` fixes any sequence problems, converts input from `.tab` to `.json` format (if required), and writes to the specified file.

```
pdp.py config tests/test_input/config/testconf.tab -v --fix_sequences tests/test_output/config/sequence_fix_conf.json
```


### 1b: `pdp.py prodigal`

The `prodigal` subcommand runs bacterial CDS prediction on the input sequences, generating one job per sequence and passing them to the requested scheduler. As the `multiprocessing` scheduler is expected always to be present, the examples below will use only this scheduler. To test the same command with an SGE-based system, use the `--scheduler SGE` option.

This subcommand should fail if the output directory already exists,

```
pdp.py prodigal -v tests/test_input/config/testconf.json tests/test_output/config/testprodigalconf.json \
                --outdir tests/test_output/prodigal
```

Using the `-f` or `--force` option should cause the Prodigal runs to go to completion, and the new config file to be written.

```
pdp.py prodigal -v tests/test_input/config/testconf.json tests/test_output/config/testprodigalconf.json \
                --outdir tests/test_output/prodigal -f
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