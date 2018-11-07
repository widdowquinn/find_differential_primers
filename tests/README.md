# README.md - tests

This directory contains test code and data for the `diagnostic_primers` module and `pdp` command-line application.

Tests are divided conceptually into the following kinds:

- tests of `pdp` script function to confirm that a command-line call works as expected
- tests for code correctness in the module (unit tests, etc.)

The tests for taking datasets in order through a `pdp` analysis are found in test scripts that have filenames `test_NN_subcmd_<action>.py`, where `NN` indicates the order in which these tests should be run. These tests take a set of 15 input bacterial genome sequences (`tests/test_input/sequences`) through a series of primer design and classification steps.

Unit tests and other tests of code correctness are found in files with names of the form `test_<activity>.py`.

## Tests that emulate command-line calls

These tests attempt to implement commands that could be run at the terminal. In general, each test creates an `argparse.Namespace` object that emulates the parsed command-line being tested.

The subdirectories under `test_input`, `test_output` and `test_targets` for these command-line tests are prefixed with `pdp_`.

The tests must be run in order, as each generates output which is necessary for the next stage in the `pdp` design/analysis workflow. The numbering of tests is deliberate, as an attempt to ensure that tests are run in the correct order.

### `test_01_subcmd_config.py`

Tests of `pdp config` subcommands, checking that config file validation, tab to JSON and JSON to tab conversion, and sequence stitching/ambiguity symbol replacement work correctly. These tests also generate/modify input data for use with later test modules, to emulate taking real data through a primer design and classification process.

Input files are in `tests/test_inputs/pdp_config`, targets are in `tests/test_targets/pdp_config`.

The `test_fix_sequences()` method is equivalent to the following command:

```bash
pdp config --disable_tqdm -v \
    --fix_sequences tests/test_input/config/seqfixed_conf.json \
    tests/test_input/config/testconf.json
```

which stitches and replaces ambiguity symbols in the input sequences referred to by the configuration file `tests/test_input/config/testconf.json`. The output configuration file `tests/test_output/pdp_config/seqfixed_conf.json` should be the same file as `tests/test_input/pdp_filter/seqfixed_conf.json`.

### `test_02_subcmd_filter.py`

Tests of `pdp filter` subcommands. Two tests generate output that is used for later workflow stage tests.

Input files are in `tests/test_inputs/pdp_filter`, targets are in `tests/test_targets/pdp_filter`.

The `test_filter_prodigal_run()` method is equivalent to the following command:

```bash
pdp filter -v --disable_tqdm --prodigal \
    --outdir tests/test_output/pdp_filter/prodigal \
    --suffix prodigal \
    tests/test_input/pdp_filter/testreducedep3conf.json \
    tests/test_output/pdp_config/prodconf.json
```

and the `test_filter_prodigaligr_run()` method is equivalent to the command:

```bash
pdp filter -v --disable_tqdm --prodigaligr \
    --outdir tests/test_output/pdp_filter/prodigaligr \
    --suffix prodigaligr \
    tests/test_input/pdp_filter/seqfixed_conf.json \
    tests/test_output/pdp_config/prodigrconf.json
```

Both commands use the same `tests/test_input/pdp_filter/seqfixed_conf.json` configuration file (which should be the same as the `tests/test_output/pdp_config/seqfixed_conf.json` output file from the earlier test `test_01_subcmd_config.py`, and describe stitched input sequences, with fixed ambiguity symbols), but place their program output in separate subdirectories under `tests/test_output/pdp_filter`. The new configuration files that they generate are placed in `tests/test_output/pdp_config`.

The output config files `tests/test_output/pdp_config/prodconf.json` and `tests/test_output/pdp_config/prodigrconf.json` are used as inputs for the next workflow stage.

### `test_03_subcmd_eprimer3.py`

Tests of `pdp eprimer3` subcommands. The tests generate output that is used for later workflow stage tests.

The `test_eprimer3_01_run()` and `test_eprimer3_02_force()` methods carry out a shortened analysis on a subset of three input files. This saves time on tests which would cause continuous integration tools like Travis CI to time out. The equivalent command-line is:

```bash
bash eprimer3 -v --disable_tqdm \
    --outdir /tests/test_output/pdp_eprimer3/subset \
    tests/test_input/pdp_eprimer3/subsetconf.json \
    tests/test_output/pdp_config/subsetep3conf.json
```

We still wish to use the complete set of `ePrimer3` results for the full input sequence set in later tests, so the results were generated manually using the commands:

```bash
bash eprimer3 -v \
    --outdir /tests/test_input/pdp_dedupe/prodigal \
    tests/test_input/pdp_eprimer3/prodconf.json \
    tests/test_input/pdp_dedupe/prod_ep3conf.json
```

for prodigal CDS-filtered regions, and

```bash
bash eprimer3 -v \
    --outdir /tests/test_input/pdp_dedupe/prodigaligr \
    tests/test_input/pdp_eprimer3/prodigrconf.json \
    tests/test_input/pdp_dedupe/prodigr_ep3conf.json
```

for intergenic regions, respectively.