# README.md - tests

This directory contains test code and data for the `diagnostic_primers` module and `pdp` command-line application.

Tests are divided conceptually into the following kinds:

- tests of `pdp` script function to confirm that a command-line call works as expected, linked into a chain emulating an example workflow using the `pdp` subcommands
- tests for code correctness in the module (unit tests, etc.)

The tests for taking datasets in order through a `pdp` analysis are found in test scripts that have filenames `test_NN_subcmd_<action>.py`, where `NN` indicates the order in which these tests should be run. These tests take a set of 15 input bacterial genome sequences (`tests/test_input/sequences`) through a series of primer design and classification steps.

Unit tests and other tests of code correctness are found in files with names of the form `test_<activity>.py`.

## Tests that emulate command-line calls

These tests attempt to implement commands that could be run at the terminal. In general, each test creates an `argparse.Namespace` object that emulates the parsed command-line being tested.

The subdirectories under `test_input`, `test_output` and `test_targets` for these command-line tests are prefixed with `pdp_`.

The tests must be run in order, as each generates output which is necessary for the next stage in the `pdp` design/analysis workflow. The numbering of tests is deliberate, as an attempt to ensure that tests are run in the correct order using `pytest`.

### `test_01_subcmd_config.py`

Tests of `pdp config` subcommands, checking that config file validation, tab to JSON and JSON to tab conversion, and sequence stitching/ambiguity symbol replacement work correctly. These tests also generate/modify input data for use with later test modules, to emulate taking real data through a primer design and classification process.

Input files are in `tests/test_inputs/pdp_config`, targets are in `tests/test_targets/pdp_config`.

The `test_fix_sequences()` method is equivalent to the following command:

```bash
pdp config --disable_tqdm -v \
    --outdir tests/test_output/pdp_config \
    --fix_sequences tests/test_output/pdp_config/seqfixed_conf.json \
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
    tests/test_output/pdp_filter/prodconf.json
```

and the `test_filter_prodigaligr_run()` method is equivalent to the command:

```bash
pdp filter -v --disable_tqdm --prodigaligr \
    --outdir tests/test_output/pdp_filter/prodigaligr \
    --suffix prodigaligr \
    tests/test_input/pdp_filter/seqfixed_conf.json \
    tests/test_output/pdp_filter/prodigrconf.json
```

Both commands use the same `tests/test_input/pdp_filter/seqfixed_conf.json` configuration file (which should be the same as the `tests/test_output/pdp_config/seqfixed_conf.json` output file from the earlier test `test_01_subcmd_config.py`, and describe stitched input sequences, with fixed ambiguity symbols), but place their program output in separate subdirectories under `tests/test_output/pdp_filter`. The new configuration files that they generate are also placed in `tests/test_output/pdp_filter`.

The output config files `tests/test_output/pdp_filter/prodconf.json` and `tests/test_output/pdp_filter/prodigrconf.json` are used as inputs for the next workflow stage.

### `test_03_subcmd_eprimer3.py`

Tests of `pdp eprimer3` subcommands. The tests generate output that is used for later workflow stage tests.

The `test_eprimer3_01_run()` and `test_eprimer3_02_force()` methods carry out a shortened analysis on a subset of three input files. This saves time on tests which would cause continuous integration tools like Travis CI to time out. The equivalent command-line is:

```bash
pdp eprimer3 -v --disable_tqdm \
    --outdir tests/test_output/pdp_eprimer3/subset \
    tests/test_input/pdp_eprimer3/subsetconf.json \
    tests/test_output/pdp_config/subsetep3conf.json
```

We still wish to use the complete set of `ePrimer3` results for the full input sequence set in later tests, so results were also generated manually using the commands:

```bash
pdp eprimer3 -v --filter \
    --numreturn 10 \
    --outdir tests/test_input/pdp_dedupe/prodigal \
    tests/test_input/pdp_eprimer3/prodconf.json \
    tests/test_input/pdp_dedupe/prod_ep3conf.json
```

for prodigal CDS-filtered regions, and

```bash
pdp eprimer3 -v --filter \
    --numreturn 10 \
    --outdir tests/test_input/pdp_dedupe/prodigaligr \
    tests/test_input/pdp_eprimer3/prodigrconf.json \
    tests/test_input/pdp_dedupe/prodigr_ep3conf.json
```

for intergenic regions, respectively (note the `--filter` option is used to enforce designing primers to the filtered sequence sets). This places two new post-primer design configuration files at `tests/test_input/pdp_dedupe/prod_ep3conf.json` and `tests/test_input/pdp_dedupe/prodigr_ep3conf.json`, with the corresponding `ePrimer3` output at `tests/test_input/pdp_dedupe/prodigal` and `tests/test_input/pdp_dedupe/prodigaligr`.

### `test_04_subcmd_dedupe.py`

Tests of `pdp dedupe` subcommands. The tests generate output used in later `pdp` workflow tests.

The two test methods `test_dedupe_prodigal_run()` and `test_dedupe_prodigaligr_run()` deduplicate primers designed to prodigal-predicted CDS and their corresponding intergenic regions, respectively. They are equivalent to the commands:

```bash
pdp dedupe -v --disable_tqdm \
    --dedupedir tests/test_output/pdp_dedupe/prodigal \
    tests/test_input/pdp_dedupe/prod_ep3conf.json \
    tests/test_output/pdp_dedupe/dedupe_prod.json
```

and

```bash
pdp dedupe -v --disable_tqdm \
    --dedupedir tests/test_output/pdp_dedupe/prodigaligr \
    tests/test_input/pdp_dedupe/prodigr_ep3conf.json \
    tests/test_output/pdp_dedupe/dedupe_prodigr.json  
```

The `prodigal` and `prodigaligr` directories produced from deduplication, and the output configuration files `dedupe_prod.json` and `dedupe_prodigr.json` are the same as those used as input for the next set of tests in `tests/test_input/test_input/pdp_blastscreen`

### `test_05_subcmd_blastscreen.py`

Tests of `pdp blastscreen` subcommands. The tests generate output used in later `pdp` workflow tests.

The two test methods `test_blastscreen_prodigal_01_run()` and `test_blastscreen_prodigaligr_run()` screen the primer sets defined in `dedupe_prod.json` and `dedupe_prodigr.json` respectively against a `BLAST+` nucleotide database of *Escherichia coli* sequences in `tests/test_input/pdp_blastscreen/blastdb`. The results are written to `test_output/pdp_blastscreen`, and the two configuration files produced (`screened_prod.json` and `screened_prodigr.json`) are used as input to the next workflow stage.

The equivalent commands are:

```bash
pdp blastscreen -v -f --disable_tqdm \
    --db=tests/test_input/pdp_blastscreen/blastdb/e_coli_screen.fna \
    --outdir=tests/test_output/pdp_blastscreen/prodigal \
    tests/test_input/pdp_blastscreen/dedupe_prod.json \
    tests/test_output/pdp_blastscreen/screened_prod.json
```

and

```bash
pdp blastscreen -v -f --disable_tqdm \
    --db=tests/test_input/pdp_blastscreen/blastdb/e_coli_screen.fna \
    --outdir=tests/test_output/pdp_blastscreen/prodigaligr \
    tests/test_input/pdp_blastscreen/dedupe_prodigr.json \
    tests/test_output/pdp_blastscreen/screened_prodigr.json
```

### `test_06_subcmd_primersearch.py`

Tests of `pdp primersearch` subcommands. The tests generate output used in later `pdp` workflow tests.

The two test methods `test_primersearch_prodigal_run()` and `test_primersearch_prodigaligr_run()` conduct primersearch cross-hybridisation tests on the primer sets defined in `screened_prod.json` and `screened_prodigr.json` against the sequences defined in the same files (and stored in `tests/test_input/sequences`). The results are written to `test_output/pdp_primersearch` and the two configuration file outputs (`primersearch_prod.json` and `primersearch_prodigr.json`) are used as input to the next workflow stage tests.

The equivalent commands are:

```bash
pdp primersearch -v --disable_tqdm \
    --outdir=tests/test_output/pdp_primersearch/prodigal \
    tests/test_input/pdp_primersearch/screened_prod.json \
    tests/test_output/pdp_primersearch/primersearch_prod.json
```

```bash
pdp primersearch -v --disable_tqdm \
    --outdir=tests/test_output/pdp_primersearch/prodigaligr \
    tests/test_input/pdp_primersearch/screened_prodigr.json \
    tests/test_output/pdp_primersearch/primersearch_prodigr.json
```

### `test_07_subcmd_classify.py`

Tests of `pdp primersearch` subcommands. The tests generate output used in later `pdp` workflow tests.

The two test methods `test_classify_prodigal_run()` and `test_classify_prodigaligr_run()` classify the outputs described in `primersearch_prod.json` and `primersearch_prodigr.json`, respectively. The outputs (diagnostic sequences, `JSON` files, and summary information) are written to `tests/test_output/pdp_classify/prodigal` and `tests/test_output/pdp_classify/prodigaligr`, respectively. The sequence outputs are used as input to the next workflow stage tests.

The equivalent commands are:

```bash
pdp classify -v -f --disable_tqdm \
    tests/test_input/pdp_classify/primersearch_prodigal.json \
    tests/test_output/pdp_classify/prodigal
```

```bash
pdp classify -v -f --disable_tqdm \
    tests/test_input/pdp_classify/primersearch_prodigaligr.json \
    tests/test_output/pdp_classify/prodigaligr
```

### `test_08_subcmd_extract.py`

Tests of `pdp extract` subcommands. The tests generate output that is not used for later workflow tests.

The two methods `test_extract_prodigal_run()` and `test_extract_progidal_align()` extract the amplicons for primers defined in `primersearch_prod.json`; the first method extracts sequences, and the second extracts and aligns them (using `MAFFT`), for the primer sets defined in `tests/test_output/pdp_classify/prodigal/Pectobacterium_primers.json`. The equivalent command-lines are:

```bash
pdp extract -v --disable_tqdm -f \
    --noalign \
    tests/test_input/pdp_extract/primersearch_prod.json \
    tests/test_output/pdp_classify/prodigal/Pectobacterium_primers.json \
    tests/test_output/pdp_extract/prodigal/noalign
```

and

```bash
pdp extract -v --disable_tqdm -f \
    tests/test_input/pdp_extract/primersearch_prod.json \
    tests/test_output/pdp_classify/prodigal/Pectobacterium_primers.json \
    tests/test_output/pdp_extract/prodigal/align
```

Similarly, the two methods `test_extract_prodigaligr_run()` and `test_extract_progidaligr_align()` extract the amplicons for primers defined in `primersearch_prodigr.json`; the first method extracts sequences, and the second extracts and aligns them (using `MAFFT`), for the primer sets defined in `tests/test_output/pdp_classify/prodigaligr/Pectobacterium_primers.json`. The equivalent command-lines are:

```bash
pdp extract -v --disable_tqdm -f \
    --noalign \
    tests/test_input/pdp_extract/primersearch_prodigr.json \
    tests/test_output/pdp_classify/prodigaligr/Pectobacterium_primers.json \
    tests/test_output/pdp_extract/prodigaligr/noalign
```

and

```bash
pdp extract -v --disable_tqdm -f \
    tests/test_input/pdp_extract/primersearch_prodigr.json \
    tests/test_output/pdp_classify/prodigaligr/Pectobacterium_primers.json \
    tests/test_output/pdp_extract/prodigaligr/align
```
