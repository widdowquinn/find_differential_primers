# README.md (`pdp`/`diagnostic_primers`)

<!-- TOC -->

- [README.md (`pdp`/`diagnostic_primers`)](#readmemd-pdpdiagnostic_primers)
    - [NOTE FOR USERS<a id="usernote"></a>](#note-for-usersa-idusernotea)
    - [NOTE FOR DEVELOPERS<a id="devnotenote"></a>](#note-for-developersa-iddevnotenotea)
    - [Overview<a id="overview"></a>](#overviewa-idoverviewa)
        - [Third-party Packages](#third-party-packages)
        - [Recent changes](#recent-changes)
    - [Walkthrough <a id="walkthrough"></a>](#walkthrough-a-idwalkthrougha)
        - [1. Producing and validating the config file](#1-producing-and-validating-the-config-file)
        - [2. Fix sequences for analysis](#2-fix-sequences-for-analysis)
        - [3. Defining CDS features on each genome (optional)](#3-defining-cds-features-on-each-genome-optional)
        - [4. Design primers to each genome in bulk](#4-design-primers-to-each-genome-in-bulk)
        - [5. Deduplicate primer sets (optional)](#5-deduplicate-primer-sets-optional)
        - [6. Screen primers against `BLASTN` database (optional)](#6-screen-primers-against-blastn-database-optional)
        - [7. Test primers against input sequences for crosshybridisation with `primersearch`](#7-test-primers-against-input-sequences-for-crosshybridisation-with-primersearch)
        - [8. Classify the primers by diagnostic capability with `classify`](#8-classify-the-primers-by-diagnostic-capability-with-classify)
    - [Usage<a id="usage"></a>](#usagea-idusagea)
        - [`pdp config`<a id="config"></a>](#pdp-configa-idconfiga)
            - [Converting between `.tab` and `JSON` format config files](#converting-between-tab-and-json-format-config-files)
            - [Validate a config file](#validate-a-config-file)
            - [Repair input sequences](#repair-input-sequences)
        - [`pdp prodigal`<a id="prodigal"></a>](#pdp-prodigala-idprodigala)
            - [Default feature prediction](#default-feature-prediction)
            - [Specify location to write `prodigal` predictions](#specify-location-to-write-prodigal-predictions)
            - [Specify the location of the `prodigal` executable](#specify-the-location-of-the-prodigal-executable)
        - [`pdp eprimer3`<a id="eprimer3"></a>](#pdp-eprimer3a-ideprimer3a)
            - [Default primer prediction](#default-primer-prediction)
            - [Change the number of predicted primers per input sequence](#change-the-number-of-predicted-primers-per-input-sequence)
            - [Change primer design parameters](#change-primer-design-parameters)
            - [Specify location to write primer prediction output](#specify-location-to-write-primer-prediction-output)
            - [Specify the location of the `eprimer3` executable](#specify-the-location-of-the-eprimer3-executable)
        - [`pdp blastscreen`<a id="blastscreen"></a>](#pdp-blastscreena-idblastscreena)
            - [Basic screening](#basic-screening)
            - [Controlling sensitivity](#controlling-sensitivity)
            - [Specify the location of the `BLAST+` executable](#specify-the-location-of-the-blast-executable)
            - [Control the number of threads used](#control-the-number-of-threads-used)
            - [Use the SGE/OGE scheduler](#use-the-sgeoge-scheduler)
        - [`pdp primersearch`<a id="primersearch"></a>](#pdp-primersearcha-idprimersearcha)
            - [Basic cross-hybridisation](#basic-cross-hybridisation)
            - [Controlling sensitivity](#controlling-sensitivity-1)
            - [Specify the location of the `primersearch` executable](#specify-the-location-of-the-primersearch-executable)
            - [Control the number of threads used](#control-the-number-of-threads-used-1)
            - [Use the SGE/OGE scheduler](#use-the-sgeoge-scheduler-1)
        - [`pdp classify`<a id="classify"></a>](#pdp-classifya-idclassifya)
            - [Basic classification](#basic-classification)
    - [FURTHER INFORMATION:](#further-information)
    - [CONTRIBUTORS](#contributors)
    - [CITATIONS](#citations)

<!-- /TOC -->

## NOTE FOR USERS<a id="usernote"></a>

The default branch for this repository is a development branch: `diagnostic_primers`. If you are looking for code to reproduce work from [Pritchard *et al.* (2012)](https://dx.doi.org/10.1371/journal.pone.0034498) or [Pritchard *et al.* (2013)](https://dx.doi.org/10.1111/j.1365-3059.2012.02678.x), please checkout the `master` branch, or download [release v0.1.3](https://github.com/widdowquinn/find_differential_primers/tree/v0.1.3).

- `diagnostic_primers`: 

[![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/diagnostic_primers/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)
[![Code Health](https://landscape.io/github/widdowquinn/find_differential_primers/diagnostic_primers/landscape.svg?style=flat)](https://landscape.io/github/widdowquinn/find_differential_primers/diagnostic_primers)
[![Build Status](https://travis-ci.org/widdowquinn/find_differential_primers.svg?branch=diagnostic_primers)](https://travis-ci.org/widdowquinn/find_differential_primers)

- `master`: 

[![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/master/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)
[![Code Health](https://landscape.io/github/widdowquinn/find_differential_primers/master/landscape.svg?style=flat)](https://landscape.io/github/widdowquinn/find_differential_primers/master)
[![Build Status](https://travis-ci.org/widdowquinn/find_differential_primers.svg?branch=master)](https://travis-ci.org/widdowquinn/find_differential_primers)

## NOTE FOR DEVELOPERS<a id="devnotenote"></a>

The default master branch for development is `diagnostic_primers`. We appreciate it when you [raise issues](https://github.com/widdowquinn/find_differential_primers/issues), or make contributions *via* [pull request](https://github.com/widdowquinn/find_differential_primers/pulls), especially if you follow the guidelines on the [wiki](https://github.com/widdowquinn/find_differential_primers/wiki).

- Current test coverage (`diagnostic_primers`): [![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/diagnostic_primers/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)

## Overview<a id="overview"></a>
This codebase performs automated finding of discriminatory (real-time) PCR or qPCR primers that distinguish between sets of provided genomes or other biological sequences of interest. Methods are provided specifically for identification of marker sequences useful for metabarcoding, that can discriminate well within a subset of bacterial genomes.

The `pdp` command-line program can be used to conduct this analysis in an interactive or scripted manner.

The `diagnostic_primers` Python model can be used within your own programs to automate primer and/or marker design.

### Third-party Packages

**This tool depends on some third-party software packages**

- [`Primer3`](https://sourceforge.net/projects/primer3/) **VERSION 1.1.4**: `Primer3` is the tool used to design primers. For compatibility with `EMBOSS`, version 1 of the software is essential.
- [`EMBOSS`](http://emboss.sourceforge.net/): This suite of tools is used to interact with `Primer3` and to perform *in silico* cross-hybridisation checks with `primersearch`.
- [`BLAST+`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download): This tool is used to screen primers against a database of off-target sequences with the `blastscreen` command.
- [`prodigal`](https://github.com/hyattpd/Prodigal): This program is used to identify candidate CDS features when using the `pdp filter` subcommand.
- [`MAFFT`](https://mafft.cbrc.jp/alignment/software/): This is required to align candidate metabarcoding marker sequences, when using the `pdp extract` subcommand

These packages can all be installed using the `conda` package manager, and the `bioconda` channel.

- [`bioconda`](https://bioconda.github.io/): a channel for the `conda` package manager, specialising in bioinformatics software.


### Recent changes

The new version of `diagnostic_primers` (formerly `find_differential_primers`) now uses a subcommand model, like the tools `git` and `subversion`. These execute the following subtasks, some or all of which may be required to perform a specific primer/marker design run.

- `pdp config`: Process/validate the configuration file and stitch input contig fragments/replace ambiguity symbols as necessary.
- `pdp filter`: Filter input genome sequences so that primers are designed only to targeted regions.
- `pdp eprimer3`/`e3`: Design amplifying primers on the input sequences
- `pdp blastscreen`/`bs`: Filter designed primers against a database of negative examples
- `pdp primersearch`/`ps`: Filter designed primers on their ability to amplify each input sequence
- `pdp classify`/`cl`: Classify designed primers by specificity for each class of input sequence
- `pdp extract`/`ex`: Extract amplicon sequences corresponding to diagnostic primer sets
- `pdp plot`/`pl`: Generate useful graphical output for interpretation of results

Each of these subcommands has specific help, accessible with `pdp <subcommand> -h` or `pdp <subcommand> --help`.

## Walkthrough <a id="walkthrough"></a>

This section describes a worked example analysis, proceeding from defining a config file to producing a diagnostic primer set result. All the files required for this analysis can be found in the subdirectory `tests/walkthrough`.

### 1. Producing and validating the config file

We begin with a small set of bacterial genomes: three *Pectobacterium* species. These are defined as `.fasta` sequences in the directory `tests/walkthrough/sequences`:

```bash
$ ls tests/walkthrough/sequences/
GCF_000011605.1.fasta	GCF_000291725.1.fasta	GCF_000749845.1.fasta
```

A basic config file defining the three genomes is provided as `tests/walkthrough/pectoconf.tab` in tab-separated tabular format (for your own data this can be generated in Excel, and saved as tab-separated text format). 

```text
# Pectobacterium genomes downloaded from GenBank/NCBI; genomovars inferred from ANIm
# Annotated Pba: genomovar 1
Pba_SCRI1043	Pectobacterium,atrosepticum_NCBI,gv1	tests/walkthrough/sequences/GCF_000011605.1.fasta	-
# Annotated Pwa: genomovars 2, 3
Pwa_CFBP_3304	Pectobacterium,wasabiae_NCBI,gv2	tests/walkthrough/sequences/GCF_000291725.1.fasta	-
# Annotated Pb	: genomovar 7
Pbe_NCPPB_2795	Pectobacterium,betavasculorum_NCBI,gv7	tests/walkthrough/sequences/GCF_000749845.1.fasta	-
```

Four columns are required: `name`; `classes` (a list of comma-separated labels); `FASTA location`; and `features`. At this point no features are defined (these would used to direct primer design to specified regions of the genome), so this column contains only the symbol `-` to mark it as being empty.

Comment lines can be included, starting with `#` as the first character. These are ignored in the analysis

To confirm that the config file is correctly-formatted, we use the `pdp config --validate` command:

```bash
$ pdp config --validate tests/walkthrough/pectoconf.tab
WARNING: Validation problems
    Pbe_NCPPB_2795 requires stitch (tests/walkthrough/sequences/GCF_000749845.1.fasta)
    Pwa_CFBP_3304 requires stitch (tests/walkthrough/sequences/GCF_000291725.1.fasta)
    Pwa_CFBP_3304 has non-N ambiguities (tests/walkthrough/sequences/GCF_000291725.1.fasta)
```

This tells us that the first two genome files are in multiple parts so must be concatenated for this analysis, and that the second file also has ambiguity base symbols that are not `N`. These must be replaced accordingly for the analysis to proceed using the third-party tools.

### 2. Fix sequences for analysis

The `pdp config` command can fix input genomes for analysis, with the `--fix_sequences` argument. This writes new sequences having the necessary changes (stitching, replacing ambiguity symbols) and generates a new config file that refers to the fixed sequences. We specify the path to this new file with the argument `--fix_sequences <NEWCONFIG>.json`: 

```bash
$ pdp config --fix_sequences tests/walkthrough/fixed.json \
                tests/walkthrough/pectoconf.tab
```

This writes corrected sequences to the `tests/walkthrough/sequences` subdirectory, and a new config file to `tests/walkthrough/fixed.json` (in `JSON` format, which is how `pdp` prefers to receive configuration files):

```bash
$ tree tests/walkthrough/
tests/walkthrough/
├── fixed.json
├── pectoconf.tab
└── sequences
    ├── GCF_000011605.1.fasta
    ├── GCF_000291725.1.fasta
    ├── GCF_000291725.1_concat.fas
    ├── GCF_000291725.1_concat_noambig.fas
    ├── GCF_000749845.1.fasta
    └── GCF_000749845.1_concat.fas
```

### 3. Defining CDS features on each genome (optional)

For prokaryotic genomes, we can use a genecaller to predict gene features with the `pdp filter --prodigal` command. This generates ba GFF file describing CDS feature locations on each genome, and FASTA files of the CDS sequences. In the primer design stage, we can take CDS locations into account to retain only primers that amplify within CDS regions (assuming greater sequence conservation in CDS regions).

To use the genecaller, we must provide an appropriate config file (the `fixed.json` config file), and the path to a new config file that will contain information about the predicted features (we'll call this `fixed_with_features.json`). We will tell `prodigal` to place the predicted gene locations in the subdirectory `tests/walkthrough/prodigal`:

```bash
pdp filter --prodigal --outdir tests/walkthrough/prodigal \
                tests/walkthrough/fixed.json \
                tests/walkthrough/fixed_with_features.json
```

The new directory containing genecaller output is created for us, as is the new config file:

```bash
$ tree tests/walkthrough/
tests/walkthrough/
├── fixed.json
├── fixed_with_features.json
├── pectoconf.tab
├── prodigal
│   ├── GCF_000011605.1.features
│   ├── GCF_000011605.1.gff
│   ├── GCF_000291725.1_concat_noambig.features
│   ├── GCF_000291725.1_concat_noambig.gff
│   ├── GCF_000749845.1_concat.features
│   └── GCF_000749845.1_concat.gff
└── sequences
[…]
```

### 4. Design primers to each genome in bulk

Using the most informative config file (`fixed_with_features.json`) we can design primers to each input genome with the EMBOSS `ePrimer3` package. At a minimum, we need to give the `pdp eprimer3` command the input config file, and the path to an output config file that will contain information about the primer description files for each genome. 

To use only the CDS regions we designed in the previous step, we need to give the `--filter` argument.

We will also use the `--outdir` argument to tell `pdp` where to put the `ePrimer3` output files:

```bash
$ pdp eprimer3 --filter --outdir tests/walkthrough/eprimer3 \
                  tests/walkthrough/fixed_with_features.json \
                  tests/walkthrough/with_primers.json
```

This places the output of `ePrimer3` into its own directory, and generates `JSON` files that describe the primers for each of the genomes.

```bash
$ tree tests/walkthrough/
tests/walkthrough/
├── eprimer3
│   ├── GCF_000011605.1.eprimer3
│   ├── GCF_000011605.1_named.eprimer3
│   ├── GCF_000011605.1_named.json
│   ├── GCF_000291725.1_concat_noambig.eprimer3
│   ├── GCF_000291725.1_concat_noambig_named.eprimer3
│   ├── GCF_000291725.1_concat_noambig_named.json
│   ├── GCF_000749845.1_concat.eprimer3
│   ├── GCF_000749845.1_concat_named.eprimer3
│   └── GCF_000749845.1_concat_named.json
├── fixed.json
├── fixed_with_features.json
├── pectoconf.tab
├── prodigal
[…]
├── sequences
[…]
└── with_primers.json
```


### 5. Deduplicate primer sets (optional)

When designing thermodynamically plausible primers to closely-related genomes, it is very likely that identical primer sets will be created on distinct genomes due to the sequence similarity between genomes. Carrying these duplicate primer sets through to later analysis stages can be a significant unnecessary computational load. To remove identical primer sets, we can use use the `pdp dedupe` command:

```bash
$ pdp dedupe --dedupedir tests/walkthrough/deduped \
                    tests/walkthrough/with_primers.json \
                    tests/walkthrough/deduped_primers.json
```

This places new `JSON` files with deduplicated primers in the directory `tests/walkthrough/deduped`, and creates a new config file in `deduped_primers.json` that points to these files.


### 6. Screen primers against `BLASTN` database (optional)

Now that primers have been designed and deduplicated, they can be screened against a `BLASTN` database to identify and filter out any primers that have potential cross-amplification. In general, we advise that this step is used not to demonstrate potential for cross-hybridisation/amplification, but to exclude primers that have *any* theoretical potential for off-target binding. That is, we recommend this process only to aid a negative screen to remove non-specific primer sets.

The `pdp blastscreen` subcommand requires the path to a suitable pre-compiled `BLASTN` database to be provided with the argument `--db` (there is a `BLASTN` *E. coli* genome database in the subdirectory `tests/walkthrough/blastdb/`). By specifiying `--outdir`, we place `BLAST` output in the `tests/walkthrough/blastn` subdirectory. The input config file we use is the `with_primers.json` file produced in the last step, having locations of the unscreened primer files, and we specify that the command writes a new config file called `screened.json` that points to a reduced set of primers, with the potentially off-target/non-specific sets excluded.

**No sequences are deleted as a result of this action.**

```bash
$ pdp blastscreen --db tests/walkthrough/blastdb/e_coli_screen.fna \
                     --outdir tests/walkthrough/blastn \
                     tests/walkthrough/deduped_primers.json \
                     tests/walkthrough/screened.json
```

This screen produces the new subdirectory `tests/walkthrough/blastn` that contains all the primer sequences in FASTA format, and the tabular output of the BLAST search. New files are added to the `eprimer3` subdirectory (with suffix `_screened`), describing the reduced sets of primers, post-screening.

```bash
$ tree tests/walkthrough/
tests/walkthrough/
├── blastdb
│   ├── e_coli_screen.fna.nhr
│   ├── e_coli_screen.fna.nin
│   └── e_coli_screen.fna.nsq
├── blastn
│   ├── GCF_000011605.1_primers.fasta
│   ├── GCF_000011605.1_primers.tab
│   ├── GCF_000291725.1_concat_noambig_primers.fasta
│   ├── GCF_000291725.1_concat_noambig_primers.tab
│   ├── GCF_000749845.1_concat_primers.fasta
│   └── GCF_000749845.1_concat_primers.tab
├── eprimer3
│   ├── GCF_000011605.1.eprimer3
│   ├── GCF_000011605.1_named.eprimer3
│   ├── GCF_000011605.1_named.json
│   ├── GCF_000011605.1_named_screened.fasta
│   ├── GCF_000011605.1_named_screened.json
│   ├── GCF_000291725.1_concat_noambig.eprimer3
│   ├── GCF_000291725.1_concat_noambig_named.eprimer3
│   ├── GCF_000291725.1_concat_noambig_named.json
│   ├── GCF_000291725.1_concat_noambig_named_screened.fasta
│   ├── GCF_000291725.1_concat_noambig_named_screened.json
│   ├── GCF_000749845.1_concat.eprimer3
│   ├── GCF_000749845.1_concat_named.eprimer3
│   ├── GCF_000749845.1_concat_named.json
│   ├── GCF_000749845.1_concat_named_screened.fasta
│   └── GCF_000749845.1_concat_named_screened.json
[…]
├── fixed.json
├── fixed_with_features.json
├── pectoconf.tab
├── prodigal
[…]
├── sequences
[…]
└── with_primers.json
```

The new configuration file can be used in the `primersearch` cross-hybridisation detection stage.

### 7. Test primers against input sequences for crosshybridisation with `primersearch`

To identify which primers might be diagnostically useful for any of the classes/labels defined in the config file, we must test whether they potentially amplify the other genomes from the input set. We use the `EMBOSS` tool `primersearch` to check whether any of the designed primers (after screening against the `BLAST` database) also have potential to amplify any of the other input sequences.

To do this we pass the appropriate `.json` config file, and specify a directory to hold `primersearch` output, as well provide the location of a new config file that will hold data about our crosshybridisation screen.

```bash
pdp primersearch \
       --outdir tests/walkthrough/primersearch \
       tests/walkthrough/screened.json \
       tests/walkthrough/primersearch.json
```

After running this command, the new directory `tests/walkthrough/primersearch` is produced, containing a new primer file (with extension `.tab` - this was needed for `primersearch`) and a `.json` file for each input sequence. In addition, there is a new `.primersearch` file for each comparison of primer sequence set against one of the input genome sequences.

```bash
$ tree tests/walkthrough
tests/walkthrough
[...]
├── primersearch
│   ├── Pba_SCRI1043_primers.tab
│   ├── Pba_SCRI1043_primersearch.json
│   ├── Pba_SCRI1043_ps_Pbe_NCPPB_2795.primersearch
│   ├── Pba_SCRI1043_ps_Pwa_CFBP_3304.primersearch
│   ├── Pbe_NCPPB_2795_primers.tab
│   ├── Pbe_NCPPB_2795_primersearch.json
│   ├── Pbe_NCPPB_2795_ps_Pba_SCRI1043.primersearch
│   ├── Pbe_NCPPB_2795_ps_Pwa_CFBP_3304.primersearch
│   ├── Pwa_CFBP_3304_primers.tab
│   ├── Pwa_CFBP_3304_primersearch.json
│   ├── Pwa_CFBP_3304_ps_Pba_SCRI1043.primersearch
│   └── Pwa_CFBP_3304_ps_Pbe_NCPPB_2795.primersearch
├── primersearch.json
[...]
```

The new `primersearch.json` config file contains information about this crosshybridisation screen, and can be used for identification and extraction of diagnostic primer sequence sets.

### 8. Classify the primers by diagnostic capability with `classify`

To extract useful information from `primersearch` output, and classify the primer sets by their ability to amplify only genomes belonging to a specific named group in the configuration file, we use the `pdp classify` subcommand. This examines the `primersearch` output and reports back diagnostic primer sets.

We pass the `.json` file produced by the `primersearch` run, and the path to a directory for the output of the `pdp classify` subcommand:

```bash
pdp classify \
       tests/walkthrough/primersearch.json \
       tests/walkthrough/classify
```

The new directory contains `.json` and `.ePrimer3` format files for each set of primers diagnostic to a given class, and summary information in `summary.tab` and `results.json` files.

```bash
$ tree tests/walkthrough/classify
tests/walkthrough/classify
├── atrosepticum_NCBI_primers.ePrimer3
├── atrosepticum_NCBI_primers.json
├── betavasculorum_NCBI_primers.ePrimer3
├── betavasculorum_NCBI_primers.json
├── gv1_primers.ePrimer3
├── gv1_primers.json
├── gv2_primers.ePrimer3
├── gv2_primers.json
├── gv7_primers.ePrimer3
├── gv7_primers.json
├── Pectobacterium_primers.ePrimer3
├── Pectobacterium_primers.json
├── results.json
├── summary.tab
├── wasabiae_NCBI_primers.ePrimer3
└── wasabiae_NCBI_primers.json
```

The `summary.tab` file contains a table showing how many primer sets were designed that are potentially specific to each of the input genomes and classes/labels:

```bash
$ cat tests/walkthrough/classify/summary.tab 
Group	NumPrimers	Primers
Pectobacterium	4	tests/walkthrough/classify/Pectobacterium_primers.json
atrosepticum_NCBI	2	tests/walkthrough/classify/atrosepticum_NCBI_primers.json
betavasculorum_NCBI	1	tests/walkthrough/classify/betavasculorum_NCBI_primers.json
gv1	2	tests/walkthrough/classify/gv1_primers.json
gv2	2	tests/walkthrough/classify/gv2_primers.json
gv7	1	tests/walkthrough/classify/gv7_primers.json
```


## Usage<a id="usage"></a>


### `pdp config`<a id="config"></a>

The `config` subcommand handles interactions with the configuration file for a primer design run. Configuration files can be provided in one of two formats:

1. `<config>.tab`: a plain text, tab-separated file descrbing the input data in a table of multiple columns [as defined below](#config_tab_format). This is intended to be an easily human-readable file, that can be prepared and edited in a spreadsheet application such as Google Sheets, or Microsoft Excel. `pdp config` will recognise `.tab` or `.conf` as a file extension.
2. `<config>.json`: a `JSON` format file describing the input data. This is not intended to be human-readable, but can be converted to and from the `.tab` format using `pdp` config. Further steps in the primer design process require that the configuration is provided in `JSON` format.

#### Converting between `.tab` and `JSON` format config files

-*1. `.tab` to JSON**

Provide the path to the output `JSON` file as an argument to `--to_json`, and the path to the `.tab` config file as input:

```bash
pdp config --to_json <OUTPUT>.json <INPUT>.tab
```

-*2. `JSON` to `tab`**

Provide the path to the output `.tab` file as an argument to `--to_tab`, and the path to the `JSON` config file as input:

```bash
pdp config --to_tab <OUTPUT>.tab <INPUT>.json
```

#### Validate a config file

`pdp` can examine the contents of a config file and determine whether it conforms to the required specification, and whether the sequences used for input require stitching, or replacement of ambiguity codons. To validate a config file, use the `--validate` flag:

```bash
$ pdp config --validate <INFILE>.tab
$ pdp config --validate <INFILE>.json
```

#### Repair input sequences

For use with this primer design tool, the input sequences must be concatenated, and cannot contain non-`N` ambiguity base symbols. `pdp` can nondestructively repair input sequences by stitching sequence fragments/contigs together, and replacing all ambiguity symbols with `N`.

```bash
pdp config --fix_sequences <REPAIRED>.json <INPUT>.[tab|json]
```

The repaired sequences are written to new files in the same directory as the input file, with one of the following suffixes:

- `_concat`: the sequence was concatenated
- `_noambig`: the sequence had ambiguity symbols replaced
- `_concat_noambig`: the sequence was concatenated, and ambiguity symbols were replaced

such that an input file `<SEQUENCE>.fas` may be repaired to generate the file `<SEQUENCE>_concat_noambig.fas` in the same directory as the original file, and a new config file pointing to the modified sequences is written to `<REPAIRED>.json`.

### `pdp prodigal`<a id="prodigal"></a>

The `prodigal` (or `prod`) subcommand runs the [`prodigal`](https://github.com/hyattpd/Prodigal) prokaryotic gene feature-calling package on the sequences listed in the passed configuration file. A new configuration file, specifying the location of the feature file for each input sequence, is written to the specified output file location.

#### Default feature prediction

`prodigal` feature prediction is run on the sequences listed in `<INPUT>.json`, and a new config file written to `<OUTPUT>.json` with the locations of the feature predictions indicated.

```bash
pdp prodigal <INPUT>.json <OUTPUT>.json
```

To overwrite existing output, pass the `-f` or `--force` argument:

```bash
pdp prodigal --force <INPUT>.json <OUTPUT>.json
```

#### Specify location to write `prodigal` predictions

By default, `pdp` writes output to the subdirectory `prodigal`. To put the feature predictions in another location, pass the directory you want to place the `prodigal` output (here, `<OUTDIR>`) as the `--outdir` argument (and use the `-f`/`--force` argument to overwrite existing output):

```bash
pdp prodigal --outdir <OUTDIR> <INPUT>.json <OUTPUT>.json
```

#### Specify the location of the `prodigal` executable

By default `pdp` will look for `prodigal` in your `$PATH`. A different executable can be specified with the `--prodigal` argument:

```bash
pdp prodigal --prodigal <PATH_TO_PRODIGAL> <INPUT>.json <OUTPUT>.json
```

### `pdp eprimer3`<a id="eprimer3"></a>

The `eprimer3` command runs primer prediction on each of the input sequences listed in the passed input configuration file. The tool used by `pdp` is the [EMBOSS `ePrimer3` package](http://bioinf.ibun.unal.edu.co/cgi-bin/emboss/help/eprimer3). A new configuration file is written describing the locations of the predicted primers.


#### Default primer prediction

Primer prediction is run on the sequences listed in `<INPUT>.json`, and the new config file written to `<OUTPUT>.json`.

```bash
pdp eprimer3 <INPUT>.json <OUTPUT>.json
```

#### Change the number of predicted primers per input sequence

By default only 10 primers are predicted per sequence. This is a choice made for speed of testing, and is unlikely to be enough to useful for designing diagnostic primers for a prokaryotic genome. Overall runtime increases exponentially with the number of primers that need to be tested for cross-hybridisation, and a suitable choice of value will depend strongly on the dataset being used. To specify the number of primers to be designed for each input sequence, use the `--numreturn` argument. For example, to design 2000 primers per input sequence, use:

```bash
pdp eprimer3 --numreturn 2000 <INPUT>.json <OUTPUT>.json
```

#### Change primer design parameters

All parameters for `eprimer3` are available to be changed at the command line. There are a large number of these arguments, and they are all described in the help text (use: `pdp eprimer3 -h`), but some useful examples are listed below:

-*Specify primer lengths**

To specify an optimal primer oligo size, and an acceptable (minimum/maximum) range of sizes, use the `--osize`, `--minsize`, `--maxsize` arguments, e.g.:

```bash
pdp eprimer3 --osize 25 --minsize 20 --maxsize 30 <INPUT>.json <OUTPUT>.json
```

-*Specify primer thermodynamics**

To specify optimal, minimum and maximum melting temperatures (Tm) for the predicted primers, use the `--opttm`, `--mintm`, and `--maxtm` arguments, e.g.:

```bash
pdp eprimer3 --opttm 65 --mintm 62 --maxtm 68 <INPUT>.json <OUTPUT>.json
```

-*Specify amplicon lengths**

To specify an optimal amplicon size, and an acceptable (minimum/maximum) range of sizes, use the `--psizeopt`, `--psizemin`, `--psizemax` arguments, e.g.:

```bash
pdp eprimer3 --psizeopt 200 --psizemin 190 --psizemax 210 <INPUT>.json <OUTPUT>.json
```

#### Specify location to write primer prediction output

By default, `pdp` writes output to the subdirectory `eprimer3`. To put the primer predictions in another location, pass the directory you want to place the output (here, `<OUTDIR>`) as the `--outdir` argument (and use the `-f`/`--force` argument to overwrite existing output):

```bash
pdp eprimer3 --outdir <OUTDIR> <INPUT>.json <OUTPUT>.json
```

#### Specify the location of the `eprimer3` executable

By default `pdp` looks for the EMBOSS `eprimer3` executable in your `$PATH`, but its location can be specified with the `--eprimer3` argument:

```
pdp eprimer3 --eprimer3 <PATH_TO_EPRIMER3> <INPUT>.json <OUTPUT>.json
```

### `pdp blastscreen`<a id="blastscreen"></a>

The `blastscreen` command screens predicted primers against a local `BLASTN` nucleotide database. Primer pairs for which at least one member produces a match in the `BLAST` database are excluded. The tool used by `pdp` is a [local `BLAST+` installation](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). `BLAST` output is written to a new directory, and a new configuration file is written describing the primer sets that pass the screen (i.e. have no matches in the database).

#### Basic screening

The `BLAST` database to screen the primers described in `<INPUT>.json` against is passed as `<BLASTDB>`. The `BLAST` results are written to the directory `<BLASTOUT>`, and the new config file written to `<OUTPUT>.json`.

```bash
pdp blastscreen --db <BLASTDB> --outdir <BLASTOUT> <INPUT>.json <OUTPUT>.json
```

The default database is `nr` and the default output directory is `blastn`.

#### Controlling sensitivity

A single argument `--maxaln` is used to control sensitivity of primer matching. This describes the maximum number of identities allowed in the `BLAST` match before the primer pair is excluded (default=`15`).

```bash
pdp blastscreen --db <BLASTDB> --outdir <BLASTOUT> --maxaln 25 <INPUT>.json <OUTPUT>.json
```

#### Specify the location of the `BLAST+` executable

The location of the `BLAST+` executable can be provided with the `--blastn` argument.

```bash
pdp blastscreen --db <BLASTDB> --outdir <BLASTOUT> --blastn <BLASTNPATH> <INPUT>.json <OUTPUT>.json
```

#### Control the number of threads used

The `BLAST` screen is parallelised on as many threads as are available, by default. The number of worker threads can be controlled with the `-w` argument.

```bash
pdp blastscreen --db <BLASTDB> --outdir <BLASTOUT> -w 4 <INPUT>.json <OUTPUT>.json
```

#### Use the SGE/OGE scheduler

The `BLAST` screen can be parallelised using a Sun Grid Engine variant, such as Son of Grid Engine, Open Grid Engine, or Univa Grid Engine. To specify this scheduler, use the `-s` argument with the value `SGE`.

```bash
pdp blastscreen --db <BLASTDB> --outdir <BLASTOUT> -s SGE <INPUT>.json <OUTPUT>.json
```

### `pdp primersearch`<a id="primersearch"></a>

The `primersearch` command performs *in silico* hybridisation of predicted primers against each of the input genomes, so that cross-hybridising primers can be identified. The tool used by `pdp` is the [EMBOSS `primersearch` tool](http://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html). `primersearch` output is written to a new directory, and a new configuration file is written describing the cross-hybridisation results.

#### Basic cross-hybridisation

The configuration file describing primers to use is passed as `<INPUT>.json`, and the path to the new output directory containing `primersearch` output is passed as `<OUTDIR>`. The path to write the new configuration file describing cross-hybridisation information is provided as `<OUTPUT>.json`.

```bash
pdp primersearch --outdir <OUTDIR> <INPUT>.json <OUTPUT>.json
```

#### Controlling sensitivity

A single argument `--mismatchpercent` is used to control the sensitivity with which `primersearch` thinks primers cross-hybridise. This describes the maximum percentage of mismatches (in the range `[0, 1]`) allowed in the primer match before `primersearch` considers that hybridisation is not possible. The default value is `0.1`.

```bash
pdp primersearch --outdir <OUTDIR> --mismatchpercent 0.25 <INPUT>.json <OUTPUT>.json
```

#### Specify the location of the `primersearch` executable

The location of the `primersearch` executable can be provided with the `--primersearch` argument.

```bash
pdp primersearch --outdir <OUTDIR> --primersearch <PSPATH> <INPUT>.json <OUTPUT>.json
```

#### Control the number of threads used

By default, cross-hybridisation is parallelised on as many threads as are available. The number of worker threads can be controlled with the `-w` argument.

```bash
pdp primersearch --outdir <OUTDIR> -w 4 <INPUT>.json <OUTPUT>.json
```

#### Use the SGE/OGE scheduler

The cross-hybridisation screen can also be parallelised using a Sun Grid Engine variant, such as Son of Grid Engine, Open Grid Engine, or Univa Grid Engine. To specify this scheduler, use the `-s` argument with the value `SGE`.

```bash
pdp primersearch --outdir <OUTDIR> -s SGE <INPUT>.json <OUTPUT>.json
```


### `pdp classify`<a id="classify"></a>

The `classify` command takes the output from the `primersearch` step, and identifies primer sets that uniquely amplify each of the target groups defined in the corresponding `.json` configuration file.

#### Basic classification

The configuration file describing primers and their PrimerSearch results is passed as `<INPUT>.json`, and the path to the new output directory that will contain the sets of primers predicted to be specific to each group, and summary information as `<OUTDIR>`.

```bash
pdp classify <INPUT>.json <OUTDIR>
```

This will produce a summary tab-separated plain text table (`summary.tab`), a `JSON` format file describing the complete set of results (`results.json`), and then a pair of `.json` and `.ePrimer3` format files for each defined group for which predicted diagnostic primers could be derived.




## FURTHER INFORMATION:
For further technical information, please read the comments contained within the top of each '*.py' file as well as the Supporting Information (['Methods S1' document](doi:10.1371/journal.pone.0034498.s006)) of [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498).

## CONTRIBUTORS
- [Leighton Pritchard](https://github.com/widdowquinn)
- [Benjamin Leopold](https://github.com/cometsong)
- [Michael Robeson](https://github.com/mikerobeson)
- [Rory McLeod](https://github.com/rory-mcleod)

## CITATIONS
Please refer to the following for methodological details:

- Pritchard L _et al._ (2012) "Alignment-Free 
Design of Highly Discriminatory Diagnostic Primer Sets for _Escherichia coli_ O104:H4 Outbreak Strains." _PLoS ONE_ **7**(4): e34498. [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498) - _Method description and application to human bacterial pathogens, sub-serotype resolution_
- Pritchard L _et al._ (2013) "Detection of phytopathogens of the genus _Dickeya_ using a PCR primer 
prediction pipeline for draft bacterial genome sequences." _Plant Pathology_, **62**, 587-596
[doi:10.1111/j.1365-3059.2012.02678.x](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-3059.2012.02678.x/full) - _Application to plant pathogens, species-level resolution_
