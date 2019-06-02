# README.md (`pdp`/`diagnostic_primers`)

<!-- TOC -->

- [Documentation](#documentation)
- [CITATIONS](#citations)
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
    - [8. Classify primers by predicted diagnostic capability with `classify`](#8-classify-primers-by-predicted-diagnostic-capability-with-classify)
    - [9. Extract and assess candidate metabarcoding markers](#9-extract-and-assess-candidate-metabarcoding-markers)
- [FURTHER INFORMATION](#further-information)
- [CONTRIBUTORS](#contributors)

<!-- /TOC -->

## Documentation

Full documentation about usage of `pdp`/`diagnostic_primers` can be found on ReadTheDocs:

- [`pdp` documentation](https://pdp.readthedocs.io/en/diagnostic_primers/index.html)

## CITATIONS

If you use `pdp`/`diagnostic_primers`, please cite one or both of the following papers:

- Pritchard L _et al._ (2012) "Alignment-Free 
Design of Highly Discriminatory Diagnostic Primer Sets for _Escherichia coli_ O104:H4 Outbreak Strains." _PLoS ONE_ **7**(4): e34498. [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498) - _Method description and application to human bacterial pathogens, sub-serotype resolution_
- Pritchard L _et al._ (2013) "Detection of phytopathogens of the genus _Dickeya_ using a PCR primer 
prediction pipeline for draft bacterial genome sequences." _Plant Pathology_, **62**, 587-596
[doi:10.1111/j.1365-3059.2012.02678.x](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-3059.2012.02678.x/full) - _Application to plant pathogens, species-level resolution_


## NOTE FOR USERS<a id="usernote"></a>

The default branch for this repository is a development branch: `diagnostic_primers`. If you are looking for code to reproduce work from [Pritchard *et al.* (2012)](https://dx.doi.org/10.1371/journal.pone.0034498) or [Pritchard *et al.* (2013)](https://dx.doi.org/10.1111/j.1365-3059.2012.02678.x), please checkout the `master` branch, or download [release v0.1.3](https://github.com/widdowquinn/find_differential_primers/tree/v0.1.3).

- `diagnostic_primers`: 

[![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/diagnostic_primers/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)
[![Build Status](https://travis-ci.org/widdowquinn/find_differential_primers.svg?branch=diagnostic_primers)](https://travis-ci.org/widdowquinn/find_differential_primers)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/bf07262b60a74d89839f241d3c61de10)](https://www.codacy.com/app/widdowquinn/find_differential_primers?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=widdowquinn/find_differential_primers&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/widdowquinn/find_differential_primers/badge)](https://www.codefactor.io/repository/github/widdowquinn/find_differential_primers)

- `master`: 

[![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/master/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)
[![Build Status](https://travis-ci.org/widdowquinn/find_differential_primers.svg?branch=master)](https://travis-ci.org/widdowquinn/find_differential_primers)

## NOTE FOR DEVELOPERS<a id="devnotenote"></a>

The default master branch for development is `diagnostic_primers`. We appreciate it when you [raise issues](https://github.com/widdowquinn/find_differential_primers/issues), or make contributions *via* [pull request](https://github.com/widdowquinn/find_differential_primers/pulls), especially if you follow the guidelines on the [wiki](https://github.com/widdowquinn/find_differential_primers/wiki). Writing tests to demonstrate/catch bugs, or to improve coverage, is especially appreciated!

- Current test coverage (`diagnostic_primers`): [![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/diagnostic_primers/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)

## Overview<a id="overview"></a>
This package automates discovery of discriminatory PCR or qPCR primers that can distinguish between arbitrary subgroups of input genomes or other biological sequences of interest. The package also aids identification of novel metabarcoding marker sequences. The package can be used in a number of ways:

- The `pdp` command-line program can be used to conduct this analysis in an interactive or scripted manner.
- The `diagnostic_primers` Python model can be used within your own programs to automate primer and/or marker design.

## Third-party Packages

**This tool depends on some third-party software packages**

- `Primer3` **VERSION 1.1.4**: `Primer3` is the tool used to design primers. For compatibility with `EMBOSS`, version 1 of the software is essential. [[webpage](https://sourceforge.net/projects/primer3/)]
- `EMBOSS`: This suite of tools is used to interact with `Primer3` and to perform *in silico* cross-hybridisation checks with `primersearch`. [[webpage](http://emboss.sourceforge.net/)]
- `BLAST+`: This tool is used to screen primers against a database of off-target sequences with the `blastscreen` command. [[webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)]
- `prodigal`: This program is used to identify candidate prokaryotic CDS features when using the `pdp filter` subcommand. [[webpage](https://github.com/hyattpd/Prodigal)]
- `MAFFT`: This is required to align candidate metabarcoding marker sequences, when using the `pdp extract` subcommand. [[webpage](https://mafft.cbrc.jp/alignment/software/)]

These packages can all be installed using the `conda` package manager, and the `bioconda` channel.

- `bioconda`: a channel for the `conda` package manager, specialising in bioinformatics software. [[webpage](https://bioconda.github.io/)]

```bash
conda install primer3=1.1.4 emboss blast prodigal mafft
```


## Recent changes

The new version of `diagnostic_primers` (formerly `find_differential_primers`) now uses a subcommand model (like the tools `git` and `subversion`). These execute the following subtasks, some or all of which may be required to perform a specific primer/marker design run. This change has been made for flexibility in composing workflows.

- `pdp config`: Process/validate the configuration file and stitch input contig fragments/replace ambiguity symbols as necessary.
- `pdp filter`: Filter input genome sequences so that primers are designed only to targeted regions.
- `pdp eprimer3`/`e3`: Design amplifying primers on the input sequences
- `pdp blastscreen`/`bs`: Filter designed primers against a database of negative examples
- `pdp primersearch`/`ps`: Filter designed primers on their ability to amplify each input sequence
- `pdp classify`/`cl`: Classify designed primers by specificity for each class of input sequence
- `pdp extract`/`ex`: Extract amplicon sequences corresponding to diagnostic primer sets

Each of these subcommands has specific help, accessible with `pdp <subcommand> -h` or `pdp <subcommand> --help`.

## Walkthrough <a id="walkthrough"></a>

This section describes a worked example analysis, proceeding from defining a config file to producing a diagnostic primer set result. All the files required for this analysis can be found in the subdirectory `tests/walkthrough`.

### 1. Producing and validating the config file

We begin with a small set of bacterial genomes: three *Pectobacterium* species. These are defined as `.fasta` sequences in the directory `tests/walkthrough/sequences`:

```bash
$ ls tests/walkthrough/sequences/
GCF_000011605.1.fasta	GCF_000291725.1.fasta	GCF_000749845.1.fasta
```

A basic config file defining the three genomes is also provided as `tests/walkthrough/pectoconf.tab` in tab-separated tabular format (you can generate equivalent files for your own data using Excel, and saving the sheet as tab-separated text format). 

```text
# Pectobacterium genomes downloaded from GenBank/NCBI; genomovars inferred from ANIm
# Annotated Pba: genomovar 1
Pba_SCRI1043	Pectobacterium,atrosepticum_NCBI,gv1	tests/walkthrough/sequences/GCF_000011605.1.fasta	-
# Annotated Pwa: genomovars 2, 3
Pwa_CFBP_3304	Pectobacterium,wasabiae_NCBI,gv2	tests/walkthrough/sequences/GCF_000291725.1.fasta	-
# Annotated Pb	: genomovar 7
Pbe_NCPPB_2795	Pectobacterium,betavasculorum_NCBI,gv7	tests/walkthrough/sequences/GCF_000749845.1.fasta	-
```

Four columns are required (no headers are read), one row per input genome:

- a unique name for the genome
- a list of comma-separated labels for this genome - this will determine specificity of predicted primer sets
- path to a `FASTA` file of the input genome sequence - this can be a draft genome in multiple fragments
- (optional) path to a `.gff` file descrigbing genome features - if no features are provided, this should be a hyphen/dash (`-`).

In the walkthrough config file no features are defined, so this column contains only the symbol `-` to mark it as being empty.

Comment lines can be included, starting with `#` as the first character. These are ignored in the analysis

To confirm that the config file is correctly-formatted, we use the `pdp config --validate` command:

```bash
$ pdp config --validate tests/walkthrough/pectoconf.tab
WARNING: Validation problems
    Pbe_NCPPB_2795 requires stitch (tests/walkthrough/sequences/GCF_000749845.1.fasta)
    Pwa_CFBP_3304 requires stitch (tests/walkthrough/sequences/GCF_000291725.1.fasta)
    Pwa_CFBP_3304 has non-N ambiguities (tests/walkthrough/sequences/GCF_000291725.1.fasta)
```

This tells us that the first two genome files are in multiple parts so must be concatenated for this analysis, and that the second file also has ambiguity base symbols that are not `N`. These must be replaced by `N` for the analysis to proceed using the third-party tools.

### 2. Fix sequences for analysis

The `pdp config` command will fix input genomes for analysis, when the `--fix_sequences` argument is supplied. In this case, `pdp` writes new sequences having the necessary changes (stitching, replacing ambiguity symbols) and generates a new config file that refers to the fixed sequences. The path to the new config file is given with the argument `--fix_sequences <NEWCONFIG>.json`: 

```bash
$ pdp config --fix_sequences tests/walkthrough/fixed.json \
                tests/walkthrough/pectoconf.tab
```

This writes corrected sequences to the `tests/walkthrough/sequences` subdirectory, and a new config file to `tests/walkthrough/fixed.json` (in `JSON`, rather than tabular format):

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

For prokaryotic genomes, we can use a genecaller to predict gene features with the `pdp filter --prodigal` command. This generates a `GFF` file describing CDS feature locations on each genome, and `FASTA` files of the CDS sequences. This can be useful because, in the primer design stage, we can take CDS locations into account to retain only primers that amplify within CDS regions.

To use the genecaller, we must provide an appropriate config file (the `fixed.json` config file), and the path to a new config file that will contain information about the predicted features (we'll call this `fixed_with_features.json`). We will tell `prodigal` to place the predicted gene locations in the subdirectory `tests/walkthrough/prodigal`:

```bash
pdp filter --prodigal --outdir tests/walkthrough/prodigal \
                tests/walkthrough/fixed.json \
                tests/walkthrough/fixed_with_features.json
```

The new directory containing genecaller output is created, as is the new config file:

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

Using the config file `fixed_with_features.json` we can design primers to each input genome with the EMBOSS `ePrimer3` package. At a minimum, we need to give the `pdp eprimer3` command the input config file, and the path to an output config file that will contain additional information about the primers designed to each genome. We can also use the `--outdir` argument to provide a path to a directory where `pdp` will put the `ePrimer3` output files.

If we wish to use only primers that amplify within the CDS regions we designed in the previous step, we need to give the `--filter` argument. Otherwise `pdp eprimer3` will ignore this information.

```bash
$ pdp eprimer3 --filter --outdir tests/walkthrough/eprimer3 \
                  tests/walkthrough/fixed_with_features.json \
                  tests/walkthrough/with_primers.json
```

This places the output of `ePrimer3` into its own directory (`.eprimer3` files), and generates `JSON` files that describe the primers for each of the genomes. Additionally, a `.bed` file is produced describing the location of each primer set on the source genome, which can be visualised using a tool such as IGV.

```bash
$ tree tests/walkthrough/
tests/walkthrough/
├── eprimer3
│   ├── GCF_000011605.1.eprimer3
│   ├── GCF_000011605.1_named.bed
│   ├── GCF_000011605.1_named.eprimer3
│   ├── GCF_000011605.1_named.json
│   ├── GCF_000291725.1_concat_noambig.eprimer3
│   ├── GCF_000291725.1_concat_noambig_named.bed
│   ├── GCF_000291725.1_concat_noambig_named.eprimer3
│   ├── GCF_000291725.1_concat_noambig_named.json
│   ├── GCF_000749845.1_concat.eprimer3
│   ├── GCF_000749845.1_concat_named.bed
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

When designing thermodynamically plausible primers to closely-related genomes, it is very likely that identical primer sets will be created on distinct genomes due to their sequence conservation. Carrying these duplicate primer sets through to later analysis stages can result in significant unnecessary computational load. To remove identical primer sets, we can use the `pdp dedupe` command:

```bash
$ pdp dedupe --dedupedir tests/walkthrough/deduped \
                    tests/walkthrough/with_primers.json \
                    tests/walkthrough/deduped_primers.json
```

This examines the primers described in the input config file (`with_primers.json`) and places equivalent sets of new `JSON` and `.bed` files, with all redundant primers removed, in the specified directory (`tests/walkthrough/deduped`), and creates a new config file (`deduped_primers.json`) that points to these files.


### 6. Screen primers against `BLASTN` database (optional)

Now that primers have been designed and deduplicated, it can be useful to screen them against a `BLASTN` database to identify and filter out any primers that have potential for off-target amplification. In general, we advise that this step is used not to demonstrate potential for target amplification, but only to exclude primers that might have potential for off-target binding, in order to remove non-specific primer sets.

The `pdp blastscreen` subcommand expects the path to a suitable pre-compiled `BLASTN` database to be provided with the argument `--db` (there is a `BLASTN` *E. coli* genome database in the subdirectory `tests/walkthrough/blastdb/`). By specifiying `--outdir`, we place `BLAST` output in the `tests/walkthrough/blastn` subdirectory. The input config file use is the `deduped_primers.json` file produced in the last step, having locations of the non-redundant primer files, and we specify that the command writes a new config file called `screened.json` that points to a reduced set of primers, with the potentially off-target/non-specific sets excluded.

**No sequences are deleted as a result of this action.**

```bash
$ pdp blastscreen --db tests/walkthrough/blastdb/e_coli_screen.fna \
                     --outdir tests/walkthrough/blastn \
                     tests/walkthrough/deduped_primers.json \
                     tests/walkthrough/screened.json
```

This screen produces the new subdirectory `tests/walkthrough/blastn` that contains all the primer sequences derived from each genome in FASTA format, and the tabular format output (`.blasttab`) of each BLAST search. New files are added to the `deduped` subdirectory (with suffix `_screened`), describing the reduced sets of primers, post-screening.

```bash
$ tree tests/walkthrough/
tests/walkthrough/
├── blastdb
│   ├── e_coli_screen.fna.nhr
│   ├── e_coli_screen.fna.nin
│   └── e_coli_screen.fna.nsq
├── blastn
│   ├── GCF_000011605.1_primers.blasttab
│   ├── GCF_000011605.1_primers.fasta
│   ├── GCF_000291725.1_concat_noambig_primers.blasttab
│   ├── GCF_000291725.1_concat_noambig_primers.fasta
│   ├── GCF_000749845.1_concat_primers.blasttab
│   └── GCF_000749845.1_concat_primers.fasta
[…]
├── deduped
│   ├── GCF_000011605.1_named_deduped.bed
│   ├── GCF_000011605.1_named_deduped.json
│   ├── GCF_000011605.1_named_deduped_screened.bed
│   ├── GCF_000011605.1_named_deduped_screened.fasta
│   ├── GCF_000011605.1_named_deduped_screened.json
│   ├── GCF_000291725.1_concat_noambig_named_deduped.bed
│   ├── GCF_000291725.1_concat_noambig_named_deduped.json
│   ├── GCF_000291725.1_concat_noambig_named_deduped_screened.bed
│   ├── GCF_000291725.1_concat_noambig_named_deduped_screened.fasta
│   ├── GCF_000291725.1_concat_noambig_named_deduped_screened.json
│   ├── GCF_000749845.1_concat_named_deduped.bed
│   ├── GCF_000749845.1_concat_named_deduped.json
│   ├── GCF_000749845.1_concat_named_deduped_screened.bed
│   ├── GCF_000749845.1_concat_named_deduped_screened.fasta
│   └── GCF_000749845.1_concat_named_deduped_screened.json
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

### 7. Test primers against input sequences for crosshybridisation with `primersearch`

To identify which primers might be diagnostically useful for any of the classes/labels defined in the config file, we test how they potentially amplify the other genomes from the input set. We use the `EMBOSS` tool `primersearch` to do this.

Using the `pdp primersearch` subcommand, we pass the appropriate `.json` config file and specify a directory to hold `primersearch` output with `--outdir`, as well provide the path to write a new config file that will hold data about this crosshybridisation screen.

```bash
pdp primersearch \
       --outdir tests/walkthrough/primersearch \
       tests/walkthrough/screened.json \
       tests/walkthrough/primersearch.json
```

After running the command, a new directory `tests/walkthrough/primersearch` is produced, containing a new primer file (`.tab` - required for `primersearch` input) and `.primersearch` and `.json` files for each genome's `primersearch` results. In addition, there are `.json` and `.bed` files describing all amplicons produced *on* a specific genome: `_amplicons.bed` and `_amplicons.json`. The `.bed` files can be visualised using a tool like `IGV`.

```bash
$ tree tests/walkthrough
tests/walkthrough
[...]
├── primersearch
│   ├── Pba_SCRI1043_amplicons.bed
│   ├── Pba_SCRI1043_amplicons.json
│   ├── Pba_SCRI1043_primers.primertab
│   ├── Pba_SCRI1043_primersearch.json
│   ├── Pba_SCRI1043_ps_Pba_SCRI1043.primersearch
│   ├── Pba_SCRI1043_ps_Pbe_NCPPB_2795.primersearch
│   ├── Pba_SCRI1043_ps_Pwa_CFBP_3304.primersearch
│   ├── Pbe_NCPPB_2795_amplicons.bed
│   ├── Pbe_NCPPB_2795_amplicons.json
│   ├── Pbe_NCPPB_2795_primers.primertab
│   ├── Pbe_NCPPB_2795_primersearch.json
│   ├── Pbe_NCPPB_2795_ps_Pba_SCRI1043.primersearch
│   ├── Pbe_NCPPB_2795_ps_Pbe_NCPPB_2795.primersearch
│   ├── Pbe_NCPPB_2795_ps_Pwa_CFBP_3304.primersearch
│   ├── Pwa_CFBP_3304_amplicons.bed
│   ├── Pwa_CFBP_3304_amplicons.json
│   ├── Pwa_CFBP_3304_primers.primertab
│   ├── Pwa_CFBP_3304_primersearch.json
│   ├── Pwa_CFBP_3304_ps_Pba_SCRI1043.primersearch
│   ├── Pwa_CFBP_3304_ps_Pbe_NCPPB_2795.primersearch
│   ├── Pwa_CFBP_3304_ps_Pwa_CFBP_3304.primersearch
│   └── target_amplicons.json
├── primersearch.json
[...]
```

The new `primersearch.json` config file contains information about this crosshybridisation screen, and can be used for identification and extraction of diagnostic primer sequence sets.

### 8. Classify primers by predicted diagnostic capability with `classify`

To classify  primer sets by their ability to amplify only genomes belonging to a specific named group/label in the configuration file, we use the `pdp classify` subcommand. This examines the `primersearch` output and reports back primer sets that amplify all genomes having a specific label/group, and only those genomes.

We pass as input the `.json` file produced by the `pdp primersearch` run, and the path to a directory for the output of the `pdp classify` subcommand:

```bash
pdp classify \
       tests/walkthrough/primersearch.json \
       tests/walkthrough/classify
```

The new directory contains `.json` and `.ePrimer3` format files for each set of primers diagnostic to a given class, and summary information for thos primers in `summary.tab` and `results.json` files.

The output directory also contains `.json` and `.bed` files describing all regions on each genome that are amplified by primers having a particular specificity. For instance, the `Pbe_NCPPB_2795_Pectobacterium_amplicons.*` files contain information about the regions of the genome `Pbe_NCPPB_2795` that are amplified by all primer sets with specificity to the label/group `Pectobacterium`.

```bash
$ tree tests/walkthrough/classify
tests/walkthrough/classify
├── Pba_SCRI1043_Pectobacterium_amplicons.bed
├── Pba_SCRI1043_Pectobacterium_amplicons.json
├── Pba_SCRI1043_atrosepticum_NCBI_amplicons.bed
├── Pba_SCRI1043_atrosepticum_NCBI_amplicons.json
├── Pba_SCRI1043_gv1_amplicons.bed
├── Pba_SCRI1043_gv1_amplicons.json
├── Pbe_NCPPB_2795_Pectobacterium_amplicons.bed
├── Pbe_NCPPB_2795_Pectobacterium_amplicons.json
├── Pbe_NCPPB_2795_betavasculorum_NCBI_amplicons.bed
├── Pbe_NCPPB_2795_betavasculorum_NCBI_amplicons.json
├── Pbe_NCPPB_2795_gv7_amplicons.bed
├── Pbe_NCPPB_2795_gv7_amplicons.json
├── Pectobacterium_primers.ePrimer3
├── Pectobacterium_primers.json
├── Pwa_CFBP_3304_Pectobacterium_amplicons.bed
├── Pwa_CFBP_3304_Pectobacterium_amplicons.json
├── Pwa_CFBP_3304_gv2_amplicons.bed
├── Pwa_CFBP_3304_gv2_amplicons.json
├── Pwa_CFBP_3304_wasabiae_NCBI_amplicons.bed
├── Pwa_CFBP_3304_wasabiae_NCBI_amplicons.json
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
atrosepticum_NCBI	1	tests/walkthrough/classify/atrosepticum_NCBI_primers.json
betavasculorum_NCBI	2	tests/walkthrough/classify/betavasculorum_NCBI_primers.json
gv1	1	tests/walkthrough/classify/gv1_primers.json
gv2	2	tests/walkthrough/classify/gv2_primers.json
gv7	2	tests/walkthrough/classify/gv7_primers.json
```

### 9. Extract and assess candidate metabarcoding markers

The candidate diagnostic primer sets designed using `pdp` can be useful for generating group-specific metabarcoding markers. Using the `pdp extract` subcommand, the region(s) predicted to be amplified by each candidate primer set (in the `primersearch` step) are extracted from the target genomes, aligned (using `MAFFT`) and then compared to each other to estimate their sequence diversity.

Where a set of primers is predicted to amplify a specific subgroup of the input genomes, and the predicted amplicons are sequence-diverse, those regions and primers are potentially useful markers for metabarcoding. `pdp extract` produces an output table of diversity measures for each primer set to assess the sequence diversity of each primer set's predicted amplicons. 

We provide `pdp extract` with paths to the `JSON` config file describing `primersearch` results, and to a `JSON` output output from the `pdp classify` step, corresponding to one of the groups/labels that the resulting primers are specific to. We also supply the path to an output directory for the command to write to.

```bash
pdp extract \
    tests/walkthrough/primersearch.json \
    tests/walkthrough/classify/Pectobacterium_primers.json \
    tests/walkthrough/extract
```

The new output directory is created (if needed), and a further subdirectory named for primer specificity. Within that directory, all predicted amplicons for each group-specific primer are written as `FASTA` files, and a corresponding `.aln` alignment file (produced using `MAFFT`). An additional `distances_summary.tab` file is also produced, describing sequence diversity measures for amplicons produced by each primer set:

```bash
$ tree tests/walkthrough/extract
tests/walkthrough/extract
└── Pectobacterium_primers
    ├── GCF_000011605.1_primer_00001.aln
    ├── GCF_000011605.1_primer_00001.fasta
    ├── GCF_000011605.1_primer_00002.aln
    ├── GCF_000011605.1_primer_00002.fasta
    ├── GCF_000011605.1_primer_00003.aln
    ├── GCF_000011605.1_primer_00003.fasta
    ├── GCF_000749845.1_concat_primer_00003.aln
    ├── GCF_000749845.1_concat_primer_00003.fasta
    └── distances_summary.tab
```

The `distances_summary` file describes, for each primer set, sequence distance (min, max, mean and stdev), and diversity/evenness measures (Shannon Index and Shannon Evenness), to aid selection of candidate markers. It also provides a count of the number of unique amplicon sequences observed, and the count of sequences that are not unique (`unique + nonunique = total count`).

```bash
$ cat tests/walkthrough/extract/Pectobacterium_primers/distances_summary.tab 
primer	dist_mean	dist_sd	dist_min	dist_max	unique	nonunique	shannon_index	shannon_evenness
GCF_000011605.1_primer_00001	0.0267	0.0153	0.0100	0.0400	3	0	1.10	1.00
GCF_000011605.1_primer_00002	0.0600	0.0000	0.0600	0.0600	3	0	1.10	1.00
GCF_000011605.1_primer_00003	0.0400	0.0100	0.0300	0.0500	3	0	1.10	1.00
GCF_000749845.1_concat_primer_00003	0.1133	0.0115	0.1000	0.1200	3	0	1.10	1.00
```




## FURTHER INFORMATION

For further technical information, please read the comments contained within the top of each '*.py' file as well as the Supporting Information (['Methods S1' document](doi:10.1371/journal.pone.0034498.s006)) of [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498).

## CONTRIBUTORS

- [Leighton Pritchard](https://github.com/widdowquinn)
- [Benjamin Leopold](https://github.com/cometsong)
- [Michael Robeson](https://github.com/mikerobeson)
- [Rory McLeod](https://github.com/rory-mcleod)
