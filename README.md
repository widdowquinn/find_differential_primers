# README.md (diagnostic_primers)

## Table of Contents

1. [Note for Users](#usernote)
1. [Note for Developers](#devnotenote)
1. [Overview](#overview)
1. [Usage](#usage)
      1. [Summary](#summary)
      2. [`pdp.py config`](#config)

## NOTE FOR USERS<a id="usernote"></a>

The default branch for this repository is a development branch: `diagnostic_primers`. If you are looking for code to reproduce work in Pritchard *et al.* (2012) or Pritchard *et al.* (2013), please checkout the `master` branch, or download [release v0.1.3](https://github.com/widdowquinn/find_differential_primers/tree/v0.1.3).

* `diagnostic_primers`: 

[![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/diagnostic_primers/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)
[![Build Status](https://travis-ci.org/widdowquinn/find_differential_primers.svg?branch=diagnostic_primers)](https://travis-ci.org/widdowquinn/find_differential_primers)

* `master`: 

[![codecov](https://codecov.io/gh/widdowquinn/find_differential_primers/branch/master/graph/badge.svg)](https://codecov.io/gh/widdowquinn/find_differential_primers)
[![Build Status](https://travis-ci.org/widdowquinn/find_differential_primers.svg?branch=master)](https://travis-ci.org/widdowquinn/find_differential_primers)

## NOTE FOR DEVELOPERS<a id="devnotenote"></a>

The default master branch for development is `diagnostic_primers`. We would appreciate contributions, especially if you follow the guidelines on the [wiki](https://github.com/widdowquinn/find_differential_primers/wiki).

* Current test coverage (`diagnostic_primers`): [https://codecov.io/gh/widdowquinn/find_differential_primers/list/diagnostic_primers](https://codecov.io/gh/widdowquinn/find_differential_primers/list/diagnostic_primers)

## Overview<a id="overview"></a>
This repository contains code for automated finding of discriminatory (real-time) PCR or qPCR primers that distinguish among genomes or other biological sequences of interest. 

## Usage<a id="usage"></a>

### Summary<a id="summary"></a>

This new version of `diagnostic_primers` (formerly `find_differential_primers`) now uses a subcommand model, like the tools `git` and `subversion`. These execute the following subtasks, some or all of which may be required in a primer design run.

* `config`: Process/validate the configuration file and stitch input contig fragments/replace ambiguity symbols as necessary.
* `prodigal`: Predict CDS locations on the input sequence
* `eprimer3`: Design amplifying primers on the input sequence
* `blastcheck`: Filter designed primers against a database of negative examples
* `primersearch`: Filter designed primers on their ability to amplify each input sequence
* `classify`: Classify designed primers by specificity for each class of input sequence

Each of these subcommands has specific help, accessible with `pdp.py <subcommand> -h` or `pdp.py <subcommand> --help`.

### `pdp.py config`<a id="config"></a>

The `config` subcommand handles interactions with the configuration file for a primer design run. Configuration files can be provided in one of two formats:

1. `<config>.tab`: a plain text, tab-separated file descrbing the input data in a table of multiple columns [as defined below](#config_tab_format). This is intended to be an easily human-readable file, that can be prepared and edited in a spreadsheet application such as Google Sheets, or Microsoft Excel. `pdp.py config` will recognise `.tab` or `.conf` as a file extension.
2. `<config>.json`: a JSON format file describing the input data. This is not intended to be human-readable, but can be converted to and from the `.tab` format using `pdp.py` config. Further steps in the primer design process require that the configuration is provided in JSON format.

#### Converting between `.tab` and JSON format config files

**1. `.tab` to JSON**

Provide the path to the output JSON file as an argument to `--to_json`, and the path to the `.tab` config file as input:

```bash
pdp.py config --to_json <OUTPUT>.json <INPUT>.tab
```

**2. JSON to `tab`**

Provide the path to the output `.tab` file as an argument to `--to_tab`, and the path to the JSON config file as input:

```bash
pdp.py config --to_tab <OUTPUT>.tab <INPUT>.json
```

## FURTHER INFORMATION:
Please read the comments contained within the top of each '*.py' file as well as the Supporting Information (['Methods S1' document](doi:10.1371/journal.pone.0034498.s006)) of [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498).

## CONTRIBUTORS
* [Leighton Pritchard](https://github.com/widdowquinn)
* [Benjamin Leopold](https://github.com/cometsong)
* [Michael Robeson](https://github.com/mikerobeson)
* [Rory McLeod](https://github.com/rory-mcleod)

## CITATIONS
Please refer to the following for methodological details:

* Pritchard L _et al._ (2012) "Alignment-Free 
Design of Highly Discriminatory Diagnostic Primer Sets for _Escherichia coli_ O104:H4 Outbreak Strains." _PLoS ONE_ **7**(4): e34498. [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498) - _Method description and application to human bacterial pathogens, sub-serotype resolution_
* Pritchard L _et al._ (2013) "Detection of phytopathogens of the genus _Dickeya_ using a PCR primer 
prediction pipeline for draft bacterial genome sequences." _Plant Pathology_, **62**, 587-596
[doi:10.1111/j.1365-3059.2012.02678.x](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-3059.2012.02678.x/full) - _Application to plant pathogens, species-level resolution_
