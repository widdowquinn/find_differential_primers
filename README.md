#README.md (find_differential_primers)

##Overview
This repository contains code for finding discriminatory primers among genomes or other biological sequences of interest. 

##DEPENDENCIES:
The following dependencies have been confirmed to work for running the 'find_differential_primers.py' pipeline, though any later version (and in some cases, some earlier versions) should work:

* **[Biopython](http://biopython.org/wiki/Download)**: (forthcoming) v1.63 - but **NOTE** you may need to clone the current Biopython repo from <https://github.com/biopython/biopython> (see below).
* **[bx-python](https://bitbucket.org/james_taylor/bx-python/src)**: current Hg (Mercurial) version works as of June 19, 2013
* **[EMBOSS](http://emboss.sourceforge.net/download/)** (for ePrimer3): v6.4.0, v6.6.0
* **[primer3](http://primer3.sourceforge.net/releases.php)**: v1.1.4 **NOTE:** primer3 version2 does not play nice with EMBOSS ePrimer3 (see e.g. [this](https://code.google.com/p/msatcommander/issues/detail?id=29), [this](https://bugs.launchpad.net/ubuntu/+source/emboss/+bug/978257) and [this](http://stackoverflow.com/questions/9866113/why-am-i-getting-this-error-in-my-primer3-eprimer3-mac-osx-build))
* **[prodigal](https://code.google.com/p/prodigal/)**  : v1.20
* **[BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)**: v2.2.22+, v2.2.28+, v2.2.29+

### NOTE ON EMBOSS/PRIMER3/BIOPYTHON DEPENDENCIES

There is a point of fragility in choice of EMBOSS, primer3, and Biopython versions, that centres around the following issues:

1. primer3 versions newer than v1.1.4 do not work with EMBOSS. This locks us, for now, into a 2008 version of primer3. As I chose to use the EMBOSS tools, I am taking their lead on this - when EMBOSS ePrimer3 changes to v2+, so will this script.
2. EMBOSS' ePrimer3 interface is not stable between version numbers. In particular, v6.5.0 changes the -otm flag to -opttm. This means that Biopython version v1.63 or lower does not use the appropriate option for EMBOSS v6.5.0. If you are using an EMBOSS version older than v6.5.0 then Biopython 1.63 should be fine. If you are using EMBOSS v6.5.0+, then note that the appropriate change has been committed at the `git` repository at <https://github.com/biopython/biopython> (as of Dec 2013), and for now you should install Biopython from the bleeding edge source. It is anticipated that this change will appear in Biopython v1.64.

#### Acceptable combinations:

* EMBOSS v6.5.0+/Biopython cloned from GitHub repository/Primer3 v1.1.4 **should work**
* EMBOSS pre-v6.5.0/Biopython pre-v1.63/Primer3 v1.1.4 **should work**
* Primer3 v2+: **will not work**

# INSTALLATION

If you have downloaded v0.1.0 or greater, and the dependencies above are satisfied, then installation should be as simple as cloning the repository:

```
$ git clone https://github.com/widdowquinn/find_differential_primers
$ cd find_differential_primers
```


then issuing:

```
$ python setup.py install
```

(or whatever variant you wish, e.g. for a home directory-local installation) from the top directory in the repository, with root permissions, if necessary.

#BASIC USE:

1. Collect all biological sequences to be distinguished into a convenient location (e.g. in same directory; this is not essential, but it simplifies things if using a `Makefile`). 
2. Construct a config file similar to the example given in `O104_primers_5.conf` or `test.conf`. This will describe each sequence by name, the classes to which it belongs, and (at least) the location of the FASTA file containing the sequence (or sequences - `find_differential_primers.py` will stitch sequences with the spacer `NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN`, if necessary).
3. If you need a BLAST database of negative screening examples, construct this with `makeblastdb` (part of **BLAST+**).
4. Run the `find_differential_primers.py` script, with suitable command-line options.

These steps are encapsulated in the accompanying `makefile` in `samples/makefile`. This file can be modified to point to your input sequence file of interest, and run by issuing `make` at the command-line. See documentation in the `makefile` for more details.


# TEST:

Change directory to `tests`, and run the script with the `test.conf` config file using default settings:

```
$ ../find_differential_primers/find_differential_primers.py -i test.conf -v
```

This should run to completion, and produce the output indicated below:

```
$ tree differential_primer_results/
differential_primer_results/
├── Erwinia_family-specific_amplicons.fas
├── Erwinia_family-specific_primers.eprimer3
├── Eta_1_99_specific_amplicons.fas
├── Eta_1_99_specific_primers.eprimer3
├── Pba_SCRI1043_specific_amplicons.fas
├── Pba_SCRI1043_specific_primers.eprimer3
├── Pca_PC1_specific_amplicons.fas
├── Pca_PC1_specific_primers.eprimer3
├── Pca_PCC21_specific_amplicons.fas
├── Pca_PCC21_specific_primers.eprimer3
├── Pectobacterium_family-specific_amplicons.fas
├── Pectobacterium_family-specific_primers.eprimer3
├── Pwa_WPP163_specific_amplicons.fas
├── Pwa_WPP163_specific_primers.eprimer3
├── atrosepticum_family-specific_amplicons.fas
├── atrosepticum_family-specific_primers.eprimer3
├── carotovorum_family-specific_amplicons.fas
├── carotovorum_family-specific_primers.eprimer3
├── differential_primer_results-families.tab
├── differential_primer_results.tab
├── tasmaniensis_family-specific_amplicons.fas
├── tasmaniensis_family-specific_primers.eprimer3
├── universal_amplicons.fas
├── universal_primers.eprimer3
├── wasabiae_family-specific_amplicons.fas
└── wasabiae_family-specific_primers.eprimer3
0 directories, 26 files
$ wc differential_primer_results/*
      51      68    2621 differential_primer_results/Erwinia_family-specific_amplicons.fas
     140     452    4296 differential_primer_results/Erwinia_family-specific_primers.eprimer3
      51      68    2621 differential_primer_results/Eta_1_99_specific_amplicons.fas
     140     469    4681 differential_primer_results/Eta_1_99_specific_primers.eprimer3
      24      32    1257 differential_primer_results/Pba_SCRI1043_specific_amplicons.fas
      68     226    2315 differential_primer_results/Pba_SCRI1043_specific_primers.eprimer3
      18      24     918 differential_primer_results/Pca_PC1_specific_amplicons.fas
      52     172    1737 differential_primer_results/Pca_PC1_specific_primers.eprimer3
      21      28    1085 differential_primer_results/Pca_PCC21_specific_amplicons.fas
      60     199    2019 differential_primer_results/Pca_PCC21_specific_primers.eprimer3
      42      56    2166 differential_primer_results/Pectobacterium_family-specific_amplicons.fas
     116     374    3583 differential_primer_results/Pectobacterium_family-specific_primers.eprimer3
      18      24     936 differential_primer_results/Pwa_WPP163_specific_amplicons.fas
      52     172    1759 differential_primer_results/Pwa_WPP163_specific_primers.eprimer3
      24      32    1257 differential_primer_results/atrosepticum_family-specific_amplicons.fas
      68     218    2138 differential_primer_results/atrosepticum_family-specific_primers.eprimer3
      24      32    1236 differential_primer_results/carotovorum_family-specific_amplicons.fas
      68     218    2107 differential_primer_results/carotovorum_family-specific_primers.eprimer3
      13      56    1150 differential_primer_results/differential_primer_results-families.tab
      15      86     933 differential_primer_results/differential_primer_results.tab
      51      68    2621 differential_primer_results/tasmaniensis_family-specific_amplicons.fas
     140     452    4301 differential_primer_results/tasmaniensis_family-specific_primers.eprimer3
       0       0       0 differential_primer_results/universal_amplicons.fas
       4      10     135 differential_primer_results/universal_primers.eprimer3
      18      24     936 differential_primer_results/wasabiae_family-specific_amplicons.fas
      52     166    1626 differential_primer_results/wasabiae_family-specific_primers.eprimer3
    1330    3726   50434 total
$ cat differential_primer_results/differential_primer_results.tab 
# Summary information table
# Generated by find_differential_primers
# Columns in the table:
# 1) Query organism ID
# 2) Query organism families
# 3) Count of organism-unique primers
# 4) Count of universal primers
# 5) Query sequence filename
# 6) Query feature filename
# 7) Query ePrimer3 primers filename
Pba_SCRI1043	Pectobacterium,atrosepticum	8	0	sequences/NC_004547.fna	sequences/NC_004547.prodigalout	sequences/NC_004547.eprimer3
Pca_PC1	Pectobacterium,carotovorum	6	1	sequences/NC_012917.fna	sequences/NC_012917.prodigalout	sequences/NC_012917.eprimer3
Pwa_WPP163	Pectobacterium,wasabiae	6	2	sequences/NC_013421.fna	sequences/NC_013421.prodigalout	sequences/NC_013421.eprimer3
Pca_PCC21	Pectobacterium,carotovorum	7	0	sequences/NC_018525.fna	sequences/NC_018525.prodigalout	sequences/NC_018525.eprimer3
Eta_1_99	Erwinia,tasmaniensis	17	0	sequences/NC_010694.fna	sequences/NC_010694.prodigalout	sequences/NC_010694.eprimer3
$ cat differential_primer_results/differential_primer_results-families.tab 
# Summary information table
# Generated by find_differential_primers
# Columns in the table:
# 1) Family
# 2) Count of family-specific primers
# 3) Family-specific primer file
# 4) Family-specific amplicon file
Erwinia	17	differential_primer_results/Erwinia_family-specific_primers.eprimer3	differential_primer_results/Erwinia_family-specific_amplicons.fas
carotovorum	8	differential_primer_results/carotovorum_family-specific_primers.eprimer3	differential_primer_results/carotovorum_family-specific_amplicons.fas
Pectobacterium	14	differential_primer_results/Pectobacterium_family-specific_primers.eprimer3	differential_primer_results/Pectobacterium_family-specific_amplicons.fas
wasabiae	6	differential_primer_results/wasabiae_family-specific_primers.eprimer3	differential_primer_results/wasabiae_family-specific_amplicons.fas
atrosepticum	8	differential_primer_results/atrosepticum_family-specific_primers.eprimer3	differential_primer_results/atrosepticum_family-specific_amplicons.fas
tasmaniensis	17	differential_primer_results/tasmaniensis_family-specific_primers.eprimer3	differential_primer_results/tasmaniensis_family-specific_amplicons.fas
```

Once you have run the tests once, you can use the `test_nocalc.conf` file to re-run the prediction/classification without having to execute `Prodigal`, `ePrimer3` or `PrimerSearch`, by running:

```
$ ../find_differential_primers/find_differential_primers.py -i test_nocalc.conf -v --noprimersearch --noprodigal --noprimer3
```

##FURTHER INFORMATION:
Please read the comments contained within the top of each '*.py' file as well as the Supporting Information (['Methods S1' document](doi:10.1371/journal.pone.0034498.s006)) of [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498).

##CONTRIBUTORS
* [Leighton Pritchard](https://github.com/widdowquinn)
* [Benjamin Leopold](https://github.com/cometsong)
* [Michael Robeson](https://github.com/mikerobeson)

##CITATIONS
Please refer to the following for methodological details:

* Pritchard L _et al._ (2012) "Alignment-Free 
Design of Highly Discriminatory Diagnostic Primer Sets for _Escherichia coli_ O104:H4 Outbreak Strains." _PLoS ONE_ **7**(4): e34498. [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498) - _Method description and application to human bacterial pathogens, sub-serotype resolution_
* Pritchard L _et al._ (2013) "Detection of phytopathogens of the genus _Dickeya_ using a PCR primer 
prediction pipeline for draft bacterial genome sequences." _Plant Pathology_, **62**, 587-596
[doi:10.1111/j.1365-3059.2012.02678.x](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-3059.2012.02678.x/full) - _Application to plant pathogens, species-level resolution_
