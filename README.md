#README.md (find_differential_primers)

##Overview
This repository contains code for finding discriminatory primers among genomes or other biological sequences of interest. 

##DEPENDENCIES:
The following Dependencies have been confirmed to work for running the 'find_differential_primers.py' pipeline, though any later version (and many earlier versions) should work:

* **[Biopython](http://biopython.org/wiki/Download)**: v1.57 onwards
* **[bx-python](https://bitbucket.org/james_taylor/bx-python/src)**: current Hg (Mercurial) version works as of June 19, 2013
* **[EMBOSS](http://emboss.sourceforge.net/download/)** (for ePrimer3): v6.30, v6.31, v6.4.0
* **[primer3](http://primer3.sourceforge.net/releases.php)**: v1.1.4
* **[prodigal](https://code.google.com/p/prodigal/)**  : v1.20
* **[BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)**: v2.2.22+ 



##CITATIONS
Please refer to the following for methodological details:

* Pritchard L _et al._ (2012) "Alignment-Free 
Design of Highly Discriminatory Diagnostic Primer Sets for _Escherichia coli_ O104:H4 Outbreak Strains." _PLoS ONE_ **7**(4): e34498. [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498) - _Method description and application to human bacterial pathogens, sub-serotype resolution_
* Pritchard L _et al._ (2013) "Detection of phytopathogens of the genus _Dickeya_ using a PCR primer 
prediction pipeline for draft bacterial genome sequences." _Plant Pathology_, **62**, 587-596
[doi:10.1111/j.1365-3059.2012.02678.x](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-3059.2012.02678.x/full) - _Application to plant pathogens, species-level resolution_


#BASIC USE:

1. Collect all biological sequences to be distinguished into a convenient location (e.g. in same directory; this is not essential, but it simplifies things if using a `Makefile`). 
2. Construct a config file similar to the example given in `O104_primers_5.conf`. This will describe each sequence by name, the classes to which it belongs, and (at least) the location of the FASTA file containing the sequence (or sequences - `find_differential_primers.py` will stitch sequences with a spacer, if necessary).
3. If you need a BLAST database of negative screening examples, construct this with `makeblastdb` (part of BLAST+).
4. Run the `find_differential_primers.py` script, with suitable command-line options.

Finally, run `find_differential_primers.py`.

##FURTHER INFORMATION:
Please read the comments contained within the top of each '*.py' file as well as the Supporting Information (['Methods S1' document](doi:10.1371/journal.pone.0034498.s006)) of [doi:10.1371/journal.pone.0034498](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0034498).

#CONTRIBUTORS
* [Leighton Pritchard](https://github.com/widdowquinn)
* [Benjamin Leopold](https://github.com/cometsong)
* [Michael Robeson](https://github.com/mikerobeson)