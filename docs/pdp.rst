.. _pdp-pdp:

==============================
``pdp`` command-line interface
==============================

The ``pdp`` command-line interface is intended to split individual stages of the complete primer/metabarcoding marker design process into logically-sensible "chunks". Each stage is invoked using the formula:

.. code-block:: bash

    pdp <SUBCOMMAND> <OPTIONS> <ARGUMENTS>

The stages are described and summarised below, in the usual order of implementation for design of qPCR primers or metabarcoding markers.

----------------------------
0. Create configuration file
----------------------------

The ``pdp`` software takes as initial input a configuration file. This is a plain text, tab-separated file containing columns describing the following data:

- The name of the sequence (for use during analysis)
- A comma-separated list of labels describing classes/groups to which the sequence belongs
- The path to the sequence's FASTA file
- (Optionally) the path to a set of features (e.g. gene annotations); if no features are supplied, this can be replaced with a hyphen: ``-``.

.. TIP::
    Any line in the configuration file that starts with a hash/octothorpe (``#``) is ignored. This is useful for commenting and/or ignoring specific sequences.

An example configuration file is included at the location ``tests/walkthrough/pectoconf.tab``. The contents are shown below (note the extensive commenting, with lines preceded by a hash/octothorpe).

.. code-block:: bash

    $ cat tests/walkthrough/pectoconf.tab
    # Pectobacterium genomes downloaded from GenBank/NCBI; genomovars inferred from ANIm
    # Annotated Pba: genomovar 1
    Pba_SCRI1043	Pectobacterium,atrosepticum_NCBI,gv1	tests/walkthrough/sequences/GCF_000011605.1.fasta	-
    # Annotated Pwa: genomovars 2, 3
    Pwa_CFBP_3304	Pectobacterium,wasabiae_NCBI,gv2	tests/walkthrough/sequences/GCF_000291725.1.fasta	-
    # Annotated Pb	: genomovar 7
    Pbe_NCPPB_2795	Pectobacterium,betavasculorum_NCBI,gv7	tests/walkthrough/sequences/GCF_000749845.1.fasta	-

The first line describing an input sequence tells us that its name is ``Pba_SCRI1043``, that it belongs to classes/groups ``Pectobacterium``, ``atrosepticum_NCBI``, and ``gv1``, and that the sequence's FASTA file can be found at ``tests/walkthrough/sequences/GCF_000011605.1.fasta``. There are no features associated with the sequence.

.. TIP::
    It can be easiest to prepare your configuration files in a spreadsheet program such as Microsoft Excel, then convert the file to ``.json`` format with ``pdp config``.


-----------------
1. ``pdp config``
-----------------

The ``pdp config`` command provides functions for validation and handling of ``pdp`` configuration files. The key functions provided by this command are:

validation
    ``pdp config`` can check the formatting and content of a configuration file, to confirm that the appropriate fields are all populated for each input sequence, and that the specified files exist on the filesystem. The program will warn if input sequences require some repair.

.. code-block:: bash

    $ pdp config --validate tests/walkthrough/pectoconf.tab
    WARNING: Validation problems
        Pbe_NCPPB_2795 requires stitch (tests/walkthrough/sequences/GCF_000749845.1.fasta)
        Pwa_CFBP_3304 requires stitch (tests/walkthrough/sequences/GCF_000291725.1.fasta)
        Pwa_CFBP_3304 has non-N ambiguities (tests/walkthrough/sequences/GCF_000291725.1.fasta)

sequence repair
    Sequences used as input to ``pdp`` may not conform to requirements of some of the third party tools, or the pipeline's other assumptions. Some stages in the pipeline cannot accommodate genomes presented as multiple sequence files, and others cannot cope with nucleotide ambiguity symbols other than ``N``. ``pdp config`` can repair input genome sequences by stitching them and replacing ambiguity codons, writing the "fixed" genome to a new file. Instead of modifying the input sequence directly, which would change the source data (a problem for reproducibility), ``pdp`` will "fix" these sequences by writing new, "stitched" and "cleaned" versions of the input sequences, then automatically update the configuration file to point to the modified files. This is done with the ``pdp config --fix_sequences`` command:

.. code-block:: bash

    pdp config --fix_sequences tests/walkthrough/fixed.json tests/walkthrough/pectoconf.tab

.. ATTENTION::
    The ``--fix_sequences`` option takes as its argument the location to write the output configuration file; the final positional argument is the path to the input configuration file.

format conversion
    The ``pdp config`` subcommand can convert configuration files between ``.tab`` and `JSON`_ format. The ``.tab`` format is easier to read and manipulate in spreadsheet software, but all the tools in the ``pdp`` pipeline require input in ``.json`` format.

.. code-block:: bash

    pdp config --to_json myconfig.json myconfig.tab
    pdp config --to_tab myconfig.tab myconfig.json


``pdp config`` output files are (except when converting to ``.tab`` format) written as `JSON`_ files. This is a machine-readable version of the configuration file, and all the other ``pdp`` tools accept ``JSON`` format configuration files, but not ``.tab`` files. As an example the ``fixed.conf`` configuration file produced in the example above is shown below. Here, the ``Pba_SCRI1043`` line is unmodified, the ``Pbe_NCPPB_2795`` is stitched (hence ``concat`` appears in the filename), and the ``Pwa_CFBP_3304`` is both stitched and has ambiguity symbols replaced (so has ``concat_noambig`` in the filename).

.. code-block:: json

    [
        {
            "features": null,
            "filestem": "GCF_000011605.1",
            "filtered_seqfile": null,
            "groups": [
                "Pectobacterium",
                "atrosepticum_NCBI",
                "gv1"
            ],
            "name": "Pba_SCRI1043",
            "primers": null,
            "primersearch": null,
            "seqfile": "tests/walkthrough/sequences/GCF_000011605.1.fasta"
        },
        {
            "features": null,
            "filestem": "GCF_000749845.1_concat",
            "filtered_seqfile": null,
            "groups": [
                "Pectobacterium",
                "betavasculorum_NCBI",
                "gv7"
            ],
            "name": "Pbe_NCPPB_2795",
            "primers": null,
            "primersearch": null,
            "seqfile": "tests/walkthrough/sequences/GCF_000749845.1_concat.fas"
        },
        {
            "features": null,
            "filestem": "GCF_000291725.1_concat_noambig",
            "filtered_seqfile": null,
            "groups": [
                "Pectobacterium",
                "gv2",
                "wasabiae_NCBI"
            ],
            "name": "Pwa_CFBP_3304",
            "primers": null,
            "primersearch": null,
            "seqfile": "tests/walkthrough/sequences/GCF_000291725.1_concat_noambig.fas"
        }
    ]

As can be seen from this file, the modified sequences are written to the ``tests/walkthrough/sequences`` subdirectory, alongside the original:

.. code-block:: bash

    $ tree tests/walkthrough/sequences/
    tests/walkthrough/sequences/
    ├── GCF_000011605.1.fasta
    ├── GCF_000291725.1.fasta
    ├── GCF_000291725.1_concat.fas
    ├── GCF_000291725.1_concat_noambig.fas
    ├── GCF_000749845.1.fasta
    └── GCF_000749845.1_concat.fas

-----------------
2. ``pdp filter``
-----------------

As an optional step in the primer/marker design pipeline, input genomes may be filtered such that primer sets are only designed to restricted sections of the input. The ``pdp filter`` subcommand has options to automate this process for the following genomic regions:

predicted coding sequences
    The ``pdp filter --prodigal`` option carries out *a priori* CDS prediction on each input genome using the `Prodigal`_ software tool. It generates a ``.gff`` file, and includes that in the output configuration file. This is suggested as an approach to maximise conserved genomic regions as sources for robust qPCR diagnostic primer design, as coding sequences are expected to be relatively conserved.

.. code-block:: bash

    pdp filter --prodigal myconfig.json cds_filtered.json

predicted intergenic regions
    Using the ``pdp filter --prodigaligr`` option will conduct CDS prediction with `Prodigal`_, but generate a ``.gff`` file describing regions *between* predicted CDS, plus a short flanking regions that extends into neighbouring predicted CDS. This is included in the output configuration file. The size of the flanking region can be controlled with the ``--flanklen`` option. This is suggested as an approach to maximise amplicon sequence variation when designing metabarcoding markers, as intergenic regions are expected to be less well-conserved.

.. code-block:: bash

    pdp filter --prodigaligr myconfig.json igr_filtered.json

consserved variable regions
    With the ``pdp filter --alnvar <GROUP>`` option, ``pdp`` will pairwise align all genomes annotated in the configuration file with group ``<GROUP>``, to identify regions conserved across all members of that group. If the ``--min_sim_error_rate <VALUE>`` option is used, then only conserved regions that have an error rate (SNP or indel) of at least ``<VALUE>``% (constrained in the range [0, 1) in *every* genome of ``<GROUP>`` will be retained for primer design. This is suggested as an approach to maximising amplicon sequence variatino for metabarcoding marker design, as a minimum level of source sequence variation is assured.

.. code-block:: bash

    pdp filter --alnvar group01 --min_sim_error_rate 0.1 myconfig.json igr_filtered.json

.. TIP::
    The ``pdp filter`` subcommand can be used with options for multiprocessing/`SGE`_-like parallelisation (see below)

-------------------
3. ``pdp eprimer3``
-------------------

.. WARNING::

    The `EMBOSS`_ ``ePrimer3`` package uses the `PRIMER3`_ primer design software, but will only work with an old version of that software: ``v1.1.4``

The ``pdp eprimer3`` command takes a ``.json`` format configuration file, and uses the  `EMBOSS`_ ``ePrimer3`` package to design primers to each input genome sequence, writing the output files to a location specified with the ``--outdir <OUTDIR>`` option. A new output ``.json`` configuration file is produced which associates primer information with the appropriate input sequence.

.. code-block:: bash

    pdp eprimer3 --outdir primers/ myconfig.json genomes_with_primers.json

In order to use the "filtered" genomes specified in an input configuration file, the ``pdp eprimer3`` subcommand must be run with the ``--filter`` option.

.. code-block:: bash

    pdp eprimer3 --filter --outdir primers/ myconfig.json genomes_with_primers.json

.. TIP::
    The ``pdp eprimer3`` subcommand can be used with options for multiprocessing/`SGE`_-like parallelisation (see below)

-----------------
4. ``pdp dedupe``
-----------------

The ``pdp eprimer3`` primer design step treats each genome independently which, especially when several related genomes are used as input, can result in design of many identical primer sets. The ``pdp dedupe`` subcommand identifies these redundant primer sets and discards them from the design process, generating a new ``.json`` configuration file that points only to non-redundant primer sets. The location of deduplicated/non-redundant primer sets generated in this process is given with the ``--dedupedir <DIRNAME>`` option, and ``BLASTN+`` output written to a location specified by the ``--outdir <OUTPUTDIR>`` option.

.. TIP::
    This step is optional, but it is typical that most primer sets in the primer design process are redundant. As some of the subsequent pipeline stages do not scale well with an increasing number of primer sets, deduplication is *highly recommended*.

.. code-block:: bash

    pdp dedupe --dedupedir deduped/ myconfig.json deduped.json

----------------------
5. ``pdp blastscreen``
----------------------

As a fast, preliminary screen of designed primers against a set of off-target sequences, the ``pdp blastscreen`` subcommand performs ``BLASTN+`` searches of designed primers against a pre-prepared nucleotide sequence database specified by the ``--db <DBPATH>`` option, where ``<DBPATH>`` is the path to the ``BLAST+`` database.

.. code-block:: bash

    pdp blastscreen --db blastdb/offtargets --outdir blastn_outputmyconfig.json screened.json

.. TIP::
    The composition of the screening database should be appropriate to your analysis/design goals. For example, if you are interested in designing primers diagnostic to all species in a particular bacterial genus, then a database comprising available genomes from sister genera may be appropriate. Alternatively, a subsampling of complete genomes from the bacterial group containing your genus of interest may be useful. However the screening database is constructed, it should represent a good range of off-target sequences that could reasonably be detected as false positives, to eliminate non-specific primer sets.

.. TIP::
    This step is optional, but it is typical that many primer sets in the primer design process have off-target matches and are not useful as diagnostic tools. As some of the subsequent pipeline stages do not scale well with an increasing number of primer sets, screening against a relevant database is *highly recommended*.

.. TIP::
    The ``pdp blastscreen`` subcommand can be used with options for multiprocessing/`SGE`_-like parallelisation (see below)

-----------------------
6. ``pdp primersearch``
-----------------------

``pdp primersearch`` carries out the most critical step in the design process: predicting which genomes produce amplicons for each of the designed primer sets. Each candidate primer set is tested in turn against all the input sequences to determine whether it has the potential to amplify that sequence. This is the most computationally-demanding step of the analysis.

The ``pdp primersearch`` command uses the `EMBOSS`_ tool ``primersearch`` to carry out *in silico* hybridisation of each of the candidate primer sets against each input (*unfiltered*) genome. The output directory into which ``primersearch`` result files are written can be specified with ``--outdir``:

.. code-block:: bash

    pdp primersearch --outdir primersearch_results myconfig.json primersearch.json

.. TIP::
    It is strongly recommended that primers are deduplicated, and an off-target pre-screen is performed using ``pdp blastscreen`` or ``pdp diamondscreen`` before carrying out this step.

.. TIP::
    The ``pdp primersearch`` subcommand can be used with options for multiprocessing/`SGE`_-like parallelisation (see below)

-------------------
7. ``pdp classify``
-------------------

Each primer set can be assessed for potential specificity using the ``pdp classify`` subcommand to determine whether it is predicted to amplify specifically genomes only from one of the groups specified in the configuration file.

.. code-block:: bash

    pdp classify myconfig.json classified_primers/

No new configuration file is produced, but new ``.json`` and ``.ePrimer3`` format files are written to the specified output directory for each of the groups specified in the configuration file. These contain accounts of the primer sets determined to be specific to that group.

.. code-block:: bash

    $ tree tests/walkthrough/classify/
    tests/walkthrough/classify/
    ├── Pectobacterium_primers.ePrimer3
    ├── Pectobacterium_primers.json
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

Two summary files are also written: ``results.json`` and ``summary.tab``. The ``summary.tab`` file is a tab-separated plain text file that describes how many primer sets were predicted to be diagnostic for each input class, and describes a path to the ``.json`` file describing their results:

.. code-block:: bash

    $ cat tests/walkthrough/classify/summary.tab
    Group   NumPrimers      Primers
    Pectobacterium  4       tests/walkthrough/classify/Pectobacterium_primers.json
    atrosepticum_NCBI       1       tests/walkthrough/classify/atrosepticum_NCBI_primers.json
    betavasculorum_NCBI     2       tests/walkthrough/classify/betavasculorum_NCBI_primers.json
    gv1     1       tests/walkthrough/classify/gv1_primers.json
    gv2     2       tests/walkthrough/classify/gv2_primers.json
    gv7     2       tests/walkthrough/classify/gv7_primers.json
    wasabiae_NCBI   2       tests/walkthrough/classify/wasabiae_NCBI_primers.json

------------------
8. ``pdp extract``
------------------

The ``pdp extract`` subcommand extracts the amplicon sequences for each of the primers in a specified ``pdp classify`` output primer file (specific for a particular genome group), aligns those sequences using ``MAFFT`` and then produces summary statistics about the sequence diversity of the amplicons.

``pdp extract`` requires the configuration file used for ``pdp classify``, the path to the ``.json`` file describing the primers specific to a particular genome group, and the path to an output directory for result files. Given the output path ``<OUTDIR>`` and the group ``group01``, the extracted sequences and summary files will be written to ``<OUTDIR>/group01``.

.. code-block:: bash

    pdp extract myconfig.json classified_primers/group01.json extracted/

The output directory will contain, for each primer set:

- a FASTA format file describing all the amplicon sequences (ending in ``.fasta``)
- a ``MAFFT``-aligned output file (ending in ``.aln``) describing the aligned sequences
- a summary tab-separated plain text file called ``distances_summary.tab``

The ``distances_summary.tab`` file is a table with one row per primer set, (and one row for headers), describing:

- primer set identifier
- mean pairwise distance between aligned sequences
- standard deviation of pairwise distance between aligned sequences
- minimum pairwise distance between aligned sequences
- maximum pairwise distance between aligned sequences
- the count of unique amplicons
- the number of non-unique amplicons
- the Shannon Index of sequence diversity (larger is more diverse)
- the Shannon Evenness of sequence diversity ([0, 1]: closer to 1 is more even)

.. TIP::
    The ``pdp extract`` subcommand can be used with options for multiprocessing/`SGE`_-like parallelisation (see below)

----------------------------------------
Multiprocessing/SGE-like parallelisation
----------------------------------------

.. ATTENTION::
    Several stages in the ``pdp`` pipeline (principally those that call third-party software tools) can take advantage of multicore systems or clusters using an `SGE`_-like scheduler. The syntax for doing this is the same for each of the stages.

multiprocessing
    By default, subcommands that can use parallelism will attempt to distribute jobs to local cores using Python's built-in ``multiprocessing`` module. This can be explicitly enabled with the ``-s multiprocessing`` option, and the number of workers controlled with the ``-w <N>`` option, to limit the total number of workers to a maximum of ``<N>``.

.. code-block:: bash

    pdp eprimer3 --outdir primers -s multiprocessing -w 4 myconfig.json eprimer3.json

`SGE`_-like schedulers
    The ``-s SGE`` option can be provided to use an `SGE`_-like scheduler (one you can invoke with ``qsub``). To cause minimal problems with queues, individual jobs are batched into job arrays, with a default array size of 10000 jobs (this can be controlled with the ``--SGEgroupsize <N>`` option). If you need to pass further arguments to SGE, this can be done with the ``--SGEargs <ARGUMENTS>`` option.

.. code-block:: bash

    pdp eprimer3 --outdir primers -s SGE --SGEgroupsize 5000 --SGEargs "-M me@domain.org -m bes" myconfig.json eprimer3.json



.. _EMBOSS: http://emboss.sourceforge.net/
.. _JSON: https://www.json.org/
.. _PRIMER3: http://primer3.sourceforge.net/
.. _Prodigal: https://github.com/hyattpd/Prodigal
.. _SGE: https://en.wikipedia.org/wiki/Oracle_Grid_Engine
