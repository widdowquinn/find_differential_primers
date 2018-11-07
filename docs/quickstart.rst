.. _pdp-quickstart:

================
QuickStart Guide
================

------------
Installation
------------

To use ``pdp`` you will need to install it on your machine (laptop, desktop, server, or cluster). Installation is easiest when using one of the two popular Python package managers:

^^^^^^^^^^^^
1. ``conda``
^^^^^^^^^^^^

``pdp`` is available through the ``bioconda`` channel of `Anaconda`_:

.. code-block:: bash

    conda install -c bioconda pdp

^^^^^^^^^^^
2. ``PyPI``
^^^^^^^^^^^

``pdp`` is available *via* the `PyPI`_ package manager for Python:

.. code-block:: bash

    pip install pdp

.. TIP::
    ``pdp`` can also be installed directly from source. More detailed installation instructions can be found on the :ref:`pdp-installation` page.


----------------------------------------------
``pdp`` Walkthrough for qPCR diagnostic design
----------------------------------------------

An example qPCR primer set design using ``pdp`` is provided as a walkthrough below. The general procedure for any ``pdp`` analysis is:

1. Create a configuration file that describes which input sequences will be used for design, the paths to the sequence files, and labels describing which diagnostic classes or groups the sequences belong to
2. Fix the input sequences (if necessary) by stitching multiple contigs together and removing ambiguity symbols
3. Design primers to each input sequence
4. Deduplicate primer designs (this reduces computation time)
5. Screen primers against a ``BLAST`` or ``DIAMOND`` database of off-target sequences
6. Perform *in silico* cross-hybridisation of primers against each input sequence
7. Classify each primer set by specificity to the input sequence group labels, to generate candidate diagnostic primer sets


.. TIP::
    To see options available for the ``pdp`` program, or any subcommand, use the ``-h`` or ``--help`` option, e.g.:

    .. code-block:: bash

        pdp -help
        pdp config -h


^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Create configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will work with three bacterial genomes in the ``tests/walkthrough/sequences`` subdirectory of the ``find_differential_primers`` repository. These genomes represent three different *Pectobacterium* species, and are provided as ``.fasta`` files containing complete genome sequences. We can see these files using the ``ls`` command:

.. code-block:: bash

    $ ls tests/walkthrough/sequences/
    GCF_000011605.1.fasta	GCF_000291725.1.fasta	GCF_000749845.1.fasta

To use these files in our analysis, we need to construct a configuration file. This is a plain text, tab-separated file containing columns describing the following data:

- The name of the sequence (for use during analysis)
- A comma-separated list of labels describing classes/groups to which the sequence belongs
- The path to the sequence's FASTA file
- (Optionally) the path to a set of features (e.g. gene annotations); if no features are supplied, this can be replaced with a hyphen: ``-``.

.. TIP::
    Any line in the configuration file that starts with a hash/octothorpe (``#``) is ignored. This is useful for commenting and/or ignoring specific sequences.

We will use a pre-prepared configuration file at the location ``tests/walkthrough/pectoconf.tab``. The contents are shown below (note the extensive commenting).

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


"""""""""""""""""""""""""""""""
Validate the configuration file
"""""""""""""""""""""""""""""""

To confirm that the configuration file can be used in the rest of the design process, use the command ``pdp config --validate`` on that file:

.. code-block:: bash

    $ pdp config --validate tests/walkthrough/pectoconf.tab
    WARNING: Validation problems
        Pbe_NCPPB_2795 requires stitch (tests/walkthrough/sequences/GCF_000749845.1.fasta)
        Pwa_CFBP_3304 requires stitch (tests/walkthrough/sequences/GCF_000291725.1.fasta)
        Pwa_CFBP_3304 has non-N ambiguities (tests/walkthrough/sequences/GCF_000291725.1.fasta)


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2. Prepare the input sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. ATTENTION::
    To generate diagnostic primers and metabarcoding markers, the input sequences must each be "stitched" so that there is only a single contiguous sequence corresponding to each input file. Also, any IUPAC ambiguity symbols (e.g. `W`, `Y`, etc.) must be replaced with `N`.

Instead of modifying the input sequence directly, which would change the source data (a problem for reproducibility), ``pdp`` will "fix" these sequences by writing new, "stitched" and "cleaned" versions of the input sequences, then automatically update the configuration file to point to the modified files. This is done with the ``pdp config --fix_sequences`` command:

.. code-block:: bash

    pdp config --fix_sequences tests/walkthrough/fixed.json tests/walkthrough/pectoconf.tab

.. ATTENTION::
    The ``--fix_sequences`` option takes as its argument the location to write the output configuration file; the final positional argument is the path to the input configuration file.

The output file is written in a different format to the input file, as a `JSON`_ file. This is a more machine-readable version of the configuration file, and all the other ``pdp`` tools accept ``JSON`` format configuration files. You don't need to be familiar with the details of the format, but for information the ``fixed.conf`` configuration file is shown below. Here, the ``Pba_SCRI1043`` is unmodified, the ``Pbe_NCPPB_2795`` is stitched (hence ``concat`` appears in the filename), and the ``Pwa_CFBP_3304`` is both stitched and has ambiguity symbols replaced (so has ``concat_noambig`` in the filename).

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

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
3. Design primers to each input sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we can design primer sets against each input sequence, using the `EMBOSS`_ package ``ePrimer3``.

.. WARNING::

    The `EMBOSS`_ ``ePrimer3`` package uses the `PRIMER3`_ primer design software, but will only work with an old version of that software: ``v1.1.4``

 To do this, we use the ``pdp eprimer3`` command, specifying with ``--outdir`` the output directory to write the predicted primers into. We also need to specify the input configuration file (the ``JSON`` file ``fixed.json`` produced in step 2), and the location of a new output ``JSON`` file (here, ``with_primers.json``) that associates primer information with the appropriate input sequence, in each case:

.. code-block:: bash

    pdp eprimer3 --outdir tests/walkthrough/eprimer3 \
        tests/walkthrough/fixed.json
        tests/walkthrough/with_primers.json

The new ``tests/walkthrough/eprimer3`` directory now contains files describing primers designed to each input sequence, and corresponding ``JSON`` files describing the primer sets.

.. code-block:: bash

    $ tree tests/walkthrough/eprimer3/
    tests/walkthrough/eprimer3/
    ├── GCF_000011605.1.eprimer3
    ├── GCF_000011605.1_named.eprimer3
    ├── GCF_000011605.1_named.json
    ├── GCF_000291725.1_concat_noambig.eprimer3
    ├── GCF_000291725.1_concat_noambig_named.eprimer3
    ├── GCF_000291725.1_concat_noambig_named.json
    ├── GCF_000749845.1_concat.eprimer3
    ├── GCF_000749845.1_concat_named.eprimer3
    └── GCF_000749845.1_concat_named.json

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
4. Deduplicate primer sets (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. ATTENTION::

    This step is recommended, but not necessary, when designing diagnostic primer sets

When designing primers to groups of closely-related genomes, it is usual to have a large number of identical primer sets that originate from different genomes. We only need to test one of these redundant primer sets to know whether it may be diagnostically useful, so we can remove duplicates with the ``pdp dedupe`` command:

.. code-block:: bash

    pdp dedupe --dedupedir tests/walkthrough/deduped \
        tests/walkthrough/with_primers.json \
        tests/walkthrough/deduped_primers.json

The complete set of nonredundant primers is written to ``tests/walkthrough/deduped``, and a new ``JSON`` configuration file recording only the deduplicated primers for each input sequence is written to ``deduped_primers.json``.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
5. Screen primers against a local sequence database (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. ATTENTION::

    This step is recommended, but not necessary, when designing diagnostic primer sets

Prescreening the primers we have just designed against a local database of off-target sequences allows us to remove primer sets that do not specifically amplify our input sequences without having to perform computationally costly *in silico* cross-hybridisation.

.. TIP::

    The composition of the screening database should be appropriate to your analysis/design goals. For example, if you are interested in designing primers diagnostic to all species in a particular bacterial genus, then a database comprising available genomes from sister genera may be appropriate. Alternatively, a subsampling of complete genomes from the bacterial group containing your genus of interest may be useful. However the screening database is constructed, it should represent a good range of off-target sequences that could reasonably be detected as false positives, to eliminate non-specific primer sets.

To use a local ``BLAST`` nucleotide database as the off-target screen, use the ``pdp blastscreen`` command, provide the location of the database with the ``--db`` argument, and the location to write output results with ``--outdir``. In the command below, we use a precompiled database of *Escherichia coli* genomes, and use the ``deduped_primers.json`` configuration file as input. A new configuration file recording only the screened primers for each sequence is written to ``screened.json``.

.. code-block:: bash

    pdp blastscreen --db tests/walkthrough/blastdb/e_coli_screen.fna \
        --outdir tests/walkthrough/blastn \
        tests/walkthrough/deduped_primers.json \
        tests/walkthrough/screened.json

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
6. Perform *in silico* cross-hybridisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the critical step in determining the predicted diagnostic specificity of the candidate primer sets.

.. TIP::

    It is strongly recommended that primers are deduplicated, and an off-target pre-screen is performed using ``pdp blastscreen`` or ``pdp diamondscreen`` before carrying out this step.

In this step, each candidate primer set is tested in turn against all the input sequences to determine whether it has the potential to amplify that sequence. This is the most computationally-demanding step of the analysis.

The ``pdp primersearch`` command uses the `EMBOSS`_ tool ``primersearch`` to carry out *in silico* hybridisation of each of the candidate primer sets. The output directory into which result files are written is specified with ``--outdir``, and here we use the configuration file of pre-screened primers ``screened.json`` as input, writing a new configuration file (that records the *in silico* hybridisation results) as ``primersearch.json``:

.. code-block:: bash

    pdp primersearch \
        --outdir tests/walkthrough/primersearch \
        tests/walkthrough/screened.json \
        tests/walkthrough/primersearch.json

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
7. Classify primer sets by specificity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The final step in determining qPCR primer set specificity is to analyse the *in silico* hybridisation results to determing which primer sets amplify exactly the members of each class/group defined in the initial configuration file. This is done using the ``pdp classify`` subcommand, which takes as arguments the configuration file written out by the *in silico* hybridisation step (``primersearch.json``), and the location of a directory into which the output will be written (here, ``tests/walkthrough/classify``):

.. code-block:: bash

    pdp classify \
        tests/walkthrough/primersearch.json \
        tests/walkthrough/classify

The output directory contains ``.json`` and ``.ePrimer3`` format files for each set of candidate primers that were determined to be specific to a class/group named in the initial configuration file, and two summary files (``results.json`` and ``summary.tab``):

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

The ``summary.tab`` file is a tab-separated plain text file that describes how many primer sets were determined to potentially be diagnostic for each input class, and describes a path to the ``JSON`` file describing their results:

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


.. _Anaconda: https://www.anaconda.com/what-is-anaconda/
.. _EMBOSS: http://emboss.sourceforge.net/
.. _JSON: https://www.json.org/
.. _PRIMER3: http://primer3.sourceforge.net/
.. _PyPI: https://pypi.python.org/pypi
