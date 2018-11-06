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

To see options available for the ``pdp`` program, use the ``-h`` or ``--help`` option:

.. code-block:: bash

    pdp -help


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

Any line in the configuration file that starts with a hash/octothorpe (``#``) is ignored. This is useful for commenting and/or ignoring specific sequences.

We will use a pre-prepared configuration file at the location ``tests/walkthrough/pectoconf.tab`` - note the extensive commenting.

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


^^^^^^^^^^^^^^^^^^^^^^^^^^
2. Fix the input sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^

To generate diagnostic primers and metabarcoding markers, the input sequences must each be "stitched" so that there is only a single contiguous sequence corresponding to each input file. Also, any IUPAC ambiguity symbols (e.g. `W`, `Y`, etc.) must be replaced with `N`.

Instead of modifying the input sequence directly, which would modify the source data, ``pdp`` will "fix" these sequences by writing new, "stitched" and "cleaned" versions of the input sequences, then automatically update the configuration file to point to the modified files. This is done with the ``pdp config --fix_sequences`` command:

.. code-block:: bash

    pdp config --fix_sequences tests/walkthrough/fixed.json tests/walkthrough/pectoconf.tab

The ``--fix_sequences`` option takes as the next argument the location to write the output configuration file; the final positional argument is the path to the input configuration file.

The output file is written in a different format to the input file, as a `JSON`_ file. This is a more machine-readable version of the configuration file. You don't need to be familiar with the details of the format, but for information the ``fixed.conf`` configuration file is shown below. Here, the ``Pba_SCRI1043`` is unmodified, the ``Pbe_NCPPB_2795`` is stitched (hence ``concat`` appears in the filename), and the ``Pwa_CFBP_3304`` is both stitched and has ambiguity symbols replaced (so has ``concat_noambig`` in the filename).

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




.. _Anaconda: https://www.anaconda.com/what-is-anaconda/
.. _JSON: https://www.json.org/
.. _PyPI: https://pypi.python.org/pypi
