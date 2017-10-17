#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""blast.py

Module that interacts with BLAST+ for diagnostic primer prediction.

The module provides functions that generate BLASTN command lines,
tuned for screening primer sequences. This module does NOT run the
commands - that is handled by one of the schedulers (multiprocessing
or SGE).

The module also has functions that apply the results of a BLASTN screen,
filtering predicted primers out of the primer definition JSON file if
they are too similar to the screening database sequences.

(c) The James Hutton Institute 2016-2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2016-2017 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import csv
import json
import os

from Bio.Blast.Applications import NcbiblastnCommandline

from . import eprimer3


def build_commands(collection, blastexe, blastdb, outdir=None):
    """Builds and returns a list of BLASTN command lines for screening

    The returned commands run BLASTN using primer sequences for each
    PDPData object in the collection as a query.

    If no output directory is provided, output will be placed under the same
    directory as the input sequence files.
    """
    clines = []

    # Create output directory if required
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    for g in collection.data:
        if outdir is None:
            stem = os.path.splitext(g.primers)[0]
        else:
            stempath = os.path.split(os.path.splitext(g.seqfile)[0])
            stem = os.path.join(outdir, stempath[-1])

        # Create a FASTA format version of the primer sequences
        fastafname = '_'.join([stem, 'primers.fasta'])
        g.write_primers(fastafname)

        clines.append(build_blastscreen_cmd(fastafname, blastexe, blastdb,
                                            outdir))
    return clines


def build_blastscreen_cmd(queryfile, blastexe, blastdb, outdir=None):
    """Build and return a BLASTN command-line.

    - queryfile         Path to the primer sequences query file
    - blastexe          Path to BLASTN executable
    - blastdb           Path to screening BLAST database (nt)
    - outdir            Path to output directory
    - force             If True, existing outdir contents will be overwritten

    Constructs a BLASTN command using the Biopython interface. This is
    intended for BLAST screening of primer FASTA sequence files from a
    GenomeData object against a database of sequences that are intended to be
    negative examples (i.e. there should not, ideally, be any matches).

    The BLASTN command uses the `blastn-short` task, as this is a search with
    short input sequence, and writes results out in tab-separated format
    (outfmt=6), with a 90% sequence identity hard cutoff, with ungapped
    matches. It reports only the top matching target sequence, as any match
    is taken to indicate a potential off-target hybridisation.

    If the output directory is not specified, the output is written to the
    same directory as the input, with the same filestem.
    """
    if outdir is None:
        stem = os.path.splitext(queryfile)[0]
    else:
        filestem = os.path.splitext(os.path.split(queryfile)[-1])[0]
        stem = os.path.join(outdir, filestem)
    return NcbiblastnCommandline(query=queryfile,
                                 cmd=blastexe,
                                 db=blastdb,
                                 out=stem + '.tab',
                                 task='blastn-short',
                                 max_target_seqs=1,
                                 outfmt=6,
                                 perc_identity=90,
                                 ungapped=True)


def apply_screen(blastfile, primerjson, maxaln=15):
    """Apply the results from a BLASTN screen to a primer JSON file.

    Loads the BLASTN .tab file, and the JSON file defining primers. Where
    one or more primers in a pair are found to make qualifying matches in
    the .tab file, that pair is removed from the set of primers.

    The string '_screened' is appended to the JSON filestem, and the
    primer set with any screening modifications is written to a new file.
    The new file path is returned, for the calling function to substitute
    in to the parent PDPCollection, if required.

    blastfile     - path to BLASTN output .tab file
    primerjson    - path to JSON file describing primers
    maxaln        - the maximum allowed alignment length

    Primer pairs where one or more sequences has alignment length greater
    than maxaln are removed from the set loaded in the JSON file
    """
    # Parse BLASTN output and identify noncompliant primers
    excluded = set()
    with open(blastfile, 'r') as bfh:
        reader = csv.reader(bfh, delimiter='\t')
        for row in reader:
            if int(row[3]) > maxaln:
                excluded.add(row[0][:-4])

    # Parse primer JSON and remove primer pairs found in excluded
    primerdata = eprimer3.load_primers(primerjson, 'json')
    primerdata = [prm for prm in primerdata if prm.name not in excluded]

    # Generate new JSON filename and write primers
    newstem = os.path.splitext(primerjson)[0] + '_screened'
    jsonpath = newstem + '.json'
    eprimer3.write_primers(primerdata, jsonpath, 'json')
    eprimer3.write_primers(primerdata, newstem + '.fasta', 'fasta')

    # Return new JSON filename
    return jsonpath
