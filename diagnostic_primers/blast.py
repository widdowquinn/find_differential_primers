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

(c) The James Hutton Institute 2016-2019
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2016-2019 The James Hutton Institute

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
import os

from Bio.Blast.Applications import NcbiblastnCommandline

from diagnostic_primers import load_primers, write_primers


def build_commands(collection, blastexe, blastdb, outdir=None, existingfiles=[]):
    """Builds and returns a list of BLASTN command lines for screening

    The returned commands run BLASTN using primer sequences for each
    PDPData object in the collection as a query. It will create a new
    output directory if outdir does not exist, and will place FASTA
    sequences for each primer in that directory.

    If no output directory is provided, output will be placed under the same
    directory as the input sequence files. Otherwise, the query primer
    sequences are placed in the specified output directory.
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
        fastafname = "_".join([stem, "primers.fasta"])
        g.write_primers(fastafname)

        cline = build_blastscreen_cmd(fastafname, blastexe, blastdb, outdir)
        if os.path.split(cline.out)[-1] not in existingfiles:
            clines.append(cline)
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
    return NcbiblastnCommandline(
        query=queryfile,
        cmd=blastexe,
        db=blastdb,
        out=stem + ".blasttab",
        task="blastn-short",
        max_target_seqs=1,
        outfmt=6,
        perc_identity=90,
        ungapped=True,
    )


def apply_screen(blastfile, primerjson, jsondir=None, maxaln=15):
    """Apply the results from a BLASTN screen to a primer JSON file.

    Loads the BLASTN .blasttab file, and the JSON file defining primers. Where
    one or more primers in a pair are found to make qualifying matches in
    the .blasttab file, that pair is removed from the set of primers.

    The string '_screened' is appended to the JSON filestem, and the
    primer set with any screening modifications is written to a new file.
    The new file path is returned, for the calling function to substitute
    in to the parent PDPCollection, if required.

    blastfile     - path to BLASTN output .blasttab file
    primerjson    - path to JSON file describing primers
    jsondir       - path to directory for screened JSON files
    maxaln        - the maximum allowed alignment length

    Primer pairs where one or more sequences has alignment length greater
    than maxaln are removed from the set loaded in the JSON file
    """
    # Parse BLASTN output and identify noncompliant primers
    excluded = set()
    with open(blastfile, "r") as bfh:
        reader = csv.reader(bfh, delimiter="\t")
        for row in reader:
            if int(row[3]) > maxaln:
                excluded.add(row[0][:-4])

    # Parse primer JSON and remove primer pairs found in excluded
    primerdata = load_primers(primerjson, "json")
    primerdata = [prm for prm in primerdata if prm.name not in excluded]

    # Generate new JSON filename and write primers
    oldpath = os.path.split(primerjson)[:-1]
    newstem = os.path.splitext(os.path.split(primerjson)[-1])[0] + "_screened"
    # If no output directory is specified, place the new .json file alongside the old
    if jsondir is not None:
        newpath = os.path.join(jsondir, newstem)
    else:
        newpath = os.path.join(*oldpath, newstem)
    write_primers(primerdata, newpath + ".json", "json")
    write_primers(primerdata, newpath + ".bed", "bed")
    write_primers(primerdata, newpath + ".fasta", "fasta")

    # Return new JSON filename
    return newpath + ".json"


def parse_blasttab(fhandle):
    """Return the passed BLAST tab output file as a list of lists.

    This is used when testing for conserved BLAST output, as the
    exact format of the BLAST result can depend on the software version.
    For instance, the locally-installed version may be BLASTN+ 2.6.0,
    which reports match identity to 3sf, and the version on TravisCI may be
    BLASTN+ 2.2.28, which reports to 2sf.

    Returning a list of lines, parsed into the appropriate data type,
    allows for direct comparison of line content independent of formatting.
    """
    retval = []
    for line in fhandle.readlines():
        splitline = line.split("\t")
        data = splitline[:2]  # First two columns are strings
        data += [float(_) for _ in splitline[2:]]  # The rest are numeric
    return retval
