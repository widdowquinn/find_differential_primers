#!/usr/bin/env python3
#
# primersearch.py
#
# Code to conduct primer crosshybridisation tests with primersearch
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016 The James Hutton Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import json
import os

from Bio.Emboss.Applications import PrimerSearchCommandline

from .eprimer3 import load_primers, write_primers


def build_commands(collection, primersearch_exe, primersearch_dir,
                   mismatchpercent):
    """Build and return a list of command-lines to run primersearch.

    collection          - PDPCollection describing analysis inputs
    primersearch_exe    - path to primersearch executable
    primersearch_dir    - path to primersearch output
    mismatchpercent     - allowed 'wobble' for primers
    """
    clines = []    # holds command lines

    # Make sure output directory exists
    os.makedirs(primersearch_dir, exist_ok=True)

    # Build command-lines for each input primer set
    # We generate a collection of target sequences from each PDPData
    # object, then loop over all PDPData objects, and build a command
    # line comparing that object's predicted primers to all other
    # target sequences.
    # Primersearch - bafflingly - doesn't accept EMBOSS' ePrimer3
    # format, so we need to write primers out as 3-column TSV
    targets = {dat.name: dat.seqfile for dat in collection.data}
    for dat in collection.data:
        # Get the primers from the JSON file and write them
        # to the output directory
        primerpath = os.path.join(primersearch_dir,
                                  '{}_primers.tab'.format(dat.name))
        primers = load_primers(dat.primers, 'json')
        write_primers(primers, primerpath, 'tsv')
        # Create a dictionary to hold target names, to be written
        # to a JSON file, and the path added to the PDPData object
        psdict = {'query': dat.name,
                  'primers': primerpath}
        for tgtname, tgtpath in targets.items():
            if dat.name != tgtname:
                # Name for output file is built from the PDPData
                # query/target object names
                outstem = os.path.join(primersearch_dir,
                                       '{}_ps_{}.primersearch'.format(dat.name,
                                                                      tgtname))
                # Add the output file to the PDPData primersearch attr
                psdict[tgtname] = outstem
                # Generate the primersearch cmd-line
                cline = build_command(primersearch_exe, primerpath,
                                      tgtpath, outstem,
                                      mismatchpercent)
                clines.append(cline)
        # Write primersearch output JSON file and add to PDPData object
        psjson = os.path.join(primersearch_dir,
                              '{}_primersearch.json'.format(dat.name))
        with open(psjson, 'w') as ofh:
            json.dump(psdict, ofh, sort_keys=True)
        dat.primersearch = psjson
    return clines


def build_command(primersearch_exe, primerfile, seqfile, filestem,
                  mismatchpercent):
    """Return a single primersearch command line.

    The returned command line queries a single primer set against a
    single genome sequence/input sequence.

    primersearch_exe      - path to primersearch
    primerfile            - path to input primer file (ePrimer3 format)
    seqfile               - path to target sequence
    filestem              - path to output file
    argdict               - optional arguments to primersearch
    mismatchpercent       - allowed 'wobble' for primers
    """
    # Build command-line with defaults
    cline = PrimerSearchCommandline()
    cline.auto = True
    cline.seqall = seqfile
    cline.infile = primerfile
    cline.outfile = filestem
    cline.mismatchpercent = mismatchpercent
    return cline
