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


class PrimerSearchRecord(object):

    """Container for single PrimerSearch record

    This will contain the entire data from a PrimerSearch record
    """

    def __init__(self, name):
        self._name = str(name)
        self._amplimers = []

    def add_amplimer(self, amplimer):
        """Add a PrimerSearchAmplimer to the list of amplimers"""
        self._amplimers.append(amplimer)

    @property
    def name(self):
        return self._name

    @property
    def amplimers(self):
        return self._amplimers[:]

    def __str__(self):
        outstr = [f"\nPrimer name {self.name}"]
        for idx, amp in enumerate(self.amplimers):
            outstr += [f"Amplimer {idx + 1}",
                       f"\tSequence: {amp.sequence}",
                       f"\tAmplimer length: {len(amp)} bp"]
        return '\n'.join(outstr) + '\n'


class PrimerSearchAmplimer(object):

    """Container for single PrimerSearch amplimer

    This will contain the entire data from a PrimerSearch record amplimer
    """

    def __init__(self, name):
        self._name = str(name)
        self._seqname = ''
        self._len = None

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._seqname

    @sequence.setter
    def sequence(self, val):
        """Name of originating sequence"""
        self._seqname = str(val)

    @property
    def length(self):
        return self._len

    @length.setter
    def length(self, val):
        self._len = int(val)

    def __len__(self):
        return self.length


def parse_output(filename):
    """Return the contents of a PrimerSearch output file

    TODO: This is a hacky wee function - Scanner/Consumer would be better
          While we're developing, this will do.
    """
    records = []
    with open(filename, 'r') as ifh:
        for line in ifh:
            if line.startswith("Primer name"):   # Start of record
                rname = line.split("Primer name")[-1].strip()
                record = PrimerSearchRecord(rname)
            if line.startswith("Amplimer"):
                aname = line.strip()
                amplimer = PrimerSearchAmplimer(aname)
                record.add_amplimer(amplimer)
            if line.strip().startswith("Sequence"):
                seqname = line.split("Sequence:")[-1].strip()
                amplimer.sequence = seqname
            if line.strip().startswith("Amplimer length"):
                alen = int(line.split("Amplimer length:")[-1].strip().split()[0])
                amplimer.length = alen
                records.append(record)
    return records
