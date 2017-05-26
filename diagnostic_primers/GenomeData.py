#!/usr/bin/env python3
#
# GenomeData.py
#
# Class to hold information about a single genome input to diagnostic design
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

import os
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SPACER = "NNNNNCATCCATTCATTAATTAATTAATGAATGAATGNNNNN"

AMBIGUITIES = re.compile('[BDHKMRSVWY]')


class GenomeData(object):

    """Container for information about an input sequence for diagnostic PCR
    primer prediction.

    LAZY INPUT CHECKS
    =================
    On instantiation, paths are checked for validity and set, but there is no
    check on file contents for whether they need to be stitched, or if
    ambiguity symbols need to be substituted.

    It is intended that the calling code first checks (if required) for
    whether stitching/ambiguity substitution is needed with, e.g.
    gd = GenomeData(...)
    if gd.needs_stitch():
        <do something>
    if gd.has_ambiguity():
        <do something>
    """

    def __init__(self, name, groups, seqfile, features, primers):
        """Instantiate a GenomeData object.

        name - identifier for the GenomeData object
        groups - list (or comma-separated string) of diagnostic groups to which
                 the object belongs
        seqfile - path to FASTA file describing the object sequence
        features - path to file describing features/excluded regions on the
                   sequence
        primers - path to file describing primers that amplify the input,
                  sequence, in ePrimer3 format
        """
        self.name = name
        self._groups = set()      # Default group: empty set
        self.groups = groups
        self._seqfile = None      # Default: None
        self.altseqfiles = []     # Holds old filenames used for the object
        self.seqfile = seqfile
        self.features = features
        self.primers = primers
        self.cmds = {}            # Command-lines used for this object

    def as_list(self):
        """Returns attributes of the object as a list."""
        return [self.name, ','.join(sorted(self.groups)),
                '-' if self.seqfile is None else self.seqfile,
                '-' if self.features is None else self.features,
                '-' if self.primers is None else self.primers]

    def stitch(self):
        """If there is more than one sequence in the file, stitch with a
        spacer sequence.

        The following actions are applied:
        - stitch sequences together into new single sequence
        - write this sequence to a new file
        - copy old filename to self.altseqfiles, replace self.seqfile with
          new filename, and delete self.seqnames
        - replace feature and primer files in this object with None, as they
          no longer relate to the input sequence
        """
        if self.needs_stitch:
            seqdata = list(SeqIO.parse(self.seqfile, 'fasta'))
            catseq = SPACER.join([str(s.seq) for s in seqdata])
            newseq = SeqRecord(Seq(catseq),
                               id='_'.join([self.name, 'concatenated']),
                               description="%s, concatenated with spacers" %
                               self.name)
            outfilename = ''.join([os.path.splitext(self.seqfile)[0],
                                   '_concat', '.fas'])
            SeqIO.write([newseq], outfilename, 'fasta')
            self.altseqfiles.append(self.seqfile)
            self.seqfile = outfilename
            if hasattr(self, '_seqnames'):
                delattr(self, '_seqnames')
            self.features = None
            self.primers = None

    def has_ambiguities(self):
        """Returns True if the sequence(s) have non-N ambiguity symbols."""
        for s in SeqIO.parse(self.seqfile, 'fasta'):
            if re.search(AMBIGUITIES, str(s.seq)):
                return True

    def replace_ambiguities(self):
        """Replace non-N ambiguity symbols in self.seqfile with N.

        The following actions are applied:
        - replace ambiguity symbols with Ns
        - write new sequence(s) to file
        - copy old filename to self.altseqfiles, replace self.seqfile with
          new filename; delete self.seqnames
        - replace feature and primer files with None, as they no longer relate
          to the input sequence
        """
        if self.has_ambiguities():
            seqdata = list(SeqIO.parse(self.seqfile, 'fasta'))
            for s in seqdata:
                s.seq = Seq(re.sub(AMBIGUITIES, 'N', str(s.seq)),
                            s.seq.alphabet)
                s.id = '_'.join([s.id, 'noambig'])
            outfilename = ''.join([os.path.splitext(self.seqfile)[0],
                                   '_noambig', '.fas'])
            SeqIO.write(seqdata, outfilename, 'fasta')
            self.altseqfiles.append(self.seqfile)
            self.seqfile = outfilename
            if hasattr(self, '_seqnames'):
                delattr(self, '_seqnames')
            self.features = None
            self.primers = None

    @property
    def name(self):
        """Identifier for the GenomeData object."""
        return self._name

    @name.setter
    def name(self, value):
        self._name = str(value)

    @property
    def groups(self):
        """Groups to which the GenomeData object belongs"""
        return self._groups

    @groups.setter
    def groups(self, value):
        """Expects a list, a set, or a comma-separated string."""
        if isinstance(value, list):
            self._groups = self._groups.union(set(value))
        elif isinstance(value, set):
            self._groups = self._groups.union(value)
        elif isinstance(value, str):
            self._groups = self._groups.union(value.split(','))
        else:
            raise TypeError("groups should be set, list or " +
                            "comma-separated str")

    @property
    def seqfile(self):
        """Path to input sequence file."""
        return self._seqfile

    @seqfile.setter
    def seqfile(self, value):
        if not os.path.isfile(value):
            raise ValueError("%s is not a valid file path" % value)
        # Keep a record of other sequence files previously used with this
        # object, if we're updating
        if self._seqfile is not None:
            self.altseqfiles.append(self._seqfile)
        self._seqfile = value

    @property
    def features(self):
        """Path to feature location file."""
        return self._features

    @features.setter
    def features(self, value):
        if value is not None:
            if not os.path.isfile(value):
                raise ValueError("%s is not a valid file path" % value)
        self._features = value

    @property
    def primers(self):
        """Path to primer file."""
        return self._primers

    @primers.setter
    def primers(self, value):
        if value is not None:
            if not os.path.isfile(value):
                raise ValueError("%s is not a valid file path" % value)
        self._primers = value

    @property
    def seqnames(self):
        """Lazily returns list of names of sequences in self.seqfile."""
        if not hasattr(self, "_seqnames"):
            self._seqnames = [s.id for s in SeqIO.parse(self.seqfile, 'fasta')]
        return self._seqnames

    @property
    def needs_stitch(self):
        """Returns True if more than one sequence in self.seqfile."""
        return len(self.seqnames) > 1
