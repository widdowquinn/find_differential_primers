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

class GenomeData(object):
    """Container for information about an input sequence for diagnostic PCR
    primer prediction.
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

    def as_list(self):
        """Returns attributes of the object as a list."""
        return [self.name, ','.join(sorted(self.groups)), 
                '-' if self.seqfile is None else self.seqfile,
                '-' if self.features is None else self.features,
                '-' if self.primers is None else self.primers]

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
            self._groups = self_groups.union(set(value))
        elif isinstance(value, set):
            self._groups = self._groups.union(value)
        elif isinstance(value, str):
            self._groups = self._groups.union(value.split(','))
        else:
            raise TypeError("groups should be set, list or comma-separated str")
    
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
