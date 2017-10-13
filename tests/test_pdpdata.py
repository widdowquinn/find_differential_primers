#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_pdpdata.py

Test instantiation, methods and attributions of the PDPData class

This test suite is intended to be run from the repository root using:

nosetests -v

(c) The James Hutton Institute 2017
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

Copyright (c) 2017 The James Hutton Institute

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

import os
import unittest

from diagnostic_primers.config import PDPData

from Bio import SeqIO
from nose.tools import assert_equal, raises


class BadName(object):

    """Dummy object that throws an error if str() is used on it."""

    def __init__(self):
        pass

    def __str__(self):
        raise ValueError("Cannot be a string")


class TestPDPData(unittest.TestCase):

    """Class defining tests of the PDPData object."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'sequences')
        self.name = "test_name"
        self.badname = BadName()
        self.groups_list = ['group1', 'group2', 'group3']
        self.groups_str = 'group1,group2,group3'
        self.groups_set = set(self.groups_list)
        self.seqfile = os.path.join(self.datadir,
                                    'GCF_000011605.1.fasta')
        self.stitchfile = os.path.join(self.datadir,
                                       'test_stitch.fasta')
        self.stitchout = os.path.join(self.datadir,
                                      'genomedata_test_concat.fas')
        self.noambigout = os.path.join(self.datadir,
                                       'genomedata_test_noambig.fas')
        self.features = None
        self.primers = None

    def test_instantiation_grouplist(self):
        """PDPData object instantiates with list of groups."""
        PDPData(self.name, self.groups_list, self.seqfile,
                self.features, self.primers)

    def test_instantiation_groupstr(self):
        """PDPData object instantiates with string of groups."""
        PDPData(self.name, self.groups_str, self.seqfile,
                self.features, self.primers)

    def test_instantiation_groupset(self):
        """PDPData object instantiates with set of groups."""
        PDPData(self.name, self.groups_set, self.seqfile,
                self.features, self.primers)

    def test_stitch(self):
        """PDPData object stitches multi-sequence input.

        The self.stitchfile path points to an input file with multiple
        sequences. We test if this needs stitching (it should), and if so
        stitch it.
        """
        gd = PDPData(self.name, self.groups_str, self.stitchfile,
                     self.features, self.primers)
        gd.stitch()  # assumes stitch is needed
        # We have to take the input filename from the GenomeData object,
        # as this will have changed if stitched/ambiguities removed.
        with open(gd.seqfile, 'r') as ifh:
            with open(self.stitchout, 'r') as ofh:
                assert_equal(ifh.read(), ofh.read())

    def test_noambig(self):
        """PDPData object replaces ambiguities in input.

        The self.stitchfile path points to an input file with multiple
        sequences, and ambiguity symbols. We test for ambiguity symbols, and
        replace them.
        """
        gd = PDPData(self.name, self.groups_str, self.stitchfile,
                     self.features, self.primers)
        if gd.has_ambiguities:
            gd.seqnames  # Forces lazy population of seqnames for testing
            gd.replace_ambiguities()
        # We have to take the input filename from the GenomeData object,
        # as this will have changed if stitched/ambiguities removed.
        with open(gd.seqfile, 'r') as ifh:
            with open(self.noambigout, 'r') as ofh:
                assert_equal(ifh.read(), ofh.read())

    def test_attributes(self):
        """PDPData attributes are accessible.

        Ensure lazily-populated attributes get assigned
        """
        gd = PDPData(self.name, self.groups_str, self.seqfile,
                     self.features, self.primers)
        assert_equal(gd.seqnames,
                     [s.id for s in SeqIO.parse(self.seqfile, 'fasta')])
        assert_equal(gd.features, self.features)

    @raises(OSError)
    def test_invalid_sequence(self):
        """PDPData errors with invalid input sequence file."""
        PDPData(self.name, self.groups_str, "seqfile.notexist",
                self.features, self.primers)

    @raises(OSError)
    def test_invalid_features(self):
        """PDPData errors with invalid feature file."""
        PDPData(self.name, self.groups_str, self.seqfile,
                "features.notexist", self.primers)

    @raises(OSError)
    def test_invalid_primers(self):
        """PDPData errors with invalid primers file."""
        PDPData(self.name, self.groups_str, self.seqfile,
                self.features, "primers.notexist")

    @raises(TypeError)
    def test_invalid_groups(self):
        """PDPData errors with invalid group type."""
        PDPData(self.name, 12345, self.seqfile,
                self.features, self.primers)

    @raises(TypeError)
    def test_invalid_name(self):
        """PDPData errors because name is not/cannot be a string."""
        PDPData(self.badname, self.groups_str, self.seqfile,
                self.features, self.primers)
