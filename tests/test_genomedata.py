#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_genomedata.py

Test instantiation, methods and attributions of the GenomeData class

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

from diagnostic_primers.GenomeData import GenomeData

from nose.tools import assert_equal, raises


class TestGenomeData(unittest.TestCase):

    """Class defining tests of the GenomeData object."""

    def setUp(self):
        """Set 'globals' for tests."""
        self.name = "test_name"
        self.groups_list = ['group1', 'group2', 'group3']
        self.groups_str = 'group1,group2,group3'
        self.seqfile = os.path.join('tests', 'test_input', 'sequences',
                                    'GCF_000011605.1.fasta')
        self.stitchfile = os.path.join('tests', 'test_input', 'sequences',
                                       'GCF_000291725.1.fasta')
        self.stitchout = os.path.join('tests', 'test_input', 'sequences',
                                      'genomedata_stitch_noambig.fas')
        self.features = None
        self.primers = None

    def test_instantiation_grouplist(self):
        """GenomeData object instantiates with list of groups."""
        gd = GenomeData(self.name, self.groups_list, self.seqfile,
                        self.features, self.primers)

    def test_instantiation_groupstr(self):
        """GenomeData object instantiates with string of groups."""
        gd = GenomeData(self.name, self.groups_str, self.seqfile,
                        self.features, self.primers)

    def test_instantiation_stitch(self):
        """GenomeData object stitches multi-sequence input.

        The self.stitchfile path points to an input file with multiple
        sequences. This is automagically stitched and has ambiguities
        removed when defined in the GenomeData object
        """
        gd = GenomeData(self.name, self.groups_str, self.stitchfile,
                        self.features, self.primers)
        # We have to take the input filename from the GenomeData object,
        # as this will have changed if stitched/ambiguities removed.
        with open(gd.seqfile, 'r') as ifh:
            with open(self.stitchout, 'r') as ofh:
                assert_equal(ifh.read(), ofh.read())

    @raises(ValueError)
    def test_invalid_features(self):
        """GenomeData errors with invalid feature file."""
        gd = GenomeData(self.name, self.groups_str, self.stitchfile,
                        "features.notexist", self.primers)

    @raises(ValueError)
    def test_invalid_primers(self):
        """GenomeData errors with invalid primers file."""
        gd = GenomeData(self.name, self.groups_str, self.stitchfile,
                        self.features, "primers.notexist")
