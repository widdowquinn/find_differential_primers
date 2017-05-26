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

from nose.tools import assert_equal


class TestGenomeData(unittest.TestCase):

    """Class defining tests of the GenomeData object."""

    def setUp(self):
        """Set 'globals' for tests."""
        self.name = "test_name"
        self.groups_list = ['group1', 'group2', 'group3']
        self.groups_str = 'group1,group2,group3'
        self.seqfile = os.path.join('tests', 'test_input', 'sequences',
                                    'GCF_000011605.1.fasta')
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
