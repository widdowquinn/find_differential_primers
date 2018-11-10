#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_parsing.py

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

from diagnostic_primers.config import PDPCollection


class TestParser(unittest.TestCase):

    """Class defining tests of config file parsing and writing."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join("tests", "test_input", "parsing")
        self.config = os.path.join(self.datadir, "testconf.json")
        self.configtest = os.path.join(self.datadir, "testout.conf")
        self.confignew = os.path.join(self.datadir, "new.conf")

    def test_parse_config_json(self):
        """Test basic JSON config file parsing.

        This test loads in a JSON format config file, and checks
        i)  the correct number of input sequences is found (8)
        ii) the eight sequence groups are correctly identified
        """
        gc = PDPCollection("test")
        gc.from_json(self.config)
        self.assertEqual(8, len(gc), msg="Expected 8 sequences, got a different number")
        self.assertEqual(
            [
                "Pectobacterium",
                "atrosepticum_NCBI",
                "betavasculorum_NCBI",
                "gv1",
                "gv2",
                "gv7",
                "wasabiae_NCBI",
            ],
            gc.groups,
            msg="Sequence groups do not match expected list",
        )
