#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_pdpcollection.py

Test instantiation, methods and attributions of the PDPCollection class

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

import json
import os
import unittest

from diagnostic_primers.config import (PDPData, PDPCollection,
                                       ConfigSyntaxError)

from nose.tools import assert_equal, raises


class TestGenomeCollection(unittest.TestCase):

    """Class defining tests of the GenomeCollection object."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'config')
        self.name = "test_collection"
        self.tabconfigfile = os.path.join(self.datadir, 'testconf.tab')
        self.jsonconfigfile = os.path.join(self.datadir, 'testconf.json')
        self.failconfig = os.path.join(self.datadir, 'broken.conf')

    def test_instantiate(self):
        """PDPCollection instantiates."""
        gc = PDPCollection(self.name)
        
    def test_load_json(self):
        """PDPCollection loads from JSON config file."""
        gc = PDPCollection(self.name)
        gc.from_json(self.jsonconfigfile)

    def test_load_tab(self):
        """PDPCollection loads from TSV config file."""
        gc = PDPCollection(self.name)
        gc.from_tab(self.tabconfigfile)
        
    @raises(ConfigSyntaxError)
    def test_load_fail_1(self):
        """PDPCollection throws error with broken TSV config."""
        gc = PDPCollection(self.name)
        gc.from_tab(self.failconfig)

    @raises(json.decoder.JSONDecodeError)
    def test_load_fail_2(self):
        """PDPCollection throws error when loading TSV config as if JSON."""
        gc = PDPCollection(self.name)
        gc.from_json(self.failconfig)

    def test_collection_data(self):
        """PDPCollection contains PDPData."""
        gc = PDPCollection(self.name)
        gc.from_json(self.jsonconfigfile)
        for item in gc.data:
            assert type(item) == PDPData
