#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_process.py

Tests for the process module

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

from diagnostic_primers import process
from diagnostic_primers.GenomeCollection import ConfigSyntaxError

from nose.tools import assert_equal, raises


class TestProcess(unittest.TestCase):

    """Class defining tests of the GenomeCollection object."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'config')
        self.name = "test_collection"
        self.configfile = os.path.join(self.datadir, 'testin.conf')
        self.failconfig = os.path.join(self.datadir, 'broken.conf')

    def test_instantiate(self):
        """Process loads GenomeCollection from config file."""
        gc = process.load_collection(self.configfile, self.name)

    @raises(ConfigSyntaxError)
    def test_load_fail(self):
        """Process throws error with broken config."""
        gc = process.load_collection(self.failconfig, self.name)
