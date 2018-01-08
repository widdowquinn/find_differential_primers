#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_tools.py

Test functions in the tools.py helper module

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

import logging
import os
import unittest

from argparse import Namespace

from nose.tools import raises

from diagnostic_primers.scripts import tools


class TestTools(unittest.TestCase):

    """Class defining tests of the tools.py module."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input')
        self.fakejsonfile = os.path.join(self.datadir, "fakefile.json")
        self.fakescript = "./fakescript"

        # Null logger for nosetests
        self.logger = logging.getLogger('TestTools logger')
        self.logger.addHandler(logging.NullHandler())

    @raises(SystemExit)
    def test_loadjson_nologger(self):
        """loading nonexistent JSON file fails."""
        args = Namespace(infilename=self.fakejsonfile)
        tools.load_config_json(args, self.logger)

    @raises(SystemExit)
    def test_runparallel_mp_fail(self):
        """running nonexistent script fails."""
        args = Namespace(scheduler="multiprocessing", workers=1,
                         verbose=False)
        clines = [' '.join([self.fakescript, "-arg1 %d"]) % val for val in range(4)]
        tools.run_parallel_jobs(clines, args, self.logger)
