#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_classify.py

Test classify subcommand for pdp.py script

This test suite is intended to be run from the repository root using:

nosetests -v

Individual test classes can be run using, e.g.:

$ nosetests -v tests/test_subcommands.py:TestConfigSubcommand

Each command CMD available at the command-line as pdp.py <CMD> is
tested in its own class (subclassing unittest.TestCase), where the
setUp() method defines input/output files, a null logger (picked up
by nosetests), and a dictionary of command lines, keyed by test name
with values that represent the command-line options.

For each test, command-line options are defined in a Namespace,
and passed as the sole argument to the appropriate subcommand
function from subcommands.py.

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

from diagnostic_primers.scripts import subcommands

from tools import assert_dirfiles_equal


class TestClassifySubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py classify subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "config")
        self.outdir = os.path.join("tests", "test_output", "classify")
        self.targetdir = os.path.join("tests", "test_targets", "classify")

        # null logger
        self.logger = logging.getLogger("TestClassifySubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command-line Namespaces
        self.argsdict = {
            "run": Namespace(
                infilename=os.path.join(self.confdir, "testclassify.json"),
                outdir=self.outdir,
                cl_force=True,
                verbose=False,
                disable_tqdm=True,
            )
        }

    def test_classify_run(self):
        """Classify command runs normally."""
        subcommands.subcmd_classify(self.argsdict["run"], self.logger)

        # Check output:
        self.logger.info("Comparing output primer sequences to targets")
        assert_dirfiles_equal(self.outdir, self.targetdir)
