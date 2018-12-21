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

(c) The James Hutton Institute 2017-2018
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

Copyright (c) 2017-2018 The James Hutton Institute

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

from argparse import Namespace

from diagnostic_primers.scripts import subcommands

from tools import PDPTestCase, modify_namespace


class TestClassifySubcommand(PDPTestCase):
    """Class defining tests of the pdp.py classify subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "pdp_classify")
        self.outdir = os.path.join("tests", "test_output", "pdp_classify")
        self.targetdir = os.path.join("tests", "test_targets", "pdp_classify")

        # null logger
        self.logger = logging.getLogger("TestClassifySubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command-line Namespaces
        self.base_namespace = Namespace(cl_force=True, verbose=True, disable_tqdm=True)

    def test_classify_prodigal_run(self):
        """Classify command runs normally.

        pdp classify -v -f --disable_tqdm \
            tests/test_input/pdp_classify/primersearch_prod.json \
            tests/test_output/pdp_classify/prodigal
        """
        subcommands.subcmd_classify(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "primersearch_prod.json"),
                    "outdir": os.path.join(self.outdir, "prodigal"),
                },
            ),
            self.logger,
        )

        # Check output:
        self.logger.info("Comparing output primer sequences to targets")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigal"),
            os.path.join(self.targetdir, "prodigal"),
        )

    def test_classify_prodigaligr_run(self):
        """Classify command runs normally.

        pdp classify -v -f --disable_tqdm \
            tests/test_input/pdp_classify/primersearch_prodigaligr.json \
            tests/test_output/pdp_classify/prodigaligr
        """
        subcommands.subcmd_classify(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(
                        self.confdir, "primersearch_prodigr.json"
                    ),
                    "outdir": os.path.join(self.outdir, "prodigaligr"),
                },
            ),
            self.logger,
        )

        # Check output:
        self.logger.info("Comparing output primer sequences to targets")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigaligr"),
            os.path.join(self.targetdir, "prodigaligr"),
        )
