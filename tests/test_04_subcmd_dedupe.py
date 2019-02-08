#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_03a_subcmd_dedupe.py

Test dedupe subcommand for pdp script

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

(c) The James Hutton Institute 2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2018 The James Hutton Institute

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
import shutil

from argparse import Namespace

from diagnostic_primers.scripts import subcommands

from tools import PDPTestCase, modify_namespace


# Defined as global so it can be seen by the TestDedupeSubcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_dedupe")


class TestDedupeSubcommand(PDPTestCase):
    """Class defining tests of the pdp.py dedupe subcommand."""

    @classmethod
    def setUpClass(TestDedupeSubcommand):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "pdp_dedupe")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdp_dedupe")

        # null logger
        self.logger = logging.getLogger("TestDedupeSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # base Namespace
        self.base_namespace = Namespace(verbose=True, disable_tqdm=True)

    def test_dedupe_prodigal_run(self):
        """dedupe command runs normally on prodigal regions.

        pdp dedupe -v --disable_tqdm \
            --dedupedir tests/test_output/pdp_dedupe/prodigal \
            tests/test_input/pdp_dedupe/prod_ep3conf.json \
            tests/test_output/pdp_dedupe/dedupe_prod.json
        """
        subcommands.subcmd_dedupe(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "prod_ep3conf.json"),
                    "outfilename": os.path.join(self.outdir, "dedupe_prod.json"),
                    "dd_dedupedir": os.path.join(self.outdir, "prodigal"),
                },
            ),
            self.logger,
        )

        # Check file contents: config
        self.logger.info("Checking output config file against target file")
        self.assertJsonEqual(
            os.path.join(self.outdir, "dedupe_prod.json"),
            os.path.join(self.targetdir, "dedupe_prod.json"),
        )
        self.logger.info("Config file checks against target correctly")

        # Check deduped sequences
        self.logger.info("Comparing output sequences/JSON to target")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigal"),
            os.path.join(self.targetdir, "prodigal"),
        )

    def test_dedupe_prodigaligr_run(self):
        """dedupe command runs normally on prodigal IGR regions.

        pdp dedupe -v --disable_tqdm \
            --dedupedir tests/test_output/pdp_dedupe/prodigaligr \
            tests/test_input/pdp_dedupe/prodigr_ep3conf.json \
            tests/test_output/pdp_dedupe/dedupe_prodigr.json
        """
        subcommands.subcmd_dedupe(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "prodigr_ep3conf.json"),
                    "outfilename": os.path.join(self.outdir, "dedupe_prodigr.json"),
                    "dd_dedupedir": os.path.join(self.outdir, "prodigaligr"),
                },
            ),
            self.logger,
        )

        # Check file contents: config
        self.logger.info("Checking output config file against target file")
        self.assertJsonEqual(
            os.path.join(self.outdir, "dedupe_prodigr.json"),
            os.path.join(self.targetdir, "dedupe_prodigr.json"),
        )
        self.logger.info("Config file checks against target correctly")

        # Check deduped sequences
        self.logger.info("Comparing output sequences/JSON to target")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigaligr"),
            os.path.join(self.targetdir, "prodigaligr"),
        )
