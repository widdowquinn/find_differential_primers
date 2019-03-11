#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_primersearch.py

Test primersearch subcommand for pdp script

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

(c) The James Hutton Institute 2017-2019
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

Copyright (c) 2017-201 The James Hutton Institute

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


# Defined as global so it can be seen by the TestPrimersearchSubcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_primersearch")


class TestPrimersearchSubcommand(PDPTestCase):
    """Class defining tests of the pdp.py primersearch subcommand."""

    @classmethod
    def setUpClass(TestPrimersearchSubcommand):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "pdp_primersearch")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdp_primersearch")
        self.ps_exe = "primersearch"
        self.mismatchpercent = 0.1  # This must be in range [0,1]
        self.scheduler = "multiprocessing"
        self.workers = None

        # Make sure output directories exist
        os.makedirs(self.outdir, exist_ok=True)

        # null logger
        self.logger = logging.getLogger("TestPrimersearchSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # base namespace
        self.base_namespace = Namespace(
            ps_exe=self.ps_exe,
            ps_force=True,
            mismatchpercent=self.mismatchpercent,
            scheduler=self.scheduler,
            workers=self.workers,
            verbose=True,
            disable_tqdm=True,
            recovery=False,
        )

    def test_primersearch_prodigal_run(self):
        """primersearch command runs normally on prodigal features.

        pdp primersearch -v -f --disable_tqdm \
            --outdir=tests/test_output/pdp_primersearch/prodigal \
            tests/test_input/pdp_primersearch/screened_prod.json \
            tests/test_output/pdp_primersearch/primersearch_prod.json
        """
        subcommands.subcmd_primersearch(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "screened_prod.json"),
                    "outfilename": os.path.join(self.outdir, "primersearch_prod.json"),
                    "ps_dir": os.path.join(self.outdir, "prodigal"),
                },
            ),
            self.logger,
        )

        # Check file contents: config
        self.assertJsonEqual(
            os.path.join(self.outdir, "primersearch_prod.json"),
            os.path.join(self.targetdir, "primersearch_prod.json"),
        )
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences.
        self.logger.info("Comparing output JSON files to targets")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigal"),
            os.path.join(self.targetdir, "prodigal"),
            filt=(".json"),
        )

    def test_primersearch_prodigaligr_run(self):
        """primersearch command runs normally on prodigal intergenic features.

        pdp primersearch -v -f --disable_tqdm \
            --outdir=tests/test_output/pdp_primersearch/prodigaligr \
            tests/test_input/pdp_primersearch/screened_prodigr.json \
            tests/test_output/pdp_primersearch/primersearch_prodigr.json """
        subcommands.subcmd_primersearch(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "screened_prodigr.json"),
                    "outfilename": os.path.join(
                        self.outdir, "primersearch_prodigr.json"
                    ),
                    "ps_dir": os.path.join(self.outdir, "prodigaligr"),
                },
            ),
            self.logger,
        )

        # Check file contents: config
        self.assertJsonEqual(
            os.path.join(self.outdir, "primersearch_prodigr.json"),
            os.path.join(self.targetdir, "primersearch_prodigr.json"),
        )
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences.
        self.logger.info("Comparing output JSON files to targets")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigaligr"),
            os.path.join(self.targetdir, "prodigaligr"),
            filt=(".json",),
        )
