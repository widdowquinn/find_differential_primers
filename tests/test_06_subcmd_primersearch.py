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

import json
import logging
import os
import unittest

from argparse import Namespace

from nose.tools import assert_equal

from diagnostic_primers.scripts import subcommands

from tools import assert_dirfiles_equal, ordered, modify_namespace


class TestPrimersearchSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py primersearch subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "pdp_primersearch")
        self.outdir = os.path.join("tests", "test_output", "pdp_primersearch")
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
        )

    def test_primersearch__prodigal_run(self):
        """primersearch command runs normally."""
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
        with open(os.path.join(self.outdir, "primersearch_prod.json")) as ofh:
            with open(os.path.join(self.targetdir, "primersearch_prod.json")) as tfh:
                assert_equal(ordered(json.load(ofh)), ordered(json.load(tfh)))
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences.
        self.logger.info("Comparing output JSON files to targets")
        assert_dirfiles_equal(
            os.path.join(self.outdir, "prodigal"),
            os.path.join(self.targetdir, "prodigal"),
            filter=(".json",),
        )

    def test_primersearch__prodigaligr_run(self):
        """primersearch command runs normally."""
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
        with open(os.path.join(self.outdir, "primersearch_prodigr.json")) as ofh:
            with open(os.path.join(self.targetdir, "primersearch_prodigr.json")) as tfh:
                assert_equal(ordered(json.load(ofh)), ordered(json.load(tfh)))
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences.
        self.logger.info("Comparing output JSON files to targets")
        assert_dirfiles_equal(
            os.path.join(self.outdir, "prodigaligr"),
            os.path.join(self.targetdir, "prodigaligr"),
            filter=(".json",),
        )
