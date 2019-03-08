#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_blastscreen.py

Test blastscreen subcommand for pdp script

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

Copyright (c) 2017-2019 The James Hutton Institute

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

import pytest

from diagnostic_primers.scripts import subcommands

from tools import PDPTestCase, modify_namespace


# Defined as global so it can be seen by the TestBlastscreenSubcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_blastscreen")


class TestBlastscreenSubcommand(PDPTestCase):
    """Class defining tests of the pdp blastscreen subcommand."""

    @classmethod
    def setUpClass(TestBlastscreenSubcommand):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.dbdir = os.path.join("tests", "test_input", "pdp_blastscreen", "blastdb")
        self.confdir = os.path.join("tests", "test_input", "pdp_blastscreen")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdp_blastscreen")
        self.blast_exe = "blastn"
        self.maxaln = 15
        self.scheduler = "multiprocessing"
        self.workers = None

        # null logger
        self.logger = logging.getLogger("TestBlastscreenSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # base Namespace
        self.base_namespace = Namespace(
            bs_db=os.path.join(self.dbdir, "e_coli_screen.fna"),
            bs_exe=self.blast_exe,
            bs_force=True,
            maxaln=self.maxaln,
            scheduler=self.scheduler,
            workers=self.workers,
            verbose=True,
            disable_tqdm=True,
            recovery=False,
        )

    def test_blastscreen_prodigal_01_run(self):
        """blastscreen of prodigal data overwrites existing folder.

        pdp blastscreen -v -f --disable_tqdm \
            --db=tests/test_input/pdp_blastscreen/blastdb/e_coli_screen.fna \
            --outdir=tests/test_output/pdp_blastscreen/prodigal \
            tests/test_input/pdp_blastscreen/dedupe_prod.json \
            tests/test_output/pdp_blastscreen/screened_prod.json
        """
        subcommands.subcmd_blastscreen(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "dedupe_prod.json"),
                    "outfilename": os.path.join(self.outdir, "screened_prod.json"),
                    "bs_dir": os.path.join(self.outdir, "prodigal"),
                    "bs_jsondir": os.path.join(self.outdir, "prodigal"),
                },
            ),
            self.logger,
        )

        # Check file contents: config
        self.logger.info("Checking output config file against target file")
        self.assertJsonEqual(
            os.path.join(self.outdir, "screened_prod.json"),
            os.path.join(self.targetdir, "screened_prod.json"),
        )
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences:
        self.logger.info("Comparing output sequences/JSON to target")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigal"),
            os.path.join(self.targetdir, "prodigal"),
        )

    def test_blastscreen_prodigaligr_run(self):
        """blastscreen of prodigaligr data overwrites existing folder.

        pdp blastscreen -v -f --disable_tqdm \
            --db=tests/test_input/pdp_blastscreen/blastdb/e_coli_screen.fna \
            --outdir=tests/test_output/pdp_blastscreen/prodigaligr \
            tests/test_input/pdp_blastscreen/dedupe_prodigr.json \
            tests/test_output/pdp_blastscreen/screened_prodigr.json
        """
        subcommands.subcmd_blastscreen(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "dedupe_prodigr.json"),
                    "outfilename": os.path.join(self.outdir, "screened_prodigr.json"),
                    "bs_dir": os.path.join(self.outdir, "prodigaligr"),
                    "bs_jsondir": os.path.join(self.outdir, "prodigaligr"),
                },
            ),
            self.logger,
        )

        # Check file contents: config
        self.logger.info("Checking output config file against target file")
        self.assertJsonEqual(
            os.path.join(self.outdir, "screened_prodigr.json"),
            os.path.join(self.targetdir, "screened_prodigr.json"),
        )
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences:
        self.logger.info("Comparing output sequences/JSON to target")
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigaligr"),
            os.path.join(self.targetdir, "prodigaligr"),
        )

    def test_blastscreen_prodigal_02_noforce(self):
        """blastscreen command does not overwrite existing folder."""
        with pytest.raises(SystemExit):
            subcommands.subcmd_blastscreen(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.confdir, "dedupe_prod.json"),
                        "outfilename": os.path.join(self.outdir, "screened_prod.json"),
                        "bs_dir": os.path.join(self.outdir, "prodigal"),
                        "bs_jsondir": os.path.join(self.outdir, "prodigal"),
                        "bs_force": False,
                    },
                ),
                self.logger,
            )

    def test_blastscreen_nodb(self):
        """blastscreen command does not run without database."""
        with pytest.raises(SystemExit):
            subcommands.subcmd_blastscreen(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.confdir, "dedupe_prod.json"),
                        "outfilename": os.path.join(self.outdir, "screened_prod.json"),
                        "bs_dir": os.path.join(self.outdir, "prodigal"),
                        "bs_jsondir": os.path.join(self.outdir, "prodigal"),
                        "bs_db": None,
                    },
                ),
                self.logger,
            )
