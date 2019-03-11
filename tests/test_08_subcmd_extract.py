#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_extract.py

Test extract subcommand for pdp script

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

from diagnostic_primers.scripts import subcommands

from tools import PDPTestCase, modify_namespace

# Defined as global so it can be seen by the TestExtractSubcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_extract")


class TestExtractSubcommand(PDPTestCase):
    """Class defining tests of the pdp.py extract subcommand."""

    @classmethod
    def setUpClass(TestExtractSubcommand):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.indir = os.path.join("tests", "test_input", "pdp_extract")
        self.outdir = OUTDIR
        self.classifydir = os.path.join("tests", "test_output", "pdp_classify")
        self.targetdir = os.path.join("tests", "test_targets", "pdp_extract")
        self.filestem = "Pectobacterium_primers"
        self.mafft_exe = "mafft"
        self.scheduler = "multiprocessing"
        self.workers = 4

        # null logger
        self.logger = logging.getLogger("TestExtractSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # base namespace
        self.base_namespace = Namespace(
            outdir=self.outdir,
            verbose=True,
            ex_force=True,
            noalign=False,
            mafft_exe=self.mafft_exe,
            scheduler=self.scheduler,
            workers=self.workers,
            disable_tqdm=True,
            ex_minamplicon=50,
            ex_maxamplicon=300,
            recovery=False,
        )

    def test_extract_prodigal_run(self):
        """Extract command runs normally on prodigal regions (no alignment).

        pdp extract -v --disable_tqdm -f \
            --noalign \
            tests/test_input/pdp_extract/primersearch_prod.json \
            tests/test_output/pdp_classify/prodigal/Pectobacterium_primers.json \
            tests/test_output/pdp_extract/prodigal
        """
        subcommands.subcmd_extract(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "primersearch_prod.json"),
                    "primerfile": os.path.join(
                        self.classifydir, "prodigal", "%s.json" % self.filestem
                    ),
                    "outdir": os.path.join(self.outdir, "prodigal", "noalign"),
                    "noalign": True,
                },
            ),
            self.logger,
        )

        # Check output:
        self.logger.info("Comparing output amplicons to targets")
        # We have to infer the output location for the extracted amplicons.
        # This is defined by the filestem of the input JSON file
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigal", "noalign", self.filestem),
            os.path.join(self.targetdir, "prodigal", "noalign", self.filestem),
        )

    def test_extract_progidal_align(self):
        """Extract command runs normally on prodigal regions (with MAFFT alignment).

        pdp extract -v --disable_tqdm -f \
            tests/test_input/pdp_extract/primersearch_prod.json \
            tests/test_output/pdp_classify/prodigal/Pectobacterium_primers.json \
            tests/test_output/pdp_extract/prodigal
        """
        subcommands.subcmd_extract(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "primersearch_prod.json"),
                    "primerfile": os.path.join(
                        self.classifydir, "prodigal", "%s.json" % self.filestem
                    ),
                    "outdir": os.path.join(self.outdir, "prodigal", "align"),
                },
            ),
            self.logger,
        )
        # Check output:
        self.logger.info("Comparing output amplicons to targets")
        # We have to infer the output location for the extracted amplicons.
        # This is defined by the filestem of the input JSON file
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigal", "align", self.filestem),
            os.path.join(self.targetdir, "prodigal", "align", self.filestem),
        )

    def test_extract_prodigaligr_run(self):
        """Extract command runs normally on prodigal IGR regions (no alignment).

        pdp extract -v --disable_tqdm -f \
            --noalign \
            tests/test_input/pdp_extract/primersearch_prodigr.json \
            tests/test_output/pdp_classify/prodigaligr/Pectobacterium_primers.json \
            tests/test_output/pdp_extract/prodigaligr
        """
        subcommands.subcmd_extract(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "primersearch_prodigr.json"),
                    "primerfile": os.path.join(
                        self.classifydir, "prodigaligr", "%s.json" % self.filestem
                    ),
                    "outdir": os.path.join(self.outdir, "prodigaligr", "noalign"),
                    "noalign": True,
                },
            ),
            self.logger,
        )

        # Check output:
        self.logger.info("Comparing output amplicons to targets")
        # We have to infer the output location for the extracted amplicons.
        # This is defined by the filestem of the input JSON file
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigaligr", "noalign", self.filestem),
            os.path.join(self.targetdir, "prodigaligr", "noalign", self.filestem),
        )

    def test_extract_progidaligr_align(self):
        """Extract command runs normally on prodigal IGR regions (with MAFFT alignment).

        pdp extract -v --disable_tqdm -f \
            tests/test_input/pdp_extract/primersearch_prodigr.json \
            tests/test_output/pdp_classify/prodigaligr/Pectobacterium_primers.json \
            tests/test_output/pdp_extract/prodigaligr
        """
        subcommands.subcmd_extract(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "primersearch_prodigr.json"),
                    "primerfile": os.path.join(
                        self.classifydir, "prodigaligr", "%s.json" % self.filestem
                    ),
                    "outdir": os.path.join(self.outdir, "prodigaligr", "align"),
                },
            ),
            self.logger,
        )
        # Check output:
        self.logger.info("Comparing output amplicons to targets")
        # We have to infer the output location for the extracted amplicons.
        # This is defined by the filestem of the input JSON file
        self.assertDirsEqual(
            os.path.join(self.outdir, "prodigaligr", "align", self.filestem),
            os.path.join(self.targetdir, "prodigaligr", "align", self.filestem),
        )
