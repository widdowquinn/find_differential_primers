#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_eprimer3.py

Test eprimer3 subcommand for pdp.py script

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

from tools import PDPTestCase, get_primer3_version, modify_namespace


# Defined as global so it can be seen by the TestEPrimer3Subcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_eprimer3")

# Available primer3 version as global so that pytest.skipif() can see it
PRIMER3_VERSION = get_primer3_version()


class TestEPrimer3Subcommand(PDPTestCase):
    """Class defining tests of the pdp.py eprimer3 subcommand."""

    @classmethod
    def setUpClass(TestEPrimer3Subcommand):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "pdp_eprimer3")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdp_eprimer3")
        self.ep3_exe = "eprimer3"
        self.scheduler = "multiprocessing"
        self.workers = None

        # null logger instance that does nothing
        self.logger = logging.getLogger("TestConfigSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # base Namespace
        self.base_namespace = Namespace(
            eprimer3_dir=self.outdir,
            eprimer3_exe=self.ep3_exe,
            eprimer3_force=True,
            scheduler=self.scheduler,
            workers=4,
            verbose=True,
            ep_hybridprobe=False,
            ep_filter=False,
            disable_tqdm=True,
            ep_numreturn=10,
            ep_osize=20,
            ep_minsize=18,
            ep_maxsize=22,
            ep_opttm=59,
            ep_mintm=58,
            ep_maxtm=60,
            ep_ogcpercent=55,
            ep_mingc=30,
            ep_maxgc=80,
            ep_psizeopt=100,
            ep_psizemin=50,
            ep_psizemax=150,
            ep_maxpolyx=3,
            ep_osizeopt=20,
            ep_ominsize=13,
            ep_omaxsize=30,
            ep_otmopt=69,
            ep_otmmin=68,
            ep_otmmax=70,
            ep_ogcopt=55,
            ep_ogcmin=30,
            ep_ogcmax=80,
            recovery=False,
        )

    @pytest.mark.skipif(PRIMER3_VERSION[0] != 1, reason="requires primer3 v1")
    def test_eprimer3_01_run(self):
        """eprimer3 subcommand recapitulates primer design for small input set.

        pdp eprimer3 -v \
            --outdir=tests/test_output/pdp_eprimer3/subset \
            tests/test_input/pdp_eprimer3/subsetconf.json \
            tests/test_output/pdp_eprimer3/subsetep3conf.json
        """
        subcommands.subcmd_eprimer3(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "subsetconf.json"),
                    "outfilename": os.path.join(self.outdir, "subsetep3conf.json"),
                    "eprimer3_dir": os.path.join(self.outdir, "subset"),
                },
            ),
            self.logger,
        )
        # Check file contents
        self.assertDirsEqual(
            os.path.join(self.outdir, "subset"), os.path.join(self.targetdir, "subset")
        )

    @pytest.mark.skipif(PRIMER3_VERSION[0] != 1, reason="requires primer3 v1")
    def test_eprimer3_02_force(self):
        """eprimer3 subcommand executes and overwrites existing output.

        This is the same test as test_eprimer3_01_run:

        pdp eprimer3 -v -f \
            --outdir=tests/test_output/pdp_eprimer3/subset \
            tests/test_input/pdp_eprimer3/subsetconf.json \
            tests/test_output/pdp_eprimer3/subsetep3conf.json \
        """
        self.test_eprimer3_01_run()

    @pytest.mark.skipif(PRIMER3_VERSION[0] != 1, reason="requires primer3 v1")
    def test_eprimer3_03_noforce(self):
        """Script exits when not forcing ePrimer3 output overwrite of existing output.

        pdp eprimer3 -v \
            --outdir=tests/test_output/pdp_eprimer3/subset \
            tests/test_input/pdp_eprimer3/subsetconf.json \
            tests/test_output/pdp_eprimer3/subsetep3conf.json
        """
        with pytest.raises(SystemExit):
            subcommands.subcmd_eprimer3(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.confdir, "subsetconf.json"),
                        "outfilename": os.path.join(self.outdir, "subsetep3conf.json"),
                        "eprimer3_dir": os.path.join(self.outdir, "subset"),
                        "eprimer3_force": False,
                    },
                ),
                self.logger,
            )

    @pytest.mark.skipif(PRIMER3_VERSION[0] != 1, reason="requires primer3 v1")
    def test_invalid_conf_file(self):
        """Script exits when ePrimer3 config file has wrong suffix.

        pdp eprimer3 -v \
            --outdir=tests/test_output/pdp_eprimer3/subset \
            tests/test_input/pdp_eprimer3/testprodigalconf.nojson \
            tests/test_output/pdp_eprimer3/ep3conf.json
        """
        with pytest.raises(SystemExit):
            subcommands.subcmd_eprimer3(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(
                            self.confdir, "testprodigalconf.nojson"
                        ),
                        "outfilename": os.path.join(self.outdir, "ep3conf.json"),
                    },
                ),
                self.logger,
            )

    @pytest.mark.skipif(PRIMER3_VERSION[0] != 1, reason="requires primer3 v1")
    def test_tsv_conf_file(self):
        """Error raised when .conf file provided for ePrimer3.

        pdp eprimer3 -v \
            --outdir=tests/test_output/pdp_eprimer3/subset \
            tests/test_input/pdp_eprimer3/testin.conf \
            tests/test_output/pdp_eprimer3/ep3conf.json
        """
        with pytest.raises(ValueError):
            subcommands.subcmd_eprimer3(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.confdir, "testin.conf"),
                        "outfilename": os.path.join(self.outdir, "ep3conf.json"),
                    },
                ),
                self.logger,
            )
