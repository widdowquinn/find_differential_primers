#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_03b_subcmd_primer3.py

Test primer3 subcommand for pdp script. These tests require primer3 v2+.

This test suite is intended to be run from the repository root using:

pytest -v

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


# Defined as global so it can be seen by the TestPrimer3Subcommand() class
# and setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_primer3")

# Available primer3 version as global so that pytest.skipif() can see it
PRIMER3_VERSION = get_primer3_version()


class TestPrimer3Subcommand(PDPTestCase):
    """Class defining tests of the pdp primer3 subcommand."""

    @classmethod
    def setUpClass(TestPrimer3Subcommand):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "pdp_primer3")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdp_primer3")
        self.p3_exe = "primer3_core"
        self.scheduler = "multiprocessing"
        self.workers = None

        # null logger instance that does nothing
        self.logger = logging.getLogger("TestPrimer3Subcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # path to thermodynamic parameters (needed for Travis/testing)
        self.therm_param_path = os.path.join(
            "tests", "test_input", "primer3", "primer3_config"
        )

        # base Namespace
        self.base_namespace = Namespace(
            primer3_dir=self.outdir,
            primer3_exe=self.p3_exe,
            primer3_force=True,
            scheduler=self.scheduler,
            workers=4,
            verbose=True,
            p3_hybridprobe=False,
            p3_filter=False,
            disable_tqdm=True,
            p3_param_path=self.therm_param_path,
            p3_numreturn=10,
            p3_osize=20,
            p3_minsize=18,
            p3_maxsize=22,
            p3_wt_lt=2,
            p3_wt_gt=2,
            p3_opttm=59,
            p3_mintm=58,
            p3_maxtm=60,
            p3_ogcpercent=55,
            p3_mingc=30,
            p3_maxgc=80,
            p3_psizeopt=100,
            p3_psizemin=50,
            p3_psizemax=150,
            p3_maxpolyx=3,
            p3_osizeopt=20,
            p3_ominsize=13,
            p3_omaxsize=30,
            p3_otmopt=69,
            p3_otmmin=68,
            p3_otmmax=70,
            p3_ogcopt=55,
            p3_ogcmin=30,
            p3_ogcmax=80,
            recovery=False,
        )

    @pytest.mark.skipif(PRIMER3_VERSION[0] < 2, reason="requires primer3 v2+")
    def test_primer3_01_run(self):
        """primer3 subcommand recapitulates primer design for small input set.

        pdp primer3 -v \
            --outdir=tests/test_output/pdp_primer3/subset \
            tests/test_input/pdp_primer3/subsetconf.json \
            tests/test_output/pdp_primer3/subsetep3conf.json
        """
        subcommands.subcmd_primer3(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "subsetconf.json"),
                    "outfilename": os.path.join(self.outdir, "subsetep3conf.json"),
                    "primer3_dir": os.path.join(self.outdir, "subset"),
                },
            ),
            self.logger,
        )
        # Check file contents
        self.assertDirsEqual(
            os.path.join(self.outdir, "subset"), os.path.join(self.targetdir, "subset")
        )

    @pytest.mark.skipif(PRIMER3_VERSION[0] < 2, reason="requires primer3 v2+")
    def test_primer3_02_force(self):
        """primer3 subcommand executes and overwrites existing output.

        This is the same test as test_primer3_01_run:


        pdp primer3 -v -f \
            --outdir=tests/test_output/pdp_primer3/subset \
            tests/test_input/pdp_primer3/subsetconf.json \
            tests/test_output/pdp_primer3/subsetep3conf.json
        """
        self.test_primer3_01_run()

    @pytest.mark.skipif(PRIMER3_VERSION[0] < 2, reason="requires primer3 v2+")
    def test_primer3_03_noforce(self):
        """Script exits when not forcing primer3 output overwrite of existing output.

        pdp primer3 -v \
            --outdir=tests/test_output/pdp_primer3/subset \
            tests/test_input/pdp_primer3/subsetconf.json \
            tests/test_output/pdp_primer3/subsetep3conf.json
        """
        with pytest.raises(SystemExit):
            subcommands.subcmd_primer3(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.confdir, "subsetconf.json"),
                        "outfilename": os.path.join(self.outdir, "subsetep3conf.json"),
                        "primer3_dir": os.path.join(self.outdir, "subset"),
                        "primer3_force": False,
                    },
                ),
                self.logger,
            )

    @pytest.mark.skipif(PRIMER3_VERSION[0] < 2, reason="requires primer3 v2+")
    def test_invalid_conf_file(self):
        """Script exits when primer3 config file has wrong suffix.

        pdp primer3 -v \
            --outdir=tests/test_output/pdp_primer3/subset \
            tests/test_input/pdp_primer3/testprodigalconf.nojson \
            tests/test_output/pdp_primer3/ep3conf.json
        """
        with pytest.raises(SystemExit):
            subcommands.subcmd_primer3(
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

    @pytest.mark.skipif(PRIMER3_VERSION[0] < 2, reason="requires primer3 v2+")
    def test_tsv_conf_file(self):
        """Error raised when .conf file provided for primer3.

        pdp primer3 -v \
            --outdir=tests/test_output/pdp_primer3/subset \
            tests/test_input/pdp_primer3/testin.conf \
            tests/test_output/pdp_primer3/ep3conf.json
        """
        with pytest.raises(ValueError):
            subcommands.subcmd_primer3(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.confdir, "testin.conf"),
                        "outfilename": os.path.join(self.outdir, "ep3conf.json"),
                    },
                ),
                self.logger,
            )
