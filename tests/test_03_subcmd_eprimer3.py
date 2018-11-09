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
import unittest

from argparse import Namespace

from nose.tools import raises

from diagnostic_primers.scripts import subcommands

from tools import assert_dirfiles_equal, modify_namespace


class TestEPrimer3Subcommand(unittest.TestCase):
    """Class defining tests of the pdp.py eprimer3 subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "pdp_eprimer3")
        self.confoutdir = os.path.join("tests", "test_output", "pdp_config")
        self.conftargets = os.path.join("tests", "test_targets", "pdp_config")
        self.outdir = os.path.join("tests", "test_output", "pdp_eprimer3")
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
        )

    def test_eprimer3_01_run(self):
        """eprimer3 subcommand recapitulates primer design for small input set."""
        subcommands.subcmd_eprimer3(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "subsetconf.json"),
                    "outfilename": os.path.join(self.confoutdir, "subsetep3conf.json"),
                    "eprimer3_dir": os.path.join(self.outdir, "subset"),
                },
            ),
            self.logger,
        )
        # Check file contents
        assert_dirfiles_equal(
            os.path.join(self.outdir, "subset"), os.path.join(self.targetdir, "subset")
        )

    def test_eprimer3_02_force(self):
        """eprimer3 subcommand executes and overwrites existing output.

        This is the same test as test_eprimer3_01_run,
        """
        self.test_eprimer3_01_run()

    @raises(SystemExit)
    def test_invalid_conf_file(self):
        """Script exits if ePrimer3 config file has wrong suffix."""
        subcommands.subcmd_eprimer3(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "testprodigalconf.nojson"),
                    "outfilename": os.path.join(self.confoutdir, "ep3conf.json"),
                },
            ),
            self.logger,
        )

    @raises(ValueError)
    def test_tsv_conf_file(self):
        """Error raised if .conf file provided for ePrimer3."""
        subcommands.subcmd_eprimer3(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "testin.conf"),
                    "outfilename": os.path.join(self.confoutdir, "ep3conf.json"),
                },
            ),
            self.logger,
        )

    @raises(SystemExit)
    def test_eprimer3_03_noforce(self):
        """Script exits if not forcing ePrimer3 output overwrite."""
        subcommands.subcmd_eprimer3(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.confdir, "subsetconf.json"),
                    "outfilename": os.path.join(self.confoutdir, "subsetep3conf.json"),
                    "eprimer3_dir": os.path.join(self.outdir, "subset"),
                    "eprimer3_force": False,
                },
            ),
            self.logger,
        )
