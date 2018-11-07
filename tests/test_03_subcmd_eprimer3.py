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

from nose.tools import raises

from diagnostic_primers.scripts import subcommands

from tools import assert_dirfiles_equal


class TestEPrimer3Subcommand(unittest.TestCase):
    """Class defining tests of the pdp.py eprimer3 subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join("tests", "test_input", "config")
        self.confoutdir = os.path.join("tests", "test_output", "config")
        self.conftargets = os.path.join("tests", "test_targets", "config")
        self.outdir = os.path.join("tests", "test_output", "eprimer3", "scriptout")
        self.targetdir = os.path.join("tests", "test_targets", "eprimer3")
        self.ep3_exe = "eprimer3"
        self.scheduler = "multiprocessing"
        self.workers = None
        self.hybridprobe = False
        # Default values for ePrimer3 run - not modified by any tests,
        # defined in parsers.py
        self.ep3_defaults = {
            "ep_numreturn": 10,
            "ep_osize": 20,
            "ep_minsize": 18,
            "ep_maxsize": 22,
            "ep_opttm": 59,
            "ep_mintm": 58,
            "ep_maxtm": 60,
            "ep_ogcpercent": 55,
            "ep_mingc": 30,
            "ep_maxgc": 80,
            "ep_psizeopt": 100,
            "ep_psizemin": 50,
            "ep_psizemax": 150,
            "ep_maxpolyx": 3,
            "ep_osizeopt": 20,
            "ep_ominsize": 13,
            "ep_omaxsize": 30,
            "ep_otmopt": 69,
            "ep_otmmin": 68,
            "ep_otmmax": 70,
            "ep_ogcopt": 55,
            "ep_ogcmin": 30,
            "ep_ogcmax": 80,
        }

        # null logger instance that does nothing
        self.logger = logging.getLogger("TestConfigSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Dictionary of command-line namespaces
        self.argsdict = {
            "run": Namespace(
                infilename=os.path.join(self.confdir, "testprodigalconf.json"),
                outfilename=os.path.join(self.confoutdir, "ep3conf.json"),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=2,
                verbose=False,
                ep_hybridprobe=self.hybridprobe,
                ep_filter=False,
                disable_tqdm=True,
                **self.ep3_defaults
            ),
            "force": Namespace(
                infilename=os.path.join(self.confdir, "testprodigalconf.json"),
                outfilename=os.path.join(self.confoutdir, "ep3conf.json"),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=2,
                verbose=False,
                ep_hybridprobe=self.hybridprobe,
                ep_filter=False,
                disable_tqdm=True,
                **self.ep3_defaults
            ),
            "notconf": Namespace(
                infilename=os.path.join(self.confdir, "testprodigalconf.nojson"),
                outfilename=os.path.join(self.confoutdir, "ep3conf.json"),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False,
                disable_tqdm=True,
                ep_filter=False,
            ),
            "notjson": Namespace(
                infilename=os.path.join(self.confdir, "testin.conf"),
                outfilename=os.path.join(self.confoutdir, "ep3conf.json"),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False,
                disable_tqdm=True,
                ep_filter=False,
            ),
            "noforce": Namespace(
                infilename=os.path.join(self.confdir, "testprodigalconf.json"),
                outfilename=os.path.join(self.confoutdir, "ep3conf.json"),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=False,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False,
                ep_hybridprobe=self.hybridprobe,
                ep_filter=False,
                disable_tqdm=True,
                **self.ep3_defaults
            ),
        }

    def test_eprimer3_01_run(self):
        """eprimer3 subcommand executes correct primer design."""
        subcommands.subcmd_eprimer3(self.argsdict["run"], self.logger)
        # Check file contents
        assert_dirfiles_equal(self.outdir, self.targetdir)

    def test_eprimer3_02_force(self):
        """eprimer3 subcommand executes correctly and overwrites output."""
        subcommands.subcmd_eprimer3(self.argsdict["force"], self.logger)
        # Check file contents
        assert_dirfiles_equal(self.outdir, self.targetdir)

    @raises(SystemExit)
    def test_invalid_conf_file(self):
        """Script exits if ePrimer3 config file has wrong suffix."""
        subcommands.subcmd_eprimer3(self.argsdict["notconf"], self.logger)

    @raises(ValueError)
    def test_tsv_conf_file(self):
        """Error raised if .conf file provided for ePrimer3."""
        subcommands.subcmd_eprimer3(self.argsdict["notjson"], self.logger)

    @raises(SystemExit)
    def test_outdir_not_forced(self):
        """Script exits if not forcing ePrimer3 output overwrite."""
        subcommands.subcmd_eprimer3(self.argsdict["noforce"], self.logger)
        # Run twice to ensure the error is thrown if tests are out of order
        subcommands.subcmd_eprimer3(self.argsdict["noforce"], self.logger)
