#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_prodigal.py

Test prodigal filtering for pdp.py script

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

import json
import logging
import os
import unittest

from argparse import Namespace

from nose.tools import assert_equal, raises

from diagnostic_primers.scripts import subcommands

from tools import assert_dirfiles_equal, ordered


class TestProdigalSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py prodigal subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join("tests", "test_input", "config")
        self.outconfdir = os.path.join("tests", "test_output", "config")
        self.outrundir = os.path.join("tests", "test_output", "prodigal")
        self.targetdir = os.path.join("tests", "test_targets", "prodigal")
        self.prodigal_exe = "prodigal"
        self.scheduler = "multiprocessing"
        self.workers = None

        # null logger instance that does nothing
        self.logger = logging.getLogger("TestConfigSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Dictionary of command-line namespaces
        self.argsdict = {
            "run": Namespace(
                infilename=os.path.join(self.datadir, "testreducedep3conf.json"),
                outfilename=os.path.join(self.outconfdir, "prodconf.json"),
                filt_prodigal=True,
                filt_prodigaligr=False,
                filt_outdir=self.outrundir,
                filt_prodigal_exe=self.prodigal_exe,
                filt_force=True,
                filt_suffix="prodigal",
                filt_spacerlen=150,
                filt_flanklen=150,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=True,
            ),
            "notconf": Namespace(
                infilename=os.path.join(self.datadir, "fixedconf.nojson"),
                outfilename=os.path.join(self.outconfdir, "prodconf.json"),
                filt_prodigal=True,
                filt_prodigaligr=False,
                filt_outdir=self.outrundir,
                filt_prodigal_exe=self.prodigal_exe,
                filt_force=True,
                filt_suffix="prodigal",
                filt_spacerlen=150,
                filt_flanklen=150,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=True,
            ),
            "notjson": Namespace(
                infilename=os.path.join(self.datadir, "testin.conf"),
                outfilename=os.path.join(self.outconfdir, "prodconf.json"),
                filt_prodigal=True,
                filt_prodigaligr=False,
                filt_outdir=self.outrundir,
                filt_prodigal_exe=self.prodigal_exe,
                filt_force=True,
                filt_suffix="prodigal",
                filt_spacerlen=150,
                filt_flanklen=150,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=True,
            ),
        }

    def test_prodigal_run(self):
        """prodigal subcommand produces correct annotation."""
        subcommands.subcmd_filter(self.argsdict["run"], self.logger)

        # Check file contents
        assert_dirfiles_equal(self.outrundir, self.targetdir)

    @raises(SystemExit)
    def test_invalid_conf_file(self):
        """Script exits if prodigal config file has wrong suffix."""
        subcommands.subcmd_filter(self.argsdict["notconf"], self.logger)

    @raises(ValueError)
    def test_tsv_conf_file(self):
        """Error raised if .conf file provided for prodigal."""
        subcommands.subcmd_filter(self.argsdict["notjson"], self.logger)
