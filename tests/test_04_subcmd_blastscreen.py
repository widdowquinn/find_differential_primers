#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_blastscreen.py

Test blastscreen subcommand for pdp.py script

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

from tools import (assert_dirfiles_equal, ordered)


class TestBlastscreenSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py blastscreen subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.dbdir = os.path.join('tests', 'test_input', 'blastdb')
        self.confdir = os.path.join('tests', 'test_input', 'config')
        self.outconfdir = os.path.join('tests', 'test_output', 'config')
        self.outdir = os.path.join('tests', 'test_output', 'blastscreen')
        self.targetdir = os.path.join('tests', 'test_targets', 'blastscreen')
        self.targetconfdir = os.path.join('tests', 'test_targets', 'config')
        self.blast_exe = 'blastn'
        self.maxaln = 15
        self.scheduler = 'multiprocessing'
        self.workers = None

        # null logger
        self.logger = logging.getLogger('TestBlastscreenSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Command-line namespaces
        self.argsdict = {
            'run':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprimer3conf.json'),
                outfilename=os.path.join(self.outconfdir, 'screened.json'),
                bs_db=os.path.join(self.dbdir, 'e_coli_screen.fna'),
                bs_exe=self.blast_exe,
                bs_force=True,
                bs_dir=self.outdir,
                maxaln=self.maxaln,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
            'noforce':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprimer3conf.json'),
                outfilename=os.path.join(self.outconfdir, 'screened.json'),
                bs_db=os.path.join(self.dbdir, 'e_coli_screen.fna'),
                bs_exe=self.blast_exe,
                bs_force=False,
                bs_dir=self.outdir,
                maxaln=self.maxaln,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
            'nodb':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprimer3conf.json'),
                outfilename=os.path.join(self.outconfdir, 'screened.json'),
                bs_db=None,
                bs_exe=self.blast_exe,
                bs_force=False,
                bs_dir=self.outdir,
                maxaln=self.maxaln,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
        }

    def test_blastscreen_run(self):
        """blastscreen command runs normally and overwrites existing folder."""
        subcommands.subcmd_blastscreen(self.argsdict['run'], self.logger)

        # Check file contents: config
        self.logger.info("Checking output config file against target file")
        with open(os.path.join(self.outconfdir, 'screened.json')) as ofh:
            with open(os.path.join(self.targetconfdir,
                                   'screened.json')) as tfh:
                assert_equal(ordered(json.load(ofh)), ordered(json.load(tfh)))
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences:
        self.logger.info("Comparing output sequences/JSON to target")
        assert_dirfiles_equal(self.outdir, self.targetdir)

    @raises(SystemExit)
    def test_blastscreen_noforce(self):
        """blastscreen command does not overwrite existing folder."""
        subcommands.subcmd_blastscreen(self.argsdict['run'], self.logger)
        subcommands.subcmd_blastscreen(self.argsdict['noforce'], self.logger)

    @raises(SystemExit)
    def test_blastscreen_nodb(self):
        """blastscreen command does not run without database."""
        subcommands.subcmd_blastscreen(self.argsdict['nodb'], self.logger)
