#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcommands.py

Test subcommand functions for pdp.py script

This test suite is intended to be run from the repository root using:

nosetests -v

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

from diagnostic_primers.scripts import subcommands

from argparse import Namespace
from nose.tools import assert_equal, raises


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


class TestConfigSubcommand(unittest.TestCase):

    """Class defining tests of the pdp.py config subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'config')
        self.outdir = os.path.join('tests', 'test_output', 'config')
        self.targetdir = os.path.join('tests', 'test_targets', 'config')
        self.tsv_to_json_fname = os.path.join(self.outdir,
                                              'tab_converted_conf.json')
        self.tsv_to_json_target = os.path.join(self.targetdir,
                                               'testconf.json')

        # null logger instance that does nothing
        self.logger = logging.getLogger('TestConfigSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Dictionary of command-line namespaces
        self.argsdict = {'validate_json_good':
                         Namespace(infilename=os.path.join(self.datadir,
                                                           'testconf.json'),
                                   verbose=True,
                                   validate=True,
                                   fix_sequences=False,
                                   to_json=False),
                         'validate_tsv_good':
                         Namespace(infilename=os.path.join(self.datadir,
                                                           'testconf.tab'),
                                   verbose=True,
                                   validate=True,
                                   fix_sequences=False,
                                   to_json=False),
                         'validate_config_bad':
                         Namespace(infilename=os.path.join(self.datadir,
                                                           'testin.conf'),
                                   verbose=True,
                                   validate=True,
                                   fix_sequences=False,
                                   to_json=False),
                         'tab_to_json':
                         Namespace(infilename=os.path.join(self.datadir,
                                                           'testconf.json'),
                                   verbose=True,
                                   validate=False,
                                   fix_sequences=False,
                                   to_json=self.tsv_to_json_fname),
                         }

    def test_validate_json_good(self):
        """config subcmd validates good JSON config file."""
        subcommands.subcmd_config(self.argsdict['validate_json_good'],
                                  self.logger)

    def test_validate_tab_good(self):
        """config subcmd validates good TSV config file."""
        subcommands.subcmd_config(self.argsdict['validate_tsv_good'],
                                  self.logger)

    def test_tsv_to_json(self):
        """config subcmd converts TSV config to JSON."""
        subcommands.subcmd_config(self.argsdict['tab_to_json'],
                                  self.logger)
        with open(self.tsv_to_json_fname, 'r') as fh1:
            with open(self.tsv_to_json_target, 'r') as fh2:
                assert_equal(ordered(json.load(fh1)),
                             ordered(json.load(fh2)))

    @raises(SystemExit)
    def test_validate_config_bad(self):
        """config subcmd validates badly-formatted config file."""
        subcommands.subcmd_config(self.argsdict['validate_config_bad'],
                                  self.logger)
