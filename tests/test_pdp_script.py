#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcommands.py

Test mains script call for pdp.py script

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

import logging
import os
import unittest

from diagnostic_primers.scripts import pdp_script

from nose.tools import raises


class TestPDPScript(unittest.TestCase):

    """Class defining tests of the main pdp.py script call."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join('tests', 'walkthrough')
        self.argsdict = {'validate': ['config', '--validate',
                                      os.path.join(
                                          self.confdir, 'pectoconf.tab'),
                                      '-v'],
                         'validate_log': ['config', '--validate',
                                          os.path.join(
                                              self.confdir, 'pectoconf.tab'),
                                          '-l',
                                          os.path.join(self.confdir,
                                                       'test.log')],
                         'failed_log': ['config', '--validate',
                                        os.path.join(
                                            self.confdir, 'pectoconf.tab'),
                                        '-l', os.path.join(self.confdir,
                                                           'notexist',
                                                           'test.log')]}

        # Null logger for nosetests
        self.logger = logging.getLogger('TestPDPScript logger')
        self.logger.addHandler(logging.NullHandler())

    def test_validate_config(self):
        """config subcmd validates known good JSON config file."""
        pdp_script.run_pdp_main(self.argsdict['validate'])

    def test_validate_and_log(self):
        """config subcmd used to test logger."""
        pdp_script.run_pdp_main(self.argsdict['validate_log'])

    @raises(SystemExit)
    def test_logfail(self):
        """should fail because log not written."""
        pdp_script.run_pdp_main(self.argsdict['failed_log'])
