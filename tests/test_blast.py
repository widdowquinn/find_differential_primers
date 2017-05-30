#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_blast.py

Test generation of BLAST command-lines, and availability of dependencies.

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

import os
import subprocess
import sys
import unittest

from diagnostic_primers import (config, blast)

from nose.tools import assert_equal


class TestCommands(unittest.TestCase):

    """Class defining tests of BLAST command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'blast')
        self.primerfile = os.path.join(self.datadir, 'primers.fasta')
        self.screendb = os.path.join(self.datadir, 'primerscreen')
        self.blastexe = 'blastn'
        self.outdir = os.path.join('tests', 'test_output', 'blast')
        self.config = os.path.join('tests', 'test_input', 'config',
                                   'testprimer3conf.json')

    def test_blastexe(self):
        """BLASTN executable exists."""
        cmd = "{0} -version".format(self.blastexe)
        result = subprocess.run(cmd, shell=sys.platform != "win32",
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=True)
        assert_equal(result.stdout[:6], b'blastn')

    def test_blastscreen_cmd(self):
        """BLASTN primer screening command builds correctly."""
        cmd = blast.build_blastscreen_cmd(self.primerfile,
                                          self.blastexe,
                                          self.screendb,
                                          self.outdir)
        testcmd = ' '.join(['blastn',
                            '-out', os.path.join(self.outdir, 'primers.tab'),
                            '-outfmt', '6',
                            '-query', self.primerfile,
                            '-db', self.screendb,
                            '-max_target_seqs', '1',
                            '-task', 'blastn-short',
                            '-perc_identity', '90',
                            '-ungapped'])
        # We must test the string representation of the NcbiblastnCommandline
        assert_equal(str(cmd), testcmd)

    def test_blastscreen_cmds(self):
        """BLASTN primer screening command lines build correctly."""
        pdpc = config.PDPCollection()
        pdpc.from_json(self.config)
        clines = blast.build_commands(pdpc, self.blastexe, self.screendb,
                                      self.outdir)
        print(clines)
