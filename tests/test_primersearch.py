#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_primersearch.py

Test generation of primersearch command-lines, primersearch file parsing and writing.

This test suite is intended to be run from the repository root using:

nosetests -v

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

import os
import subprocess
import sys

from diagnostic_primers import primersearch, config

from tools import PDPTestCase


class TestCommands(PDPTestCase):

    """Class defining tests of primersearch command lines."""

    def setUp(self):
        """Set parameters for tests."""
        self.ps_exe = "primersearch"
        self.inconf = os.path.join("tests", "test_input", "primersearch", "testps.json")
        self.outconf = os.path.join(
            "tests", "test_output", "primersearch", "psconf.json"
        )
        self.targetconf = os.path.join(
            "tests", "test_targets", "primersearch", "primersearch.json"
        )
        self.outdir = os.path.join("tests", "test_output", "primersearch")
        self.mismatchpercent = 10

    def test_primersearch_exe(self):
        """primersearch executable exists and runs."""
        cmd = "{} --version".format(self.ps_exe)
        result = subprocess.run(
            cmd,
            shell=sys.platform != "win32",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        # EMBOSS writes information out to STDERR
        self.assertEqual(result.stderr[:6], b"EMBOSS")

    def test_primersearch_cmd(self):
        """primersearch command builds correctly."""
        cmd = primersearch.build_command(
            self.ps_exe,
            "testfile.ep3",
            "testfile.fas",
            "query_ps_subject.primersearch",
            self.mismatchpercent,
        )
        target = " ".join(
            [
                self.ps_exe,
                "-auto",
                "-outfile=query_ps_subject.primersearch",
                "-seqall=testfile.fas",
                "-infile=testfile.ep3",
                "-mismatchpercent=10",
            ]
        )
        self.assertEqual(str(cmd), target)

    def test_primersearch_cmds(self):
        """primersearch command creation completes with no errors."""
        pdpc = config.PDPCollection()
        pdpc.from_json(self.inconf)
        primersearch.build_commands(
            pdpc, self.ps_exe, self.outdir, self.mismatchpercent
        )
