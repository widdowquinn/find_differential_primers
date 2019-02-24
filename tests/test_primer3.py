#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_primer3.py

Test generation of Primer3 (v2+) command-lines, primer file parsing and writing.

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

import os
import subprocess
import shlex
import shutil

import pytest

from diagnostic_primers import config, primer3, load_primers, write_primers

from tools import PDPTestCase, get_primer3_version, modify_namespace


# Defined as global so it can be seen by the TestCommands() and TestParsing() classes
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "primer3")

# Available primer3 version as global so that pytest.skipif() can see it
PRIMER3_VERSION = get_primer3_version()


class TestCommands(PDPTestCase):
    """Class defining tests of ePrimer3 command-line generation."""

    @classmethod
    def setUpClass(TestCommands):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)
        os.makedirs(OUTDIR, exist_ok=True)

    def setUp(self):
        """Set parameters for tests."""
        self.primer3_exe = "primer3_core"
        self.datadir = os.path.join("tests", "test_input", "primer3")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "primer3")
        self.config = os.path.join(
            "tests", "test_input", "eprimer3", "testprodigalconf.json"
        )
        self.seqfile = os.path.join(self.datadir, "GCF_000740965.1_concat.fas")
        self.therm_param_path = os.path.join(self.datadir, "primer3_config")
        # Default values for ePrimer3 run - not modified by any tests,
        # defined in parsers.py
        self.p3_defaults = {
            "p3_param_path": self.therm_param_path,
            "p3_hybridprobe": False,
            "p3_numreturn": 10,
            "p3_osize": 20,
            "p3_minsize": 18,
            "p3_maxsize": 22,
            "p3_wt_lt": 1,
            "p3_wt_gt": 1,
            "p3_opttm": 59,
            "p3_mintm": 58,
            "p3_maxtm": 60,
            "p3_ogcpercent": 55,
            "p3_mingc": 30,
            "p3_maxgc": 80,
            "p3_psizeopt": 100,
            "p3_psizemin": 50,
            "p3_psizemax": 150,
            "p3_maxpolyx": 3,
            "p3_osizeopt": 20,
            "p3_ominsize": 13,
            "p3_omaxsize": 30,
            "p3_otmopt": 69,
            "p3_otmmin": 68,
            "p3_otmmax": 70,
            "p3_ogcopt": 55,
            "p3_ogcmin": 30,
            "p3_ogcmax": 80,
            "p3_filter": False,
        }

    @pytest.mark.skipif(PRIMER3_VERSION[0] < 2, reason="requires primer3 v2+")
    def test_primer3_exe(self):
        """Primer3 executable exists and runs, and is version 2 or greater."""
        cmd = [shlex.quote(self.primer3_exe), "--version"]
        # check=False necessary because primer3_core exits with 255 and the
        # informative licence/help statement
        subprocess.run(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )  #  nosec
        #  Check version
        self.assertGreaterEqual(2, PRIMER3_VERSION[0])

    @pytest.mark.skipif(PRIMER3_VERSION[0] < 2, reason="requires primer3 v2+")
    def test_primer3_build_and_run_cmd(self):
        """Primer3 command builds correctly."""
        stem = os.path.join(self.outdir, "GCF_000740965.1_concat")
        cmd = primer3.build_command(
            self.primer3_exe, "Pba_21A", self.seqfile, stem, self.p3_defaults
        )
        target = "primer3_core -output tests/test_output/primer3/GCF_000740965.1_concat.primer3 tests/test_output/primer3/GCF_000740965.1_concat.boulder"
        self.assertEqual(" ".join(cmd.cline), target)
        subprocess.run(
            cmd.cline,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertDirsEqual(self.outdir, self.targetdir)
