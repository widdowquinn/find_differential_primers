#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_eprimer3.py

Test generation of ePrimer3 command-lines, primer file parsing and writing.

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

from Bio.Emboss import Primer3

from diagnostic_primers import config, eprimer3, load_primers, write_primers

from tools import PDPTestCase


# Defined as global so it can be seen by the TestCommands() and TestParsing() classes
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "eprimer3")


class TestCommands(PDPTestCase):
    """Class defining tests of ePrimer3 command-line generation."""

    @classmethod
    def setUpClass(TestCommands):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.ep3_exe = "eprimer3"
        self.datadir = os.path.join("tests", "test_input", "eprimer3")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "eprimer3")
        self.config = os.path.join(
            "tests", "test_input", "eprimer3", "testprodigalconf.json"
        )
        self.seqfile = os.path.join(self.datadir, "GCF_000011605.1.fasta")
        self.existingfiles = []
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
            "ep_filter": False,
            "recovery": False,
        }

    def test_eprimer3_exe(self):
        """ePrimer3 executable exists and runs."""
        cmd = [shlex.quote(self.ep3_exe), "--version"]
        pipe = subprocess.run(
            cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )  #  nosec
        # EMBOSS writes information out to STDERR
        self.assertEqual(pipe.stderr[:6], b"EMBOSS")

    def test_eprimer3_cmd(self):
        """ePrimer3 primer creation command builds correctly."""
        filestem = os.path.join(self.outdir, os.path.split(self.seqfile)[-1])
        cmd = eprimer3.build_command(
            self.ep3_exe, self.seqfile, filestem, self.ep3_defaults
        )
        target = " ".join(
            [
                "eprimer3 -auto",
                "-outfile="
                + "tests/test_output/eprimer3/GCF_000011605.1.fasta.eprimer3",
                "-sequence=tests/test_input/eprimer3/GCF_000011605.1.fasta",
                "-numreturn=10",
                "-osize=20",
                "-minsize=18",
                "-maxsize=22",
                "-opttm=59",
                "-mintm=58",
                "-maxtm=60",
                "-ogcpercent=55",
                "-mingc=30",
                "-maxgc=80",
                "-maxpolyx=3",
                "-psizeopt=100",
                "-prange=50-150",
                "-osizeopt=20",
                "-ominsize=13",
                "-omaxsize=30",
                "-otmopt=69",
                "-otmmin=68",
                "-otmmax=70",
                "-ogcopt=55",
                "-ogcmin=30",
                "-ogcmax=80",
            ]
        )
        self.assertEqual(str(cmd), target)

    def test_eprimer3_cmds(self):
        """ePrimer3 primer creation commands build with no errors."""
        pdpc = config.PDPCollection()
        pdpc.from_json(self.config)
        eprimer3.build_commands(
            pdpc, self.ep3_exe, self.outdir, self.existingfiles, self.ep3_defaults
        )


class TestParsing(PDPTestCase):
    """Class defining tests of primer file parsing."""

    @classmethod
    def setUpClass(TestParsing):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join("tests", "test_input", "eprimer3")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "eprimer3")
        # The three paths below should point to the same data in three
        # different formats
        self.ep3primerfile = os.path.join(self.datadir, "GCF_000011605.1.eprimer3")
        self.ep3extprimerfile_in = os.path.join(
            self.datadir, "GCF_000011605.1_named.eprimer3"
        )
        self.ep3extprimerfile_tgt = os.path.join(
            self.targetdir, "GCF_000011605.1_named.eprimer3"
        )
        self.jsonprimerfile = os.path.join(self.datadir, "GCF_000011605.1_named.json")
        self.fastaprimerfile = os.path.join(self.datadir, "GCF_000011605.1_named.fasta")
        # The target files below should contain the same data as the three
        # files above (but not the same format)
        with open(self.ep3primerfile, "r") as tfh1:  # bare ePrimer3
            self.ep3primertargets = Primer3.read(tfh1).primers
        with open(self.ep3primerfile, "r") as tfh2:  # named ePrimer3
            self.namedprimertargets = Primer3.read(tfh2).primers
            # We need to add a name to each primer, to emulate the named input
            # (Biopython's parser does not record this)
            for idx, primer in enumerate(self.namedprimertargets, 1):
                stem = os.path.splitext(os.path.split(self.ep3primerfile)[-1])[0]
                primer.name = "%s_primer_%05d" % (stem, idx)

    def test_load_primers_eprimer3(self):
        """ePrimer3 format primers load correctly."""
        primers = load_primers(self.ep3primerfile, fmt="eprimer3", noname=True)
        for primer1, primer2 in zip(primers, self.ep3primertargets):
            self.assertDictEqual(primer1.__dict__, primer2.__dict__)

    def test_load_primers_eprimer3extended(self):
        """ePrimer3 extended format primers load without error."""
        primers = load_primers(self.ep3extprimerfile_in, fmt="eprimer3", noname=True)
        for primer1, primer2 in zip(primers, self.ep3primertargets):
            self.assertDictEqual(primer1.__dict__, primer2.__dict__)

    def test_load_primers_json(self):
        """JSON format primers load without error."""
        primers = load_primers(self.jsonprimerfile, fmt="json")
        for primer1, primer2 in zip(primers, self.namedprimertargets):
            self.assertDictEqual(primer1.__dict__, primer2.__dict__)

    def test_write_primers_eprimer3(self):
        """parse primers and write in ePrimer3 format."""
        primers = load_primers(self.ep3primerfile, fmt="eprimer3")
        outfname = os.path.join(self.outdir, "test_write_primers.eprimer3")
        write_primers(primers, outfname, fmt="eprimer3")
        self.assertEprimer3Equal(outfname, self.ep3extprimerfile_tgt)

    def test_write_primers_json(self):
        """parse primers and write in JSON format."""
        primers = load_primers(self.ep3primerfile, fmt="ep3")
        outfname = os.path.join(self.outdir, "test_write_primers.json")
        write_primers(primers, outfname, fmt="json")
        self.assertJsonEqual(outfname, self.jsonprimerfile)

    def test_write_primers_fasta(self):
        """parse primers and write in FASTA format."""
        primers = load_primers(self.jsonprimerfile, fmt="json")
        outfname = os.path.join(self.outdir, "test_write_primers.fasta")
        write_primers(primers, outfname, fmt="fasta")
        self.assertFilesEqual(outfname, self.fastaprimerfile)
