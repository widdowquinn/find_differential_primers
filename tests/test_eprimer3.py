#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_eprimer3.py

Test generation of ePrimer3 command-lines, primer file parsing and writing.

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
import os
import subprocess
import sys
import unittest

from Bio.Emboss import Primer3
from nose.tools import assert_equal, raises, nottest

from diagnostic_primers import (eprimer3, config)


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    elif isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    elif isinstance(obj, Primer3.Primers):
        return sorted(ordered(x) for x in obj.__dict__.items())
    else:
        return obj


class TestCommands(unittest.TestCase):

    """Class defining tests of ePrimer3 command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        self.ep3_exe = 'eprimer3'
        self.datadir = os.path.join('tests', 'test_input', 'sequences')
        self.outdir = os.path.join('tests', 'test_output', 'eprimer3')
        self.targetdir = os.path.join('tests', 'test_targets', 'eprimer3')
        self.config = os.path.join('tests', 'test_input', 'config',
                                   'testprodigalconf.json')
        self.seqfile = os.path.join(self.datadir, "GCF_000011605.1.fasta")
        # Default values for ePrimer3 run - not modified by any tests,
        # defined in parsers.py
        self.ep3_defaults = {'ep_numreturn': 10,
                             'ep_osize': 20,
                             'ep_minsize': 18,
                             'ep_maxsize': 22,
                             'ep_opttm': 59,
                             'ep_mintm': 58,
                             'ep_maxtm': 60,
                             'ep_ogcpercent': 55,
                             'ep_mingc': 30,
                             'ep_maxgc': 80,
                             'ep_psizeopt': 100,
                             'ep_psizemin': 50,
                             'ep_psizemax': 150,
                             'ep_maxpolyx': 3,
                             'ep_osizeopt': 20,
                             'ep_ominsize': 13,
                             'ep_omaxsize': 30,
                             'ep_otmopt': 69,
                             'ep_otmmin': 68,
                             'ep_otmmax': 70,
                             'ep_ogcopt': 55,
                             'ep_ogcmin': 30,
                             'ep_ogcmax': 80,
                             }

    def test_eprimer3_exe(self):
        """ePrimer3 executable exists and runs."""
        cmd = "{0} --help".format(self.ep3_exe)
        result = subprocess.run(cmd, shell=sys.platform != "win32",
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=True)
        print(cmd)
        # EMBOSS writes information out to STDERR
        assert_equal(result.stderr[:6], b'EMBOSS')

    def test_eprimer3_cmd(self):
        """ePrimer3 primer creation command builds correctly."""
        filestem = os.path.join(self.outdir,
                                os.path.split(self.seqfile)[-1])
        cmd = eprimer3.build_command(self.ep3_exe, self.seqfile,
                                     filestem, self.ep3_defaults)
        target = ' '.join(["eprimer3 -auto",
                           "-outfile=tests/test_output/eprimer3/GCF_000011605.1.fasta.eprimer3",
                           "-sequence=tests/test_input/sequences/GCF_000011605.1.fasta",
                           "-numreturn=10", "-osize=20", "-minsize=18", "-maxsize=22",
                           "-opttm=59", "-mintm=58", "-maxtm=60", "-ogcpercent=55",
                           "-mingc=30", "-maxgc=80", "-maxpolyx=3", "-psizeopt=100",
                           "-prange=50-150", "-osizeopt=20", "-ominsize=13",
                           "-omaxsize=30", "-otmopt=69", "-otmmin=68", "-otmmax=70",
                           "-ogcopt=55", "-ogcmin=30", "-ogcmax=80"])
        assert_equal(str(cmd), target)

    def test_eprimer3_cmds(self):
        """ePrimer3 primer creation commands build with no errors."""
        pdpc = config.PDPCollection()
        pdpc.from_json(self.config)
        clines = eprimer3.build_commands(pdpc, self.ep3_exe, self.outdir,
                                         self.ep3_defaults)


class TestParsing(unittest.TestCase):

    """Class defining tests of primer file parsing."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'eprimer3')
        self.outdir = os.path.join('tests', 'test_output', 'eprimer3')
        self.targetdir = os.path.join('tests', 'test_targets', 'eprimer3')
        # The three paths below should point to the same data in three
        # different formats
        self.ep3primerfile = os.path.join(self.datadir, "GCF_000011605.1.eprimer3")
        self.ep3extprimerfile = os.path.join(self.datadir,
                                             "GCF_000011605.1_named.eprimer3")
        self.jsonprimerfile = os.path.join(self.datadir,
                                           "GCF_000011605.1_named.json")
        self.fastaprimerfile = os.path.join(self.datadir,
                                            "GCF_000011605.1_named.fasta")
        # The target files below should contain the same data as the three
        # files above (but not the same format)
        with open(self.ep3primerfile, 'r') as tfh1:     # bare ePrimer3
            self.ep3primertargets = Primer3.read(tfh1).primers
        with open(self.ep3primerfile, 'r') as tfh2:  # named ePrimer3
            self.namedprimertargets = Primer3.read(tfh2).primers
            # We need to add a name to each primer, to emulate the named input
            # (Biopython's parser does not record this)
            for idx, primer in enumerate(self.namedprimertargets, 1):
                stem = os.path.splitext(os.path.split(self.ep3primerfile)[-1])[0]
                primer.name = "%s_primer_%05d" % (stem, idx)

    def test_load_primers_eprimer3(self):
        """ePrimer3 format primers load correctly."""
        primers = eprimer3.load_primers(self.ep3primerfile,
                                        format="eprimer3",
                                        noname=True)
        for primer1, primer2 in zip(primers, self.ep3primertargets):
            assert_equal(ordered(primer1), ordered(primer2))

    def test_load_primers_eprimer3extended(self):
        """ePrimer3 extended format primers load without error."""
        primers = eprimer3.load_primers(self.ep3extprimerfile,
                                        format="eprimer3",
                                        noname=True)
        for primer1, primer2 in zip(primers, self.ep3primertargets):
            assert_equal(ordered(primer1), ordered(primer2))

    def test_load_primers_json(self):
        """JSON format primers load without error."""
        primers = eprimer3.load_primers(self.jsonprimerfile,
                                        format="json")
        for primer1, primer2 in zip(primers, self.namedprimertargets):
            assert_equal(ordered(primer1), ordered(primer2))

    def test_write_primers_eprimer3(self):
        """parse primers and write in ePrimer3 format."""
        primers = eprimer3.load_primers(self.ep3primerfile,
                                        format="eprimer3")
        outfname = os.path.join(self.outdir, "test_write_primers.eprimer3")
        eprimer3.write_primers(primers, outfname, format="eprimer3")
        with open(outfname, 'r') as wfh:
            with open(self.ep3extprimerfile, 'r') as tfh:
                # We need to skip the first comment line as the file
                # paths are different
                assert_equal(wfh.readlines()[1:],
                             tfh.readlines()[1:])

    def test_write_primers_json(self):
        """parse primers and write in JSON format."""
        primers = eprimer3.load_primers(self.ep3primerfile,
                                        format="ep3")
        outfname = os.path.join(self.outdir, "test_write_primers.json")
        eprimer3.write_primers(primers, outfname, format="json")
        with open(outfname, 'r') as wfh:
            with open(self.jsonprimerfile, 'r') as tfh:
                assert_equal(ordered(json.load(wfh)),
                             ordered(json.load(tfh)))

    def test_write_primers_fasta(self):
        """parse primers and write in FASTA format."""
        primers = eprimer3.load_primers(self.jsonprimerfile,
                                        format="json")
        outfname = os.path.join(self.outdir, "test_write_primers.fasta")
        eprimer3.write_primers(primers, outfname, format="fasta")
        with open(outfname, 'r') as wfh:
            with open(self.fastaprimerfile, 'r') as tfh:
                assert_equal(wfh.read(), tfh.read())
