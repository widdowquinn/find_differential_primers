#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_pdpdata.py

Test instantiation, methods and attributions of the PDPData class

This test suite is intended to be run from the repository root using:

nosetests -v

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
import shutil

from collections import namedtuple

import pytest

from Bio import SeqIO

from diagnostic_primers.config import PDPData
from tools import PDPTestCase


class BadName:
    """Dummy object that throws an error if str() is used on it."""

    def __init__(self):
        pass

    def __str__(self):
        raise ValueError("Cannot be a string")


# Defined as global so it can be seen by the TestPDPData() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdpdata")


# Convenience struct for TestPDPData groups
PDPTestGroups = namedtuple("TestGroups", "list str set")

# Convenience struct for TestPDPData directory paths
PDPTestDirs = namedtuple("TestDirs", "outdir tgtdir datadir")

# Convenience struct for TestPDPData test sequence files
PDPTestSeqFiles = namedtuple(
    "TestSeqFiles", "seqfile filtered_seqfile stitchfile stitchtgt noambigtgt"
)

# Convenience struct for TestPDPData test data files
PDPTestDataFiles = namedtuple(
    "TestDataFiles", "features primers primersearch target_amplicons"
)


class TestPDPData(PDPTestCase):
    """Class defining tests of the PDPData object."""

    @classmethod
    def setUpClass(cls):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.testdirs = PDPTestDirs(
            OUTDIR,
            os.path.join("tests", "test_targets", "pdpdata"),
            os.path.join("tests", "test_input", "pdpdata"),
        )
        os.makedirs(self.testdirs.outdir, exist_ok=True)
        self.name = "test_name"
        self.testgroups = PDPTestGroups(
            ["group1", "group2", "group3"],
            "group1,group2,group3",
            set(["group1", "group2", "group3"]),
        )
        self.testseqfiles = PDPTestSeqFiles(
            os.path.join(self.testdirs.datadir, "GCF_000011605.1.fasta"),
            os.path.join(self.testdirs.datadir, "GCF_000011605.1_prodigal.fasta"),
            os.path.join(self.testdirs.datadir, "test_stitch.fasta"),
            os.path.join(self.testdirs.tgtdir, "test_stitch_concat.fas"),
            os.path.join(self.testdirs.tgtdir, "test_stitch_noambig.fas"),
        )
        self.testdatafiles = PDPTestDataFiles(None, None, None, None)

    def test_instantiation_grouplist(self):
        """PDPData object instantiates with list of groups."""
        PDPData(
            self.name,
            self.testgroups.list,
            self.testseqfiles.seqfile,
            self.testseqfiles.filtered_seqfile,
            self.testdatafiles.features,
            self.testdatafiles.primers,
            self.testdatafiles.primersearch,
            self.testdatafiles.target_amplicons,
        )

    def test_instantiation_groupstr(self):
        """PDPData object instantiates with string of groups."""
        PDPData(
            self.name,
            self.testgroups.str,
            self.testseqfiles.seqfile,
            self.testseqfiles.filtered_seqfile,
            self.testdatafiles.features,
            self.testdatafiles.primers,
            self.testdatafiles.primersearch,
            self.testdatafiles.target_amplicons,
        )

    def test_instantiation_groupset(self):
        """PDPData object instantiates with set of groups."""
        PDPData(
            self.name,
            self.testgroups.set,
            self.testseqfiles.seqfile,
            self.testseqfiles.filtered_seqfile,
            self.testdatafiles.features,
            self.testdatafiles.primers,
            self.testdatafiles.primersearch,
            self.testdatafiles.target_amplicons,
        )

    def test_stitch(self):
        """PDPData object stitches multi-sequence input.

        The self.stitchfile path points to an input file with multiple
        sequences. We test if this needs stitching (it should), and if so
        stitch it.
        """
        gdata = PDPData(
            self.name,
            self.testgroups.str,
            self.testseqfiles.stitchfile,
            self.testseqfiles.filtered_seqfile,
            self.testdatafiles.features,
            self.testdatafiles.primers,
            self.testdatafiles.primersearch,
            self.testdatafiles.target_amplicons,
        )
        gdata.stitch(outdir=self.testdirs.outdir)  # assumes stitch is needed
        # We have to take the input filename from the GenomeData object,
        # as this will have changed if stitched/ambiguities removed.
        self.assertFilesEqual(gdata.seqfile, self.testseqfiles.stitchtgt)

    def test_noambig(self):
        """PDPData object replaces ambiguities in input.

        The self.stitchfile path points to an input file with multiple
        sequences, and ambiguity symbols. We test for ambiguity symbols, and
        replace them.
        """
        gdata = PDPData(
            self.name,
            self.testgroups.str,
            self.testseqfiles.stitchfile,
            self.testseqfiles.filtered_seqfile,
            self.testdatafiles.features,
            self.testdatafiles.primers,
            self.testdatafiles.primersearch,
            self.testdatafiles.target_amplicons,
        )
        if gdata.has_ambiguities:
            gdata.replace_ambiguities(outdir=self.testdirs.outdir)
        # We have to take the input filename from the GenomeData object,
        # as this will have changed if stitched/ambiguities removed.
        self.assertFilesEqual(gdata.seqfile, self.testseqfiles.noambigtgt)

    def test_attributes(self):
        """PDPData attributes are accessible.

        Ensure lazily-populated attributes get assigned
        """
        gdata = PDPData(
            self.name,
            self.testgroups.str,
            self.testseqfiles.seqfile,
            self.testseqfiles.filtered_seqfile,
            self.testdatafiles.features,
            self.testdatafiles.primers,
            self.testdatafiles.primersearch,
            self.testdatafiles.target_amplicons,
        )
        self.assertEqual(
            gdata.seqnames,
            [s.id for s in SeqIO.parse(self.testseqfiles.seqfile, "fasta")],
        )
        self.assertEqual(gdata.features, self.testdatafiles.features)

    def test_invalid_sequence(self):
        """PDPData errors with invalid input sequence file."""
        with pytest.raises(OSError):
            PDPData(
                self.name,
                self.testgroups.str,
                "seqfile.notexist",
                self.testseqfiles.filtered_seqfile,
                self.testdatafiles.features,
                self.testdatafiles.primers,
                self.testdatafiles.primersearch,
                self.testdatafiles.target_amplicons,
            )

    def test_invalid_filtered_sequence(self):
        """PDPData errors with invalid input filtered sequence file."""
        with pytest.raises(OSError):
            PDPData(
                self.name,
                self.testgroups.str,
                self.testseqfiles.seqfile,
                "filtered_seqfile.notexist",
                self.testdatafiles.features,
                self.testdatafiles.primers,
                self.testdatafiles.primersearch,
                self.testdatafiles.target_amplicons,
            )

    def test_invalid_features(self):
        """PDPData errors with invalid feature file."""
        with pytest.raises(OSError):
            PDPData(
                self.name,
                self.testgroups.str,
                self.testseqfiles.seqfile,
                self.testseqfiles.filtered_seqfile,
                "features.notexist",
                self.testdatafiles.primers,
                self.testdatafiles.primersearch,
                self.testdatafiles.target_amplicons,
            )

    def test_invalid_primers(self):
        """PDPData errors with invalid primers file."""
        with pytest.raises(OSError):
            PDPData(
                self.name,
                self.testgroups.str,
                self.testseqfiles.seqfile,
                self.testseqfiles.filtered_seqfile,
                self.testdatafiles.features,
                "primers.notexist",
                self.testdatafiles.primersearch,
                self.testdatafiles.target_amplicons,
            )

    def test_invalid_groups(self):
        """PDPData errors with invalid group type."""
        with pytest.raises(TypeError):
            PDPData(
                self.name,
                12345,
                self.testseqfiles.seqfile,
                self.testseqfiles.filtered_seqfile,
                self.testdatafiles.features,
                self.testdatafiles.primers,
                self.testdatafiles.primersearch,
                self.testdatafiles.target_amplicons,
            )

    def test_invalid_name(self):
        """PDPData errors because name is not/cannot be a string."""
        with pytest.raises(TypeError):
            PDPData(
                BadName(),
                self.testgroups.str,
                self.testseqfiles.seqfile,
                self.testseqfiles.filtered_seqfile,
                self.testdatafiles.features,
                self.testdatafiles.primers,
                self.testdatafiles.primersearch,
                self.testdatafiles.target_amplicons,
            )
