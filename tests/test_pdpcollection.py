#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_pdpcollection.py

Test instantiation, methods and attributions of the PDPCollection class

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

import json
import os
import shutil

from diagnostic_primers.config import (
    PDPData,
    PDPCollection,
    PDPCollectionException,
    PDPEncoder,
    ConfigSyntaxError,
)

import pytest

from tools import PDPTestCase

# Defined as global so it can be seen by the TestBlastscreenSubcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdpcollection")


class TestGenomeCollection(PDPTestCase):
    """Class defining tests of the GenomeCollection object."""

    @classmethod
    def setUpClass(TestGenomeCollection):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join("tests", "test_input", "pdpcollection")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdpcollection")
        os.makedirs(self.outdir, exist_ok=True)
        self.name = "test_collection"
        self.tabconfigfile = os.path.join(self.datadir, "testconf.tab")
        self.jsonconfigfile = os.path.join(self.datadir, "testconf.json")
        self.jsonoutfile = os.path.join(self.outdir, "testconf.json")
        self.taboutfile = os.path.join(self.outdir, "testconf.tab")
        self.tabtargetfile = os.path.join(self.targetdir, "testconf.tab")
        self.failconfig = os.path.join(self.datadir, "broken.conf")
        self.failmissingprimerfile = os.path.join(
            self.datadir, "missingprimerfile.json"
        )

    def test_instantiate(self):
        """PDPCollection instantiates."""
        PDPCollection(self.name)

    def test_load_json(self):
        """PDPCollection loads from JSON config file."""
        gc = PDPCollection(self.name)
        gc.from_json(self.jsonconfigfile)
        self.assertEqual(
            len(gc), 8, msg="Wrong count of sequences in {}".format(self.jsonconfigfile)
        )

    def test_write_json(self):
        """PDPCollection writes JSON config file."""
        gc = PDPCollection(self.name)
        gc.from_json(self.jsonconfigfile)
        gc.write_json(self.jsonoutfile)
        self.assertJsonEqual(self.jsonconfigfile, self.jsonoutfile)

    def test_write_tab(self):
        """PDPCollection writes TSV config file."""
        gc = PDPCollection(self.name)
        gc.from_json(self.jsonconfigfile)
        gc.write_tab(self.taboutfile)
        self.assertFilesEqual(self.taboutfile, self.tabtargetfile)

    def test_json_encoder_not_pdpd(self):
        """PDPEncoder writes non-PDPData object correctly."""
        data = [[1, 2, 3], "string"]
        data_json = json.dumps(data, sort_keys=True, cls=PDPEncoder)
        # String representations from JSON encoders use double-quotes,
        # Python representations use single-quotes
        self.assertEqual(data_json.replace('"', "'"), str(data))

    def test_load_tab(self):
        """PDPCollection loads from TSV config file."""
        gc = PDPCollection(self.name)
        gc.from_tab(self.tabconfigfile)

    def test_load_fail_1(self):
        """PDPCollection throws error with broken TSV config."""
        with pytest.raises(ConfigSyntaxError):
            gc = PDPCollection(self.name)
            gc.from_tab(self.failconfig)

    def test_load_fail_2(self):
        """PDPCollection throws error when loading TSV config as if JSON."""
        with pytest.raises(json.decoder.JSONDecodeError):
            gc = PDPCollection(self.name)
            gc.from_json(self.failconfig)

    def test_collection_data(self):
        """PDPCollection contains PDPData."""
        gc = PDPCollection(self.name)
        gc.from_json(self.jsonconfigfile)
        for item in gc.data:
            self.assertIsInstance(item, PDPData)

    def test_missing_primerfile(self):
        """PDPCollection throws error when primer file not defined."""
        with pytest.raises(PDPCollectionException):
            gc = PDPCollection(self.name)
            gc.from_json(self.failmissingprimerfile)
