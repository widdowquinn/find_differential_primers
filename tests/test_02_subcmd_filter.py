#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_02_subcmd_filter.py

Test prodigal filtering for pdp script

This test suite is intended to be run from the repository root using:

nosetests -v

Individual test classes can be run using, e.g.:

$ nosetests -v tests/test_subcommands.py:TestConfigSubcommand

Each command CMD available at the command-line as pdp <CMD> is
tested in its own class (subclassing unittest.TestCase), where the
setUp() method defines input/output files, a null logger (picked up
by nosetests), and a dictionary of command lines, keyed by test name
with values that represent the command-line options.

For each test, command-line options are defined in a Namespace,
and passed as the sole argument to the appropriate subcommand
function from subcommands.py.

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

import logging
import os
import shutil

from argparse import Namespace
from collections import namedtuple

import pytest

from diagnostic_primers.scripts import subcommands
from diagnostic_primers.scripts.subcommands.subcmd_filter import PDPFilterException

from tools import PDPTestCase, modify_namespace

# Defined as global so it can be seen by the TestFilterSubcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_filter")

# Convenience struct for passing arguments to run filter commands
Params = namedtuple("Params", "infname outfname outdir tgtdir suffix filt")

# Convenience struct for scheduling
Scheduling = namedtuple("Scheduling", "scheduler workers")


class TestFilterSubcommand(PDPTestCase):
    """Class defining tests of the pdp filter subcommand."""

    @classmethod
    def setUpClass(cls):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join("tests", "test_input", "pdp_filter")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdp_filter")
        self.prodigal_exe = "prodigal"
        self.scheduling = Scheduling("multiprocessing", None)

        # null logger instance that does nothing
        self.logger = logging.getLogger("TestFilterSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # base Namespace
        self.base_namespace = Namespace(
            filt_prodigal=False,
            filt_prodigaligr=False,
            filt_alnvar=None,
            filt_outdir=self.outdir,
            filt_prodigal_exe=self.prodigal_exe,
            filt_force=True,
            filt_suffix="prodigal",
            filt_spacerlen=150,
            filt_flanklen=150,
            filt_minserate=0.05,
            filt_minsecount=0,
            scheduler=self.scheduling.scheduler,
            workers=self.scheduling.workers,
            verbose=True,
            disable_tqdm=True,
            recovery=False,
            nucmer_exe="nucmer",
            deltafilter_exe="delta-filter",
            maxmatch=False,
            jobprefix="test_02_subcmd_filter",
        )

    def filter_run(self, params):
        """Runs a pdp filter command with passed settings

        :param params:    Params namedtuple defining command parameters

        Checks the output and target directories for equality
        """
        filt_ns = modify_namespace(
            self.base_namespace,
            {
                "infilename": params.infname,
                "outfilename": params.outfname,
                "filt_outdir": params.outdir,
                "filt_suffix": params.suffix,
            },
        )

        # Decide whether we're using alnvar, prodigal or prodigaligr filters
        if params.filt == "prodigal":
            filt_ns = modify_namespace(filt_ns, {"filt_prodigal": True})
        elif params.filt == "prodigaligr":
            filt_ns = modify_namespace(filt_ns, {"filt_prodigaligr": True})
        elif params.filt == "alnvar":
            filt_ns = modify_namespace(filt_ns, {"filt_alnvar": "atrosepticum_NCBI"})
        subcommands.subcmd_filter(filt_ns, self.logger)

        # Check file contents
        self.assertDirsEqual(params.outdir, params.tgtdir)

    def test_filter_prodigal_run(self):
        """filter subcommand produces correct annotation with --prodigal.

        pdp filter -v --disable_tqdm --prodigal \
            --outdir tests/test_output/pdp_filter/prodigal \
            --suffix prodigal \
            tests/test_input/pdp_filter/seqfixed_conf.json \
            tests/test_output/pdp_filter/prodconf.json
        """
        suffix = "prodigal"
        self.filter_run(
            Params(
                os.path.join(self.datadir, "seqfixed_conf.json"),
                os.path.join(self.outdir, "prodconf.json"),
                os.path.join(self.outdir, suffix),
                os.path.join(self.targetdir, suffix),
                suffix,
                "prodigal",
            )
        )

    def test_filter_prodigaligr_run(self):
        """filter subcommand produces correct annotation with --prodigaligr.

        pdp filter -v --disable_tqdm --prodigal \
            --outdir tests/test_output/pdp_filter/prodigaligr \
            --suffix prodigaligr \
            tests/test_input/pdp_filter/seqfixed_conf.json \
            tests/test_output/pdp_filter/prodigrconf.json
        """
        suffix = "prodigaligr"
        self.filter_run(
            Params(
                os.path.join(self.datadir, "seqfixed_conf.json"),
                os.path.join(self.outdir, "prodigrconf.json"),
                os.path.join(self.outdir, suffix),
                os.path.join(self.targetdir, suffix),
                suffix,
                "prodigaligr",
            )
        )

    def test_filter_alnvar_run(self):
        """filter subcommand produces correct annotation with --alnvar.

        pdp filter -v --disable_tqdm --alnvar betavasculorum_NCBI\
            --outdir tests/test_output/pdp_filter/alnvar \
            --min_sim_error_rate 0.005 \
            --suffix alnvar \
            tests/test_input/pdp_filter/seqfixed_conf.json \
            tests/test_output/pdp_filter/alnvarconf.json
        """
        suffix = "alnvar"
        self.filter_run(
            Params(
                os.path.join(self.datadir, "seqfixed_conf.json"),
                os.path.join(self.outdir, "alnvarconf.json"),
                os.path.join(self.outdir, suffix),
                os.path.join(self.targetdir, suffix),
                suffix,
                "alnvar",
            )
        )

    def test_invalid_conf_file(self):
        """Script exits when filter config file has wrong suffix.

        pdp filter -v --disable_tqdm --prodigal \
            --outdir tests/test_output/pdp_filter \
            --suffix prodigal \
            tests/test_input/pdp_filter/fixedconf.nojson \
            tests / test_output / pdp_config / prodconf.json """
        with pytest.raises(PDPFilterException):
            subcommands.subcmd_filter(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.datadir, "fixedconf.nojson"),
                        "outfilename": os.path.join(self.outdir, "prodconf.json"),
                        "filt_outdir": os.path.join(self.outdir, "prodigal"),
                        "filt_prodigal": True,
                    },
                ),
                self.logger,
            )

    def test_tsv_conf_file(self):
        """Error raised when tab .conf file provided for filter.

        pdp filter -v --disable_tqdm --prodigal \
            --outdir tests/test_output/pdp_filter \
            --suffix prodigal \
            tests/test_input/pdp_filter/testin.conf \
            tests / test_output / pdp_config / prodconf.json """
        with pytest.raises(PDPFilterException):
            subcommands.subcmd_filter(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.datadir, "testin.conf"),
                        "outfilename": os.path.join(self.outdir, "prodconf.json"),
                        "filt_prodigal": True,
                    },
                ),
                self.logger,
            )
