#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_01_subcmd_config.py

Test config subcommand for pdp script

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

import logging
import os
import shutil

from argparse import Namespace

import pytest

from diagnostic_primers.scripts import subcommands

from tools import PDPTestCase, modify_namespace

# Defined as global so it can be seen by the TestConfigSubcommand() class
# setUpClass() classmethod.
OUTDIR = os.path.join("tests", "test_output", "pdp_config")


class TestConfigSubcommand(PDPTestCase):
    """Class defining tests of the pdp config subcommand."""

    @classmethod
    def setUpClass(TestConfigSubcommand):
        # Clean up old output directory
        if os.path.isdir(OUTDIR):
            shutil.rmtree(OUTDIR)

    def setUp(self):
        """Set parameters for tests."""
        self.indir = os.path.join("tests", "test_input", "pdp_config")
        self.outdir = OUTDIR
        self.targetdir = os.path.join("tests", "test_targets", "pdp_config")

        # null logger instance that does nothing
        self.logger = logging.getLogger("TestConfigSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # base Namespace
        self.base_namespace = Namespace(
            outdir=self.outdir,
            verbose=True,  # verbose on
            disable_tqdm=True,  # turn tqdm progress bar off
            validate=False,  # no other options selected
            fix_sequences=False,
            to_json=False,
            to_tab=False,
        )

    def test_validate_json_good(self):
        """config subcmd validates known good JSON config file.

        pdp config -v --disable_tqdm --validate \
            tests/test_input/config/testconf.json \
            --outdir tests/test_output/pdp_config
        """
        subcommands.subcmd_config(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "testconf.json"),
                    "validate": True,
                },
            ),
            self.logger,
        )

    def test_validate_tab_good(self):
        """config subcmd validates known good TSV config file.

        pdp config -v --disable_tqdm --validate \
            tests/test_input/config/testconf.tab \
            --outdir tests/test_output/pdp_config
        """
        subcommands.subcmd_config(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "testconf.tab"),
                    "validate": True,
                },
            ),
            self.logger,
        )

    def test_tsv_to_json(self):
        """config subcmd converts TSV config to JSON.

        pdp config -v --disable_tqdm \
            --to_json tests/test_output/pdp_config/json_converted_conf.json \
            tests/test_input/config/testconf.tab \
            --outdir tests/test_output/pdp_config
        """
        conf_out = os.path.join(self.outdir, "tab_converted_conf.json")
        conf_tgt = os.path.join(self.targetdir, "tab_converted_conf.json")
        subcommands.subcmd_config(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "testconf.tab"),
                    "to_json": os.path.join(conf_out),
                },
            ),
            self.logger,
        )
        self.assertJsonEqual(conf_out, conf_tgt)

    def test_json_to_tsv(self):
        """config subcmd converts JSON config to TSV.

        pdp config -v --disable_tqdm \
            --to_tab tests/test_output/pdp_config/json_converted_conf.tab \
            tests/test_input/config/testconf.json \
            --outdir tests/test_output/pdp_config
        """
        conf_out = os.path.join(self.outdir, "json_converted_conf.tab")
        conf_tgt = os.path.join(self.targetdir, "json_converted_conf.tab")
        subcommands.subcmd_config(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "testconf.json"),
                    "to_tab": conf_out,
                },
            ),
            self.logger,
        )
        self.assertFilesEqual(conf_out, conf_tgt)

    def test_fix_sequences(self):
        """config subcmd fixes sequences and writes JSON/TSV.

        pdp config -v --disable_tqdm \
            --fix_sequences tests/test_output/pdp_config/seqfixed_conf.json \
            tests/test_input/config/testconf.json \
            --outdir tests/test_output/pdp_config
        """
        # Fix sequences
        conf_out = os.path.join(self.outdir, "seqfixed_conf.json")
        conf_tgt = os.path.join(self.targetdir, "seqfixed_conf.json")
        subcommands.subcmd_config(
            modify_namespace(
                self.base_namespace,
                {
                    "infilename": os.path.join(self.indir, "testconf.json"),
                    "fix_sequences": conf_out,
                },
            ),
            self.logger,
        )
        self.assertJsonEqual(conf_out, conf_tgt)
        # Write equivalent tabular output
        tabconf_out = os.path.join(self.outdir, "seqfixed_conf.tab")
        tabconf_tgt = os.path.join(self.targetdir, "seqfixed_conf.tab")
        subcommands.subcmd_config(
            modify_namespace(
                self.base_namespace, {"infilename": conf_out, "to_tab": tabconf_out}
            ),
            self.logger,
        )
        self.assertFilesEqual(tabconf_out, tabconf_tgt)

    def test_validate_config_bad(self):
        """config subcmd errors on validating badly-formatted config file.

        pdp config -v --disable_tqdm --validate \
            tests/test_input/config/testin.conf \
            --outdir tests/test_output/pdp_config
        """
        with pytest.raises(SystemExit):
            subcommands.subcmd_config(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.indir, "testin.conf"),
                        "validate": True,
                    },
                ),
                self.logger,
            )

    def test_validate_config_bad_suffix(self):
        """config subcmd errors with wrong file extension for config.

        pdp config -v --disable_tqdm \
            --fix_sequences tests/test_output/pdp_config/seqfixed_conf.json \
            tests/test_input/config/testconf.notjson \
            --outdir tests/test_output/pdp_config
        """
        with pytest.raises(SystemExit):
            subcommands.subcmd_config(
                modify_namespace(
                    self.base_namespace,
                    {
                        "infilename": os.path.join(self.indir, "testconf.notjson"),
                        "fix_sequences": os.path.join(
                            self.outdir, "seqfixed_conf.json"
                        ),
                    },
                ),
                self.logger,
            )
