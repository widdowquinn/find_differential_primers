#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcommands.py

Test subcommand functions for pdp.py script

This test suite is intended to be run from the repository root using:

nosetests -v

Individual test classes can be run using, e.g.:

$ nosetests -v tests/test_subcommands.py:TestConfigSubcommand

Each command CMD available at the command-line as pdp.py <CMD> is
tested in its own class (subclassing unittest.TestCase), where the
setUp() method defines input/output files, a null logger (picked up
by nosetests), and a dictionary of command lines, keyed by test name
with values that represent the command-line options.

For each test, command-line options are defined in a Namespace,
and passed as the sole argument to the appropriate subcommand
function from subcommands.py.

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
import logging
import os
import unittest

from diagnostic_primers import (blast, )
from diagnostic_primers.scripts import subcommands

from argparse import Namespace
from nose.tools import assert_equal, raises


def ordered(obj):
    """Return ordered version of the passed object

    Dictionaries are not ordered in all Python versions, and the
    implementation of sort_keys() in the the JSON library seems
    erratic in terms of effect
    """
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def assert_dirfiles_equal(dir1, dir2, listonly=False):
    """Assert that dir1 and dir2 contain the same files and they are the same

    - dir1         path to directory of files for comparison
    - dir2         path to directory of files for comparison
    - listonly     Boolean, if True then don't compare file contents

    The list of files present in dir1 and dir2 are compared, and then
    the contents of each file are compared.

    The test is performed differently depending on file extension:

    - .json        test equality of ordered dictionary
    - .fasta       test that file streams are equal
    - .tab         test that file streams are equal
    - .ePrimer3    test that each line after the header is equal
    """
    # Skip hidden files
    dir1files = [_ for _ in os.listdir(dir1) if not _.startswith('.')]
    dir2files = [_ for _ in os.listdir(dir2) if not _.startswith('.')]
    assert_equal(sorted(dir1files), sorted(dir2files),
                 msg="%s and %s have differing contents" % (dir1, dir2))
    if not listonly:
        for fname in dir1files:
            msg = "%s not equal in both directories" % fname
            with open(os.path.join(dir1, fname)) as ofh:
                with open(os.path.join(dir2, fname)) as tfh:
                    if os.path.splitext(fname)[-1] == '.json':
                        # Order the parsed JSON before comparing
                        assert_equal(
                            ordered(json.load(ofh)), ordered(json.load(tfh)),
                            msg=msg)
                    elif os.path.splitext(fname)[-1].lower() == '.eprimer3':
                        # Skip first line
                        assert_equal(ofh.readlines()[
                                     1:], tfh.readlines()[1:], msg=msg)
                    elif os.path.splitext(fname)[-1].lower() == '.blasttab':
                        # Depending on BLAST version, columns may be formatted
                        # differently, so we can't compare characters only
                        data1 = blast.parse_blasttab(ofh)
                        data2 = blast.parse_blasttab(tfh)
                        for line1, line2 in zip(data1, data2):
                            assert_equal(line1, line2)
                    else:  #  Compare the file contents directly
                        assert_equal(ofh.read(), tfh.read(), msg=msg)


class TestConfigSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py config subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'config')
        self.outdir = os.path.join('tests', 'test_output', 'config')
        self.targetdir = os.path.join('tests', 'test_targets', 'config')
        self.tsv_to_json_fname = os.path.join(self.outdir,
                                              'tab_converted_conf.json')
        self.tsv_to_json_target = os.path.join(self.targetdir,
                                               'tab_converted_conf.json')
        self.json_to_tsv_fname = os.path.join(self.outdir,
                                              'json_converted_conf.tab')
        self.json_to_tsv_target = os.path.join(self.targetdir,
                                               'json_converted_conf.tab')
        self.fixed_fname = os.path.join(self.outdir, 'seqfixed_conf.json')
        self.fixed_target = os.path.join(self.targetdir, 'seqfixed_conf.json')

        # null logger instance that does nothing
        self.logger = logging.getLogger('TestConfigSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Dictionary of command-line namespaces
        self.argsdict = {
            'validate_json_good':
            Namespace(
                infilename=os.path.join(self.datadir, 'testconf.json'),
                verbose=True,
                validate=True,
                fix_sequences=False,
                to_json=False,
                to_tab=False),
            'validate_tsv_good':
            Namespace(
                infilename=os.path.join(self.datadir, 'testconf.tab'),
                verbose=True,
                validate=True,
                fix_sequences=False,
                to_json=False,
                to_tab=False),
            'validate_config_bad':
            Namespace(
                infilename=os.path.join(self.datadir, 'testin.conf'),
                verbose=True,
                validate=True,
                fix_sequences=False,
                to_json=False,
                to_tab=False),
            'tab_to_json':
            Namespace(
                infilename=os.path.join(self.datadir, 'testconf.json'),
                verbose=True,
                validate=False,
                fix_sequences=False,
                to_json=self.tsv_to_json_fname,
                to_tab=False),
            'json_to_tab':
            Namespace(
                infilename=os.path.join(self.datadir, 'testconf.tab'),
                verbose=True,
                validate=False,
                fix_sequences=False,
                to_json=False,
                to_tab=self.json_to_tsv_fname),
            'fix_sequences':
            Namespace(
                infilename=os.path.join(self.datadir, 'testconf.json'),
                verbose=True,
                validate=False,
                fix_sequences=self.fixed_fname,
                to_json=False,
                to_tab=False),
            'notconf':
            Namespace(
                infilename=os.path.join(self.datadir, 'testconf.notjson'),
                verbose=True,
                validate=False,
                fix_sequences=self.fixed_fname,
                to_json=False,
                to_tab=False),
        }

    def test_validate_json_good(self):
        """config subcmd validates known good JSON config file."""
        subcommands.subcmd_config(self.argsdict['validate_json_good'],
                                  self.logger)

    def test_validate_tab_good(self):
        """config subcmd validates known good TSV config file."""
        subcommands.subcmd_config(self.argsdict['validate_tsv_good'],
                                  self.logger)

    def test_tsv_to_json(self):
        """config subcmd converts TSV config to JSON."""
        subcommands.subcmd_config(self.argsdict['tab_to_json'], self.logger)
        with open(self.tsv_to_json_fname, 'r') as fh1:
            with open(self.tsv_to_json_target, 'r') as fh2:
                assert_equal(ordered(json.load(fh1)), ordered(json.load(fh2)))

    def test_json_to_tsv(self):
        """config subcmd converts JSON config to TSV."""
        subcommands.subcmd_config(self.argsdict['json_to_tab'], self.logger)
        with open(self.json_to_tsv_fname, 'r') as fh1:
            with open(self.json_to_tsv_target, 'r') as fh2:
                assert_equal(fh1.read(), fh2.read())

    def test_fix_sequences(self):
        """config subcmd fixes sequences and writes JSON."""
        subcommands.subcmd_config(self.argsdict['fix_sequences'], self.logger)

        # Output JSON is correct
        with open(self.fixed_fname, 'r') as fh1:
            with open(self.fixed_target, 'r') as fh2:
                assert_equal(ordered(json.load(fh1)), ordered(json.load(fh2)))

    @raises(SystemExit)
    def test_validate_config_bad(self):
        """config subcmd errors on validating badly-formatted config file."""
        subcommands.subcmd_config(self.argsdict['validate_config_bad'],
                                  self.logger)

    @raises(SystemExit)
    def test_validate_config_bad_suffix(self):
        """config subcmd errors with wrong file extension for config."""
        subcommands.subcmd_config(self.argsdict['notconf'], self.logger)


class TestProdigalSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py prodigal subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.datadir = os.path.join('tests', 'test_input', 'config')
        self.outconfdir = os.path.join('tests', 'test_output', 'config')
        self.outrundir = os.path.join('tests', 'test_output', 'prodigal')
        self.targetdir = os.path.join('tests', 'test_targets', 'prodigal')
        self.prodigal_exe = 'prodigal'
        self.scheduler = 'multiprocessing'
        self.workers = None

        # null logger instance that does nothing
        self.logger = logging.getLogger('TestConfigSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Dictionary of command-line namespaces
        self.argsdict = {
            'run':
            Namespace(
                infilename=os.path.join(self.datadir,
                                        'testreducedep3conf.json'),
                outfilename=os.path.join(self.outconfdir, 'prodconf.json'),
                prodigaldir=self.outrundir,
                prodigal_exe=self.prodigal_exe,
                prodigalforce=True,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=True),
            'notconf':
            Namespace(
                infilename=os.path.join(self.datadir, 'fixedconf.nojson'),
                outfilename=os.path.join(self.outconfdir, 'prodconf.json'),
                prodigaldir=self.outrundir,
                prodigal_exe=self.prodigal_exe,
                prodigalforce=True,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=True),
            'notjson':
            Namespace(
                infilename=os.path.join(self.datadir, 'testin.conf'),
                outfilename=os.path.join(self.outconfdir, 'prodconf.json'),
                prodigaldir=self.outrundir,
                prodigal_exe=self.prodigal_exe,
                prodigalforce=True,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=True),
            'noforce':
            Namespace(
                infilename=os.path.join(self.datadir, 'fixedconf.json'),
                outfilename=os.path.join(self.outconfdir, 'prodconf.json'),
                prodigaldir=self.outrundir,
                prodigal_exe=self.prodigal_exe,
                prodigalforce=False,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=True),
        }

    def test_prodigal_run(self):
        """prodigal subcommand produces correct annotation."""
        subcommands.subcmd_prodigal(self.argsdict['run'], self.logger)

        # Check file contents
        assert_dirfiles_equal(self.outrundir, self.targetdir)

    @raises(SystemExit)
    def test_invalid_conf_file(self):
        """Script exits if prodigal config file has wrong suffix."""
        subcommands.subcmd_prodigal(self.argsdict['notconf'], self.logger)

    @raises(ValueError)
    def test_tsv_conf_file(self):
        """Error raised if .conf file provided for prodigal."""
        subcommands.subcmd_prodigal(self.argsdict['notjson'], self.logger)

    @raises(SystemExit)
    def test_outdir_not_forced(self):
        """Script exits if not forcing prodigal output overwrite."""
        subcommands.subcmd_prodigal(self.argsdict['noforce'], self.logger)
        # We don't enforce test order, so run twice to force the error
        subcommands.subcmd_prodigal(self.argsdict['noforce'], self.logger)


class TestEPrimer3Subcommand(unittest.TestCase):
    """Class defining tests of the pdp.py eprimer3 subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join('tests', 'test_input', 'config')
        self.confoutdir = os.path.join('tests', 'test_output', 'config')
        self.conftargets = os.path.join('tests', 'test_targets', 'config')
        self.outdir = os.path.join('tests', 'test_output', 'eprimer3',
                                   'scriptout')
        self.targetdir = os.path.join('tests', 'test_targets', 'eprimer3')
        self.ep3_exe = 'eprimer3'
        self.scheduler = 'multiprocessing'
        self.workers = None
        self.hybridprobe = False
        # Default values for ePrimer3 run - not modified by any tests,
        # defined in parsers.py
        self.ep3_defaults = {
            'ep_numreturn': 10,
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

        # null logger instance that does nothing
        self.logger = logging.getLogger('TestConfigSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Dictionary of command-line namespaces
        self.argsdict = {
            'run':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprodigalconf.json'),
                outfilename=os.path.join(self.confoutdir, 'ep3conf.json'),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=2,
                verbose=False,
                ep_hybridprobe=self.hybridprobe,
                **self.ep3_defaults),
            'force':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprodigalconf.json'),
                outfilename=os.path.join(self.confoutdir, 'ep3conf.json'),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=2,
                verbose=False,
                ep_hybridprobe=self.hybridprobe,
                **self.ep3_defaults),
            'notconf':
            Namespace(
                infilename=os.path.join(self.confdir,
                                        'testprodigalconf.nojson'),
                outfilename=os.path.join(self.confoutdir, 'ep3conf.json'),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
            'notjson':
            Namespace(
                infilename=os.path.join(self.confdir, 'testin.conf'),
                outfilename=os.path.join(self.confoutdir, 'ep3conf.json'),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=True,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
            'noforce':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprodigalconf.json'),
                outfilename=os.path.join(self.confoutdir, 'ep3conf.json'),
                eprimer3_dir=self.outdir,
                eprimer3_exe=self.ep3_exe,
                eprimer3_force=False,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False,
                ep_hybridprobe=self.hybridprobe,
                **self.ep3_defaults),
        }

    def test_eprimer3_run(self):
        """eprimer3 subcommand executes correct primer design."""
        subcommands.subcmd_eprimer3(self.argsdict['run'], self.logger)
        # Check file contents
        assert_dirfiles_equal(self.outdir, self.targetdir)

    def test_eprimer3_force(self):
        """eprimer3 subcommand executes correctly and overwrites output."""
        subcommands.subcmd_eprimer3(self.argsdict['force'], self.logger)
        # Check file contents
        assert_dirfiles_equal(self.outdir, self.targetdir)

    @raises(SystemExit)
    def test_invalid_conf_file(self):
        """Script exits if ePrimer3 config file has wrong suffix."""
        subcommands.subcmd_eprimer3(self.argsdict['notconf'], self.logger)

    @raises(ValueError)
    def test_tsv_conf_file(self):
        """Error raised if .conf file provided for ePrimer3."""
        subcommands.subcmd_eprimer3(self.argsdict['notjson'], self.logger)

    @raises(SystemExit)
    def test_outdir_not_forced(self):
        """Script exits if not forcing ePrimer3 output overwrite."""
        subcommands.subcmd_eprimer3(self.argsdict['noforce'], self.logger)
        # Run twice to ensure the error is thrown if tests are out of order
        subcommands.subcmd_eprimer3(self.argsdict['noforce'], self.logger)


class TestBlastscreenSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py blastscreen subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.dbdir = os.path.join('tests', 'test_input', 'blastdb')
        self.confdir = os.path.join('tests', 'test_input', 'config')
        self.outconfdir = os.path.join('tests', 'test_output', 'config')
        self.outdir = os.path.join('tests', 'test_output', 'blastscreen')
        self.targetdir = os.path.join('tests', 'test_targets', 'blastscreen')
        self.targetconfdir = os.path.join('tests', 'test_targets', 'config')
        self.blast_exe = 'blastn'
        self.maxaln = 15
        self.scheduler = 'multiprocessing'
        self.workers = None

        # null logger
        self.logger = logging.getLogger('TestBlastscreenSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Command-line namespaces
        self.argsdict = {
            'run':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprimer3conf.json'),
                outfilename=os.path.join(self.outconfdir, 'screened.json'),
                bs_db=os.path.join(self.dbdir, 'e_coli_screen.fna'),
                bs_exe=self.blast_exe,
                bs_force=True,
                bs_dir=self.outdir,
                maxaln=self.maxaln,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
            'noforce':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprimer3conf.json'),
                outfilename=os.path.join(self.outconfdir, 'screened.json'),
                bs_db=os.path.join(self.dbdir, 'e_coli_screen.fna'),
                bs_exe=self.blast_exe,
                bs_force=False,
                bs_dir=self.outdir,
                maxaln=self.maxaln,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
            'nodb':
            Namespace(
                infilename=os.path.join(self.confdir, 'testprimer3conf.json'),
                outfilename=os.path.join(self.outconfdir, 'screened.json'),
                bs_db=None,
                bs_exe=self.blast_exe,
                bs_force=False,
                bs_dir=self.outdir,
                maxaln=self.maxaln,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False),
        }

    def test_blastscreen_run(self):
        """blastscreen command runs normally and overwrites existing folder."""
        subcommands.subcmd_blastscreen(self.argsdict['run'], self.logger)

        # Check file contents: config
        self.logger.info("Checking output config file against target file")
        with open(os.path.join(self.outconfdir, 'screened.json')) as ofh:
            with open(os.path.join(self.targetconfdir,
                                   'screened.json')) as tfh:
                assert_equal(ordered(json.load(ofh)), ordered(json.load(tfh)))
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences:
        self.logger.info("Comparing output sequences/JSON to target")
        assert_dirfiles_equal(self.outdir, self.targetdir)

    @raises(SystemExit)
    def test_blastscreen_noforce(self):
        """blastscreen command does not overwrite existing folder."""
        subcommands.subcmd_blastscreen(self.argsdict['run'], self.logger)
        subcommands.subcmd_blastscreen(self.argsdict['noforce'], self.logger)

    @raises(SystemExit)
    def test_blastscreen_nodb(self):
        """blastscreen command does not run without database."""
        subcommands.subcmd_blastscreen(self.argsdict['nodb'], self.logger)


class TestPrimersearchSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py primersearch subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join('tests', 'test_input', 'config')
        self.outconfdir = os.path.join('tests', 'test_output', 'config')
        self.outdir = os.path.join('tests', 'test_output', 'primerscreen_cmd')
        self.targetdir = os.path.join('tests', 'test_targets',
                                      'primerscreen_cmd')
        self.targetconfdir = os.path.join('tests', 'test_targets', 'config')
        self.confname = 'test_primersearch_cmd.json'
        self.ps_exe = 'primersearch'
        self.mismatchpercent = 0.1  # This must be in range [0,1]
        self.scheduler = 'multiprocessing'
        self.workers = None

        # null logger
        self.logger = logging.getLogger('TestPrimersearchSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Command-line namespaces
        self.argsdict = {
            'run':
            Namespace(
                infilename=os.path.join(self.confdir, self.confname),
                outfilename=os.path.join(self.outconfdir, self.confname),
                ps_exe=self.ps_exe,
                ps_dir=self.outdir,
                ps_force=True,
                mismatchpercent=self.mismatchpercent,
                scheduler=self.scheduler,
                workers=self.workers,
                verbose=False)
        }

    def test_primersearch_run(self):
        """primersearch command runs normally."""
        subcommands.subcmd_primersearch(self.argsdict['run'], self.logger)

        # Check file contents: config
        self.logger.info("Checking output config file against target file")
        with open(os.path.join(self.outconfdir, self.confname)) as ofh:
            with open(
                    os.path.join(self.targetconfdir, self.confname)) as tfh:
                assert_equal(ordered(json.load(ofh)), ordered(json.load(tfh)))
        self.logger.info("Config file checks against target correctly")

        # Check filtered sequences:
        self.logger.info("Comparing output sequences/JSON to target")
        assert_dirfiles_equal(self.outdir, self.targetdir)


class TestClassifySubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py classify subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join('tests', 'test_input', 'config')
        self.outdir = os.path.join('tests', 'test_output', 'classify')
        self.targetdir = os.path.join('tests', 'test_targets', 'classify')

        # null logger
        self.logger = logging.getLogger('TestClassifySubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Command-line Namespaces
        self.argsdict = {
            'run':
            Namespace(
                infilename=os.path.join(self.confdir, 'testclassify.json'),
                outdir=self.outdir,
                cl_force=True,
                verbose=False)
        }

    def test_classify_run(self):
        """Classify command runs normally."""
        subcommands.subcmd_classify(self.argsdict['run'], self.logger)

        # Check output:
        self.logger.info("Comparing output primer sequences to targets")
        assert_dirfiles_equal(self.outdir, self.targetdir)


class TestExtractSubcommand(unittest.TestCase):
    """Class defining tests of the pdp.py extract subcommand."""

    def setUp(self):
        """Set parameters for tests."""
        self.confdir = os.path.join('tests', 'test_input', 'config')
        self.indir = os.path.join('tests', 'test_input', 'extract')
        self.outdir = os.path.join('tests', 'test_output', 'extract')
        self.alignoutdir = os.path.join('tests', 'test_output',
                                        'extract_mafft')
        self.targetdir = os.path.join('tests', 'test_targets', 'extract')
        self.aligntargetdir = os.path.join('tests', 'test_targets',
                                           'extract_mafft')
        self.filestem = "Pectobacterium_primers"
        self.mafft_exe = "mafft"

        # null logger
        self.logger = logging.getLogger('TestExtractSubcommand logger')
        self.logger.addHandler(logging.NullHandler())

        # Command-line Namespaces
        self.argsdict = {
            'run':
            Namespace(
                infilename=os.path.join(self.confdir, 'testclassify.json'),
                primerfile=os.path.join(self.indir, '%s.json' % self.filestem),
                outdir=self.outdir,
                verbose=False,
                ex_force=True,
                noalign=True,
                mafft_exe=self.mafft_exe),
            'align':
            Namespace(
                infilename=os.path.join(self.confdir, 'testclassify.json'),
                primerfile=os.path.join(self.indir, '%s.json' % self.filestem),
                outdir=self.alignoutdir,
                verbose=False,
                ex_force=True,
                noalign=False,
                mafft_exe=self.mafft_exe),
        }

    def test_extract_run(self):
        """Extract command runs normally (no alignment)."""
        args = self.argsdict['run']
        subcommands.subcmd_extract(args, self.logger)

        # Check output:
        self.logger.info("Comparing output amplicons to targets")
        # We have to infer the output location for the extracted amplicons.
        # This is defined by the filestem of the input JSON file
        outputdir = os.path.join(args.outdir, self.filestem)
        targetdir = os.path.join(self.targetdir, self.filestem)
        assert_dirfiles_equal(outputdir, targetdir)

    def test_extract_align(self):
        """Extract command runs normally (with MAFFT alignment)."""
        args = self.argsdict['align']
        subcommands.subcmd_extract(args, self.logger)

        # Check output:
        self.logger.info("Comparing output amplicons to targets")
        # We have to infer the output location for the extracted amplicons.
        # This is defined by the filestem of the input JSON file
        outputdir = os.path.join(args.outdir, self.filestem)
        targetdir = os.path.join(self.aligntargetdir, self.filestem)
        assert_dirfiles_equal(outputdir, targetdir)
