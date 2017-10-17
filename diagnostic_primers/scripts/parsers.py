#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""parsers.py

Provides functions to generate subcommand parsers for the pdp.py script

- config:        interact with config files
- prodigal:      run prodigal bacterial gene prediction
- eprimer3:      predict primers with EMBOSS ePrimer3
- primersearch:  run in-silico hybridisation with EMBOSS PrimerSearch
- blastscreen:   screen primers with BLASTN
- classify:      classify primers

(c) The James Hutton Institute 2017

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

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

import sys

from argparse import ArgumentParser

from . import subcommands


# Build common parser for all subcommands
def build_common_parser():
    """Returns the common argument parser for the script.

    This parser implements options that are common to all subcommands.
    """
    parser_common = ArgumentParser(add_help=False)
    parser_common.add_argument('infilename',
                               help='path to configuration file')
    parser_common.add_argument('-l', '--logfile', dest='logfile',
                               action='store', default=None,
                               help='logfile location')
    parser_common.add_argument('-v', '--verbose', action='store_true',
                               dest='verbose', default=False,
                               help='report progress to log')
    return parser_common


# Build common parser for commands that use a scheduler
def build_scheduler_parser():
    """Returns a parser with options for a scheduler.

    This parser implements options that are common to commands that use
    a scheduler.
    """
    parser_scheduler = ArgumentParser(add_help=False)
    parser_scheduler.add_argument('-s', '--scheduler', dest='scheduler',
                                  action='store', default='multiprocessing',
                                  help='Job scheduler [multiprocessing|SGE]')
    parser_scheduler.add_argument('-w', '--workers', dest='workers',
                                  action='store', default=None, type=int,
                                  help='Number of parallel workers to use')
    return parser_scheduler


# Build subcommand parsers
# Each subcommand gets its own parser, and the entry point to the subcommand
# is defined with parser.set_defaults().
# See https://docs.python.org/3.6/library/argparse.html#sub-commands
def build_parser_config(subparsers, parents=None):
    """Add parser for `config` subcommand to subparsers

    This parser implements options for processing configuration files.
    """
    parser = subparsers.add_parser('config', parents=parents)
    parser.add_argument('--validate', action='store_true',
                        dest='validate', default=False,
                        help='Validate config file, then exit')
    parser.add_argument('--fix_sequences', action='store',
                        dest='fix_sequences', default=None,
                        help='Fix config file sequences and write new JSON ' +
                        'config file')
    parser.add_argument('--to_json', action='store',
                        dest='to_json', default=None,
                        help='Convert .tab config file to JSON and write')
    parser.add_argument('--to_tab', action='store',
                        dest='to_tab', default=None,
                        help='Convert JSON config file to .tab and write')
    parser.set_defaults(func=subcommands.subcmd_config)


def build_parser_prodigal(subparsers, parents=None):
    """Add parser for `prodigal` subcommand to subparsers

    This parser implements options for controlling the prodigal bacterial
    gene prediction tool.
    """
    parser = subparsers.add_parser('prodigal', aliases=['prod'],
                                   parents=parents)
    parser.add_argument('outfilename',
                        help='Path to write new configuration file')
    parser.add_argument('--prodigal', dest='prodigal_exe',
                        action='store', default='prodigal',
                        help='path to Prodigal executable')
    parser.add_argument('--outdir', dest='prodigaldir',
                        action='store', default='prodigal',
                        help='path to directory for Prodigal output')
    parser.add_argument('-f', '--force', dest='prodigalforce',
                        action='store_true', default=False,
                        help='Allow overwrite in Prodigal output directory')
    parser.set_defaults(func=subcommands.subcmd_prodigal)


def build_parser_eprimer3(subparsers, parents=None):
    """Add parser for `eprimer3` subcommand to subparsers

    This parser implements options for controlling primer creation with the
    EMBOSS ePrimer3 tool.
    """
    parser = subparsers.add_parser('eprimer3', aliases=['e3'],
                                   parents=parents)
    # Primer prediction options - subcommand eprimer3
    parser.add_argument('outfilename',
                        help='Path to write new configuration file')
    parser.add_argument('--eprimer3', dest='eprimer3_exe',
                        action='store', default='eprimer3',
                        help='path to ePrimer3 executable')
    parser.add_argument('--outdir', dest='eprimer3_dir',
                        action='store', default='eprimer3',
                        help='path to directory for ePrimer3 output')
    parser.add_argument('-f', '--force', dest='eprimer3_force',
                        action='store_true', default=False,
                        help='Overwrite old ePrimer3 output')
    parser.add_argument('--numreturn', dest='ep_numreturn',
                        action='store', default=10, type=int,
                        help='number of primers to return')
    parser.add_argument('--osize', dest='ep_osize',
                        action='store', default=20, type=int,
                        help='optimal size for primer oligo')
    parser.add_argument('--minsize', dest='ep_minsize',
                        action='store', default=18, type=int,
                        help='minimum size for primer oligo')
    parser.add_argument('--maxsize', dest='ep_maxsize',
                        action='store', default=22, type=int,
                        help='maximum size for primer oligo')
    parser.add_argument('--opttm', dest='ep_opttm',
                        action='store', default=59, type=int,
                        help='optimal Tm for primer oligo')
    parser.add_argument('--mintm', dest='ep_mintm',
                        action='store', default=58, type=int,
                        help='minimum Tm for primer oligo')
    parser.add_argument('--maxtm', dest='ep_maxtm',
                        action='store', default=60, type=int,
                        help='maximum Tm for primer oligo')
    parser.add_argument('--ogcpercent', dest='ep_ogcpercent',
                        action='store', default=55, type=int,
                        help='optimal %%GC for primer oligo')
    parser.add_argument('--mingcpercent', dest='ep_mingc',
                        action='store', default=30, type=int,
                        help='minimum %%GC for primer oligo')
    parser.add_argument('--maxgcpercent', dest='ep_maxgc',
                        action='store', default=80, type=int,
                        help='maximum %%GC for primer oligo')
    parser.add_argument('--psizeopt', dest='ep_psizeopt',
                        action='store', default=100, type=int,
                        help='optimal size of amplified region')
    parser.add_argument('--psizemin', dest='ep_psizemin',
                        action='store', default=50, type=int,
                        help='minimum size of amplified region')
    parser.add_argument('--psizemax', dest='ep_psizemax',
                        action='store', default=150, type=int,
                        help='maximum size of amplified region')
    parser.add_argument('--maxpolyx', dest='ep_maxpolyx',
                        action='store', default=3, type=int,
                        help='maximum run of repeated nucleotides ' +
                        'in primer')
    parser.add_argument('--hybridprobe', dest='ep_hybridprobe',
                        action='store_true', default=False,
                        help='design a reporter oligo')
    parser.add_argument('--oligoosize', dest='ep_osizeopt',
                        action='store', default=20, type=int,
                        help='optimal size for internal oligo')
    parser.add_argument('--oligominsize', dest='ep_ominsize',
                        action='store', default=13, type=int,
                        help='minimum size for internal oligo')
    parser.add_argument('--oligomaxsize', dest='ep_omaxsize',
                        action='store', default=30, type=int,
                        help='maximum size for internal oligo')
    parser.add_argument('--oligootm', dest='ep_otmopt',
                        action='store', default=69, type=int,
                        help='optimal Tm for internal oligo')
    parser.add_argument('--oligomintm', dest='ep_otmmin',
                        action='store', default=68, type=int,
                        help='minimum Tm for internal oligo')
    parser.add_argument('--oligomaxtm', dest='ep_otmmax',
                        action='store', default=70, type=int,
                        help='maximum Tm for internal oligo')
    parser.add_argument('--oligoogcpercent', dest='ep_ogcopt',
                        action='store', default=55, type=int,
                        help='optimal %%GC for internal oligo')
    parser.add_argument('--oligomingcpercent', dest='ep_ogcmin',
                        action='store', default=30, type=int,
                        help='minimum %%GC for internal oligo')
    parser.add_argument('--oligomaxgcpercent', dest='ep_ogcmax',
                        action='store', default=80, type=int,
                        help='maximum %%GC for internal oligo')
# Commented out until Biopython command line code catches up with new
# EMBOSS options
#    parser.add_argument("--oligomaxpolyx", dest="ep_opolymax",
#                                 action="store",
#                                 default=3, type=int,
#                                 help="maximum run of repeated nucleotides " +
#                                 "in internal primer")
    parser.set_defaults(func=subcommands.subcmd_eprimer3)


def build_parser_blastscreen(subparsers, parents=None):
    """Add parser for `blastscreen` command to subparsers

    This parser controls options for screening predicted primers against a
    negative control database with BLASTN.
    """
    parser = subparsers.add_parser('blastscreen', aliases=['bs'],
                                   parents=parents)
    # BLASTN screen options - subcommand blastscreen
    parser.add_argument('outfilename',
                        help='Path to write new configuration file')
    parser.add_argument('--blastn', dest='bs_exe',
                        action='store', default='blastn',
                        help='path to BLASTN+ executable')
    parser.add_argument('--db', dest='bs_db',
                        action='store', default=None,
                        help='path to BLASTN+ database')
    parser.add_argument('--maxaln', dest='maxaln',
                        action='store', default=15, type=int,
                        help='exclude primers with longer alignment')
    parser.add_argument('--outdir', dest='bs_dir',
                        action='store', default='blastn',
                        help='path to directory for BLASTN+ output')
    parser.add_argument('-f', '--force', dest='bs_force',
                        action='store_true', default=False,
                        help='Overwrite old BLASTN+ output')
    parser.set_defaults(func=subcommands.subcmd_blastscreen)


def build_parser_primersearch(subparsers, parents=None):
    """Add parser for `primersearch` command to subparsers

    This parser controls options for in silico hybridisation of primers against
    positive and negative examples using the EMBOSS PrimerSearch tool.
    """
    parser = subparsers.add_parser('primersearch', aliases=['ps'],
                                   parents=parents)
    # primersearch options - subcommand primersearch
    parser.add_argument('outfilename',
                        help='Path to write new configuration file')
    parser.add_argument('--primersearch', dest='ps_exe',
                        action='store', default='primersearch',
                        help='path to primersearch executable')
    parser.add_argument('--outdir', dest='ps_dir',
                        action='store', default='primersearch',
                        help='path to directory for primersearch output')
    parser.add_argument('-f', '--force', dest='ps_force',
                        action='store_true', default=False,
                        help='Overwrite old primersearch output')
    parser.add_argument('--mismatchpercent', '-m', dest='mismatchpercent',
                        action='store', type=float, default=0.1,
                        help='Allowed percentage primer mismatch')
    parser.set_defaults(func=subcommands.subcmd_primersearch)


def build_parser_classify(subparsers, parents=None):
    """Add parser for `classify` command to subparsers

    This parser controls options for classifying predicted primer sets.
    """
    parser = subparsers.add_parser('classify', aliases=['cl'],
                                   parents=parents)
    parser.set_defaults(func=subcommands.subcmd_classify)


# Process command-line
def parse_cmdline(args=None):
    """Parse command-line arguments for script.

    The script offers a single main parser, with subcommands for the actions:

    process - process/check input data (stitch fragments/fix sequence problems)
    prodigal - run Prodigal to predict CDS features, and update .conf file
    eprimer3 - run ePrimer3 to predict primer pairs, and update .conf file
    check_blast - check/filter designed primers against negative example
                  BLAST database
    primersearch - check/filter designed primers against complete genome
                   negative examples
    classify - classify designed primers against input genome/classes
    """
    # Main parent parser
    parser_main = ArgumentParser(prog='pdp.py')
    subparsers = parser_main.add_subparsers(title='subcommands',
                                            description='valid subcommands',
                                            help='additional help')

    # Common parser to be included with all the subcommand parsers
    parser_common = build_common_parser()
    parser_scheduler = build_scheduler_parser()

    # Add subcommand parsers to the main parser's subparsers
    build_parser_config(subparsers, parents=[parser_common])
    build_parser_prodigal(subparsers, parents=[
                          parser_common, parser_scheduler])
    build_parser_eprimer3(subparsers, parents=[
                          parser_common, parser_scheduler])
    build_parser_blastscreen(subparsers, parents=[
                             parser_common, parser_scheduler])
    build_parser_primersearch(
        subparsers, parents=[parser_common, parser_scheduler])
    build_parser_classify(subparsers, parents=[parser_common])

    # Parse arguments
    if args is None:
        args = sys.argv[1:]
    return parser_main.parse_args(args)
