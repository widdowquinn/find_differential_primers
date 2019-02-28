# -*- coding: utf-8 -*-
"""Module providing parser definitions for find_differential_primers.py

(c) The James Hutton Institute 2019

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

Copyright (c) 2019 The James Hutton Institute
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

import multiprocessing
import os

from argparse import ArgumentParser


# Default output directory
OUTDIR = "differential_primer_results"


def build_io_parser():
    """Return file IO command-line parser for find_differential_primers.py"""
    io_parser = ArgumentParser(add_help=False)
    io_parser.add_argument(
        "-i",
        "--infile",
        dest="fdp_infilename",
        action="store",
        help="path to configuration file",
        default=None,
    )
    io_parser.add_argument(
        "-o",
        "--outdir",
        dest="fdp_outdir",
        action="store",
        help="path to directory for output files",
        default=OUTDIR,
    )
    io_parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="report verbose output",
        default=False,
    )
    io_parser.add_argument(
        "--debug",
        action="store_true",
        dest="fdp_debug",
        help="deprecated: does nothing",
        default=False,
    )
    io_parser.add_argument(
        "--disable_tqdm",
        action="store_true",
        dest="fdp_disable_tqdm",
        default=False,
        help="turn off tqdm progress bar",
    )
    io_parser.add_argument(
        "--clean",
        action="store_true",
        dest="fdp_clean",
        help="clean up old output files before running",
        default=False,
    )
    io_parser.add_argument(
        "--cleanonly",
        action="store_true",
        dest="fdp_cleanonly",
        help="clean up old output files and exit",
        default=False,
    )
    return io_parser


def build_log_parser():
    """Return logging command-line parser for find_differential_primers.py"""
    log_parser = ArgumentParser(add_help=False)
    log_parser.add_argument(
        "-l",
        "--logfile",
        dest="fdp_logfile",
        action="store",
        default="fdp",
        help="filestem for recording logfiles",
    )
    log_parser.add_argument(
        "--log_dir",
        dest="fdp_log_dir",
        action="store",
        help="path to common directory for log files",
        default=OUTDIR,
    )
    log_parser.add_argument(
        "--keep_logs",
        action="store_true",
        dest="fdp_keep_logs",
        help="deprecated: does nothing",
        default=False,
    )
    return log_parser


def build_prodigal_parser():
    """Return prodigal command-line parser for find_differential_primers.py"""
    prod_parser = ArgumentParser(add_help=False)
    prod_parser.add_argument(
        "--prodigal",
        dest="fdp_prodigal_exe",
        action="store",
        help="path to Prodigal executable",
        default="prodigal",
    )
    prod_parser.add_argument(
        "--nocds",
        dest="fdp_nocds",
        action="store_true",
        help="deprecated: does nothing",
        default=False,
    )
    prod_parser.add_argument(
        "--noprodigal",
        dest="fdp_noprodigal",
        action="store_true",
        help="skip Prodigal CDS prediction step",
        default=False,
    )
    return prod_parser


def build_eprimer3_parser():
    """Return ePrimer3 command-line parser for find_differential_primers.py"""
    ep3_parser = ArgumentParser(add_help=False)
    ep3_parser.add_argument(
        "--eprimer3",
        dest="fdp_eprimer3_exe",
        action="store",
        help="path to EMBOSS ePrimer3 executable",
        default="eprimer3",
    )
    ep3_parser.add_argument(
        "--numreturn",
        dest="fdp_numreturn",
        action="store",
        help="number of primers to be reported by ePrimer3",
        default=20,
        type=int,
    )
    ep3_parser.add_argument(
        "--hybridprobe",
        dest="fdp_hybridprobe",
        action="store_true",
        help="generate internal oligo as a hybridisation probe",
        default=False,
    )
    ep3_parser.add_argument(
        "--filtergc3prime",
        dest="fdp_filtergc3prime",
        action="store_true",
        help="allow no more than two GC at the 3` end of each primer",
        default=False,
    )
    ep3_parser.add_argument(
        "--osize",
        dest="fdp_osize",
        action="store",
        help="optimal size for primer oligo",
        default=20,
        type=int,
    )
    ep3_parser.add_argument(
        "--minsize",
        dest="fdp_minsize",
        action="store",
        help="minimum size for primer oligo",
        default=18,
        type=int,
    )
    ep3_parser.add_argument(
        "--maxsize",
        dest="fdp_maxsize",
        action="store",
        help="maximum size for primer oligo",
        default=22,
        type=int,
    )
    ep3_parser.add_argument(
        "--otm",
        dest="fdp_otm",
        action="store",
        help="optimal melting temperature for primer oligo",
        default=59,
        type=int,
    )
    ep3_parser.add_argument(
        "--mintm",
        dest="fdp_mintm",
        action="store",
        help="minimum melting temperature for primer oligo",
        default=58,
        type=int,
    )
    ep3_parser.add_argument(
        "--maxtm",
        dest="fdp_maxtm",
        action="store",
        help="maximum melting temperature for primer oligo",
        default=60,
        type=int,
    )
    ep3_parser.add_argument(
        "--ogcpercent",
        dest="fdp_ogcpercent",
        action="store",
        help="optimal percent GC for primer oligo",
        default=55,
        type=int,
    )
    ep3_parser.add_argument(
        "--mingc",
        dest="fdp_mingc",
        action="store",
        help="minimum percent GC for primer oligo",
        default=30,
        type=int,
    )
    ep3_parser.add_argument(
        "--maxgc",
        dest="fdp_maxgc",
        action="store",
        help="maximum percent GC for primer oligo",
        default=80,
        type=int,
    )
    ep3_parser.add_argument(
        "--psizeopt",
        dest="fdp_psizeopt",
        action="store",
        help="optimal size for amplified region",
        default=100,
        type=int,
    )
    ep3_parser.add_argument(
        "--psizemin",
        dest="fdp_psizemin",
        action="store",
        help="minimum size for amplified region",
        default=50,
        type=int,
    )
    ep3_parser.add_argument(
        "--psizemax",
        dest="fdp_psizemax",
        action="store",
        help="maximum size for amplified region",
        default=150,
        type=int,
    )
    ep3_parser.add_argument(
        "--maxpolyx",
        dest="fdp_maxpolyx",
        action="store",
        help="maximum run of repeated nucleotides in primer",
        default=3,
        type=int,
    )
    ep3_parser.add_argument(
        "--mismatchpercent",
        dest="fdp_mismatchpercent",
        action="store",
        help="allowed percentage mismatch in primersearch",
        default=10,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligoosize",
        dest="fdp_oligoosize",
        action="store",
        help="optimal size for internal oligo",
        default=20,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligominsize",
        dest="fdp_oligominsize",
        action="store",
        help="minimum size for internal oligo",
        default=13,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligomaxsize",
        dest="fdp_oligomaxsize",
        action="store",
        help="maximum size for internal oligo",
        default=30,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligootm",
        dest="fdp_oligootm",
        action="store",
        help="optimal melting temperature for internal oligo",
        default=69,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligomintm",
        dest="fdp_oligomintm",
        action="store",
        help="minimum melting temperature for internal oligo",
        default=68,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligomaxtm",
        dest="fdp_oligomaxtm",
        action="store",
        help="maximum melting temperature for internal oligo",
        default=70,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligoogcpercent",
        dest="fdp_oligoogcpercent",
        action="store",
        help="optimal percent GC for internal oligo",
        default=55,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligomingc",
        dest="fdp_oligomingc",
        action="store",
        help="minimum percent GC for internal oligo",
        default=30,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligomaxgc",
        dest="fdp_oligomaxgc",
        action="store",
        help="maximum percent GC for internal oligo",
        default=80,
        type=int,
    )
    ep3_parser.add_argument(
        "--oligomaxpolyx",
        dest="fdp_oligomaxpolyx",
        action="store",
        help="maximum run of repeated nt in internal oligo",
        default=3,
        type=int,
    )
    return ep3_parser


def build_blast_parser():
    """Return blast command-line parser for find_differential_primers.py"""
    blast_parser = ArgumentParser(add_help=False)
    blast_parser.add_argument(
        "--blast_exe",
        dest="fdp_blast_exe",
        action="store",
        help="path to BLASTN/ executable",
        default="blastn",
    )
    blast_parser.add_argument(
        "--blastdb",
        dest="fdp_blastdb",
        action="store",
        help="path to off-target BLAST database",
        default=None,
    )
    blast_parser.add_argument(
        "--useblast",
        dest="fdp_useblast",
        action="store_true",
        help="deprecated: does nothing",
        default=False,
    )
    return blast_parser


def build_primersearch_parser():
    """Return primersearch command-line parser for find_differential_primers.py"""
    ps_parser = ArgumentParser(add_help=False)
    ps_parser.add_argument(
        "--noprimersearch",
        dest="fdp_noprimersearch",
        action="store_true",
        help="do not carry out PrimerSearch step",
        default=False,
    )
    return ps_parser


def build_classify_parser():
    """Return classify command-line parser for find_differential_primers.py"""
    classify_parser = ArgumentParser(add_help=False)
    classify_parser.add_argument(
        "--noclassify",
        dest="fdp_noclassify",
        action="store_true",
        help="do not carry out primer classification step",
        default=False,
    )
    return classify_parser


def build_scheduling_parser():
    """Return scheduling command-line parser for find_differential_primers.py"""
    scheduling_parser = ArgumentParser(add_help=False)
    scheduling_parser.add_argument(
        "--cpus",
        dest="fdp_cpus",
        action="store",
        help="number of CPUs to use in multiprocessing",
        default=multiprocessing.cpu_count(),
        type=int,
    )
    scheduling_parser.add_argument(
        "--sge",
        dest="fdp_sge",
        action="store_true",
        help="use SGE job scheduler",
        default=False,
    )
    return scheduling_parser
