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

import os

from argparse import ArgumentParser


def build_io_parser():
    """Return file IO command-line parser for find_differential_primers.py"""
    io_parser = ArgumentParser(add_help=False)
    io_parser.add_argument(
        "-i",
        "--infile",
        dest="infilename",
        action="store",
        help="path to configuration file",
        default=None,
    )
    io_parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        action="store",
        help="path to directory for output files",
        default="differential_primer_results",
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
        dest="debug",
        help="deprecated: does nothing",
        default=False,
    )
    io_parser.add_argument(
        "--disable_tqdm",
        action="store_true",
        dest="disable_tqdm",
        default=False,
        help="turn off tqdm progress bar",
    )
    return io_parser


def build_log_parser():
    """Return logging command-line parser for find_differential_primers.py"""
    log_parser = ArgumentParser(add_help=False)
    log_parser.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        action="store",
        default="fdp",
        help="filename stem for recording logfiles",
    )
    log_parser.add_argument(
        "--log_dir",
        dest="log_dir",
        action="store",
        help="path to common directory for log files",
        default=os.curdir,
    )
    log_parser.add_argument(
        "--keep_logs",
        action="store_true",
        dest="keep_logs",
        help="deprecated: does nothing",
        default=False,
    )
    return log_parser
