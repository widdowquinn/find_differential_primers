# -*- coding: utf-8 -*-
"""Parser for pdp extract subcommand

(c) The James Hutton Institute 2017-2019

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

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

from diagnostic_primers.scripts import subcommands
from argparse import ArgumentDefaultsHelpFormatter

def build(subparsers, parents=None):
    """Add parser for `extract` command to subparsers

    This parser controls options for extracting amplicons and other
    information, given a set of primers and corresponding genomes/
    primersearch output.

    The parser requires input:

    - the JSON file describing the primers
    - a JSON config file including primersearch output
    """
    parser = subparsers.add_parser("extract", aliases=["ex"], parents=parents,
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("primerfile", help="Path to the JSON format primer file")
    parser.add_argument("outdir", help="Path to directory for output")
    parser.add_argument(
        "-f",
        "--force",
        dest="ex_force",
        action="store_true",
        default=False,
        help="Overwrite old extract output",
    )
    parser.add_argument(
        "--mafft",
        dest="mafft_exe",
        action="store",
        default="mafft",
        help="Path to MAFFT executable",
    )
    parser.add_argument(
        "--noalign",
        dest="noalign",
        action="store_true",
        default=False,
        help="Suppress amplicon alignment",
    )
    parser.add_argument(
        "--minamplicon",
        dest="ex_minamplicon",
        action="store",
        type=int,
        default=50,
        help="Smallest amplicon size to process",
    )
    parser.add_argument(
        "--maxamplicon",
        dest="ex_maxamplicon",
        action="store",
        type=int,
        default=300,
        help="Longest amplicon size to process",
    )
    parser.set_defaults(func=subcommands.subcmd_extract)
