# -*- coding: utf-8 -*-
"""Parser for blastscreen subcommand

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

from argparse import ArgumentDefaultsHelpFormatter
from diagnostic_primers.scripts import subcommands


def build(subparsers, parents=None):
    """Add parser for `blastscreen` command to subparsers

    subparsers           ArgumentParser.subparser
    parents              additional parser objects

    This parser controls options for screening predicted primers against a
    negative control database with BLASTN.
    """
    parser = subparsers.add_parser("blastscreen", aliases=["bs"], parents=parents,
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    # BLASTN screen options - subcommand blastscreen
    parser.add_argument("outfilename", help="Path to write new configuration file")
    parser.add_argument(
        "--blastn",
        dest="bs_exe",
        action="store",
        default="blastn",
        help="path to BLASTN+ executable",
    )
    parser.add_argument(
        "--db",
        dest="bs_db",
        action="store",
        default="nr",
        help="path to BLASTN+ database",
    )
    parser.add_argument(
        "--maxaln",
        dest="maxaln",
        action="store",
        default=15,
        type=int,
        help="exclude primers with longer alignment",
    )
    parser.add_argument(
        "--outdir",
        dest="bs_dir",
        action="store",
        default="blastn",
        help="path to directory for BLASTN+ output",
    )
    parser.add_argument(
        "--jsondir",
        dest="bs_jsondir",
        action="store",
        default=None,
        help="path to directory for screened JSON output",
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="bs_force",
        action="store_true",
        default=False,
        help="Overwrite old BLASTN+ output",
    )
    parser.set_defaults(func=subcommands.subcmd_blastscreen)
