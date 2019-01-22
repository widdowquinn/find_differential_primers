# -*- coding: utf-8 -*-
"""Parser for classify subcommand

(c) The James Hutton Institute 2017-2019

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


def build(subparsers, parents=None):
    """Add parser for `classify` command to subparsers

    subparsers           ArgumentParser.subparser
    parents              additional parser objects

    This parser controls options for classifying predicted primer sets.
    """
    parser = subparsers.add_parser("classify", aliases=["cl"], parents=parents)
    parser.add_argument("outdir", help="Path to directory for output")
    parser.add_argument(
        "-f",
        "--force",
        dest="cl_force",
        action="store_true",
        default=False,
        help="Overwrite old classifier output",
    )
    parser.add_argument(
        "--minamplicon",
        dest="cl_minamplicon",
        action="store",
        type=int,
        default=50,
        help="Smallest amplicon size to accept as cross-hybridisation",
    )
    parser.add_argument(
        "--maxamplicon",
        dest="cl_maxamplicon",
        action="store",
        type=int,
        default=300,
        help="Longest amplicon size to accept as cross-hybridisation",
    )
    parser.set_defaults(func=subcommands.subcmd_classify)
