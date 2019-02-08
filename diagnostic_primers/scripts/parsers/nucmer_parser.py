# -*- coding: utf-8 -*-
"""Parser for nucmer pdp options

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

from argparse import ArgumentParser


# Build plugin parser for commands using nucmer
def build():
    """Returns a parser with options for running nucmer.

    The parser implements options common to commands that need to run
    the nucmer sequence alignment tool.
    """
    parser_nucmer = ArgumentParser(add_help=False)
    parser_nucmer.add_argument(
        "--nucmer_exe",
        dest="nucmer_exe",
        action="store",
        default="nucmer",
        type=str,
        help="path to nucmer executable",
    )
    parser_nucmer.add_argument(
        "--deltafilter_exe",
        dest="deltafilter_exe",
        action="store",
        default="delta-filter",
        type=str,
        help="path to nucmer executable",
    )
    parser_nucmer.add_argument(
        "--maxmatch",
        dest="maxmatch",
        action="store_true",
        default=False,
        help="path to nucmer executable",
    )
    return parser_nucmer
