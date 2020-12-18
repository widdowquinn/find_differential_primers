# -*- coding: utf-8 -*-
"""Parser for common pdp options

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
from argparse import ArgumentDefaultsHelpFormatter


# Build common parser options for all subcommands
def build():
    """Return a common argument parser for pdp.

    This parser implements options that are common to all subcommands.
    """
    parser_common = ArgumentParser(add_help=False,
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    parser_common.add_argument("infilename", help="path to configuration file")
    parser_common.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        action="store",
        default=None,
        help="logfile location",
    )
    parser_common.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="report progress to log",
    )
    parser_common.add_argument(
        "--disable_tqdm",
        action="store_true",
        dest="disable_tqdm",
        default=False,
        help="turn off tqdm progress bar",
    )
    return parser_common
