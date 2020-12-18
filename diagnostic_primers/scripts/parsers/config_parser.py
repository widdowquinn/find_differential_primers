# -*- coding: utf-8 -*-
"""Parser for pdp config subcommand

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
    """Add parser for `config` subcommand to subparsers

    This parser implements options for processing configuration files.
    """
    parser = subparsers.add_parser("config", parents=parents,
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "--outdir",
        action="store",
        dest="outdir",
        default=None,
        help="Path to directory to place output files",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        dest="validate",
        default=False,
        help="Validate config file, then exit",
    )
    parser.add_argument(
        "--fix_sequences",
        action="store",
        dest="fix_sequences",
        default=None,
        help="Fix config file sequences and write new JSON " + "config file",
    )
    parser.add_argument(
        "--to_json",
        action="store",
        dest="to_json",
        default=None,
        help="Convert .tab config file to JSON and write",
    )
    parser.add_argument(
        "--to_tab",
        action="store",
        dest="to_tab",
        default=None,
        help="Convert JSON config file to .tab and write",
    )
    parser.set_defaults(func=subcommands.subcmd_config)
