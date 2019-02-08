# -*- coding: utf-8 -*-
"""Parser for pdp dedupe subcommand

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


def build(subparsers, parents=None):
    """Add parser for `ddedupe` command to subparsers

    This parser controls options for deduplicating primers in an input config
    JSON file.
    """
    parser = subparsers.add_parser("dedupe", aliases=["dd"], parents=parents)
    # dedupe options - subcommand dedupe
    parser.add_argument("outfilename", help="Path to write new configuration file")
    parser.add_argument(
        "-f",
        "--force",
        dest="dd_force",
        action="store_true",
        default=False,
        help="Overwrite old config file target",
    )
    parser.add_argument(
        "--dedupedir",
        dest="dd_dedupedir",
        action="store",
        default=None,
        help="Output directory for deduplicated primer JSON files",
    )
    parser.set_defaults(func=subcommands.subcmd_dedupe)
