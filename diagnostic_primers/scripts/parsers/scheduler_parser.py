# -*- coding: utf-8 -*-
"""Parser for common pdp scheduler options

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


# Build common parser for commands that use a scheduler
def build():
    """Returns a parser with options for a scheduler.

    This parser implements options that are common to commands that use
    a scheduler.
    """
    parser_scheduler = ArgumentParser(add_help=False)
    parser_scheduler.add_argument(
        "-s",
        "--scheduler",
        dest="scheduler",
        action="store",
        default="multiprocessing",
        type=str,
        help="Job scheduler [multiprocessing|SGE]",
    )
    parser_scheduler.add_argument(
        "-w",
        "--workers",
        dest="workers",
        action="store",
        default=None,
        type=int,
        help="Number of parallel workers to use",
    )
    parser_scheduler.add_argument(
        "--SGEgroupsize",
        dest="sgegroupsize",
        action="store",
        default=10000,
        type=int,
        help="Number of jobs to place in an SGE array group " "(default 10000)",
    )
    parser_scheduler.add_argument(
        "--SGEargs",
        dest="sgeargs",
        action="store",
        default=None,
        type=str,
        help="Additional arguments for qsub",
    )
    parser_scheduler.add_argument(
        "--jobprefix",
        dest="jobprefix",
        action="store",
        default="pdp",
        type=str,
        help="prefix for scheduled jobs",
    )
    parser_scheduler.add_argument(
        "--recovery",
        dest="recovery",
        action="store_true",
        default=False,
        help="skip scheduled third-party tool calls and reuse output",
    )
    return parser_scheduler
