#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pdp_script.py

Implements the pdp script for finding differential primers

(c) The James Hutton Institute 2017-2018

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

Copyright (c) 2017 The James Hutton Institute
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

import sys
import time

from diagnostic_primers import __version__
from diagnostic_primers.scripts import parsers
from diagnostic_primers.scripts.logger import build_logger


def run_pdp_main(argv=None, logger=None):
    """Main process for pdp script"""
    # If we need to (i.e. a namespace isn't passed), parse the command-line
    if argv is None:
        args = parsers.parse_cmdline()
    else:
        args = parsers.parse_cmdline(argv)

    # Catch execution with no arguments
    if len(sys.argv) == 1:
        sys.stderr.write("pdp version: {0}\n".format(__version__))
        return 0

    # Set up logging
    time0 = time.time()
    if logger is None:
        logger = build_logger("pdp", args)

    # Run the subcommand
    returnval = args.func(args, logger)
    logger.info("Completed. Time taken: %.3f", (time.time() - time0))
    return returnval
