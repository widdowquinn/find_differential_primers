#!/usr/bin/env/python3
# -*- coding: utf-8 -*-
"""find_differential_primers.py

Reproduces the behaviour/interface of the old find_differential_primers.py script

(c) The James Hutton Institute 2019

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

import os
import sys
import time

from diagnostic_primers import __version__
from diagnostic_primers.scripts import parsers
from diagnostic_primers.scripts import subcommands
from diagnostic_primers.scripts.logger import build_logger


def run_fdp_main(argv=None, logger=None):
    """Main process for the find_differential_primers.py script

    :param argv:  arguments for the script (Namespace preferred)
    :param logger:  Python logging object

    This entry point attempts to replicate the soon to be deprecated
    find_differential_primers.py interface using the updated pdp package.
    """
    # If no namespace is passed as argv, parse the command-line arguments
    if argv is None:
        args = parsers.parse_fdp()
    else:
        args = parsers.parse_fdp(argv)

    # Catch execution with no supplied arguments
    if len(sys.argv) == 1:
        sys.stderr.write(
            "pdp (find_differential_primers.py) version: {0}".format(__version__)
        )

    # Set up logging
    time0 = time.time()
    # Set up actual path to logfile for FDP
    args.logfile = os.path.join(args.fdp_log_dir, args.fdp_logfile + ".log")
    if logger is None:
        logger = build_logger("find_differential_primers.py", args)

    # Call and return the script function
    returnval = subcommands.subcmd_fdp(args, logger)
    logger.info("Completed. Time taken: %.3f", (time.time() - time0))
    return returnval
