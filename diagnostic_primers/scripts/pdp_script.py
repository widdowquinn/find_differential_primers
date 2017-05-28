#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pdp_script.py

Implements the pdp script for finding differential primers

(c) The James Hutton Institute 2017

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

import logging
import logging.handlers
import os
import shutil
import sys
import time
import traceback

from . import parsers, tools


def run_parallel_jobs(clines):
    """Run the passed command-lines in parallel."""
    logger.info('Running jobs using scheduler: %s' % args.scheduler)
    # Pass lines to scheduler and run
    if args.scheduler == 'multiprocessing':
        retvals = multiprocessing.run(clines, workers=args.workers,
                                      verbose=args.verbose)
        if sum([r.returncode for r in retvals]):
            logger.error('At least one run has problems (exiting).')
            for retval in retvals:
                if retval.returncode != 0:
                    logger.error('Failing command: %s' % retval.args)
                    logger.error('Failing stderr:\n %s' % retval.stderr)
            sys.exit(1)
        else:
            logger.info('Runs completed without error.')
    elif args.scheduler == 'SGE':
        joblist = [sge_jobs.Job("pdp_%06d" % idx, cmd) for idx, cmd in
                   enumerate(clines)]
        sge.run_dependency_graph(joblist, verbose=True, logger=logger)
    else:
        raise ValueError('Scheduler must be one of ' +
                         '[multiprocessing|SGE], got %s' % args.scheduler)


def log_clines(clines):
    """Log command-lines, one per line."""
    logger.info('...%d commands returned:\n%s' %
                (len(clines), '\n'.join(['\t%s' % c for c in clines])))


def run_pdp_main(namespace=None):
    """Main process for pdp.py script"""
    # If we need to (a namespace isn't passed), parse the command-line
    if namespace is None:
        args = parsers.parse_cmdline()
    else:
        args = namespace

    # Set up logging
    logger = logging.getLogger('pdp.py: %s' % time.asctime())
    t0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error('Could not open %s for logging' %
                         args.logfile)
            logger.error(last_exception())
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info('Processed arguments: %s' % args)
    logger.info('command-line: %s' % ' '.join(sys.argv))

    # Run the subcommand
    return args.func(args, logger)
