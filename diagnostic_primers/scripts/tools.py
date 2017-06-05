#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tools.py

Provides useful functions for the pdp.py script

- config:        interact with config files
- prodigal:      run prodigal bacterial gene prediction
- eprimer3:      predict primers with EMBOSS ePrimer3
- primersearch:  run in-silico hybridisation with EMBOSS PrimerSearch
- blastscreen:   screen primers with BLASTN
- classify:      classify primers

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

import sys
import traceback

from diagnostic_primers import (multiprocessing, sge, sge_jobs, config)


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Load PDPCollection from .tab file
def load_config_tab(args, logger):
    """Load tab format config to PDPCollection."""
    pdpc = config.PDPCollection()
    try:
        pdpc.from_tab(args.infilename)
    except:
        logger.error('Could not read config file %s (exiting)',
                     args.infilename)
        logger.error(last_exception())
        raise SystemExit(1)
    return pdpc


# Load PDPCollection from .json file
def load_config_json(args, logger):
    """Load JSON format config to PDPCollection."""
    pdpc = config.PDPCollection()
    try:
        pdpc.from_json(args.infilename)
    except:
        logger.error('Could not read config file %s (exiting)',
                     args.infilename)
        logger.error(last_exception())
        raise SystemExit(1)
    return pdpc


# Report a list of command lines to a logger, in pretty format
def log_clines(clines, logger):
    """Log command-lines, one per line."""
    logger.info('...%d commands returned:\n%s' %
                (len(clines), '\n'.join(['  %s' % c for c in clines])))


# Pass jobs to the appropriate scheduler
def run_parallel_jobs(clines, args, logger):
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
            raise SystemExit(1)
        else:
            logger.info('Runs completed without error.')
    elif args.scheduler == 'SGE':
        joblist = [sge_jobs.Job("pdp_%06d" % idx, cmd) for idx, cmd in
                   enumerate(clines)]
        sge.run_dependency_graph(joblist, verbose=True, logger=logger)
    else:
        raise ValueError('Scheduler must be one of ' +
                         '[multiprocessing|SGE], got %s' % args.scheduler)
