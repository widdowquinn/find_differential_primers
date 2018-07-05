#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_prodigal.py

Provides the prodigal subcommand for pdp.py

(c) The James Hutton Institute 2017-18

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

Copyright (c) 2017-18 The James Hutton Institute
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

from diagnostic_primers import prodigal
from tqdm import tqdm

from ..tools import (create_output_directory, load_config_json, log_clines,
                     run_parallel_jobs)


def subcmd_prodigal(args, logger):
    """Run Prodigal to predict bacterial CDS on the input sequences."""
    # Determine config input type
    configtype = os.path.splitext(args.infilename)[-1][1:]
    if configtype not in ('tab', 'json', 'conf'):
        logger.error(
            "Expected config file to end in .conf, .json or .tab " +
            "got %s (exiting)", configtype)
        raise SystemExit(1)

    if configtype in ('tab', 'conf'):
        raise ValueError("prodigal subcommand requires JSON config file")
    coll = load_config_json(args, logger)

    # Check if output exists and if we should overwrite
    create_output_directory(args.prodigaldir, args.prodigalforce, logger)

    # Build command-lines for Prodigal and run
    logger.info('Building Prodigal command lines...')
    clines = prodigal.build_commands(coll, args.prodigal_exe, args.prodigaldir)
    log_clines(clines, logger)
    run_parallel_jobs(clines, args, logger)

    # Add Prodigal output files to the GenomeData objects and write
    # the config file
    for gcc in tqdm(coll.data):
        gcc.features = gcc.cmds['prodigal'].split()[-1].strip()
        logger.info('%s feature file:\t%s' % (gcc.name, gcc.features))
    logger.info('Writing new config file to %s' % args.outfilename)
    coll.write_json(args.outfilename)
    return 0
