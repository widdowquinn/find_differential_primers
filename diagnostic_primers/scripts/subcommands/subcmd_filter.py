#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_filter.py

Provides the filter subcommand for pdp.py

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

from diagnostic_primers import (
    prodigal, )
from tqdm import tqdm

from ..tools import (create_output_directory, load_config_json, log_clines,
                     run_parallel_jobs,)


def subcmd_filter(args, logger):
    """Run filter on input sequences, add filter and filtered_seqfile to config file."""
    # Exactly one of the following filter modes must be selected
    filtermodes = [args.filt_prodigal, args.filt_prodigaligr]
    if sum(filtermodes) != 1:
        logger.error("Invalid filter options chosen (exiting)")
        if sum(filtermodes) == 0:
            logger.error("No filter modes chosen.")
        if sum(filtermodes) > 0:
            logger.error("More than one filter mode chosen")
        raise SystemExit(1)

    # Determine config input type
    configtype = os.path.splitext(args.infilename)[-1][1:]
    if configtype not in ('tab', 'json', 'conf'):
        logger.error(
            "Expected config file to end in .conf, .json or .tab " +
            "got %s (exiting)", configtype)
        raise SystemExit(1)

    if configtype in ('tab', 'conf'):
        raise ValueError(
            "filter subcommand requires JSON config file, please convert the input"
        )

    # Load config file
    logger.info("Reading data from %s", args.infilename)
    coll = load_config_json(args, logger)

    # Check if output exists and if we should overwrite
    create_output_directory(args.filt_outdir, args.filt_force, logger)

    # Both the --prodigal and --prodigaligr modes require prodigal genecaller output
    # Build command-lines for Prodigal and run
    if args.filt_prodigal or args.filt_prodigaligr:
        logger.info('Building Prodigal command lines...')
        clines = prodigal.build_commands(coll, args.filt_prodigal_exe,
                                         args.filt_outdir)
        log_clines(clines, logger)
        run_parallel_jobs(clines, args, logger)

        # If the filter is prodigal, add Prodigal output files to the GenomeData
        # objects and write the config file; if the filter is prodigaligr, then
        # identify the intergenic loci, and create a new GFF feature file of
        # intergenic regions, with a buffer sequence into the flanking genes
        if args.filt_prodigal:
            logger.info("Collecting Prodigal prediction output")
            for gcc in tqdm(coll.data):
                gcc.features = gcc.cmds['prodigal'].split()[-1].strip()
                logger.info('%s feature file:\t%s', gcc.name, gcc.features)
            logger.info('Writing new config file to %s', args.outfilename)
        elif args.filt_prodigaligr:
            logger.info("Calculating intergenic regions from Prodigal output")
            for gcc in tqdm(coll.data):
                prodigalout = gcc.cmds['prodigal'].split()[-1].strip()
                bedpath = os.path.splitext(prodigalout)[0] + '_igr.gff'
                logger.info("Prodigal input: %s; BED output: %s", prodigalout, bedpath)
                prodigal.generate_igr(prodigalout, gcc.seqfile, bedpath)
                gcc.features = bedpath
                logger.info('%s feature file:\t%s', gcc.name, gcc.features)

    # Compile a new genome from the gcc.features file
    logger.info(
        "Compiling filtered sequences from primer design target features")
    for gcc in tqdm(coll.data):
        logger.info("Reading primer design feature data from %s", gcc.features)
        stem, ext = os.path.splitext(gcc.seqfile)
        filtered_path = stem + '_' + args.filt_suffix + ext
        logger.info("Writing filtered sequence for primer design to %s",
                    filtered_path)
        gcc.create_filtered_genome(filtered_path, args.filt_spacerlen,
                                   args.filt_suffix)

    # Write updated config file
    logger.info("Writing new config file to %s", args.outfilename)
    coll.write_json(args.outfilename)

    return 0
