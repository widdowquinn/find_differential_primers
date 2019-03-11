#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_primersearch.py

Provides the primersearch subcommand for pdp.py

(c) The James Hutton Institute 2017-19

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

Copyright (c) 2017-19 The James Hutton Institute
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

from diagnostic_primers import primersearch
from diagnostic_primers.scripts.tools import (
    collect_existing_output,
    create_output_directory,
    load_config_json,
    log_clines,
    run_parallel_jobs,
)


def subcmd_primersearch(args, logger):
    """Perform in silico hybridisation with EMBOSS PrimerSearch."""
    # Does output already exist, and should we overwrite?
    create_output_directory(args.ps_dir, args.ps_force, logger)

    # Get config file data
    coll = load_config_json(args, logger)

    # If we are in recovery mode, we are salvaging output from a previous
    # run, and do not necessarily need to rerun all the jobs. In this case,
    # we prepare a list of output files we want to recover from the results
    # in the output directory.
    existingfiles = []
    if args.recovery:
        logger.warning("Entering recovery mode")
        logger.info(
            "\tIn this mode, existing comparison output from %s is reused", args.ps_dir
        )
        existingfiles = collect_existing_output(args.ps_dir, "primersearch", args)
        logger.info(
            "Existing files found:\n\t%s", "\n\t".join([_ for _ in existingfiles])
        )

    # Construct command lines for primersearch
    logger.info("Building primersearch command-lines...")
    mismatchpercent = int(100 * args.mismatchpercent)  # for EMBOSS
    clines = primersearch.build_commands(
        coll, args.ps_exe, args.ps_dir, mismatchpercent, existingfiles
    )
    if len(clines):
        pretty_clines = [str(c).replace(" -", " \\\n          -") for c in clines]
        log_clines(pretty_clines, logger)
        run_parallel_jobs(clines, args, logger)
    else:
        logger.warning(
            "No primersearch jobs were scheduled (you may see this if the --recovery option is active)"
        )

    # Load PrimerSearch output and generate .json/.bed files of amplimers
    # (regions on each target genome amplified by a primer)
    logger.info("Identifying target amplicoms")
    amplimers = primersearch.load_collection_amplicons(coll)
    amplimerpath = os.path.join(args.ps_dir, "target_amplicons.json")
    logger.info("Writing all target amplicons to %s", amplimerpath)
    amplimers.write_json(amplimerpath)
    # Subdivide the amplimers into a new PDPGenomeAmplicons object - one per
    # input genome, and write bed/JSON files accordingly
    logger.info("Writing individual amplicon files for each target")
    for obj in amplimers.split_on_targets():
        jsonpath = os.path.join(args.ps_dir, "{}_amplicons.json".format(obj.targets[0]))
        logger.info("\tWorking with target %s", obj.targets[0])
        obj.write_json(jsonpath)
        obj.write_bed(args.ps_dir)
        # Add the JSON file to the appropriate entry in the collection
        coll[obj.name].target_amplicons = jsonpath

    # Write new config file, and exit
    logger.info("Writing new config file to %s", args.outfilename)
    coll.write_json(args.outfilename)
    return 0
