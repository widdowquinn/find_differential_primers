#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_blastscreen.py

Provides the blastscreen subcommand for pdp

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


from tqdm import tqdm

from diagnostic_primers import blast
from diagnostic_primers.scripts.tools import (
    collect_existing_output,
    create_output_directory,
    load_config_json,
    log_clines,
    run_parallel_jobs,
)


def subcmd_blastscreen(args, logger):
    """Screen primer sequences against a local BLAST database."""
    # We can't go on if there's no BLASTN+ database to check against
    if args.bs_db is None:
        logger.error(
            "No BLASTN+ database provided, cannot perform " + "screen (exiting)"
        )
        raise SystemExit(1)

    # Check if output exists and if we should overwrite
    create_output_directory(args.bs_dir, args.bs_force, logger)

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
            "\tIn this mode, existing comparison output from %s is reused", args.bs_dir
        )
        existingfiles = collect_existing_output(args.bs_dir, "blastscreen", args)
        logger.info(
            "Existing files found:\n\t%s", "\n\t".join([_ for _ in existingfiles])
        )

    # Run BLASTN search with primer sequences
    logger.info("Building BLASTN screen command-lines...")
    clines = blast.build_commands(
        coll, args.bs_exe, args.bs_db, args.bs_dir, existingfiles
    )
    if len(clines):
        pretty_clines = [str(c).replace(" -", " \\\n          -") for c in clines]
        log_clines(pretty_clines, logger)
        run_parallel_jobs(clines, args, logger)
        logger.info("BLASTN+ search complete")
    else:
        logger.warning(
            "No BLASTN+ jobs were scheduled (you may see this if the --recovery option is active)"
        )

    # Amend primer JSON files to remove screened primers
    for blastout, indata in tqdm(
        zip([cline.out for cline in clines], coll.data),
        desc="removing screened primers",
        disable=args.disable_tqdm,
    ):
        logger.info(
            "Amending primer file %s with results from %s", indata.primers, blastout
        )
        newprimers = blast.apply_screen(
            blastout, indata.primers, jsondir=args.bs_jsondir, maxaln=args.maxaln
        )
        logger.info("Screened primers placed in %s", newprimers)
        indata.primers = newprimers

    # Write new config file post-BLASTN screen
    logger.info("Writing new config file to %s", args.outfilename)
    coll.write_json(args.outfilename)

    return 0
