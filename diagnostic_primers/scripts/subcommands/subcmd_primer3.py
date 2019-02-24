#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_primer3.py

Provides the primer3 subcommand for pdp

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

Copyright (c) 2019 The James Hutton Institute

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

from tqdm import tqdm

from diagnostic_primers import primer3, load_primers, write_primers
from diagnostic_primers.scripts.tools import (
    create_output_directory,
    load_config_json,
    log_clines,
    run_parallel_jobs,
)


def subcmd_primer3(args, logger):
    """Run Primer3 (v2+) to design primers for each input sequence."""
    # Determine config input type
    configtype = os.path.splitext(args.infilename)[-1][1:]
    if configtype not in ("tab", "json", "conf"):
        logger.error(
            "Expected config file to end in .conf, .json or .tab " + "got %s (exiting)",
            configtype,
        )
        raise SystemExit(1)

    if configtype in ("tab", "conf"):
        raise ValueError("primer3 subcommand requires JSON config file")
    coll = load_config_json(args, logger)

    # Check if output exists and if we should overwrite
    create_output_directory(args.primer3_dir, args.primer3_force, logger)

    # Build command-lines for primer3 and run
    # This will write 'bare' primer3 files, with unnamed primer pairs
    logger.info("Building primer3 command lines...")
    clines = primer3.build_commands(
        coll, args.primer3_exe, args.primer3_dir, vars(args)
    )
    logger.info(
        "Created input files for Primer3 (v2+):\n\t",
        "\n\t".join([_.infile for _ in clines]),
    )
    pretty_clines = [str(c).replace(" -", " \\\n          -") for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs([" ".join(_.cline) for _ in clines], args, logger)

    # Parse Primer3 output and write out as JSON
    pbar = tqdm(coll.data, desc="writing primer sets", disable=args.disable_tqdm)
    for gcc in pbar:
        p3file = gcc.cmds["Primer3"].outfile
        primers = load_primers(p3file, fmt="primer3")
        # Add source genome to each primer
        for primer in primers:
            primer.source = gcc.seqfile
            primer.sourcename = gcc.name
        # Write named ePrimer3
        outfname = os.path.splitext(p3file)[0] + "_named.eprimer3"
        write_primers(primers, outfname, fmt="ep3")
        # Write named BED
        outfname = os.path.splitext(p3file)[0] + "_named.bed"
        pbar.set_description("Writing: %s" % outfname)
        write_primers(primers, outfname, fmt="bed")
        # Write named JSON (the reference description)
        outfname = os.path.splitext(p3file)[0] + "_named.json"
        write_primers(primers, outfname, fmt="json")
        gcc.primers = outfname

    logger.info("Writing new config file to %s" % args.outfilename)
    coll.write_json(args.outfilename)
    return 0
