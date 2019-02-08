#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_config.py

Provides the config subcommand for pdp.py

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

import os

from tqdm import tqdm

from diagnostic_primers.scripts.tools import load_config_tab, load_config_json


def ensure_path_to(fname):
    """If the path to the passed filename doesn't exist, it is created."""
    absfname = os.path.abspath(fname)
    pathto = os.path.dirname(absfname)
    os.makedirs(pathto, exist_ok=True)
    return absfname


def subcmd_config(args, logger):
    """Run `config` subcommand operations.

    Available subcommands:

    - to_json: convert .tab format config file to JSON and write out
    - to_tab: convert JSON format config file to .tab and write out
    - validate: report on the status of sequence files/classes
    - fix_sequences: stitch multiple sequences together, convert ambiguity
                     symbols to Ns, and write out a new JSON config file

    All subcommands should take either .tab or .json config files,
    distinguished by file extension
    """
    # Determine input config file type
    configtype = os.path.splitext(args.infilename)[-1][1:]
    if configtype not in ("tab", "json", "conf"):
        logger.error(
            "Expected config file to end in .conf, .json or .tab " + "got %s (exiting)",
            configtype,
        )
        raise SystemExit(1)

    if configtype in ("tab", "conf"):
        coll = load_config_tab(args, logger)
    elif configtype in ("json",):
        coll = load_config_json(args, logger)

    if args.outdir is not None:
        if not os.path.isdir(args.outdir):
            try:
                logger.info("Attempting to create output directory %s", args.outdir)
                os.makedirs(args.outdir)
            except Exception:
                logger.info(
                    "Could not create output directory %s (exiting)", args.outdir
                )
                raise SystemExit(1)

    # Do sequences need to be stitched or their ambiguities replaced?
    # If --validate is active, we report only and do not modify.
    # Note that the underlying code doesn't require us to check whether
    # the sequence files need stitching or symbols replacing, and
    # we could more efficiently use
    # for g in gc.data:
    #     g.stitch()
    #     g.replace_ambiguities()
    # but we're being helpfully verbose.
    logger.info(
        "Checking whether input sequences require stitching, "
        + "or have non-N ambiguities."
    )
    problems = ["Validation problems"]  # Holds messages about problem files
    pbar = tqdm(
        coll.data, desc="identifying validation problems", disable=args.disable_tqdm
    )
    for gcc in pbar:
        if gcc.needs_stitch:
            msg = "%s requires stitch" % gcc.name
            problems.append("%s (%s)" % (msg, gcc.seqfile))
            if args.fix_sequences:
                gcc.stitch(outdir=args.outdir)
        if gcc.has_ambiguities:
            msg = "%s has non-N ambiguities" % gcc.name
            problems.append("%s (%s)" % (msg, gcc.seqfile))
            if args.fix_sequences:
                gcc.replace_ambiguities(outdir=args.outdir)

    # If we were not fixing sequences, report problems
    if not args.fix_sequences:
        if len(problems) > 1:
            logger.warning("\n    ".join(problems))

    # Write post-processing config file and exit
    if args.to_json:
        logger.info("Writing JSON config file to %s", ensure_path_to(args.to_json))
        coll.write_json(args.to_json)
    elif args.to_tab:
        logger.info("Writing .tab file to %s", ensure_path_to(args.to_tab))
        coll.write_tab(args.to_tab)
    elif args.fix_sequences:
        logger.info(
            "Writing JSON config file to %s", ensure_path_to(args.fix_sequences)
        )
        coll.write_json(args.fix_sequences)
    return 0
