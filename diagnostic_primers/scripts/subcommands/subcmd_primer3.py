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

from diagnostic_primers import primer3
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
    pretty_clines = [str(c).replace(" -", " \\\n          -") for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs(clines, args, logger)
