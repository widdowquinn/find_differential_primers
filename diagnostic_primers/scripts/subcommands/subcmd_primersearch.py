#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_primersearch.py

Provides the primersearch subcommand for pdp.py

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

from diagnostic_primers import primersearch

from ..tools import (create_output_directory, load_config_json, log_clines,
                     run_parallel_jobs)


def subcmd_primersearch(args, logger):
    """Perform in silico hybridisation with EMBOSS PrimerSearch."""
    # Does output already exist, and should we overwrite?
    create_output_directory(args.ps_dir, args.ps_force, logger)

    # Get config file data
    coll = load_config_json(args, logger)

    # Construct command lines for primersearch
    logger.info("Building primersearch command-lines...")
    mismatchpercent = int(100 * args.mismatchpercent)  # for EMBOSS
    clines = primersearch.build_commands(coll, args.ps_exe, args.ps_dir,
                                         mismatchpercent)
    pretty_clines = [str(c).replace(' -', ' \\\n          -') for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs(clines, args, logger)

    # Write new config file, and exit
    logger.info('Writing new config file to %s', args.outfilename)
    coll.write_json(args.outfilename)
    return 0
