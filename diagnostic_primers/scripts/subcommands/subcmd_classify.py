#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_classify.py

Provides the classify subcommand for pdp.py

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

from tqdm import tqdm

from diagnostic_primers import classify

from ..tools import create_output_directory, has_primersearch, load_config_json


def subcmd_classify(args, logger):
    """Perform classification of predicted primers."""
    logger.info("Classifying primers for specificity")

    # Create output directory, if needed
    create_output_directory(args.outdir, args.cl_force, logger)

    # Load the JSON config file (post-primersearch)
    coll = load_config_json(args, logger)

    # Test whether the collection has primersearch output
    if not has_primersearch(coll):
        logger.error(
            " ".join(
                [
                    "To use the classify subcommand, the JSON file",
                    "must contain links to primersearch data.",
                    "(exiting)",
                ]
            )
        )
        raise SystemExit(1)
    logger.info("All input genomes have linked path to PrimerSearch data:")
    pbar = tqdm(coll.data, disable=args.disable_tqdm)
    for genome in pbar:
        pbar.set_description("%s -> %s" % (genome.name, genome.primersearch))

    # Obtain classification of all primer sets linked from config file, and
    # report to logger
    results = classify.classify_primers(coll)
    logger.info(
        "Identified primers specific to groups:\n\t%s", "\n\t".join(results.groups)
    )
    for group in results.groups:
        logger.info(
            "Primers specific to %s:\n\t%s",
            group,
            "\n\t".join([primer.name for primer in results.diagnostic_primer(group)]),
        )

    # Write diagnostic primer outputs to the output directory
    classify.write_results(results, os.path.join(args.outdir, "results.json"))
    classify.write_results(
        results, os.path.join(args.outdir, "summary.tab"), fmt="summary"
    )

    return 0
