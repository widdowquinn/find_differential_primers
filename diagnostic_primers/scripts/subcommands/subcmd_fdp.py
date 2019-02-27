#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_fdp.py

Provides find_differential_primers.py back-compatibility

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

from argparse import Namespace

from diagnostic_primers.scripts import subcommands
from diagnostic_primers.scripts.logger import build_logger


def subcmd_fdp(args, logger):
    """Provides find_differential_primers.py script back-compatibility

    :param args:  Namespace with command-line arguments
    :param logger:  logging.logger object

    This function aims to emulate, using the current diagnostic_primers package
    and modules, the functionality of the old find_differential_primers script.

    The general approach is to create a new Namespace for each pdp subcommand,
    using arguments taken from the passed Namespace in args, and call the
    corresponding subcommand.

    We use the passed logger object as the master logger, and also generate
    a new logger for each step in the pdp pipeline.
    """
    # Get the logger filestem
    logfilestem = os.path.splitext(args.logfile)[0]
    logger.info("Using parent filestem for log files: %s", logfilestem)

    # Generate namespace and logger for pdp config step
    logger.info("Preparing pdp config step")
    config_namespace = Namespace(
        outdir=args.fdp_outdir,
        verbose=args.verbose,
        disable_tqdm=args.fdp_disable_tqdm,
        validate=False,
        fix_sequences=os.path.join(args.fdp_outdir, "fixed_config.json"),
        infilename=args.fdp_infilename,
        logfile=logfilestem + "_config.log",
        to_json=False,
        to_tab=False,
    )
    logger.info(
        "pdp config namespace:\n\t%s",
        "\n\t".join(
            ["{}: {}".format(_[0], _[1]) for _ in config_namespace.__dict__.items()]
        ),
    )
    config_logger = build_logger(
        "CONFIG (find_differential_primers.py)", config_namespace
    )
    logger.info("Running pdp config step")
    subcommands.subcmd_config(config_namespace, config_logger)
