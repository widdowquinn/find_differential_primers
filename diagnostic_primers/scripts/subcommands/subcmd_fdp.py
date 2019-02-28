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

from diagnostic_primers import PDPException
from diagnostic_primers.scripts import subcommands
from diagnostic_primers.scripts.logger import build_logger
from diagnostic_primers.scripts.tools import get_primer3_version


class PDPFDPException(PDPException):
    """General exception when running find_differential_primers.py script"""

    def __init__(self, msg="Error using find_differential_primers.py"):
        PDPException.__init__(self, msg)


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
    # Check if we have the right primer3 version (this script only runs with
    # EMBOSS and primer3 v1.1.4)
    p3_version = get_primer3_version()
    if p3_version[0] > 1:
        raise PDPFDPException(
            "Primer3 v1.1.4 required, version {} found".format(
                ".".join([str(_) for _ in p3_version])
            )
        )

    # Get the logger filestem
    args.logfilestem = os.path.splitext(args.logfile)[0]
    logger.info("Using parent filestem for log files: %s", args.logfilestem)

    # Set scheduler
    if args.fdp_sge:
        args.scheduler = "SGE"
    else:
        args.scheduler = "multiprocessing"

    # Run each step of the pipeline
    args.infilename = run_config(args, logger)  # pdp config
    if args.fdp_noprodigal:  # pdp filter
        logger.info("Skipping pdp filter (prodigal) step")
    else:
        args.infilename = run_filter(args, logger)
    args.infilename = run_eprimer3(args, logger)
    args.infilename = run_dedupe(args, logger)


def run_config(args, logger):
    """Run pdp config step with find_differential_primers.py script

    :param args:  Namespace with command-line arguments
    :param logger:  logging.logger object
    """
    # Generate namespace and logger for pdp config step
    logger.info("Preparing pdp config step")
    outconfpath = os.path.join(args.fdp_outdir, "fixed_config.json")
    config_namespace = Namespace(
        outdir=args.fdp_outdir,
        verbose=args.verbose,
        disable_tqdm=args.fdp_disable_tqdm,
        validate=False,
        fix_sequences=outconfpath,
        infilename=args.fdp_infilename,
        logfile=args.logfilestem + "_config.log",
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
    return outconfpath


def run_filter(args, logger):
    """Run pdp filter step with find_differential_primers.py script

    :param args:  Namespace with command-line arguments
    :param logger:  logging.logger object
    """
    # Generate namespace and logger for pdp filter step
    logger.info("Preparing pdp filter step")
    outconfpath = os.path.join(args.fdp_outdir, "filtered.json")
    filter_namespace = Namespace(
        infilename=args.infilename,
        outfilename=outconfpath,
        filt_prodigal=True,
        filt_prodigaligr=False,
        filt_alnvar=None,
        filt_outdir=args.fdp_outdir,
        filt_prodigal_exe=args.fdp_prodigal_exe,
        filt_force=True,
        filt_suffix="prodigal",
        filt_spacerlen=150,
        filt_flanklen=150,
        scheduler=args.scheduler,
        workers=args.fdp_cpus,
        verbose=args.verbose,
        disable_tqdm=args.fdp_disable_tqdm,
        logfile=args.logfilestem + "_filter.log",
    )
    logger.info(
        "pdp filter namespace:\n\t%s",
        "\n\t".join(
            ["{}: {}".format(_[0], _[1]) for _ in filter_namespace.__dict__.items()]
        ),
    )
    filter_logger = build_logger(
        "FILTER (find_differential_primers.py)", filter_namespace
    )
    logger.info("Running pdp filter step")
    subcommands.subcmd_filter(filter_namespace, filter_logger)
    return outconfpath


def run_eprimer3(args, logger):
    """Run pdp eprimer3 step with find_differential_primers.py script

    :param args:  Namespace with command-line arguments
    :param logger:  logging.logger object
    """
    #  Generate namespace and logger for pdp eprimer3 step
    logger.info("Preparing pdp eprimer3 step")
    outconfpath = os.path.join(args.fdp_outdir, "eprimer3.json")
    eprimer3_namespace = Namespace(
        infilename=args.infilename,
        outfilename=outconfpath,
        eprimer3_dir=args.fdp_outdir,
        eprimer3_exe=args.fdp_eprimer3_exe,
        eprimer3_force=True,
        scheduler=args.scheduler,
        workers=args.fdp_cpus,
        verbose=args.verbose,
        ep_hybridprobe=args.fdp_hybridprobe,
        ep_filter=args.fdp_filtergc3prime,
        disable_tqdm=args.fdp_disable_tqdm,
        ep_numreturn=args.fdp_numreturn,
        ep_osize=args.fdp_osize,
        ep_minsize=args.fdp_minsize,
        ep_maxsize=args.fdp_maxsize,
        ep_opttm=args.fdp_otm,
        ep_mintm=args.fdp_mintm,
        ep_maxtm=args.fdp_maxtm,
        ep_ogcpercent=args.fdp_ogcpercent,
        ep_mingc=args.fdp_mingc,
        ep_maxgc=args.fdp_maxgc,
        ep_psizeopt=args.fdp_psizeopt,
        ep_psizemin=args.fdp_psizemin,
        ep_psizemax=args.fdp_psizemax,
        ep_maxpolyx=args.fdp_maxpolyx,
        ep_osizeopt=args.fdp_oligoosize,
        ep_ominsize=args.fdp_oligominsize,
        ep_omaxsize=args.fdp_oligomaxsize,
        ep_otmopt=args.fdp_oligootm,
        ep_otmmin=args.fdp_oligomintm,
        ep_otmmax=args.fdp_oligomaxtm,
        ep_ogcopt=args.fdp_oligoogcpercent,
        ep_ogcmin=args.fdp_oligomingc,
        ep_ogcmax=args.fdp_oligomaxgc,
        logfile=args.logfilestem + "_eprimer3.log",
    )
    logger.info(
        "pdp eprimer3 namespace:\n\t%s",
        "\n\t".join(
            ["{}: {}".format(_[0], _[1]) for _ in eprimer3_namespace.__dict__.items()]
        ),
    )
    eprimer3_logger = build_logger(
        "EPRIMER3 (find_differential_primers.py)", eprimer3_namespace
    )
    logger.info("Running pdp eprimer3 step")
    subcommands.subcmd_eprimer3(eprimer3_namespace, eprimer3_logger)
    return outconfpath


def run_dedupe(args, logger):
    """Run pdp dedupe step with find_differential_primers.py script

    :param args:  Namespace with command-line arguments
    :param logger:  logging.logger object
    """
    #  Generate namespace and logger for pdp dedupe step
    logger.info("Preparing pdp dedupe step")
    outconfpath = os.path.join(args.fdp_outdir, "deduped.json")
    dedupe_namespace = Namespace(
        verbose=args.verbose,
        disable_tqdm=args.fdp_disable_tqdm,
        infilename=args.infilename,
        outfilename=outconfpath,
        dd_dedupedir=args.fdp_outdir,
        logfile=args.logfilestem + "_dedupe.log",
    )
    logger.info(
        "pdp dedupe namespace:\n\t%s",
        "\n\t".join(
            ["{}: {}".format(_[0], _[1]) for _ in dedupe_namespace.__dict__.items()]
        ),
    )
    dedupe_logger = build_logger(
        "DEDUPE (find_differential_primers.py)", dedupe_namespace
    )
    logger.info("Running pdp dedupe step")
    subcommands.subcmd_dedupe(dedupe_namespace, dedupe_logger)
    return outconfpath
