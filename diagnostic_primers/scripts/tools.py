#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tools.py

Provides useful functions for the pdp.py script

- config:        interact with config files
- prodigal:      run prodigal bacterial gene prediction
- eprimer3:      predict primers with EMBOSS ePrimer3
- primersearch:  run in-silico hybridisation with EMBOSS PrimerSearch
- blastscreen:   screen primers with BLASTN
- classify:      classify primers

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
import sys
import traceback

from diagnostic_primers import config, multiprocessing, sge, sge_jobs, PDPException


class PDPScriptError(PDPException):
    """Exception thrown when script fails."""

    def __init__(self, msg="Error in `pdp` script"):
        PDPException.__init__(self, msg)


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))


# Load PDPCollection from .tab file
def load_config_tab(args, logger):
    """Load tab format config to PDPCollection."""
    pdpc = config.PDPCollection()
    try:
        pdpc.from_tab(args.infilename)
    except config.ConfigSyntaxError:
        logger.error("Could not read config file %s (exiting)", args.infilename)
        logger.error(last_exception())
        raise SystemExit(1)
    return pdpc


# Load PDPCollection from .json file
def load_config_json(args, logger):
    """Load JSON format config to PDPCollection."""
    pdpc = config.PDPCollection()
    try:
        logger.info("Loading config from %s", args.infilename)
        pdpc.from_json(args.infilename)
    except config.ConfigSyntaxError:
        logger.error("Could not read config file %s (exiting)", args.infilename)
        logger.error(last_exception())
        raise SystemExit(1)
    except FileNotFoundError:
        logger.error("Config file %s does not exist (exiting)", args.infilename)
        logger.error(last_exception)
        raise SystemExit(1)
    return pdpc


# Report a list of command lines to a logger, in pretty format
def log_clines(clines, logger):
    """Log command-lines, one per line."""
    logger.info(
        "...%d commands returned:\n%s"
        % (len(clines), "\n".join(["  %s" % c for c in clines]))
    )


# Pass jobs to the appropriate scheduler
def run_parallel_jobs(clines, args, logger):
    """Run the passed command-lines in parallel."""
    logger.info("Running jobs using scheduler: %s" % args.scheduler)
    # Pass lines to scheduler and run
    if args.scheduler == "multiprocessing":
        retvals = multiprocessing.run(clines, workers=args.workers)
        if retvals != 0:
            logger.error("At least one run has problems (exiting).")
            raise SystemExit(1)
        else:
            logger.info("Runs completed without error.")
    elif args.scheduler == "SGE":
        joblist = [
            sge_jobs.Job("pdp_%06d" % idx, cmd) for idx, cmd in enumerate(clines)
        ]
        sge.run_dependency_graph(joblist, logger=logger)
    else:
        raise ValueError(
            "Scheduler must be one of "
            + "[multiprocessing|SGE], got %s" % args.scheduler
        )


# Test whether the passed PDPCollection has primersearch output linked
def has_primersearch(coll):
    """Returns True if the passed PDPCollection has primersearch output

    - coll      PDPCollection describing genomes for a run
    """
    if None in [genome.primersearch for genome in coll.data]:
        return False
    return True


# Create an output directory, if it doesn't already exist
def create_output_directory(outdirname, force, logger):
    """Create output directory, respecting cmd-line arguments

    - outdirname          path to output directory
    - force               force creation True/False
    - logger              logger for program

    Attempts to create a directory, but will not attempt to overwrite
    unless the command-line argument to force is set True
    """
    if os.path.exists(outdirname) and not force:
        logger.error(
            "Output directory %s exists - not overwriting (exiting)", outdirname
        )
        raise SystemExit(1)
    if force:
        logger.warning(
            "Output directory %s exists: forcing use and potential overwrite",
            outdirname,
        )
    else:
        logger.info("Creating output directory %s", outdirname)
    os.makedirs(outdirname, exist_ok=True)


def chunk(iterable, size):
    for i in range(0, len(iterable), size):
        yield iterable[i : i + size]


# Collect existing output files when in recovery mode
def collect_existing_output(dirpath, step, args):
    """Returns a list of existing output files at dirpath

    :param dirpath:       path to existing output directory
    :param step:          the pipeline step being run
    :param args:          command-line arguments for the run
    """
    # The suffixes dict has step:file_suffix pairs. Step refers to the pdp
    # pipeline step, and the file_suffix the the suffix of files that are
    # written to the output at the end of a called third-party application
    # job. So, the prodigal step calls prodigal, which writes a .gff file.
    #  When in --recovery mode, the presence of a .gff file is taken to
    #  indicate that the third-party tool ran to completion.
    suffixes = {
        "eprimer3": ".eprimer3",
        "primer3": ".primer3",
        "prodigal": ".gff",
        "alnvar": ".filter",
        "blastscreen": ".blasttab",
        "primersearch": ".primersearch",
        "extract": ".aln",
    }
    try:
        existingfiles = [
            fname
            for fname in os.listdir(dirpath)
            if os.path.splitext(fname)[-1] == suffixes[step]
        ]
    except KeyError:
        raise PDPScriptError(
            "PDP step {} not recognised when collecting output".format(step)
        )
    return existingfiles
