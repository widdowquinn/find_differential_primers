#!/usr/bin/env python3
#
# predict_diagnostic_primers.py
#
# This script aids prediction of diagnostic PCR primers from nucleotide 
# sequence files, where sequence files are associated with group/class
# identifiers
# 
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016 The James Hutton Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import logging
import os
import shutil
import sys
import time
import traceback

from argparse import ArgumentParser

from diagnostic_primers import multiprocessing, process, prodigal

# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))

# Process command-line
def parse_cmdline(args):
    """Parse command-line arguments for script.

    The script offers a single main parser, with subcommands for the actions:

    process - process/check input data (stitch fragments/fix sequence problems)
    prodigal - run Prodigal to predict CDS features, and update .conf file
    eprimer3 - run ePrimer3 to predict primer pairs, and update .conf file
    check_blast - check/filter designed primers against negative example
                  BLAST database
    primersearch - check/filter designed primers against complete genome
                   negative examples
    classify - classify designed primers against input genome/classes
    """
    # Main parent parser
    parser_main = ArgumentParser(prog="pdp.py")
    subparsers = parser_main.add_subparsers(title="subcommands",
                                            description="valid subcommands",
                                            help="additional help")

    # A 'common' parser, with the shared commands for all subcommands
    parser_common = ArgumentParser(add_help=False)
    parser_common.add_argument("infilename",
                               help="path to configuration file")
    parser_common.add_argument("-l", "--logfile", dest="logfile",
                               action="store", default=None,
                               help="logfile location")
    parser_common.add_argument("-v", "--verbose", action="store_true",
                               dest="verbose", default=False,
                               help="report progress to log")
    parser_common.add_argument("-s", "--scheduler", dest="scheduler",
                               action="store", default="multiprocessing",
                               help="Job scheduler [multiprocessing|SGE]")
    parser_common.add_argument("-w", "--workers", dest="workers",
                               action="store", default=None, type=int,
                               help="Number of parallel workers to use")

    # Subcommand parsers
    parser_process = subparsers.add_parser('process', aliases=['pr'],
                                           parents=[parser_common])
    parser_prodigal = subparsers.add_parser('prodigal',
                                            parents=[parser_common])
    parser_eprimer3 = subparsers.add_parser('eprimer3', aliases=['e3'],
                                            parents=[parser_common])
    parser_blastcheck = subparsers.add_parser('blastcheck', aliases=['bc'],
                                              parents=[parser_common])
    parser_primersearch = subparsers.add_parser('primersearch',
                                                aliases=['ps'],
                                                parents=[parser_common])
    parser_classify = subparsers.add_parser('classify', aliases=['cl'],
                                            parents=[parser_common])
    
    # Config file processing options - subcommand process
    parser_process.add_argument("outfilename",
                               help="Path to write new configuration file")
    parser_process.add_argument('--validate', action="store_true",
                                dest='validate', default=False,
                                help="Validate config file, then exit")

    # CDS prediction options - subcommand prodigal
    parser_prodigal.add_argument("outfilename",
                               help="Path to write new configuration file")
    parser_prodigal.add_argument("--prodigal", dest="prodigal_exe",
                                 action="store", default="prodigal",
                                 help='path to Prodigal executable')
    parser_prodigal.add_argument("--outdir", dest="prodigaldir",
                                 action="store", default="prodigal",
                                 help="path to directory for Prodigal output")
    parser_prodigal.add_argument("-f", "--force", dest="prodigalforce",
                                 action="store_true", default=False,
                                 help="Overwrite old Prodigal output")

    # Primer prediction options - subcommand eprimer3
    parser_eprimer3.add_argument("outfilename",
                               help="Path to write new configuration file")
    parser_eprimer3.add_argument("--eprimer3", dest="eprimer3_exe",
                                 action="store", default="eprimer3",
                                 help='path to ePrimer3 executable')
    parser_eprimer3.add_argument("--outdir", dest="eprimer3dir",
                                 action="store", default="eprimer3",
                                 help="path to directory for ePrimer3 output")
    parser_eprimer3.add_argument("-f", "--force", dest="eprimer3force",
                                 action="store_true", default=False,
                                 help="Overwrite old ePrimer3 output")

    # Parse arguments
    return parser_main.parse_args()


def load_config_file():
    """Load config file and return a GenomeCollection."""
    try:
        gc = process.load_collection(args.infilename, name="pdp.py")
    except:
        logger.error("Could not parse config file %s (exiting)" % 
                     args.infilename)
        logger.error(last_exception())
        sys.exit(1)
    logger.info("Parsed config file %s: %d sequences in %d groups" %
                (args.infilename, len(gc), len(gc.groups())))
    logger.info("Diagnostic groups:\n%s" %
                '\n'.join(["\t%s" % g for g in gc.groups()]))
    return gc


def run_parallel_jobs(clines):
    """Run the passed command-lines in parallel."""
    logger.info("Running jobs using scheduler: %s" % args.scheduler)
    # Pass lines to scheduler and run
    if args.scheduler == 'multiprocessing':
        retval = multiprocessing.run(clines)
        if retval:
            logger.error("At least one Prodigal run has problems (exiting).")
            sys.exit(1)
        else:
            logger.info("Prodigal runs completed without error.")
    elif args.scheduler == 'SGE':
        sge.run(clines)
    else:
        raise ValueError("Scheduler must be one of " +
                         "[multiprocessing|SGE], got %s" % args.scheduler)


###
# Run as script
if __name__ == '__main__':

    # Parse command-line
    assert len(sys.argv) > 1, "pdp.py requires a valid subcommand"
    args = parse_cmdline(sys.argv)
    subcmd = sys.argv[1]

    # Set up logging
    logger = logging.getLogger('pdp.py: %s' % time.asctime())
    t0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         args.logfile)
            logger.error(last_exception())
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info("Processed arguments: %s" % args)
    logger.info("command-line: %s" % ' '.join(sys.argv))

    # PROCESS
    # If we're running the process operation, the goal is to parse the input
    # config file, verify that the config file can be read, report some basic
    # information back (if verbose), and clean up sequence data by stitching it
    # and replacing ambiguity symbols with 'N'.
    # If sequence data needs to be stitched, or symbols replaced, then new
    # sequence files are produced and written (if --validate is not in
    # operation)
    # A new config file, pointing to the revised files, is written out (if
    # --validate is not in operation).
    if subcmd == 'process':
        gc = load_config_file()
        
        # Do sequences need to be stitched or their ambiguities replaced?
        # If --validate is active, we report only and do not modify.
        # Note that the underlying code doesn't require us to check whether
        # the sequence files need stitching or symbols replacing, and
        # we could more efficiently use
        # for g in gc.data:
        #     g.stitch()
        #     g.replace_ambiguities()
        # but we're being helpfully verbose.
        logger.info("Checking whether input sequences require stitching, " +\
                    "or have non-N ambiguities.")
        for g in gc.data:
            if g.needs_stitch:
                logger.info("%s requires stitch" % g.name)
                if not args.validate:
                    g.stitch()
            else:
                logger.info("%s does not require stitch" % g.name)
            if g.has_ambiguities():
                logger.info("%s contains non-N ambiguities" % g.name)
                if not args.validate:
                    g.replace_ambiguities()
            else:
                logger.info("%s does not contain non-N ambiguities" % g.name)
            logger.info("Sequence file: %s" % g.seqfile)

        # Write post-processing config file and exit
        if not args.validate:
            logger.info("Writing processed config file to %s" %
                        args.outfilename)
            gc.write(args.outfilename)
        sys.exit(0)


    # PRODIGAL
    # The prodigal subcommand is used if the user wants to run Prodigal to
    # predict CDS on the input sequences.
    if subcmd == 'prodigal':
        gc = load_config_file()

        # Build command-lines for Prodigal and run
        logger.info("Building Prodigal command lines...")
        if args.prodigalforce:
            logger.warning("Forcing Prodigal to run. This may overwrite " +\
                           "existing output.")
        else:
            logger.info("Prodigal will fail if output directory exists")
        clines = prodigal.build_commands(gc, args.prodigal_exe,
                                         args.prodigaldir, args.prodigalforce)
        logger.info("...%d commands returned:\n%s" %
                    (len(clines), '\n'.join(['\t%s' % c for c in clines])))
        run_parallel_jobs(clines)

        # Add Prodigal output files to the GenomeData objects and write
        # the config file
        for g in gc.data:
            g.features = g.cmds['prodigal'].split()[-1].strip()
            logger.info("%s feature file:\t%s" % (g.name, g.features))
        logger.info("Writing new config file to %s" % args.outfilename)
        gc.write(args.outfilename)
        sys.exit(0)
            
        

    # Exit as if all is well
    sys.exit(0)
