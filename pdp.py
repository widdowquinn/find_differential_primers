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

from diagnostic_primers import multiprocessing, process, prodigal, eprimer3,\
    sge, sge_jobs


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
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
    parser_main = ArgumentParser(prog='pdp.py')
    subparsers = parser_main.add_subparsers(title='subcommands',
                                            description='valid subcommands',
                                            help='additional help')

    # A 'common' parser, with the shared commands for all subcommands
    parser_common = ArgumentParser(add_help=False)
    parser_common.add_argument('infilename',
                               help='path to configuration file')
    parser_common.add_argument('-l', '--logfile', dest='logfile',
                               action='store', default=None,
                               help='logfile location')
    parser_common.add_argument('-v', '--verbose', action='store_true',
                               dest='verbose', default=False,
                               help='report progress to log')
    parser_common.add_argument('-s', '--scheduler', dest='scheduler',
                               action='store', default='multiprocessing',
                               help='Job scheduler [multiprocessing|SGE]')
    parser_common.add_argument('-w', '--workers', dest='workers',
                               action='store', default=None, type=int,
                               help='Number of parallel workers to use')

    # Subcommand parsers
    parser_process = subparsers.add_parser('process',
                                           parents=[parser_common])
    parser_prodigal = subparsers.add_parser('prodigal', aliases=['prod'],
                                            parents=[parser_common])
    parser_eprimer3 = subparsers.add_parser('eprimer3', aliases=['e3'],
                                            parents=[parser_common])
    parser_blastscreen = subparsers.add_parser('blastscreen', aliases=['bs'],
                                               parents=[parser_common])
    parser_primersearch = subparsers.add_parser('primersearch',
                                                aliases=['ps'],
                                                parents=[parser_common])
    parser_classify = subparsers.add_parser('classify', aliases=['cl'],
                                            parents=[parser_common])

    # Config file processing options - subcommand process
    parser_process.add_argument('outfilename',
                                help='Path to write new configuration file')
    parser_process.add_argument('--validate', action='store_true',
                                dest='validate', default=False,
                                help='Validate config file, then exit')

    # CDS prediction options - subcommand prodigal
    parser_prodigal.add_argument('outfilename',
                                 help='Path to write new configuration file')
    parser_prodigal.add_argument('--prodigal', dest='prodigal_exe',
                                 action='store', default='prodigal',
                                 help='path to Prodigal executable')
    parser_prodigal.add_argument('--outdir', dest='prodigaldir',
                                 action='store', default='prodigal',
                                 help='path to directory for Prodigal output')
    parser_prodigal.add_argument('-f', '--force', dest='prodigalforce',
                                 action='store_true', default=False,
                                 help='Overwrite old Prodigal output')

    # Primer prediction options - subcommand eprimer3
    parser_eprimer3.add_argument('outfilename',
                                 help='Path to write new configuration file')
    parser_eprimer3.add_argument('--eprimer3', dest='eprimer3_exe',
                                 action='store', default='eprimer3',
                                 help='path to ePrimer3 executable')
    parser_eprimer3.add_argument('--outdir', dest='eprimer3_dir',
                                 action='store', default='eprimer3',
                                 help='path to directory for ePrimer3 output')
    parser_eprimer3.add_argument('-f', '--force', dest='eprimer3_force',
                                 action='store_true', default=False,
                                 help='Overwrite old ePrimer3 output')
    parser_eprimer3.add_argument('--numreturn', dest='ep_numreturn',
                                 action='store', default=10, type=int,
                                 help='number of primers to return')
    parser_eprimer3.add_argument('--osize', dest='ep_osize',
                                 action='store', default=20, type=int,
                                 help='optimal size for primer oligo')
    parser_eprimer3.add_argument('--minsize', dest='ep_minsize',
                                 action='store', default=18, type=int,
                                 help='minimum size for primer oligo')
    parser_eprimer3.add_argument('--maxsize', dest='ep_maxsize',
                                 action='store', default=22, type=int,
                                 help='maximum size for primer oligo')
    parser_eprimer3.add_argument('--opttm', dest='ep_opttm',
                                 action='store', default=59, type=int,
                                 help='optimal Tm for primer oligo')
    parser_eprimer3.add_argument('--mintm', dest='ep_mintm',
                                 action='store', default=58, type=int,
                                 help='minimum Tm for primer oligo')
    parser_eprimer3.add_argument('--maxtm', dest='ep_maxtm',
                                 action='store', default=60, type=int,
                                 help='maximum Tm for primer oligo')
    parser_eprimer3.add_argument('--ogcpercent', dest='ep_ogcpercent',
                                 action='store', default=55, type=int,
                                 help='optimal %%GC for primer oligo')
    parser_eprimer3.add_argument('--mingcpercent', dest='ep_mingc',
                                 action='store', default=30, type=int,
                                 help='minimum %%GC for primer oligo')
    parser_eprimer3.add_argument('--maxgcpercent', dest='ep_maxgc',
                                 action='store', default=80, type=int,
                                 help='maximum %%GC for primer oligo')
    parser_eprimer3.add_argument('--psizeopt', dest='ep_psizeopt',
                                 action='store', default=100, type=int,
                                 help='optimal size of amplified region')
    parser_eprimer3.add_argument('--psizemin', dest='ep_psizemin',
                                 action='store', default=50, type=int,
                                 help='minimum size of amplified region')
    parser_eprimer3.add_argument('--psizemax', dest='ep_psizemax',
                                 action='store', default=150, type=int,
                                 help='maximum size of amplified region')
    parser_eprimer3.add_argument('--maxpolyx', dest='ep_maxpolyx',
                                 action='store', default=3, type=int,
                                 help='maximum run of repeated nucleotides ' +
                                 'in primer')
    parser_eprimer3.add_argument('--hybridprobe', dest='ep_hybridprobe',
                                 action='store_true', default=False,
                                 help='design a reporter oligo')
    parser_eprimer3.add_argument('--oligoosize', dest='ep_osizeopt',
                                 action='store', default=20, type=int,
                                 help='optimal size for internal oligo')
    parser_eprimer3.add_argument('--oligominsize', dest='ep_ominsize',
                                 action='store', default=13, type=int,
                                 help='minimum size for internal oligo')
    parser_eprimer3.add_argument('--oligomaxsize', dest='ep_omaxsize',
                                 action='store', default=30, type=int,
                                 help='maximum size for internal oligo')
    parser_eprimer3.add_argument('--oligootm', dest='ep_otmopt',
                                 action='store', default=69, type=int,
                                 help='optimal Tm for internal oligo')
    parser_eprimer3.add_argument('--oligomintm', dest='ep_otmmin',
                                 action='store', default=68, type=int,
                                 help='minimum Tm for internal oligo')
    parser_eprimer3.add_argument('--oligomaxtm', dest='ep_otmmax',
                                 action='store', default=70, type=int,
                                 help='maximum Tm for internal oligo')
    parser_eprimer3.add_argument('--oligoogcpercent', dest='ep_ogcopt',
                                 action='store', default=55, type=int,
                                 help='optimal %%GC for internal oligo')
    parser_eprimer3.add_argument('--oligomingcpercent', dest='ep_ogcmin',
                                 action='store', default=30, type=int,
                                 help='minimum %%GC for internal oligo')
    parser_eprimer3.add_argument('--oligomaxgcpercent', dest='ep_ogcmax',
                                 action='store', default=80, type=int,
                                 help='maximum %%GC for internal oligo')
# Commented out until Biopython command line code catches up with new
# EMBOSS options
#    parser_eprimer3.add_argument("--oligomaxpolyx", dest="ep_opolymax",
#                                 action="store",
#                                 default=3, type=int,
#                                 help="maximum run of repeated nucleotides " +
#                                 "in internal primer")

    # Parse arguments
    return parser_main.parse_args()


def load_config_file():
    """Load config file and return a GenomeCollection."""
    try:
        gc = process.load_collection(args.infilename, name='pdp.py')
    except:
        logger.error('Could not parse config file %s (exiting)' %
                     args.infilename)
        logger.error(last_exception())
        sys.exit(1)
    logger.info('Parsed config file %s: %d sequences in %d groups' %
                (args.infilename, len(gc), len(gc.groups())))
    logger.info('Diagnostic groups:\n%s' %
                '\n'.join(['\t%s' % g for g in gc.groups()]))
    return gc


def run_parallel_jobs(clines):
    """Run the passed command-lines in parallel."""
    logger.info('Running jobs using scheduler: %s' % args.scheduler)
    # Pass lines to scheduler and run
    if args.scheduler == 'multiprocessing':
        retvals = multiprocessing.run(clines, workers=args.workers,
                                      verbose=args.verbose)
        if sum([r.returncode for r in retvals]):
            logger.error('At least one run has problems (exiting).')
            for retval in retvals:
                if retval.returncode != 0:
                    logger.error('Failing command: %s' % retval.args)
                    logger.error('Failing stderr:\n %s' % retval.stderr)
            sys.exit(1)
        else:
            logger.info('Runs completed without error.')
    elif args.scheduler == 'SGE':
        joblist = [sge_jobs.Job("pdp_%06d" % idx, cmd) for idx, cmd in
                   enumerate(clines)]
        sge.run_dependency_graph(joblist, verbose=True, logger=logger)
    else:
        raise ValueError('Scheduler must be one of ' +
                         '[multiprocessing|SGE], got %s' % args.scheduler)


def log_clines(clines):
    """Log command-lines, one per line."""
    logger.info('...%d commands returned:\n%s' %
                (len(clines), '\n'.join(['\t%s' % c for c in clines])))


def subcmd_process():
    """Run the process subcommand operations. Here, the goal is to parse the
    input config file, verify that the config file can be read, report some
    basic information back (if verbose), and clean up sequence data by
    stitching multiple contigs/sequences together and replacing ambiguity
    symbols with 'N'.

    If sequence data needs to be stitched, or symbols replaced, then new
    sequence files are produced and written (if --validate is not in
    operation)

    A new config file, pointing to the revised files, is written out (if
    --validate is not in operation).

    """
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
    logger.info('Checking whether input sequences require stitching, ' +
                'or have non-N ambiguities.')
    for g in gc.data:
        if g.needs_stitch:
            logger.info('%s requires stitch' % g.name)
            if not args.validate:
                g.stitch()
        else:
            logger.info('%s does not require stitch' % g.name)
        if g.has_ambiguities():
            logger.info('%s contains non-N ambiguities' % g.name)
            if not args.validate:
                g.replace_ambiguities()
        else:
            logger.info('%s does not contain non-N ambiguities' % g.name)
        logger.info('Sequence file: %s' % g.seqfile)

    # Write post-processing config file and exit
    if not args.validate:
        logger.info('Writing processed config file to %s' %
                    args.outfilename)
        gc.write(args.outfilename)
    sys.exit(0)


def subcmd_prodigal():
    """Run Prodigal to predict bacterial CDS on the input sequences."""
    gc = load_config_file()

    # Build command-lines for Prodigal and run
    logger.info('Building Prodigal command lines...')
    if args.prodigalforce:
        logger.warning('Forcing Prodigal to run. This may overwrite ' +
                       'existing output.')
    else:
        logger.info('Prodigal will fail if output directory exists')
    clines = prodigal.build_commands(gc, args.prodigal_exe,
                                     args.prodigaldir, args.prodigalforce)
    log_clines(clines)
    run_parallel_jobs(clines)

    # Add Prodigal output files to the GenomeData objects and write
    # the config file
    for g in gc.data:
        g.features = g.cmds['prodigal'].split()[-1].strip()
        logger.info('%s feature file:\t%s' % (g.name, g.features))
    logger.info('Writing new config file to %s' % args.outfilename)
    gc.write(args.outfilename)
    sys.exit(0)


def subcmd_eprimer3():
    """Run ePrimer3 to design primers for each input sequence."""
    gc = load_config_file()

    # Build command-lines for ePrimer3 and run
    logger.info('Building ePrimer3 command lines...')
    if args.eprimer3_force:
        logger.warning('Forcing ePrimer3 to run. This may overwrite ' +
                       'existing output.')
    else:
        logger.info('ePrimer3 may fail if output directory exists')
    clines = eprimer3.build_commands(gc, args.eprimer3_exe,
                                     args.eprimer3_dir,
                                     args.eprimer3_force, vars(args))
    pretty_clines = [str(c).replace(' -', '\\\n\t\t-') for c in clines]
    log_clines(pretty_clines)
    run_parallel_jobs(clines)

    # Load ePrimer3 data for each input sequence, and write JSON representation
    # Record JSON representation in the object
    for g in gc.data:
        logger.info("Loading/naming primers in %s" % g.cmds['ePrimer3'].outfile)
        primers, outfname = eprimer3.load_primers(g.cmds['ePrimer3'].outfile)
        logger.info('Named primers ePrimer3 file created:\t%s' % outfname)
        outfname = os.path.splitext(outfname)[0] + '.json'
        logger.info('Writing primers to %s' % outfname)
        eprimer3.primers_to_json(primers, outfname)
        g.primers = outfname

    logger.info('Writing new config file to %s' % args.outfilename)
    gc.write(args.outfilename)
    sys.exit(0)


def subcmd_blastscreen():
    """Screen primer sequences against a local BLAST database."""
    gc = load_config_file()

    for g in gc.data:
        logger.info("Loading primer data from %s" % g.primers)
        fastafname = eprimer3.json_to_fasta(g.primers)
        logger.info("Wrote primer FASTA sequences to %s" % fastafname)


###
# Run as script
if __name__ == '__main__':

    # Parse command-line
    assert len(sys.argv) > 1, 'pdp.py requires a valid subcommand'
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
            logger.error('Could not open %s for logging' %
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
    logger.info('Processed arguments: %s' % args)
    logger.info('command-line: %s' % ' '.join(sys.argv))

    # TODO: turn the if statements below into a distribution dictionary
    subcmds = {'process': subcmd_process,
               'prodigal': subcmd_prodigal, 'prod': subcmd_prodigal,
               'eprimer3': subcmd_eprimer3, 'e3': subcmd_eprimer3,
               'blastscreen': subcmd_blastscreen, 'bs': subcmd_blastscreen}
    try:
        subcmds[subcmd]()
    except KeyError:
        raise NotImplementedError

    # Exit as if all is well
    sys.exit(0)
