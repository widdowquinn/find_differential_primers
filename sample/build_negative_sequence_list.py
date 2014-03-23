#!/usr/bin/env python
#
# build_negative_list.py
#
# Script that takes as input a FASTA file, and a list of sequence identifiers,
# and returns a list of sequence identifiers which are those in the FASTA file,
# but not in the initial list of sequence identifiers.
#
# (c) The James Hutton Institute 2013
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
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#=============
# IMPORTS

import logging
import logging.handlers
import sys
import time
import traceback

try:
    from Bio import SeqIO
except ImportError:
    print "Biopython required for script, but not found (exiting)"
    sys.exit(1)
from optparse import OptionParser

#=============
# FUNCTIONS

# Parse command-line
def parse_cmdline():
    """ Parse command-line arguments for the script
    """
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfile", dest="outfilename",
                      action="store", default=None,
                      help="Output filename")
    parser.add_option("-f", "--fastafile", dest="fastafilename",
                      action="store", default=None,
                      help="Input FASTA filename")
    parser.add_option("-i", "--listfile", dest="listfilename",
                      action="store", default=None,
                      help="Input ID list filename")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    parser.add_option("-l", "--logfile", dest="logfile",
                      action="store", default=None,
                      help="Logfile location")
    return parser.parse_args()

# Report last exception as string
def last_exception():
    """ Returns last exception as a string, for use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))

# Return a set of IDs from the input filename
def get_fasta_ids():
    """ Returns a set of strings. These are the sequence IDs from the
        input FASTA file.
    """
    try:
        seqids = set([s.id for s in
                      SeqIO.parse(options.fastafilename, 'fasta')])
    except IOError:
        logger.error("Could not extract sequence IDs from %s",
                     options.fastafilename)
        logger.error(last_exception())
        sys.exit(1)
    logger.info("Recovered %d sequence IDs", len(seqids))
    return seqids

# Return the sequence IDs that are in the input FASTA file, but not the
# input list
def get_negative_set():
    """ Get a set of sequence IDs from the input FASTA file, and the list of
        input sequenc IDs, then return the IDs that are in the FASTA file,
        but not the input list.
    """
    fastaset = get_fasta_ids()
    try:
        listids = set([l.strip() for l in open(options.listfilename, 'rU') if
                       len(l.strip()) and not l.startswith('#')])
    except IOError:
        logger.error("Could not extract sequence IDs from %s",
                     options.listfilename)
        logger.error(last_exception())
        sys.exit(1)
    logger.info("Recovered %d IDs from the input list", len(listids))
    negset = fastaset - listids
    logger.info("Identified %d IDs for the negative set", len(negset))
    return negset


#=============
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    # options are all options - no arguments
    options, args = parse_cmdline()

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    # err_handler points to sys.stderr
    # err_handler_file points to a logfile, if named
    logger = logging.getLogger('build_negative_sequence_list.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = \
                  logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.logfile is not None:
        try:
            logstream = open(options.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except IOError:
            logger.error("Could not open %s for logging",
                         options.logfile)
            sys.exit(1)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    logger.info('# build_negative_sequence_list.py logfile')
    logger.info('# Run: %s', time.asctime())

    # Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    # Have we got an input and output file? If not, exit.
    if options.fastafilename is None:
        logger.error("No input FASTA file name (exiting)")
        sys.exit(1)
    logger.info("Input FASTA file: %s", options.fastafilename)
    if options.listfilename is None:
        logger.error("No input seqID list file name (exiting)")
        sys.exit(1)
    logger.info("Input list file: %s", options.listfilename)
    if options.outfilename is None:
        logger.error("No output filename (exiting)")
        sys.exit(1)
    logger.info("Output file: %s", options.outfilename)

    # Get the set of sequence IDs that are in the FASTA file,
    # but not the input list
    neg_id_set = get_negative_set()

    # Write the negative ID set to file
    try:
        with open(options.outfilename, 'w') as fh:
            fh.write('\n'.join(list(neg_id_set)))
    except IOError:
        logger.error("Could not write negative ID set to %s",
                     options.outfilename)
        logger.error(last_exception())
        sys.exit(1)
