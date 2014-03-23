#!/usr/bin/env python
#
# stitch_six_frame_stops.py
#
# Takes an input (multiple) FASTA sequence file, and replaces all runs of
# N with the sequence NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN, which
# contains start and stop codons in all frames.  All the sequences in the
# input file are then stitched together with the same sequence.
#
# Overall, the effect is to replace all regions of base uncertainty with
# the insert sequence, forcing stops and starts for gene-calling, to avoid
# chimeras or frame-shift errors due to placement of Ns.
#
# We also convert other ambiguity codons to Ns
#
# This script is intended for use in assembly pipelines, where contigs are
# provided in the correct (or, at least, an acceptable) order.
#
# (c) The James Hutton Institute 2011
# Authors: Leighton Pritchard, Benjamin Leopold, Michael Robeson
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

###
# IMPORTS

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from optparse import OptionParser

import logging
import logging.handlers
import re
import sys

SEPARATOR = 'NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN'

###
# FUNCTIONS


# Parse command-line
def parse_cmdline():
    """ Parse command-line arguments
    """
    usage = "usage: %prog [options] <organism_type> <infile>"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfile", dest="outfilename",
                      action="store", default=None,
                      help="Output filename")
    parser.add_option("-i", "--infile", dest="infilename",
                      action="store", default=None,
                      help="Input filename")
    parser.add_option("-n", "--nrun", dest="nrun",
                      action="store", default=2,
                      help="Minimum number of Ns in substituted run " +\
                           "(shorter runs not substituted)")
    parser.add_option("--noambiguity", dest="noambiguity",
                      action="store_true", default=False,
                      help="Do not replace IUPAC ambiguity codes other than " +\
                           "N with N")
    parser.add_option("--id", dest="seqid",
                      action="store", default=None,
                      help="ID/Accession for the output stitched sequence")
    parser.add_option("--nodesc", dest="nodesc",
                      action="store_true", default=False,
                      help="Suppress sequence description string")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    return parser.parse_args()


# Replace runs of N in each sequence with the separator
def stitch_ns(sequences, nrun):
    """ Loop over each input sequence in sequences, and replace any
        occurrences of N longer than nrun with the separator sequence.
    """
    nrun = int(nrun)              # Ensure we're dealing with an int
    new_sequences = []
    for seq in sequences:
        seqdata = str(seq.seq)
        ncount = seqdata.count('n') + seqdata.count('N')
        logger.info("%d Ns in %s", ncount, seq.id)
        if 0 == ncount:
            logger.info("Keeping unmodified %s", seq.id)
            new_sequences.append(seq)
            continue
        new_seq, repcount = re.subn('[nN]{%d,}' % nrun, SEPARATOR, seqdata)
        logger.info("Replaced %d runs of N runs >= %d in %s",
                    repcount, nrun, seq.id)
        new_seqrecord = SeqRecord(Seq(new_seq.upper()),
                                  id=''.join([seq.id, "_N-replaced"]),
                                  name=''.join([seq.name, "_N-replaced"]),
                                  description=seq.description)
        logger.info("New SeqRecord created:\n%s", new_seqrecord)
        new_sequences.append(new_seqrecord)
    return new_sequences


# Stitch passed sequences together with the separator
def stitch_seqs(sequences, seqid, nodesc, noambiguity):
    """ Stitch the passed sequences together, giving the result the passed ID,
        but a description that derives from each of the passed sequences
    """
    new_seq = SEPARATOR.join([str(s.seq) for s in sequences])
    if noambiguity:
        new_seq, repcount = re.subn('[^acgtACGTnN]', 'N', new_seq)
        logger.info("Substituted %d ambiguity bases", repcount)
    new_id = seqid
    new_name = seqid + "_stitched"
    if nodesc:
        new_desc = ''
    else:
        new_desc = '+'.join([s.id for s in sequences])
    stitched_seq = SeqRecord(Seq(new_seq), id=new_id, name=new_name,
                             description=new_desc)
    logger.info("Created stitched sequence (len:%d):\n%s",
                len(stitched_seq), stitched_seq)
    return stitched_seq

###
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    # options are options, arguments are the .sff files
    options, args = parse_cmdline()

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    logger = logging.getLogger('stitch_six_frame_stops.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = \
                  logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    # Load data from file (STDIN if no filename provided)
    if options.infilename is None:
        inhandle = sys.stdin
        logger.info("Using STDIN for input")
    else:
        inhandle = open(options.infilename, 'rU')
    data = list(SeqIO.parse(inhandle, 'fasta'))
    inhandle.close()

    # Stitch individual sequences
    data = stitch_ns(data, options.nrun)

    if options.seqid is None:
        options.seqid = '_'.join(data[0].description.split(' ')[1:])

    # Stitch all sequences together
    stitchedseq = stitch_seqs(data, options.seqid, options.nodesc,
                              options.noambiguity)

    # Write the stitched sequence to file (or STDOUT if no filename provided)
    if options.outfilename is None:
        outhandle = sys.stdout
    else:
        outhandle = open(options.outfilename, 'w')
    SeqIO.write([stitchedseq], outhandle, 'fasta')
    outhandle.close()
