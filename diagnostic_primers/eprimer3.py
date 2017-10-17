#!/usr/bin/env python3
#
# eprimer3.py
#
# Code to conduct primer prediction with ePrimer3
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

import json
import os

from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Emboss import Primer3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class PrimersEncoder(json.JSONEncoder):

    """JSON encoder for Primer3.Primers objects."""

    def default(self, obj):
        if not isinstance(obj, Primer3.Primers):
            return super(PrimersEncoder, self).default(obj)

        # Convert complex Primer3.Primers object to serialisable dictionary
        # and return
        return obj.__dict__


def build_commands(collection, eprimer3_exe, eprimer3_dir,
                   argdict=None):
    """Builds and returns a list of command-lines to run ePrimer3

    The commands will run on each sequence in the passed GenomeCollection.
    """
    clines = []  # Holds command-lines

    # Ensure output directory exists
    os.makedirs(eprimer3_dir, exist_ok=True)

    for g in collection.data:
        if eprimer3_dir is None:
            stem = os.path.splitext(g.seqfile)[0]
        else:
            stempath = os.path.split(os.path.splitext(g.seqfile)[0])
            stem = os.path.join(eprimer3_dir, stempath[-1])
        cline = build_command(eprimer3_exe, g.seqfile, stem, argdict)
        g.cmds['ePrimer3'] = cline
        clines.append(cline)
    return clines


def build_command(eprimer3_exe, seqfile, filestem, argdict=None):
    """Builds and returns ePrimer3 command line.

    The ePrimer3 command uses the Biopython interface
    """
    cline = Primer3Commandline(cmd=eprimer3_exe)
    cline.sequence = seqfile
    cline.auto = True
    cline.outfile = filestem + '.eprimer3'
    if argdict is not None:
        prange = [0, 200]
        args = [(a[3:], v) for a, v in argdict.items() if
                a.startswith('ep_')]
        for arg, val in args:
            if 'psizemin' == arg:
                prange[0] = val
            elif 'psizemax' == arg:
                prange[1] = val
            else:
                setattr(cline, arg, val)
    setattr(cline, 'prange', '%d-%d' % tuple(prange))
    return cline


def load_primers(infname, fmt='eprimer3', noname=False):
    """Load primers from a file.

    The function can load JSON or ePrimer3 files - ePrimer3 by default.
    """
    if fmt in ('ep3', 'eprimer3'):
        return __load_primers_eprimer3(infname, noname)
    elif fmt in ('json', ):
        return __load_primers_json(infname)


def __load_primers_eprimer3(infname, noname=False):
    """Loads and names primers from the passed ePrimer3 file. Returns JSON.

    noname  - Boolean flag. Unless set to true, each primer receives
              a unique name based on the input filename.
    """
    # Load primers with Biopython. This does not respect a 'name' attribute,
    # as the bare ePrimer3 files don't provide names for primer sets.
    with open(infname, 'r') as primerfh:
        ep3primers = Primer3.read(primerfh).primers

    # Create unique identifier for each primer set, based on the input
    # filename
    if not noname:
        for idx, primer in enumerate(ep3primers, 1):
            stem = os.path.splitext(os.path.split(infname)[-1])[0]
            primer.name = "%s_primer_%05d" % (stem, idx)

    # Return primers
    return ep3primers


def __load_primers_json(infname):
    """Loads and returns primers from a JSON file."""
    primers = []
    with open(infname, 'r') as primerfh:
        for pdata in json.load(primerfh):
            primer = Primer3.Primers()
            for k, v in pdata.items():
                setattr(primer, k, v)
            primers.append(primer)
    return primers


def write_primers(primers, outfilename, fmt='fasta'):
    """Write Primer3.Primers to file.

    primers      - collection of Biopython primer objects
    outfilename  - path to output file
    format       - sequence format to write
    """
    if fmt in ('json',):
        __write_primers_json(primers, outfilename)
    elif fmt in ('ep3', 'eprimer3'):
        __write_primers_eprimer3(primers, outfilename)
    elif fmt in ('tsv'):
        __write_primers_tsv(primers, outfilename)
    else:
        __write_primers_seqio(primers, outfilename, fmt)


def __write_primers_seqio(primers, outfilename, fmt):
    """Write primers to file, using SeqIO.

    Returns the number of records written
    """
    seqrecords = []

    for primer in primers:
        seqrecords.append(SeqRecord(Seq(primer.forward_seq),
                                    id=primer.name + '_fwd',
                                    description=''))
        seqrecords.append(SeqRecord(Seq(primer.reverse_seq),
                                    id=primer.name + '_rev',
                                    description=''))
        if len(primer.internal_seq):  # This is '' id no oligo
            seqrecords.append(SeqRecord(Seq(primer.internal_seq),
                                        id=primer.name + '_int',
                                        description=''))

    return SeqIO.write(seqrecords, outfilename, fmt)


def __write_primers_tsv(primers, outfname):
    """Write primers to file in three-column tab-separated format.

    This is required for primersearch input with EMBOSS
    """
    # Don't use more than one newline at the end of the header, or
    # else primersearch treats the blank line as a primer!
    header = '\n'.join(['# EPRIMER3 PRIMERS %s' % outfname,
                        '# Name       FWD        REV']) + '\n'
    with open(outfname, 'w') as outfh:
        outfh.write(header)
        for primer in primers:
            outfh.write('\t'.join([primer.name,
                                   primer.forward_seq,
                                   primer.reverse_seq]) + '\n')


def __write_primers_eprimer3(primers, outfname):
    """Write Primer3 primer objects in ePrimer3 format (Extended)."""
    header = '\n'.join(["# EPRIMER3 PRIMERS %s " % outfname,
                        "#                      Start  Len   Tm     " +
                        "GC%   Sequence"]) + '\n'

    with open(outfname, 'w') as outfh:
        outfh.write(header)
        for idx, primer in enumerate(primers, 1):
            outfh.write("# %s\n" % primer.name)
            outfh.write("%-4d PRODUCT SIZE    %d\n" % (idx, primer.size))
            outfh.write("     FORWARD PRIMER  %-9d  %-3d  %.02f  %.02f  %s\n" %
                        (primer.forward_start, primer.forward_length,
                         primer.forward_tm, primer.forward_gc,
                         primer.forward_seq))
            outfh.write("     REVERSE PRIMER  %-9d  %-3d  %.02f  %.02f  %s\n" %
                        (primer.reverse_start, primer.reverse_length,
                         primer.reverse_tm, primer.reverse_gc,
                         primer.reverse_seq))
            if hasattr(primer, 'internal_start'):
                outfh.write("     INTERNAL OLIGO  " +
                            "%-9d  %-3d  %.02f  %.02f  %s\n" %
                            (primer.internal_start, primer.internal_length,
                             primer.internal_tm, primer.internal_gc,
                             primer.internal_seq))
            outfh.write('\n' * 3)


def __write_primers_json(primers, outfname):
    """Write Primer3 primer objects in JSON format."""
    with open(outfname, 'w') as ofh:
        json.dump(primers, ofh, cls=PrimersEncoder)
