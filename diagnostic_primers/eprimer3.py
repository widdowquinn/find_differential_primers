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

import errno
import os

from Bio.Emboss.Applications import Primer3Commandline
from Bio.Emboss import Primer3


def build_commands(collection, eprimer3_exe, eprimer3_dir=None, force=False,
                   argdict=None):
    """Builds and returns a list of command-lines to run ePrimer3 on each
    sequence in the passed GenomeCollection, using the Biopython interface.
    """
    clines = []
    for g in collection.data:
        if eprimer3_dir is None:
            stem = os.path.splitext(g.seqfile)[0]
        else:
            stempath = os.path.split(os.path.splitext(g.seqfile)[0])
            stemdir = os.path.join(*stempath[:-1], eprimer3_dir)
            os.makedirs(stemdir, exist_ok=force)  # Python 3.2+ only
            stem = os.path.join(stemdir, stempath[-1])
        cline = Primer3Commandline(cmd=eprimer3_exe)
        cline.sequence = g.seqfile
        cline.auto = True
        cline.outfile = stem + '.eprimer3'
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
        g.cmds['ePrimer3'] = cline
        clines.append(cline)
    return clines


def eprimer3_to_fasta(infname, outfname):
    """Converts ePrimer3 format files to FASTA, with one sequence per primer
    oligo/internal oligo.
    """
    with open(infname, 'r') as infh:
        primers = Primer3.read(infh).primers

    seqrecords = []
    for primer in primers:
        print(vars(primer))


def add_primer_names(infname):
    """Adds a name to each of the ePrimer3 primer sets in the passed file.
    Writes a new file (with the suffix '_named'), and returns the filename.
    """
    with open(infname, 'r') as primerfh:
        primers = Primer3.read(primerfh).primers

    # Add a name to each primer set, based on input filename
    for idx, primer in enumerate(primers, 1):
        stem = os.path.splitext(os.path.split(infname)[-1])[0]
        primer.name = "%s_primer_%05d" % (stem, idx)

    # Write named primers to output file
    outfname = os.path.splitext(infname)[0] + '_named.eprimer3'
    write_eprimer3(primers, outfname)
    return outfname


def write_eprimer3(primers, outfname):
    """Write the Primer3 primer objects to the named file, in
    Primer3-compatible form.
    """
    header = '\n'.join(["# EPRIMER3 PRIMERS %s " % outfname,
                        "#                      Start  Len   Tm     " +
                        "GC%   Sequence"]) + '\n'

    with open(outfname, 'w') as outfh:
        outfh.write(header)
        for idx, primer in enumerate(primers, 1):
            outfh.write("# %s\n" % primer.name)
            outfh.write("%-4d PRODUCT SIZE: %d\n" % (idx, primer.size))
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
