#!/usr/bin/env python
#
# fasta_summary.py
#
# Quick script that reports the number of sequences in a multiple FASTA
# file, and the total length of those sequences
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

from Bio import SeqIO
import sys

infilenames = sys.argv[1:]

for infilename in infilenames:
    data = [len(s) for s in SeqIO.parse(infilename, 'fasta')]
    largedata = [d for d in data if d > 1000]
    median = sorted(data)[len(data) / 2]
    sortedlist = sorted(data, reverse=True)
    tot = 0
    while tot < sum(data) / 2:
        n50 = sortedlist.pop()
        tot += n50
    sys.stdout.write("Filename: %s" % infilename)
    sys.stdout.write("Number of sequences: %d" % len(data))
    sys.stdout.write("Median: %d" % median)
    sys.stdout.write("N50 (%d): %d" % (tot, n50))
    sys.stdout.write("Total sequence length: %d" % sum(data))
    if len(largedata):
        sys.stdout.write("Number of large (>1000bp) sequences: %d" %
                         len(largedata))
        largemedian = sorted(largedata)[len(largedata) / 2]
        lsortedlist = sorted(largedata, reverse=True)
        largetot = 0
        while largetot < sum(data) / 2:
            try:
                largen50 = lsortedlist.pop()
            except IndexError:              # Large contigs don't reach N50
                sys.stdout.write("Large contigs don't cover 50% of " +
                                 "assembled sequence")
                largen50 = 0
                break
            largetot += largen50
        sys.stdout.write("Median (>1000bp): %d" % largemedian)
        sys.stdout.write("N50 (>1000bp) (%d): %d" % (largetot, largen50))
    else:
        sys.stdout.write("No large (>1kbp) contigs!")

    print
