#!/usr/bin/env python3
#
# prodigal.py
#
# Code to conduct CDS prediction with Prodigal
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


def build_commands(collection, prodigal_exe, prodigal_dir=None):
    """Builds and returns a list of prodigal command-lines

    The returned commands will run Prodigal on each sequence in the passed
    PDPCollection.

    If no output directory is provided, output will be placed in the same
    directory as the input sequence files.
    """
    clines = []  # Holds command-lines

    # Create the output directory, if needed
    os.makedirs(prodigal_dir, exist_ok=True)

    for g in collection.data:
        if prodigal_dir is None:
            stem = os.path.splitext(g.seqfile)[0]
        else:
            stempath = os.path.split(os.path.splitext(g.seqfile)[0])
            stem = os.path.join(prodigal_dir, stempath[-1])
        ftfile = stem + '.features'
        outfile = stem + '.gff'
        cline = ' \\\n          '.join([prodigal_exe,
                                        '-m -c -f gff',
                                        '-a %s' % ftfile,
                                        '-i %s' % g.seqfile,
                                        '-o %s' % outfile])
        g.cmds['prodigal'] = cline
        clines.append(cline)
    return clines
