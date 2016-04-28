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
            print(stem)
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
