#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""eprimer3.py

Code to conduct primer prediction with ePrimer3

(c) The James Hutton Institute 2016-2019
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

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

Copyright (c) 2016-2019 The James Hutton Institute

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

from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline


def build_commands(collection, eprimer3_exe, eprimer3_dir, argdict=None):
    """Builds and returns a list of command-lines to run ePrimer3

    The commands will run on each sequence in the passed PDPCollection.
    """
    clines = []  # Holds command-lines

    # Ensure output directory exists
    os.makedirs(eprimer3_dir, exist_ok=True)

    for g in collection.data:
        stempath = os.path.split(os.path.splitext(g.seqfile)[0])
        stem = os.path.join(eprimer3_dir, stempath[-1])
        # Are we using the filter region information for this design?
        # Two ways we won't use the filter region: the --filter argument/
        # ep_filter argument is not set; the --filter argument/ep_filter
        # argument *is* set, but there's no filtered genome sequence file.
        if argdict["ep_filter"] and g.filtered_seqfile is not None:
            seqfile = g.filtered_seqfile
        else:
            seqfile = g.seqfile
        cline = build_command(eprimer3_exe, seqfile, stem, argdict)
        g.cmds["ePrimer3"] = cline
        clines.append(cline)
    return clines


def build_command(eprimer3_exe, seqfile, filestem, argdict=None):
    """Builds and returns ePrimer3 command line.

    The ePrimer3 command uses the Biopython interface
    """
    cline = Primer3Commandline(cmd=eprimer3_exe)
    cline.sequence = seqfile
    cline.auto = True
    cline.outfile = filestem + ".eprimer3"
    if argdict is not None:
        prange = [0, 200]
        args = [(a[3:], v) for a, v in argdict.items() if a.startswith("ep_")]
        for arg, val in args:
            if "psizemin" == arg:
                prange[0] = val
            elif "psizemax" == arg:
                prange[1] = val
            else:
                setattr(cline, arg, val)
    setattr(cline, "prange", "%d-%d" % tuple(prange))
    return cline
