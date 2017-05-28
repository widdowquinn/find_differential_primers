#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""blast.py

Code for interaction with BLAST+ for diagnostic primer prediction

(c) The James Hutton Institute 2016-2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2016-2017 The James Hutton Institute

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

from Bio.Blast.Applications import NcbiblastnCommandline


def build_commands(collection, blast_db, blast_exe, blast_dir,
                   force=False):
    """Builds and returns a list of BLAST command-lines to run on each
    sequence in the passed GenomeCollection, using the Biopython interface.
    """
    clines = []
    for g in collection.data:
        if blast_dir is None:
            stem = os.path.splitext(g.seqfile)[0]
        else:
            stempath = os.path.split(os.path.splitext(g.seqfile)[0])
            stemdir = os.path.join(*stempath[:-1], blast_dir)
            os.makedirs(stemdir, exist_ok=force)  # Python 3.2+ only
            stem = os.path.join(stemdir, stempath[-1])
        cline = NcbiblastnCommandline(query=g.fastafname,
                                      db=blast_db,
                                      out=stem + '.tab',
                                      task='blastn-short',
                                      max_target_seqs=1,
                                      outfmt=6,
                                      perc_identity=90,
                                      ungapped=True)
        g.cmds['blastscreen'] = cline
        clines.append(cline)
    return clines
