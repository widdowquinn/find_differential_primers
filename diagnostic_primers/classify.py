#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""classify.py

Code to classify primer sets by predicted specificity

(c) The James Hutton Institute 2018

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

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

Copyright (c) 2018 The James Hutton Institute
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

import json


def classify_primers(coll, min_amplicon=50, max_amplicon=300):
    """Classifies each of the primer sets referred to in the passed collection

    - coll      PDPCollection descibing the genomes in the run, with links
                to the predicted primer sets, and their matches to other
                genomes post-primersearch
    - min_amplicon   The minimum length of an amplicon that could be
                     considered a false positive
    - max_amplicon   The maximum length of an amplicon that could be
                     considered a false positive

    First, the collection data is parsed, and a dictionary generated, keyed
    by group, with values being a set of genome names belonging to the
    group.

    For each genome in the passed PDPCollection, the path to corresponding
    PrimerSearch JSON file is followed, and the file parsed to obtain the
    path to the relevant PrimerSearch output.

    Each PrimerSearch output file is parsed and a dictionary populated.
    The dict is keyed by primer ID, with value being a set of the names
    of genomes where an amplicon is theoretically produced (filtered for
    amplicon length).

    Using set logic, each primer's set of potential amplified genomes is
    compared against the sets of specific groups to determine if there is
    an exact match.

    The primers that amplify exactly those genomes which are members of one
    of the defined classes are returned as a PDPDiagnosticPrimers object that
    is a collection of Primer3.Primers objects.
    """
    raise NotImplementedError
