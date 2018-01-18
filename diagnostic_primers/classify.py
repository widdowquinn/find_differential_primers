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

from collections import defaultdict

from Bio.Emboss import PrimerSearch

from .eprimer3 import load_primers
from .primersearch import parse_output


class PDPDiagnosticPrimers(object):

    """Collection of diagnostic primers."""

    def __init__(self, name):
        self.name = str(name)
        self._groups = defaultdict(set)   # Groups with diagnostic primers
        self._primers = dict()

    def add_diagnostic_primer(self, primer, group):
        """Add a diagnostic primer set, specifying the group

        - primer     Eprimer3 object describing the primer set
        - group      the group it's specific to
        """
        self._groups[group].add(primer)
        self._primers[primer.name] = primer

    def diagnostic_primer(self, group):
        """Returns list of primers diagnostic for the passed group"""
        return self._groups[group]

    @property
    def groups(self):
        return sorted(list(self._groups.keys()))


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
    # Parse passed collection and generate dictionary keyed by all groups,
    # with values a set of names of members of those groups
    names = set()                  # set of all genome names
    groups = defaultdict(set)      # group name: set of genome names
    for genome in coll.data:
        for group in genome.groups:
            groups[group].add(genome.name)
            names.add(genome.name)
    groupdata = [(members, name) for (name, members) in
                 groups.items()]

    # Create dictionary to hold primer cross-hybridisation targets keyed by
    # primer name with value a set of all genome targets
    crosshyb = defaultdict(set)

    # Parse the collection and follow the linked primersearch JSON file
    primers = {}
    for genome in coll.data:
        # All primers amplify their own source genome. Load the list
        # of primers and populate the crosshyb dictionary
        for primer in load_primers(genome.primers, fmt="json"):
            crosshyb[primer.name].add(genome.name)
            primers[primer.name] = primer

        # Load the data for primersearch cross-hybridisation, and populate
        # the crosshyb dictionary
        with open(genome.primersearch, 'r') as ifh:
            psdata = json.load(ifh)

            # Each key other than "query" and "primers" is the name of
            # the genome being tested against, and has a PrimerSearch
            # output file
            crosshybnames = [_ for _ in psdata.keys() if _ not in
                             ('primers', 'query')]
            for name in crosshybnames:
                data = parse_output(psdata[name])
                for primer in data:
                    lentest = [max_amplicon > len(amplimer) > min_amplicon for
                               amplimer in primer.amplimers]
                    if sum(lentest):
                        crosshyb[primer.name].add(name)
    crosshybdata = [(targets, primer) for (primer, targets) in
                    crosshyb.items()]

    # To determine group-specific primer sets, we loop through the group
    # data, and compare their members to each of the targets. As we find
    # specific primer sets, we remove them from the pool.
    results = PDPDiagnosticPrimers(coll.name)
    for members, group in groupdata:
        for targets, primer in crosshybdata:
            if members == targets:  # Primers are specific
                results.add_diagnostic_primer(primers[primer], group)

    return results


def write_results(outdir):
    """Writes files describing PDPDiagnosticPrimers object data to outdir

    - outdir     path to directory containing output file results
    """
    raise NotImplementedError
