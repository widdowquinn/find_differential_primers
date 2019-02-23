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
DD2 5DA,
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
import os

from collections import defaultdict

from Bio.Emboss.Primer3 import Primers

from diagnostic_primers import load_primers, PrimersEncoder, write_primers
from diagnostic_primers.primersearch import parse_output


class PDPDiagnosticPrimersEncoder(json.JSONEncoder):

    """JSON encoder for PDPDiagnosticPrimers objects"""

    def default(self, obj):
        if isinstance(obj, Primers):
            encoder = PrimersEncoder()
            return encoder.default(obj)
        if not isinstance(obj, PDPDiagnosticPrimers):
            return json.JSONEncoder.default(self, obj)

        # Convert PDPDiagnosticPrimers object to serialisable dictionary
        # and return
        return obj.__dict__


class PDPDiagnosticPrimers(object):

    """Collection of diagnostic primers."""

    def __init__(self, name):
        self.name = str(name)
        self._groups = defaultdict(list)  # Groups with diagnostic primers
        self._primers = dict()

    def add_diagnostic_primer(self, primer, group):
        """Add a diagnostic primer set, specifying the group

        - primer     Eprimer3 object describing the primer set
        - group      the group it's specific to
        """
        self._groups[group].append(primer)
        self._primers[primer.name] = primer

    def diagnostic_primer(self, group):
        """Returns list of primers diagnostic for the passed group"""
        return self._groups[group]

    @property
    def groups(self):
        return sorted(list(self._groups.keys()))

    @property
    def primers(self):
        return sorted(list(self._primers.keys()))


def classify_primers(coll, min_amplicon, max_amplicon):
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

    For each genome in the passed PDPCollection, the path to the
    PrimerSearch JSON file corresponding to the primers generated from that
    genome is followed, and the file parsed to obtain the
    path to the relevant PrimerSearch output.

    Each PrimerSearch output file is parsed and a dictionary
    populated:
    - dictionary keyed by primer ID, with value being a set of the names
    of genomes where an amplicon is theoretically produced (filtered for
    amplicon length).

    Using set logic, each primer's set of potential amplified genomes is
    compared against the sets of specific groups to determine if there is
    an exact match.

    The primers that amplify exactly those genomes which are members of one
    of the defined classes are returned as a PDPDiagnosticPrimers object that
    is a collection of Primer3.Primers objects.
    """
    # Parse passed collection and generate local dictionary with a key for
    # each group in the collection. This will take with values that are a
    # set of names of members of those groups.
    names = set()  # set of all genome names
    groups = defaultdict(set)  # group name: set of genome names
    for genome in coll.data:
        for group in genome.groups:
            groups[group].add(genome.name)
            names.add(genome.name)
    groupdata = [(members, name) for (name, members) in groups.items()]

    # Create dictionary to hold primer cross-hybridisation targets
    # - keyed by primer name with value a set of all genome targets
    crosshyb = defaultdict(set)

    # Parse the collection and follow the linked primersearch JSON file
    # - populate a dictionary of paths to each input genome, keyed by name
    genomepaths = {}  # dictionary of paths to each genome file, keyed by name
    for item in coll.data:
        genomepaths[item.name] = item.seqfile
    primers = {}
    for genome in coll.data:
        # All primers amplify their own source genome. Load the list
        # of primers and populate the crosshyb dictionaries
        for primer in load_primers(genome.primers, fmt="json"):
            crosshyb[primer.name].add(genome.name)
            primers[primer.name] = primer

        # Load the data for primersearch cross-hybridisation, and populate
        # the crosshyb dictionary
        with open(genome.primersearch, "r") as ifh:
            psdata = json.load(ifh)

            # Each key other than "query" and "primers" is the name of
            # the genome being tested against, and has a PrimerSearch
            # output file
            # crosshybnames is a list of target genome names
            crosshybnames = [_ for _ in psdata.keys() if _ not in ("primers", "query")]
            for name in crosshybnames:
                data = parse_output(
                    psdata[name], genomepaths[name]
                )  # primersearch results
                for primer in data:
                    for amplimer in primer.amplimers:
                        if max_amplicon > len(amplimer) > min_amplicon:
                            crosshyb[primer.name].add(name)
    crosshybdata = [(targets, primer) for (primer, targets) in crosshyb.items()]

    # To determine group-specific primer sets, we loop through the group
    # data, and compare their members to each of the targets. As we find
    # specific primer sets, we remove them from the pool.
    results = PDPDiagnosticPrimers(coll.name)
    for members, group in groupdata:
        for targets, primer in crosshybdata:
            if members == targets:  # Primers are specific
                results.add_diagnostic_primer(primers[primer], group)

    return results


def write_results(results, outfilename, fmt="json"):
    """Writes files describing PDPDiagnosticPrimers object data to outdir

    - results       PDPDiagnosticPrimers object
    - outfilename   Path to output file/directory
    - fmt           Result format

    Several output formats may be written:

    json
    ====
    JSON representation of the complete PDPDiagnosticPrimers object

    summary
    =======
    writes primers, and also tab-separated plain text table with the columns:

    Group          (specifically-amplified group)
    NumPrimers     (number of primer sets specifically amplifying that group)
    Primers        (path to file describing specific primers)

    primers
    =======
    ePrimer3 and JSON format files describing specific primers, one per group
    """
    funcs = {
        "json": __write_results_json,
        "summary": __write_results_summary,
        "primers": __write_results_primers,
    }

    funcs[fmt](results, outfilename)


def __write_results_json(results, outfilename):
    """Write PDPDiagnosticPrimers JSON representation

    - results       PDPDiagnosticPrimers object
    - outfilename   path to output file
    """
    with open(outfilename, "w") as ofh:
        json.dump(results, ofh, cls=PDPDiagnosticPrimersEncoder)


def __write_results_primers(results, outdir):
    """Write JSON/ePrimer3 files describing diagnostic primers

    - results       PDPDiagnosticPrimers object
    - outdir        path to directory for output files
    """
    for group in results.groups:
        outstem = os.path.join(outdir, "%s_primers" % group)
        write_primers(results.diagnostic_primer(group), outstem + ".json", "json")
        write_primers(results.diagnostic_primer(group), outstem + ".ePrimer3", "ep3")
        # We don't write .bed files using write_primers(), as this reports the primer
        # sets in the context of the genomes they're defined from. However, for
        # classify (and extract), we want locations on each genome they amplify from
        # write_primers(results.diagnostic_primer(group), outstem + ".bed", "bed")


def __write_results_summary(results, outfilename):
    """Write summary table of diagnostic primer results"""
    # Write the diagnostic primer sets first
    outdir = os.path.split(outfilename)[0]
    __write_results_primers(results, outdir)

    # Write the summary table
    outstr = ["\t".join(["Group", "NumPrimers", "Primers"])]
    for group in results.groups:
        outstr.append(
            "\t".join(
                [
                    group,
                    str(len(results.diagnostic_primer(group))),
                    os.path.join(outdir, "%s_primers.json" % group),
                ]
            )
        )
    with open(outfilename, "w") as ofh:
        ofh.write("\n".join(outstr) + "\n")
