#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""extract.py

Code to extract amplicons corresponding to primer sets from source genomes

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
import math
import statistics

from collections import defaultdict, namedtuple

from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

from diagnostic_primers import load_primers
from diagnostic_primers.primersearch import parse_output


class PDPAmpliconError(Exception):
    """Custom exception for handling amplicons"""

    def __init__(self, message):
        super(PDPAmpliconError, self).__init__(message)


class PDPAmplicon(object):
    """Data about a primer amplicon

    This is the primer itself, and the amplified sequence from a
    single target genome.
    """

    def __init__(self, name, primer=None, primersearch=None, amplimer=None, seq=None):
        """Initialise object.

        - name       name for the amplicon
        - primer     primer object
        - psresult   primersearch result object (self-amplifiers don't have it)
        - amplimer   amplified region of genome (Biopython Seq)
        """
        self._name = str(name)
        self._primer = None
        self._psresult = None
        self._amplimer = None
        self._seq = None
        self._primer_indexed = None
        self.primer = primer
        self.psresult = primersearch
        self.amplimer = amplimer
        self.seq = seq

    @property
    def name(self):
        return self._name

    @property
    def primer(self):
        return self._primer

    @primer.setter
    def primer(self, val):
        """Primer object for the amplimer"""
        self._primer = val

    @property
    def primersearch(self):
        return self._psresult

    @primersearch.setter
    def primersearch(self, val):
        """PrimerSearch object describing the matches for the primer"""
        self._psresult = val

    @property
    def amplimer(self):
        return self._amplimer

    @amplimer.setter
    def amplimer(self, val):
        """PrimerSearchAmplimer describing the amplimer used"""
        self._amplimer = val

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, val):
        """Seq object describing amplified sequence

        TODO: Assert that the sequence corresponds to the defined primers
        """
        self._seq = val
        if self._seq is not None:
            self._seq.id = self.name
            self._seq.description = "Predicted diagnostic amplicon"

    def __len__(self):
        """Return length of amplified sequence"""
        return len(self.seq)


class PDPAmpliconCollection(object):
    """Collection of PDPAmplicon objects

    Provides methods to operate on a collection of PDPAmplicons.
    """

    def __init__(self, name):
        self._name = str(name)
        self._amplicons = {}  # Amplicons stored, keyed by name
        self._primers = set()

    def new_amplicon(self, name, primer, primersearch, amplimer, seq):
        """Create and return a new PDPAmplicon object

        - name            Identifier for the amplicon (unique)
        - primer          ePrimer3.Primer
        - primersearch    PrimerSearchRecord
        - amplimer        PrimerSearchAmplimer
        - seq             Biopython Seq
        """
        if name in self._amplicons:  # Name must be unique
            raise PDPAmpliconError("New amplicon name must be unique: %s exists" % name)
        self._amplicons[name] = PDPAmplicon(name, primer, primersearch, amplimer, seq)
        self.__index_by_primer()  # Build primer-indexed dict of amplicons
        self._primers.add(primer)
        return self._amplicons[name]

    def get_primer_amplicon_sequences(self, primer_name):
        """Returns a list of amplicon sequences for named primer

        - primer_name       Name of the primer we want sequences for
        """
        return [_.seq for _ in self._primer_indexed[primer_name]]

    def write_amplicon_sequences(self, pname, fname):
        """Write all amplicon sequences as FASTA to passed location

        - pname       Primer for which to write amplicons
        - fname             Path to write FASTA file
        """
        with open(fname, "w") as ofh:
            seqdata = self.get_primer_amplicon_sequences(pname)
            # Order sequence data for consistent output (aids testing)
            seqdata = [_[1] for _ in sorted([(seq.id, seq) for seq in seqdata])]
            SeqIO.write(seqdata, ofh, "fasta")

    def __iter__(self):
        """Iterate over amplicons in the collection"""
        for _ in self._amplicons.values():
            yield _

    def __getitem__(self, key):
        """Return an amplicon by name"""
        return self._amplicons[key]

    def __len__(self):
        """Return number of amplicons in the collection"""
        return len(self._amplicons)

    def __index_by_primer(self):
        """Create internal index of amplicons keyed by primer name"""
        self._primer_indexed = defaultdict(set)
        for _ in self._amplicons.values():
            self._primer_indexed[_.primer.name].add(_)

    @property
    def names(self):
        """List of names of amplicons"""
        return list(self._amplicons.keys())

    @property
    def primers(self):
        """List of primer sets in the collection"""
        return list(self._primers)

    @property
    def primer_names(self):
        """List of primer names in the collection"""
        return [_.name for _ in self._primers]

    @property
    def primer_amplicons(self):
        """set of amplicons for a named primer"""
        return self._primer_indexed


def extract_amplicons(
    name, primer, pdpcoll, min_amplicon, max_amplicon, seq_cache=None
):
    """Return PDPAmpliconCollection corresponding to primers in the passed file

    - name        identifier for this action
    - primer      Primer3.Primers object
    - pdpcoll     PDPCollection containing information about the primer
                  and target genome sources (primersearch, seqfile,
                  filestem)
    - seq_cache   a directory of potential target genomes, cached locally to
                  save on file IO
    """
    # Make dictionaries of each config entry and source genome path by name
    namedict = {_.name: _ for _ in pdpcoll.data}
    genomepaths = {_.name: _.seqfile for _ in pdpcoll.data}

    # For each primer, we identify the corresponding config file entry
    # We're going to extract each primer amplicon individually. We know
    # at this point that the primer is predicted specific only for
    # targets of a particular class, so we can scrape all the relevant
    # primersearch files. We cache those files as we see them, to save on
    # file IO, in a dictionary keyed by primersearch output filename.
    #
    # We also cache genome data, and this is returned by the function,
    # for reuse and reduced fileIO.
    #
    # The complete primer list for each source genome is also cached.
    psoutput_cache = {}  # primersearch output
    if seq_cache is None:  # genomes (may be reused)
        seq_cache = {}
    sourceprimer_cache = {}  # source genome primers
    amplicons = PDPAmpliconCollection(name)

    stem = primer.name.split("_primer_")[0]
    source_data = namedict[primer.sourcename]  # the source sequence for the primers
    seq_cache[source_data.name] = SeqIO.read(source_data.seqfile, "fasta")

    # Cache the source genome primer information
    if stem not in sourceprimer_cache:
        sourceprimer_cache[stem] = {
            _.name: _ for _ in load_primers(source_data.primers, fmt="json")
        }

    # Get primersearch output for each primer as we encounter it
    with open(source_data.primersearch, "r") as ifh:
        psdata = json.load(ifh)
        targets = [_ for _ in psdata.keys() if _ not in ("primers", "query")]

        # Examine each target for the primer
        for target in targets:

            # Cache primersearch output for the target
            if psdata[target] not in psoutput_cache:
                psoutput_cache[psdata[target]] = {
                    _.name: _ for _ in parse_output(psdata[target], genomepaths[target])
                }
            # TODO: turn output from parse_output into an indexable object,
            #       so we can use primer names to get results, rather than
            #       hacking that, as we do above.

            # Cache genome data for the target (read each file once, saves time)
            if target not in seq_cache:
                seq_cache[target] = SeqIO.read(namedict[target].seqfile, "fasta")
            target_genome = seq_cache[target]

            # psresult holds the primersearch result - we create an
            # amplicon for each amplimer in the psresult
            # If this primer isn't in the set that amplifies the target,
            # continue to the next target
            if primer.name not in psoutput_cache[psdata[target]]:
                continue
            psresult = psoutput_cache[psdata[target]][primer.name]
            for ampidx, amplimer in enumerate(psresult.amplimers):
                coords = (amplimer.forward_start, amplimer.reverse_end)
                # Extract the genome sequence
                # We have to account here for forward/reverse primer
                # amplification wrt target genome sequence. We want
                # all the amplimer sequences to be identically-stranded
                # for downstream alignments so, if the forward/reverse
                # primer sequences don't match between the primer sets and
                # the PrimerSearch results, we flip the sequence here.
                seq = target_genome[min(coords) - 1 : max(coords)]
                if primer.forward_seq != amplimer.forward_seq:
                    seq = seq.reverse_complement()
                if max_amplicon > len(seq) > min_amplicon:
                    amplicons.new_amplicon(
                        "_".join([primer.name, target, str(ampidx + 1)]),
                        primer,
                        psresult,
                        amplimer,
                        seq,
                    )

    return amplicons, seq_cache


# Results object for returning distance calculations
DistanceResults = namedtuple(
    "DistanceResults",
    "matrix distances mean sd min max unique nonunique shannon evenness",
)


def calculate_distance(aln, calculator="identity"):
    """Report distance measures for the passed nucleotide AlignIO object

    - aln           A Bio.AlignIO.MultipleSeqAlignment
    - calculator    The metric to use when calculating distance

    If the alignment contains only a single sequence, return zeroes for all fields,
    and a None for the distance matrix
    """
    if len(aln) == 1:  # We can't calculate distances or a matrix
        raise PDPAmpliconError(
            "Alignment contains a single sequence: cannot calculate distances"
        )
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(aln)
    # Flatten the DistanceMatrix's matrix, discarding the diagonal
    # The DistanceMatrix is a lower-triangular matrix, so the last item
    # in each row is on the diagonal
    # This gives a list of all pairwise distances
    distances = [_ for sublist in dm.matrix for _ in sublist[:-1]]
    # The number of unique amplicons is found by taking the length
    # of the set comprehension of sequences in the alignment
    unique = len({str(_.seq) for _ in aln})
    nonunique = len(aln) - unique
    shannon, evenness = shannon_index(aln)
    # The standard deviation calculation throws an error if there's only one
    # distance, which occurs when there are only two sequences in the
    # alignment
    return DistanceResults(
        dm,
        distances,
        statistics.mean(distances),
        statistics.stdev(distances) if len(aln) > 2 else 0,
        min(distances),
        max(distances),
        unique,
        nonunique,
        shannon,
        evenness,
    )


def shannon_index(aln):
    """Returns the Shannon index and evenness of a sequence alignment

    aln         A Bio.AlignIO.MultipleSeqAlignment

    Shannon index is a measure of the sequence diversity of an alignment. It accounts
    for abundance and evenness of the number of unique sequences present.

    Suppose our alignment is composed of N sequences, but only S distinct sequences.
    For each distinct sequence i, let the proportion of sequences N that are sequence i
    be p_i. Then the Shannon index H is:

    $H = - \\sum_i^S p_i \\ln p_i$



    The evenness of the distribution can then be calculated as

    $E_H = H / \\ln(S)$

    E_H is constrained between 1 (completely even) and 0 (uneven distribution).
    """
    # Make dictionary of unique sequence counts
    countdict = defaultdict(int)
    alnsize = len(aln)
    for seq in aln:
        countdict[hash(str(seq.seq))] += 1

    # Calculate Shannon index
    sindex = 0
    for _, count in countdict.items():
        p_i = count / alnsize
        sindex -= p_i * math.log(p_i)

    # Calculate evenness
    try:
        seven = sindex / math.log(len(countdict))
    except ZeroDivisionError:
        seven = 0

    return sindex, seven
