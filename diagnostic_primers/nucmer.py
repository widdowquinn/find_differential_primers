#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""nucmer.py

Provides helper functions to run nucmer from pdp

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

import os

from collections import namedtuple
from itertools import permutations

from diagnostic_primers import PDPException
from diagnostic_primers.sge_jobs import Job


# convenience factories for nucmer results
# Information useful for context of nucmer comparison output
NucmerOutput = namedtuple(
    "NucmerOutput", "query subject out_delta out_filter cmd_nucmer cmd_delta"
)
# Parsed nucmer .delta file data
NucmerDelta = namedtuple(
    "NucmerDelta", "deltafile queryfile subjectfile query_intervals"
)


# Exceptions for nucmer processing
class PDPNucmerException(PDPException):
    """Exception thrown during nucmer run or processing"""

    def __init__(self, msg="Error in `nucmer` processing"):
        PDPException.__init__(self, msg)


# Class to describe a nucmer command, for consistency with Biopython commands
class NucmerCommand(object):
    """Command-line for nucmer"""

    def __init__(self, cline, infile, outfile):
        self.cline = cline
        self.infile = infile
        self.outfile = outfile

    def __str__(self):
        return " ".join(self.cline)


# Class to describe a delta-filter command, for consistency with Biopython commands
class DeltaFilterCommand(object):
    """Command-line for delta-filter"""

    def __init__(self, cline, infile, outfile):
        self.cline = cline
        self.infile = infile
        self.outfile = outfile

    def __str__(self):
        return " ".join(self.cline)


class DeltaData(object):
    """Class to hold MUMmer/nucmer output "delta" data

    This is required because the ordering within files differs depending on MUMmer
    build, for the same version (as evidenced by differences between OSX and Linux
    builds), and a means of testing for equality of outputs is necessary.

    The output file structure and format is described at
    http://mummer.sourceforge.net/manual/#nucmeroutput

    Each file is represented as:

    - header: first line of the .delta file, naming the two input comparison files; stored
              as a tuple (path1, path2), returned as the combined string; the individual
              files are stored as self._query and self._subject
    - program: name of the MUMmer program that produced the output
    - query: path to the query sequence file
    - subject: path to the subject sequence file
    """

    def __init__(self, name, handle=None):
        self.name = name
        self._metadata = None
        self._comparisons = []
        if handle is not None:
            self.from_delta(handle)

    def from_delta(self, handle):
        """Populate the object from the passed .delta or .filter filehandle"""
        parser = DeltaIterator(handle)
        for element in parser:
            if isinstance(element, DeltaMetadata):
                self._metadata = element
            if isinstance(element, DeltaComparison):
                self._comparisons.append(element)

    @property
    def comparisons(self):
        """Comparisons in the .delta file"""
        return self._comparisons

    @property
    def metadata(self):
        """Metadata from the .delta file"""
        return self._metadata

    @property
    def reference(self):
        """The reference file for the MUMmer comparison"""
        return self._metadata.reference

    @property
    def program(self):
        """The MUMmer program used for the comparison"""
        return self._metadata.program

    @property
    def query(self):
        """The query file for the MUMmer comparison"""
        return self._metadata.query

    def __eq__(self, other):
        # We do not enforce equality of metadata, as the full path to both query and reference is
        # written in the .delta file, and we care only about the alignment data, and the program
        # that was used.
        if not isinstance(other, DeltaData):
            return False
        return (self.program == other.program) and (
            self._comparisons == other._comparisons
        )

    def __len__(self):
        return len(self._comparisons)

    def __str__(self):
        """Return the object in .delta format output"""
        outstr = os.linesep.join(
            [str(self._metadata)] + [str(_) for _ in self._comparisons]
        )
        return outstr


class DeltaMetadata(object):
    """Represents the metadata header for a MUMmer .delta file"""

    def __init__(self):
        self.reference = None
        self.query = None
        self.program = None

    def __eq__(self, other):
        if not isinstance(other, DeltaMetadata):
            return False
        return (self.reference, self.query, self.program) == (
            other.reference,
            other.query,
            other.program,
        )

    def __str__(self):
        return "{} {}{}{}".format(self.reference, self.query, os.linesep, self.program)


class DeltaComparison(object):
    """Represents a comparison between two sequences in a .delta file"""

    def __init__(self, header, alignments):
        self.header = header
        self.alignments = alignments

    def add_alignment(self, aln):
        """Add passed alignment to this object

        :param aln:  DeltaAlignment object
        """
        self.alignments.append(aln)

    def __eq__(self, other):
        if not isinstance(other, DeltaComparison):
            return False
        return (self.header == other.header) and (
            sorted(self.alignments) == sorted(other.alignments)
        )

    def __len__(self):
        return len(self.alignments)

    def __str__(self):
        outstr = os.linesep.join([str(self.header)] + [str(_) for _ in self.alignments])
        return outstr


class DeltaHeader(object):
    """Represents a single sequence comparison header from a MUMmer .delta file"""

    def __init__(self, reference, query, reflen, querylen):
        self.reference = reference
        self.query = query
        self.referencelen = int(reflen)
        self.querylen = int(querylen)

    def __eq__(self, other):
        if not isinstance(other, DeltaHeader):
            return False
        return (self.reference, self.query, self.referencelen, self.querylen) == (
            other.reference,
            other.query,
            other.referencelen,
            other.querylen,
        )

    def __str__(self):
        return ">{} {} {} {}".format(
            self.reference, self.query, self.referencelen, self.querylen
        )


class DeltaAlignment(object):
    """Represents a single alignment region and scores for a pairwise comparison"""

    def __init__(self, refstart, refend, qrystart, qryend, errs, simerrs, stops):
        self.refstart = int(refstart)
        self.refend = int(refend)
        self.querystart = int(qrystart)
        self.queryend = int(qryend)
        self.errs = int(errs)
        self.simerrs = int(simerrs)
        self.stops = int(stops)
        self.indels = []

    def __lt__(self, other):
        return (self.refstart, self.refend, self.querystart, self.queryend) < (
            other.refstart,
            other.refend,
            other.querystart,
            other.queryend,
        )

    def __eq__(self, other):
        return (self.refstart, self.refend, self.querystart, self.queryend) == (
            other.refstart,
            other.refend,
            other.querystart,
            other.queryend,
        )

    def __str__(self):
        outstr = [
            "{} {} {} {} {} {} {}".format(
                self.refstart,
                self.refend,
                self.querystart,
                self.queryend,
                self.errs,
                self.simerrs,
                self.stops,
            )
        ] + [str(_) for _ in self.indels]
        return os.linesep.join(outstr)


class DeltaIterator(object):
    """Iterator for MUMmer .delta files. Returns a stream of DeltaMetadata,
    DeltaComparison and DeltaAlignment objects when iterated over a filehandle

    The .delta file structure and format is described at
    http://mummer.sourceforge.net/manual/#nucmeroutput
    """

    def __init__(self, handle):
        """Instantiate the class with the passed filehandle"""
        self._handle = handle
        self._metadata = None  # metadata for a .delta file
        self._header = None  # header information for a pairwise comparison
        self._comparison = None  # current comparison region

    def __iter__(self):
        """Iterate over elements of the .delta file as DeltaHeader and DeltaAlignment objects"""
        return iter(self.__next__, None)

    def __next__(self):
        """Parse the next element from the .delta file"""
        # Parse .delta file metadata
        if self._metadata is None:
            self._metadata = DeltaMetadata()
            self._metadata.reference, self._metadata.query = (
                self._handle.readline().strip().split()
            )
            self._metadata.program = self._handle.readline().strip()
            return self._metadata

        # Parse remaining lines into a DeltaHeader for each comparison, and corresponding
        # DeltaAlignments
        line = self._handle.readline()
        while line:
            # If we're at the start of an alignment, create a new DeltaAlignment
            if line.startswith(">"):
                if self._comparison is not None:
                    return self._comparison
                self._header = DeltaHeader(*(line[1:].split()))
                self._comparison = DeltaComparison(self._header, [])
            # Populate the current pairwise alignment with each individual alignment
            else:
                alndata = line.rstrip().split()
                if len(alndata) > 1:  # alignment header
                    alignment = DeltaAlignment(*alndata)
                elif alndata[0] == "0":
                    alignment.indels.append(alndata[0])
                    self._comparison.add_alignment(alignment)
                else:
                    alignment.indels.append(alndata[0])
            # Get the next line and return the final comparison if we're at the end of file
            line = self._handle.readline()
            if not line:
                return self._comparison


# Generate list of Job objects, one per NUCmer run
def generate_nucmer_jobs(
    groupdata, outdir, nucmer_exe, filter_exe, maxmatch=False, jobprefix="PDPNUCmer"
):
    """Return list of Jobs describing NUCmer command-lines for PDP.

    - groupdata - iterable of PDPData objects
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files, generating Jobs describing NUCmer command lines
    for each pairwise comparison.
    """
    nucmerdata = generate_nucmer_commands(
        groupdata, outdir, nucmer_exe, filter_exe, maxmatch
    )
    joblist = []
    for idx, ndata in enumerate(nucmerdata):
        njob = Job("%s_%06d-n" % (jobprefix, idx), ndata.cmd_nucmer)  # nucmer job
        fjob = Job("%s_%06d-f" % (jobprefix, idx), ndata.cmd_delta)  # filter job
        fjob.add_dependency(njob)
        joblist.append((fjob, ndata))
    return joblist


# Generate list of NUCmer pairwise comparison command lines from
# passed sequence filenames
def generate_nucmer_commands(groupdata, outdir, nucmer_exe, filter_exe, maxmatch=False):
    """Return list of NUCmer command-lines for PDP.

    The first element returned is a list of NUCmer commands, and the
    second a corresponding list of delta_filter_wrapper.py commands.
    The NUCmer commands should each be run before the corresponding
    delta-filter command.

    - groupdata - iterable of PDPData objects
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files generating NUCmer command lines for each
    pairwise comparison.
    """
    ndata = []
    comparisons = permutations(groupdata, 2)
    for (query, subject) in comparisons:
        nucmerdata = construct_nucmer_cmdline(
            query, subject, outdir, nucmer_exe, filter_exe, maxmatch
        )
        ndata.append(nucmerdata)
    return ndata


# Generate single NUCmer pairwise comparison command line from pair of
# input filenames
def construct_nucmer_cmdline(
    query, subject, outdir, nucmer_exe, filter_exe, maxmatch=False
):
    """Return a tuple of corresponding NUCmer and delta-filter commands

    The split into a tuple was made necessary by changes to SGE/OGE.
    The delta-filter command must now be run as a dependency of the NUCmer
    command, and be wrapped in a Python script to capture STDOUT.

    NOTE: This command-line writes output data to a subdirectory of the passed
    outdir, called "nucmer_output".

    - query - query PDPData object
    - subject - subject PDPData object
    - outdir - path to output directory
    - maxmatch - Boolean flag indicating whether to use NUCmer's -maxmatch
    option. If not, the -mum option is used instead
    """
    outprefix = os.path.join(
        outdir,
        "%s_vs_%s"
        % (
            os.path.splitext(os.path.split(query.seqfile)[-1])[0],
            os.path.splitext(os.path.split(subject.seqfile)[-1])[0],
        ),
    )
    if maxmatch:
        mode = "--maxmatch"
    else:
        mode = "--mum"

    # output filepaths
    out_delta = outprefix + ".delta"
    out_filter = outprefix + ".filter"

    # command-line objects to return
    nucmercmd = NucmerCommand(
        [nucmer_exe, mode, "-p", outprefix, query.seqfile, subject.seqfile],
        query.seqfile,
        out_delta,
    )
    filtercmd = DeltaFilterCommand(
        ["delta_filter_wrapper.py ", filter_exe, "-1", out_delta, out_filter],
        out_delta,
        out_filter,
    )

    return NucmerOutput(
        query=query,
        subject=subject,
        out_delta=out_delta,
        out_filter=out_filter,
        cmd_nucmer=nucmercmd,
        cmd_delta=filtercmd,
    )


def parse_delta_query_regions(fname, min_sim_errors=0, min_err_rate=0):
    """Return NucmerDelta object nucmer aligned regions on query

    fname             path to the input .delta file
    min_sim_errors    skip aligned regions that contain fewer than the passed
                      count of similarity errors
    min_err_rate      skip aligned regions that have a similarity error rate
                      per base less than the passed value

    Extracts a list of intervals from a nucmer .delta output file. The regions
    are defined as (start, end, sim_error) tuples, where start and end refer to
    start locations on the query genome, and sim_error refers to the count of
    similarity errors.

    Returns a NucmerDelta object that includes the path to the .delta file, and
    paths to the query and subject sequence files

    The query filestem is used as the identifier for the interval, enabling
    interval calculations with the returned values
    """
    intervals = []
    with open(fname, "r") as dfh:
        # First line is paths to query and subject files
        qpath, spath = dfh.readline().strip().split(" ")
        stem = os.path.splitext(os.path.split(qpath)[-1])[0]  # query filestem
        # Remaining lines are either:
        # - NUCMER (announcing which alignment package was used)
        # - >desc1 desc2 len1 len2   description lines and lengths of input files
        # - seven integers (start/end positions, errors, and stops)
        # - signed integer (symbols between indels)
        # We loop over these, parsing only the lines with seven integers
        for line in [_.strip().split() for _ in dfh.readlines()]:
            if len(line) == 7:
                start, end, sim_errors = (int(line[0]), int(line[1]), int(line[5]))
                if sim_errors < min_sim_errors or (
                    sim_errors / (end - start) < min_err_rate
                ):
                    continue
                intervals.append((stem, start, end))
    return NucmerDelta(
        deltafile=fname, queryfile=qpath, subjectfile=spath, query_intervals=intervals
    )
