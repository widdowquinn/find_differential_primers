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

import os

from collections import namedtuple
from itertools import permutations

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
        njob = Job("%s_%06d-n" % (jobprefix, idx), ndata.cmd_nucmer)
        fjob = Job("%s_%06d-f" % (jobprefix, idx), ndata.cmd_delta)
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
    nucmercmd = [nucmer_exe, mode, "-p", outprefix, query.seqfile, subject.seqfile]
    out_delta = outprefix + ".delta"
    out_filter = outprefix + ".filter"
    filtercmd = ["delta_filter_wrapper.py ", filter_exe, "-1", out_delta, out_filter]
    return NucmerOutput(
        query=query,
        subject=subject,
        out_delta=out_delta,
        out_filter=out_filter,
        cmd_nucmer=nucmercmd,
        cmd_delta=filtercmd,
    )


def parse_delta(fname, no_exact=False):
    """Return a NucmerDelta object describing nucmer results in fname

    fname         path to the input .delta file
    no_exact      skip aligned regions that contain no similarity errors

    Extracts a list of intervals from a nucmer .delta output file. The regions
    are defined as (start, end, sim_error) tuples, where start and end refer to
    start locations on the query genome, and sim_error refers to the count of
    similarity errors.

    Returns a NucmerDelta object that includes the path to the .delta file, and
    paths to the query and subject sequence files
    """
    intervals = []
    with open(fname, "r") as dfh:
        # First line is paths to query and subject files
        qpath, spath = dfh.readline().strip().split(" ")
        stem = os.path.splitext(os.path.split(spath)[-1])[0]  # subject filestem
        # Remaining lines are either:
        # - NUCMER (announcing which alignment package was used)
        # - >desc1 desc2 len1 len2   description lines and lengths of input files
        # - seven integers (start/end positions, errors, and stops)
        # - signed integer (symbols between indels)
        # We loop over these, parsing only the lines with seven integers
        for line in [_.strip().split() for _ in dfh.readlines()]:
            if len(line) == 7:
                start, end, sim_errors = (int(line[0]), int(line[1]), int(line[5]))
                if no_exact and (sim_errors == 0):
                    continue
                intervals.append((stem, start, end))
    return NucmerDelta(
        deltafile=fname, queryfile=qpath, subjectfile=spath, query_intervals=intervals
    )
