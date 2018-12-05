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

from itertools import permutations

from diagnostic_primers.sge_jobs import Job


# Generate list of Job objects, one per NUCmer run
def generate_nucmer_jobs(
    filenames, outdir, nucmer_exe, filter_exe, maxmatch=False, jobprefix="PDPNUCmer"
):
    """Return list of Jobs describing NUCmer command-lines for PDP.

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files, generating Jobs describing NUCmer command lines
    for each pairwise comparison.
    """
    ncmds, fcmds = generate_nucmer_commands(
        filenames, outdir, nucmer_exe, filter_exe, maxmatch
    )
    joblist = []
    for idx, ncmd in enumerate(ncmds):
        njob = Job("%s_%06d-n" % (jobprefix, idx), ncmd)
        fjob = Job("%s_%06d-f" % (jobprefix, idx), fcmds[idx])
        fjob.add_dependency(njob)
        joblist.append(fjob)
    return joblist


# Generate list of NUCmer pairwise comparison command lines from
# passed sequence filenames
def generate_nucmer_commands(filenames, outdir, nucmer_exe, filter_exe, maxmatch=False):
    """Return list of NUCmer command-lines for PDP.

    The first element returned is a list of NUCmer commands, and the
    second a corresponding list of delta_filter_wrapper.py commands.
    The NUCmer commands should each be run before the corresponding
    delta-filter command.

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files generating NUCmer command lines for each
    pairwise comparison.
    """
    nucmer_cmdlines, delta_filter_cmdlines = [], []
    comparisons = permutations(filenames, 2)
    for (fname1, fname2) in comparisons:
        ncmd, dcmd = construct_nucmer_cmdline(
            fname1, fname2, outdir, nucmer_exe, filter_exe, maxmatch
        )
        nucmer_cmdlines.append(ncmd)
        delta_filter_cmdlines.append(dcmd)
    return (nucmer_cmdlines, delta_filter_cmdlines)


# Generate single NUCmer pairwise comparison command line from pair of
# input filenames
def construct_nucmer_cmdline(
    fname1, fname2, outdir, nucmer_exe, filter_exe, maxmatch=False
):
    """Return a tuple of corresponding NUCmer and delta-filter commands

    The split into a tuple was made necessary by changes to SGE/OGE.
    The delta-filter command must now be run as a dependency of the NUCmer
    command, and be wrapped in a Python script to capture STDOUT.

    NOTE: This command-line writes output data to a subdirectory of the passed
    outdir, called "nucmer_output".

    - fname1 - query FASTA filepath
    - fname2 - subject FASTA filepath
    - outdir - path to output directory
    - maxmatch - Boolean flag indicating whether to use NUCmer's -maxmatch
    option. If not, the -mum option is used instead
    """
    outprefix = os.path.join(
        outdir,
        "%s_vs_%s"
        % (
            os.path.splitext(os.path.split(fname1)[-1])[0],
            os.path.splitext(os.path.split(fname2)[-1])[0],
        ),
    )
    if maxmatch:
        mode = "--maxmatch"
    else:
        mode = "--mum"
    nucmercmd = "{0} {1} -p {2} {3} {4}".format(
        nucmer_exe, mode, outprefix, fname1, fname2
    )
    filtercmd = "delta_filter_wrapper.py " + "{0} -1 {1} {2}".format(
        filter_exe, outprefix + ".delta", outprefix + ".filter"
    )
    return (nucmercmd, filtercmd)
