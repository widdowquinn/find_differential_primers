#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_extract.py

Provides the extract subcommand for pdp.py

(c) The James Hutton Institute 2017-18

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

Copyright (c) 2017-18 The James Hutton Institute
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
import subprocess

from Bio import (AlignIO, SeqIO)

from diagnostic_primers import (eprimer3, extract)

from ..tools import (create_output_directory, load_config_json,
                     run_parallel_jobs)


def subcmd_extract(args, logger):
    """Extract amplicons corresponding to primer sets."""
    logger.info("Extracting amplicons for primer set %s", args.primerfile)
    logger.info("PrimerSearch and genome information provided by %s",
                args.infilename)
    if not args.noalign:
        logger.info("MAFFT executable for alignment: %s", args.mafft_exe)

    # Create output directory, if needed
    task_name = os.path.splitext(os.path.split(args.primerfile)[-1])[0]
    outdir = os.path.join(args.outdir, task_name)
    create_output_directory(outdir, args.ex_force, logger)

    # Load the config file and extract the amplicons for each primer set
    # in turn. Put the amplicons into a .fasta file and record the location
    # for each primer set
    primers = eprimer3.load_primers(args.primerfile, fmt='json')
    coll = load_config_json(args, logger)
    logger.info("Extracting amplicons from source genomes")
    amplicon_fasta = {}
    seq_cache = {}
    for primer in primers:
        amplicons, seq_cache = extract.extract_amplicons(
            task_name, primer, coll, seq_cache=seq_cache)

        # Write the amplicons and primers to FASTA, and record the path against
        # the primer name
        for pname in amplicons.primer_names:
            seqoutfname = os.path.join(outdir, pname + ".fasta")
            logger.info("Writing amplified sequences for %s to %s", pname,
                        seqoutfname)
            amplicons.write_amplicon_sequences(pname, seqoutfname)
            amplicon_fasta[pname] = seqoutfname

    # Align the sequences with MAFFT
    amplicon_alnfiles = {}
    if not args.noalign:
        clines = []
        logger.info("Compiling MAFFT alignment commands")
        for pname, fname in amplicon_fasta.items():
            alnoutfname = os.path.join(outdir, pname + ".aln")
            amplicon_alnfiles[pname] = alnoutfname
            # MAFFT is run with --quiet flag to suppress verbiage in STDERR
            cline = "{} --quiet {} > {}".format(args.mafft_exe, fname,
                                                alnoutfname)
            clines.append(cline)
        # Pass command-lines to the appropriate scheduler
        logger.info("Aligning amplicons with MAFFT")
        run_parallel_jobs(clines, args, logger)
    else:
        # If we're not aligning, reuse the FASTA files
        amplicon_alnfiles = amplicon_fasta

    # Calculate distance matrix information
    distances = {}
    logger.info("Calculating distance matrices")
    for pname, fname in amplicon_alnfiles.items():
        aln = AlignIO.read(open(fname), 'fasta')
        distances[pname] = extract.calculate_distance(aln)

    # Write distance information to summary file
    distoutfname = os.path.join(outdir, "distances_summary.tab")
    logger.info("Writing distance metric summaries to %s", distoutfname)
    with open(distoutfname, "w") as ofh:
        ofh.write("\t".join([
            "primer", "dist_mean", "dist_sd", "dist_min", "dist_max", "unique",
            "nonunique"
        ]) + "\n")
        # Note: ordered output for the table
        for pname, result in sorted(distances.items()):
            ofh.write('\t'.join([
                pname,
                "%0.4f" % result.mean,
                "%0.4f" % result.sd,
                "%0.4f" % result.min,
                "%0.4f" % result.max,
                "%d" % result.unique,
                "%d" % result.nonunique
            ]) + "\n")

    return 0
