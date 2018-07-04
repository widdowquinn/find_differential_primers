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

from ..tools import (create_output_directory, load_config_json)


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
    # in turn
    primers = eprimer3.load_primers(args.primerfile, fmt='json')
    coll = load_config_json(args, logger)
    logger.info("Extracting amplicons from source genomes")
    # TODO: make this call for a single primer, and cache sequence
    #       and other (primersearch, primers) data here
    distances = {}
    for primer in primers:
        amplicons = extract.extract_amplicons(task_name, primer, coll)

        # Write the amplicons and primers to suitable output files
        # TODO: put this into extract.py as a function write_amplicon_sequences()
        for pname in amplicons.primer_names:
            seqoutfname = os.path.join(outdir, pname + ".fasta")
            logger.info("Writing amplified sequences for %s to %s", pname,
                        seqoutfname)
            with open(seqoutfname, "w") as ofh:
                seqdata = amplicons.get_primer_amplicon_sequences(pname)
                # Order sequence data for consistent output (aids testing)
                seqdata = [
                    _[1] for _ in sorted([(seq.id, seq) for seq in seqdata])
                ]
                SeqIO.write(seqdata, ofh, 'fasta')

            # Align the sequences with MAFFT
            # TODO: Neaten this up so we're multiprocessing it and catching output/
            #       errors
            if not args.noalign:
                alnoutfname = os.path.join(outdir, pname + ".aln")
                logger.info("Aligning amplicons with MAFFT and writing to %s",
                            alnoutfname)
                # MAFFT is run with --quiet flag to suppress verbiage in STDERR
                result = subprocess.run(
                    [args.mafft_exe, "--quiet", seqoutfname],
                    stdout=subprocess.PIPE)
                if result.returncode:  # MAFFT failed
                    logger.error(
                        "There was an error aligning %s with MAFFT (exiting)",
                        seqoutfname)
                    raise SystemExit(1)
                with open(alnoutfname, "w") as ofh:
                    ofh.write(result.stdout.decode('utf-8'))
            else:
                alnoutfname = seqoutfname

            # Calculate distance matrix information
            logger.info("Calculating distance matrices")
            aln = AlignIO.read(open(alnoutfname),
                               'fasta')  # TODO: avoid file IO
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
