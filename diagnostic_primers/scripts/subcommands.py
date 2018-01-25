#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcommands.py

Provides subcommand functions for the pdp.py script

- config:        interact with config files
- prodigal:      run prodigal bacterial gene prediction
- eprimer3:      predict primers with EMBOSS ePrimer3
- primersearch:  run in-silico hybridisation with EMBOSS PrimerSearch
- blastscreen:   screen primers with BLASTN
- classify:      classify primers

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

from Bio import (SeqIO, AlignIO)

from diagnostic_primers import (prodigal, eprimer3, blast, primersearch,
                                classify, extract)

from .tools import (log_clines, run_parallel_jobs, load_config_tab,
                    load_config_json, has_primersearch,
                    create_output_directory)


def subcmd_config(args, logger):
    """Run `config` subcommand operations.

    Available subcommands:

    - to_json: convert .tab format config file to JSON and write out
    - to_tab: convert JSON format config file to .tab and write out
    - validate: report on the status of sequence files/classes
    - fix_sequences: stitch multiple sequences together, convert ambiguity
                     symbols to Ns, and write out a new JSON config file

    All subcommands should take either .tab or .json config files,
    distinguished by file extension
    """
    # Determine input config file type
    configtype = os.path.splitext(args.infilename)[-1][1:]
    if configtype not in ('tab', 'json', 'conf'):
        logger.error("Expected config file to end in .conf, .json or .tab " +
                     "got %s (exiting)", configtype)
        raise SystemExit(1)

    if configtype in ('tab', 'conf'):
        coll = load_config_tab(args, logger)
    elif configtype in ('json', ):
        coll = load_config_json(args, logger)

    # Do sequences need to be stitched or their ambiguities replaced?
    # If --validate is active, we report only and do not modify.
    # Note that the underlying code doesn't require us to check whether
    # the sequence files need stitching or symbols replacing, and
    # we could more efficiently use
    # for g in gc.data:
    #     g.stitch()
    #     g.replace_ambiguities()
    # but we're being helpfully verbose.
    logger.info('Checking whether input sequences require stitching, ' +
                'or have non-N ambiguities.')
    problems = ["Validation problems"]  # Holds messages about problem files
    for gcc in coll.data:
        if gcc.needs_stitch:
            msg = '%s requires stitch' % gcc.name
            logger.info(msg)
            problems.append('%s (%s)' % (msg, gcc.seqfile))
            if args.fix_sequences:
                gcc.stitch()
        else:
            logger.info('%s does not require stitch', gcc.name)
        if gcc.has_ambiguities:
            msg = '%s has non-N ambiguities' % gcc.name
            logger.info(msg)
            problems.append('%s (%s)' % (msg, gcc.seqfile))
            if args.fix_sequences:
                gcc.replace_ambiguities()
        else:
            logger.info('%s does not contain non-N ambiguities', gcc.name)
        logger.info('Sequence file: %s', gcc.seqfile)

    # If we were not fixing sequences, report problems
    if not args.fix_sequences:
        if len(problems) > 1:
            logger.warning('\n    '.join(problems))

    # Write post-processing config file and exit
    if args.to_json:
        logger.info('Writing JSON config file to %s', args.to_json)
        coll.write_json(args.to_json)
    elif args.to_tab:
        logger.info('Writing .tab file to %s', args.to_tab)
        coll.write_tab(args.to_tab)
    elif args.fix_sequences:
        logger.info('Writing JSON config file to %s', args.fix_sequences)
        coll.write_json(args.fix_sequences)
    return 0


def subcmd_prodigal(args, logger):
    """Run Prodigal to predict bacterial CDS on the input sequences."""
    # Determine config input type
    configtype = os.path.splitext(args.infilename)[-1][1:]
    if configtype not in ('tab', 'json', 'conf'):
        logger.error("Expected config file to end in .conf, .json or .tab " +
                     "got %s (exiting)", configtype)
        raise SystemExit(1)

    if configtype in ('tab', 'conf'):
        raise ValueError("prodigal subcommand requires JSON config file")
    coll = load_config_json(args, logger)

    # Check if output exists and if we should overwrite
    create_output_directory(args.prodigaldir, args.prodigalforce, logger)

    # Build command-lines for Prodigal and run
    logger.info('Building Prodigal command lines...')
    clines = prodigal.build_commands(coll, args.prodigal_exe, args.prodigaldir)
    log_clines(clines, logger)
    run_parallel_jobs(clines, args, logger)

    # Add Prodigal output files to the GenomeData objects and write
    # the config file
    for gcc in coll.data:
        gcc.features = gcc.cmds['prodigal'].split()[-1].strip()
        logger.info('%s feature file:\t%s' % (gcc.name, gcc.features))
    logger.info('Writing new config file to %s' % args.outfilename)
    coll.write_json(args.outfilename)
    return 0


def subcmd_eprimer3(args, logger):
    """Run ePrimer3 to design primers for each input sequence."""
    # Determine config input type
    configtype = os.path.splitext(args.infilename)[-1][1:]
    if configtype not in ('tab', 'json', 'conf'):
        logger.error("Expected config file to end in .conf, .json or .tab " +
                     "got %s (exiting)", configtype)
        raise SystemExit(1)

    if configtype in ('tab', 'conf'):
        raise ValueError("eprimer3 subcommand requires JSON config file")
    coll = load_config_json(args, logger)

    # Check if output exists and if we should overwrite
    create_output_directory(args.eprimer3_dir, args.eprimer3_force, logger)

    # Build command-lines for ePrimer3 and run
    # This will write 'bare' ePrimer3 files, with unnamed primer pairs
    logger.info('Building ePrimer3 command lines...')
    clines = eprimer3.build_commands(coll, args.eprimer3_exe,
                                     args.eprimer3_dir, vars(args))
    pretty_clines = [str(c).replace(' -', ' \\\n          -') for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs(clines, args, logger)

    # Load bare ePrimer3 data for each input sequence, and write JSON
    # representation with named primer sets
    # Record the path to the JSON representation in the PDPData object
    for gcc in coll.data:
        ep3file = gcc.cmds['ePrimer3'].outfile
        logger.info("Loading primers from ePrimer3 output %s", ep3file)
        primers = eprimer3.load_primers(ep3file, fmt='eprimer3')
        # Write named ePrimer3
        outfname = os.path.splitext(ep3file)[0] + '_named.eprimer3'
        logger.info('Writing named primer sequences to %s' % outfname)
        eprimer3.write_primers(primers, outfname, fmt='ep3')
        # Write named JSON
        outfname = os.path.splitext(ep3file)[0] + '_named.json'
        logger.info('Writing primer JSON sequences to %s' % outfname)
        eprimer3.write_primers(primers, outfname, fmt='json')
        gcc.primers = outfname

    logger.info('Writing new config file to %s' % args.outfilename)
    coll.write_json(args.outfilename)
    return 0


def subcmd_primersearch(args, logger):
    """Perform in silico hybridisation with EMBOSS PrimerSearch."""
    # Does output already exist, and should we overwrite?
    create_output_directory(args.ps_dir, args.ps_force, logger)

    # Get config file data
    coll = load_config_json(args, logger)

    # Construct command lines for primersearch
    logger.info("Building primersearch command-lines...")
    mismatchpercent = int(100 * args.mismatchpercent)  # for EMBOSS
    clines = primersearch.build_commands(coll, args.ps_exe, args.ps_dir,
                                         mismatchpercent)
    pretty_clines = [str(c).replace(' -', ' \\\n          -') for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs(clines, args, logger)

    # Write new config file, and exit
    logger.info('Writing new config file to %s', args.outfilename)
    coll.write_json(args.outfilename)
    return 0


def subcmd_blastscreen(args, logger):
    """Screen primer sequences against a local BLAST database."""
    # We can't go on if there's no BLASTN+ database to check against
    if args.bs_db is None:
        logger.error("No BLASTN+ database provided, cannot perform " +
                     "screen (exiting)")
        raise SystemExit(1)

    # Check if output exists and if we should overwrite
    create_output_directory(args.bs_dir, args.bs_force, logger)

    # Get config file data
    coll = load_config_json(args, logger)

    # Run BLASTN search with primer sequences
    logger.info("Building BLASTN screen command-lines...")
    clines = blast.build_commands(coll, args.bs_exe, args.bs_db, args.bs_dir)
    pretty_clines = [str(c).replace(' -', ' \\\n          -') for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs(clines, args, logger)

    logger.info("BLASTN+ search complete")

    # Amend primer JSON files to remove screened primers
    for blastout, indata in zip([cline.out for cline in clines], coll.data):
        logger.info("Amending primer file %s with results from %s",
                    indata.primers, blastout)
        newprimers = blast.apply_screen(blastout, indata.primers, args.maxaln)
        logger.info("Screened primers placed in %s", newprimers)
        indata.primers = newprimers

    # Write new config file post-BLASTN screen
    logger.info('Writing new config file to %s', args.outfilename)
    coll.write_json(args.outfilename)
    return 0


def subcmd_classify(args, logger):
    """Perform classification of predicted primers."""
    args = args
    logger = logger
    logger.info("Classifying primers for specificity")

    # Create output directory, if needed
    create_output_directory(args.outdir, args.cl_force, logger)

    # Load the JSON config file (post-primersearch)
    coll = load_config_json(args, logger)

    # Test whether the collection has primersearch output
    if not has_primersearch(coll):
        logger.error(' '.join([
            "To use the classify subcommand, the JSON file",
            "must contain links to primersearch data.", "(exiting)"
        ]))
        raise SystemExit(1)
    logger.info("All input genomes have linked path to PrimerSearch data:")
    for genome in coll.data:
        logger.info("\t%s:\t%s", genome.name, genome.primersearch)

    # Obtain classification of all primer sets linked from config file, and
    # report to logger
    results = classify.classify_primers(coll)
    logger.info("Identified primers specific to groups:\n\t%s", '\n\t'.join(
        results.groups))
    for group in results.groups:
        logger.info("Primers specific to %s:\n\t%s", group, '\n\t'.join(
            [primer.name for primer in results.diagnostic_primer(group)]))

    # Write diagnostic primer outputs to the output directory
    classify.write_results(results, os.path.join(args.outdir, 'results.json'))
    classify.write_results(
        results, os.path.join(args.outdir, 'summary.tab'), fmt='summary')


def subcmd_extract(args, logger):
    """Extract amplicons corresponding to primer sets."""
    args = args
    logger = logger
    logger.info("Extracting amplicons for primer set %s", args.primerfile)
    logger.info("PrimerSearch and genome information provided by %s",
                args.infilename)
    if not args.noalign:
        logger.info("MAFFT executable for alignment: %s", args.mafft_exe)

    # Create output directory, if needed
    task_name = os.path.splitext(os.path.split(args.primerfile)[-1])[0]
    outdir = os.path.join(args.outdir, task_name)
    create_output_directory(outdir, args.ex_force, logger)

    # Load the config file and extract the amplicons
    primers = eprimer3.load_primers(args.primerfile, fmt='json')
    coll = load_config_json(args, logger)
    logger.info("Extracting amplicons from source genomes")
    # TODO: make this call for a single primer, and cache sequence
    #       and other (primersearch, primers) data here
    amplicons = extract.extract_amplicons(task_name, primers, coll)

    # Write the amplicons and primers to suitable output files
    # TODO: put this into extract.py as a function write_amplicon_sequences()
    distances = {}
    for pname in amplicons.primer_names:
        seqoutfname = os.path.join(outdir, pname + ".fasta")
        logger.info("Writing amplified sequences for %s to %s", pname,
                    seqoutfname)
        with open(seqoutfname, "w") as ofh:
            SeqIO.write(
                amplicons.get_primer_amplicon_sequences(pname), ofh, 'fasta')

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
                ofh.write(result.stdout.decode("utf-8"))
        else:
            alnoutfname = seqoutfname

        # Calculate distance matrix information
        logger.info("Calculating distance matrices")
        aln = AlignIO.read(open(alnoutfname), 'fasta')  # TODO: avoid file IO
        distances[pname] = extract.calculate_distance(aln)

    # Write distance information to summary file
    distoutfname = os.path.join(outdir, "distances_summary.tab")
    logger.info("Writing distance metric summaries to %s", distoutfname)
    with open(distoutfname, "w") as ofh:
        ofh.write("\t".join([
            "Primer", "Mean distance", "SD", "Min distance", "Max distance"
        ]) + "\n")
        for pname, result in distances.items():
            ofh.write('\t'.join([
                pname,
                "%0.4f" % result.mean,
                "%0.4f" % result.sd,
                "%0.4f" % result.min,
                "%0.4f" % result.max
            ]) + "\n")
