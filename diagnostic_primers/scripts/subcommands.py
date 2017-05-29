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

(c) The James Hutton Institute 2017

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

Copyright (c) 2017 The James Hutton Institute
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

from diagnostic_primers import (prodigal, eprimer3, blast)

from .tools import (load_config_file, log_clines, run_parallel_jobs)


def subcmd_config(args, logger):
    """Run the `config` subcommand operations. Here, the goal is to parse the
    input config file, verify that the config file can be read, report some
    basic information back (if verbose), and clean up sequence data by
    stitching multiple contigs/sequences together and replacing ambiguity
    symbols with 'N'.

    If sequence data needs to be stitched, or symbols replaced, then new
    sequence files are produced and written (if --validate is not in
    operation)

    A new config file, pointing to the revised files, is written out (if
    --validate is not in operation).

    """
    coll = load_config_file(args, logger)

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
    for gcc in coll.data:
        if gcc.needs_stitch:
            logger.info('%s requires stitch' % gcc.name)
            if not args.validate:
                gcc.stitch()
        else:
            logger.info('%s does not require stitch' % gcc.name)
        if gcc.has_ambiguities():
            logger.info('%s contains non-N ambiguities' % gcc.name)
            if not args.validate:
                gcc.replace_ambiguities()
        else:
            logger.info('%s does not contain non-N ambiguities' % gcc.name)
        logger.info('Sequence file: %s' % gcc.seqfile)

    # Write post-processing config file and exit
    if not args.validate:
        logger.info('Writing processed config file to %s' %
                    args.outfilename)
        coll.write(args.outfilename)
    return 0


def subcmd_prodigal(args, logger):
    """Run Prodigal to predict bacterial CDS on the input sequences."""
    coll = load_config_file(args, logger)

    # Build command-lines for Prodigal and run
    logger.info('Building Prodigal command lines...')
    if args.prodigalforce:
        logger.warning('Forcing Prodigal to run. This may overwrite ' +
                       'existing output.')
    else:
        logger.info('Prodigal will fail if output directory exists')
    clines = prodigal.build_commands(coll, args.prodigal_exe,
                                     args.prodigaldir, args.prodigalforce)
    log_clines(clines, logger)
    run_parallel_jobs(clines, args, logger)

    # Add Prodigal output files to the GenomeData objects and write
    # the config file
    for gcc in coll.data:
        gcc.features = gcc.cmds['prodigal'].split()[-1].strip()
        logger.info('%s feature file:\t%s' % (gcc.name, gcc.features))
    logger.info('Writing new config file to %s' % args.outfilename)
    coll.write(args.outfilename)
    return 0


def subcmd_eprimer3(args, logger):
    """Run ePrimer3 to design primers for each input sequence."""
    coll = load_config_file(args, logger)

    # Build command-lines for ePrimer3 and run
    logger.info('Building ePrimer3 command lines...')
    if args.eprimer3_force:
        logger.warning('Forcing ePrimer3 to run. This may overwrite ' +
                       'existing output.')
    else:
        logger.info('ePrimer3 may fail if output directory exists')
    clines = eprimer3.build_commands(coll, args.eprimer3_exe,
                                     args.eprimer3_dir,
                                     args.eprimer3_force, vars(args))
    pretty_clines = [str(c).replace(' -', '\\\n\t\t-') for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs(clines, args, logger)

    # Load ePrimer3 data for each input sequence, and write JSON representation
    # Record JSON representation in the object
    for gcc in coll.data:
        logger.info("Loading/naming primers in %s" %
                    gcc.cmds['ePrimer3'].outfile)
        primers, outfname = eprimer3.load_primers(gcc.cmds['ePrimer3'].outfile)
        logger.info('Named primers ePrimer3 file created:\t%s' % outfname)
        outfname = os.path.splitext(outfname)[0] + '.json'
        logger.info('Writing primers to %s' % outfname)
        eprimer3.primers_to_json(primers, outfname)
        gcc.primers = outfname

    logger.info('Writing new config file to %s' % args.outfilename)
    coll.write(args.outfilename)
    return 0


def subcmd_primersearch(args, logger):
    """Perform in silico hybridisation with EMBOSS PrimerSearch."""
    args = args
    logger = logger
    raise NotImplementedError


def subcmd_blastscreen(args, logger):
    """Screen primer sequences against a local BLAST database."""
    # We can't go on if there's no BLASTN+ database to check against
    if args.bs_db is None:
        logger.error("No BLASTN+ database provided, cannot perform " +
                     "screen (exiting)")
        return 1

    coll = load_config_file(args, logger)

    for gcc in coll.data:
        logger.info("Loading primer data from %s" % gcc.primers)
        gcc.fastafname = eprimer3.json_to_fasta(gcc.primers)
        logger.info("Wrote primer FASTA sequences to %s" % gcc.fastafname)

    logger.info("Building BLASTN screen command-lines...")
    if args.bs_force:
        logger.warning("Forcing BLASTN screen to run. This may " +
                       "overwrite existing output.")
    else:
        logger.info("BLASTN screen may fail if output directory exists.")
    clines = blast.build_commands(coll,
                                  args.bs_db, args.bs_exe,
                                  args.bs_dir, args.bs_force)
    pretty_clines = [str(c).replace(' -', '\\\n\t\t-') for c in clines]
    log_clines(pretty_clines, logger)
    run_parallel_jobs(clines, args, logger)

    logger.info("BLASTN+ screen complete")
    return 0


def subcmd_classify(args, logger):
    """Perform classification of predicted primers."""
    args = args
    logger = logger
    raise NotImplementedError