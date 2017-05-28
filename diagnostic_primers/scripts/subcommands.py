import sys

from diagnostic_primers import (multiprocessing, process, prodigal,
                                eprimer3, sge, sge_jobs, blast)

from .tools import (last_exception, load_config_file, log_clines)


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
    gc = load_config_file(args, logger)

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
    for g in gc.data:
        if g.needs_stitch:
            logger.info('%s requires stitch' % g.name)
            if not args.validate:
                g.stitch()
        else:
            logger.info('%s does not require stitch' % g.name)
        if g.has_ambiguities():
            logger.info('%s contains non-N ambiguities' % g.name)
            if not args.validate:
                g.replace_ambiguities()
        else:
            logger.info('%s does not contain non-N ambiguities' % g.name)
        logger.info('Sequence file: %s' % g.seqfile)

    # Write post-processing config file and exit
    if not args.validate:
        logger.info('Writing processed config file to %s' %
                    args.outfilename)
        gc.write(args.outfilename)
    return 0


def subcmd_prodigal(args, logger):
    """Run Prodigal to predict bacterial CDS on the input sequences."""
    gc = load_config_file(args, logger)

    # Build command-lines for Prodigal and run
    logger.info('Building Prodigal command lines...')
    if args.prodigalforce:
        logger.warning('Forcing Prodigal to run. This may overwrite ' +
                       'existing output.')
    else:
        logger.info('Prodigal will fail if output directory exists')
    clines = prodigal.build_commands(gc, args.prodigal_exe,
                                     args.prodigaldir, args.prodigalforce)
    log_clines(clines)
    run_parallel_jobs(clines)

    # Add Prodigal output files to the GenomeData objects and write
    # the config file
    for g in gc.data:
        g.features = g.cmds['prodigal'].split()[-1].strip()
        logger.info('%s feature file:\t%s' % (g.name, g.features))
    logger.info('Writing new config file to %s' % args.outfilename)
    gc.write(args.outfilename)
    return 0


def subcmd_eprimer3(args, logger):
    """Run ePrimer3 to design primers for each input sequence."""
    gc = load_config_file(args, logger)

    # Build command-lines for ePrimer3 and run
    logger.info('Building ePrimer3 command lines...')
    if args.eprimer3_force:
        logger.warning('Forcing ePrimer3 to run. This may overwrite ' +
                       'existing output.')
    else:
        logger.info('ePrimer3 may fail if output directory exists')
    clines = eprimer3.build_commands(gc, args.eprimer3_exe,
                                     args.eprimer3_dir,
                                     args.eprimer3_force, vars(args))
    pretty_clines = [str(c).replace(' -', '\\\n\t\t-') for c in clines]
    log_clines(pretty_clines)
    run_parallel_jobs(clines)

    # Load ePrimer3 data for each input sequence, and write JSON representation
    # Record JSON representation in the object
    for g in gc.data:
        logger.info("Loading/naming primers in %s" %
                    g.cmds['ePrimer3'].outfile)
        primers, outfname = eprimer3.load_primers(g.cmds['ePrimer3'].outfile)
        logger.info('Named primers ePrimer3 file created:\t%s' % outfname)
        outfname = os.path.splitext(outfname)[0] + '.json'
        logger.info('Writing primers to %s' % outfname)
        eprimer3.primers_to_json(primers, outfname)
        g.primers = outfname

    logger.info('Writing new config file to %s' % args.outfilename)
    gc.write(args.outfilename)
    return 0


def subcmd_primersearch(args, logger):
    """Perform in silico hybridisation with EMBOSS PrimerSearch."""
    raise NotImplementedError


def subcmd_blastscreen(args, logger):
    """Screen primer sequences against a local BLAST database."""
    # We can't go on if there's no BLASTN+ database to check against
    if args.bs_db is None:
        logger.error("No BLASTN+ database provided, cannot perform " +
                     "screen (exiting)")
        return 1

    gc = load_config_file()

    for g in gc.data:
        logger.info("Loading primer data from %s" % g.primers)
        g.fastafname = eprimer3.json_to_fasta(g.primers)
        logger.info("Wrote primer FASTA sequences to %s" % g.fastafname)

    logger.info("Building BLASTN screen command-lines...")
    if args.bs_force:
        logger.warning("Forcing BLASTN screen to run. This may " +
                       "overwrite existing output.")
    else:
        logger.info("BLASTN screen may fail if output directory exists.")
    clines = blast.build_commands(gc,
                                  args.bs_db, args.bs_exe,
                                  args.bs_dir, args.bs_force)
    pretty_clines = [str(c).replace(' -', '\\\n\t\t-') for c in clines]
    log_clines(pretty_clines)
    run_parallel_jobs(clines)

    logger.info("BLASTN+ screen complete")
    return 0


def subcmd_classify(args, logger):
    """Perform classification of predicted primers."""
    raise NotImplementedError
