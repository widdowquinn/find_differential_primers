import sys
import traceback

from diagnostic_primers import (multiprocessing, process, prodigal,
                                eprimer3, sge, sge_jobs, blast)


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Load a config file and record details in a logger
def load_config_file(args, logger):
    """Load config file and return a GenomeCollection."""
    try:
        gc = process.load_collection(args.infilename, name='pdp.py')
    except:
        logger.error('Could not parse config file %s (exiting)' %
                     args.infilename)
        logger.error(last_exception())
        raise SystemExit(1)
    logger.info('Parsed config file %s: %d sequences in %d groups' %
                (args.infilename, len(gc), len(gc.groups())))
    logger.info('Diagnostic groups:\n%s' %
                '\n'.join(['\t%s' % g for g in gc.groups()]))
    return gc


# Report a list of command lines to a logger, in pretty format
def log_clines(clines, logger):
    """Log command-lines, one per line."""
    logger.info('...%d commands returned:\n%s' %
                (len(clines), '\n'.join(['\t%s' % c for c in clines])))
