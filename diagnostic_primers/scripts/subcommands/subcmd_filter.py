#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_filter.py

Provides the filter subcommand for pdp.py

(c) The James Hutton Institute 2017-19

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

Copyright (c) 2017-19 The James Hutton Institute
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

from pybedtools import BedTool
from tqdm import tqdm

from diagnostic_primers import PDPException, multiprocessing, prodigal, sge
from diagnostic_primers.nucmer import generate_nucmer_jobs, parse_delta_query_regions
from diagnostic_primers.scripts.tools import (
    collect_existing_output,
    create_output_directory,
    load_config_json,
    log_clines,
    run_parallel_jobs,
)


# Exceptions for filter subcommand
class PDPFilterException(PDPException):
    """Exception thrown by `pdp filter` subcommand"""

    def __init__(self, msg="Error in `pdp filter` subcommand"):
        PDPException.__init__(self, msg)


def subcmd_filter(args, logger):
    """Run filter on input sequences, add filter and filtered_seqfile to config file.

    The subcommand is invoked with `pdp filter <ARGUMENT>` where ARGUMENT is one of a
    restricted set of options:

    - prodigal
    - prodigaligr
    - alnvar

    The general action taken is to calculate (depending on ARGUMENT) a set of features
    for each (eligible) input genome that define valid areas to which primers may be
    defined. This set is written to file and stored in the `features` attr of each
    corresponding PDPData object.

    Using this set of features, a new "filtered" genome is generated, which is the
    extracted set of features, joined by a spacer as a single FASTA file. This is
    written to disk and the path recorded in the `filtered_seqfile` slot of the
    corresponding PDPData object.

    These data are written to the output config JSON file.

    In the following `pdp eprimer3` subcommand, using the `--filter` option will
    cause the tool to use the `filtered_seqfile` sequence for primer design, rather
    than the usual `seqfile` describing the complete genome.

    As later steps in the tool to calculate cross-hybridisation query all primer
    designs against the originating complete genome, as well as comparator
    genomes, the initial location on the `filtered_seqfile` is unimportant.
    """
    # Exactly one of the following filter modes must be selected
    check_filtermodes(
        logger, args.filt_prodigal, args.filt_prodigaligr, args.filt_alnvar is not None
    )

    # Check config input type is valid
    check_config_extension(args.infilename, logger)

    # Load config file
    logger.info("Reading data from %s", args.infilename)
    coll = load_config_json(args, logger)

    # Check if output directory exists and if we should overwrite
    create_output_directory(args.filt_outdir, args.filt_force, logger)

    # Both the --prodigal and --prodigaligr modes require prodigal genecaller output
    # Build command-lines for Prodigal and run
    if args.filt_prodigal or args.filt_prodigaligr:
        # If we are in recovery mode, we are salvaging output from a previous
        # run, and do not necessarily need to rerun all the jobs. In this case,
        # we prepare a list of output files we want to recover from the results
        # in the output directory.
        existingfiles = []
        if args.recovery:
            logger.warning("Entering recovery mode")
            logger.info(
                "\tIn this mode, existing comparison output from %s is reused",
                args.filt_outdir,
            )
            existingfiles = collect_existing_output(args.filt_outdir, "prodigal", args)
            logger.info(
                "Existing files found:\n\t%s", "\n\t".join([_ for _ in existingfiles])
            )

        logger.info("Building Prodigal command lines...")
        clines = prodigal.build_commands(
            coll, args.filt_prodigal_exe, existingfiles, args.filt_outdir
        )
        if len(clines):
            log_clines(clines, logger)
            run_parallel_jobs(clines, args, logger)
        else:
            logger.warning(
                "No prodigal jobs were scheduled (you may see this if the --recovery option is active)"
            )

        # If the filter is prodigal, add Prodigal output files to the GenomeData
        # objects and write the config file; if the filter is prodigaligr, then
        # identify the intergenic loci, and create a new GFF feature file of
        # intergenic regions, with a buffer sequence into the flanking genes
        if args.filt_prodigal:  # Use Prodigal features
            logger.info("Collecting Prodigal prediction output")
            pbar = tqdm(
                coll.data, desc="collecting prodigal output", disable=args.disable_tqdm
            )
            for gcc in pbar:
                gcc.features = gcc.cmds["prodigal"].outfile.split()[-1].strip()
            logger.info("Writing new config file to %s", args.outfilename)
        elif args.filt_prodigaligr:  # Use intergenic regions
            logger.info("Calculating intergenic regions from Prodigal output")
            pbar = tqdm(
                coll.data,
                desc="calculating intergenic regions",
                disable=args.disable_tqdm,
            )
            for gcc in pbar:
                prodigalout = gcc.cmds["prodigal"].outfile.split()[-1].strip()
                bedpath = os.path.splitext(prodigalout)[0] + "_igr.bed"
                prodigal.generate_igr(prodigalout, gcc.seqfile, bedpath)
                gcc.features = bedpath

    # The --alnvar mode requires an all-vs-all comparison of all genomes having the
    # specified class. For instance `pdp filter --alnvar gv01` should carry out
    # all-vs-all nucmer comparisons for the source genomes classified as `gv01`.
    # Once calculated, these are the basis for identification of regions on each of
    # the genomes that are common to all members of the class (here, `gv01`) by
    # calculating the intersections of all alignments on each query genome.
    # The regions are filtered to remove those where the match is exact across all
    # compared genomes.
    if args.filt_alnvar:
        # check the passed argument is an existing class
        filterclass = check_filterclass(args.filt_alnvar, coll, logger)
        logger.info(
            "Attempting to target primer design for class %s to specific and conserved variable regions",
            filterclass,
        )
        # We use a slightly modified location for output
        nucmerdir = os.path.join(args.filt_outdir, "nucmer_output")

        # If we are in recovery mode, we are salvaging output from a previous
        # run, and do not necessarily need to rerun all the jobs. In this case,
        # we prepare a list of output files we want to recover from the results
        # in the output directory.
        existingfiles = []  # Holds list of existing output files
        if args.recovery:
            logger.warning("Entering recovery mode")
            logger.info(
                "\tIn this mode, existing comparison output from %s is reused",
                nucmerdir,
            )
            existingfiles = collect_existing_output(nucmerdir, "alnvar", args)
            logger.info(
                "Existing files found:\n\t%s", "\n\t".join([_ for _ in existingfiles])
            )

        # Collect PDPData objects belonging to the prescribed class
        groupdata = coll.get_groupmembers(filterclass)
        logger.info(
            "Recovered %d PDPData objects with group %s:", len(groupdata), filterclass
        )
        for dataset in groupdata:
            logger.info("\t%s", dataset.name)
        # Carry out nucmer pairwise comparisons for all PDPData objects in the group
        # 1. create directory to hold output
        logger.info("Creating output directory for comparisons: %s", nucmerdir)
        os.makedirs(nucmerdir, exist_ok=True)
        # 2. create and run nucmer comparisons for each of the PDPData objects in the class
        logger.info("Running nucmer pairwise comparisons of group genomes")
        nucmerdata = run_nucmer_comparisons(
            groupdata, nucmerdir, existingfiles, args, logger
        )
        # 3. process alignment files for each of the PDPData objects
        logger.info("Processing nucmer alignment files for the group genomes")
        alndata = process_nucmer_comparisons(groupdata, nucmerdata, args, logger)
        for genome, intervals in alndata:
            bedpath = os.path.join(args.filt_outdir, "%s_alnvar.bed" % genome.name)
            logger.info(
                "Writing interval file to %s and modifying %s PDPData object",
                bedpath,
                genome.name,
            )
            intervals.saveas(bedpath)
            genome.features = bedpath

    # Compile a new genome from the gcc.features file
    logger.info("Compiling filtered sequences from primer design target features")
    logger.info(
        "Spacer length: %d, flanking length: %d",
        args.filt_spacerlen,
        args.filt_flanklen,
    )
    pbar = tqdm(coll.data, desc="compiling filtered genomes", disable=args.disable_tqdm)
    for gcc in [_ for _ in pbar if _.features is not None]:
        stem, ext = os.path.splitext(gcc.seqfile)
        if args.filt_outdir is None:
            filtered_path = stem + "_" + args.filt_suffix + ext
        else:
            filtered_path = (
                os.path.join(args.filt_outdir, os.path.split(stem)[-1])
                + "_"
                + args.filt_suffix
                + ext
            )
        gcc.create_filtered_genome(
            filtered_path, args.filt_spacerlen, args.filt_suffix, args.filt_flanklen
        )

    # Write updated config file
    logger.info("Writing new config file to %s", args.outfilename)
    coll.write_json(args.outfilename)

    return 0


def check_filtermodes(logger, *filtermodes):
    """Raise an exception if the filter modes are invalid

    The modes argument is a variable number of True/False or 1/0 values
    corresponding to mutually exclusive choices in the pdp filter parser
    """
    if sum(filtermodes) != 1:
        logger.error("Invalid filter options chosen (exiting)")
        if sum(filtermodes) == 0:
            logger.error("No filter modes chosen.")
        if sum(filtermodes) > 0:
            logger.error("More than one filter mode chosen")
        raise PDPFilterException("Incorrect filter mode specified")


def check_config_extension(fname, logger):
    """Raise an exception if the config file type is not valid

    We expect a JSON file with extension .json, but do not check for
    valid formatting.
    """
    configtype = os.path.splitext(fname)[-1]
    # check extension
    if configtype not in (".tab", ".json", ".conf"):
        logger.error(
            "Expected config file to end in .conf, .json or .tab " + "got %s (exiting)",
            configtype,
        )
        raise PDPFilterException("Could not read configuration file")
    # raise issue if type is not correct
    if configtype in (".tab", ".conf"):
        raise PDPFilterException(
            "filter subcommand requires JSON config file, please convert input"
        )


def check_filterclass(group, coll, logger):
    """Return group or raise exception if not found in the collection"""
    if group not in coll.groups:
        logger.error(
            "Class %s is not defined in the configuration file (exiting)", group
        )
        raise PDPFilterException("--alnvar class is not present in the dataset")
    return group


def run_nucmer_comparisons(groupdata, outdir, existingfiles, args, logger):
    """Run nucmer alignments and eturn iterable of NucmerOutput objects

    :param groupdata:         iterable of PDPData objects
    :param outdir:            parent directory for nucmer output
    :param existingfiles:     iterable of existing nucmer results files to skip generation
    :param args:              pdp filter command-line arguments
    :param logger:            a logging object

    For each PDPData object in groupdata, generate a Job and a NucmerOutput
    object, then pass the jobs to the appropriate scheduler.

    Return an iterable of NucmerOutput objects, one per nucmer alignment.
    """
    logger.info("Composing nucmer command-lines into jobs:")
    nucmer_jobs = generate_nucmer_jobs(
        groupdata,
        outdir,
        args.nucmer_exe,
        args.deltafilter_exe,
        args.maxmatch,
        args.jobprefix,
    )
    jobs, nucmerdata = zip(*nucmer_jobs)

    # Avoid jobs with existing output
    runjobs = [
        _ for _ in jobs if os.path.split(_.command.outfile)[-1] not in existingfiles
    ]
    if len(runjobs) == 0:
        logger.warning(
            "No nucmer jobs were scheduled (you may see this if the --recovery option is active)"
        )
    else:
        logger.info("Running the following jobs:")
        for job in runjobs:
            logger.info("\t%s", job.name)
        logger.info("Running jobs with scheduler: %s", args.scheduler)
        if args.scheduler == "multiprocessing":
            multiprocessing.run_dependency_graph(runjobs, args.workers, logger)
        elif args.scheduler == "SGE":
            sge.run_dependency_graph(
                runjobs,
                logger,
                jgprefix=args.jobprefix,
                sgegroupsize=args.sgegroupsize,
                sgeargs=args.sgeargs,
            )
        else:
            logger.error("Scheduler %s not recognised (exiting)", args.scheduler)
            raise PDPFilterException("Scheduler not recognised by PDP")
    return nucmerdata


def process_nucmer_comparisons(groupdata, nucmerdata, args, logger):
    """Return intervals on each genome describing regions aligned to all group members

    groupdata         iterable of PDPData objects
    nucmerdata        iterable of nucmer.NucmerOutput namedtuples
    args              pdp filter command-line arguments
    logger            a logging object

    For each PDPData object in groupdata, identify comparisons in nucmerdata where
    that object is the query, and process the corresponding output file into a set
    of intervals. When processing the output file, do not use aligned intervals
    where there are no similarity errors.

    When all comparisons are collected for a query, identify the intervals
    corresponding to the intersection of all alignments to the query.

    For each PDPData object, return the collection of all alignment intervals,
    plus the intersection.
    """
    # Collect results into a single list for return
    intervals = []
    for genome in groupdata:
        logger.info("Processing %s nucmer files", genome.name)
        nucmer_out = [_ for _ in nucmerdata if _.query == genome]
        nucmer_results = [
            parse_delta_query_regions(
                _.out_delta,
                min_sim_errors=args.filt_minsecount,
                min_err_rate=args.filt_minserate,
            )
            for _ in nucmer_out
        ]
        nucmer_intervals = [BedTool(_.query_intervals) for _ in nucmer_results]
        if len(nucmer_intervals) == 1:
            common_regions = nucmer_intervals[0].sort().merge()
        else:
            common_regions = (
                nucmer_intervals[0].intersect(nucmer_intervals[1:]).sort().merge()
            )
        intervals.append((genome, common_regions.sort()))

    logger.info("Common regions identified:")
    for genome, regions in intervals:
        logger.info(
            "\t%s: regions: %d, total length: %d",
            genome.name,
            len(regions),
            regions.total_coverage(),
        )

    return intervals


def recursive_intersection(bedtools, current=None):
    """Return a BedTool describing the recursive intersection of BedTool objects

    bedtools behaviour is

    a.intersect(b).intersect(c) != a.intersect([b, c])

    The first term identifies regions in a & b & c, the second identifies
    regions in a & b | a & c.

    To obtain the first behaviour, we recursively apply .intersect() to
    each of the elements in the passed iterable of bedtools, to identify
    regions common to all three (or more) sets.

    The result of these intersections may contain duplicate intervals,
    so the final result is sorted and merged.

    The returned BedTool object describes regions common to all passed
    BedTools, relative to the first BedTool in the iterable
    """
    if len(bedtools) == 0:
        return current.sort().merge()
    if current is None:
        current = bedtools.pop()
    else:
        current = current.intersect(bedtools.pop())
    return recursive_intersection(bedtools, current)


def chained_intersection(bedtools):
    """Return a BedTool describing the chained intersection of BedTool objects

    bedtools behaviour is

    a.intersect(b).intersect(c) != a.intersect([b, c])

    The first term identifies regions in a & b & c, the second identifies
    regions in a & b | a & c.

    To obtain the first behaviour, we apply .intersect() to
    each of the elements in the passed iterable of bedtools, starting
    at the last, then working to the first, to identify
    regions common to all three (or more) sets.

    The result of these intersections may contain duplicate intervals,
    so the final result is sorted and merged.

    The returned BedTool object describes regions common to all passed
    BedTools, relative to the first BedTool in the iterable
    """
    current = bedtools.pop()
    while len(bedtools):
        current = current.sort().merge().intersect(bedtools.pop())
    return current.sort().merge()
