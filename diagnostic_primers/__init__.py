# -*- coding: utf-8 -*-
__version__ = "0.2.0.dev"

import json
import os
import re
import sys
import traceback

from Bio import SeqIO
from Bio.Emboss import Primer3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Define base PDPException
class PDPException(Exception):
    """General exception for PDP/find_differential_primers"""

    def __init__(self, msg="Error in pdp/diagnostic_primers"):
        Exception.__init__(self, msg)


# Report last exception as string
def last_exception():
    """Return last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))


# Code to read/write primers
class PrimersEncoder(json.JSONEncoder):
    """JSON encoder for Primer3.Primers objects."""

    def default(self, obj):
        if not isinstance(obj, Primer3.Primers):
            return json.JSONEncoder.default(self, obj)

        # Convert complex Primer3.Primers object to serialisable dictionary
        # and return
        return obj.__dict__


def load_primers(infname, fmt="eprimer3", noname=False):
    """Load primers from a file.

    :param infname: path to input file describing primers
    :param fmt: expected format of primer input file
    :param noname: Boolean flag. If not set to true, each primer receives
              a unique name based on the input filename.


    The function can load JSON or ePrimer3 files - ePrimer3 by default.
    """
    if fmt in ("ep3", "eprimer3"):
        return __load_primers_eprimer3(infname, noname)
    elif fmt in ("p3", "primer3"):
        return __load_primers_primer3(infname, noname)
    elif fmt in ("json",):
        return __load_primers_json(infname)


def __load_primers_eprimer3(infname, noname=False):
    """Loads and names primers from the passed ePrimer3 file.

    :param infname: path to input file describing primers in ePrimer3 format
    :param noname: Boolean flag. If not set to true, each primer receives
              a unique name based on the input filename.

    Primers are returned as a list of Primer3.Primers objects
    """
    # Load primers with Biopython. This does not respect a 'name' attribute,
    # as the bare ePrimer3 files don't provide names for primer sets.
    with open(infname, "r") as primerfh:
        ep3primers = Primer3.read(primerfh).primers

    # Create unique identifier for each primer set, based on the input
    # filename
    if not noname:
        for idx, primer in enumerate(ep3primers, 1):
            stem = os.path.splitext(os.path.split(infname)[-1])[0]
            primer.name = "%s_primer_%05d" % (stem, idx)

    # Return primers
    return ep3primers


def __load_primers_json(infname):
    """Loads and returns primers from a JSON file.

    :param infname: path to input file describing primers in JSON format

    Primers are returned as a list of Primer3.Primers objects
    """
    primers = []
    with open(infname, "r") as primerfh:
        for pdata in json.load(primerfh):
            primer = Primer3.Primers()
            for k, v in pdata.items():
                setattr(primer, k, v)
            primers.append(primer)
    return primers


def __parse_line(line, primer):
    """Return Primer3.Primers object, modified by the Primer3 output line

    :param line:  single line from Primer3 (v2+) output
    :param primer:  Primer3.Primers object

    """
    key, val = line.split("=")
    info = key.split("_", 3)
    if info[1] == "LEFT":
        if len(info) == 3:
            start, prodlen = val.split(",")
            primer.forward_start = int(start)
            primer.forward_length = int(prodlen)
        if info[-1] == "SEQUENCE":
            primer.forward_seq = str(val)
        if info[-1] == "TM":
            primer.forward_tm = float(val)
        if info[-1] == "GC_PERCENT":
            primer.forward_gc = float(val)
    elif info[1] == "RIGHT":
        if len(info) == 3:
            start, prodlen = val.split(",")
            primer.reverse_start = int(start)
            primer.reverse_length = int(prodlen)
        if info[-1] == "SEQUENCE":
            primer.reverse_seq = str(val)
        if info[-1] == "TM":
            primer.reverse_tm = float(val)
        if info[-1] == "GC_PERCENT":
            primer.reverse_gc = float(val)
    elif info[1] == "PAIR":
        if info[-1] == "PRODUCT_SIZE":
            primer.size = int(val)
    return primer


def __load_primers_primer3(infname, noname=False):
    """Loads and names primers from the passed Primer3 (v2+) file.

    :param infname: path to input file describing primers in ePrimer3 format
    :param noname: Boolean flag. If not set to true, each primer receives
              a unique name based on the input filename.

    Primers are returned as a list of Primer3.Primers objects
    """
    # Load primers with Biopython. This does not respect a 'name' attribute,
    # as the bare Primer3 files don't provide names for primer sets.
    stem = os.path.splitext(os.path.split(infname)[-1])[0]
    p3primers = []
    with open(infname, "r") as pfh:
        primer = None
        for line in [_.strip() for _ in pfh.readlines() if len(_.strip())]:
            if line == "=":  # end of file
                p3primers.append(primer)
                primer = None
                break
            # Instantiate new primer on PRIMER_PAIR_[0-9]_PENALTY
            primeridx = re.search("(?<=PRIMER_PAIR_)[0-9]*(?=_PENALTY)", line)
            if primeridx is not None:
                if primer is not None:  #  Record last primer
                    p3primers.append(primer)
                primer = Primer3.Primers()
                if not noname:
                    primer.name = "{}_primer_{:05d}".format(
                        stem, int(primeridx.group()) + 1
                    )
            # If there's a primer object to populate, read the PRIMER_LEFT,
            # PRIMER_RIGHT and PRIMER_PAIR data, matching regexes to an
            # attribute dictionary
            elif primer is not None:
                primer = __parse_line(line, primer)

    return p3primers


def write_primers(primers, outfilename, fmt="fasta"):
    """Write Primer3.Primers to file.

    :param primers:  collection of Biopython primer objects
    :param outfilename:  path to output file
    :param fmt:  sequence format to write

    TODO: distribution dictionary
    """
    # Ensure output directory exists
    outdir = os.path.join(*os.path.split(outfilename)[:-1])
    os.makedirs(outdir, exist_ok=True)

    # Order primers before writing
    primers = [_[1] for _ in sorted([(primer.name, primer) for primer in primers])]
    if fmt in ("json",):
        __write_primers_json(primers, outfilename)
    elif fmt in ("ep3", "eprimer3"):
        __write_primers_eprimer3(primers, outfilename)
    elif fmt in ("tsv", "tab"):
        __write_primers_tsv(primers, outfilename)
    elif fmt in ("bed"):
        __write_primers_bed(primers, outfilename)
    else:
        __write_primers_seqio(primers, outfilename, fmt)


def __write_primers_seqio(primers, outfilename, fmt):
    """Write primers to file, using SeqIO.

    :param primers:  iterable of Primer3.Primers objects
    :param outfilename:  path to output file
    :param fmt:  string describing format for primer file output

    Returns the number of records written
    """
    seqrecords = []

    for primer in primers:
        seqrecords.append(
            SeqRecord(Seq(primer.forward_seq), id=primer.name + "_fwd", description="")
        )
        seqrecords.append(
            SeqRecord(Seq(primer.reverse_seq), id=primer.name + "_rev", description="")
        )
        if len(primer.internal_seq):  # This is '' id no oligo
            seqrecords.append(
                SeqRecord(
                    Seq(primer.internal_seq), id=primer.name + "_int", description=""
                )
            )

    return SeqIO.write(seqrecords, outfilename, fmt)


def __write_primers_tsv(primers, outfname):
    """Write primers to file in three-column tab-separated format.

    :param primers:  iterable of Primer3.Primers objects
    :param outfname:  path to output file

    This is required for primersearch input with EMBOSS
    """
    # Don't use more than one newline at the end of the header, or
    # else primersearch treats the blank line as a primer!
    header = (
        "\n".join(["# EPRIMER3 PRIMERS %s" % outfname, "# Name       FWD        REV"])
        + "\n"
    )
    with open(outfname, "w") as outfh:
        outfh.write(header)
        for primer in primers:
            outfh.write(
                "\t".join([primer.name, primer.forward_seq, primer.reverse_seq]) + "\n"
            )


def __write_primers_eprimer3(primers, outfname):
    """Write Primer3 primer objects in ePrimer3 format (Extended).

    :param primers:  iterable of Primer3.Primers objects
    :param outfname:  path to output file
    """
    header = (
        "\n".join(
            [
                "# EPRIMER3 PRIMERS %s " % outfname,
                "#                      Start  Len   Tm     " + "GC%   Sequence",
            ]
        )
        + "\n"
    )

    with open(outfname, "w") as outfh:
        outfh.write(header)
        for idx, primer in enumerate(primers, 1):
            outfh.write("# %s\n" % primer.name)
            outfh.write("%-4d PRODUCT SIZE: %d\n" % (idx, primer.size))
            outfh.write(
                "     FORWARD PRIMER  %-9d  %-3d  %.02f  %.02f  %s\n"
                % (
                    primer.forward_start,
                    primer.forward_length,
                    primer.forward_tm,
                    primer.forward_gc,
                    primer.forward_seq,
                )
            )
            outfh.write(
                "     REVERSE PRIMER  %-9d  %-3d  %.02f  %.02f  %s\n"
                % (
                    primer.reverse_start,
                    primer.reverse_length,
                    primer.reverse_tm,
                    primer.reverse_gc,
                    primer.reverse_seq,
                )
            )
            if hasattr(primer, "internal_start"):
                outfh.write(
                    "     INTERNAL OLIGO  "
                    + "%-9d  %-3d  %.02f  %.02f  %s\n"
                    % (
                        primer.internal_start,
                        primer.internal_length,
                        primer.internal_tm,
                        primer.internal_gc,
                        primer.internal_seq,
                    )
                )
            outfh.write("\n" * 3)


def __write_primers_json(primers, outfname):
    """Write Primer3 primer objects in JSON format.

    :param primers:  iterable of Primer3.Primers objects
    :param outfname:  path to output file
    """
    with open(outfname, "w") as ofh:
        json.dump(primers, ofh, cls=PrimersEncoder)


def __write_primers_bed(primers, outfname):
    """Write Primer3 primer objects in BED format.

    :param primers:  iterable of Primer3.Primers objects
    :param outfname:  path to output file
    """
    with open(outfname, "w") as outfh:
        sourceids = {}
        for primer in primers:
            try:
                source_id = sourceids[primer.source]  #  lazily acquire primer source ID
            except KeyError:
                source_id = sourceids.setdefault(
                    primer.source, load_fasta_id(primer.source)
                )
            outfh.write(
                "{}\t{}\t{}\t{}\n".format(
                    source_id,
                    primer.forward_start,
                    primer.reverse_start - 1,
                    primer.name,
                )
            )


def load_fasta_id(fname):
    """Return the identifier from the passed FASTA file.

    :param fname:  path to input FASTA file
    """
    with open(fname, "r") as ifh:
        return SeqIO.read(ifh, "fasta").id
