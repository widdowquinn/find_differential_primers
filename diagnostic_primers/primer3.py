#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""primer3.py

Code to conduct primer prediction with Primer3 (v2+)

(c) The James Hutton Institute 2016-2019
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

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

Copyright (c) 2016-2019 The James Hutton Institute

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

import json
import os
import shlex

from collections import namedtuple

from Bio import SeqIO

from diagnostic_primers import PDPException


# Define PDPPrimer3Exception
class PDPPrimer3Exception(PDPException):
    """General exception for interacting with Primer3 from PDP"""

    def __init__(self, msg="Error using Primer3"):
        PDPException.__init__(self, msg)


class Primer3Command(object):
    """Command-line for Primer3"""

    def __init__(self, cline, infile, outfile):
        self.cline = cline
        self.infile = infile
        self.outfile = outfile

    def __str__(self):
        return " ".join(self.cline)


def build_commands(collection, primer3_exe, primer3_dir, existingfiles, argdict=None):
    """Builds and returns a list of command-lines to run Primer3 (v2+)

    :param collection: PDPCollection object describing a set of input files
                    for primer design
    :param primer3_exe: path to primer3_core (v2+) executable
    :param primer3_dir: path to directory for writing input BoulderIO file,
                    and collecting primer3 output
    :param existingfiles:  iterable of existing output files to be reused
    :param argdict: dictionary of arguments from parser

    The commands will run on each sequence in the passed PDPCollection.

    For each input file in the PDPCollection, we need to generate an
    input file in BoulderIO format that will be provided as input
    to the primer3_core executable; we also need to provide an output
    directory to write the results to (the input file will be
    written here, also).
    """
    clines = []  # Holds command-lines

    # Ensure output directory exists
    os.makedirs(primer3_dir, exist_ok=True)

    for g in collection.data:
        stempath = os.path.split(os.path.splitext(g.seqfile)[0])
        stem = os.path.join(primer3_dir, stempath[-1])
        # Are we using the filter region information for this design?
        # Two ways we won't use the filter region: the --filter argument/
        # ep_filter argument is not set; the --filter argument/ep_filter
        # argument *is* set, but there's no filtered genome sequence file.
        if argdict["p3_filter"] and g.filtered_seqfile is not None:
            seqfile = g.filtered_seqfile
        else:
            seqfile = g.seqfile
        cline = build_command(primer3_exe, g.name, seqfile, stem, argdict)
        g.cmds["Primer3"] = cline
        if os.path.split(cline.outfile)[-1] not in existingfiles:
            clines.append(cline)
    return clines


def build_command(primer3_exe, seqname, seqfile, stem, argdict):
    """Build and return Primer3 (v2+) command line.

    :param primer3_exe: path to primer3_core executable (v2+)
    :param seqname: name to identify the sequence
    :param seqfile: path to input sequence file for primer design
    :param stem: filestem (with path) for input BoulderIO file and
                primer3 output
    :param argdict: dictionary of arguments from parser

    The command line will be structured as follows

    primer3_core \
        -output [path to output file] \
        [input file in BoulderIO format]

    The BoulderIO input file needs to be constructed, containing the
    sequence data, metadata and primer3 settings, for each run. This
    is written to the output directory in argdict["primer3_dir"],
    as is the output from primer3_core
    """
    # Create the input file
    try:
        infname = build_input_file(seqname, seqfile, stem, argdict)
    except Exception:
        raise PDPPrimer3Exception(
            "Error creating BoulderIO input file for {}".format(seqfile)
        )

    # Define path to output file, and return completed command-line
    ofname = stem + ".primer3"
    cline = [shlex.quote(_) for _ in [primer3_exe, "-output", ofname, infname]]
    return Primer3Command(cline, infname, ofname)


def build_input_file(seqname, seqfile, stem, argdict):
    """Return path to BoulderIO input file for Primer3

    :param seqname: name to identify the sequence
    :param seqfile: path to input sequence file for primer design
    :param stem: filestem (with path) for input BoulderIO file and
                primer3 output
    :param argdict: dictionary of arguments from parser

    Constructs a new Primer3 BoulderIO-format input file from the
    passed sequence file and argument dictionary, and writes it to
    the location indicated by stem, with file extension .boulder
    """
    ofname = stem + ".boulder"
    with open(ofname, "w") as ofh:
        # Define sequence name and template
        ofh.write("SEQUENCE_ID={}\n".format(seqname))
        inseq = str(SeqIO.read(seqfile, "fasta").seq)
        ofh.write("SEQUENCE_TEMPLATE={}\n".format(inseq))

        # Define primer design job
        ofh.write("PRIMER_TASK=generic\n")
        ofh.write("PRIMER_PICK_LEFT_PRIMER=1\n")
        ofh.write(
            "PRIMER_PICK_INTERNAL_OLIGO={}\n".format(
                1 if argdict["p3_hybridprobe"] else 0
            )
        )
        ofh.write("PRIMER_PICK_RIGHT_PRIMER=1\n")

        # Define number of primers to return
        ofh.write("PRIMER_NUM_RETURN={}\n".format(argdict["p3_numreturn"]))

        # Define primer sequence size
        ofh.write("PRIMER_OPT_SIZE={}\n".format(argdict["p3_osize"]))
        ofh.write("PRIMER_MIN_SIZE={}\n".format(argdict["p3_minsize"]))
        ofh.write("PRIMER_MAX_SIZE={}\n".format(argdict["p3_maxsize"]))

        # We need to penalise non-optimal size primers. This is done automatically
        # with ePrimer3, but not with Primer3 (v2+)
        ofh.write("PRIMER_PAIR_WT_PRODUCT_SIZE_LT={}\n".format(argdict["p3_wt_lt"]))
        ofh.write("PRIMER_PAIR_WT_PRODUCT_SIZE_GT={}\n".format(argdict["p3_wt_gt"]))

        # Define product size ranges
        ofh.write("PRIMER_PRODUCT_OPT_SIZE={}\n".format(argdict["p3_psizeopt"]))
        ofh.write(
            "PRIMER_PRODUCT_SIZE_RANGE={}-{}\n".format(
                argdict["p3_psizemin"], argdict["p3_psizemax"]
            )
        )

        # Define primer TM ranges
        ofh.write("PRIMER_OPT_TM={}\n".format(argdict["p3_opttm"]))
        ofh.write("PRIMER_MIN_TM={}\n".format(argdict["p3_mintm"]))
        ofh.write("PRIMER_MAX_TM={}\n".format(argdict["p3_maxtm"]))

        # Define primer GC content ranges
        ofh.write("PRIMER_OPT_GC_PERCENT={}\n".format(argdict["p3_ogcpercent"]))
        ofh.write("PRIMER_MIN_GC_PERCENT={}\n".format(argdict["p3_mingc"]))
        ofh.write("PRIMER_MAX_GC_PERCENT={}\n".format(argdict["p3_maxgc"]))

        # Define internal oligo size
        ofh.write("PRIMER_INTERNAL_OPT_SIZE={}\n".format(argdict["p3_osizeopt"]))
        ofh.write("PRIMER_INTERNAL_MIN_SIZE={}\n".format(argdict["p3_ominsize"]))
        ofh.write("PRIMER_INTERNAL_MAX_SIZE={}\n".format(argdict["p3_omaxsize"]))

        # Define internal oligo TM ranges
        ofh.write("PRIMER_INTERNAL_OPT_TM={}\n".format(argdict["p3_otmopt"]))
        ofh.write("PRIMER_INTERNAL_MIN_TM={}\n".format(argdict["p3_otmmin"]))
        ofh.write("PRIMER_INTERNAL_MAX_TM={}\n".format(argdict["p3_otmmax"]))

        # Define internal oligo GC content ranges
        ofh.write("PRIMER_INTERNAL_OPT_GC_PERCENT={}\n".format(argdict["p3_ogcopt"]))
        ofh.write("PRIMER_INTERNAL_MIN_GC_PERCENT={}\n".format(argdict["p3_ogcmin"]))
        ofh.write("PRIMER_INTERNAL_MAX_GC_PERCENT={}\n".format(argdict["p3_ogcmax"]))

        # Define maximum acceptable mononucleotide repeat
        ofh.write("PRIMER_MAX_POLY_X={}\n".format(argdict["p3_maxpolyx"]))

        # If an alternative path to thermodynamic parameters is provided, use it
        if argdict["p3_param_path"] is not None:
            if argdict["p3_param_path"][-1] != "/":
                argdict["p3_param_path"] += "/"
            ofh.write(
                "PRIMER_THERMODYNAMIC_PARAMETERS_PATH={}\n".format(
                    argdict["p3_param_path"]
                )
            )

        # Terminate input file
        ofh.write("=\n")

    return ofname
