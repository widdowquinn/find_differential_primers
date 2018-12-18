# -*- coding: utf-8 -*-
"""Parser for pdp dedupe subcommand

(c) The James Hutton Institute 2017-2019

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

Copyright (c) 2017-2019 The James Hutton Institute
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

from diagnostic_primers.scripts import subcommands


def build(subparsers, parents=None):
    """Add parser for `filter` subcommand to subparsers

    This parser implements options for adding a filter to each genome to
    restrict primer design locations.
    """
    parser = subparsers.add_parser("filter", aliases=["filt"], parents=parents)
    parser.add_argument("outfilename", help="Path to write new configuration file")
    parser.add_argument(
        "--prodigal",
        dest="filt_prodigal",
        action="store_true",
        default=False,
        help="use prodigal to predict CDS, restrict primer design to these regions",
    )
    parser.add_argument(
        "--prodigaligr",
        dest="filt_prodigaligr",
        action="store_true",
        default=False,
        help="use prodigal to predict CDS, restrict primer design to intergenic regions",
    )
    parser.add_argument(
        "--alnvar",
        dest="filt_alnvar",
        action="store",
        default=None,
        help="use mummer to align genomes in the specified class, restrict primer design to aligned regions with sequence variability",
    )
    parser.add_argument(
        "--prodigal_exe",
        dest="filt_prodigal_exe",
        action="store",
        default="prodigal",
        help="path to prodigal executable",
    )
    parser.add_argument(
        "--outdir",
        dest="filt_outdir",
        action="store",
        default="filter_output",
        help="path to directory for filter program output",
    )
    parser.add_argument(
        "--suffix",
        dest="filt_suffix",
        action="store",
        type=str,
        default="filtered",
        help="suffix for filtered sequence file",
    )
    parser.add_argument(
        "--spacerlen",
        dest="filt_spacerlen",
        action="store",
        type=int,
        default=150,
        help="length of N spacer between regions",
    )
    parser.add_argument(
        "--flanklen",
        dest="filt_flanklen",
        action="store",
        type=int,
        default=0,
        help="length of feature flanking region to use when designing probes",
    )
    parser.add_argument(
        "--min_sim_error_count",
        dest="filt_minsecount",
        action="store",
        type=int,
        default=0,
        help="minimum number of similarity errors for retaining nucmer alignment with --alnvar",
    )
    parser.add_argument(
        "--min_sim_error_rate",
        dest="filt_minserate",
        action="store",
        type=float,
        default=0,
        help="minimum rate of similarity errors for retaining nucmer alignment with --alnvar",
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="filt_force",
        action="store_true",
        default=False,
        help="allow overwriting of filter output directory",
    )
    parser.set_defaults(func=subcommands.subcmd_filter)
