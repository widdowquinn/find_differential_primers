# -*- coding: utf-8 -*-
"""Parser for pdp primer3 subcommand

(c) The James Hutton Institute 2019

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

Copyright (c) 2019 The James Hutton Institute
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
    """Add parser for `primer3` subcommand to subparsers

    This parser implements options for controlling primer creation with the
    popular Primer3 tool (v2+).

    Initially, we're using the same command-line options as for the eprimer3
    subcommand.
    """
    parser = subparsers.add_parser("primer3", aliases=["p3"], parents=parents)
    # Primer prediction options - subcommand eprimer3
    parser.add_argument("outfilename", help="Path to write new configuration file")
    parser.add_argument(
        "--primer3",
        dest="primer3_exe",
        action="store",
        default="primer3_core",
        help="path to Primer3 (v2+) executable",
    )
    parser.add_argument(
        "--therm_param_path",
        dest="p3_param_path",
        action="store",
        default=None,
        help="PRIMER_THERMODYNAMIC_PARAMETERS_PATH setting",
    )
    parser.add_argument(
        "--outdir",
        dest="primer3_dir",
        action="store",
        default="primer3",
        help="path to directory for Primer3 (v2+) output",
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="primer3_force",
        action="store_true",
        default=False,
        help="Overwrite old ePrimer3 output",
    )
    parser.add_argument(
        "--filter",
        dest="p3_filter",
        action="store_true",
        default=False,
        help="use the filtered_seqfile to design primer sets",
    )
    parser.add_argument(
        "--numreturn",
        dest="p3_numreturn",
        action="store",
        default=10,
        type=int,
        help="number of primers to return",
    )
    parser.add_argument(
        "--osize",
        dest="p3_osize",
        action="store",
        default=20,
        type=int,
        help="optimal size for primer oligo",
    )
    parser.add_argument(
        "--minsize",
        dest="p3_minsize",
        action="store",
        default=18,
        type=int,
        help="minimum size for primer oligo",
    )
    parser.add_argument(
        "--maxsize",
        dest="p3_maxsize",
        action="store",
        default=22,
        type=int,
        help="maximum size for primer oligo",
    )
    parser.add_argument(
        "--wt_lt",
        dest="p3_wt_lt",
        action="store",
        default=1,
        type=float,
        help="penalty for products less than optimal size",
    )
    parser.add_argument(
        "--wt_gt",
        dest="p3_wt_gt",
        action="store",
        default=1,
        type=float,
        help="penalty for products greater than optimal size",
    )
    parser.add_argument(
        "--opttm",
        dest="p3_opttm",
        action="store",
        default=59,
        type=int,
        help="optimal Tm for primer oligo",
    )
    parser.add_argument(
        "--mintm",
        dest="p3_mintm",
        action="store",
        default=58,
        type=int,
        help="minimum Tm for primer oligo",
    )
    parser.add_argument(
        "--maxtm",
        dest="p3_maxtm",
        action="store",
        default=60,
        type=int,
        help="maximum Tm for primer oligo",
    )
    parser.add_argument(
        "--ogcpercent",
        dest="p3_ogcpercent",
        action="store",
        default=55,
        type=int,
        help="optimal %%GC for primer oligo",
    )
    parser.add_argument(
        "--mingcpercent",
        dest="p3_mingc",
        action="store",
        default=30,
        type=int,
        help="minimum %%GC for primer oligo",
    )
    parser.add_argument(
        "--maxgcpercent",
        dest="p3_maxgc",
        action="store",
        default=80,
        type=int,
        help="maximum %%GC for primer oligo",
    )
    parser.add_argument(
        "--psizeopt",
        dest="p3_psizeopt",
        action="store",
        default=100,
        type=int,
        help="optimal size of amplified region",
    )
    parser.add_argument(
        "--psizemin",
        dest="p3_psizemin",
        action="store",
        default=50,
        type=int,
        help="minimum size of amplified region",
    )
    parser.add_argument(
        "--psizemax",
        dest="p3_psizemax",
        action="store",
        default=150,
        type=int,
        help="maximum size of amplified region",
    )
    parser.add_argument(
        "--maxpolyx",
        dest="p3_maxpolyx",
        action="store",
        default=3,
        type=int,
        help="maximum run of repeated nucleotides " + "in primer",
    )
    parser.add_argument(
        "--hybridprobe",
        dest="p3_hybridprobe",
        action="store_true",
        default=False,
        help="design a reporter oligo",
    )
    parser.add_argument(
        "--oligoosize",
        dest="p3_osizeopt",
        action="store",
        default=20,
        type=int,
        help="optimal size for internal oligo",
    )
    parser.add_argument(
        "--oligominsize",
        dest="p3_ominsize",
        action="store",
        default=13,
        type=int,
        help="minimum size for internal oligo",
    )
    parser.add_argument(
        "--oligomaxsize",
        dest="p3_omaxsize",
        action="store",
        default=30,
        type=int,
        help="maximum size for internal oligo",
    )
    parser.add_argument(
        "--oligootm",
        dest="p3_otmopt",
        action="store",
        default=69,
        type=int,
        help="optimal Tm for internal oligo",
    )
    parser.add_argument(
        "--oligomintm",
        dest="p3_otmmin",
        action="store",
        default=68,
        type=int,
        help="minimum Tm for internal oligo",
    )
    parser.add_argument(
        "--oligomaxtm",
        dest="p3_otmmax",
        action="store",
        default=70,
        type=int,
        help="maximum Tm for internal oligo",
    )
    parser.add_argument(
        "--oligoogcpercent",
        dest="p3_ogcopt",
        action="store",
        default=55,
        type=int,
        help="optimal %%GC for internal oligo",
    )
    parser.add_argument(
        "--oligomingcpercent",
        dest="p3_ogcmin",
        action="store",
        default=30,
        type=int,
        help="minimum %%GC for internal oligo",
    )
    parser.add_argument(
        "--oligomaxgcpercent",
        dest="p3_ogcmax",
        action="store",
        default=80,
        type=int,
        help="maximum %%GC for internal oligo",
    )
    parser.set_defaults(func=subcommands.subcmd_primer3)
