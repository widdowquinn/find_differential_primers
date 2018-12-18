# -*- coding: utf-8 -*-
"""Parser for pdp eprimer3 subcommand

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
    """Add parser for `eprimer3` subcommand to subparsers

    This parser implements options for controlling primer creation with the
    EMBOSS ePrimer3 tool.
    """
    parser = subparsers.add_parser("eprimer3", aliases=["e3"], parents=parents)
    # Primer prediction options - subcommand eprimer3
    parser.add_argument("outfilename", help="Path to write new configuration file")
    parser.add_argument(
        "--eprimer3",
        dest="eprimer3_exe",
        action="store",
        default="eprimer3",
        help="path to ePrimer3 executable",
    )
    parser.add_argument(
        "--outdir",
        dest="eprimer3_dir",
        action="store",
        default="eprimer3",
        help="path to directory for ePrimer3 output",
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="eprimer3_force",
        action="store_true",
        default=False,
        help="Overwrite old ePrimer3 output",
    )
    parser.add_argument(
        "--filter",
        dest="ep_filter",
        action="store_true",
        default=False,
        help="use the filtered_seqfile to design primer sets",
    )
    parser.add_argument(
        "--numreturn",
        dest="ep_numreturn",
        action="store",
        default=10,
        type=int,
        help="number of primers to return",
    )
    parser.add_argument(
        "--osize",
        dest="ep_osize",
        action="store",
        default=20,
        type=int,
        help="optimal size for primer oligo",
    )
    parser.add_argument(
        "--minsize",
        dest="ep_minsize",
        action="store",
        default=18,
        type=int,
        help="minimum size for primer oligo",
    )
    parser.add_argument(
        "--maxsize",
        dest="ep_maxsize",
        action="store",
        default=22,
        type=int,
        help="maximum size for primer oligo",
    )
    parser.add_argument(
        "--opttm",
        dest="ep_opttm",
        action="store",
        default=59,
        type=int,
        help="optimal Tm for primer oligo",
    )
    parser.add_argument(
        "--mintm",
        dest="ep_mintm",
        action="store",
        default=58,
        type=int,
        help="minimum Tm for primer oligo",
    )
    parser.add_argument(
        "--maxtm",
        dest="ep_maxtm",
        action="store",
        default=60,
        type=int,
        help="maximum Tm for primer oligo",
    )
    parser.add_argument(
        "--ogcpercent",
        dest="ep_ogcpercent",
        action="store",
        default=55,
        type=int,
        help="optimal %%GC for primer oligo",
    )
    parser.add_argument(
        "--mingcpercent",
        dest="ep_mingc",
        action="store",
        default=30,
        type=int,
        help="minimum %%GC for primer oligo",
    )
    parser.add_argument(
        "--maxgcpercent",
        dest="ep_maxgc",
        action="store",
        default=80,
        type=int,
        help="maximum %%GC for primer oligo",
    )
    parser.add_argument(
        "--psizeopt",
        dest="ep_psizeopt",
        action="store",
        default=100,
        type=int,
        help="optimal size of amplified region",
    )
    parser.add_argument(
        "--psizemin",
        dest="ep_psizemin",
        action="store",
        default=50,
        type=int,
        help="minimum size of amplified region",
    )
    parser.add_argument(
        "--psizemax",
        dest="ep_psizemax",
        action="store",
        default=150,
        type=int,
        help="maximum size of amplified region",
    )
    parser.add_argument(
        "--maxpolyx",
        dest="ep_maxpolyx",
        action="store",
        default=3,
        type=int,
        help="maximum run of repeated nucleotides " + "in primer",
    )
    parser.add_argument(
        "--hybridprobe",
        dest="ep_hybridprobe",
        action="store_true",
        default=False,
        help="design a reporter oligo",
    )
    parser.add_argument(
        "--oligoosize",
        dest="ep_osizeopt",
        action="store",
        default=20,
        type=int,
        help="optimal size for internal oligo",
    )
    parser.add_argument(
        "--oligominsize",
        dest="ep_ominsize",
        action="store",
        default=13,
        type=int,
        help="minimum size for internal oligo",
    )
    parser.add_argument(
        "--oligomaxsize",
        dest="ep_omaxsize",
        action="store",
        default=30,
        type=int,
        help="maximum size for internal oligo",
    )
    parser.add_argument(
        "--oligootm",
        dest="ep_otmopt",
        action="store",
        default=69,
        type=int,
        help="optimal Tm for internal oligo",
    )
    parser.add_argument(
        "--oligomintm",
        dest="ep_otmmin",
        action="store",
        default=68,
        type=int,
        help="minimum Tm for internal oligo",
    )
    parser.add_argument(
        "--oligomaxtm",
        dest="ep_otmmax",
        action="store",
        default=70,
        type=int,
        help="maximum Tm for internal oligo",
    )
    parser.add_argument(
        "--oligoogcpercent",
        dest="ep_ogcopt",
        action="store",
        default=55,
        type=int,
        help="optimal %%GC for internal oligo",
    )
    parser.add_argument(
        "--oligomingcpercent",
        dest="ep_ogcmin",
        action="store",
        default=30,
        type=int,
        help="minimum %%GC for internal oligo",
    )
    parser.add_argument(
        "--oligomaxgcpercent",
        dest="ep_ogcmax",
        action="store",
        default=80,
        type=int,
        help="maximum %%GC for internal oligo",
    )
    # Commented out until Biopython command line code catches up with new
    # EMBOSS options
    #    parser.add_argument("--oligomaxpolyx", dest="ep_opolymax",
    #                                 action="store",
    #                                 default=3, type=int,
    #                                 help="maximum run of repeated nucleotides " +
    #                                 "in internal primer")
    parser.set_defaults(func=subcommands.subcmd_eprimer3)
