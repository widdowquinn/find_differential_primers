#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""subcmd_dedupe.py

Provides the dedupe subcommand for pdp

(c) The James Hutton Institute 2018-2019

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

Copyright (c) 2018-2019 The James Hutton Institute
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

from tqdm import tqdm

from diagnostic_primers import load_primers, write_primers
from diagnostic_primers.scripts.tools import load_config_json


def ensure_path_to(fname):
    """If the path to the passed filename doesn't exist, it is created."""
    absfname = os.path.abspath(fname)
    pathto = os.path.dirname(absfname)
    os.makedirs(pathto, exist_ok=True)
    return absfname


def subcmd_dedupe(args, logger):
    """Remove duplicated primer sets from a primer config file

    This requires us to

    - read in the config file and identify primer sets
    - hash/uniquely identify primer sets (as left/right primer sequences)
    - create a new config file with only one representative primer per
      unique primer set, recording the synonyms
    """
    logger.info("Deduplicating primers")

    # Load the JSON config file (any stage post-ePrimer3)
    coll = load_config_json(args, logger)

    # Inspect each primer set in turn, and key by tuple of (FWD, REV) primer sequence
    # If the key has been seen before
    seen = set()
    removed = 0
    kept = 0
    pbar = tqdm(coll.data, desc="deduplicating primers", disable=args.disable_tqdm)
    for cdata in pbar:
        primers = load_primers(cdata.primers, "json")
        outpfname = os.path.splitext(cdata.primers)[0] + "_deduped"
        if args.dd_dedupedir is not None:
            outpfname = os.path.join(args.dd_dedupedir, os.path.split(outpfname)[-1])
            ensure_path_to(outpfname)
        nonredundant = []
        for primer in primers:
            key = (primer.forward_seq, primer.reverse_seq)
            if key in seen:
                # Remove primer from primer list and write out new file when we're done
                removed += 1
            else:
                kept += 1
                nonredundant.append(primer)
                seen.add(key)
        # write deduplicated primers
        write_primers(nonredundant, outpfname + ".json", "json")
        write_primers(nonredundant, outpfname + ".bed", "bed")
        cdata.primers = (
            outpfname + ".json"
        )  # update PDPCollection with new primer location
    logger.info("%d primer sets were duplicated (%d kept)", removed, kept)
    logger.info("Writing deduped config file to %s", args.outfilename)
    coll.write_json(args.outfilename)  # write updated PDPCollection

    return 0
