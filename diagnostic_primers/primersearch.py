# -*- coding: utf-8 -*-
"""primersearch.py

Code to conduct primer crosshybridisation tests with primersearch

(c) The James Hutton Institute 2016-2019

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
import re

from collections import defaultdict, namedtuple

from Bio import SeqIO
from Bio.Emboss.Applications import PrimerSearchCommandline
from pybedtools import BedTool

from diagnostic_primers import load_primers, write_primers


# Convenience struct for amplimer
Amplimer = namedtuple("Amplimer", "seq start end")


def build_commands(collection, primersearch_dir, mismatchpercent, existingfiles):
    """Build and return a list of command-lines to run primersearch.

    :param collection:  PDPCollection describing analysis inputs
    :param primersearch_dir:  path to primersearch output
    :param mismatchpercent:  allowed 'wobble' for primers
    """
    clines = []  # holds command lines

    # Make sure output directory exists
    os.makedirs(primersearch_dir, exist_ok=True)

    # Build command-lines for each input primer set
    # We generate a collection of target sequences from each PDPData
    # object, then loop over all PDPData objects, and build a command
    # line comparing that object's predicted primers to all other
    # target sequences.
    # Primersearch - bafflingly - doesn't accept EMBOSS' ePrimer3
    # format, so we need to write primers out as 3-column TSV
    for dat in collection.data:
        # Get the primers from the JSON file and write them
        # to the output directory
        primerpath = os.path.join(
            primersearch_dir, "{}_primers.primertab".format(dat.name)
        )
        primers = load_primers(dat.primers, "json")
        write_primers(primers, primerpath, "tsv")
        # Create a dictionary to hold target names, to be written
        # to a JSON file, and the path added to the PDPData object
        psdict = {"query": dat.name, "primers": primerpath}
        for _ in collection.data:
            # if dat.name != tgtname:  # Don't query against self-genome
            # Name for output file is built from the PDPData
            # query/target object names
            outstem = os.path.join(
                primersearch_dir, "{}_ps_{}.primersearch".format(dat.name, _.name)
            )
            # Add the output file to the PDPData primersearch attr
            psdict[_.name] = outstem
            # Generate the primersearch cmd-line
            cline = build_command(primerpath, _.seqfile, outstem, mismatchpercent)
            # If there are no primers, attempting to run the primersearch
            # command will give a non-zero exit code, and fail tests. In this
            # instance we write a blank primersearch "output file" and don't
            # append the command-line to the list which will be returned.
            if not primers:
                with open(outstem, "w") as ofh:
                    ofh.write("")
            elif not os.path.split(outstem)[-1] in existingfiles:
                clines.append(cline)
        # Write primersearch output JSON file and add to PDPData object
        psjson = os.path.join(primersearch_dir, f"{dat.name}_primersearch.json")
        with open(psjson, "w") as ofh:
            json.dump(psdict, ofh, sort_keys=True)
        dat.primersearch = psjson
    return clines


def build_command(primerfile, seqfile, filestem, mismatchpercent):
    """Return a single primersearch command line.

    The returned command line queries a single primer set against a
    single genome sequence/input sequence.

    :param primerfile:          path to input primer file (ePrimer3 format)
    :param seqfile:             path to target sequence
    :param filestem:            path to output file
    :param argdict:             optional arguments to primersearch
    :param mismatchpercent:     allowed 'wobble' for primers
    """
    # Build command-line with defaults
    cline = PrimerSearchCommandline()
    cline.auto = True
    cline.seqall = seqfile
    cline.infile = primerfile
    cline.outfile = filestem
    cline.mismatchpercent = mismatchpercent
    return cline


class PrimerSearchRecord:
    """Container for single PrimerSearch record

    This will contain the entire data from a PrimerSearch record
    """

    def __init__(self, name):
        self._name = str(name)
        self._amplimers = []

    def add_amplimer(self, amplimer):
        """Add a PrimerSearchAmplimer to the list of amplimers"""
        self._amplimers.append(amplimer)

    @property
    def name(self):
        """Get name of PrimerSearchRecord"""
        return self._name

    @property
    def amplimers(self):
        """Get amplimers from PrimerSearchRecord"""
        return self._amplimers[:]

    def __str__(self):
        outstr = ["\nPrimer name %s" % self.name]
        for idx, amp in enumerate(self.amplimers):
            outstr += [
                "Amplimer %d" % (idx + 1),
                "\tSequence: %s" % amp.sequence,
                "\tAmplimer length: %d bp" % len(amp),
            ]
        return "\n".join(outstr) + "\n"


class PrimerSearchAmplimer:
    """Container for single PrimerSearch amplimer

    This will contain the entire data from a PrimerSearch record amplimer

    TODO: Tidy the forward/reverse start and maybe preprocess
    """

    def __init__(self, name):
        self._name = str(name)
        self._seqname = ""
        self._len = None
        self.primer_name = None
        self.target = None
        self.target_fasta_id = None
        self._fwd = None  # namedtuple of forward primer
        self._rev = None  # namedtuple of reverse primer

    @property
    def name(self):
        """Get name of PrimerSearchAmplimer"""
        return self._name

    @property
    def fwd(self):
        """Forward primer"""
        return self._fwd

    @fwd.setter
    def fwd(self, val):
        """Set forward primer"""
        self._fwd = val

    @property
    def rev(self):
        """Forward primer"""
        return self._rev

    @rev.setter
    def rev(self, val):
        """Set forward primer"""
        self._rev = val

    @property
    def sequence(self):
        """Get PrimerSearchAmplimer sequence name"""
        return self._seqname

    @sequence.setter
    def sequence(self, val):
        """Set name of originating sequence"""
        self._seqname = str(val)

    @property
    def length(self):
        """Get length of PrimerSearchAmplimer"""
        return self._len

    @length.setter
    def length(self, val):
        """Set length of PrimerSearchAmplimer"""
        self._len = int(val)

    def __len__(self):
        return self.length

    def __str__(self):
        return "Amplimer: {}, length: {}".format(self.name, self.length)


class AmplimersEncoder(json.JSONEncoder):
    """JSON encoder for PrimerSearchAmplimer objects."""

    def default(self, o):  # pylint: disable=E0202
        if not isinstance(o, PrimerSearchAmplimer):
            return json.JSONEncoder.default(self, o)

        # Convert complex PrimerSearchAmplimer object to serialisable dictionary
        # and return
        return o.__dict__


class PDPGenomeAmplicons:
    """Collection of primersearch amplimers.

    The genomes dictionary contains a set of names of amplimers from that genome,
    keyed by genome name. The amplimer dictionary contains PrimerSearchAmplimer
    objects, keyed by name.
    """

    def __init__(self, name):
        self.name = str(name)
        self._targets = defaultdict(list)

    def add_amplimer(self, amplimer, targetname):
        """Add an amplimer to the collection, grouped by genome name

        amplimer          PrimerSearchAmplimer object
        target            name of the target amplified genome
        """
        self._targets[targetname].append(amplimer)

    def from_json(self, filename):
        """Load data from JSON format file

        filename         path to JSON data file

        This overwrites the objects current data
        """
        with open(filename, "r") as ifh:
            data = json.load(ifh)
        for target, amplimers in data.items():
            for ampdata in amplimers:
                newamp = PrimerSearchAmplimer(ampdata["_name"])
                for key, val in ampdata.items():
                    if key == "_fwd":  # forward amplimer
                        newamp.fwd = Amplimer(*val)
                    elif key == "_rev":  # reverse amplimer
                        newamp.rev = Amplimer(*val)
                    else:
                        newamp.__setattr__(key, val)
                self.add_amplimer(newamp, target)

    def filter_primers(self, primers):
        """Return PDPGenomeAmplicons object containing only amplimers in primers"""
        obj = PDPGenomeAmplicons(self.name + "_filtered")
        for target in self._targets:
            for amplimer in self._targets[target]:
                if amplimer.primer_name in primers:
                    obj.add_amplimer(amplimer, target)
        return obj

    def get_target_amplimers(self, target):
        """Return the amplimer collection for the named target."""
        return self._targets[target]

    def split_on_targets(self):
        """Return a list of new PDPGenomeAmplicons objects for each target."""
        split_list = []
        for key in self._targets:
            obj = PDPGenomeAmplicons(key)
            obj.targets = {key: self._targets[key]}
            split_list.append(obj)
        return split_list

    def write_json(self, outfilename):
        """Write the object to a JSON format file

        outfilename       path to output JSON file

        Writes serialised objects to JSON file
        """
        with open(outfilename, "w") as ofh:
            json.dump(self._targets, ofh, sort_keys=True, cls=PDPGenomeAmpliconsEncoder)

    def write_bed(self, outdir):
        """Write one BED file per target describing each genomes amplicons

        outdir            path to directory for BED output
        """
        for target, amplimers in self._targets.items():
            ofpath = os.path.join(outdir, target + "_amplicons.bed")
            # Create list of (target, start, end, primer_name) tuples
            regions = BedTool(
                [
                    (amp.target_fasta_id, amp.fwd.start, amp.rev.end, amp.primer_name)
                    for amp in amplimers
                ]
            )
            regions.saveas(ofpath)

    def write_target_bed(self, target, outfname):
        """Write single BED file for one target.

        target         name of target genome to write out
        outfname       path to output file
        """
        amplimers = self._targets[target]
        regions = BedTool(
            [
                (amp.target_fasta_id, amp.fwd.start, amp.rev.end, amp.primer_name)
                for amp in amplimers
            ]
        )
        regions.saveas(outfname)

    @property
    def targets(self):
        """Return names of amplimers"""
        return sorted(list(self._targets.keys()))

    @targets.setter
    def targets(self, tgtdict):
        """Set targets dictionary"""
        self._targets = tgtdict


class PDPGenomeAmpliconsEncoder(json.JSONEncoder):

    """JSON encoder for PDPGenomeAmplicons objects"""

    def default(self, o):  # pylint: disable=E0202
        if isinstance(o, PrimerSearchAmplimer):
            encoder = AmplimersEncoder()
            return encoder.default(o)
        if not isinstance(o, PDPGenomeAmplicons):
            return super(PDPGenomeAmpliconsEncoder, self).default(o)

        # Convert PDPDiagnosticPrimers object to serialisable dictionary
        # and return
        return o.__dict__


def parse_output(filename, genomepath):
    """Return the contents of a PrimerSearch output file

    filename         path to the primersearch output file
    genomepath       path to the target genome from primersearch analysis

    TODO: This is a hacky wee function - Scanner/Consumer would be better
          While we're developing, this will do. But we could cope with a
          more complete model of the data.
    """
    records = []
    with open(genomepath, "r") as ifh:
        target = SeqIO.read(ifh, "fasta")
    with open(filename, "r") as ifh:
        record = None
        for line in ifh:
            if line.startswith("Primer name"):  # Start of record
                if record is not None:
                    records.append(record)
                rname = line.split("Primer name")[-1].strip()
                record = PrimerSearchRecord(rname)
            if line.startswith("Amplimer"):
                aname = line.strip()
                amplimer = PrimerSearchAmplimer(aname)
                amplimer.target_fasta_id = target.id
                amplimer.primer_name = rname
                amplimer.target = genomepath
                record.add_amplimer(amplimer)
            if line.strip().startswith("Sequence"):
                seqname = line.split("Sequence:")[-1].strip()
                amplimer.sequence = seqname
            if line.strip().startswith("Amplimer length"):
                alen = int(line.split("Amplimer length:")[-1].strip().split()[0])
                amplimer.length = alen
            # PrimerSearch output records matches in the direction of strandedness
            # of the target genome, not the direction of "forward" or "reverse"
            # primers. We deal with this elsewhere to preserve the PrimerSearch
            # output file data.
            if "forward strand" in line:
                fseq = line.strip().split()[0]
                fstart = int(re.search(r"(?<=at )[0-9]*", line).group())
                amplimer.fwd = Amplimer(fseq, fstart, fstart + len(fseq))
                # amplimer.forward_seq = line.strip().split()[0]
                # amplimer.forward_start = int(re.search(r"(?<=at )[0-9]*", line).group())
                # amplimer.forward_end = amplimer.forward_start + len(
                #     amplimer.forward_seq
                # )
            if "reverse strand" in line:
                rseq = line.strip().split()[0]
                rend = (
                    len(target) - int(re.search(r"(?<=at \[)[0-9]*", line).group()) + 1
                )
                amplimer.rev = Amplimer(rseq, rend - len(rseq), rend)
                # amplimer.reverse_seq = line.strip().split()[0]
                # amplimer.reverse_end = (
                #     len(target) - int(re.search(r"(?<=at \[)[0-9]*", line).group()) + 1
                # )
                # amplimer.reverse_start = amplimer.reverse_end - len(
                #     amplimer.reverse_seq
                # )
        if record is not None:
            records.append(record)
    return records


def load_collection_amplicons(coll):
    """Return PDPGenomeAmplicons object for a passed PDPCollection

    coll         PDPCollection object with primersearch output

    Parses the primersearch output for all input genomes in a collection,
    returning a PDPGenomeAmplicons object describing the target genomes and
    the amplimers that are amplified from them.
    """
    # For each input genome...
    # - populate a dictionary of paths to each input genome, keyed by name
    genomepaths = {}  # dictionary of paths to each genome file, keyed by name
    for item in coll.data:
        genomepaths[item.name] = item.seqfile

    # - get the (forward) PrimerSearch results of its child primers amplifying the
    #   other genomes
    targetamplicons = PDPGenomeAmplicons(coll.name)
    for item in coll.data:
        # Open the 'forward' PrimerSearch output
        with open(item.primersearch, "r") as psfh:
            psdata = json.load(psfh)
            # Populate collection of amplimers for each target genome
            for targetname in [
                _ for _ in psdata.keys() if _ not in ("primers", "query")
            ]:
                psresult = parse_output(psdata[targetname], genomepaths[targetname])
                for amplicon in psresult:
                    for amplimer in amplicon.amplimers:
                        targetamplicons.add_amplimer(amplimer, targetname)
    return targetamplicons
