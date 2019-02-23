#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""config.py

Provides PDPCollection and PDPData classes

(c) The James Hutton Institute 2017-2018

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

Copyright (c) 2017-2018 The James Hutton Institute
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

import csv
import json
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqRecord import SeqRecord

from pybedtools import BedTool

from diagnostic_primers import PDPException


# Exception of syntax error in config file
class ConfigSyntaxError(Exception):
    """Custom exception for parsing config files."""

    def __init__(self, message):
        Exception.__init__(self, message)


# Exception for PDP data object
class PDPCollectionException(PDPException):
    """Exception thrown for problems with PDPCollection objects"""

    def __init__(self, msg="Problem with PDPCollection"):
        PDPException.__init__(self, msg)


class PDPEncoder(json.JSONEncoder):
    """JSON encoder for PDP objects."""

    def default(self, obj):
        if not isinstance(obj, PDPData):
            return json.JSONEncoder.default(self, obj)

        # Convert complex PDPData object to serialisable dictionary and return
        objdict = {
            "target_amplicons": obj.target_amplicons,
            "name": obj.name,
            "groups": obj.groups,
            "seqfile": obj.seqfile,
            "filtered_seqfile": obj.filtered_seqfile,
            "features": obj.features,
            "primers": obj.primers,
            "primersearch": obj.primersearch,
            "filestem": obj.filestem,
        }

        return objdict


# Class that contains PDPData objects and interfaces with config files
class PDPCollection(object):
    """Container for PDPData objects and config file interface."""

    def __init__(self, name="pdp.py"):
        self.name = str(name)
        self._data = {}  # PDPData objects, keyed by .name attribure

    def from_tab(self, filename):
        """Load data from tab-format config file.

        Tab-format config files are expected to contain three or four columns
        only:

        <name>\t<groups>\t<sequence>
        <name>\t<groups>\t<sequence>\t<features>

        describing the unique name for an input, the groups to which it
        belongs, and the sequence file corresponding.

        The fourth column is optional and describes regions that can be
        excluded from primer design, or that primers will be restricted to.

        Comment lines begin with a hash and are ignored.
        """
        with open(filename, newline="") as ifh:
            reader = csv.reader(ifh, delimiter="\t")
            for row in reader:
                if len(row) and not row[0].startswith("#"):
                    self.add_data(*self.__parse_row(row))

    def from_json(self, filename):
        """Load data from JSON format config file.

        JSON format config files describe arrays of input genomes/sequences,
        where each array entry is a JSON object with the following keys:

        'name', 'groups', 'seqfile', 'features', 'primers'

        These are used directly to populate the collection's PDPData objects.
        """
        with open(filename, "r") as ifh:
            data = json.load(ifh)
        try:
            for item in data:
                # Not all config files have all data, so we make some fields
                # conditional
                if "primersearch" in item:
                    primersearch_val = item["primersearch"]
                else:
                    primersearch_val = None
                if "filtered_seqfile" in item:
                    filtered_seqval = item["filtered_seqfile"]
                else:
                    filtered_seqval = None
                self.add_data(
                    item["name"],
                    item["groups"],
                    item["seqfile"],
                    filtered_seqval,
                    item["features"],
                    item["primers"],
                    primersearch_val,
                    item["target_amplicons"],
                )
        except KeyError:
            raise PDPCollectionException(
                "Expected JSON field is missing, for instantiation"
            )

    def add_data(
        self,
        name=None,
        groups=None,
        seqfile=None,
        filtered_seqfile=None,
        features=None,
        primers=None,
        primersearch=None,
        target_amplicons=None,
    ):
        """Create a new PDPData object from passed info and add to collection.

        name         -    unique identifier for object
        groups       -    list of groups to which the object belongs
        seqfile      -    path to sequence file
        filtered_seqfile - path to sequence file composed of regions in features
        features     -    path to regions for inclusion/exclusion
        primers      -    path to primers in JSON format
        primersearch -    path to primersearch results in JSON format
        target_amplicons - path to target_amplicons JSON file
        """
        self._data[name] = PDPData(
            name,
            groups,
            seqfile,
            filtered_seqfile,
            features,
            primers,
            primersearch,
            target_amplicons,
        )

    def write_json(self, outfilename):
        """Write the Collection data contents to JSON format config file.

        outfilename  -    path to JSON config file

        Writes an array of serialised PDPData objects to JSON configuration
        file.
        """
        with open(outfilename, "w") as ofh:
            json.dump(self.data, ofh, sort_keys=True, cls=PDPEncoder)

    def write_tab(self, outfilename):
        """Write the Collection data contents to .tab format config file.

        outfilename -     path to .tab config file

        Writes a table of tab-separated columns describing the collection
        """
        with open(outfilename, "w") as ofh:
            for inseq in self.data:
                ofh.write(
                    "\t".join(
                        [
                            inseq.name,
                            ",".join(inseq.groups),
                            inseq.seqfile,
                            inseq.features or "-",
                        ]
                    )
                    + "\n"
                )

    @staticmethod
    def __parse_row(row):
        """Parse row from a tab-format config file to list.

        For tab-format configs, we expect three or four columns. We throw
        an error for anything else.
        """
        if len(row) < 3 or len(row) > 4:
            raise ConfigSyntaxError(
                "Row must contain 3 or 4 columns, "
                + "got %d at %s" % (len(row), "\t".join(row))
            )

        # Column 2 should be a comma-separated list of strings, defining
        # groups to which an input belongs:
        groups = [e.strip() for e in row[1].split(",")]
        row[1] = groups

        # Substitute None for '-' placeholders
        return [None if e == "-" else e for e in row]

    @property
    def data(self):
        """List of contained PDPData objects."""
        return [v for (k, v) in sorted(self._data.items())]

    @property
    def groups(self):
        """Groups present in the contained PDPData objects"""
        groups = set()
        for d in self.data:
            groups = groups.union(set(d.groups))
        return sorted(list(groups))

    def get_groupmembers(self, val):
        """PDPData objects having the passed group identity"""
        if val not in self.groups:
            raise PDPCollectionException("Group not found in PDPCollection")
        return [d for d in self.data if val in d.groups]

    def __getitem__(self, key):
        """Return the PDPData object with the named key"""
        return self._data[key]

    def __len__(self):
        """Return number of items in collection."""
        return len(self._data)


# Class defining paths to data for inputs and methods to operate on it
class PDPData(object):
    """Container for input sequence data and operations on that data."""

    def __init__(
        self,
        name,
        groups,
        seqfile,
        filtered_seqfile,
        features,
        primers,
        primersearch,
        target_amplicons,
    ):
        self._name = ""  # Set up private attributes
        self._groups = set()
        self._seqfile = None
        self._filtered_seqfile = None
        self._features = None
        self._primers = None
        self._primersearch = None
        self._filestem = None
        self._target_amplicons = None
        self.cmds = {}  # command-lines used to generate this object
        self.name = name  # Populate attributes
        self.groups = groups
        self.seqfile = seqfile
        self.filtered_seqfile = filtered_seqfile
        self.features = features
        self.primers = primers
        self.primersearch = primersearch
        if target_amplicons is not None:
            self.target_amplicons = target_amplicons
        # Useful values
        self.spacer = "NNNNNCATCCATTCATTAATTAATTAATGAATGAATGNNNNN"
        self.ambiguities = re.compile("[BDHKMRSVWY]")

    def stitch(self, outdir=None):
        """Stitch sequences in the sequence file, if necessary

        The following actions are applied:
        - stitch sequences together into new single sequence
        - write this sequence to a new file
        - replace self.seqfile with new filename, and delete self.seqnames
        - replace feature and primer files in this object with None, as they
          no longer relate to the input sequence
        """
        # Extract relevant file paths
        if outdir is None:
            try:
                outdir = os.path.join(os.path.split(self.seqfile)[:-1])
            except TypeError:  # Raised if file is in current directory
                outdir = ""
        outstem = os.path.splitext(os.path.split(self.seqfile)[-1])[0]
        if self.needs_stitch:
            outfilename = os.path.join(outdir, "{}_concat.fas".format(outstem))
            seqdata = list(SeqIO.parse(self.seqfile, "fasta"))
            catseq = self.spacer.join([str(s.seq) for s in seqdata])
            newseq = SeqRecord(
                Seq(catseq),
                id="_".join([self.name, "concatenated"]),
                description="%s, concatenated with spacers" % self.name,
            )
            SeqIO.write([newseq], outfilename, "fasta")
            self.seqfile = outfilename
            if hasattr(self, "_seqnames"):
                delattr(self, "_seqnames")
            self.features = None
            self.primers = None

    def replace_ambiguities(self, outdir=None):
        """Replace non-N ambiguity symbols in self.seqfile with N, if needed.

        The following actions are applied:
        - replace ambiguity symbols with Ns
        - write new sequence(s) to file
        - replace self.seqfile with new filename; delete self.seqnames
        - replace feature and primer files with None, as they no longer relate
          to the input sequence
        """
        # Extract relevant file paths
        if outdir is None:
            try:
                outdir = os.path.join(os.path.split(self.seqfile)[:-1])
            except TypeError:  # Raised if file is in current directory
                outdir = ""
        outstem = os.path.splitext(os.path.split(self.seqfile)[-1])[0]
        if self.has_ambiguities:
            outfilename = os.path.join(outdir, "{}_noambig.fas".format(outstem))
            seqdata = list(SeqIO.parse(self.seqfile, "fasta"))
            for s in seqdata:
                s.seq = Seq(re.sub(self.ambiguities, "N", str(s.seq)), s.seq.alphabet)
                s.id = "_".join([s.id, "noambig"])
            SeqIO.write(seqdata, outfilename, "fasta")
            self.seqfile = outfilename
            if hasattr(self, "_seqnames"):
                delattr(self, "_seqnames")
            self.features = None
            self.primers = None

    def create_filtered_genome(self, filteredpath, spacerlen, suffix, flanklen=0):
        """Create a new 'filtered_seqfile' of regions specified in self.features.

        - Load feature data
        - Extract regions from self.seqfile
        - Concatenate regions with the appropriate spacer
        - Write the concatenated sequence out and update self.filtered_seqfile

        filteredpath      - path to write filtered_seqfile
        spacerlen         - length of run of Ns to use as spacer between regions
        suffix            - string to add to filtered sequence ID/description etc.
        flanklen          - length of flanking region to include for features (only
                            meant to be used for intergenic regions)
        """
        # The GFF file in self.features should contain a single record with a
        # number of features.
        # We treat each feature as a region of the genome to which it is valid
        # to design a primer, and extract the sequence from self.seqfile,
        # compiling it in a list of sequences.
        # When we are done, we concatenate those sequences with a spacer.
        seqdata = SeqIO.read(self.seqfile, format="fasta")
        regions = list()
        spacer = "N" * spacerlen  # Can't concatenate Seq objects in Biopython yet
        for idx, feature in enumerate(BedTool(self.features)):
            ftr = SeqFeature(
                FeatureLocation(
                    ExactPosition(max(0, feature.start - flanklen)),
                    ExactPosition(min(len(seqdata), feature.end + flanklen)),
                ),
                type="misc",
                id="IGR_{}".format(idx),
            )
            regions.append(ftr.extract(seqdata).seq)
        filtered_seqdata = SeqRecord(
            Seq(spacer.join([str(region) for region in regions])),
            id="_".join([seqdata.id, suffix]),
            name="_".join([seqdata.name, suffix]),
            description=seqdata.description + ", filtered and concatenated",
        )
        SeqIO.write([filtered_seqdata], filteredpath, "fasta")
        self.filtered_seqfile = filteredpath

    def write_primers(self, outfilename, fmt="fasta"):
        """Write the primers for this object to file.

        outfilename  - path to output file
        fmt -       sequence format to write

        The output file format is controlled by Biopython's formatting
        """
        if self.primers is None:
            raise ValueError("No primer file is defined for this object")

        with open(self.primers, "r") as infh:
            primers = json.load(infh)

        seqrecords = []

        for primer in primers:
            seqrecords.append(
                SeqRecord(
                    Seq(primer["forward_seq"]),
                    id=primer["name"] + "_fwd",
                    description="",
                )
            )
            seqrecords.append(
                SeqRecord(
                    Seq(primer["reverse_seq"]),
                    id=primer["name"] + "_rev",
                    description="",
                )
            )
            if len(primer["internal_seq"]):  # This is '' id no oligo
                seqrecords.append(
                    SeqRecord(
                        Seq(primer["internal_seq"]),
                        id=primer["name"] + "_int",
                        description="",
                    )
                )

        return SeqIO.write(seqrecords, outfilename, fmt)

    @property
    def name(self):
        """Identifier for the GenomeData object."""
        return self._name

    @name.setter
    def name(self, value):
        try:
            self._name = str(value)
        except (TypeError, ValueError):
            raise TypeError("PDPData name should be a string")

    @property
    def groups(self):
        """Groups to which the GenomeData object belongs"""
        return sorted(list(self._groups))

    @groups.setter
    def groups(self, value):
        """Expects a list, a set, or a comma-separated string."""
        if isinstance(value, list):
            self._groups = self._groups.union(set(value))
        elif isinstance(value, set):
            self._groups = self._groups.union(value)
        elif isinstance(value, str):
            self._groups = self._groups.union(value.split(","))
        else:
            raise TypeError(
                "PDPData groups should be set, list or " + "comma-separated str"
            )

    @property
    def seqfile(self):
        """Path to input sequence file."""
        return self._seqfile

    @seqfile.setter
    def seqfile(self, value):
        if not os.path.isfile(value):
            raise OSError("%s is not a valid file path" % value)
        self._seqfile = value
        self._filestem = os.path.splitext(os.path.split(self._seqfile)[-1])[0]

    @property
    def target_amplicons(self):
        """Path to target_amplicons file."""
        return self._target_amplicons

    @target_amplicons.setter
    def target_amplicons(self, value):
        if not os.path.isfile(value):
            raise OSError("%s is not a valid file path" % value)
        self._target_amplicons = value
        self._filestem = os.path.splitext(os.path.split(self._target_amplicons)[-1])[0]

    @property
    def filtered_seqfile(self):
        """Path to filtered input sequence file."""
        return self._filtered_seqfile

    @filtered_seqfile.setter
    def filtered_seqfile(self, value):
        if value is not None:
            if not os.path.isfile(value):
                raise OSError("%s is not a valid file path" % value)
            self._filtered_seqfile = value

    @property
    def filestem(self):
        """Filestem for input sequence file."""
        return self._filestem

    @property
    def features(self):
        """Path to feature location file."""
        return self._features

    @features.setter
    def features(self, value):
        if value is not None:
            if not os.path.isfile(value):
                raise OSError("%s is not a valid file path" % value)
        self._features = value

    @property
    def primers(self):
        """Path to primer file."""
        return self._primers

    @primers.setter
    def primers(self, value):
        if value is not None:
            if not os.path.isfile(value):
                raise OSError("%s is not a valid file path" % value)
        self._primers = value

    @property
    def primersearch(self):
        """Path to primersearch file."""
        return self._primersearch

    @primersearch.setter
    def primersearch(self, value):
        if value is not None:
            if not os.path.isfile(value):
                raise OSError("%s is not a valid file path" % value)
        self._primersearch = value

    @property
    def seqnames(self):
        """Lazily returns list of names of sequences in self.seqfile."""
        if not hasattr(self, "_seqnames"):
            self._seqnames = [s.id for s in SeqIO.parse(self.seqfile, "fasta")]
        return self._seqnames

    @property
    def needs_stitch(self):
        """Returns True if more than one sequence in self.seqfile."""
        return len(self.seqnames) > 1

    @property
    def has_ambiguities(self):
        """Returns True if the sequence(s) have non-N ambiguity symbols."""
        for s in SeqIO.parse(self.seqfile, "fasta"):
            if re.search(self.ambiguities, str(s.seq)):
                return True
