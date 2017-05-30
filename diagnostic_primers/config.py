#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""config.py

Provides PDPCollection and PDPData classes

(c) The James Hutton Institute 2017

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

Copyright (c) 2017 The James Hutton Institute
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
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class ConfigSyntaxError(Exception):
    """Custom exception for parsing config files."""
    def __init__(self, message):
        super(ConfigSyntaxError, self).__init__(message)


class PDPEncoder(json.JSONEncoder):

    """JSON encoder for PDP objects."""

    def default(self, obj):
        if not isinstance(obj, PDPData):
            return super(PDPEncoder, self).default(obj)

        # Convert complex PDPData object to serialisable dictionary and return
        objdict = {'name': obj.name,
                   'groups': obj.groups,
                   'seqfile': obj.seqfile,
                   'features': obj.features,
                   'primers': obj.primers}

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
        with open(filename, newline='') as ifh:
            reader = csv.reader(ifh, delimiter='\t')
            for row in reader:
                if not row[0].startswith('#'):
                    self.add_data(*self.__parse_row(row))

    def from_json(self, filename):
        """Load data from JSON format config file.

        JSON format config files describe arrays of input genomes/sequences,
        where each array entry is a JSON object with the following keys:

        'name', 'groups', 'seqfile', 'features', 'primers'

        These are used directly to populate the collection's PDPData objects.
        """
        with open(filename, 'r') as ifh:
            data = json.load(ifh)
        for input in data:
            self.add_data(input['name'], input['groups'], input['seqfile'],
                          input['features'], input['primers'])

    def add_data(self, name=None, groups=None, seqfile=None, features=None,
                 primers=None):
        """Create a new PDPData object from passed info and add to collection.

        name     -    unique identifier for object
        groups   -    list of groups to which the object belongs
        seqfile  -    path to sequence file
        features -    path to regions for inclusion/exclusion
        primers  -    path to primers in JSON format
        """
        self._data[name] = PDPData(name, groups, seqfile, features, primers)

    def write_json(self, outfilename):
        """Write the Collection data contents to JSON format config file.

        outfilename  -    path to JSON config file

        Writes an array of serialised PDPData objects to JSON configuration
        file.
        """
        with open(outfilename, 'w') as ofh:
            json.dump(self.data, ofh, cls=PDPEncoder)

    def __parse_row(self, row):
        """Parse row from a tab-format config file to list.

        For tab-format configs, we expect three or four columns. We throw
        an error for anything else.
        """
        if len(row) < 3 or len(row) > 4:
            raise ConfigSyntaxError("Row must contain 3 or 4 columns, " +
                                    "got %d at %s" %
                                    (len(row), "\t".join(row)))

        # Column 2 should be a comma-separated list of strings, defining
        # groups to which an input belongs:
        groups = [e.strip() for e in row[1].split(',')]
        row[1] = groups

        # Substitute None for '-' placeholders
        return [None if e == '-' else e for e in row]

    def __len__(self):
        """Return number of items in collection."""
        return len(self._data)

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


# Class defining paths to data for inputs and methods to operate on it
class PDPData(object):

    """Container for input sequence data and operations on that data."""

    def __init__(self, name, groups, seqfile, features, primers):
        self._name = ""         # Set up private attributes
        self._groups = set()
        self._seqfile = None
        self._features = None
        self._primers = None
        self.cmds = {}           # command-lines used to generate this object
        self.name = name        # Populate attributes
        self.groups = groups
        self.seqfile = seqfile
        self.features = features
        self.primers = primers
        # Useful values
        self.spacer = "NNNNNCATCCATTCATTAATTAATTAATGAATGAATGNNNNN"
        self.ambiguities = re.compile('[BDHKMRSVWY]')

    def stitch(self):
        """Stitch sequences in the sequence file, if necessary

        The following actions are applied:
        - stitch sequences together into new single sequence
        - write this sequence to a new file
        - replace self.seqfile with new filename, and delete self.seqnames
        - replace feature and primer files in this object with None, as they
          no longer relate to the input sequence
        """
        if self.needs_stitch:
            seqdata = list(SeqIO.parse(self.seqfile, 'fasta'))
            catseq = self.spacer.join([str(s.seq) for s in seqdata])
            newseq = SeqRecord(Seq(catseq),
                               id='_'.join([self.name, 'concatenated']),
                               description="%s, concatenated with spacers" %
                               self.name)
            outfilename = ''.join([os.path.splitext(self.seqfile)[0],
                                   '_concat', '.fas'])
            SeqIO.write([newseq], outfilename, 'fasta')
            self.seqfile = outfilename
            if hasattr(self, '_seqnames'):
                delattr(self, '_seqnames')
            self.features = None
            self.primers = None

    def replace_ambiguities(self):
        """Replace non-N ambiguity symbols in self.seqfile with N, if needed.

        The following actions are applied:
        - replace ambiguity symbols with Ns
        - write new sequence(s) to file
        - replace self.seqfile with new filename; delete self.seqnames
        - replace feature and primer files with None, as they no longer relate
          to the input sequence
        """
        if self.has_ambiguities:
            seqdata = list(SeqIO.parse(self.seqfile, 'fasta'))
            for s in seqdata:
                s.seq = Seq(re.sub(self.ambiguities, 'N', str(s.seq)),
                            s.seq.alphabet)
                s.id = '_'.join([s.id, 'noambig'])
            outfilename = ''.join([os.path.splitext(self.seqfile)[0],
                                   '_noambig', '.fas'])
            SeqIO.write(seqdata, outfilename, 'fasta')
            self.seqfile = outfilename
            if hasattr(self, '_seqnames'):
                delattr(self, '_seqnames')
            self.features = None
            self.primers = None

    @property
    def name(self):
        """Identifier for the GenomeData object."""
        return self._name

    @name.setter
    def name(self, value):
        try:
            self._name = str(value)
        except:
            raise TypeError("PDPData name should be a string")

    @property
    def groups(self):
        """Groups to which the GenomeData object belongs"""
        return list(self._groups)

    @groups.setter
    def groups(self, value):
        """Expects a list, a set, or a comma-separated string."""
        if isinstance(value, list):
            self._groups = self._groups.union(set(value))
        elif isinstance(value, set):
            self._groups = self._groups.union(value)
        elif isinstance(value, str):
            self._groups = self._groups.union(value.split(','))
        else:
            raise TypeError("PDPData groups should be set, list or " +
                            "comma-separated str")

    @property
    def seqfile(self):
        """Path to input sequence file."""
        return self._seqfile

    @seqfile.setter
    def seqfile(self, value):
        if not os.path.isfile(value):
            raise OSError("%s is not a valid file path" % value)
        self._seqfile = value

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
    def seqnames(self):
        """Lazily returns list of names of sequences in self.seqfile."""
        if not hasattr(self, "_seqnames"):
            self._seqnames = [s.id for s in SeqIO.parse(self.seqfile, 'fasta')]
        return self._seqnames

    @property
    def needs_stitch(self):
        """Returns True if more than one sequence in self.seqfile."""
        return len(self.seqnames) > 1

    @property
    def has_ambiguities(self):
        """Returns True if the sequence(s) have non-N ambiguity symbols."""
        for s in SeqIO.parse(self.seqfile, 'fasta'):
            if re.search(self.ambiguities, str(s.seq)):
                return True
