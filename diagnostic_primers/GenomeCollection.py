#!/usr/bin/env python3
#
# GenomeCollection.py
#
# Class to hold information about a single genome input to diagnostic design
# 
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016 The James Hutton Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import csv

from .GenomeData import GenomeData

# Custom exception for parsing config file
class ConfigSyntaxError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


# Class representing a collection of several GenomeData objects for analysis
class GenomeCollection(object):
    """Container for one or more GenomeData objects."""
    def __init__(self, name, config_file=None):
        """Instantiates GenomeCollection.

        name - an identifier for the collection
        """
        self.name = str(name)
        self._data = {}  # Hold GenomeData objects in dictionary, keyed by ID
        if config_file is not None:
            self.read_config(config_file)

    def add_gdata(self, name=None, groups=None, seqfile=None, features=None,
                  primers=None, relpath=None):
        """Create a new GenomeData object from the passed data, and add to
        self._data.

        This method expects sufficient information to instantiate a GenomeData
        object, and ignores everything in *cols:

        name - identifier for GenomeData object
        groups - diagnostic groups to which object belongs
        """
        self._data[name] = GenomeData(name, groups, seqfile, features,
                                      primers)

    def read_config(self, filename):
        """Compiles a GenomeCollection from the passed config file.

        filename - config filename
        relpath - relative path from calling script/code to config file
        """
        with open(filename, newline='') as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row[0].startswith('#'):
                    self.add_gdata(*self.__parse_row(row))

    def write(self, filename):
        """Writes a new config file from the data in the collection."""
        with open(filename, 'w', newline='') as o:
            writer = csv.writer(o, delimiter='\t')
            writer.writerows([v.as_list() for (k, v) in 
                              sorted(self._data.items())])

    def groups(self):
        """Returns sorted list of groups represented in collection."""
        groups = set()
        for d in self._data.values():
            groups = groups.union(d.groups)
        return sorted(list(groups))

                
    def __parse_row(self, row):
        """Parse a config file row. Returns parsed row as list."""
        if not len(row) == 5:
            raise ConfigSyntaxError("Row must contain five columns, " +\
                                    "got %d columns at %s" %
                                    (len(row), "\t".join(row)))
        # Replace '-' placeholders with None
        return [None if e == '-' else e for e in row]                

    def __len__(self):
        """Return number of items in collection."""
        return len(self._data)

    @property
    def data(self):
        """List of contained GenomeData objects."""
        return [v for (k, v) in sorted(self._data.items())]
