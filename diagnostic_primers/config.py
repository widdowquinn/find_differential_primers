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


# Custom exception for parsing config files
class ConfigSyntaxError(Exception):
    def __init__(self, message):
        super(ConfigSyntaxError, self).__init__(message)


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

        Comment lines begin with '#' and are ignored.
        """
        with open(filename, newline='') as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if not row[0].startswith('#'):
                    self.add_data(*self.__parse_row(row))

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

    def __parse_row(self, row):
        """Parse row from a tab-format config file to list.

        For tab-format configs, we expect three or four columns. We throw
        an error for anything else.
        """
        if len(row) < 3 or len(row) > 4:
            raise ConfigSyntaxError("Row must contain 3 or 4 columns, " +
                                    "got %d at %s" %
                                    (len(row), "\t".join(row)))
        # If we need to process/validate/modify row data, do it here
        return list(row)

    @property
    def data(self):
        """List of contained PDPData objects."""
        return [v for (k, v) in sorted(self._data.items())]
