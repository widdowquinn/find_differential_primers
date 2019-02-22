#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tools.py

This module provides helper functions for tests.

(c) The James Hutton Institute 2017-2019
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

import copy
import json
import os
import re
import subprocess
import unittest

from diagnostic_primers import blast


class PDPFileEqualityTests(object):
    """Tests for equality of filetypes used in PDP.

    Each test defines a comparison for a specific filetype, with contents that are
    expected to be equal in some way.
    """

    def assertJsonEqual(self, json1, json2):
        """Assert that two passed JSON files are equal.

        As we can't always be sure that JSON elements are in the same order in otherwise
        equal files, we compare ordered components.
        """
        with open(json1, "r") as fh1:
            with open(json2, "r") as fh2:
                self.assertEqual(
                    ordered(json.load(fh1)),
                    ordered(json.load(fh2)),
                    msg="JSON files {} and {} do not contain equal contents".format(
                        json1, json2
                    ),
                )

    def assertEprimer3Equal(self, fname1, fname2):
        """Assert that two passed ePrimer3 output files are equal.

        This is a standard file comparison, skipping the first line.
        """
        with open(fname1, "r") as fh1:
            with open(fname2, "r") as fh2:
                fdata1 = fh1.readlines()[1:]
                fdata2 = fh2.readlines()[1:]
                self.assertEqual(
                    fdata1,
                    fdata2,
                    msg="ePrimer3 files {} and {} are not equivalent".format(
                        fname1, fname2
                    ),
                )

    def assertBlasttabEqual(self, fname1, fname2):
        """Assert that two passed BLAST+ .tab output files contain the same data.

        This is not a simple comparison, as we can't rely on the same ordering,
        so we parse the files and compare objects.
        """
        with open(fname1, "r") as fh1:
            with open(fname2, "r") as fh2:
                data1 = blast.parse_blasttab(fh1)
                data2 = blast.parse_blasttab(fh2)
                for line1, line2 in zip(data1, data2):
                    self.assertEqual(line1, line2)

    def assertFilesEqual(self, fname1, fname2):
        """Assert that the two passed files have the same contents."""
        with open(fname1, "r") as fh1:
            with open(fname2, "r") as fh2:
                self.assertEqual(
                    fh1.read(),
                    fh2.read(),
                    msg="Files {} and {} do not have the same contents".format(
                        fname1, fname2
                    ),
                )


class PDPTestCase(unittest.TestCase, PDPFileEqualityTests):
    """Specific PDP unit tests."""

    def assertDirsEqual(self, dir1, dir2, filt=None):
        """Assert that two passed directories have the same contents.

        Files to be compared can be restricted using the filter argument. For
        instance:

        assertDirsEqual(d1, d2, filter=".tab") will only compare files with
        the tab extension.

        Directories are compared recursively.
        """
        # List directories and skip hidden files
        dir1files = ordered([_ for _ in os.listdir(dir1) if not _.startswith(".")])
        dir2files = ordered([_ for _ in os.listdir(dir2) if not _.startswith(".")])
        self.assertEqual(
            dir1files,
            dir2files,
            msg="{} and {} do not have same file listings".format(dir1, dir2),
        )
        #  Compare contents of directories; descend through directories, but
        # filter file extensions if needed
        if filt is not None:
            dir1files = [
                _
                for _ in dir1files
                if (os.path.isdir(_) is False) and (os.path.splitext(_)[-1] == filt)
            ]
        for fpath in dir1files:
            if os.path.isdir(fpath):  # Compare dictionaries
                self.assertDirsEqual(
                    os.path.join(dir1, fpath), os.path.join(dir2, fpath)
                )
            else:  # Compare files
                ext = os.path.splitext(fpath)[-1]
                fname1 = os.path.join(dir1, fpath)
                fname2 = os.path.join(dir2, fpath)
                if ext.lower() == ".json":  # Compare JSON files
                    self.assertJsonEqual(fname1, fname2)
                elif ext.lower() == ".blasttab":  # Compare BLAST+ .tab output
                    self.assertBlasttabEqual(fname1, fname2)
                elif ext.lower() == ".eprimer3":  # Compare ePrimer3 output
                    self.assertEprimer3Equal(fname1, fname2)
                else:  # Compare standard files
                    self.assertFilesEqual(fname1, fname2)


# Define primer3 version as global variable
def get_primer3_version():
    """Return primer3_core version as a tuple of ints.

    Assumes primer3_core is in the $PATH
    """
    output = subprocess.run(
        ["primer3_core", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    version_line = [
        _ for _ in output.stderr.split(b"\n") if _.startswith(b"This is primer3")
    ][0].decode("utf-8")
    return tuple(
        [
            int(_)
            for _ in re.search("(?<=release ).*(?=\\))", version_line)
            .group()
            .split(".")
        ]
    )


def ordered(obj):
    """Return ordered version of the passed object

    Dictionaries are not ordered in all Python versions, and the
    implementation of sort_keys() in the the JSON library seems
    erratic in terms of effect
    """
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def modify_namespace(namespace, args):
    """Modify the specified arguments in the passed Namespace.

    namespace       argparse.Namespace object
    args            dict of argument: value pairs

    For most command-line tests, we define a base argparse.Namespace object, then
    change a few arguments. This function takes the base namespace and a dictionary
    of argument: value pairs, and returns the modified namespace.
    """
    new_namespace = copy.deepcopy(namespace)
    for argname, argval in args.items():
        setattr(new_namespace, argname, argval)
    return new_namespace
