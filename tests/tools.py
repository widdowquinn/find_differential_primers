#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tools.py

This module provides helper functions for tests.

(c) The James Hutton Institute 2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

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

import json
import os

from nose.tools import (
    assert_equal, )

from diagnostic_primers import (
    blast, )


def assert_dirfiles_equal(dir1, dir2, listonly=False, filter=None):
    """Assert that dir1 and dir2 contain the same files and they are the same

    - dir1         path to directory of files for comparison
    - dir2         path to directory of files for comparison
    - listonly     Boolean, if True then don't compare file contents
    - filter       Iterable of file extensions to compare

    The list of files present in dir1 and dir2 are compared, and then
    the contents of each file are compared.

    The test is performed differently depending on file extension:

    - .json        test equality of ordered dictionary
    - .blasttab    test that all the columns are equal
    - .ePrimer3    test that each line after the header is equal
    - otherwise test for equality of the file contents directly
    """
    # Skip hidden files
    dir1files = [_ for _ in os.listdir(dir1) if not _.startswith('.')]
    dir2files = [_ for _ in os.listdir(dir2) if not _.startswith('.')]
    assert_equal(
        sorted(dir1files),
        sorted(dir2files),
        msg="%s and %s have differing contents" % (dir1, dir2))
    if filter is not None:
        dir1files = [_ for _ in dir1files if os.path.splitext(_)[-1] in filter]
    if not listonly:
        for fname in dir1files:
            msg = "%s not equal in both directories (%s, %s)" % (fname, dir1,
                                                                 dir2)
            # TODO: make this a distribution dictionary
            with open(os.path.join(dir1, fname)) as ofh:
                with open(os.path.join(dir2, fname)) as tfh:
                    if os.path.splitext(fname)[-1] == '.json':
                        # Order the parsed JSON before comparing
                        assert_equal(
                            ordered(json.load(ofh)),
                            ordered(json.load(tfh)),
                            msg=msg)
                    elif os.path.splitext(fname)[-1].lower() == '.eprimer3':
                        # Skip first line
                        assert_equal(
                            ofh.readlines()[1:], tfh.readlines()[1:], msg=msg)
                    elif os.path.splitext(fname)[-1].lower() == '.blasttab':
                        # Depending on BLAST version, columns may be formatted
                        # differently, so we can't compare characters only
                        data1 = blast.parse_blasttab(ofh)
                        data2 = blast.parse_blasttab(tfh)
                        for line1, line2 in zip(data1, data2):
                            assert_equal(line1, line2)
                    else:  #  Compare the file contents directly
                        assert_equal(ofh.read(), tfh.read(), msg=msg)


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
