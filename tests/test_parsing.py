#!/usr/bin/env python3
#
# test_parsing.py
#
# Tests for config file parsing
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

import os

from diagnostic_primers.GenomeCollection import GenomeCollection

from nose.tools import assert_equal


TESTDATA = "test_data"


def test_parse_config():
    """Test basic parsing of config file."""
    inconf = os.path.join(TESTDATA, "testin.conf")     # input config file
    gc = GenomeCollection("test", config_file=inconf)
    assert_equal(16, len(gc))
    assert_equal(['Pectobacterium', 'atrosepticum_NCBI', 
                  'betavasculorum_NCBI', 'gv1', 'gv2', 'gv3', 'gv7',
                  'wasabiae_NCBI'],
                 gc.groups())

def test_parse_and_write():
    """Test basic parsing/generation of config file."""
    inconf = os.path.join(TESTDATA, "testin.conf")     # input config file
    testconf = os.path.join(TESTDATA, "testout.conf")  # output config file
    newconf = os.path.join(TESTDATA, "new.conf")       # generated config file
    gc = GenomeCollection("test", config_file=inconf)
    gc.write(newconf)
    assert_equal(open(newconf, 'r').read(), open(testconf, 'r').read())




# Run as script
if __name__ == '__main__':
    import inspect
    import test_parsing
    functions = [o[0] for o in inspect.getmembers(test_parsing) if
                 inspect.isfunction(o[1])]
    for fn in functions:
        print("\nFunction called: {}()".format(fn))
        locals()[fn]()
