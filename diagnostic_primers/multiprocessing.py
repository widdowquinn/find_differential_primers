#!/usr/bin/env python3
#
# prodigal.py
#
# Code to conduct CDS prediction with Prodigal
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

"""Code to run a set of command-line jobs using multiprocessing.

For parallelisation on multi-core desktop/laptop systems, etc. we use
Python's multiprocessing module to distribute command-line jobs.
"""

import multiprocessing
import subprocess
import sys

# Run a set of command lines using multiprocessing
def run(cmdlines, workers=None, verbose=False):
    """Distributes passed command-line jobs using multiprocessing.

    - cmdlines - an iterable of command line strings

    Returns CUMRETVAL, a sum of exit codes from each job that was run. If
    all goes well, we should have CUMRETVAL==0. Anything else and the calling
    function should act accordingly.
    """
    # Keep track of return values for this pool, reset to zero
    global CUMRETVAL
    # Run jobs
    # If workers is None or greater than the number of cores available,
    # it will be set to the maximum number of cores
    pool = multiprocessing.Pool(processes=workers)
    results = [pool.apply_async(subprocess.call, (str(cline), ),
                                {'shell': sys.platform != "win32"})
               for cline in cmdlines]
    pool.close()        # Run jobs
    pool.join()         # Collect output
    return [r.get() for r in results]
