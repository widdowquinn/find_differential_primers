#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""multiprocessing.py

Code to run a set of command-line jobs using multiprocessing.

For parallelisation on multi-core desktop/laptop systems, etc. we use
Python's multiprocessing module to distribute command-line jobs.

(c) The James Hutton Institute 2016-2019
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

import multiprocessing
import subprocess
import sys

CUMRETVAL = 0


# Run a job dependency graph with multiprocessing
def run_dependency_graph(jobgraph, workers=None, logger=None):
    """Create and run pools of jobs based on the passed jobgraph.

    :param jobgraph:  list of jobs, which may have dependencies.
    :param verbose:  flag for multiprocessing verbosity
    :param logger:  Logging.Logger (optional)

    The strategy here is to loop over each job in the list of jobs (jobgraph),
    and create/populate a series of Sets of commands, to be run in
    reverse order with multiprocessing_run as asynchronous pools.
    """
    cmdsets = []
    for job in jobgraph:
        cmdsets = populate_cmdsets(job, cmdsets, depth=1)

    # Put command sets in reverse order, and submit to multiprocessing_run
    cmdsets.reverse()
    cumretval = 0
    for cmdset in cmdsets:
        if logger:  # Try to be informative, if the logger module is being used
            logger.info("Command pool now running:")
            for cmd in cmdset:
                logger.info(cmd)
        cumretval += run(cmdset, workers)
        if logger:  # Try to be informative, if the logger module is being used
            logger.info("Command pool done.")
    return cumretval


def populate_cmdsets(job, cmdsets, depth):
    """Create list of jobsets at different depths of dependency tree.

    :param job:  Job object
    :param cmdsets:  sets of command-lines, one at each depth
    :param depth:  depth of dependencies to begin at

    This is a recursive function (is there something quicker in the itertools
    module?) that descends each 'root' job in turn, populating each
    """
    if len(cmdsets) < depth:
        cmdsets.append(set())
    if isinstance(job.command, list):  # Typical command representation
        cmdsets[depth - 1].add(" ".join(job.command))
    else:  # Command object
        cmdsets[depth - 1].add(str(job.command))
    if len(job.dependencies) == 0:
        return cmdsets
    for j in job.dependencies:
        cmdsets = populate_cmdsets(j, cmdsets, depth + 1)
    return cmdsets


# Run a set of command lines using multiprocessing
def run(cmdlines, workers=None):
    """Distributes passed command-line jobs using multiprocessing.

    :param cmdlines:  an iterable of command line strings
    :param workers:  number of CPUS to use

    Returns CompletedProcess objects for each command. These provide access
    to the return code, stdout and stderr, along with the arguments that
    launched the process (in this case the full command-line).
    """
    # Run jobs
    # If workers is None or greater than the number of cores available,
    # it will be set to the maximum number of cores
    # stderr is piped to DEVNULL because piping to PIPE could block other
    # threads (seen on the JHI development machine as a result of an old
    # Prodigal version). We may want to revisit this to capture the output
    # of processes in a Manager.
    # The command-lines in this package may be provided as any object whose
    # __str__() attribute returns the command-line as a string.
    pool = multiprocessing.Pool(processes=workers)
    results = [
        pool.apply_async(
            subprocess.run,
            (str(cline),),
            {
                "shell": sys.platform != "win32",
                "stdout": subprocess.PIPE,
                "stderr": subprocess.PIPE,
            },
        )
        for cline in cmdlines
    ]
    pool.close()  # Run jobs
    pool.join()  # Collect output
    return sum([r.get().returncode for r in results])
