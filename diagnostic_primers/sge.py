#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sge_jobs.py

Code to run a set of command-line jobs using SGE/Grid Engine

For parallelisation on multi-node system, we use some custom code to submit
jobs.

(c) The James Hutton Institute 2013-2019
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

Copyright (c) 2013-2018 The James Hutton Institute

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

import itertools
import os
import shlex
import subprocess

from collections import defaultdict

from diagnostic_primers.sge_jobs import JobGroup

QSUB_DEFAULT = "qsub"

JGPREFIX = "pdp"


def split_seq(iterable, size):
    """Splits a passed iterable into chunks of a given size.

    :param iterable:  iterable object to split
    :param size:  the size of chunks to return
    """
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))


# Convert joblist into jobgroups
def compile_jobgroups_from_joblist(joblist, jgprefix, sgegroupsize):
    """Return list of jobgroups, rather than list of jobs.

    :param joblist:  list of Job objects to be run
    :param jgprefix:  prefix string for this set of jobs
    :param sgegroupsize:  the number of individual Jobs to group together
    """
    jobcmds = defaultdict(list)
    for job in joblist:
        if not isinstance(job.command, list):
            cline = str(job.command).split()
            jobcmds[cline[0]].append(" ".join(cline))
        else:
            jobcmds[job.command[0]].append(" ".join(job.command))
    jobgroups = []
    for cmds in list(jobcmds.items()):
        # Break arglist up into batches of sgegroupsize (default: 10,000)
        sublists = split_seq(cmds[1], sgegroupsize)
        count = 0
        for sublist in sublists:
            count += 1
            sge_jobcmdlist = ['"%s"' % jc for jc in sublist]
            jobgroups.append(
                JobGroup(
                    "%s_%d" % (jgprefix, count),
                    "$cmds",
                    arguments={"cmds": sge_jobcmdlist},
                )
            )
    return jobgroups


# Run a job dependency graph, with SGE
def run_dependency_graph(
    jobgraph, logger=None, jgprefix="PDP_SGE_JG", sgegroupsize=10000, sgeargs=None
):
    """Create and runs SGE scripts for jobs based on passed jobgraph.

    :param jobgraph:  JobGraph of Job objects
    :param verbose:  flag for verbosity
    :param logger:  Logging.Logger object
    :param jgprefix:  prefix string for submitted jobs
    :param sgegroupsize:  the maximum size for an array job submission
    :param sgeargs:  additional arguments to qsub

    The strategy here is to loop over each job in the dependency graph
    and, because we expect a single main delta-filter (wrapped) job,
    with a single nucmer dependency for each analysis, we can split
    the dependency graph into two lists of corresponding jobs, and
    run the corresponding nucmer jobs before the delta-filter jobs.
    """
    jobs_main = []  # Can be run first, before deps
    jobs_deps = []  # Depend on the main jobs

    # Try to be informative by telling the user what jobs will run
    dep_count = 0  # how many dependencies are there
    if logger:
        logger.info("Jobs to run with scheduler")
        for job in jobgraph:
            logger.info("{0}: {1}".format(job.name, job.command))
            jobs_main.append(job)
            if len(job.dependencies):
                dep_count += len(job.dependencies)
                for dep in job.dependencies:
                    logger.info("\t[^ depends on: %s (%s)]", dep.name, dep.command)
                    jobs_deps.append(dep)
    logger.info("There are %d job dependencies" % dep_count)
    # Clear dependencies in main group
    for job in jobs_main:
        job.dependencies = []

    # We can use an array (or series of arrays) to schedule our jobs.
    # This cuts down on social problems with long job lists choking up
    # the queue.
    # We split the main and dependent jobs into separate JobGroups.
    # These JobGroups are paired, in order
    logger.info("Compiling main and dependent jobs into separate JobGroups")
    maingroups = compile_jobgroups_from_joblist(
        jobs_main, jgprefix + "_main", sgegroupsize
    )
    depgroups = compile_jobgroups_from_joblist(
        jobs_deps, jgprefix + "_deps", sgegroupsize
    )

    # Assign dependencies to jobgroups
    for mg, dg in zip(maingroups, depgroups):
        mg.add_dependency(dg)
    jobgroups = maingroups + depgroups

    # Send jobs to scheduler
    logger.info("Running jobs with scheduler...")
    logger.info("Jobs passed to scheduler in order:")
    for job in jobgroups:
        logger.info("\t%s" % job.name)
    build_and_submit_jobs(os.curdir, jobgroups, sgeargs)
    logger.info("Waiting for SGE-submitted jobs to finish (polling)")
    for job in jobgroups:
        job.wait()


def populate_jobset(job, jobset, depth):
    """ Return a set of jobs, flattening the dependency tree

    :param job:  Job to add to the set of jobs
    :param jobset:  set of Jobs
    :param depth:  depth of the dependency tree to add jobs to

    The returned set of Jobs contains Jobs at different depths of the
    dependency tree, retaining dependencies as strings, not Jobs.
    """
    jobset.add(job)
    if len(job.dependencies) == 0:
        return jobset
    for j in job.dependencies:
        jobset = populate_jobset(j, jobset, depth + 1)
    return jobset


def build_directories(root_dir):
    """Constructs SGE subdirectories at passed location

    :param root_dir:  Path to the top-level directory for creation of subdirectories

    Directories output, stderr, stdout, and jobs are created in the
    passed root directory
    """
    # If the root directory doesn't exist, create it
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    # Create subdirectories
    directories = [
        os.path.join(root_dir, subdir)
        for subdir in ("output", "stderr", "stdout", "jobs")
    ]
    for dirname in directories:
        os.makedirs(dirname, exist_ok=True)


def build_job_scripts(root_dir, jobs):
    """Constructs the script for each passed Job in the jobs iterable

    :param root_dir:  path to output directory
    """
    # Loop over the job list, creating each job script in turn, and then adding
    # scriptPath to the Job object
    for job in jobs:
        scriptPath = os.path.join(root_dir, "jobs", job.name)
        with open(scriptPath, "w") as scriptFile:
            scriptFile.write("#!/bin/sh\n#$ -S /bin/bash\n%s\n" % job.script)
        job.scriptPath = scriptPath


def extract_submittable_jobs(waiting):
    """Returns a list of jobs from pending list that can be submitted

    :param waiting: List of Job objects
    """
    submittable = set()  # Holds jobs that are able to be submitted
    # Loop over each job, and check all the subjobs in that job's dependency
    # list.  If there are any, and all of these have been submitted, then
    # append the job to the list of submittable jobs.
    for job in waiting:
        unsatisfied = sum([(subjob.submitted is False) for subjob in job.dependencies])
        if 0 == unsatisfied:
            submittable.add(job)
    return list(submittable)


def submit_safe_jobs(root_dir, jobs, sgeargs=None):
    """Submit the passed list of jobs to the Grid Engine server, using the passed
    directory as the root for scheduler output.

    - root_dir      Path to output directory
    - jobs          Iterable of Job objects
    """
    # Loop over each job, constructing SGE command-line
    for job in jobs:
        job.out = os.path.join(root_dir, "stdout")
        job.err = os.path.join(root_dir, "stderr")

        # Add job name, current working directory, SGE stdout and stderr
        # directories to the SGE command line
        args = " -N %s " % (job.name)
        args += " -cwd "
        args += " -o %s -e %s " % (job.out, job.err)

        # If a queue is specified, add this to the SGE command line
        # if job.queue is not None and job.queue in local_queues:
        #    args += local_queues[job.queue]
        #    #args += "-q %s " % job.queue

        # If the job is actually a JobGroup, add the task numbering argument
        if isinstance(job, JobGroup):
            args += "-t 1:%d " % (job.tasks)

        # If there are dependencies for this job, hold the job until they are
        # complete
        if len(job.dependencies) > 0:
            args += "-hold_jid "
            for dep in job.dependencies:
                args += dep.name + ","
            args = args[:-1]

        # Build the qsub SGE commandline (passing local environment)
        qsubcmd = "%s -V %s %s" % (QSUB_DEFAULT, args, job.scriptPath)
        if sgeargs is not None:
            qsubcmd = "%s %s" % (qsubcmd, sgeargs)
        args = [shlex.quote(_) for _ in qsubcmd.split()]
        subprocess.run(args)  # nosec
        job.submitted = True  # Set the job's submitted flag to True


def submit_jobs(root_dir, jobs, sgeargs=None):
    """ Submit each of the passed jobs to the SGE server, using the passed
    directory as root for SGE output.

    - root_dir       Path to output directory
    - jobs           List of Job objects
    """
    waiting = list(jobs)  # List of jobs still to be done
    # Loop over the list of pending jobs, while there still are any
    while len(waiting) > 0:
        # extract submittable jobs
        submittable = extract_submittable_jobs(waiting)
        # run those jobs
        submit_safe_jobs(root_dir, submittable, sgeargs)
        # remove those from the waiting list
        for job in submittable:
            waiting.remove(job)


def build_and_submit_jobs(root_dir, jobs, sgeargs=None):
    """Submits the passed iterable of Job objects to SGE, placing SGE's output in
    the passed root directory

    - root_dir   Root directory for SGE and job output

    - jobs       List of Job objects, describing each job to be submitted
    """
    # If the passed set of jobs is not a list, turn it into one.
    # This makes use of a single JobGroup a little more intutitive
    if not isinstance(jobs, list):
        jobs = [jobs]

    # Build and submit the passed jobs
    build_directories(root_dir)  # build all necessary directories
    build_job_scripts(root_dir, jobs)  # build job scripts
    submit_jobs(root_dir, jobs, sgeargs)  # submit the jobs to SGE
