#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""setup.py

Setup configuration file for use with setuptools

(c) The James Hutton Institute 2018-2019
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

Copyright (c) 2018-2019 The James Hutton Institute

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

try:
    import distribute_setup

    distribute_setup.use_setuptools()
except ImportError:
    pass

import setuptools

import os
import sys
import re

# Get long description from README.md
with open("README.md", "r") as dfh:
    long_description = dfh.read()

# parse version from package/module without importing or evaluating the code
with open("diagnostic_primers/__init__.py") as fh:
    for line in fh:
        m = re.search(r'^__version__ = "(?P<version>[^"]+)"$', line)
        if m:
            version = m.group("version")
            break

if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: diagnostic_primers requires Python3 (exiting)\n")
    sys.exit(1)

setuptools.setup(
    name="diagnostic_primers",
    version=version,
    author="Leighton Pritchard",
    author_email="leighton.pritchard@hutton.ac.uk",
    description=(
        "diagnostic_primers is a module providing tools for "
        "diagnostic qPCR and metabarcoding primer design.",
    ),
    long_description=long_description,
    license="MIT",
    keywords="bioinformatics PCR qPCR primers script",
    platforms="Posix; MacOS X",
    url="https://github.com/widdowquinn/find_differential_primers",
    download_url="https://github.com/widdowquinn/find_differential_primers/releases",
    entry_points={
        "console_scripts": ["pdp = diagnostic_primers.scripts.pdp_script:run_pdp_main"]
    },
    scripts=[
        "pdp.py",
        os.path.join("bin", "pdp_mafft_wrapper.py"),
        os.path.join("bin", "delta_filter_wrapper.py"),
    ],
    packages=[
        "diagnostic_primers",
        "diagnostic_primers/scripts",
        "diagnostic_primers/scripts/subcommands",
    ],
    install_requires=["biopython", "pandas", "plotly", "joblib", "tqdm", "pybedtools"],
    package_data={},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
