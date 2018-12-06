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
    scripts=["pdp.py", os.path.join("bin", "pdp_mafft_wrapper.py")],
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
