# try using distribute or setuptools or distutils.
try:
    import distribute_setup
    distribute_setup.use_setuptools()
except:
    pass

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


import sys
import re

# parse version from package/module without importing or evaluating the code
with open('find_differential_primers/find_differential_primers.py') as fh:
    for line in fh:
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group('version')
            break

if sys.version_info <= (2, 6):
    sys.stderr.write("ERROR: pyADHoRe requires Python Version 2.7 or " +
                     "above...exiting.\n")
    sys.exit(1)

setup(
    name = "find_differential_primers.py",
    version = version,
    author = "Leighton Pritchard",
    author_email = "leighton.pritchard@hutton.ac.uk",
    description = "find_differential_primers.py is a script providing a " +\
    "PCR primer prediction pipeline.",
    license = "GPLv3",
    keywords = "bioinformatics PCR qPCR primers script",
    platforms = "Posix; MacOS X",
    url = "https://github.com/widdowquinn/find_differential_primers",
    download_url = \
    "https://github.com/widdowquinn/find_differential_primers/releases",
    scripts = ['find_differential_primers/find_differential_primers.py'],
    packages = [],
    install_requires = ['biopython', 'bx-python'],
    package_data = {},
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPLv3 License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )
