#!/usr/bin/env python
# find_differential_primers.py
#
# A Python script that identifies pairs of forward and reverse primers which
# are capable of amplifying either individual organisms, or a particular
# family of organisms, from a set of genome sequences.  Primers are expected
# to be located within CDS features, in an attempt to maximise sequence
# stability of the primers.
#
# The script reads from a configuration file containing sequence names and,
# at a minimum, the location of a complete genome sequence.  Optionally, the
# configuration file may also indicate:
# -  the location of a GenBank file containing CDS feature locations,
#    or an equivalent output file from the Prodigal genefinder
#    (http://compbio.ornl.gov/prodigal/)
# -  the locations on the genome, and sequences of, primers predicted in
#    EMBOSS ePrimer3 output format
#    (http://emboss.bioinformatics.nl/cgi-bin/emboss/help/eprimer3)
#
# The first step of the script, if no primer file is specified, is to use
# the sequence file as the basis for a call to EMBOSS ePrimer3
# (http://emboss.bioinformatics.nl/cgi-bin/emboss/help/eprimer3), which must
# be installed and either on the $PATH, or its location specified at the
# command line.  This will generate an output file with the same stem as the
# sequence file, but with the extension '.eprimer3'.  Some ePrimer3 settings,
# such as the number of primers to find, are command-line options.
#
# If no CDS feature file is specified, and the --noCDS flag is not set,
# the script will attempt first to use Prodigal
# (http://compbio.ornl.gov/prodigal/) to predict CDS locations, placing the
# output in the same directory as the sequence source.  If Prodigal cannot be
# found, a warning will be given, and the script will proceed as if the
# --noCDS flag is set.  If this flag is set, then all primers are carried
# through to a query with the EMBOSS PrimerSearch package
# (http://emboss.bioinformatics.nl/cgi-bin/emboss/help/primersearch) against
# all other sequences in the dataset.  If the flag is not set, then all
# primers that are not located within a CDS feature are excluded from the
# PrimerSearch input.  To enable this, the PrimerSearch input is written to
# an intermediate file with the same stem as the input sequence, but the
# extension '.primers'.
#
# A run of PrimerSearch is carried out with every set of primers against
# all other sequences in the dataset.  The output of this search is written to
# a file with the following naming convention:
# <query>_primers_vs_<target>.primersearch
# Where <query> is the name given to the query sequence in the config file, and
# <target> is the name given to the target sequence in the config file.  This
# step is not carried out if the --noprimersearch flag is set.  When this flag
# is set, the script will look for the corresponding PrimerSearch output in
# the same directory as the sequence file, and will report an error if it is
# not present.
#
# Finally, the script uses the PrimerSearch results to identify primers that
# are unique to each query sequence, and to each family named in the config
# file.  These are reported in files with the following naming convention:
# <query>_specific_primers.eprimer3
# <family>_specific_primers.primers
# We use ePrimer3 format for the family-specific primers, even though the
# start and end positions are meaningless, as they will amplify at different
# sites in each family member.  However, the source sequence is indicated in a
# comment line, and the primer sequences and T_m/GC% values should be the same,
# regardless.
# Primers that are universal to all sequences in the sample are written in
# ePrimer3 format to the file:
# universal_primers.eprimer3
# This file has the same caveats as the family-specific file above.
#
# Copyright (C) 2011 Leighton Pritchard
# Further modifications (C) 2012,2013 Leighton Pritchard, (C) 2012 Benjamin Leopold, 
# and (C) 2013 Michael Robeson
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
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# TODO:
#
# 1) Make each of the callout routines create a list of command-line calls, so
#    that these can be passed either to SGE or run via multiprocessing
# 2) Create a second parallelising callout function that runs SGE command-lines
#    and make this safe for asynchronous passage to the next stage of callouts
#    (maybe that should be 2a)
# 3) Make the input configuration file (optionally) XML
# 4) allow option file holding all command line options instead/alongside
#    actual cmdline args


###
# IMPORTS

from Bio import SeqIO                           # Parsing biological sequence data
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML                   # BLAST XML parser
from Bio.Emboss.Applications import Primer3Commandline, PrimerSearchCommandline
from Bio.Emboss import Primer3, PrimerSearch    # Primer3, PrimerSearch parsers
from Bio.GenBank import _FeatureConsumer        # For parsing GenBank locations
from Bio.Seq import Seq                         # Represents a sequence
from Bio.SeqRecord import SeqRecord             # Represents an annotated record
from Bio.SeqFeature import SeqFeature           # Represents an annotated record
from bx.intervals.cluster import ClusterTree    # Interval tree building
from collections import defaultdict             # Syntactic sugar
from optparse import OptionParser               # Cmd-line parsing
import multiprocessing
import os
import subprocess
import sys
import time
import re


###
# CLASSES

# Class describing an organism's genome, and associated data.
class GenomeData:
    """ Describes an organism's genome, and has attributes:

        name   - short, unique (not enforced) identification string
        family - string indicating a family membership
        seqfilename    - location of representative genome sequence file
        ftfilename     - location of GBK/Prodigal feature file
        primerfilename - location of ePrimer3 format primers file

        primers     - dictionary collection of Bio.Emboss.Primer3.Primer
                      objects, keyed by primer name

        Exposed methods are:

    """
    def __init__(self, name, family=None, seqfilename=None, ftfilename=None,
                 primerfilename=None, primersearchfilename=None):
        """ Expects at minimum a name to identify the organism.  Optionally
            filenames describing the location of sequence, feature, and
            primer data may be specified, along with a family classification.

            name   - short, unique (not enforced) identification string
            family - string indicating a family membership
            seqfilename    - location of representative genome sequence file
            ftfilename     - location of GBK/Prodigal feature file
            primerfilename - location of ePrimer3 format primers file
            primersearchfilename - location of PrimerSearch format primers file

            Rather hackily, passing '-' to any of the keyword arguments also
            sets them to None; this is to aid in config file parsing, and
            is a wee bit ugly.
        """
        self.name = name                       # Short identifier
        self.family = family if family != '-' else None
        self.seqfilename = seqfilename if seqfilename != '-' else None
        self.ftfilename = ftfilename if ftfilename != '-' else None
        self.primerfilename = primerfilename if primerfilename != '-' \
            else None
        self.primersearchfilename = primersearchfilename if\
            primersearchfilename != '-' else None
        self.primers = {}              # Dict of Primer objects, keyed by name
        self.load_sequence()

    def load_sequence(self):
        """ Load the sequence defined in self.seqfile into memory.  We
            assume it's FASTA format.  This can then be used to calculate
            amplicons when loading primers in.
        """
        if self.seqfilename is not None:
            try:
                self.sequence = SeqIO.read(open(self.seqfilename, 'rU'), 'fasta')
            except ValueError, data:
                print "Loading sequence file %s failed" % self.seqfilename
                print "Error message: ", data
                sys.exit(1)

    def write_primers(self, verbose=False):
        """ Write the primer pairs in self.primers out to file in an
            appropriate format for PrimerSearch.  If the filename is not
            already defined, the filestem of the
            source sequencefile is used for the output file, with the
            extension '.primers'.
            The method returns the number of lines written.
        """
        # Define output filename, if not already defined
        if self.primersearchfilename is None:
            self.primersearchfilename = os.path.splitext(self.seqfilename)[0] + \
                '.primers'
        if verbose:
            t0 = time.time()
            print "Writing primers to file %s ..." % self.primersearchfilename
        # Open handle and write data
        outfh = open(self.primersearchfilename, 'w')
        print >> outfh, "# Primers for %s" % self.name
        print >> outfh, "# Automatically generated by find_differential_primers"
        for primers in self.primers.values():
            print >> outfh, "%s\t%s\t%s" %\
                (primers.name, primers.forward_seq,
                 primers.reverse_seq)
        if not len(self.primers):
            print >> sys.stderr, "WARNING: no primers written to %s!" % \
                self.primersearchfilename
        # Being tidy
        outfh.close()
        print_verbose("... wrote %d primers to %s (%.3fs)" %\
                (len(self.primers), self.primersearchfilename, time.time() - t0))

    def get_unique_primers(self, cds_overlap=False,
                           gc3primevalid=False, oligovalid=False,
                           blastfilter=False, single_product=False,
                           verbose=False):
        """ Returns a list of primers that have the .amplifies_organism
            attribute, but where this is an empty set.
            If cds_overlap is True, then this list is restricted to those
            primers whose .cds_overlap attribute is also True
        """
        return self.get_primers_amplify_count(0, cds_overlap, gc3primevalid,
                                              oligovalid, blastfilter, single_product,
                                              verbose)

    def get_family_unique_primers(self, family_members, cds_overlap=False,
                                  gc3primevalid=False, oligovalid=False,
                                  blastfilter=False, single_product=False,
                                  verbose=False):
        """ Returns a list of primers that have the .amplifies_organism
            attribute, and where the set of organisms passed in family_members
            is the same as that in .amplifies organism, with the addition of
            self.name.
            If cds_overlap is True, then this list is restricted to those
            primers whose .cds_overlap attribute is also True
        """
        primerlist = []
        for p in self.primers.values():
            if family_members == set([self.name]).union(p.amplifies_organism):
                primerlist.append(p)
        print_verbose("[%s] %d family primers" % (self.name,
                                              len(primerlist)))
        if cds_overlap:
            primerlist = [p for p in primerlist if p.cds_overlap]
            print_verbose("[%s] %d primers after CDS filter" % \
                    (self.name, len(primerlist)))
        if gc3primevalid:
            primerlist = [p for p in primerlist if p.gc3primevalid]
            print_verbose("[%s] %d primers after GC 3` filter" % \
                    (self.name, len(primerlist)))
        if oligovalid:
            primerlist = [p for p in primerlist if p.oligovalid]
            print_verbose("[%s] %d primers after oligo filter" % \
                    (self.name, len(primerlist)))
        if blastfilter:
            primerlist = [p for p in primerlist if p.blastpass]
            print_verbose("[%s] %d primers after BLAST filter" % \
                    (self.name, len(primerlist)))
        if single_product:
            primerlist = [p for p in primerlist if p.negative_control_amplimers == 1]
            print_verbose("[%s] %d primers after single_product filter" % \
                    (self.name, len(primerlist)))
        print_verbose("[%s] returning %d primers" % (self.name, len(primerlist)))
        return primerlist

    def get_primers_amplify_count(self, count, cds_overlap=False,
                                  gc3primevalid=False, oligovalid=False,
                                  blastfilter=False, single_product=False,
                                  verbose=False):
        """ Returns a list of primers that have the .amplifies_organism
            attribute and the length of this set is equal to the passed count.
            If cds_overlap is True, then this list is restricted to those
            primers whose .cds_overlap attribute is also True
        """
        primerlist = [p for p in self.primers.values() if \
                          count == len(p.amplifies_organism)]
        print_verbose("[%s] %d family primers that amplify %d orgs" % \
                (self.name, len(primerlist), count))
        if cds_overlap:
            primerlist = [p for p in primerlist if p.cds_overlap]
            print_verbose("[%s] %d primers after CDS filter" % \
                    (self.name, len(primerlist)))
        if gc3primevalid:
            primerlist = [p for p in primerlist if p.gc3primevalid]
            print_verbose("[%s] %d primers after GC 3` filter" % \
                    (self.name, len(primerlist)))
        if oligovalid:
            primerlist = [p for p in primerlist if p.oligovalid]
            print_verbose("[%s] %d primers after oligo filter" % \
                    (self.name, len(primerlist)))
        if blastfilter:
            primerlist = [p for p in primerlist if p.blastpass]
            print_verbose("[%s] %d primers after BLAST filter" % \
                    (self.name, len(primerlist)))
        if single_product:
            primerlist = [p for p in primerlist if p.negative_control_amplimers == 1]
            print_verbose("[%s] %d primers after single_product filter" % \
                    (self.name, len(primerlist)))
        print_verbose("[%s] returning %d primers" % (self.name, len(primerlist)))
        return primerlist

    def __str__(self):
        """ Pretty string description of object contents
        """
        outstr = ['GenomeData object: %s' % self.name]
        outstr.append('Family: %s' % self.family)
        outstr.append('Sequence file: %s' % self.seqfilename)
        outstr.append('Feature file: %s' % self.ftfilename)
        outstr.append('Primers file: %s' % self.primerfilename)
        outstr.append('PrimerSearch file: %s' % self.primersearchfilename)
        outstr.append('Primers: %d' % len(self.primers))
        if len(self.primers):
            outstr.append('Primers overlapping CDS: %d' %\
                              len([p for p in self.primers.values()\
                                       if p.cds_overlap]))
        return os.linesep.join(outstr) + os.linesep


###
# FUNCTIONS

# Parse command-line options
def parse_cmdline(args):
    """ Parse command line, accepting args obtained from sys.argv
    """
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-i", "--infile", dest="filename", action="store",
                      help="location of configuration file",
                      default=None)
    parser.add_option("--nocds", dest="nocds", action="store_true",
                      help="do not restrict primer prediction to CDS",
                      default=False)
    parser.add_option("--noprodigal", dest="noprodigal", action="store_true",
                      help="do not carry out Prodigal prediction step",
                      default=False)
    parser.add_option("--noprimer3", dest="noprimer3", action="store_true",
                      help="do not carry out ePrimer3 prediction step",
                      default=False)
    parser.add_option("--noprimersearch", dest="noprimersearch",
                      action="store_true",
                      help="do not carry out PrimerSearch step",
                      default=False)
    parser.add_option("--noclassify", dest="noclassify",
                      action="store_true",
                      help="do not carry out primer classification step",
                      default=False)
    parser.add_option("--single_product", dest="single_product",
                      action="store",
                      help="location of FASTA sequence file containing sequences from which " +\
                           "a sequence-specific primer must amplify exactly one product.",
                      default=None)
    parser.add_option("--filtergc3prime", dest="filtergc3prime",
                      action="store_true",
                      help="allow no more than two GC at the 3` end of primers",
                      default=False)
    parser.add_option("--prodigal", dest="prodigal_exe", action="store",
                      help="location of Prodigal executable",
                      default="prodigal")
    parser.add_option("--eprimer3", dest="eprimer3_exe", action="store",
                      help="location of EMBOSS eprimer3 executable",
                      default="eprimer3")
    parser.add_option("--numreturn", dest="numreturn", action="store",
                      help="number of primers to find",
                      default=20, type="int")
    parser.add_option("--osize", dest="osize", action="store",
                      help="optimal size for primer oligo",
                      default=20, type="int")
    parser.add_option("--minsize", dest="minsize", action="store",
                      help="minimum size for primer oligo",
                      default=18, type="int")
    parser.add_option("--maxsize", dest="maxsize", action="store",
                      help="maximum size for primer oligo",
                      default=22, type="int")
    parser.add_option("--otm", dest="otm", action="store",
                      help="optimal melting temperature for primer oligo",
                      default=59, type="int")
    parser.add_option("--mintm", dest="mintm", action="store",
                      help="minimum melting temperature for primer oligo",
                      default=58, type="int")
    parser.add_option("--maxtm", dest="maxtm", action="store",
                      help="maximum melting temperature for primer oligo",
                      default=60, type="int")
    parser.add_option("--ogcpercent", dest="ogcpercent", action="store",
                      help="optimal %GC for primer oligo",
                      default=55, type="int")
    parser.add_option("--mingc", dest="mingc", action="store",
                      help="minimum %GC for primer oligo",
                      default=30, type="int")
    parser.add_option("--maxgc", dest="maxgc", action="store",
                      help="maximum %GC for primer oligo",
                      default=80, type="int")
    parser.add_option("--psizeopt", dest="psizeopt", action="store",
                      help="optimal size for amplified region",
                      default=100, type="int")
    parser.add_option("--psizemin", dest="psizemin", action="store",
                      help="minimum size for amplified region",
                      default=50, type="int")
    parser.add_option("--psizemax", dest="psizemax", action="store",
                      help="maximum size for amplified region",
                      default=150, type="int")
    parser.add_option("--maxpolyx", dest="maxpolyx", action="store",
                      help="maximum run of repeated nucleotides in primer",
                      default=3, type="int")
    parser.add_option("--mismatchpercent", dest="mismatchpercent",
                      action="store",
                      help="allowed percentage mismatch in primersearch",
                      default=10, type="int")
    parser.add_option("--hybridprobe", dest="hybridprobe", action="store_true",
                      help="generate internal oligo as a hybridisation probe",
                      default=False)
    parser.add_option("--oligoosize", dest="oligoosize", action="store",
                      help="optimal size for internal oligo",
                      default=20, type="int")
    parser.add_option("--oligominsize", dest="oligominsize", action="store",
                      help="minimum size for internal oligo",
                      default=13, type="int")
    parser.add_option("--oligomaxsize", dest="oligomaxsize", action="store",
                      help="maximum size for internal oligo",
                      default=30, type="int")
    parser.add_option("--oligootm", dest="oligootm", action="store",
                      help="optimal melting temperature for internal oligo",
                      default=69, type="int")
    parser.add_option("--oligomintm", dest="oligomintm", action="store",
                      help="minimum melting temperature for internal oligo",
                      default=68, type="int")
    parser.add_option("--oligomaxtm", dest="oligomaxtm", action="store",
                      help="maximum melting temperature for internal oligo",
                      default=70, type="int")
    parser.add_option("--oligoogcpercent", dest="oligoogcpercent",
                      action="store",
                      help="optimal %GC for internal oligo",
                      default=55, type="int")
    parser.add_option("--oligomingc", dest="oligomingc", action="store",
                      help="minimum %GC for internal oligo",
                      default=30, type="int")
    parser.add_option("--oligomaxgc", dest="oligomaxgc", action="store",
                      help="maximum %GC for internal oligo",
                      default=80, type="int")
    parser.add_option("--oligomaxpolyx", dest="oligomaxpolyx", action="store",
                      help="maximum run of repeated nt in internal oligo",
                      default=3, type="int")
    parser.add_option("--cpus", dest="cpus", action="store",
                      help="number of CPUs to use in multiprocessing",
                      default=multiprocessing.cpu_count(), type="int")
    parser.add_option("--sge", dest="sge", action="store_true",
                      help="use SGE job scheduler",
                      default=False)
    parser.add_option("-o", "--outdir", dest="outdir", action="store",
                      help="directory for output files",
                      default="differential_primer_results")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="report progress to stdout",
                      default=False)
    parser.add_option("--clean", action="store_true", dest="clean",
                      help="clean up old output files before running",
                      default=False)
    parser.add_option("--cleanonly", action="store_true", dest="cleanonly",
                      help="clean up old output files and exit",
                      default=False)
    parser.add_option("--blast_exe", dest="blast_exe", action="store",
                      help="location of BLASTN/BLASTALL executable",
                      default="blastn")
    parser.add_option("--blastdb", dest="blastdb", action="store",
                      help="location of BLAST database",
                      default=None)
    parser.add_option("--useblast", dest="useblast", action="store_true",
                      help="use existing BLAST results",
                      default=False)
    parser.add_option("--debug", action="store_true", dest="debug",
                      help="report extra progress to stdout for debugging",
                      default=False)
    parser.add_option("--keep_logs", action="store_true", dest="keep_logs",
                      help="store log files from each process",
                      default=False)
    parser.add_option("--log_dir", action="store", dest="log_dir",
                      help="store kept log files in this directory",
                      default=None)
    (options, args) = parser.parse_args()
    return (options, args)


# Create a list of GenomeData objects corresponding to config file entries
def create_gd_from_config(filename, verbose):
    """ Parses data from a configuration file into a list of GenomeData
        objects.
        Returns a list of GenomeData objects.

        Each line of the config file describes a single genome.
        The config file format is six tab-separated columns, where columns
        may be separated by multiple tabs.  'Empty' data values are indicated
        by the '-' symbol, and these are converted into None objects in
        parsing.
        Comment lines start with '#', as in Python.
        The five columns are:
        1) Genome name
        2) Genome family
        3) Location of FASTA format sequence data
        4) Location of GENBANK/PRODIGAL format feature data
        5) Location of EPRIMER3 format primer data
        6) Location of PRIMERSEARCH input format primer data

        The data would, of course, be better presented as an XML file, but it
        might be useful to maintain both tab- and XML-formatted approaches to
        facilitate human construction as well as computational.
    """
    if verbose:
        t0 = time.time()
        print "Creating list of genomes from config file %s ..." % filename
    gdlist = []                                   # Hold GenomeData objects
    # Ignore blank lines and comments...
    for line in [l.strip() for l in open(filename, 'rU')\
                     if l.strip() and not l.startswith('#')]:
        # Split data and create new GenomeData object, adding it to the list
        data = [e.strip() for e in line.strip().split('\t') if e.strip()]
        name, family, sfile, ffile, pfile, psfile = tuple(data)
        gdlist.append(GenomeData(name, family, sfile, ffile, pfile, psfile))
        print_verbose("... created GenomeData object for %s ..." % name)
        print_verbose(gdlist[-1])
    print_verbose("... created %d GenomeData objects (%.3fs)" % (len(gdlist),
                                                             time.time() - t0))
    return gdlist


# Check whether each GenomeData object has multiple sequence and, if so,
# concatenate them sensibly, resetting feature and primer file locations to
# None
def check_single_sequence(gdlist, verbose):
    """ Loops over the GenomeData objects in the passed list and, where the
        sequence file contains multiple sequences, concatenates them into
        a single sequence using a spacer that facilitates gene-finding.  As
        this process changes feature and primer locations, the ftfilename and
        primerfilename attributes are reset to None, and these are
        recalculated later on in the script, where necessary.
    """
    if verbose:
        t0 = time.time()
        print "Checking for multiple sequences ..."
    for gd in gdlist:
        # Verify that the sequence file contains a single sequence
        seqdata = [s for s in SeqIO.parse(open(gd.seqfilename, 'rU'), 'fasta')]
        if len(seqdata) != 1:
            print_verbose("... %s describes multiple sequences ..." %\
                gd.seqfilename)
            gd.seqfilename = concatenate_sequences(gd, verbose)  # Concatenate
            print_verbose("... clearing feature and primer file locations ...")
            gd.ftfilename, gd.primerfilename, gd.primersearchfilename = \
                None, None, None
    print_verbose("... checked %d GenomeData objects (%.3fs)" % (len(gdlist),
                                                             time.time() - t0))
#    return gdlist


# Concatenate multiple fragments of a genome to a single file
def concatenate_sequences(gd, verbose=False):
    """ Takes a GenomeData object and concatenates sequences with the spacer
        sequence NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN (this contains
        start and stop codons in all frames, to cap individual sequences).
        We write this data out to a new file

        For filename convention, we just add '_concatenated' to the end
        of the sequence filestem, and use the '.fas' extension.
    """
    # Spacer contains start and stop codons in all six frames
    spacer = 'NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN'
    if verbose:
        t0 = time.time()
        print "Concatenating sequences from %s ..." % gd.seqfilename
    newseq = SeqRecord(Seq(spacer.join([s.seq.data for s in \
                                        SeqIO.parse(open(gd.seqfilename, 'rU'),
                                                    'fasta')])),
                       id=gd.name + "_concatenated",
                       description="%s, concatenated with spacers" %\
                           gd.name)
    outfilename = os.path.splitext(gd.seqfilename)[0] + '_concatenated' +\
        '.fas'
    SeqIO.write([newseq], open(outfilename, 'w'), 'fasta')
    print_verbose("... wrote concatenated data to %s (%.3fs)" % (outfilename,
                                                             time.time() - t0))
    return outfilename


# Check for each GenomeData object in a passed list, the existence of
# the feature file, and create one using Prodigal if it doesn't exist already
def check_ftfilenames(gdlist, prodigal_exe, poolsize, sge, verbose):
    """ Loop over the GenomeData objects in gdlist and, where no feature file
        is specified, add the GenomeData object to the list of
        packets to be processed in parallel by Prodigal using multiprocessing.
    """
    print_verbose("Checking and predicting features for GenomeData files ...")
    # We split the GenomeData objects into those with, and without,
    # defined feature files, but we don't test the validity of the files
    # that were predefined, here.
    gds_with_ft = [gd for gd in gdlist if
                   (gd.ftfilename is not None and \
                        os.path.isfile(gd.ftfilename))]
    gds_no_ft = [gd for gd in gdlist if \
                     (gd.ftfilename is None or \
                          not os.path.isfile(gd.ftfilename))]
    # Predict features for those GenomeData objects with no feature file
    print_verbose("... %d GenomeData objects have no feature file ..." %\
            len(gds_no_ft))
    print_verbose("... running %d Prodigal jobs to predict CDS ..." %\
            len(gds_no_ft))
    # Create a list of command-line tuples, for Prodigal
    # gene prediction applied to each GenomeData object in gds_no_ft.
    clines = []
    for gd in gds_no_ft:
        gd.ftfilename = os.path.splitext(gd.seqfilename)[0] + '.prodigalout'
        seqfilename = os.path.splitext(gd.seqfilename)[0] + '.features'
        cline = "%s -a %s < %s > %s" % (prodigal_exe, seqfilename,
                                        gd.seqfilename, gd.ftfilename)
        clines.append(cline + log_output(gd.name + ".prodigal"))
    print_verbose("... Prodigal jobs to run:")
    print_verbose_list(clines)
    # Depending on the type of parallelisation required, these command-lines
    # are either run locally via multiprocessing, or passed out to SGE
    if not sge:
        multiprocessing_run(clines, poolsize, verbose)
    else:
        sge_run(clines, verbose)


# Check whether GenomeData objects have a valid primer definition file
def check_primers(gdlist, verbose):
    """ Loop over GenomeData objects in the passed gdlist and, if they have
        a defined primerfilename attribute, attempt to parse it.  If this
        is successful, do nothing.  If it fails, set the primerfilename
        attribute to None.
    """
    if verbose:
        t0 = time.time()
        print "Checking ePrimer3 output files ..."
    for gd in [g for g in gdlist if g.primerfilename]:
        try:
            Primer3.read(open(gd.primerfilename, 'rU'))
            print_verbose("... %s primer file %s OK ..." % (gd.name,
                                                        gd.primerfilename))
        except:
            print_verbose("... %s primer file %s not OK ..." % (gd.name,
                                                            gd.primerfilename))
            gd.primerfilename = None


# Check for each GenomeData object in a passed list, the existence of
# the ePrimer3 file, and create one using ePrimer3 if it doesn't exist already
def predict_primers(gdlist, eprimer3_exe, poolsize, numreturn,
                    osize, minsize, maxsize,
                    otm, mintm, maxtm,
                    ogcpercent, mingc, maxgc,
                    psizeopt, psizemin, psizemax,
                    maxpolyx,
                    hybridprobe, oligoosize, oligominsize, oligomaxsize,
                    oligootm, oligomintm, oligomaxtm, oligoogcpercent,
                    oligomingc, oligomaxgc, oligomaxpolyx,
                    sge, verbose):
    """ Loop over the GenomeData objects in gdlist and, where no primer file
        is specified, add the GenomeData object to the list of
        packets to be processed in parallel by Prodigal using multiprocessing.
    """
    if verbose:
        t0 = time.time()
        print "Checking and predicting primers for GenomeData files ..."
    input_count = len(gdlist)  # For sanity later
    # We need to split the GenomeData objects into those with, and without,
    # defined primer files, but we don't test the validity of these files
    gds_with_primers = [gd for gd in gdlist if gd.primerfilename is not None]
    gds_no_primers = [gd for gd in gdlist if gd.primerfilename is None]
    # Predict primers for those GenomeData objects with no primer file
    print_verbose("... %d GenomeData objects have no primer file ..." %\
            len(gds_no_primers))
    print_verbose("... running %d ePrimer3 jobs to predict CDS ..." %\
            len(gds_no_primers))
    # Create command-lines to run ePrimer3
    clines = []
    for gd in gds_no_primers:
        # Create ePrimer3 command-line.
        cline = Primer3Commandline(cmd=eprimer3_exe)
        cline.sequence = gd.seqfilename
        cline.auto = True
        cline.osize = "%d" % osize            # Optimal primer size
        cline.minsize = "%d" % minsize        # Min primer size
        cline.maxsize = "%d" % maxsize        # Max primer size
        cline.otm = "%d" % otm                # Optimal primer Tm
        cline.mintm = "%d" % mintm            # Min primer Tm
        cline.maxtm = "%d" % maxtm            # Max primer Tm
        cline.ogcpercent = "%d" % ogcpercent  # Optimal primer %GC
        cline.mingc = "%d" % mingc            # Min primer %GC
        cline.maxgc = "%d" % maxgc            # Max primer %GC
        cline.psizeopt = "%d" % psizeopt      # Optimal product size
        cline.maxpolyx = "%d" % maxpolyx      # Longest polyX run in primer
        cline.prange = "%d-%d" % (psizemin, psizemax)  # Allowed product sizes
        cline.numreturn = "%d" % numreturn    # Number of primers to predict
        cline.hybridprobe = hybridprobe       # Predict internal oligo?
        cline.osizeopt = "%d" % oligoosize    # Internal oligo parameters;
        cline.ominsize = "%d" % oligominsize  # We use EMBOSS v6 parameter
        cline.omaxsize = "%d" % oligomaxsize  # names, here.
        cline.otmopt = "%d" % oligootm
        cline.otmmin = "%d" % oligomintm
        cline.otmmax = "%d" % oligomaxtm
        cline.ogcopt = "%d" % oligoogcpercent
        cline.ogcmin = "%d" % oligomingc
        cline.ogcmax = "%d" % oligomaxgc
        cline.opolyxmax = "%d" % oligomaxpolyx
        cline.outfile = os.path.splitext(gd.seqfilename)[0] + '.eprimer3'
        gd.primerfilename = cline.outfile
        clines.append(str(cline) + log_output(gd.name + ".eprimer3"))
    print_verbose("... ePrimer3 jobs to run:")
    print_verbose_list(clines)
    # Parallelise jobs
    if not sge:
        multiprocessing_run(clines, poolsize, verbose)
    else:
        sge_run(clines, verbose)


# Load primers from ePrimer3 files into each GenomeData object
def load_primers(gdlist, minlength, nocds=True, filtergc3prime=False,
                 hybridprobe=False, verbose=False):
    """ Load primer data from an ePrimer3 output file into a dictionary of
        Bio.Emboss.Primer3.Primer objects (keyed by primer name) in a
        GenomeData object, for each such object in the passed list.
        Each primer object is given a new ad hoc attribute 'cds_overlap' which
        takes a Boolean, indicating whether the primer is found wholly within
        a CDS defined in the GenomeData object's feature file; this status
        is determined using an interval tree approach.
    """
    if verbose:
        t0 = time.time()
        print "Loading primers, %sfiltering on CDS overlap" %\
            ('not ' if nocds else '')
    # Load in the primers, assigning False to a new, ad hoc attribute called
    # cds_overlap in each
    for gd in gdlist:
        print_verbose("... loading primers into %s from %s ..." %\
                (gd.name, gd.primerfilename))
        try:
            os.path.isfile(gd.primerfilename)
        except TypeError:
            raise IOError("Primer file %s does not exist." % gd.primerfilename)
        primers = Primer3.read(open(gd.primerfilename, 'rU')).primers
        # Add primer pairs to the gd.primers dictionary
        primercount = 0
        for primer in primers:
            primercount += 1
            primer.cds_overlap = False           # default state
            primer.name = "%s_primer_%04d" % (gd.name, primercount)
            primer.amplifies_organism = set()    # Organisms amplified
            primer.amplifies_family = set()      # Organism families amplified
            primer.gc3primevalid = True          # Passes GC 3` test
            primer.oligovalid = True             # Oligo passes filter
            primer.blastpass = True              # Primers pass BLAST screen
            gd.primers.setdefault(primer.name, primer)
            primer.amplicon = \
                gd.sequence[primer.forward_start - 1:\
                                primer.reverse_start - 1 + primer.reverse_length]
            primer.amplicon.description = primer.name
        print_verbose("... loaded %d primers into %s ..." %\
                (len(gd.primers), gd.name))
        # Now that the primers are in the GenomeData object, we can filter
        # them on location, if necessary
        if not nocds:
            filter_primers(gd, minlength, verbose)
        # We also filter primers on the basis of GC presence at the 3` end
        if filtergc3prime:
            filter_primers_gc_3prime(gd, verbose)
        # Filter primers on the basis of internal oligo characteristics
        if hybridprobe:
            filter_primers_oligo(gd, verbose)


# Filter primers in a passed gd object on the basis of CDS features
def filter_primers(gd, minlength, verbose):
    """ Takes a passed GenomeData object, and the minimum size of an amplified
        region, and then uses a ClusterTree to find clusters of CDS and
        primer regions that overlap by this minimum size.
        There is a possibility that, by stacking primer regions, some of
        the reported overlapping primers may in fact not overlap CDS regions
        directly, so this function may overreport primers.
    """
    # Load in the feature data.  This is done using either SeqIO for
    # files with the .gbk extension, or an ad hoc parser for
    # .prodigalout prediction files
    if verbose:
        t0 = time.time()
        print "Loading feature data from %s ..." % gd.ftfilename
    if os.path.splitext(gd.ftfilename)[-1] == '.gbk':  # GenBank
        seqrecord = [r for r in SeqIO.parse(open(gd.ftfilename, 'rU'),
                                            'genbank')]
    elif os.path.splitext(gd.ftfilename)[-1] == '.prodigalout':
        seqrecord = parse_prodigal_features(gd.ftfilename, verbose)
    else:
        raise IOError, "Expected .gbk or .prodigalout file extension"
    print_verbose("... loaded %d features ..." % len(seqrecord.features))
    # Use a ClusterTree as an interval tree to identify those
    # primers that overlap with features.  By setting the minimum overlap to
    # the minimum size for a primer region, we ensure that we capture every
    # primer that overlaps a CDS feature by this amount, but we may also
    # extend beyond the CDS by stacking primers, in principle.
    print_verbose("... adding CDS feature locations to ClusterTree ...")
    ct = ClusterTree(-minlength, 2)
    # Loop over CDS features and add them to the tree with ID '-1'.  This
    # allows us to easily separate the features from primers when reviewing
    # clusters.
    for ft in [f for f in seqrecord.features if f.type == 'CDS']:
        ct.insert(ft.location.nofuzzy_start, ft.location.nofuzzy_end, -1)
    # ClusterTree requires us to identify elements on the tree by integers,
    # so we have to relate each primer added to an integer in a temporary
    # list of the gd.primers values
    print_verbose("... adding primer locations to cluster tree ...")
    aux = gd.primers.values()
    for i in range(len(gd.primers)):
        ct.insert(aux[i].forward_start,
                  aux[i].reverse_start + aux[i].reverse_length,
                  i)
    # Now we find the overlapping regions, extracting all element ids that are
    # not -1.  These are the indices for aux, and we modify the gd.cds_overlap
    # attribute directly
    print_verbose("... finding overlapping primers ...")
    overlap_primer_ids = set()                         # CDS overlap primers
    for (s, e, ids) in ct.getregions():
        primer_ids = set([i for i in ids if i != -1])  # get non-feature ids
        overlap_primer_ids = overlap_primer_ids.union(primer_ids)
    print_verbose("... %d primers overlap CDS features (%.3fs)" %\
            (len(overlap_primer_ids), time.time() - t0))
    for i in overlap_primer_ids:
        aux[i].cds_overlap = True


# Filter primers on the basis of GC content at 3` end
def filter_primers_gc_3prime(gd, verbose):
    """ Loops over the primer pairs in the passed GenomeData object and,
        if either primer has more than 2 G+C in the last five nucleotides,
        sets the .gc3primevalid flag to False.
    """
    if verbose:
        t0 = time.time()
        print "Filtering %s primers on 3` GC content ..." % gd.name
    invalidcount = 0
    for primer in gd.primers.values():
        fseq, rseq = primer.forward_seq[-5:], primer.reverse_seq[-5:]
        if (fseq.count('C') + fseq.count('G') > 2) or \
                (rseq.count('C') + fseq.count('G') > 2):
            primer.gc3primevalid = False
            invalidcount += 1
    print_verbose("... %d primers failed (%.3fs)" % (invalidcount, time.time() - t0))


# Filter primers on the basis of internal oligo characteristics
def filter_primers_oligo(gd, verbose):
    """ Loops over the primer pairs in the passed GenomeData object and,
        mark the primer.oligovalid as False if the internal oligo corresponds
        to any of the following criteria:
        - G at 5` end or 3` end
        - two or more counts of 'CC'
        - G in second position at 5` end
    """
    if verbose:
        t0 = time.time()
        print "Filtering %s primers on internal oligo characteristics ..." %\
            gd.name
    invalidcount = 0
    for primer in gd.primers.values():
        if (primer.oligo.seq.startswith('G') or \
                primer.oligo.seq.endswith('G') or \
                primer.oligo.seq[1:-1].count('CC') > 1 or \
                primer.oligo.seq[1] == 'G'):
            primer.oligovalid = False
            invalidcount += 1
    print_verbose("... %d primers failed (%.3fs)" % (invalidcount, time.time() - t0))


# Screen passed GenomeData primers against BLAST database
def blast_screen(gdlist, blast_exe, blastdb, cpus, sge, verbose):
    """ The BLAST screen takes three stages.  Firstly we construct a FASTA
        sequence file containing all primer forward and reverse sequences,
        for all primers in each GenomeData object of the list.
        We then use the local BLAST+ (not legacy BLAST) interface to BLASTN to
        query the named database with the input file.  The multiprocessing
        of BLASTN is handled by either our multiprocessing threading approach,
        or by SGE; we don't use the built-in threading of BLAST so that we
        retain flexibility when moving to SGE.  It's a small modification to
        revert to using the BLAST multithreading code.  The output file is
        named according to the GenomeData object.
        The final step is to parse the BLAST output, and label the primers
        that make hits as not having passed the BLAST filter.
    """
    build_blast_input(gdlist, verbose)
    run_blast(gdlist, blast_exe, blastdb, cpus, sge, verbose)
    parse_blast(gdlist, cpus, verbose)


# Write BLAST input files for each GenomeData object
def build_blast_input(gdlist, verbose):
    """ Loops over each GenomeData object in the list, and writes forward
        and reverse primer sequences out in FASTA format to a file with
        filename derived from the GenomeData object name.
    """
    if verbose:
        t0 = time.time()
        print "Writing files for BLAST input ..."
    for gd in gdlist:
        gd.blastinfilename = os.path.join(os.path.split(gd.seqfilename)[0],
                                          "%s_BLAST_input.fas" % gd.name)
        seqrecords = []
        for name, primer in gd.primers.items():
            seqrecords.append(SeqRecord(Seq(primer.forward_seq), 
                                        id = name + '_forward'))
            seqrecords.append(SeqRecord(Seq(primer.reverse_seq), 
                                        id = name + '_reverse'))
        print_verbose("... writing %s ..." % gd.blastinfilename)
        SeqIO.write(seqrecords,
                    open(gd.blastinfilename, 'w'),
                    'fasta')
    print_verbose("... done (%.3fs)" % (time.time() - t0))


# Run BLAST screen for each GenomeData object
def run_blast(gdlist, blast_exe, blastdb, poolsize, sge, verbose):
    """ Loop over the GenomeData objects in the passed list, and run a
        suitable BLASTN query with the primer sequences, writing to a file
        with name derived from the GenomeData object, in XML format.
    """
    if verbose:
        t0 = time.time()
        print "Compiling BLASTN command-lines ..."
    clines = []
    for gd in gdlist:
        gd.blastoutfilename = os.path.join(os.path.split(gd.seqfilename)[0],
                                          "%s_BLAST_output.xml" % gd.name)
        cline = NcbiblastnCommandline(query=gd.blastinfilename,
                                            db=blastdb,
                                            task='blastn',  # default: MEGABLAST
                                            out=gd.blastoutfilename,
                                            num_alignments=1,
                                            num_descriptions=1,
                                            outfmt=5,
                                            perc_identity=90,
                                            ungapped=True)
        clines.append(str(cline) + log_output(gd.name + ".blastn"))
    print_verbose("... BLASTN+ jobs to run:")
    print_verbose_list(clines)
    if not sge:
        multiprocessing_run(clines, poolsize, verbose)
    else:
        sge_run(clines, verbose)


# Parse BLAST output for each GenomeData object
def parse_blast(gdlist, poolsize, verbose):
    """ Loop over the GenomeData objects in the passed list, and parse the
        BLAST XML output indicated in the .blastoutfilename attribute.
        For each query that makes a suitable match, mark the appropriate
        primer's .blastpass attribute as False
    """
    if verbose:
        t0 = time.time()
        print "Parsing BLASTN output with multiprocessing ..."
    # Here I'm cheating a bit and using multiprocessing directly so that
    # we can speed up the parsing process a bit
    pool = multiprocessing.Pool(processes=poolsize)
    pool_results = [pool.apply_async(process_blastxml,
                                     (gd.blastoutfilename, gd.name, verbose)) \
                        for gd in gdlist]
    pool.close()
    pool.join()
    # Process the results returned from the BLAST searches.  Create a
    # dictionary of GenomeData objects, keyed by name, and loop over the
    # result sets, setting .blastpass attributes for the primers as we go
    gddict = {}
    [gddict.setdefault(gd.name, gd) for gd in gdlist]
    failcount = 0
    for r in [r.get() for r in pool_results]:
        for name in r:
            gd = gddict[name.split('_primer_')[0]]
            gd.primers[name].blastpass = False
            failcount += 1
    print_verbose("... %d primers failed BLAST screen ..." % failcount)
    print_verbose("... multiprocessing BLAST parsing complete (%.3fs)" % \
            (time.time() - t0))


# BLAST XML parsing function for multiprocessing
def process_blastxml(filename, name, verbose):
    """ Takes a BLAST output file, and a process name as input.  Returns
        a set of query sequence names that make a suitably strong hit to
        the database.
        We are using the database as a screen, so *any* hit that passes
        our criteria will do; BLAST+ reports the hits in quality order, so
        we only need to see this top hit.
        We care if the screening match is identical for at least 90% of
        the query, and we're using ungapped alignments, so we check
        the alignment HSP identities against the length of the query.
    """
    if verbose:
        t0 = time.time()
        print "[process name: %s] Parsing BLAST XML ..." % name
    # List to hold queries that hit the database
    matching_primers = set()
    recordcount = 0
    # Parse the file
    try:
        for record in NCBIXML.parse(open(filename, 'rU')):
            recordcount += 1         # Increment our count of matches
            # We check whether the number of identities in the alignment is
            # greater than our (arbitrary) 90% cutoff.  If so, we add the
            # query name to our set of failing/matching primers
            if len(record.alignments):
                identities = float(record.alignments[0].hsps[0].identities) / \
                    float(record.query_letters)
                if 0.9 <= identities:
                    matching_primers.add('_'.join(record.query.split('_')[:-1]))
        print_verbose("[process name: %s] Parsed %d records" % (name,
                                                                recordcount))
    except:
        print_verbose("[process name: %s] Error reading BLAST XML file" % name)
    print_verbose("[process name: %s] Time spent in process: (%.3fs)" % \
            (name, time.time() - t0))
    # Return the list of matching primers
    return matching_primers


# A function for parsing features from Prodigal output
def parse_prodigal_features(filename, verbose):
    """ Parse Prodigal 'GenBank' output.
        We try to emulate SeqIO.read() SeqRecord output as much as possible,
        but the information provided by Prodigal is limited to feature type
        and location, on a single line.
        Amended: Newer versions of Prodigal write closer match to GenBank
        format, and thus if the first line matches "DEFINITION" we use SeqIO.
        RE-amended: Latest version of Prodigal is still not good enough for SeqIO,
        so a new function is created to parse line-by-line.
    """
    record = SeqRecord(None)       # record gets a dummy sequence
    # Open filehandle and parse contents
    handle = open(filename, 'rU')
    # init feature list from file parsing
    record.features = seqrecord_parse(handle)
    #record.features = seqrecord_parse_seqio(handle, 'genbank')  # to be used later (see NOTE)
    return record


# Parse record features from the lines of prodigal or genbank format file
def seqrecord_parse(filehandle):
    features = []
    for line in filehandle:                             # for each line in file
        if (re.search("CDS", line)):                    # check if line has CDS as type
            data = [e.strip() for e in line.split()]    # split elements
            f = gb_string_to_feature(data[-1])          # last element is 'SeqFeature'
            f.type = data[0]                            # first elem (CDS) is type
            features.append(f)                          # append to feature list
    #print_debug("Number of features: %d" % len(features))
    return features                                     # return feature list


# Parse record features from sequence file, using SeqIO
def seqrecord_parse_seqio(filehandle, seqformat):
# NOTE: Latest version of prodigal output is *closer* to GenBank format
#       but not close enough for SeqIO to find the genome.features
#       Thus: this function NOT USED (until potential update to prodigal
#       or SeqIO).
    features = []
    seqrecord = list(SeqIO.parse(filehandle, seqformat))
    #print_debug("record[0].id: %s" % seqrecord[0].id)

    for r in seqrecord:
        print_debug("record seq: [%s]..." % r.seq[0:12])
        #print_debug("record info: " + r)
        features.append(r.features)

    #print_debug("Seqrecord Features: %d" % len(features))
    return features


# Code (admittedly hacky) from Brad Chapman to parse a GenBank command line
def gb_string_to_feature(content, use_fuzziness=True):
    """Convert a GenBank location string into a SeqFeature.
    """
    consumer = _FeatureConsumer(use_fuzziness)
    consumer._cur_feature = SeqFeature()
    consumer.location(content)
    return consumer._cur_feature


# Run PrimerSearch all-against-all on a list of GenomeData objects
def primersearch(gdlist, poolsize, mismatchpercent, sge, verbose):
    """ Loop over the GenomeData objects in the passed list, and construct
        command lines for an all-against-all PrimerSearch run.
        Output files are of the format
        <query name>_vs_<target name>.primersearch
        Where <query name> and <target name> are the gd.name attributes of
        the source and target GenomeData objects, respectively.
        The output file goes in the same location as the source sequence
        file.
    """
    if verbose:
        t0 = time.time()
        print "Constructing all-against-all PrimerSearch runs " +\
            "for %d objects ..." % len(gdlist)
    # Create list of command-lines
    clines = []
    for query_gd in gdlist:
        query_gd.primersearch_output = []
        for target_gd in gdlist:
            if query_gd != target_gd:
                # Location of PrimerSearch output
                outdir = os.path.split(query_gd.seqfilename)[0]
                outfilename = os.path.join(outdir, "%s_vs_%s.primersearch" %\
                                               (query_gd.name, target_gd.name))
                query_gd.primersearch_output.append(outfilename)
                # Create command-line
                cline = PrimerSearchCommandline()
                cline.auto = True
                cline.seqall = target_gd.seqfilename
                cline.infile = query_gd.primersearchfilename
                cline.outfile = outfilename
                cline.mismatchpercent = mismatchpercent
                clines.append(str(cline) + log_output(os.path.basename(outfilename)))
    print_verbose("... PrimerSearch jobs to run: ...")
    print_verbose_list(clines)
    # Parallelise jobs
    if not sge:
        multiprocessing_run(clines, poolsize, verbose)
    else:
        sge_run(clines, verbose)


# Load in existing PrimerSearch output
def load_existing_primersearch_results(gdlist, verbose):
    """ Associates PrimerSearch output files with each GenomeData object
        and returns a list of (name, filename) tuples for all GenomeData
        objects
    """
    if verbose:
        t0 = time.time()
        print "Locating existing PrimerSearch input files ..."
    primersearch_results = []
    for gd in gdlist:
        gd.primersearch_output = []
        filedir = os.path.split(gd.seqfilename)[0]
        primersearch_files = [f for f in os.listdir(filedir) if \
                       os.path.splitext(f)[-1] == '.primersearch' and \
                                  f.startswith(gd.name)]
        for filename in primersearch_files:
            print_verbose("... found %s for %s ..." % (filename, gd.name))
            gd.primersearch_output.append(os.path.join(filedir,
                                                       filename))
    print_verbose("... found %d PrimerSearch input files (%.3fs)" % \
            (len(primersearch_results), time.time() - t0))


# Run primersearch to find whether and where the predicted primers amplify our negative target
# (the one we expect exactly one match to)
def find_negative_target_products(gdlist, filename, mismatchpercent, cpus, sge, verbose):
    """ We run primersearch using the predicted primers as queries, with the passed filename
        as the target sequence.  We exploit multiprocessing, and use the prescribed number of
        CPUs.  Happily, primersearch accepts multiple sequence FASTA files.
    """
    if verbose:
        t0 = time.time()
        print "Constructing negative control PrimerSearch runs " +\
            "for %d objects ..." % len(gdlist)
    # Create list of command-lines
    clines = []
    for query_gd in gdlist:
        query_gd.primersearch_output = []
        outdir = os.path.split(query_gd.seqfilename)[0]  # Primersearch output goes here
        outfilename = os.path.join(outdir, "%s_negative_control.primersearch" %\
                                       query_gd.name)
        query_gd.primersearch_output.append(outfilename)
        # Create command-line
        cline = PrimerSearchCommandline()
        cline.auto = True
        cline.seqall = filename
        cline.infile = query_gd.primersearchfilename
        cline.outfile = outfilename
        cline.mismatchpercent = mismatchpercent
        clines.append(str(cline) + log_output(os.path.basename(outfilename)))
    print_verbose("... PrimerSearch jobs to run: ...")
    print_verbose_list(clines)
    # Parallelise jobs and run
    if not sge:
        multiprocessing_run(clines, cpus, verbose)
    else:
        sge_run(clines, verbose)


# Classify the primers in a list of GenomeData objects according to the
# other sequences that they amplify
def classify_primers(gdlist, single_product, verbose):
    """ Takes a list of GenomeData objects and loops over the primersearch
        results, loading in the primersearch results and applying them to the
        associated query GenomeData object.
        If a primer is reported, by PrimerSearch, to amplify a region of the
        target genome, two changes are made to the corresponding Primer
        object in the amplifies_object and amplifies_family ad hoc attributes,
        with the target name and family, respectively, being added to those
        sets.
    """
    if verbose:
        t0 = time.time()
        print "Classifying primers by PrimerSearch results ..."
    # Convenience dictionary, keying each GenomeData object by name
    gddict = {}
    [gddict.setdefault(gd.name, gd) for gd in gdlist]
    # Parse the PrimerSearch output, updating the primer contents of the
    # appropriate GenomeData object, for each set of results
    for gd in gdlist:
        print_verbose("... GenomeData for %s ..." % gd.name)
        for filename in gd.primersearch_output:
            print_verbose("... processing %s ..." % filename)
            # Identify the target organism
            targetname = \
             os.path.splitext(os.path.split(filename)[-1])[0].split('_vs_')[-1]
            # Only classify amplimers to sequences in the gdlist dataset
            # This avoids problems with recording counts of matches to sequences
            # that we're not considering, artifically lowering the specificity
            # counts.
            if targetname in gddict:
                # Load the contents of the PrimerSearch output
                psdata = PrimerSearch.read(open(filename, 'rU'))
                # We loop over each primer in psdata and, if the primer has a
                # length this indicates that it amplifies the target.  When
                # this is the case we add the organism name and the family
                # name to the appropriate primer in the query GenomeData object
                for pname, pdata in psdata.amplifiers.items():
                    if len(pdata):
                        # Primer amplifies
                        gd.primers[pname].amplifies_organism.add(targetname)
                        gd.primers[pname].amplifies_family.add(gddict[targetname].family)
            # Consider the negative control primersearch output
            elif 'negative_control' in filename:
                # Load PrimerSearch data
                psdata = PrimerSearch.read(open(filename, 'rU'))
                # We loop over each primer, and find the number of amplimers.  We note the
                # number of amplimers as an attribute of the primer
                for pname, pdata in psdata.amplifiers.items():
                    gd.primers[pname].negative_control_amplimers = len(pdata)
                    print_verbose("Found %d amplimers in negative control" % len(pdata))
        print_verbose("... processed %d Primersearch results for %s ..." % \
                (len(gd.primersearch_output), gd.name))
    print_verbose("... processed PrimerSearch results (%.3fs)" % (time.time() - t0))


# Write analysis data to files
def write_report(gdlist, nocds, gc3primevalid, hybridprobe,
                 blastfilter, single_product, outdir, verbose):
    """ Write a tab-separated table of information to the passed
        filename, summarising the distribution of unique, family-unique,
        and universal (for this set) primers amongst the GenomeData objects
        in gdlist. Also write out to this file the locations of the files
        containing the data used to generate the information.
        In addition, write out the following files in ePrimer3 format:
        i) <query_name>_specific.eprimer3 - unique primers for each query
           GenomeData object
        ii) <family>_specific.eprimer3 - unique primers for each family in
           the GenomeData set
        iii) universal_primers.eprimer3 - primers that amplify all members of
           the GenomeData set
    """
    if verbose:
        t0 = time.time()
        print "Creating summary output ..."
    # First we need to generate a dictionary of GenomeData object names, keyed
    # by family
    families = defaultdict(set)
    for gd in gdlist:
        families[gd.family].add(gd.name)
    # Rectify nocds flag
    cds_overlap = not nocds
    # Check whether output directory exists and, if not, create it
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    # Open output file, and write header
    outfh = open(os.path.join(outdir, 'differential_primer_results.tab'), 'w')
    print >> outfh, os.linesep.join(["# Summary information table",
                                     "# Generated by find_differential_primers",
                                     "# Columns in the table:",
                                     "# 1) Query organism ID",
                                     "# 2) Query organism family",
                                     "# 3) Count of organism-unique primers",
                                     "# 4) Count of family-unique primers",
                                     "# 5) Count of universal primers",
                                     "# 6) Query sequence filename",
                                     "# 7) Query feature filename",
                                     "# 8) Query ePrimer3 primers filename"])
    # Write data for each GenomeData object
    other_org_count = len(gdlist) - 1  # Amplifications for 'universal' set
    # We store 'universal' primers in their own list, and family-specific
    # primers in a dicitonary, keyed by family
    all_universal_primers = []
    family_specific_primers = defaultdict(list)
    # Loop over each GenomeData object and populate family-specific and
    # universal primer collections, as well as organism-specific and
    # summary information
    for gd in gdlist:
        print_verbose_list([
                "... writing data for %s ..." % gd.name,
                "... cds_overlap: %s ..." % cds_overlap,
                "... gc3primevalid: %s ..." % gc3primevalid,
                "... oligovalid: %s ..." % hybridprobe,
                "... blastpass: %s ..." % blastfilter,
                "... single_product %s ..." % (single_product is not None),
                "... retrieving primer pairs ...",
                # Get the unique, family-specific and universal primers
                "... finding strain-specific primers for %s ..." % gd.name
                ])

        unique_primers = gd.get_unique_primers(cds_overlap, gc3primevalid,
                                               hybridprobe, blastfilter,
                                               single_product,
                                               verbose)
        print_verbose("... finding family-specific primers for %s ..." % gd.name)
        family_unique_primers = \
            gd.get_family_unique_primers(families[gd.family], cds_overlap,
                                         gc3primevalid, hybridprobe,
                                         blastfilter, single_product,
                                         verbose)
        family_specific_primers[gd.family] += family_unique_primers
        print_verbose("... finding universal primers for %s ..." % gd.name)
        universal_primers = \
            gd.get_primers_amplify_count(other_org_count, cds_overlap,
                                         gc3primevalid, hybridprobe,
                                         blastfilter, single_product,
                                         verbose)
        all_universal_primers += universal_primers
        # Write summary data to file
        print >> outfh, '\t'.join([gd.name, gd.family,
                                   str(len(unique_primers)),
                                   str(len(family_unique_primers)),
                                   str(len(universal_primers)),
                                   str(gd.seqfilename),
                                   str(gd.ftfilename),
                                   str(gd.primerfilename)])
        # Write organism-specific primers to file
        write_eprimer3(unique_primers,
                       os.path.join(outdir, "%s_specific_primers.eprimer3" %\
                                        gd.name),
                       gd.seqfilename, verbose=verbose)
        # Write organism-specific amplicons to file
        SeqIO.write([p.amplicon for p in unique_primers],
                    os.path.join(outdir, "%s_specific_amplicons.fas" % gd.name),
                    'fasta')
    # Write universal primers to file
    write_eprimer3(universal_primers,
                   os.path.join(outdir, "universal_primers.eprimer3"),
                   '', append=True, verbose=verbose)
    # Write organism-specific amplicons to file
    SeqIO.write([p.amplicon for p in universal_primers],
                open(os.path.join(outdir, "universal_amplicons.fas"), 'w'),
                'fasta')
    # Write family-specific primers to file
    for family, primers in family_specific_primers.items():
        write_eprimer3(primers,
                       os.path.join(outdir,
                                    "%s_family-specific_primers.eprimer3" % family),
                       '', verbose=verbose)
        # Write family-specific amplicons to file
        SeqIO.write([p.amplicon for p in primers],
                    open(os.path.join(outdir, "%s_family-specific_amplicons.fas" % \
                                          family), 'w'),
                    'fasta')
    # Being tidy...
    outfh.close()
    print_verbose("... data written (%.3fs)" % (time.time() - t0))


# Write ePrimer3 format primer file
def write_eprimer3(primers, filename, sourcefilename, append=False, verbose=False):
    """
    """
    if verbose:
        t0 = time.time()
        print "Writing %d primer pairs to %s ..." % (len(primers), filename)
    # Open file
    filemode = 'a' if append else 'w'     # Do we append or write anew?
    outfh = open(filename, 'w')
    # Write header
    print >> outfh, os.linesep.join([
                "# EPRIMER3 PRIMERS %s " % filename,
                "#                      Start  Len   Tm     GC%   Sequence",
                os.linesep])
    primercount = 0
    for p in primers:
        primercount += 1
        print >> outfh, "# %s %s" % (p.name, sourcefilename)
        print >> outfh, "%-4d PRODUCT SIZE: %d" % (primercount,
                                                    p.size)
        print >> outfh, "     FORWARD PRIMER  %-9d  %-3d  %.02f  %.02f  %s" %\
            (p.forward_start, p.forward_length, p.forward_tm, p.forward_gc,
             p.forward_seq)
        print >> outfh, "     REVERSE PRIMER  %-9d  %-3d  %.02f  %.02f  %s" %\
            (p.reverse_start, p.reverse_length, p.reverse_tm, p.reverse_gc,
             p.reverse_seq)
        if hasattr(p, 'internal_start'):
            print >> outfh, "     INTERNAL OLIGO  %-9d  %-3d  %.02f  %.02f  %s" %\
                (p.internal_start, p.internal_length, p.internal_tm, p.internal_gc,
                 p.internal_seq)
        print >> outfh, os.linesep * 2
    # Be tidy
    outfh.close()


# Run the passed list of command-lines using a multiprocessing.Pool
def multiprocessing_run(clines, poolsize, verbose):
    """ We create a multiprocessing Pool to handle command-lines  We
        pass the (unique) GenomeData object name, and the location of the
        sequence file.  The called function returns the GenomeData name and the
        corresponding location of the generated feature file.  The GenomeData
        objects are stored in a temporary dictionary, keyed by gd.name, to allow
        association of the results of the asynchronous pool jobs with the
        correct GenomeData object
    """
    if verbose:
        t0 = time.time()
        print "Running %d jobs with multiprocessing ..." % \
            len(clines)
    pool = multiprocessing.Pool(processes=poolsize)      # create process pool
    completed = []
    if verbose:
        callback_fn = multiprocessing_callback
    else:
        callback_fn = completed.append
    pool_outputs = [pool.apply_async(subprocess.call,
                                     (str(cline), ),
                                     {'stderr': subprocess.PIPE,
                                      'shell': sys.platform != "win32"},
                                     callback=callback_fn) \
                        for cline in clines]
    pool.close()      # Run jobs
    pool.join()
    print_verbose_list(completed)
    print_verbose("... all multiprocessing jobs ended (%.3fs)" % (time.time() - t0))


# Add a multiprocessing callback function here
def multiprocessing_callback(val):
    """ A verbose callback function for multiprocessing runs.  It uses the
        return value to indicate run completion or failure.  Failure is
        indicated by a nonzero return from the multiprocessing call.
    """
    if 0 == val:
        print >> sys.stderr, \
            "... multiprocessing run completed (status: %s) ..." % val
    else:
        print >> sys.stderr, \
            "... problem with multiprocessing run (status: %s) ..." % val


# Clean output for each GenomeData object in the passed list
def clean_output(gdlist, verbose):
    """ Remove .eprimer3, .primers, .prodigalout, and .primersearch files
        from the same directory as the sequence file for each passed PrimerSearch
        object
    """
    if verbose:
        t0 = time.time()
        print "Cleaning up output files for GenomeData objects ..."
    # Loop over each GenomeData object, and remove each output file
    for gd in gdlist:
        seqdir = os.path.split(gd.seqfilename)[0]
        for filename in [f for f in os.listdir(seqdir) \
                             if os.path.splitext(f)[-1] in \
                             ['.eprimer3', 'primers', '.prodigalout',
                              '.primersearch', '.xml']]:
            abspath = os.path.join(seqdir, filename)
            print_verbose("... deleting %s ..." % abspath)
            os.remove(abspath)     # You can never go back after this point
    print_verbose("... done (%.3fs)" % (time.time() - t0))


# Print if --debug was in command line.
def print_debug(string):
    """ check options.debug, print if exists
    """
    if options.debug:
        #print string
        print "DEBUG: ", string


# construct str to concat on end of cline if option.keep_logs is set
def log_output(filename):
    """ predefine file extension and stream to print to.
        if log_dir exists, join it to filename
        else output to base filename.
    """
    log_extension = ".log"
    log_out_handle = " 2> "
    if options.keep_logs and options.log_dir:
        return log_out_handle + os.path.join(options.log_dir, filename) + log_extension
    elif options.keep_logs:
        return log_out_handle + filename + log_extension
    else:
        return ""


# Print if --verbose was in command line.
def print_verbose(string):
    # check options.verbose, print if exists
    if options.verbose:
        print string


# Print list elements if --verbose was in command line.
def print_verbose_list(*lists):
    if options.verbose:
        for lst in lists:
            print "\n".join(str(e) for e in lst)


###
# SCRIPT
if __name__ == '__main__':
    # Parse cmd-line
    options, args = parse_cmdline(sys.argv)
    print_verbose("Script called with options: %s" % ' '.join(sys.argv[1:]))
    # Create our GenomeData objects.  If there is no configuration file
    # specified, raise an error and exit.  Otherwise we end up with a list
    # of GenomeData objects that are populated only with the data from the
    # config file
    if options.filename is None:
        raise IOError("No configuration file specified")
    gdlist = create_gd_from_config(options.filename, options.verbose)
    # If the user wants to clean the directory before starting, do so
    if options.clean or options.cleanonly:
        clean_output(gdlist, options.verbose)
    if options.cleanonly:
        sys.exit(0)
    # It is possible that the sequence file for a GenomeData object might
    # be a multi-sequence file describing scaffolds or contigs.  We create a
    # concatenated sequence to facilitate further analyses, if this is the
    # case.  Where a sequence needs to be concatenated, this will affect the
    # placement of features and/or primers, so any specified files are
    # reset to None
    check_single_sequence(gdlist, options.verbose)
    # We need to check the existence of a prescribed feature file and, if
    # there is not one, create it.  We don't bother if the --nocds flag is set.
    if not (options.nocds or options.noprodigal):
        print_verbose("--nocds option not set: Checking existence of features...")
        check_ftfilenames(gdlist, options.prodigal_exe, options.cpus,
                          options.sge, options.verbose)
    elif options.nocds:
        print_verbose("--nocds option set: Not checking or creating feature files")
    else:
        print_verbose("--noprodigal option set: Not predicting new CDS")
    # We need to check for the existence of primer sequences for the organism
    # and, if they do not exist, create them using ePrimer3.  If the
    # --noprimer3 flag is set, we do not create new primers, but even if the
    # --noprimersearch flag is set, we still need to check whether the
    # primer files are valid
    if not options.noprimer3:
        print_verbose("--noprimer3 flag not set: Predicting new primers")
        check_primers(gdlist, options.verbose)
        predict_primers(gdlist, options.eprimer3_exe,
                        options.cpus, options.numreturn,
                        options.osize, options.minsize, options.maxsize,
                        options.otm, options.mintm, options.maxtm,
                        options.ogcpercent, options.mingc, options.maxgc,
                        options.psizeopt, options.psizemin, options.psizemax,
                        options.maxpolyx,
                        options.hybridprobe, options.oligoosize,
                        options.oligominsize, options.oligomaxsize,
                        options.oligootm, options.oligomintm,
                        options.oligomaxtm,
                        options.oligoogcpercent, options.oligomingc,
                        options.oligomaxgc, options.oligomaxpolyx,
                        options.sge, options.verbose)
    else:
        print_verbose("--noprimer3 flag set: Not predicting new primers")
    # With a set of primers designed for the organism, we can load them into
    # the GenomeData object, filtering for those present only in the CDS,
    # if required.  This step is necessary, whether or not a new ePrimer3
    # prediction is made.  We also filter on GC content at the primer 3' end,
    # if required.
    print_verbose("Loading primers...")
    load_primers(gdlist, options.psizemin, options.nocds,
                 options.filtergc3prime, options.hybridprobe, options.verbose)
    # At this point, we can check our primers against a prescribed BLAST
    # database.  How we filter these depends on the user's preference.
    # We screen against BLAST here so that we can flag an attribute on
    # each primer to say whether or not it passed the BLAST screen.
    if options.blastdb and not options.useblast:
        print_verbose("--blastdb options set: BLAST screening primers...")
        blast_screen(gdlist, options.blast_exe, options.blastdb,
                     options.cpus, options.sge,
                     options.verbose)
    elif options.useblast:
        print_verbose("--useblast options set: using existing BLAST results...")
    else:
        print_verbose("No BLAST options set, not BLAST screening primers...")
    # Having a set of (potentially CDS-filtered) primers for each organism,
    # we then scan these primers against each of the other organisms in the
    # set, using the EMBOSS PrimerSearch package
    # (http://embossgui.sourceforge.net/demo/manual/primersearch.html)
    # Now we have all the data we need to run PrimerSearch in an all-against-all
    # manner, so make a cup of tea, put your feet up, and do the comparisons
    # with EMBOSS PrimerSearch
    # (http://embossgui.sourceforge.net/demo/manual/primersearch.html)
    if options.noprimersearch:
        print_verbose("--noprimersearch flag set: Not running PrimerSearch")
        # Load the appropriate primersearch output files for each
        # GenomeData object
        load_existing_primersearch_results(gdlist, options.verbose)
    else:
        print_verbose("--noprimersearch flag not set: Running PrimerSearch")
        # We write input for PrimerSearch ignoring all the filters; this lets
        # us turn off PrimerSearch and rerun the analysis with alternative filter
        # settings
        for gd in gdlist:
            gd.write_primers(options.verbose)
        # Run PrimerSearch
        primersearch(gdlist, options.cpus, options.mismatchpercent,
                     options.sge, options.verbose)
    # If the --single_product option is specified, we load in the sequence file to which the
    # passed argument refers, and filter the primer sequences on the basis of how many
    # amplification products are produced from these sequences.  We expect exactly one
    # amplification product per primer set, if it's not degenerate on the target sequence
    # (note that this filter is meaningless for family-specific primers)
    if options.single_product:
        find_negative_target_products(gdlist, options.single_product, options.mismatchpercent,
                                      options.cpus, options.sge,
                                      options.verbose)
    # Now we classify the primer sets according to which sequences they amplify
    if not options.noclassify:
        print_verbose("Classifying primers and writing output files ...")
        # Classify the primers in each GenomeData object according to
        # the organisms and families that they amplify, using the
        # PrimerSearch results.
        classify_primers(gdlist, options.single_product, options.verbose)
        # All the data has been loaded and processed, so we can now create our
        # plaintext summary report of the number of unique, family-unique and
        # universal primers in each of the organisms
        write_report(gdlist, options.nocds, options.filtergc3prime,
                     options.hybridprobe,
                     (options.blastdb is not None or options.useblast),
                     options.single_product,
                     options.outdir, options.verbose)

