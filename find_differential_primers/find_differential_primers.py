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
# (c) The James Hutton Institute 2011
# Authors: Leighton Pritchard, Benjamin Leopold, Michael Robeson
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

# script version
# should match r"^__version__ = '(?P<version>[^']+)'$" for setup.py
__version__ = '0.1.2'


###
# IMPORTS

import logging
import logging.handlers
import multiprocessing
import os
import subprocess
import sys
import time
import traceback
import re

from collections import defaultdict             # Syntactic sugar
from optparse import OptionParser               # Cmd-line parsing

try:
    from Bio import SeqIO   # Parsing biological sequence data
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML                   # BLAST XML parser
    from Bio.Emboss.Applications import Primer3Commandline, \
        PrimerSearchCommandline
    from Bio.Emboss import Primer3, PrimerSearch  # EMBOSS parsers
    from Bio.GenBank import _FeatureConsumer      # For GenBank locations
    from Bio.Seq import Seq                       # Represents a sequence
    from Bio.SeqRecord import SeqRecord     # Represents annotated record
    from Bio.SeqFeature import SeqFeature   # Represents annotated record
except ImportError:
    sys.stderr.write("Biopython required for script, but not found (exiting)")
    sys.exit(1)

try:
    from bx.intervals.cluster import ClusterTree    # Interval tree building
except ImportError:
    sys.stderr.write("bx-python required for script, but not found (exiting)")
    sys.exit(1)


###
# CLASSES

# Class describing an organism's genome, and associated data.
class GenomeData(object):
    """ Describes an organism's genome, and has attributes:

        name   - short, unique (not enforced) identification string
        families - string indicating family memberships
        seqfilename    - location of representative genome sequence file
        ftfilename     - location of GBK/Prodigal feature file
        primerfilename - location of ePrimer3 format primers file

        primers     - dictionary collection of Bio.Emboss.Primer3.Primer
                      objects, keyed by primer name

        Exposed methods are:

    """
    def __init__(self, name, families=None, seqfilename=None,
                 ftfilename=None, primerfilename=None,
                 primersearchfilename=None):
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
        self.name = name  # Short identifier
        self.families = families.split(',') if families != '-' else None
        self.seqfilename = seqfilename if seqfilename != '-' else None
        self.ftfilename = ftfilename if ftfilename != '-' else None
        self.primerfilename = primerfilename if primerfilename != '-' \
            else None
        self.primersearchfilename = primersearchfilename if\
            primersearchfilename != '-' else None
        self.primers = {}  # Dict of Primer objects, keyed by name
        self.sequence = None  # Will hold genome sequence
        self.load_sequence()

    def load_sequence(self):
        """ Load the sequence defined in self.seqfile into memory.  We
            assume it's FASTA format.  This can then be used to calculate
            amplicons when loading primers in.
        """
        if self.seqfilename is not None:
            try:
                self.sequence = SeqIO.read(open(self.seqfilename, 'rU'),
                                           'fasta')
            except ValueError:
                logger.error("Loading sequence file %s failed",
                             self.seqfilename)
                logger.error(last_exception())
                sys.exit(1)

    def write_primers(self):
        """ Write the primer pairs in self.primers out to file in an
            appropriate format for PrimerSearch.  If the filename is not
            already defined, the filestem of the
            source sequencefile is used for the output file, with the
            extension '.primers'.
            The method returns the number of lines written.
        """
        # Define output filename, if not already defined
        if self.primersearchfilename is None:
            self.primersearchfilename = \
                os.path.splitext(self.seqfilename)[0] + '.primers'
        time_start = time.time()
        logger.info("Writing primers to file %s ...",
                    self.primersearchfilename)
        # Open handle and write data
        outfh = open(self.primersearchfilename, 'w')
        outfh.write("# Primers for %s\n" % self.name)
        outfh.write("# Automatically generated by find_differential_primers\n")
        for primers in self.primers.values():
            outfh.write("%s\t%s\t%s\n" %
                        (primers.name, primers.forward_seq,
                         primers.reverse_seq))
        if not len(self.primers):
            logger.warning("WARNING: no primers written to %s!",
                           self.primersearchfilename)
        # Being tidy
        outfh.close()
        logger.info("... wrote %d primers to %s (%.3fs)",
                    len(self.primers),
                    self.primersearchfilename, time.time() - time_start)

    def get_unique_primers(self, cds_overlap=False,
                           oligovalid=False,
                           blastfilter=False):
        """ Returns a list of primers that have the .amplifies_organism
            attribute, but where this is an empty set.
            If cds_overlap is True, then this list is restricted to those
            primers whose .cds_overlap attribute is also True
        """
        return self.get_primers_amplify_count(0, cds_overlap, blastfilter)

    def get_family_unique_primers(self, family_members, cds_overlap=False,
                                  blastfilter=False):
        """ Returns a list of primers that have the .amplifies_organism
            attribute, and where the set of organisms passed in family_members
            is the same as that in .amplifies_organism, with the addition of
            self.name.
            If cds_overlap is True, then this list is restricted to those
            primers whose .cds_overlap attribute is also True
        """
        primerlist = []
        for primer in self.primers.values():
            if family_members == \
               set([self.name]).union(primer.amplifies_organism):
                primerlist.append(primer)
        logger.info("[%s] %d family primers",
                    self.name, len(primerlist))
        if cds_overlap:
            primerlist = [p for p in primerlist if p.cds_overlap]
            logger.info("[%s] %d primers after CDS filter",
                        self.name, len(primerlist))
        if options.filtergc3prime:
            primerlist = [p for p in primerlist if p.gc3primevalid]
            logger.info("[%s] %d primers after GC 3` filter",
                        self.name, len(primerlist))
        if options.hybridprobe:
            primerlist = [p for p in primerlist if p.oligovalid]
            logger.info("[%s] %d primers after oligo filter",
                        self.name, len(primerlist))
        if blastfilter:
            primerlist = [p for p in primerlist if p.blastpass]
            logger.info("[%s] %d primers after BLAST filter",
                        self.name, len(primerlist))
        if options.single_product:
            primerlist = [p for p in primerlist if
                          p.negative_control_amplimers == 1]
            logger.info("[%s] %d primers after single_product filter",
                        self.name, len(primerlist))
        logger.info("[%s] returning %d primers",
                    self.name, len(primerlist))
        return primerlist

    def get_primers_amplify_count(self, count, cds_overlap=False,
                                  blastfilter=False):
        """ Returns a list of primers that have the .amplifies_organism
            attribute and the length of this set is equal to the passed count.
            If cds_overlap is True, then this list is restricted to those
            primers whose .cds_overlap attribute is also True
        """
        primerlist = [p for p in self.primers.values() if
                      count == len(p.amplifies_organism)]
        logger.info("[%s] %d family primers that amplify %d orgs",
                    self.name, len(primerlist), count)
        if cds_overlap:
            primerlist = [p for p in primerlist if p.cds_overlap]
            logger.info("[%s] %d primers after CDS filter",
                        self.name, len(primerlist))
        if options.filtergc3prime:
            primerlist = [p for p in primerlist if p.gc3primevalid]
            logger.info("[%s] %d primers after GC 3` filter",
                        self.name, len(primerlist))
        if options.hybridprobe:
            primerlist = [p for p in primerlist if p.oligovalid]
            logger.info("[%s] %d primers after oligo filter",
                        self.name, len(primerlist))
        if blastfilter:
            primerlist = [p for p in primerlist if p.blastpass]
            logger.info("[%s] %d primers after BLAST filter",
                        self.name, len(primerlist))
        if options.single_product:
            primerlist = [p for p in primerlist if
                          p.negative_control_amplimers == 1]
            logger.info("[%s] %d primers after single_product filter",
                        self.name, len(primerlist))
        logger.info("[%s] returning %d primers",
                    self.name, len(primerlist))
        return primerlist

    # Filter primers on the basis of CDS feature overlap
    def filter_primers(self, psizemin):
        """ Takes the minimum size of an amplified
            region, and then uses a ClusterTree to find clusters of CDS and
            primer regions that overlap by this minimum size.
            There is a possibility that, by stacking primer regions, some of
            the reported overlapping primers may in fact not overlap CDS
            regions directly, so this method may overreport primers.

            - psizemin (int): minimum size of an amplified region
        """
        # Load in the feature data.  This is done using either SeqIO for
        # files with the .gbk extension, or an ad hoc parser for
        # .prodigalout prediction files
        time_start = time.time()
        logger.info("Loading feature data from %s ...", self.ftfilename)
        if os.path.splitext(self.ftfilename)[-1] == '.gbk':  # GenBank
            seqrecord = [r for r in SeqIO.parse(open(self.ftfilename, 'rU'),
                                                'genbank')]
        elif os.path.splitext(self.ftfilename)[-1] == '.prodigalout':
            seqrecord = parse_prodigal_features(self.ftfilename)
        else:
            raise IOError("Expected .gbk or .prodigalout file extension")
        logger.info("... loaded %d features ...", len(seqrecord.features))

        # Use a ClusterTree as an interval tree to identify those
        # primers that overlap with features.  By setting the minimum overlap
        # to the minimum size for a primer region, we ensure that we capture
        # every primer that overlaps a CDS feature by this amount, but we may
        # also extend beyond the CDS by stacking primers, in principle.
        logger.info("... adding CDS feature locations to ClusterTree ...")
        ctree = ClusterTree(-psizemin, 2)

        # Loop over CDS features and add them to the tree with ID '-1'.  This
        # allows us to easily separate the features from primers when
        # reviewing clusters.
        for feature in [f for f in seqrecord.features if f.type == 'CDS']:
            ctree.insert(feature.location.nofuzzy_start,
                         feature.location.nofuzzy_end, -1)

        # ClusterTree requires us to identify elements on the tree by integers,
        # so we have to relate each primer added to an integer in a temporary
        # list of the self.primers values
        logger.info("... adding primer locations to cluster tree ...")
        aux = {}
        for i, primer in enumerate(self.primers.values()):
            ctree.insert(primer.forward_start,
                         primer.reverse_start + primer.reverse_length, i)
            aux[i] = primer

        # Now we find the overlapping regions, extracting all element ids
        # that are not -1.  These are the indices for aux, and we modify the
        # self.cds_overlap attribute directly
        logger.info("... finding overlapping primers ...")
        overlap_primer_ids = set()  # CDS overlap primers
        for (s, e, ids) in ctree.getregions():
            primer_ids = set([i for i in ids if i != -1])  # get non-ft ids
            overlap_primer_ids = overlap_primer_ids.union(primer_ids)
        logger.info("... %d primers overlap CDS features (%.3fs)",
                    len(overlap_primer_ids), time.time() - time_start)
        for i in overlap_primer_ids:
            aux[i].cds_overlap = True

    # Filter primers on the basis of internal oligo characteristics
    def filter_primers_oligo(self):
        """ Loops over the primer pairs in this GenomeData object and
            mark primer.oligovalid as False if the internal oligo corresponds
            to any of the following criteria:
            - G at 5` end or 3` end
            - two or more counts of 'CC'
            - G in second position at 5` end
        """
        time_start = time.time()
        logger.info("Filtering %s primers on internal oligo...",
                    self.name)
        invalidcount = 0
        for primer in self.primers.values():
            primer.oligovalid = not(primer.internal_seq.startswith('G')
                                    or primer.internal_seq.endswith('G')
                                    or primer.internal_seq[1:-1].count('CC')
                                    > 1 or primer.internal_seq[1] == 'G')
            if not primer.oligovalid:
                invalidcount += 1
            #if (primer.internal_seq.startswith('G') or
            #    primer.internal_seq.endswith('G') or
            #    primer.internal_seq[1:-1].count('CC') > 1 or
            #    primer.internal_seq[1] == 'G'):
            #    primer.oligovalid = False
            #    invalidcount += 1
        logger.info("... %d primers failed (%.3fs)", invalidcount,
                    time.time() - time_start)

    # Filter primers on the basis of GC content at 3` end
    def filter_primers_gc_3prime(self):
        """ Loops over the primer pairs in the passed GenomeData object and,
            if either primer has more than 2 G+C in the last five nucleotides,
            sets the .gc3primevalid flag to False.
        """
        time_start = time.time()
        logger.info("Filtering %s primers on 3` GC content ...", self.name)
        invalidcount = 0
        for primer in self.primers.values():
            fseq, rseq = primer.forward_seq[-5:], primer.reverse_seq[-5:]
            if (fseq.count('C') + fseq.count('G') > 2) or \
                    (rseq.count('C') + fseq.count('G') > 2):
                primer.gc3primevalid = False
                invalidcount += 1
        logger.info("... %d primers failed (%.3fs)", invalidcount,
                    time.time() - time_start)

    # Concatenate multiple fragments of a genome to a single file
    def concatenate_sequences(self):
        """ Takes a GenomeData object and concatenates sequences with the
            spacer sequence NNNNNCATCCATTCATTAATTAATTAATGAATGAATGNNNNN (this
            contains start and stop codons in all frames, to cap individual
            sequences). We write this data out to a new file

            For filename convention, we just add '_concatenated' to the end
            of the sequence filestem, and use the '.fas' extension.
        """
        # Spacer contains start and stop codons in all six frames
        spacer = 'NNNNNCATCCATTCATTAATTAATTAATGAATGAATGNNNNN'
        time_start = time.time()
        logger.info("Concatenating sequences from %s ...", self.seqfilename)
        newseq = SeqRecord(Seq(spacer.join([s.seq.data for s in
                                            SeqIO.parse(open(self.seqfilename,
                                                             'rU'),
                                                        'fasta')])),
                           id=self.name + "_concatenated",
                           description="%s, concatenated with spacers" %
                           self.name)
        outfilename = ''.join([os.path.splitext(self.seqfilename)[0],
                               '_concatenated', '.fas'])
        SeqIO.write([newseq], open(outfilename, 'w'), 'fasta')
        logger.info("... wrote concatenated data to %s (%.3fs)",
                    outfilename, time.time() - time_start)
        return outfilename

    def __str__(self):
        """ Pretty string description of object contents
        """
        outstr = ['GenomeData object: %s' % self.name]
        outstr.append('Families: %s' % list(self.families))
        outstr.append('Sequence file: %s' % self.seqfilename)
        outstr.append('Feature file: %s' % self.ftfilename)
        outstr.append('Primers file: %s' % self.primerfilename)
        outstr.append('PrimerSearch file: %s' % self.primersearchfilename)
        outstr.append('Primers: %d' % len(self.primers))
        if len(self.primers):
            outstr.append('Primers overlapping CDS: %d' %
                          len([p for p in self.primers.values() if
                               p.cds_overlap]))
        return os.linesep.join(outstr) + os.linesep


###
# FUNCTIONS

# Parse command-line options
def parse_cmdline():
    """ Parse command line, accepting args obtained from sys.argv
    """
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-i", "--infile", dest="filename", action="store",
                      help="location of configuration file",
                      default=None)
    parser.add_option("-o", "--outdir", dest="outdir", action="store",
                      help="directory for output files",
                      default="differential_primer_results")
    parser.add_option("--numreturn", dest="numreturn", action="store",
                      help="number of primers to find",
                      default=20, type="int")
    parser.add_option("--hybridprobe", dest="hybridprobe", action="store_true",
                      help="generate internal oligo as a hybridisation probe",
                      default=False)
    parser.add_option("--filtergc3prime", dest="filtergc3prime",
                      action="store_true",
                      help="allow no more than two GC at the 3` " +
                           "end of primers",
                      default=False)
    parser.add_option("--single_product", dest="single_product",
                      action="store",
                      help="location of FASTA sequence file containing " +
                           "sequences from which a sequence-specific " +
                           "primer must amplify exactly one product.",
                      default=None)
    parser.add_option("--prodigal", dest="prodigal_exe", action="store",
                      help="location of Prodigal executable",
                      default="prodigal")
    parser.add_option("--eprimer3", dest="eprimer3_exe", action="store",
                      help="location of EMBOSS eprimer3 executable",
                      default="eprimer3")
    parser.add_option("--blast_exe", dest="blast_exe", action="store",
                      help="location of BLASTN/BLASTALL executable",
                      default="blastn")
    parser.add_option("--blastdb", dest="blastdb", action="store",
                      help="location of BLAST database",
                      default=None)
    parser.add_option("--useblast", dest="useblast", action="store_true",
                      help="use existing BLAST results",
                      default=False)
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
    parser.add_option("--clean", action="store_true", dest="clean",
                      help="clean up old output files before running",
                      default=False)
    parser.add_option("--cleanonly", action="store_true", dest="cleanonly",
                      help="clean up old output files and exit",
                      default=False)
    parser.add_option("-l", "--logfile", dest="logfile",
                      action="store", default=None,
                      help="script logfile location")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="report progress to log",
                      default=False)
    parser.add_option("--debug", action="store_true", dest="debug",
                      help="report extra progress to log for debugging",
                      default=False)
    parser.add_option("--keep_logs", action="store_true", dest="keep_logs",
                      help="store log files from each process",
                      default=False)
    parser.add_option("--log_dir", action="store", dest="log_dir",
                      help="store called process log files in this directory",
                      default=None)
    (optsparsed, argsparsed) = parser.parse_args()
    return (optsparsed, argsparsed, parser)


# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Create a list of GenomeData objects corresponding to config file entries
def create_gd_from_config(filename):
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
    time_start = time.time()
    logger.info("Creating list of genomes from config file %s ...", filename)
    gd_list = []                                   # Hold GenomeData objects
    # Ignore blank lines and comments...
    for line in [l.strip() for l in open(filename, 'rU')
                 if l.strip() and not l.startswith('#')]:
        # Split data and create new GenomeData object, adding it to the list
        data = [e.strip() for e in line.strip().split('\t') if e.strip()]
        name, family, sfile, ffile, pfile, psfile = tuple(data)
        gd_list.append(GenomeData(name, family, sfile, ffile, pfile, psfile))
        logger.info("... created GenomeData object for %s ...", name)
        logger.info(gd_list[-1])
    logger.info("... created %d GenomeData objects (%.3fs)",
                len(gd_list), time.time() - time_start)
    return gd_list


# Check whether each GenomeData object has multiple sequence and, if so,
# concatenate them sensibly, resetting feature and primer file locations to
# None
def check_single_sequence(gd_list):
    """ Loops over the GenomeData objects in the passed list and, where the
        sequence file contains multiple sequences, concatenates them into
        a single sequence using a spacer that facilitates gene-finding.  As
        this process changes feature and primer locations, the ftfilename and
        primerfilename attributes are reset to None, and these are
        recalculated later on in the script, where necessary.
    """
    time_start = time.time()
    logger.info("Checking for multiple sequences ...")
    for gd_obj in gd_list:
        # Verify that the sequence file contains a single sequence
        seqdata = [s for s in SeqIO.parse(open(gd_obj.seqfilename, 'rU'),
                                          'fasta')]
        if len(seqdata) != 1:
            logger.info("... %s describes multiple sequences ...",
                        gd_obj.seqfilename)
            gd_obj.seqfilename = gd_obj.concatenate_sequences()  # Concatenate
            logger.info("... clearing feature and primer file locations ...")
            (gd_obj.ftfilename, gd_obj.primerfilename,
             gd_obj.primersearchfilename) = \
                (None, None, None)
    logger.info("... checked %d GenomeData objects (%.3fs)",
                len(gd_list), time.time() - time_start)


# Check for each GenomeData object in a passed list, the existence of
# the feature file, and create one using Prodigal if it doesn't exist already
def check_ftfilenames(gd_list):
    """ Loop over the GenomeData objects in gdlist and, where no feature file
        is specified, add the GenomeData object to the list of
        packets to be processed in parallel by Prodigal using multiprocessing.
    """
    logger.info("Checking and predicting features for GenomeData files ...")

    # We split the GenomeData objects into those with, and without,
    # defined feature files, but we don't test the validity of the files
    # that were predefined, here.
    # We don't use the objects with features here, though
    #gds_with_ft = [gd_obj for gd_obj in gd_list if
    #               (gd_obj.ftfilename is not None and
    #                os.path.isfile(gd_obj.ftfilename))]
    gds_no_ft = [gd_obj for gd_obj in gd_list if
                 (gd_obj.ftfilename is None or
                  not os.path.isfile(gd_obj.ftfilename))]

    # Predict features for those GenomeData objects with no feature file
    logger.info("... %d GenomeData objects have no feature file ...",
                len(gds_no_ft))
    logger.info("... running %d Prodigal jobs to predict CDS ...",
                len(gds_no_ft))

    # Create a list of command-line tuples, for Prodigal
    # gene prediction applied to each GenomeData object in gds_no_ft.
    clines = []
    for gd_obj in gds_no_ft:
        gd_obj.ftfilename = os.path.splitext(gd_obj.seqfilename)[0] +\
            '.prodigalout'
        seqfilename = os.path.splitext(gd_obj.seqfilename)[0] + '.features'
        cline = "%s -a %s < %s > %s" % (options.prodigal_exe, seqfilename,
                                        gd_obj.seqfilename, gd_obj.ftfilename)
        clines.append(cline + log_output(gd_obj.name + ".prodigal"))
    logger.info("... Prodigal jobs to run:")
    logger.info("Running:\n" + "\n".join(clines))

    # Depending on the type of parallelisation required, these command-lines
    # are either run locally via multiprocessing, or passed out to SGE
    if not options.sge:
        multiprocessing_run(clines)
    else:
        sge_run(clines)


# Check whether GenomeData objects have a valid primer definition file
def check_primers(gd_list):
    """ Loop over GenomeData objects in the passed gdlist and, if they have
        a defined primerfilename attribute, attempt to parse it.  If this
        is successful, do nothing.  If it fails, set the primerfilename
        attribute to None.
    """
    logger.info("Checking ePrimer3 output files ...")
    for gd_obj in [g for g in gd_list if g.primerfilename]:
        try:
            Primer3.read(open(gd_obj.primerfilename, 'rU'))
            logger.info("... %s primer file %s OK ...",
                        gd_obj.name, gd_obj.primerfilename)
        except IOError:
            logger.info("... %s primer file %s not OK ...",
                        gd_obj.name, gd_obj.primerfilename)
            gd_obj.primerfilename = None


# Check for each GenomeData object in a passed list, the existence of
# the ePrimer3 file, and create one using ePrimer3 if it doesn't exist already
def predict_primers(gd_list, emboss_version):
    """ Loop over the GenomeData objects in gdlist and, where no primer file
        is specified, add the GenomeData object to the list of
        packets to be processed in parallel by Prodigal using multiprocessing.
    """
    logger.info("Checking and predicting primers for GenomeData files ...")
    # We need to split the GenomeData objects into those with, and without,
    # defined primer files, but we don't test the validity of these files
    # We don't use the gds
    #gds_with_primers = [g for g in gd_list if g.primerfilename is not None]
    gds_no_primers = [g for g in gd_list if g.primerfilename is None]
    # Predict primers for those GenomeData objects with no primer file
    logger.info("... %d GenomeData objects have no primer file ...",
                len(gds_no_primers))
    logger.info("... running %d ePrimer3 jobs to predict CDS ...",
                len(gds_no_primers))
    # Create command-lines to run ePrimer3
    clines = []
    for gd_obj in gds_no_primers:
        # Create ePrimer3 command-line.
        cline = Primer3Commandline(cmd=options.eprimer3_exe)
        cline.sequence = gd_obj.seqfilename
        cline.auto = True
        cline.osize = "%d" % options.osize            # Optimal primer size
        cline.minsize = "%d" % options.minsize        # Min primer size
        cline.maxsize = "%d" % options.maxsize        # Max primer size
        # Optimal primer Tm option dependent on EMBOSS version
        if float('.'.join(emboss_version.split('.')[:2])) >= 6.5:
            cline.opttm = "%d" % options.otm              # Optimal primer Tm
        else:
            cline.otm = "%d" % options.otm
        cline.mintm = "%d" % options.mintm            # Min primer Tm
        cline.maxtm = "%d" % options.maxtm            # Max primer Tm
        cline.ogcpercent = "%d" % options.ogcpercent  # Optimal primer %GC
        cline.mingc = "%d" % options.mingc            # Min primer %GC
        cline.maxgc = "%d" % options.maxgc            # Max primer %GC
        cline.psizeopt = "%d" % options.psizeopt      # Optimal product size
        # Longest polyX run in primer
        cline.maxpolyx = "%d" % options.maxpolyx
        # Allowed product sizes
        cline.prange = "%d-%d" % (options.psizemin, options.psizemax)
        # Number of primers to predict
        cline.numreturn = "%d" % options.numreturn
        cline.hybridprobe = options.hybridprobe  # Predict internal oligo?
        # Internal oligo parameters;
        cline.osizeopt = "%d" % options.oligoosize
        # We use EMBOSS v6 parameter names, here.
        cline.ominsize = "%d" % options.oligominsize
        cline.omaxsize = "%d" % options.oligomaxsize
        cline.otmopt = "%d" % options.oligootm
        cline.otmmin = "%d" % options.oligomintm
        cline.otmmax = "%d" % options.oligomaxtm
        cline.ogcopt = "%d" % options.oligoogcpercent
        cline.ogcmin = "%d" % options.oligomingc
        cline.ogcmax = "%d" % options.oligomaxgc
        cline.opolyxmax = "%d" % options.oligomaxpolyx
        cline.outfile = os.path.splitext(gd_obj.seqfilename)[0] + '.eprimer3'
        gd_obj.primerfilename = cline.outfile
        clines.append(str(cline) + log_output(gd_obj.name + ".eprimer3"))
    logger.info("... ePrimer3 jobs to run:")
    logger.info("Running:\n" + '\n'.join(clines))
    # Parallelise jobs
    if not options.sge:
        multiprocessing_run(clines)
    else:
        sge_run(clines)


# Load primers from ePrimer3 files into each GenomeData object
def load_primers(gd_list):
    """ Load primer data from an ePrimer3 output file into a dictionary of
        Bio.Emboss.Primer3.Primer objects (keyed by primer name) in a
        GenomeData object, for each such object in the passed list.
        Each primer object is given a new ad hoc attribute 'cds_overlap' which
        takes a Boolean, indicating whether the primer is found wholly within
        a CDS defined in the GenomeData object's feature file; this status
        is determined using an interval tree approach.
    """
    logger.info("Loading primers, %sfiltering on CDS overlap",
                'not ' if options.nocds else '')
    # Load in the primers, assigning False to a new, ad hoc attribute called
    # cds_overlap in each
    for gd_obj in gd_list:
        logger.info("... loading primers into %s from %s ...",
                    gd_obj.name, gd_obj.primerfilename)
        try:
            os.path.isfile(gd_obj.primerfilename)
        except TypeError:
            raise IOError("Primer file %s does not exist." %
                          gd_obj.primerfilename)
        primers = Primer3.read(open(gd_obj.primerfilename, 'rU')).primers
        # Add primer pairs to the gd.primers dictionary
        primercount = 0
        for primer in primers:
            primercount += 1
            primer.cds_overlap = False           # default state
            primer.name = "%s_primer_%04d" % (gd_obj.name, primercount)
            primer.amplifies_organism = set()    # Organisms amplified
            primer.amplifies_family = set()      # Organism families amplified
            primer.gc3primevalid = True          # Passes GC 3` test
            primer.oligovalid = True             # Oligo passes filter
            primer.blastpass = True              # Primers pass BLAST screen
            gd_obj.primers.setdefault(primer.name, primer)
            primer.amplicon = \
                gd_obj.sequence[primer.forward_start - 1:
                                primer.reverse_start - 1 +
                                primer.reverse_length]
            primer.amplicon.description = primer.name
        logger.info("... loaded %d primers into %s ...",
                    len(gd_obj.primers), gd_obj.name)
        # Now that the primers are in the GenomeData object, we can filter
        # them on location, if necessary
        if not options.nocds:
            gd_obj.filter_primers(options.psizemin)
        # We also filter primers on the basis of GC presence at the 3` end
        if options.filtergc3prime:
            gd_obj.filter_primers_gc_3prime()
        # Filter primers on the basis of internal oligo characteristics
        if options.hybridprobe:
            gd_obj.filter_primers_oligo()


# Screen passed GenomeData primers against BLAST database
def blast_screen(gd_list):
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
    build_blast_input(gd_list)
    run_blast(gd_list)
    parse_blast(gd_list)


# Write BLAST input files for each GenomeData object
def build_blast_input(gd_list):
    """ Loops over each GenomeData object in the list, and writes forward
        and reverse primer sequences out in FASTA format to a file with
        filename derived from the GenomeData object name.
    """
    time_start = time.time()
    logger.info("Writing files for BLAST input ...")
    for gd_obj in gd_list:
        gd_obj.blastinfilename =\
            os.path.join(os.path.split(gd_obj.seqfilename)[0],
                         "%s_BLAST_input.fas" % gd_obj.name)
        seqrecords = []
        for name, primer in gd_obj.primers.items():
            seqrecords.append(SeqRecord(Seq(primer.forward_seq),
                                        id=name + '_forward'))
            seqrecords.append(SeqRecord(Seq(primer.reverse_seq),
                                        id=name + '_reverse'))
        logger.info("... writing %s ...", gd_obj.blastinfilename)
        SeqIO.write(seqrecords,
                    open(gd_obj.blastinfilename, 'w'),
                    'fasta')
    logger.info("... done (%.3fs)", time.time() - time_start)


# Run BLAST screen for each GenomeData object
def run_blast(gd_list):
    """ Loop over the GenomeData objects in the passed list, and run a
        suitable BLASTN query with the primer sequences, writing to a file
        with name derived from the GenomeData object, in XML format.
    """
    logger.info("Compiling BLASTN command-lines ...")
    clines = []
    for gd_obj in gd_list:
        gd_obj.blastoutfilename =\
            os.path.join(os.path.split(gd_obj.seqfilename)[0],
                         "%s_BLAST_output.xml" % gd_obj.name)
        cline = NcbiblastnCommandline(query=gd_obj.blastinfilename,
                                      db=options.blastdb,
                                      task='blastn',  # default: MEGABLAST
                                      out=gd_obj.blastoutfilename,
                                      num_alignments=1,
                                      num_descriptions=1,
                                      outfmt=5,
                                      perc_identity=90,
                                      ungapped=True)
        clines.append(str(cline) + log_output(gd_obj.name + ".blastn"))
    logger.info("... BLASTN+ jobs to run:")
    logger.info("Running:\n" + '\n'.join(clines))
    if not options.sge:
        multiprocessing_run(clines)
    else:
        sge_run(clines)


# Parse BLAST output for each GenomeData object
def parse_blast(gd_list):
    """ Loop over the GenomeData objects in the passed list, and parse the
        BLAST XML output indicated in the .blastoutfilename attribute.
        For each query that makes a suitable match, mark the appropriate
        primer's .blastpass attribute as False
    """
    time_start = time.time()
    logger.info("Parsing BLASTN output with multiprocessing ...")
    # Here I'm cheating a bit and using multiprocessing directly so that
    # we can speed up the parsing process a bit
    pool = multiprocessing.Pool(processes=options.cpus)
    pool_results = [pool.apply_async(process_blastxml,
                                     (g.blastoutfilename, g.name))
                    for g in gd_list]
    pool.close()
    pool.join()
    # Process the results returned from the BLAST searches.  Create a
    # dictionary of GenomeData objects, keyed by name, and loop over the
    # result sets, setting .blastpass attributes for the primers as we go
    gddict = {}
    for gd_obj in gd_list:
        gddict.setdefault(gd_obj.name, gd_obj)
    failcount = 0
    for result in [r.get() for r in pool_results]:
        for name in result:
            gd_obj = gddict[name.split('_primer_')[0]]
            gd_obj.primers[name].blastpass = False
            failcount += 1
    logger.info("... %d primers failed BLAST screen ...", failcount)
    logger.info("... multiprocessing BLAST parsing complete (%.3fs)",
                time.time() - time_start)


# BLAST XML parsing function for multiprocessing
def process_blastxml(filename, name):
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
    time_start = time.time()
    logger.info("[process name: %s] Parsing BLAST XML ...", name)
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
                    matching_primers.add('_'.join(
                        record.query.split('_')[:-1]))
        logger.info("[process name: %s] Parsed %d records",
                    name, recordcount)
    except IOError:
        logger.info("[process name: %s] Error reading BLAST XML file", name)
    logger.info("[process name: %s] Time spent in process: (%.3fs)",
                name, time.time() - time_start)
    # Return the list of matching primers
    return matching_primers


# A function for parsing features from Prodigal output
def parse_prodigal_features(filename):
    """ Parse Prodigal 'GenBank' output.
        We try to emulate SeqIO.read() SeqRecord output as much as possible,
        but the information provided by Prodigal is limited to feature type
        and location, on a single line.
        Amended: Newer versions of Prodigal write closer match to GenBank
        format, and thus if the first line matches "DEFINITION" we use SeqIO.
        RE-amended: Latest version of Prodigal is still not good enough for
        SeqIO, so a new function is created to parse line-by-line.
    """
    record = SeqRecord(None)       # record gets a dummy sequence
    # Open filehandle and parse contents
    handle = open(filename, 'rU')
    # init feature list from file parsing
    record.features = seqrecord_parse(handle)
    return record


# Parse record features from the lines of prodigal or genbank format file
def seqrecord_parse(filehandle):
    """ Parses the head lines of CDS features from a Prodigal or GenBank
        file.
        This is still necessary, as Prodigal's GenBank output is not
        SeqIO.read() friendly.
    """
    features = []
    for line in filehandle:
        if re.search("CDS", line):
            data = [e.strip() for e in line.split()]
            feature = gb_string_to_feature(data[-1])
            feature.type = data[0]
            features.append(feature)
    return features


# Parse record features from sequence file, using SeqIO
def seqrecord_parse_seqio(filehandle, seqformat):
    """ NOTE: Latest version of prodigal output is *closer* to GenBank format
        but not close enough for SeqIO to find the genome.features
        Thus: this function NOT USED (until potential update to prodigal
        or SeqIO).
    """
    features = []
    seqrecord = list(SeqIO.parse(filehandle, seqformat))
    for record in seqrecord:
        logger.debug("record seq: [%s]...", record.seq[0:12])
        features.append(record.features)
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
def primersearch(gd_list):
    """ Loop over the GenomeData objects in the passed list, and construct
        command lines for an all-against-all PrimerSearch run.
        Output files are of the format
        <query name>_vs_<target name>.primersearch
        Where <query name> and <target name> are the gd.name attributes of
        the source and target GenomeData objects, respectively.
        The output file goes in the same location as the source sequence
        file.
    """
    logger.info("Constructing all-against-all PrimerSearch runs " +
                "for %d objects ...", len(gd_list))
    # Create list of command-lines
    clines = []
    for query_gd in gd_list:
        query_gd.primersearch_output = []
        for target_gd in gd_list:
            if query_gd != target_gd:
                # Location of PrimerSearch output
                outdir = os.path.split(query_gd.seqfilename)[0]
                outfilename = os.path.join(outdir, "%s_vs_%s.primersearch" %
                                           (query_gd.name, target_gd.name))
                query_gd.primersearch_output.append(outfilename)
                # Create command-line
                cline = PrimerSearchCommandline()
                cline.auto = True
                cline.seqall = target_gd.seqfilename
                cline.infile = query_gd.primersearchfilename
                cline.outfile = outfilename
                cline.mismatchpercent = options.mismatchpercent
                clines.append(str(cline) +
                              log_output(os.path.basename(outfilename)))
    logger.info("... PrimerSearch jobs to run: ...")
    logger.info("Running:\n" + '\n'.join(clines))
    # Parallelise jobs
    if not options.sge:
        multiprocessing_run(clines)
    else:
        sge_run(clines)


# Load in existing PrimerSearch output
def load_existing_primersearch_results(gd_list):
    """ Associates PrimerSearch output files with each GenomeData object
        and returns a list of (name, filename) tuples for all GenomeData
        objects
    """
    time_start = time.time()
    logger.info("Locating existing PrimerSearch input files ...")
    primersearch_results = []
    for gd_obj in gd_list:
        gd_obj.primersearch_output = []
        filedir = os.path.split(gd_obj.seqfilename)[0]
        primersearch_files = [f for f in os.listdir(filedir) if
                              os.path.splitext(f)[-1] == '.primersearch' and
                              f.startswith(gd_obj.name)]
        for filename in primersearch_files:
            logger.info("... found %s for %s ...", filename, gd_obj.name)
            gd_obj.primersearch_output.append(os.path.join(filedir,
                                                           filename))
    logger.info("... found %d PrimerSearch input files (%.3fs)",
                len(primersearch_results), time.time() - time_start)


# Run primersearch to find whether and where the predicted primers amplify
# our negative target (the one we expect exactly one match to)
def find_negative_target_products(gd_list):
    """ We run primersearch using the predicted primers as queries, with
        options.single_product as the target sequence.  We exploit
        multiprocessing, and use the prescribed number of
        CPUs.  Happily, primersearch accepts multiple sequence FASTA files.
    """
    logger.info("Constructing negative control PrimerSearch runs " +
                "for %d objects ...", len(gd_list))
    # Create list of command-lines
    clines = []
    for query_gd in gd_list:
        query_gd.primersearch_output = []
        outdir = os.path.split(query_gd.seqfilename)[0]
        outfilename = os.path.join(outdir, "%s_negative_control.primersearch" %
                                   query_gd.name)
        query_gd.primersearch_output.append(outfilename)
        # Create command-line
        cline = PrimerSearchCommandline()
        cline.auto = True
        cline.seqall = options.single_product
        cline.infile = query_gd.primersearchfilename
        cline.outfile = outfilename
        cline.mismatchpercent = options.mismatchpercent
        clines.append(str(cline) + log_output(os.path.basename(outfilename)))
    logger.info("... PrimerSearch jobs to run: ...")
    logger.info("Running:\n" + '\n'.join(clines))
    # Parallelise jobs and run
    if not options.sge:
        multiprocessing_run(clines)
    else:
        sge_run(clines)


# Classify the primers in a list of GenomeData objects according to the
# other sequences that they amplify
def classify_primers(gd_list):
    """ Takes a list of GenomeData objects and loops over the primersearch
        results, loading in the primersearch results and applying them to the
        associated query GenomeData object.
        If a primer is reported, by PrimerSearch, to amplify a region of the
        target genome, two changes are made to the corresponding Primer
        object in the amplifies_object and amplifies_family ad hoc attributes,
        with the target name and family, respectively, being added to those
        sets.
    """
    time_start = time.time()
    logger.info("Classifying primers by PrimerSearch results ...")
    # Convenience dictionary, keying each GenomeData object by name
    # We need to load this fully before checking the PrimerSearch output
    # below.
    gddict = {}
    for gd_obj in gd_list:
        gddict.setdefault(gd_obj.name, gd_obj)
    # Parse the PrimerSearch output, updating the primer contents of the
    # appropriate GenomeData object, for each set of results
    for gd_obj in gd_list:
        logger.info("... GenomeData for %s ...", gd_obj.name)
        for filename in gd_obj.primersearch_output:
            logger.info("... processing %s ...", filename)
            # Identify the target organism
            targetname = \
                os.path.splitext(os.path.split(
                    filename)[-1])[0].split('_vs_')[-1]
            # Only classify amplimers to sequences in the gdlist dataset
            # This avoids problems with recording counts of matches to
            # sequences that we're not considering, artifically lowering the
            # specificity counts.
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
                        gd_obj.primers[pname].amplifies_organism.add(
                            targetname)
                        for family in gddict[targetname].families:
                            gd_obj.primers[pname].amplifies_family.add(family)
            # Consider the negative control primersearch output
            elif 'negative_control' in filename:
                # Load PrimerSearch data
                psdata = PrimerSearch.read(open(filename, 'rU'))
                # We loop over each primer, and find the number of amplimers.
                # We note the number of amplimers as an attribute of the primer
                for pname, pdata in psdata.amplifiers.items():
                    gd_obj.primers[pname].negative_control_amplimers =\
                        len(pdata)
                    logger.info("Found %d amplimers in negative control",
                                len(pdata))
        logger.info("... processed %d Primersearch results for %s ...",
                    len(gd_obj.primersearch_output), gd_obj.name)
    logger.info("... processed PrimerSearch results (%.3fs)",
                time.time() - time_start)


# Write analysis data to files
def write_report(gd_list, blastfilter):
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
    time_start = time.time()
    logger.info("Creating summary output ...")
    # First we need to generate a dictionary of GenomeData object names, keyed
    # by family
    families = defaultdict(set)
    for gd_obj in gd_list:
        for family in gd_obj.families:
            families[family].add(gd_obj.name)
    # Rectify nocds flag
    cds_overlap = not options.nocds
    # Check whether output directory exists and, if not, create it
    if not os.path.isdir(options.outdir):
        os.mkdir(options.outdir)
    # Open output file, and write header
    outfh = open(os.path.join(options.outdir,
                              'differential_primer_results.tab'), 'w')
    outfh.write(os.linesep.join([
        "# Summary information table",
        "# Generated by find_differential_primers",
        "# Columns in the table:",
        "# 1) Query organism ID",
        "# 2) Query organism families",
        "# 3) Count of organism-unique primers",
        "# 4) Count of universal primers",
        "# 5) Query sequence filename",
        "# 6) Query feature filename",
        "# 7) Query ePrimer3 primers filename"]) + '\n')
    # Write data for each GenomeData object
    other_org_count = len(gd_list) - 1  # Amplifications for 'universal' set
    # We store 'universal' primers in their own list, and family-specific
    # primers in a dicitonary, keyed by family
    all_universal_primers = []
    family_specific_primers = defaultdict(list)
    # Loop over each GenomeData object and populate family-specific and
    # universal primer collections, as well as organism-specific and
    # summary information
    for gd_obj in gd_list:
        logger.info('\n'.join([
            "... writing data for %s ..." % gd_obj.name,
            "... cds_overlap: %s ..." % cds_overlap,
            "... gc3primevalid: %s ..." % options.filtergc3prime,
            "... oligovalid: %s ..." % options.hybridprobe,
            "... blastpass: %s ..." % blastfilter,
            "... single_product %s ..." % (options.single_product is
                                           not None),
            "... retrieving primer pairs ...",
            "... finding strain-specific primers for %s ..." % gd_obj.name
        ]))
        unique_primers = gd_obj.get_unique_primers(cds_overlap, blastfilter)
        logger.info("... finding family-specific primers for %s ...",
                    gd_obj.name)
        family_unique_primers = {}
        for family in gd_obj.families:
            logger.info("Checking family: %s" % family)
            logger.info("families[%s]: %s" % (family, families[family]))
            family_unique_primers[family] = \
                gd_obj.get_family_unique_primers(families[family], cds_overlap,
                                                 blastfilter)
            family_specific_primers[family] += family_unique_primers[family]
            logger.info("family_unique_primers[%s]: %d" % 
                        (family, len(family_unique_primers[family])))
            logger.info("family_specific_primers[%s]: %d" % 
                        (family, len(family_specific_primers[family])))
        logger.info("... finding universal primers for %s ...", gd_obj.name)
        universal_primers = \
            gd_obj.get_primers_amplify_count(other_org_count, cds_overlap,
                                             blastfilter)
        all_universal_primers.extend(universal_primers)
        # Write summary data to file
        outfh.write('\t'.join([gd_obj.name, ','.join(gd_obj.families),
                               str(len(unique_primers)),
                               str(len(universal_primers)),
                               str(gd_obj.seqfilename),
                               str(gd_obj.ftfilename),
                               str(gd_obj.primerfilename)]) + '\n')
        # Write organism-specific primers to file
        write_eprimer3(unique_primers,
                       os.path.join(options.outdir,
                                    "%s_specific_primers.eprimer3" %
                                    gd_obj.name), gd_obj.seqfilename)
        # Write organism-specific amplicons to file
        SeqIO.write([p.amplicon for p in unique_primers],
                    os.path.join(options.outdir,
                                 "%s_specific_amplicons.fas" % gd_obj.name),
                    'fasta')
    outfh.close()

    # Write universal primers to file
    write_eprimer3(universal_primers,
                   os.path.join(options.outdir, "universal_primers.eprimer3"),
                   '', append=True)

    # Write organism-specific amplicons to file
    SeqIO.write([p.amplicon for p in universal_primers],
                open(os.path.join(options.outdir,
                                  "universal_amplicons.fas"), 'w'),
                'fasta')

    # Write family-specific primers to files
    outfh = open(os.path.join(options.outdir,
                              'differential_primer_results-families.tab'),
                 'w')
    outfh.write(os.linesep.join([
        "# Summary information table",
        "# Generated by find_differential_primers",
        "# Columns in the table:",
        "# 1) Family",
        "# 2) Count of family-specific primers",
        "# 3) Family-specific primer file",
        "# 4) Family-specific amplicon file"]) + '\n')
    for family, primers in family_specific_primers.items():
        outstr = [family, str(len(primers))]
        fname = os.path.join(options.outdir,
                             "%s_family-specific_primers.eprimer3" %
                             family)
        write_eprimer3(primers, fname, '')
        outstr.append(fname)
        # Write family-specific amplicons to file
        fname = os.path.join(options.outdir,
                             "%s_family-specific_amplicons.fas" %
                             family)
        SeqIO.write([p.amplicon for p in primers], open(fname, 'w'), 'fasta')
        outstr.append(fname)
        outfh.write('\t'.join(outstr) + '\n')
    # Being tidy...
    outfh.close()
    logger.info("... data written (%.3fs)", time.time() - time_start)


# Write ePrimer3 format primer file
def write_eprimer3(primers, filename, sourcefilename, append=False):
    """ Write the passed primers to the passed file, in ePrimer3-
        compatible format.
    """
    logger.info("Writing %d primer pairs to %s ...", len(primers), filename)
    # Open file
    filemode = 'a' if append else 'w'     # Do we append or write anew?
    outfh = open(filename, filemode)
    # Write header
    outfh.write(os.linesep.join([
        "# EPRIMER3 PRIMERS %s " % filename,
        "#                      Start  Len   Tm     GC%   Sequence",
        os.linesep]) + '\n')
    primercount = 0
    for primer in primers:
        primercount += 1
        outfh.write("# %s %s\n" % (primer.name, sourcefilename))
        outfh.write("%-4d PRODUCT SIZE: %d\n" % (primercount, primer.size))
        outfh.write("     FORWARD PRIMER  %-9d  %-3d  %.02f  %.02f  %s\n" %
                    (primer.forward_start, primer.forward_length,
                     primer.forward_tm, primer.forward_gc,
                     primer.forward_seq))
        outfh.write("     REVERSE PRIMER  %-9d  %-3d  %.02f  %.02f  %s\n" %
                    (primer.reverse_start, primer.reverse_length,
                     primer.reverse_tm, primer.reverse_gc,
                     primer.reverse_seq))
        if hasattr(primer, 'internal_start'):
            outfh.write("     INTERNAL OLIGO  %-9d  %-3d  %.02f  %.02f  %s\n" %
                        (primer.internal_start, primer.internal_length,
                         primer.internal_tm, primer.internal_gc,
                         primer.internal_seq))
        outfh.write(os.linesep * 3)
    # Be tidy
    outfh.close()


# Run the passed list of command-lines using a multiprocessing.Pool
def multiprocessing_run(clines):
    """ We create a multiprocessing Pool to handle command-lines  We
        pass the (unique) GenomeData object name, and the location of the
        sequence file.  The called function returns the GenomeData name and the
        corresponding location of the generated feature file.  The GenomeData
        objects are stored in a temporary dictionary, keyed by gd.name, to
        allow association of the results of the asynchronous pool jobs with the
        correct GenomeData object
    """
    time_start = time.time()
    logger.info("Running %d jobs with multiprocessing ...",
                len(clines))
    pool = multiprocessing.Pool(processes=options.cpus)  # create process pool
    completed = []
    if options.verbose:
        callback_fn = multiprocessing_callback
    else:
        callback_fn = completed.append
    for cline in clines:
        pool.apply_async(subprocess.call,
                         (str(cline), ),
                         {'stderr': subprocess.PIPE,
                          'shell': sys.platform != "win32"},
                         callback=callback_fn)
    pool.close()      # Run jobs
    pool.join()
    logger.info("Completed:\n" + '\n'.join([str(e) for e in completed]))
    logger.info("... all multiprocessing jobs ended (%.3fs)",
                time.time() - time_start)


# Add a multiprocessing callback function here
def multiprocessing_callback(val):
    """ A verbose callback function for multiprocessing runs.  It uses the
        return value to indicate run completion or failure.  Failure is
        indicated by a nonzero return from the multiprocessing call.
    """
    if 0 == val:
        logger.info("... multiprocessing run completed (status: %s) ...", val)
    else:
        logger.error("... problem with multiprocessing run (status: %s) ...",
                     val)


# Clean output for each GenomeData object in the passed list
def clean_output(gd_list):
    """ Remove .eprimer3, .primers, .prodigalout, and .primersearch files
        from the same directory as the sequence file for each passed
        PrimerSearch object
    """
    time_start = time.time()
    logger.info("Cleaning up output files for GenomeData objects ...")
    # Loop over each GenomeData object, and remove each output file
    for gd_obj in gd_list:
        seqdir = os.path.split(gd_obj.seqfilename)[0]
        for filename in [f for f in os.listdir(seqdir)
                         if os.path.splitext(f)[-1] in
                         ['.eprimer3', 'primers', '.prodigalout',
                          '.primersearch', '.xml']]:
            abspath = os.path.join(seqdir, filename)
            logger.info("... deleting %s ...", abspath)
            os.remove(abspath)     # You can never go back after this point
    logger.info("... done (%.3fs)", time.time() - time_start)


# construct str to concat on end of cline if option.keep_logs is set
def log_output(filename):
    """ predefine file extension and stream to print to.
        if log_dir exists, join it to filename
        else output to base filename.
    """
    log_extension = ".log"
    log_out_handle = " 2> "
    if options.keep_logs and options.log_dir:
        return log_out_handle + os.path.join(options.log_dir, filename) +\
            log_extension
    elif options.keep_logs:
        return log_out_handle + filename + log_extension
    else:
        return ""


# run list of command-line jobs with SGE
def sge_run(*args):
    """ Function intended to compile a passed list of command lines, and
        run them on SGE.
    """
    raise NotImplementedError


###
# SCRIPT
if __name__ == '__main__':
    # Parse cmd-line
    options, arguments, optparser = parse_cmdline()

    # Set up logging, and modify loglevel according to whether we need
    # verbosity or not
    # err_handler points to sys.stderr
    # err_handler_file points to a logfile, if named
    logger = logging.getLogger('find_differential_primers.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.logfile is not None:
        try:
            logstream = open(options.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except IOError:
            logger.error("Could not open %s for logging",
                         options.logfile)
            sys.exit(1)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    logger.info('# find_differential_primers.py logfile')
    logger.info('# Run: %s', time.asctime())

    # Report arguments, if verbose
    logger.info(options)
    logger.info(arguments)

    # Create our GenomeData objects.  If there is no configuration file
    # specified, raise an error and exit.  Otherwise we end up with a list
    # of GenomeData objects that are populated only with the data from the
    # config file
    if options.filename is None:
        optparser.print_help()
        raise IOError("No configuration file specified")
    gdlist = create_gd_from_config(options.filename)

    # If the user wants to clean the directory before starting, do so
    if options.clean or options.cleanonly:
        clean_output(gdlist)
    if options.cleanonly:
        sys.exit(0)

    # It is possible that the sequence file for a GenomeData object might
    # be a multi-sequence file describing scaffolds or contigs.  We create a
    # concatenated sequence to facilitate further analyses, if this is the
    # case.  Where a sequence needs to be concatenated, this will affect the
    # placement of features and/or primers, so any specified files are
    # reset to None
    check_single_sequence(gdlist)

    # What EMBOSS version is available? This is important as the ePrimer3
    # command-line changes in v6.6.0, which is awkward for the Biopython
    # interface.
    embossversion = \
        subprocess.check_output("embossversion",
                                stderr=subprocess.PIPE,
                                shell=sys.platform != "win32").strip()
    logger.info("EMBOSS version reported as: %s", embossversion)

    # We need to check the existence of a prescribed feature file and, if
    # there is not one, create it.  We don't bother if the --nocds flag is set.
    if not (options.nocds or options.noprodigal):
        logger.info("--nocds option not set: " +
                    "Checking existence of features...")
        check_ftfilenames(gdlist)
    elif options.nocds:
        logger.warning("--nocds option set: Not checking or " +
                       "creating feature files")
    else:
        logger.warning("--noprodigal option set: Not predicting new CDS")

    # We need to check for the existence of primer sequences for the organism
    # and, if they do not exist, create them using ePrimer3.  If the
    # --noprimer3 flag is set, we do not create new primers, but even if the
    # --noprimersearch flag is set, we still need to check whether the
    # primer files are valid
    if not options.noprimer3:
        logger.info("--noprimer3 flag not set: Predicting new primers")
        check_primers(gdlist)
        predict_primers(gdlist, embossversion)
    else:
        logger.warning("--noprimer3 flag set: Not predicting new primers")

    # With a set of primers designed for the organism, we can load them into
    # the GenomeData object, filtering for those present only in the CDS,
    # if required.  This step is necessary, whether or not a new ePrimer3
    # prediction is made.  We also filter on GC content at the primer 3' end,
    # if required.
    logger.info("Loading primers...")
    load_primers(gdlist)

    # At this point, we can check our primers against a prescribed BLAST
    # database.  How we filter these depends on the user's preference.
    # We screen against BLAST here so that we can flag an attribute on
    # each primer to say whether or not it passed the BLAST screen.
    if options.blastdb and not options.useblast:
        logger.info("--blastdb options set: BLAST screening primers...")
        blast_screen(gdlist)
    elif options.useblast:
        logger.warning("--useblast option set: " +
                       "using existing BLAST results...")
    else:
        logger.warning("No BLAST options set, not BLAST screening primers...")
    # Having a set of (potentially CDS-filtered) primers for each organism,
    # we then scan these primers against each of the other organisms in the
    # set, using the EMBOSS PrimerSearch package
    # (http://embossgui.sourceforge.net/demo/manual/primersearch.html)
    # Now we have all the data we need to run PrimerSearch in an all-vs-all
    # manner, so make a cup of tea, put your feet up, and do the comparisons
    # with EMBOSS PrimerSearch
    # (http://embossgui.sourceforge.net/demo/manual/primersearch.html)
    if options.noprimersearch:
        logger.warning("--noprimersearch flag set: Not running PrimerSearch")
        # Load the appropriate primersearch output files for each
        # GenomeData object
        load_existing_primersearch_results(gdlist)
    else:
        logger.info("--noprimersearch flag not set: Running PrimerSearch")
        # We write input for PrimerSearch ignoring all the filters; this lets
        # us turn off PrimerSearch and rerun the analysis with alternative
        # filter settings
        for gd in gdlist:
            gd.write_primers()
        # Run PrimerSearch
        primersearch(gdlist)
    # If the --single_product option is specified, we load in the sequence
    # file to which the passed argument refers, and filter the primer
    # sequences on the basis of how many amplification products are produced
    # from these sequences.  We expect exactly one amplification product per
    # primer set, if it's not degenerate on the target sequence
    # (note that this filter is meaningless for family-specific primers)
    if options.single_product:
        find_negative_target_products(gdlist)
        logger.info("--blastdb options set: BLAST screening primers...")
        blast_screen(gdlist)

    # Now we classify the primer sets according to which sequences they amplify
    if not options.noclassify:
        logger.info("Classifying primers and writing output files ...")
        # Classify the primers in each GenomeData object according to
        # the organisms and families that they amplify, using the
        # PrimerSearch results.
        classify_primers(gdlist)
        # All the data has been loaded and processed, so we can now create our
        # plaintext summary report of the number of unique, family-unique and
        # universal primers in each of the organisms
        write_report(gdlist, (options.blastdb is not None or options.useblast))
