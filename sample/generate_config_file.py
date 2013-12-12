# generate_config_file.py
#
# Generate a config file for the primer prediction script, from the
# downloaded published E. coli sequences (in sequences), and the O104 sequences
# (in O104_sequences)
#
# All GenBank
#
# These sequences come from links at
# http://pathogenomics.bham.ac.uk/blog/2011/06/ehec-genome-assembly/:
#
# Nick Loman's assembly of the released BGI reads (TY2482):
# http://static.xbase.ac.uk/files/results/nick/TY2482/TY2482.fasta.txt
#
# ERA7's copies of the released BGI assemblies (TY2482):
# https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/blob/master/strains/TY2482/seqProject/BGI/assemblies/BGI/Escherichia_coli_TY-2482.chromosome.20110616.fa.gz
# https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/tree/master/strains/TY2482/seqProject/BGI/assemblies/BGI
# This contains nine strains, each under the ehec...a0bcca3/strains/<STRAIN>/
# seqProject/<X>/assemblies/<X>/ directory
# I moved each set of contigs up to O104_sequences
# The latest BGI assemblies were concatenated with
# cat Escherichia_coli_TY-2482.*.fa > TY-2482-BGI.fas
#
#
# NCBI's assembly (LB226692)
# http://www.ncbi.nlm.nih.gov/nuccore/AFOB00000000
# http://static.xbase.ac.uk/files/results/nick/TY2482/German_vs_TWEC.gbk
# This file was in .gbk format, and needed to be converted to FASTA prior
# to analysis.  Since the 'chromosome' was entirely reconstructed with Ns, I
# had to split these before generating the FASTA file.
#
# We need to stitch these files with the six-frame separator
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

### GLOBALS

conf_file_header = """
# Configuration file for a script to find differential primers in the Dickeya
# genomes.
#
# This file defines the following data, in tab-separated format:
# Column 1: Organism abbreviation
# Column 2: Family/group for the organism, comma-separated for multiple groups
# Column 3: Location of sequence data in FASTA format
# Column 4: Location of GenBank file describing features ('-' if none)
# Column 5: Location of ePrimer3 primer definitions ('-' if none)
# Column 6: Location of PrimerSearch input format primer definitions ('-' if none)
#
# BLANK COLUMNS ARE IGNORED! USE '-' IF THERE IS NO DATA.
#
# Copyright Leighton Pritchard 2011
#
# E. coli O104 primer prediction 20110627

"""

### DEFINE LIST OF O104 SEQUENCES AND CLASS

o104seqs = [['GOS1', 'O104', 'O104_sequences/GOS1_contigs_stitched.fasta', '-', '-', '-'],
            ['GOS2', 'O104', 'O104_sequences/GOS2_contigs_stitched.fasta', '-', '-', '-'],
            ['H112180280', 'O104', 'O104_sequences/Sample280_contigs_stitched.fasta', '-', '-', '-'],
            ['H112240282', 'O104', 'O104_sequences/Sample282_contigs_stitched.fasta', '-', '-', '-'],
            ['H112240283', 'O104', 'O104_sequences/Sample283_contigs_stitched.fasta', '-', '-', '-'],
            ['H112240540', 'O104', 'O104_sequences/Sample540_contigs_stitched.fasta', '-', '-', '-'],
            ['H112240541', 'O104', 'O104_sequences/Sample541_contigs_stitched.fasta', '-', '-', '-'],
            ['LB226692-NCBI', 'O104', 'O104_sequences/LB226692_contigs_stitched.fasta', '-', '-', '-'],
            ['LB226692-Muenster', 'O104', 'O104_sequences/wgs.AFOB.1_stitched.fasta', '-', '-', '-'],
            ['TY2482-BGI1', 'O104', 'O104_sequences/TY-2482-BGI_stitched.fasta', '-', '-', '-'],
            ['TY2482-LOMAN', 'O104', 'O104_sequences/TY2482.fasta_stitched.fasta', '-', '-', '-']
            ]

### AUTOMATICALLY GENERATE LIST OF NON-O104 SEQUENCES

import os
from Bio import SeqIO

non_o104seqs = []
for filename in [f for f in os.listdir('sequences') if 'fixed.fasta' in f]:
    seqid = SeqIO.read(os.path.join('sequences', filename), \
                           'fasta').description.split(' ', 1)[-1]
    non_o104seqs.append(['_'.join(seqid.split()), 'not_O104', os.path.join('sequences', filename),
                         '-', '-', '-'])

# Add the feature locations and eprimer3 data locations to the config file
for datalist in [non_o104seqs, o104seqs]:
    for line in datalist:
        line[3] = line[2].split('.')[0] + '.prodigalout'
        line[4] = line[2].split('.')[0] + '.eprimer3'

### WRITE OUT CONFIG FILE

outfhandle = open('O104_primers_5.conf', 'w')
print >> outfhandle, conf_file_header
print >> outfhandle, '\n'.join(['\t'.join(e) for e in o104seqs])
print >> outfhandle, '\n'.join(['\t'.join(e) for e in non_o104seqs])
outfhandle.close()
