#
# Example makefile for find_differential_primers.py
#
# This makefile serves the following use-case:
#
# 0. Amplification of a specific subset of sequences from a larger set of
#    sequences, e.g. a gene complement
#
# Subject to the following requirements
#
# 1. There is a subset of predicted/annotated genes to which you wish to 
#    design individual diagnostic primers.
# 2. You have a FASTA file containing all sequences - e.g. a predicted
#    gene complement - those you do and do not want to amplify.
# 3. You have a list of sequence IDs for those sequences that you wish
#    to amplify.
# 4. All sequences are prefixed with the same string.
#
# The makefile uses MEME's fasta-fetch to extract the positive 
# sequence set, and builds a suitable configuration file for 
# find_differential_primers, to design primers to these sequences.
# It then constructs the set of negative sequences (those you don't
# want to amplify) as the remaining sequences from the input set,
# and builds a BLAST+ database for screening primers.
# It then runs find_differential_primers.py with suitable parameters.
# Finally, the files describing sequence-specific primers, the 
# corresponding amplicons, and the final table of results, are 
# separated off, and compressed for email or other distribution.
#
# TO RUN:
# At the command-line, run with ‘make’
#  $ make
#
# DEPENDENCIES:
# - The Python script find_differential_primers.py must be installed and 
#   visible in the $PATH, with all its dependencies (BioPython, EMBOSS,
#   BLAST+, etc.) satisfied. The script may be downloaded from 
#   https://github.com/widdowquinn/find_differential_primers
# - MEME must be installed, with the script fasta-fetch being visible
#   in the $PATH.
#   Note that, at least up until version 4.9.1 of MEME, the patch in 
#   accompanying file fasta-fetch.diff must be applied to the version
#   of fasta-fetch in the $PATH.
# - The companion Python scripts build_config_file.py, and 
#   build_negative_sequence_list.py must be present in the running
#   directory.
#
# (c) The James Hutton Institute 2013
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

########
# VARIABLES
#
# NOTE: These variables may be extracted to an additional file, 
#       e.g. parameters.txt, and included into this file with an
#       include directive, e.g. "include parameters.txt". This
#       would allow for defining multiple primer design runs, 
#       but avoiding the need to edit this makefile directly
#
# $(fasta_all_sequences) holds the location of the input multiple
#   FASTA sequence file that describes all input positive and negative
#   (screening) nucleotide sequences, e.g. a whole gene complement
fasta_all_sequences = phytophthora_infestans_t30-4_1_genes.fasta
# $(positive_sequence_list) holds the location of a plaintext file,
#   containing one sequence ID per line, where the sequence IDs are
#   in $(fasta_all_sequences), and define the positive sequence input
#   set.
positive_sequence_list = infestans_subset_list.txt
# $(prefix) holds a string that is the prefix to all input sequences,
#   but is only critical for the input positive sequence set
prefix = PITG_
# $(dir_input) holds the location of a directory that will contain 
#   individual FASTA files for each 'positive' sequence for which a 
#   primer set will be designed
dir_input = input_sequences
# $(dir_output) holds the location of a directory that will take 
#   the complete output from find_differential_primers
dir_output = candidate_primers
# $(dir_results) holds the location of a directory that will take 
#   a subset of the complete output: the results table, and sequence-
#   specific primers
dir_results = candidate_sequence-specific_results
# $(conf_filename) holds the location of the configuration file that 
#   will be generated for this analysis
conf_filename = candidate_primers.conf


all: data config primers

clean: 
# Clean up old input/output sequence directories, intermediate and config files
	rm -rf $(dir_input)
	rm -rf $(dir_output)
	rm -rf $(dir_results)
	rm -f $(conf_filename)
	rm -f negative_sequence_set.fasta.n??
	rm -f $(fasta_all_sequences).index

data: $(fasta_all_sequences).index
# Extract subset of input FASTA sequences, on the basis of the input gene list
# of 'positive' sequences. These are moved to a new $(dir_input) directory, for 
# neatness. 
# This uses MEME's fasta-fetch script.
# All sequences should have a sequence ID that begins with $(prefix).
	fasta-fetch $(fasta_all_sequences) -f $(positive_sequence_list) -s fasta
	mkdir -p $(dir_input)
	mv $(prefix)* $(dir_input)

$(fasta_all_sequences).index: $(fasta_all_sequences)
# Create index of input FASTA sequences, using MEME's fasta-make-index script.
# We use -f to force rewriting of this index.
	fasta-make-index $(fasta_all_sequences) -f

config: $(conf_filename)
# Dummy build for the config file

$(conf_filename):
# Build a configuration file for find_differential_primers.py from the
# files in the input sequence directory
	python build_config_file.py -i $(dir_input) -o $(conf_filename) -f -v -l build_config_file.log

negative_sequence_set.fasta: negative_sequence_set.txt
# Build a FASTA file of sequences from whole genome genecalls not in the input
# sequence list
	fasta-fetch $(fasta_all_sequences) -f negative_sequence_set.txt > negative_sequence_set.fasta

negative_sequence_set.txt:
# Build a list of sequences from the whole genome genecalls that aren't in 
# the input sequence list
	python build_negative_sequence_list.py -f $(fasta_all_sequences) \
		-i $(positive_sequence_list) -o negative_sequence_set.txt -v -l build_negative_ids.log

primers: blastdb
# Run the find_differential_primers.py script, using the generate config file.
# This is a decent set of thermodynamic parameters, but it's not universal,
# so may need to be changed for specific applications
	find_differential_primers.py -i $(conf_filename) \
		--nocds --noprodigal --numreturn 50 --osize=20 --minsize=18 --maxsize=25 \
		--otm=60 --mintm=57 --maxtm=63 --psizeopt=150 --psizemin=100 --psizemax=200 \
		-o $(dir_output) --blastdb=negative_sequence_set.fasta -v > find_differential_primers.log
# Copy the key results into their own subdirectory and zip, for delivery to
# interested parties. We only copy the sequence-specific results, not 
# family-specific.
	rm -rf $(dir_results)
	mkdir -p $(dir_results)
	cp $(dir_output)/*_specific_* $(dir_results)
	cp $(dir_output)/differential_primer_results.tab $(dir_results)
	zip -rv $(dir_results).zip $(dir_results)

blastdb: negative_sequence_set.fasta
# Build a BLAST database of negative examples for screening
	makeblastdb -in negative_sequence_set.fasta -dbtype nucl -title "Negative sequences for screening"
