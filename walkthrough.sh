#!/usr/bin/env bash
#
# Description: carries out walkthrough as described in `pdp` documentation
# Usage: walkthrough.sh

# Validate config file
pdp config --validate tests/walkthrough/pectoconf.tab

# Fix input sequences
pdp config --fix_sequences tests/walkthrough/fixed.json tests/walkthrough/pectoconf.tab

# Design primers
pdp eprimer3 -f \
    --outdir tests/walkthrough/eprimer3 \
    tests/walkthrough/fixed.json \
    tests/walkthrough/with_primers.json

# Remove redundant primers
pdp dedupe -f \
    --dedupedir tests/walkthrough/deduped \
    tests/walkthrough/with_primers.json \
    tests/walkthrough/deduped_primers.json

# Screen primers against BLAST db
pdp blastscreen -f \
    --db tests/walkthrough/blastdb/e_coli_screen.fna \
    --outdir tests/walkthrough/blastn \
    tests/walkthrough/deduped_primers.json \
    tests/walkthrough/screened.json

# Cross-hybridise primers against input genomes
pdp primersearch -f \
    --outdir tests/walkthrough/primersearch \
    tests/walkthrough/screened.json \
    tests/walkthrough/primersearch.json

# Classify primers
pdp classify -f \
    tests/walkthrough/primersearch.json \
    tests/walkthrough/classify        

# Display results
cat tests/walkthrough/classify/summary.tab
echo "\n"