#!/usr/bin/env bash
#
# Description: carries out walkthrough as described in `pdp` documentation
# Usage: walkthrough.sh

# Clean walkthrough output
OUTDIR=tests/walkthrough
rm ${OUTDIR}/*.json
rm -rf ${OUTDIR}/blastn ${OUTDIR}/classify ${OUTDIR}/deduped ${OUTDIR}/eprimer3 ${OUTDIR}/primersearch ${OUTDIR}/prodigal ${OUTDIR}/prodigaligr

# Validate config file
pdp config --validate ${OUTDIR}/pectoconf.tab

# Fix input sequences
pdp config --fix_sequences ${OUTDIR}/fixed.json \
    --outdir ${OUTDIR}/config \
    ${OUTDIR}/pectoconf.tab

# Filter on CDS (not in text walkthrough)
pdp filter -f \
    --prodigal \
    --outdir ${OUTDIR}/prodigal \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/filter_prodigal.json

# Filter on IGR (not in text walkthrough)
pdp filter -f \
    --prodigaligr \
    --outdir ${OUTDIR}/prodigaligr \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/filter_prodigaligr.json

# Design primers
pdp eprimer3 -f \
    --outdir ${OUTDIR}/eprimer3 \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/with_primers.json

# Remove redundant primers
pdp dedupe -f \
    --dedupedir ${OUTDIR}/deduped \
    ${OUTDIR}/with_primers.json \
    ${OUTDIR}/deduped_primers.json

# Screen primers against BLAST db
pdp blastscreen -f \
    --db ${OUTDIR}/blastdb/e_coli_screen.fna \
    --outdir ${OUTDIR}/blastn \
    ${OUTDIR}/deduped_primers.json \
    ${OUTDIR}/screened.json

# Cross-hybridise primers against input genomes
pdp primersearch -f \
    --outdir ${OUTDIR}/primersearch \
    ${OUTDIR}/screened.json \
    ${OUTDIR}/primersearch.json

# Classify primers
pdp classify -f \
    ${OUTDIR}/primersearch.json \
    ${OUTDIR}/classify        

# Extract amplicons
pdp extract -f \
    ${OUTDIR}/primersearch.json \
    ${OUTDIR}/classify/atrosepticum_NCBI_primers.json \
    ${OUTDIR}/extract   

# Display results
cat ${OUTDIR}/classify/summary.tab
echo "\n"