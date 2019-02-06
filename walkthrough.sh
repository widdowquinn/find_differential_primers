#!/usr/bin/env bash
#
# Description: carries out walkthrough as described in `pdp` documentation
# Usage: walkthrough.sh

# 1. Clean walkthrough output
OUTDIR=tests/walkthrough
rm ${OUTDIR}/*.json
rm -rf ${OUTDIR}/blastn* ${OUTDIR}/classify* ${OUTDIR}/deduped* ${OUTDIR}/eprimer3* ${OUTDIR}/primersearch* ${OUTDIR}/prodigal ${OUTDIR}/prodigaligr

# 2. Standard workflow (no filtering)
# Validate config file
pdp config --validate ${OUTDIR}/pectoconf.tab

# Fix input sequences
pdp config --fix_sequences ${OUTDIR}/fixed.json \
    --outdir ${OUTDIR}/config \
    ${OUTDIR}/pectoconf.tab

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

# Screen deduped primers against BLAST db
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
printf "\n"
cat ${OUTDIR}/extract/atrosepticum_NCBI_primers/distances_summary.tab
printf "\n"


# 3. Only use primers that are located in CDS regions
# Here, we can start from the fixed.json file created above
pdp filter -f \
    --prodigal \
    --outdir ${OUTDIR}/prodigal \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/filter_prodigal.json

pdp eprimer3 -f \
    --filter \
    --outdir ${OUTDIR}/eprimer3 \
    ${OUTDIR}/filter_prodigal.json \
    ${OUTDIR}/with_cds_primers.json

pdp dedupe -f \
    --dedupedir ${OUTDIR}/deduped_cds \
    ${OUTDIR}/with_cds_primers.json \
    ${OUTDIR}/deduped_cds_primers.json

pdp blastscreen -f \
    --db ${OUTDIR}/blastdb/e_coli_screen.fna \
    --outdir ${OUTDIR}/blastn_cds \
    ${OUTDIR}/deduped_cds_primers.json \
    ${OUTDIR}/screened_cds.json

pdp primersearch -f \
    --outdir ${OUTDIR}/primersearch_cds \
    ${OUTDIR}/screened_cds.json \
    ${OUTDIR}/primersearch_cds.json

pdp classify -f \
    ${OUTDIR}/primersearch_cds.json \
    ${OUTDIR}/classify_cds

pdp extract -f \
    ${OUTDIR}/primersearch_cds.json \
    ${OUTDIR}/classify_cds/atrosepticum_NCBI_primers.json \
    ${OUTDIR}/extract_cds

cat ${OUTDIR}/classify_cds/summary.tab
printf "\n"
cat ${OUTDIR}/extract_cds/atrosepticum_NCBI_primers/distances_summary.tab
printf "\n"



# 4. Only use primers that are located in intergenic regions
# Again, we can start from the fixed.json file created above
pdp filter -f \
    --prodigaligr \
    --outdir ${OUTDIR}/prodigaligr \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/filter_prodigaligr.json

pdp eprimer3 -f \
    --filter \
    --outdir ${OUTDIR}/eprimer3 \
    ${OUTDIR}/filter_prodigaligr.json \
    ${OUTDIR}/with_igr_primers.json

pdp dedupe -f \
    --dedupedir ${OUTDIR}/deduped_igr \
    ${OUTDIR}/with_igr_primers.json \
    ${OUTDIR}/deduped_igr_primers.json

pdp blastscreen -f \
    --db ${OUTDIR}/blastdb/e_coli_screen.fna \
    --outdir ${OUTDIR}/blastn_igr \
    ${OUTDIR}/deduped_igr_primers.json \
    ${OUTDIR}/screened_igr.json

pdp primersearch -f \
    --outdir ${OUTDIR}/primersearch_igr \
    ${OUTDIR}/screened_igr.json \
    ${OUTDIR}/primersearch_igr.json

pdp classify -f \
    ${OUTDIR}/primersearch_igr.json \
    ${OUTDIR}/classify_igr

pdp extract -f \
    ${OUTDIR}/primersearch_igr.json \
    ${OUTDIR}/classify_igr/atrosepticum_NCBI_primers.json \
    ${OUTDIR}/extract_igr

cat ${OUTDIR}/classify_igr/summary.tab
printf "\n"
cat ${OUTDIR}/extract_cds/atrosepticum_NCBI_primers/distances_summary.tab
printf "\n"
