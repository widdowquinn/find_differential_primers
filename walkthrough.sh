#!/usr/bin/env bash
#
# Description: carries out walkthrough as described in `pdp` documentation,
# and provides alternative workflows, where primer design locations are filtered.
# Usage: walkthrough.sh

# Which primer design tool are we using
PRIMERTOOL=primer3
# Uncomment below to use eprimer3 instead
# PRIMERTOOL=eprimer3

# THERMODYNAMIC PARAMETERS FOR PRIMER3
PRIMER3PARAMS=tests/test_input/primer3/primer3_config

# 1. Clean walkthrough output
OUTDIR=tests/walkthrough
DATADIR=tests/walkthrough_data
rm -rf ${OUTDIR}
mkdir -p ${OUTDIR}

# 2. Standard workflow (no filtering)
# Validate config file
pdp config --validate ${DATADIR}/pectoconf.tab

# Fix input sequences
pdp config --fix_sequences ${OUTDIR}/fixed.json \
    --outdir ${OUTDIR}/config \
    ${DATADIR}/pectoconf.tab

# Design primers
mkdir -p ${OUTDIR}/${PRIMERTOOL}
pdp ${PRIMERTOOL} -f \
    --outdir ${OUTDIR}/${PRIMERTOOL} \
    --therm_param_path ${PRIMER3PARAMS} \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/with_primers.json

# Remove redundant primers
pdp dedupe -f \
    --dedupedir ${OUTDIR}/deduped \
    ${OUTDIR}/with_primers.json \
    ${OUTDIR}/deduped_primers.json

# Screen deduped primers against BLAST db
pdp blastscreen -f \
    --db ${DATADIR}/blastdb/e_coli_screen.fna \
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
    ${OUTDIR}/classify/Pectobacterium_primers.json \
    ${OUTDIR}/extract

# Display results
cat ${OUTDIR}/classify/summary.tab
printf "\n"
cat ${OUTDIR}/extract/Pectobacterium_primers/distances_summary.tab
printf "\n"


# 3. Only use primers that are located in CDS regions
# Here, we can start from the fixed.json file created above
pdp filter -f \
    --prodigal \
    --outdir ${OUTDIR}/prodigalcds \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/filter_prodigal.json

pdp ${PRIMERTOOL} -f \
    --filter \
    --outdir ${OUTDIR}/${PRIMERTOOL}_cds \
    --therm_param_path ${PRIMER3PARAMS} \
    ${OUTDIR}/filter_prodigal.json \
    ${OUTDIR}/with_cds_primers.json

pdp dedupe -f \
    --dedupedir ${OUTDIR}/deduped_cds \
    ${OUTDIR}/with_cds_primers.json \
    ${OUTDIR}/deduped_cds_primers.json

pdp blastscreen -f \
    --db ${DATADIR}/blastdb/e_coli_screen.fna \
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
    ${OUTDIR}/classify_cds/Pectobacterium_primers.json \
    ${OUTDIR}/extract_cds

cat ${OUTDIR}/classify_cds/summary.tab
printf "\n"
cat ${OUTDIR}/extract_cds/Pectobacterium_primers/distances_summary.tab
printf "\n"



# 4. Only use primers that are located in intergenic regions
# Again, we can start from the fixed.json file created above
pdp filter -f \
    --prodigaligr \
    --outdir ${OUTDIR}/prodigaligr \
    ${OUTDIR}/fixed.json \
    ${OUTDIR}/filter_prodigaligr.json

pdp ${PRIMERTOOL} -f \
    --filter \
    --outdir ${OUTDIR}/${PRIMERTOOL}_igr \
    --therm_param_path ${PRIMER3PARAMS} \
    ${OUTDIR}/filter_prodigaligr.json \
    ${OUTDIR}/with_igr_primers.json

pdp dedupe -f \
    --dedupedir ${OUTDIR}/deduped_igr \
    ${OUTDIR}/with_igr_primers.json \
    ${OUTDIR}/deduped_igr_primers.json

pdp blastscreen -f \
    --db ${DATADIR}/blastdb/e_coli_screen.fna \
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
    ${OUTDIR}/classify_igr/Pectobacterium_primers.json \
    ${OUTDIR}/extract_igr

cat ${OUTDIR}/classify_igr/summary.tab
printf "\n"
cat ${OUTDIR}/extract_igr/Pectobacterium_primers/distances_summary.tab
printf "\n"
