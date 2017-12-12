#!/usr/bin/env bash
#
# Script to generate example files for PrimerSearch tests

pdp.py config --validate ps_prep.conf -v
pdp.py config --to_json ps_prep.json ps_prep.conf -v
pdp.py prodigal --outdir prodigal ps_prep.json ps_features.json -v
pdp.py eprimer3 --outdir eprimer3 ps_features.json ps_primers.json -v
pdp.py blastscreen --db ../walkthrough/blastdb/e_coli_screen.fna --outdir blastn ps_primers.json ps_screened.json -v
pdp.py primersearch --outdir primersearch ps_screened.json ps_primersearch.json -v
