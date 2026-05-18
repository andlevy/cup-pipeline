#!/bin/bash

# =============================================================================
# nf-core/sarek 3.8.1 -- Test Run Script
# =============================================================================
# Runs the built-in test profile to verify the pipeline, containers, and
# environment are all working correctly before running on real data.
#
# Usage: bash run_sarek_test.sh
# =============================================================================

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline

# --- Environment setup -------------------------------------------------------

export NXF_SINGULARITY_CACHEDIR=${BASE}/Singularity-cache/nf-sarek/sarek-3.8.1

# --- Paths -------------------------------------------------------------------

PIPELINE=${BASE}/CUP_WES_Pipeline/sarek-3.8.1/3_8_1/main.nf
OUTDIR=${BASE}/Results/sarek_test

# --- Run ---------------------------------------------------------------------

nextflow run ${PIPELINE} \
    -profile test,singularity \
    --outdir ${OUTDIR}