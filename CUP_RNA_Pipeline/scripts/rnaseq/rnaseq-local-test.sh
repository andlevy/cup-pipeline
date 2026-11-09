#!/bin/bash
# =============================================================================
# nf-core/rnaseq 3.23.0 -- Test Run Script
# =============================================================================
# Runs the built-in test profile with synthetic data to verify the pipeline,
# Singularity containers, and environment are all working correctly before
# running on real data. Completes in ~15 minutes.
#
# Note: requires brief internet access to fetch synthetic test data from
# nf-core's S3 bucket. Do not use -offline here.
# =============================================================================

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline

source /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/etc/profile.d/conda.sh
conda activate nfcore

export TMPDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export NXF_TEMP=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export NXF_SINGULARITY_CACHEDIR=${BASE}/Singularity-cache/nf-rnaseq/rnaseq-3.23.0

PIPELINE=${BASE}/CUP_RNA_Pipeline/rnaseq-3.23.0/3_23_0/main.nf
OUTDIR=${BASE}/Results/rnaseq_test

nextflow run ${PIPELINE} \
    -profile test,singularity \
    --outdir ${OUTDIR}
