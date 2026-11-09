#!/bin/bash
# nf-core/rnaseq 3.23.0 -- One-time pipeline and container download
# Run once with internet access. After this, all runs are fully offline.

set -e

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline

source /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/etc/profile.d/conda.sh
conda activate nfcore

export NXF_SINGULARITY_CACHEDIR=${BASE}/Singularity-cache/nf-rnaseq/rnaseq-3.23.0

mkdir -p ${NXF_SINGULARITY_CACHEDIR}

nf-core pipelines download rnaseq \
    -r 3.23.0 \
    --container-system singularity \
    --compress none \
    --force \
    --outdir ${BASE}/CUP_RNA_Pipeline/rnaseq-3.23.0

echo ""
echo "Pipeline downloaded to: ${BASE}/CUP_RNA_Pipeline/rnaseq-3.23.0/3_23_0/main.nf"
echo "Singularity images at:  ${NXF_SINGULARITY_CACHEDIR}"
