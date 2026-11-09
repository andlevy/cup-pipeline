#!/bin/bash
# =============================================================================
# CUP RNA Pipeline -- Production Run Script
# nf-core/rnaseq 3.23.0
# =============================================================================
# Tumour RNA-seq: STAR alignment + Salmon quantification
# Reference: Ensembl GRCh38 v115 (local)
# Version: v0.1.0
# =============================================================================

set -e

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline
REFS=${BASE}/References

source /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/etc/profile.d/conda.sh
conda activate nfcore

mkdir -p /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export TMPDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export NXF_TEMP=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export JAVA_OPTS="-Djava.io.tmpdir=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp"
export NXF_SINGULARITY_CACHEDIR=${BASE}/Singularity-cache/nf-rnaseq/rnaseq-3.23.0
export NXF_OFFLINE=true

# --- Update this for each run ------------------------------------------------

POOL=Pool_3493

# -----------------------------------------------------------------------------

SAMPLESHEET=${BASE}/Data/ICR-Test/${POOL}/samplesheet.csv
OUTDIR=${BASE}/Results/${POOL}
PIPELINE=${BASE}/CUP_RNA_Pipeline/rnaseq-3.23.0/3_23_0/main.nf

FASTA=${REFS}/ensembl_115/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
GTF=${REFS}/ensembl_115/Homo_sapiens.GRCh38.115.gtf.gz

echo "=================================================="
echo " CUP RNA Pipeline"
echo " Pool:    ${POOL}"
echo " Outdir:  ${OUTDIR}"
echo " Started: $(date)"
echo "=================================================="
echo ""

nextflow run ${PIPELINE} \
    -profile singularity \
    -c ${BASE}/CUP_RNA_Pipeline/configs/local.config \
    -resume \
    --input ${SAMPLESHEET} \
    --outdir ${OUTDIR} \
    --fasta ${FASTA} \
    --gtf ${GTF} \
    --aligner star_salmon \
    --save_align_intermeds \
    --trimmer trimgalore \
    -offline

echo ""
echo "=================================================="
echo " Pipeline complete: $(date)"
echo "=================================================="
