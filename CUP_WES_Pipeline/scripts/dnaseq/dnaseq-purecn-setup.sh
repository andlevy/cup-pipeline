#!/bin/bash
# =============================================================================
# PureCN -- Setup Script
# =============================================================================
# Generates the PureCN interval file from the KAPA HyperExome V2 BED file
# and the GRCh38 reference FASTA. This only needs to be run once.
#
# Usage: bash dnaseq-purecn-setup.sh
# Version: v0.1.0
# =============================================================================

set -e

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline
REFS=${BASE}/References

# --- Activate PureCN conda environment ---------------------------------------

source /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/etc/profile.d/conda.sh
conda activate purecn

# --- Paths -------------------------------------------------------------------

BED=${REFS}/KAPA-HyperExomeV2/HyperExomeV2_hg38/HyperExomeV2_capture_targets.bed
FASTA=${REFS}/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
OUTDIR=${REFS}/purecn
OUTFILE=${OUTDIR}/HyperExomeV2_hg38_intervals.txt

# --- Create output directory -------------------------------------------------

mkdir -p ${OUTDIR}

# --- Get PureCN IntervalFile.R script location -------------------------------

INTERVAL_SCRIPT=$(Rscript -e "cat(system.file('extdata', 'IntervalFile.R', package='PureCN'))")

echo "=============================================="
echo " PureCN Setup"
echo " Interval script: ${INTERVAL_SCRIPT}"
echo " BED file:        ${BED}"
echo " FASTA:           ${FASTA}"
echo " Output:          ${OUTFILE}"
echo " Started:         $(date)"
echo "=============================================="
echo ""

# --- Generate interval file --------------------------------------------------

Rscript ${INTERVAL_SCRIPT} \
    --infile  ${BED} \
    --fasta   ${FASTA} \
    --outfile ${OUTFILE} \
    --genome  hg38 \
    --force

echo ""
echo "=============================================="
echo " PureCN setup complete: $(date)"
echo " Interval file: ${OUTFILE}"
echo "=============================================="
