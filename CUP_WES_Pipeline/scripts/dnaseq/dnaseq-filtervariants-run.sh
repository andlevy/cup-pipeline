#!/bin/bash
# =============================================================================
# CUP WES Pipeline -- Somatic Variant Filter Runner
# =============================================================================
# Version: v0.2.0 Andrew Levy
# =============================================================================

set -e

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline
POOL=Pool_3827

VCF=${BASE}/Results/${POOL}/annotation/mutect2/S1/S1.mutect2.filtered_VEP.ann.vcf.gz
COSMIC=${BASE}/References/cosmic/Census_all.csv
OUTDIR=${BASE}/Results/${POOL}/annotation/filtered
SCRIPT=${BASE}/CUP_WES_Pipeline/scripts/dnaseq/dnaseq-filtervariants.py

# Mutect2 filter flags to exclude
# Remove any from this list to allow those variants through
EXCLUDE_FILTERS="germline panel_of_normals strand_bias base_qual weak_evidence clustered_events haplotype slippage multiallelic"

echo "=============================================="
echo " CUP WES -- Somatic Variant Filter"
echo " Pool:    ${POOL}"
echo " Started: $(date)"
echo "=============================================="
echo ""

python3 ${SCRIPT} \
    --vcf            ${VCF} \
    --outdir         ${OUTDIR} \
    --min_vaf        0.05 \
    --min_depth      20 \
    --exclude_filters ${EXCLUDE_FILTERS} \
    --no_consequence_filter \
    --report

echo ""
echo "Completed: $(date)"
