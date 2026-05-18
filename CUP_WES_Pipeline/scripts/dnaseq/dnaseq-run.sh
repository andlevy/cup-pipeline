#!/bin/bash
#Version: v0.1.2

# Set up paths
parent_dir="/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/nf-icr-testdatasets"
data_dir="$parent_dir/data"
scripts_dir="$parent_dir/nf-scripts"
output_dir="$parent_dir/nf-results/test-trim"
references_dir="$scripts_dir/references"

# Select the pool to process
pool=3827

# Run the Nextflow pipeline
pool_dir="$data_dir/Pool_$pool"
sample_pth="$pool_dir/samplesheet.csv"
output_pth="$output_dir/Pool_$pool"

cd "$pool_dir"
# Run the Nextflow pipeline
nextflow run nf-core/sarek \
    -profile docker \
    --input "$sample_pth" \
    --outdir "$output_pth" \
    --intervals "$references_dir/KAPA-HyperExomeV2/HyperExomeV2_hg38/HyperExomeV2_primary_targets.bed" \
    --pon "$references_dir/PoN/1000g_pon.hg38.vcf.gz" \
    --aligner bwa-mem2 \
    --wes \
    --trim_fastq \
    --trim_nextseq \
    --use_gatk_spark markduplicates,baserecalibrator \
    --tools mutect2,freebayes,lofreq \
    --filter_vcfs \
    --normalize_vcfs \
    --snv_consensus_calling \
    --consensus_min_count 2 \
    -resume \
    -r "3.8.1"\



