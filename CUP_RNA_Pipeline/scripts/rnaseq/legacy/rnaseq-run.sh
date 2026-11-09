#!/bin/bash

# Set up paths
parent_dir="/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/nf-icr-testdatasets"
data_dir="$parent_dir/data"
scripts_dir="$parent_dir/nf-scripts"
output_dir="$parent_dir/nf-results"
references_dir="$parent_dir/references/ensembl_115"

# Select the pool to process
pool=3493

# Run the Nextflow pipeline
pool_dir="$data_dir/Pool_$pool"
sample_pth="$pool_dir/samplesheet.csv"
output_pth="$output_dir/Pool_$pool"

cd "$pool_dir"
nextflow run nf-core/rnaseq \
    --input "$sample_pth" \
    --outdir "$output_pth" \
    --genome GRCh38 \
    --save_align_intermeds \
    -r "3.23.0" \
    -profile docker \
    -resume


# --gtf "$references_dir/Homo_sapiens.GRCh38.115.gtf.gz" \
#     --fasta "$references_dir/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz" \
 