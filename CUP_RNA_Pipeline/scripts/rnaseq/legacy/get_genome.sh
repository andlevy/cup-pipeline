parent_dir="/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/nf-icr-testdatasets"
references_dir="$parent_dir/references/ensembl_115"

mkdir -p "$references_dir"
cd "$references_dir"

VERSION=115
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz