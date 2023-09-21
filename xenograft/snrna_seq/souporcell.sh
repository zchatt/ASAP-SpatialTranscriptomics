#!/bin/bash
#PBS -P FFbigdata
#PBS -N souporcell
#PBS -l select=1:ncpus=8:mem=16GB
#PBS -l walltime=4:00:00

# NOTES
# 1. This cannot be run on the Artemis as the cluster does not enable OverlayFS for Singularity containers. Please see 
# 2. The singularity container works on "gbphc" cluster, hence have performed analysis there.

# pull container
# singularity pull shub://wheaton5/souporcell

# # # download test files
# cd /data/zac/asap_snrnaseq/test_out
# wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398235/A.merged.bam.1 -O A.merged.bam
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560245/suppl/GSM2560245_barcodes.tsv.gz
# gunzip GSM2560245_barcodes.tsv.gz
# wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
# tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz

# 1. run test files - yes all working as expected
soup_or_cell_sif="/home/zac/souporcell/souporcell_latest.sif"
num_clusters=4
num_threads_to_use=8
output_dir_name="/data/zac/asap_snrnaseq/souporcell/test_out"
reference_fasta="/data/zac/asap_snrnaseq/reference/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
barcodes_tsv="/data/zac/asap_snrnaseq/souporcell/test_out/GSM2560245_barcodes.tsv"
bam_file="/data/zac/asap_snrnaseq/souporcell/test_out/A.merged.bam"
export SINGULARITY_BIND="/data/zac/"


cd $output_dir_name

singularity exec $soup_or_cell_sif souporcell_pipeline.py \
	-i $bam_file \
	-b $barcodes_tsv \
	-f $reference_fasta \
	-t $num_threads_to_use \
	-o $output_dir_name \
	-k $num_clusters

# 2. run cellranger local .bam file
soup_or_cell_sif="/home/zac/souporcell/souporcell_latest.sif"
num_clusters=4
num_threads_to_use=24
output_dir_name="/data/zac/asap_snrnaseq/souporcell/run_HJ377DRX2_12w"
#reference_fasta="/data/zac/asap_snrnaseq/souporcell/test12w_out/ref/genome.fa"
reference_fasta="/data/zac/asap_snrnaseq/reference/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
barcodes_tsv="/data/zac/asap_snrnaseq/snrnaseq/run_HJ377DRX2_12w/outs/filtered_feature_bc_matrix/barcodes.tsv"
bam_file="/data/zac/asap_snrnaseq/snrnaseq/run_HJ377DRX2_12w/outs/possorted_genome_bam.bam"
export SINGULARITY_BIND="/data/zac/"

cd $output_dir_name

singularity exec $soup_or_cell_sif souporcell_pipeline.py \
	-i $bam_file \
	-b $barcodes_tsv \
	-f $reference_fasta \
	-t $num_threads_to_use \
	-o $output_dir_name \
	-k $num_clusters

# 3. run cellranger local .bam file with variants
soup_or_cell_sif="/home/zac/souporcell/souporcell_latest.sif"
num_clusters=4
num_threads_to_use=24
output_dir_name="/data/zac/asap_snrnaseq/souporcell/run_HJ377DRX2_12w_vcf"
#reference_fasta="/data/zac/asap_snrnaseq/souporcell/test12w_out/ref/genome.fa"
reference_fasta="/data/zac/asap_snrnaseq/reference/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
barcodes_tsv="/data/zac/asap_snrnaseq/snrnaseq/run_HJ377DRX2_12w/outs/filtered_feature_bc_matrix/barcodes.tsv"
bam_file="/data/zac/asap_snrnaseq/snrnaseq/run_HJ377DRX2_12w/outs/possorted_genome_bam.bam"
vcf="/data/zac/asap_snrnaseq/souporcell/merged_ppmi.feb.1.2015_RM35.vcf"
export SINGULARITY_BIND="/data/zac/"

mkdir -p $output_dir_name
cd $output_dir_name

singularity exec $soup_or_cell_sif souporcell_pipeline.py \
	-i $bam_file \
	-b $barcodes_tsv \
	-f $reference_fasta \
	-t $num_threads_to_use \
	-o $output_dir_name \
	-k $num_clusters \
	--known_genotypes $vcf





merged_ppmi.feb.1.2015_RM35.vcf



## Potential method to run on USyd artemis cluster
# Just in case anyone else is trying to run this on a cluster that doesnt enable OverlayFS for Singularity. 
# Can potentially get around this issue by adding the required paths to the .sif file directly by following advice here - https://stackoverflow.com/questions/67851786/edit-runscript-of-singularity-sif-container-after-building
# # 1. convert a SIF file to a (writable) sandbox which is indeed a directory:
# singularity build --sandbox souporcell_sandbox souporcell_latest.sif

# # 2. may need to make writeable eg. sudo singularity shell --writable souporcell_sandbox
# # 3. edit the souporcell_sandbox/environment to include the required paths
# # 4. convert the sandbox back to a (new) SIF file:
# sudo singularity build souporcell_latest_custom.sif souporcell_sandbox
# # note will need to use "module load singularity/3.7.0"



AAACATTGCATGGT-1
CGAACGAAAGCC
AAGGCAGCAGCAAATCCTTCGT



