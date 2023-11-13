#!/bin/bash
#PBS -P FFbigdata
#PBS -N cellranger
#PBS -l select=1:ncpus=8:mem=128GB
#PBS -l walltime=48:00:00

#################################
# ## USYD HPC - note cellranger is not compatible with CentOS6.9.
#################################
# module load cellranger
# source $CellRangerEnvFile
# # # input
# run_dir="/project/RDS-SMS-FFbigdata-RW/spatial_transcriptomics/xenograft/scrna_seq/run_human"
# fastq_files="/project/RDS-SMS-FFbigdata-RW/spatial_transcriptomics/xenograft/scrna_seq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF221012416_HJ377DRX2/12w"
# ID="HJ377DRX2_12w"

# # # get pre-built human reference
# # cd /project/RDS-SMS-FFbigdata-RW/spatial_transcriptomics/xenograft/scrna_seq/run_human/reference
# # wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# # tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz

# # run cellranger
# cd $run_dir
# cellranger count --id=run_${ID} \
#    --fastqs=$fastq_files \
#    --transcriptome="/project/RDS-SMS-FFbigdata-RW/spatial_transcriptomics/xenograft/scrna_seq/run_human/reference/refdata-gex-GRCh38-2020-A"
#################################

#################################
# ## GPHPC - Ubuntu 22.04 LTS  - TESTING
#################################
# # download and install
# cd /home/zac/lib
# curl -o cellranger-7.1.0.tar.xz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.xz?Expires=1694687368&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIueHoiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2OTQ2ODczNjh9fX1dfQ__&Signature=EvbKm6WmZlx2gjueaKNkIimf~iH6D0vhOUUWvt8HcggZnS4ZXF5hVSM913fciuQ3CQE0tKyY8m9n0if36ubMgvHBEXYzKZ9B8Bk2BR62SYzxPZFRslBazufcC3IBVh3~1L14-RxPVme1enTBr8xTjIug2lQEUsZmK8IgnR4QhW92uGGy1WFZA~NK4YWT2MONf7zlytyuHU7Ys1GJcnurkWCOO0oHtrcuJUUdSY5bKA9NMz5fjlE1JspkZp8LIvoKOR9RiHrIE9L3wtSBvgy-5zVvEoRRXxLk8akGYVscCaLE--0AW~bU~mptyIn-a4j8qyME5oSOB8RcNY9ovu8RGw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
# tar xf cellranger-7.1.0.tar.xz
# # export path
# export PATH=/home/zac/lib/cellranger-7.1.0:$PATH
# # perform check
# cellranger sitecheck > sitecheck.txt
# cellranger upload zacchatterton@sydney.edu.au sitecheck.txt
# # verify install 
# cellranger testrun --id=tiny

# # human reference data
# cd /data/zac/asap_snrnaseq/reference
# curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# tar xvf refdata-gex-GRCh38-2020-A.tar.gz

## index reference
# cellranger mkref --genome=GRCh38_and_mRatBN7_281023 --fasta=GRCh38_and_mRatBN7/fasta/genome.fa --genes=GRCh38_and_mRatBN7/genes/genes.gtf.gz

# # running with test data  - yes working
# # # input
# run_dir="/data/zac/asap_snrnaseq/snrnaseq/test_12w"
# fastq_files="/data/zac/asap_snrnaseq/snrnaseq/test_12w"
# ID="HJ377DRX2_12w"
# export PATH=/home/zac/lib/cellranger-7.1.0:$PATH

# # run cellranger
# cd $run_dir
# cellranger count --id=run_${ID} \
#    --fastqs=$fastq_files \
#    --transcriptome="/data/zac/asap_snrnaseq/reference/refdata-gex-GRCh38-2020-A"

#################################
# ## GPHPC - Ubuntu 22.04 LTS - REAL SAMPLES
#################################
# # running with single sample - yes working with human genome, repeating with barnyard genome from 10X (human+mouse) and our own (human+rat)
export PATH=/home/zac/lib/cellranger-7.1.0:$PATH
run_dir="/data/zac/asap_snrnaseq/snrnaseq/"
cd $run_dir

# # 12w_HJ377DRX2 - human/rat - screen -r 1010991
fastq_files="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF221012416_HJ377DRX2/12w"
ID="HJ377DRX2_12w_hr"
transcriptome_ref="/data/zac/asap_snrnaseq/reference/GRCh38_and_mRatBN7_281023"

cellranger count --id=run_${ID} --localcores 20 --localmem 64 \
   --fastqs=$fastq_files \
   --transcriptome=$transcriptome_ref

# # 6w_HJ377DRX2 - human/rat
fastq_files="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF221012416_HJ377DRX2/6w"
ID="HJ377DRX2_6w_hr"

cellranger count --id=run_${ID} --localcores 20 --localmem 64 \
   --fastqs=$fastq_files \
   --transcriptome=$transcriptome_ref

# # 6w_HLNGJDRX2 - human/rat
fastq_files="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/6w"
ID="HLNGJDRX2_6w_hr"

cellranger count --id=run_${ID} --localcores 20 --localmem 64 \
   --fastqs=$fastq_files \
   --transcriptome=$transcriptome_ref

# # 1y_HLNGJDRX2 - human/mouse - DONE
fastq_files="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/1y"
ID="HLNGJDRX2_1y_hm"
transcriptome_ref="/data/zac/asap_snrnaseq/reference/refdata-gex-GRCh38-and-mm10-2020-A"

cellranger count --id=run_${ID} --localcores 20 --localmem 64 \
   --fastqs=$fastq_files \
   --transcriptome=$transcriptome_ref

# # moved output files to local for downstream R analysis
# rsync -av zac@10.65.82.197:/data/zac/asap_snrnaseq/snrnaseq/run_HJ377DRX2_12w_hr/ /Users/zacc/USyd/ASAP/snrna_seq/cellranger_out/run_HJ377DRX2_12w_hr --exclude outs/*bam
# rsync -av zac@10.65.82.197:/data/zac/asap_snrnaseq/snrnaseq/run_HJ377DRX2_6w_hr/ /Users/zacc/USyd/ASAP/snrna_seq/cellranger_out/run_HJ377DRX2_6w_hr --exclude outs/*bam
# rsync -av zac@10.65.82.197:/data/zac/asap_snrnaseq/snrnaseq/run_HLNGJDRX2_6w_hr/ /Users/zacc/USyd/ASAP/snrna_seq/cellranger_out/run_HLNGJDRX2_6w_hr --exclude outs/*bam
# rsync -av zac@10.65.82.197:/data/zac/asap_snrnaseq/snrnaseq/run_HLNGJDRX2_1y_hm/ /Users/zacc/USyd/ASAP/snrna_seq/cellranger_out/run_HLNGJDRX2_1y_hm --exclude outs/*bam





















