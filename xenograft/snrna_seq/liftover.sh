#!/bin/bash
#PBS -P FFbigdata
#PBS -N liftover
#PBS -l select=1:ncpus=8:mem=200GB
#PBS -l walltime=4:00:00

# Note 1. The ppmi.feb.1.2015.vcf.gz was aligned against the reference human genome (UCSC hg19). 
# Lifover hg19 to hg38 variants
module load picard
module load samtools

# # get hg38 reference
# cd /scratch/RDS-SMS-FFbigdata-RW/zacc/souporcell
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# # create dictionary
# picard CreateSequenceDictionary \
# R=hg38.fa.gz \
# O=hg38.dict

# Liftover
cd /scratch/RDS-SMS-FFbigdata-RW/zacc/souporcell

picard LiftoverVcf \
I=/project/RDS-SMS-FFbigdata-RW/spatial_transcriptomics/xenograft/ppmi/ppmi.feb.1.2015.vcf.gz \
O=/scratch/RDS-SMS-FFbigdata-RW/zacc/souporcell/ppmi.feb.1.2015_liftedhg19tohg38.vcf \
CHAIN=/project/RDS-SMS-FFbigdata-RW/spatial_transcriptomics/xenograft/scrna_seq/souporcell/hg19ToHg38.over.chain \
REJECT=/scratch/RDS-SMS-FFbigdata-RW/zacc/souporcell/rejected_variants.vcf \
R=hg38.fa.gz \
MAX_RECORDS_IN_RAM=100000 \
WARN_ON_MISSING_CONTIG=true \
DISABLE_SORT=true