#!/bin/bash

# Project Details:
# Line details: /Users/zacc/USyd/ASAP_snRNAseq/ASAPcells_workflows.xlsx
# Line used in initial snRNAseq experiment: rm3.5 pitx3-egfp 

# Notes: 
# The RM3.5 is the control iPSC line from the fetal fibroblast parental line and isn’t part of PPMI.
# Ryan provided GATK best practices vcf output for the RM3.5 line. There might be an issue with this but Dario can give you the Dragen output when he can off ICA if need be.

# Inputs
PPMI_WES_vcf="/Users/zacc/USyd/ASAP/snrna_seq/ppmi.feb.1.2015.vcf.gz" # "chr1 format" (707183 lines)
RM35_WGS_vcf="/Users/zacc/USyd/ASAP/snrna_seq/7_21_RM35.g.vcf.gz" # "chr1 format" (177658469 lines)
fasta_ref="/Users/zacc/USyd/ASAP/snrna_seq/genome.fa"
working_dir="/Users/zacc/USyd/ASAP/snrna_seq"
#vcftools="/Users/zacc/vcftools_0.1.13/bin/vcftools"

# 1. check sample name
cd $working_dir
bcftools query -l $RM35_WGS_vcf # 7_21_RM35
bcftools query -l $PPMI_WES_vcf # 645 PPMI_SI IDs

# 2. reduce to only overlapping WES variants (40488 variants)
bcftools query -f'%CHROM\t%POS\n' $PPMI_WES_vcf > exome_variants.txt
bcftools view -T exome_variants.txt $RM35_WGS_vcf > ${RM35_WGS_vcf%%.vcf.gz}_wesvars.vcf #(43901 lines, 40488 variants)
bcftools query -f'%CHROM\t%POS\n' ${RM35_WGS_vcf%%.vcf.gz}_wesvars.vcf > overlap_exome_variants.txt
bcftools view -T overlap_exome_variants.txt $PPMI_WES_vcf > ${PPMI_WES_vcf%%.vcf.gz}_wesvars.vcf 

# re-assign
PPMI_WES_vcf=${PPMI_WES_vcf%%.vcf.gz}_wesvars.vcf
RM35_WGS_vcf=${RM35_WGS_vcf%%.vcf.gz}_wesvars.vcf

# 2. normalize against reference - as the .vcf do not always contain the same REF allele, we need to fix (exclude)
bcftools norm --check-ref x --fasta-ref $fasta_ref $RM35_WGS_vcf > ${RM35_WGS_vcf%%.vcf.gz}_norm.vcf
bcftools norm --check-ref x --fasta-ref $fasta_ref $PPMI_WES_vcf > ${PPMI_WES_vcf%%.vcf.gz}_norm.vcf

# 3. resort
bcftools sort ${RM35_WGS_vcf%%.vcf.gz}_norm.vcf -o ${RM35_WGS_vcf%%.vcf.gz}_normsort.vcf
bcftools sort ${PPMI_WES_vcf%%.vcf.gz}_norm.vcf -o ${PPMI_WES_vcf%%.vcf.gz}_normsort.vcf

# 4. compress 
bcftools view -Oz -o ${RM35_WGS_vcf%%.vcf.gz}_normsort.vcf.gz ${RM35_WGS_vcf%%.vcf.gz}_normsort.vcf
htsfile ${RM35_WGS_vcf%%.vcf.gz}_normsort.vcf.gz
bcftools index ${RM35_WGS_vcf%%.vcf.gz}_normsort.vcf.gz

bcftools view -Oz -o ${PPMI_WES_vcf%%.vcf.gz}_normsort.vcf.gz ${PPMI_WES_vcf%%.vcf.gz}_normsort.vcf
htsfile ${PPMI_WES_vcf%%.vcf.gz}_normsort.vcf.gz
bcftools index ${PPMI_WES_vcf%%.vcf.gz}_normsort.vcf.gz

# 5. merge .vcf files (40577 variants)
bcftools merge ${RM35_WGS_vcf%%.vcf.gz}_normsort.vcf.gz ${PPMI_WES_vcf%%.vcf.gz}_normsort.vcf.gz > merged_ppmi.feb.1.2015_RM35.vcf

