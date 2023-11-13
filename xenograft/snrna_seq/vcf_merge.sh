#!/bin/bash

# Project Details:
# Line details: /Users/zacc/USyd/ASAP_snRNAseq/ASAPcells_workflows.xlsx
# Line used in initial snRNAseq experiment: rm3.5 pitx3-egfp 

# Inputs
PPMI_WES_vcf="ppmi.feb.1.2015_liftedhg19tohg38.vcf" # "chr1 format" (707183 lines)
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

# 3. re-sort
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
# check file
grep -v "^#" merged_ppmi.feb.1.2015_RM35.vcf | wc -l
bcftools query -l merged_ppmi.feb.1.2015_RM35.vcf
# zip 
gzip merged_ppmi.feb.1.2015_RM35.vcf

# 5b. cleanup temporary files
rm ${RM35_WGS_vcf%%.vcf.gz}_normsort.vcf.gz ${PPMI_WES_vcf%%.vcf.gz}_normsort.vcf.gz
rm *wesvars.vcf *_norm.vcf *_normsort.vcf *normsort.vcf.gz
rm exome_variants.txt

# 6. correct chromosome nomenclature (from 'chr1', to '1')
# rm chr_name_conv.txt
# echo -e '\t' "chr1\t1" >> chr_name_conv.txt
# echo -e '\t' "chr10\t10" >> chr_name_conv.txt
# echo -e '\t' "chr11\t11" >> chr_name_conv.txt
# echo -e '\t' "chr12\t12" >> chr_name_conv.txt
# echo -e '\t' "chr13\t13" >> chr_name_conv.txt
# echo -e '\t' "chr14\t14" >> chr_name_conv.txt
# echo -e '\t' "chr15\t15" >> chr_name_conv.txt
# echo -e '\t' "chr16\t16" >> chr_name_conv.txt
# echo -e '\t' "chr17\t17" >> chr_name_conv.txt
# echo -e '\t' "chr18\t18" >> chr_name_conv.txt
# echo -e '\t' "chr19\t19" >> chr_name_conv.txt
# echo -e '\t' "chr2\t2" >> chr_name_conv.txt
# echo -e '\t' "chr20\t20" >> chr_name_conv.txt
# echo -e '\t' "chr21\t21" >> chr_name_conv.txt
# echo -e '\t' "chr22\t22" >> chr_name_conv.txt
# echo -e '\t' "chr3\t3" >> chr_name_conv.txt
# echo -e '\t' "chr4\t4" >> chr_name_conv.txt
# echo -e '\t' "chr5\t5" >> chr_name_conv.txt
# echo -e '\t' "chr6\t6" >> chr_name_conv.txt
# echo -e '\t' "chr7\t7" >> chr_name_conv.txt
# echo -e '\t' "chr8\t8" >> chr_name_conv.txt
# echo -e '\t' "chr9\t9" >> chr_name_conv.txt
# echo -e '\t' "chrMT\tMT" >> chr_name_conv.txt
# echo -e '\t' "chrX\tX" >> chr_name_conv.txt
# echo -e '\t' "chrY\tY" >> chr_name_conv.txt

# bcftools annotate --rename-chrs chr_name_conv.txt merged_ppmi.feb.1.2015_RM35.vcf > merged_ppmi.feb.1.2015_RM35_2.vcf


