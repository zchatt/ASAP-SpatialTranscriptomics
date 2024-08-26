## ASAP-snrnaseq

This repository contains analysis scripts used by The University of Sydney collaborators to analyse single-nuclei RNA sequencing (snRNAseq) as part of [Aligning Science Accross Parkinson's](https://parkinsonsroadmap.org/#).

NOTE: This repo currently in development and likely to change frequently 
 

## snRNAseq - genomics
### Overview
The following scripts are used to perform variant calling on 10X snRNAseq data from xenograph models. The aim is to establish a method that can deconvolute pooled single-nuclei from unique donors to 1) reduce the 10X snRNAseq library prep costs and/or 2) pool (vallage approach) donor lines when constructing xenograph models.

1. vcf_check.sh - Format and merging of .vcf files
2. souporcell.sh - Variant calling of snRNAseq data

## snRNAseq - transcriptomics
### Overview

1. cellranger.sh - Alignment and gene expression counting of snRNAseq data
2. GRCh38_and_mRatBN7-2023-A_build.sh - Building of merged human and rat reference genome
3. snrnaseq_lowlevel.Rmd - R scripts to define human cells population, DA neuron populations and perform DEG.

### License
This project is licensed under the MIT License.

### Acknowledgments
We would like to acknowledge the developers of related packages that have been utilized in this repository.

Please refer to the individual script files for more detailed information about each script's functionality, input requirements, and output formats.
