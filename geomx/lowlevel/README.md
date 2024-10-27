## ASAP-SpatialTranscriptomics - Nanostring GeoMx processing pipeline

This repository contains analysis scripts used by The University of Sydney collaborators to analyse spatial transcriptomics as part of [Aligning Science Accross Parkinson's](https://parkinsonsroadmap.org/#).

NOTE: This directory contains scripts used in the processing of Nanostring GeoMx data for the following datasets;
1. Spatial Transcriptomics data (GeoMx) of midbrain dopamine cells in control and PD subjects [10.5281/zenodo.13626106](10.5281/zenodo.13626106)
2. Spatial Transcriptomics data (GeoMx) of midbrain tissue in control and PD subjects [10.5281/zenodo.13626167](10.5281/zenodo.13626167)
3. Spatial Transcriptomics data (GeoMx) of locus coeruleus dopamine cells in control and PD subjects [10.5281/zenodo.13626177](10.5281/zenodo.13626177)

## Low Level Data Processing of Nanostring GeoMx
### Overview
The scripts utilize several popular packages for analyzing spatial transcriptomic data including [Seurat](https://satijalab.org/seurat/) and [GeomxTools](https://github.com/Nanostring-Biostats/GeomxTools). 


1. FASTQ files are uploaded to Illumina Basespace using Illumina BaseSpace Sequence Hub command line interface (CLI), [basespace_CLI.sh](/geomx/lowlevel/basespace_CLI.sh), and DCC files files were generated using the [GeoMx NGS Pipeline](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/nanostring-geomxr-ngs-pipeline.html). 

2. [quality_control_gx.R](/geomx/lowlevel/quality_control_gx.R) - Scripts for reading and performing quality control of GeoMx data including ROI and RNA filtering and generation of QC plots.

3. [normalize.batch.cluster_gx.R](/geomx/lowlevel/normalize.batch.cluster_gx.R) - Scripts for the Normalization (please refer to [normalization_review_gx.R](/geomx/lowlevel/normalization_review_gx.R) for selection of normalization process). The scripts also perform batch correction and clustering.


### License
This research was funded in whole or in part by Aligning Science Across Parkinson’s (ASAP-020529) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF). For the purpose of open access, the author has applied a CC BY 4.0 public copyright license to all Author Accepted Manuscripts arising from this submission.
