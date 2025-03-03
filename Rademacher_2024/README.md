## ASAP-SpatialTranscriptomics

This repository contains analysis scripts used by The University of Sydney collaborators to analyse spatial transcriptomics as part of [Aligning Science Accross Parkinson's](https://parkinsonsroadmap.org/#).

NOTE: The "Rademacher_2024" directory contains scripts used in the analysis of "Chronic hyperactivation of midbrain dopamine neurons causes preferential dopamine neuron degeneration"

## Data Processing and Analysis for Nanostring GeoMx
### Overview
The scripts utilize several popular packages for analyzing spatial transcriptomic data including [Seurat](https://satijalab.org/seurat/) and [GeomxTools](https://github.com/Nanostring-Biostats/GeomxTools). The scripts were used to perform various analyses on GeoMx (GeoMx Hu WTA) data generated from human brain tissue including data preprocessing, quality control, normalisation and assessment of directional changes of identifed mouse DEGs in human tissue.

### The analysis takes the following inputs:

[DCC files](https://zenodo.org/records/10499187): These are raw data files generated using the [GeoMx NGS Pipeline](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/nanostring-geomxr-ngs-pipeline.html). 

OR

[Seurat object](https://zenodo.org/records/10499187): An R object containing processed counts for 9012 features across 73 ROIs (n=18 samples) 

### The analysis should be performed in the following order:

1. [quality_control_gx.R](/Rademacher_2024/R/quality_control_gx.R) - Scripts for reading and performing quality control of GeoMx data including ROI and RNA filtering and generation of QC plots.

2. [normalize.batch.cluster_gx.R](/Rademacher_2024/R/normalize.batch.cluster_gx.R) - Scripts for the Normalization (please refer to [normalization_review_vis.R](/Rademacher_2024/R/normalization_review_gx.R) for selection of normalization process). The scripts also perform batch correction and clustering.

3. [validate_DEG.R](/Rademacher_2024/R/validate_DEG.R) - Scripts used for human validation of DEGs identified by chemogenetic (DREADD) mouse model to chronically hyperactivate of DA neurons.

4. [neuroestimator_run.R](/Rademacher_2024/R/neuroestimator_run.R) - Scripts used to run NEUROeSTIMator (Bahl et al., 2024)(https://www.nature.com/articles/s41467-023-44503-5) on Visium mouse datasets.

5. [neuroestimator_compare.R](/Rademacher_2024/R/neuroestimator_compare.R) - Scripts used to normalize and compare NEUROeSTIMator results generated from Visium assessment of mouse.

4. [neuroestimator_run.R](/R/neuroestimator_run.R) - Scripts used to run NEUROeSTIMator (Bahl et al., 2024)(https://www.nature.com/articles/s41467-023-44503-5) on Visium mouse datasets.

5. [neuroestimator_compare.R](/R/neuroestimator_compare.R) - Scripts used to normalize and compare NEUROeSTIMator results generated from Visium assessment of mouse.

### Usage
To run the script, follow these steps:

1. Install the required dependencies.

2. Download the [data](https://zenodo.org/records/10499187) and place them in the appropriate directories.

3. Open the R script and modify the file paths and parameters as needed.

4. Run each script using an R environment or an integrated development environment (IDE).

5. Review the generated results and output files for further analysis.


### License
This research was funded in whole or in part by Aligning Science Across Parkinson’s (ASAP-020529) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF). For the purpose of open access, the author has applied a CC BY 4.0 public copyright license to all Author Accepted Manuscripts arising from this submission.
