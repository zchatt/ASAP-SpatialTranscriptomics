## ASAP-SpatialTranscriptomics

This repository contains analysis scripts used by The University of Sydney collaborators to analyse spatial transcriptomics as part of [Aligning Science Accross Parkinson's](https://parkinsonsroadmap.org/#).

NOTE: This repo currently in development and likely to change frequently 
 
## Data Processing and Analysis for 10X Visium
### Overview
The scripts utilize several popular packages for analyzing spatial transcriptomic data including [Seurat](https://satijalab.org/seurat/), [Giotto](https://giottosuite.readthedocs.io/en/latest/index.html) and [spacxr](https://github.com/dmcable/spacexr)

The scripts are designed to perform various analyses on Visium data, including data preprocessing, dimension reduction, clustering, cell-type enrichment and cell-type differentially expressed genes (ctDEG's) analysis.

### The analysis takes the following inputs:

Visium data matrix: A matrix containing gene expression values for spatially resolved spots.
Brain tissue image: An image representing the spatial coordinates of the spots in the Visium dataset.

### The analysis should be performed in the following order:

1. [qc.filt.norm_vis.R](/visium/qc.filt.norm_vis.R) - Scripts for reading in and performing quality control of Visium data including manual annotate spots of concern, enabling the removal of unwanted spots from the analysis. The scripts also perform spot and RNA filtering and Normalization (please refer to [normalization_review_vis.R](/visium/normalization_review_vis.R) for selection of normalization process). These scripts should be run per Visium window.

2. [batch.cluster_vis.R](/visium/batch.cluster_vis.R) - Scripts for joining multiple datasets (Visium windows), batch correction and clustering batch correction and clustering of multiple Visium windows.

3. [PAGE_vis.R](/visium/PAGE_vis.R) - Scripts for performing PAGE analysis that calculates an enrichment score based on the fold change of cell type marker genes that have been curated [here](/celltype_markers/format_cell_markers.R)

4. [ctDEG_vis.R](/visium/ctDEG_vis.R) - Scripts for identifying ctDEG's using the [RCTD & C-SIDE](https://github.com/dmcable/spacexr) technique. 

### Usage
To run the script, follow these steps:

1. Install the required dependencies.

2. Prepare your Visium data matrix and brain tissue image and place them in the appropriate directories.

3. Open the R script and modify the file paths and parameters as needed.

4. Run each script using an R environment or an integrated development environment (IDE).

5. Review the generated results and output files for further analysis.


## Data Processing and Analysis for Nanostring GeoMx
### Overview
The scripts utilize several popular packages for analyzing spatial transcriptomic data including [Seurat](https://satijalab.org/seurat/), [Giotto](https://giottosuite.readthedocs.io/en/latest/index.html) and [spacxr](https://github.com/dmcable/spacexr)

The scripts are designed to perform various analyses on Visium data, including data preprocessing, dimension reduction, clustering, cell-type enrichment and cell-type differentially expressed genes (ctDEG's) analysis.

### The analysis takes the following inputs:

DCC files: These are generated using the [GeoMx NGS Pipeline](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/nanostring-geomxr-ngs-pipeline.html). Please refer to [basespace_CLI.sh](/geomx/basespace_CLI.sh) for helpful Illumina BaseSpace Sequence Hub command line interface (CLI) arguments.
Annotation file : Contains ROI metadata

### The analysis should be performed in the following order:

1. [quality_control_gx.R](/geomx/quality_control_gx.R) - Scripts for reading in and performing quality control of GeoMx data including ROI and RNA filtering and generation of QC plots.

2. [normalize.batch.cluster_gx.R](/geomx/normalize.batch.cluster_gx.R) - Scripts for the Normalization (please refer to [normalization_review_vis.R](/geomx/normalization_review_gx.R) for selection of normalization process). The scripts also perform batch correction and clustering.

3. [PAGE_gx.R](/geomx/PAGE_gx.R) - Scripts for performing PAGE analysis that calculates an enrichment score based on the fold change of cell type marker genes that have been curated [here](/celltype_markers/format_cell_markers.R)

4. [ctDEG_gx.R](/geomx/ctDEG_gx.R) - Scripts for identifying ctDEG's using the [RCTD & C-SIDE](https://github.com/dmcable/spacexr) technique. 

### Usage
To run the script, follow these steps:

1. Install the required dependencies.

2. Prepare your Visium data matrix and brain tissue image and place them in the appropriate directories.

3. Open the R script and modify the file paths and parameters as needed.

4. Run each script using an R environment or an integrated development environment (IDE).

5. Review the generated results and output files for further analysis.


## Xenograft
This is a collection of scripts to perform ST and scnRNAseq analysis in xenografts.

1. [human_probe_mm_vis.R](/xenograft/human_probe_mm_vis.R) - Scripts to identify Visium human probes (v2.0) that are capable of quantifying mouse transcripts from the mouse tissue of Xenografts run on the Visium 10X human array. Note; probes were aligned using [align_species.pbs](/xenograft/align_species.pbs). The probe set [Visium_Human_mouseapplicable_Transcriptome_Probe_Set_v2.0_GRCh38-2023-A.csv](/xenograft/Visium_Human_mouseapplicable_Transcriptome_Probe_Set_v2.0_GRCh38-2023-A.csv) can be passed to "spaceranger count" function.

2. [snrna_seq](/xenograft/snrna_seq) - Currently a collection of scripts for snrnaseq xenograft alignment, counting and variant calling. Please refer to [README](/xenograft/snrna_seq/README.md).

### License
This project is licensed under the MIT License.

### Acknowledgments
We would like to acknowledge the developers of Giotto and other related packages that have been utilized in this repository.

Please refer to the individual script files for more detailed information about each script's functionality, input requirements, and output formats.
