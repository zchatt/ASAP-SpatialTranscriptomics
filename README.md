## ASAP-SpatialTranscriptomics

## This repository contains analysis scripts used by The University of Sydney collaborators in the analysis of spatial transcriptomics as part of [Aligning Science Accross Parkinson's](https://parkinsonsroadmap.org/#).

 
## Visium_giotto.R: Data Processing and Analysis 
### Overview
This script utilize [Giotto](https://giottosuite.readthedocs.io/en/latest/index.html), a popular R package for analyzing spatial transcriptomic data. The scripts are designed to perform various analyses on Visium data, including data preprocessing, dimension reduction, clustering, and cell-type enrichment.

### The script takes the following inputs:

Visium data matrix: A matrix containing gene expression values for spatially resolved spots.
Brain tissue image: An image representing the spatial coordinates of the spots in the Visium dataset.

### The script performs the following steps:

1. Creating a Giotto image: The script generates a Giotto image object that combines the Visium data matrix and brain tissue image.

2. Manual annotation: Users can manually annotate spots of concern, enabling the removal of unwanted spots from the analysis.

3. RNA filtering: The script filters RNA based on expression thresholds and identifies RNA occurring in a minimum number of spots.

4. Spot filtering: The script filters spots based on the minimum number of RNA features detected.

5. Dimension reduction: Principal Component Analysis (PCA) is performed to reduce the dimensionality of the data.

6. Clustering: The reduced data is clustered using t-distributed Stochastic Neighbor Embedding (tSNE) and Uniform Manifold Approximation and Projection (UMAP) algorithms.

7. Cell-type enrichment: The script performs cell-type enrichment analysis using gene markers for the following cell types: "Astrocytes-1", "Astrocytes-2", "CALB1_CALCR", "CALB1_CRYM_CCDC68", "CALB1_GEM", "CALB1_PPP1R17", "CALB1_RBP4", "CALB1_TRHR", "CD4_NaiveLike", "CD8_EarlyActiv", "CD8_EffectorMemory", "CD8_NaiveLike", "CD8_Tex", "CD8_Tpex", "Endothelial", "GABA neurons", "Microglia", "NE", "ODC-1", "ODC-2", "ODC-3", "OPC", "SOX6_AGTR1", "SOX6_DDT", "SOX6_GFRA2", "SOX6_PART1", "Tfh", "Th1", and "Treg".

### Usage
To run the script, follow these steps:

1. Install the required dependencies specified in the requirements.txt file.

2. Prepare your Visium data matrix and brain tissue image and place them in the appropriate directories.

3. Open the R script and modify the file paths and parameters as needed.

4. Run the script using an R environment or an integrated development environment (IDE).

5. Review the generated results and output files for further analysis.

## Geomx_giotto.R: GeoMx Data Processing and Analysis
### Overview
This scrips is designed to analyze GeoMx spatial transcriptomics data produced by NanoString consisting spatial "Areas of Interest" (AOIs).

The script performs similar steps to , with some modifications to accommodate the GeoMx data:

1. Creating a Giotto image: The script generates a Giotto image object that combines the GeoMx data matrix and relevant metadata.

2. Manual annotation: Users can manually annotate AOIs of concern, enabling the removal of unwanted AOIs from the analysis.

3. RNA filtering: The script filters RNA based on expression thresholds and identifies RNA occurring in a minimum number of AOIs.

4. AOI filtering: The script filters AOIs based on the minimum number of RNA features detected.

5. Dimension reduction: Principal Component Analysis (PCA) is performed to reduce the dimensionality of the data.

6. Clustering: The reduced data is clustered using t-distributed Stochastic Neighbor Embedding (tSNE) and Uniform Manifold Approximation and Projection (UMAP) algorithms.

7. Cell-type enrichment: The script performs cell-type enrichment analysis using gene markers for specific cell types, as described above.

### Usage
To run the script, follow these steps:

1. Install the required dependencies specified in the requirements.txt file.

2. Prepare your GeoMx data matrix and relevant metadata, and place them in the appropriate directories.

3. Open the R script and modify the file paths and parameters as needed.

4. Run the script using an R environment or an integrated development environment (IDE).

5. Review the generated results and output files for further analysis.

### License
This project is licensed under the MIT License.

### Acknowledgments
We would like to acknowledge the developers of Giotto and other related packages that have been utilized in this repository.

Please refer to the individual script files for more detailed information about each script's functionality, input requirements, and output formats.
