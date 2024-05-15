# scripts for running RCTD and CSIDE method to identify ctDEG

library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(plyr)
library(data.table)
library(msigdbr)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(EnhancedVolcano)
library(viridis)
library(rstatix)
library(ComplexHeatmap)
library(readxl)
source("/Users/zacc/github_repo/spacexr/R/CSIDE_plots.R")
register(SerialParam())

############################################################################################
#### Inputs
############################################################################################
results_folder = "/sers/zacc/USyd/spatial_transcriptomics/analysis/visium/vis_130524"
setwd(results_folder)

### ST datasets ###

# i) visium original (n=2)
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/SN_210823/SN_210823_merged.seurat.Rdata")
st_obj <- subset(x = brain.merge, subset = tissue == "SN")


# ii) visium 051324 - Sandy Controls
working_dir = "/sers/zacc/USyd/spatial_transcriptomics/analysis/visium/vis_130524"
seurat_obj <- LoadSeuratRds("/Volumes/research-data/PRJ-ASAPbioin/Bioinformics/datasets/seurat_objects/midbrain-controls-10x/controls_midbrain.harmony.rds")
setwd(working_dir)


#### single-cell data ###
## Kamath dataset ## 

# cell_ranger_data = "/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger" # single-cell cell ranger output
# filenames <- list.files("/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger", pattern="*UMAP.tsv", full.names=TRUE) # read in metadata / umap of cell IDs
# cell_id <- do.call(rbind, lapply(filenames,fread))
# #cell_id <- read.table(file = '/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger/da_UMAP.tsv', sep = '\t', header = TRUE) # read in metadata / umap of cell IDs
# cell_id <- cell_id[!cell_id$NAME == "TYPE",]
# # #single-cell data
# sc_obj <- Read10X(cell_ranger_data) # single-cell cell ranger output
# # subset for cell-types
# sc_obj <- sc_obj[,cell_id$NAME]
# cell_types <- setNames(cell_id$Cell_Type, cell_id$NAME)
# #cell_types <- as.factor(cell_id$Cell_Type) # convert to factor data type
# cell_types <- as.factor(cell_types)
# nUMI <- colSums(sc_obj)


## merged dataset ## 
load("/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/merged_kam.sil.web_seurat.Rdata")

# pre-processing
# i) taking excitatory, inhibitory and non-neurons from Kamath and DA neuron populations from Siletti re-code
toMatch <- c("Ex_","Inh_","Astro_","Endo_","Ependyma_","Macro_","MG_","Olig_","OPC_",
             names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Siletti"])))
group1 <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
sc_obj <- merge.combined.sct[,group1]

# ii) remove cell-types with majority of cells within the Periaqueductal gray matter - Abaurre reference
sc_obj <- sc_obj[,!sc_obj@meta.data$cell_type_merge %in% c("CALB1_NPW_ONECUT1")]

# iii) joining layers
DefaultAssay(sc_obj) <- 'RNA'
sc_obj <- JoinLayers(sc_obj)


# # remove cell-types with <25 cells
# cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
# sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# # get count data
# DefaultAssay(sc_obj) <- "SCT"
# counts_sc <- sc_obj[["SCT"]]$counts
# # set cell types and quant nUMI
# cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
# cell_types <- as.factor(cell_types)
# nUMIsc <- colSums(counts_sc)


############################################################################################
###### Part 1:  Format spatial data
############################################################################################

# load in counts matrix
counts <- st_obj@assays$SCT@counts
# load in coordinates
offset = 500
kth_image <- st_obj@images
num_elements <- length(kth_image)
coordinates_list <- vector("list", length = num_elements)
for (i in seq_along(kth_image)) {
  if ("coordinates" %in% slotNames(kth_image[[i]])) {
    coordinates_list[[i]] <- kth_image[[i]]@coordinates
    if (i > 1){
      coordinates_list[[i]]["imagecol"] <-  coordinates_list[[i]]["imagecol"] + max(coordinates_list[[i-1]]["imagecol"]) + offset
    }
  }
}
combined_coordinates <- do.call(rbind, coordinates_list)
coords = combined_coordinates[,c("imagecol","imagerow")]

# inputs for RCTD
nUMI <- colSums(counts) # In this case, total counts per spot
puck <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 


############################################################################################
###### Part 2:  Format reference (sc data)
############################################################################################

# remove cell-types with < 25 cells
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]

# get count data and restrict to only genes within ST
counts_sc <- sc_obj[["RNA"]]$counts
counts_sc <- counts_sc[row.names(counts_sc) %in% row.names(puck@counts),]

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)

## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)


############################################################################################
###### Part 3: run RCTD to obtain cell proportions
############################################################################################

# run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add data to ST metadata
st_obj@meta.data <- cbind(st_obj@meta.data,normalize_weights(myRCTD@results$weights))

# plot meta
ratio = 1.9
meta_to_plot <- c("CALB1_CBLN4_PAX5","CALB1_CRYM","CALB1_CRYM_CALCR","CALB1_CRYM_CCDC68",
                  "CALB1_NEUROD6_PPP1R17", "CALB1_PRLHR_TRHR","CALB1_SEMA3D_RSPO3","CALB1_VIP_NPPC",
                  "GAD2_CALCRL_KCNK13","GAD2_EBF2_NPSR1","SOX6_AGTR1_NOX4", "SOX6_GFRA2_TBC1D8B","SOX6_SMOC1_LPL")

for (i in 1:length(meta_to_plot)){
  plot2 <- SpatialFeaturePlot(st_obj, features = meta_to_plot[i], images = "slice1") + theme(legend.position = "right", aspect.ratio = ratio)
  ggsave(plot2,file=paste0("spatplot_rctd_",meta_to_plot[i],".pdf"))
}

# plot gene expression


