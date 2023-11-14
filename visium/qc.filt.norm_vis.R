# Scripts for reading in and performing quality control of Visium data
# Modified from  https://giottosuite.readthedocs.io/en/latest/subsections/datasets/mouse_visium_brain.html - Thank you Giotto Team!

# libraries
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/giotto_env.R")
library(Giotto)
library(GiottoData)
library(glmGamPoi)
library(readxl)
library(corrplot)
library(edgeR)
library(ssizeRNA)
library(ggplot2)
library(magick)
library(MASS)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(readxl)
library(stringr)
library(Seurat)

############################################################################################
#### Inputs
############################################################################################

# run_name = "V52Y16-079-A1"
# results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/V52Y16-079-A1'
# data_path = '/Users/zacc/USyd/spatial_transcriptomics/data/SANPIN_VisiumFFPE_Cytassist_results/221011/V52Y16-079-A1/VISIUM/V52Y16-079-A1/outs'
# ctr_image = 'tissue_lowres_image.png' # name of images file within '/spatial' directory of visium folder

# run_name = "V52Y16-079-B1"
# results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/V52Y16-079-B1'
# data_path = '/Users/zacc/USyd/spatial_transcriptomics/data/SANPIN_VisiumFFPE_Cytassist_results/221011/V52Y16-079-B1/VISIUM/V52Y16-079-B1/outs'
# ctr_image = 'tissue_lowres_image.png' # name of images file within '/spatial' directory of visium folder
# 

run_name = "S11W1"
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/S11W1'
data_path = '/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/GIMR_GWCCG_220131_SANPIN_10X_Visium.analyses/analyses/S11W1/spaceranger_count_cytassist/S11W1/outs'
ctr_image = 'tissue_lowres_image.png' # name of images file within '/spatial' directory of visium folder


## thresholds
expression_threshold = 1
feat_det_in_min_cells = 50
min_det_feats_per_cell = 500

# source
#source("/Users/zacc/github_repo/Giotto/R/general_help.R")
#source("/Users/zacc/github_repo/Giotto/R/utilities.R")

############################################################################################
###### Part 1: Create Giotto Object 
############################################################################################
setwd(results_folder)
# define Giotto instructions
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = my_python_path)
## Create Giotto object
visium_brain = createGiottoVisiumObject(visium_dir = data_path,
                                        expr_data = 'filter',
                                        png_name = ctr_image,
                                        gene_column_index = 2,
                                        instructions = instrs)
dim(pDataDT(visium_brain))

## show plot
my_spatPlot <- spatPlot2D(gobject = visium_brain, cell_color = 'in_tissue', point_size = 1,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image', save_param = list(save_name = "qc_1_2d.intissue"))

############################################################################################
###### Part 2: Manually define polygon coordinates 
############################################################################################
# This step is highly specific to the sample and needs to be performed interactively.
# The goal is to define polygon coordinates that;
# i) classify per tissue spots, if multiple tissues run on single visium array
# ii) classify spots for removal due to technical issues eg. tissue overlaps.

# #uncomment and run interactively
# my_polygon_coordinates_tissue <- plotInteractivePolygons(my_spatPlot)
# my_polygon_coordinates_remove <- plotInteractivePolygons(my_spatPlot)
# save(my_polygon_coordinates_tissue,my_polygon_coordinates_remove,file=paste0(run_name,"interact_polygons.rds"))
load(paste0(run_name,"interact_polygons.rds"))

## Transform the coordinates into a Giotto polygon object
my_giotto_polygons_tissue <- createGiottoPolygonsFromDfr(my_polygon_coordinates_tissue, name = 'tissue')
my_giotto_polygons_remove <- createGiottoPolygonsFromDfr(my_polygon_coordinates_remove, name = 'remove')

## Add the polygons to the Giotto object
visium_brain <- addGiottoPolygons(gobject = visium_brain, gpolygons = list(my_giotto_polygons_tissue,my_giotto_polygons_remove))

## Add corresponding polygon IDs to cell metadata
visium_brain <- addPolygonCells(visium_brain,polygon_name = 'tissue')
visium_brain <- addPolygonCells(visium_brain,polygon_name = 'remove')

## Remove spots identified for removal & keep spots in tissue
metadata = pDataDT(visium_brain)
remove_barcodes = metadata[remove == "no_polygon"]$cell_ID
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_brain = subsetGiotto(visium_brain, cell_ids = remove_barcodes)
visium_brain = subsetGiotto(visium_brain, cell_ids = in_tissue_barcodes)
dim(pDataDT(visium_brain))

## show plot
spatPlot2D(gobject = visium_brain, cell_color = 'in_tissue', point_size = 1,
                          cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
                          show_image = T, image_name = 'image', save_param = list(save_name = "qc_2_2d.keep.intissue"))

############################################################################################
###### Part 3: Filter spots 
############################################################################################
## filter
visium_brain <- filterGiotto(gobject = visium_brain,
                          expression_threshold = expression_threshold,
                          feat_det_in_min_cells = feat_det_in_min_cells,
                          min_det_feats_per_cell = min_det_feats_per_cell,
                          expression_values = c('raw'),
                          verbose = T)
dim(pDataDT(visium_brain))
dim(fDataDT(visium_brain))

## show plot
spatPlot2D(gobject = visium_brain, cell_color = 'in_tissue', point_size = 1,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image', save_param = list(save_name = "qc_3_2d.filter.intissue"))

## add raw gene & cell statistics
visium_brain <- addStatistics(gobject = visium_brain, expression_values = "raw")

## plot raw expression values
spatPlot2D(gobject = visium_brain, point_alpha = 0.8,point_size = 1.5,
           cell_color = 'nr_feats', color_as_factor = F,
           show_image = T, image_name = 'image', save_param = list(save_name = "qc_4_2d.raw.exp"))


############################################################################################
#### Part 4: Normalization using Seurat SCTtransform 
############################################################################################
# Please refer to "ASAP-SpatialTranscriptomics/geomx/normalization_review_vis.R" for methods in selecting SCTtransform
# read in data using seurat - default is filtered
brain <- Load10X_Spatial(data_path)

# filter for QC'd spots and features
brain <- subset(brain, cells = pDataDT(visium_brain)$cell_ID , features = fDataDT(visium_brain)$feat_ID )

# Normalization with sctransform residuals for all genes. Note; not using v2 regularization as was found to 
# have a higher correlation with normalized data and number of UMIs and hence not used.
brain <- SCTransform(brain, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

# Standard log normalization for comparison
brain <- NormalizeData(brain, verbose = FALSE, assay = "Spatial")

# plots
# Computes the correlation of the log normalized data and sctransform residuals with the
# number of UMIs
brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
p1 <- GroupCorrelationPlot(brain, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(-1,1)
p2 <- GroupCorrelationPlot(brain, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(-1,1)

ggsave("qc_4_scatter.normeffect.png",
       ggarrange(p1,p2, ncol=2,nrow=2) + bgcolor("white"),
       device = "png")

############################################################################################
#### Part 5: Add metadata
############################################################################################
# extract giotto metadata
metadata = as.data.frame(pDataDT(visium_brain))
row.names(metadata) <- metadata$cell_ID
metadata <- metadata[,!colnames(metadata) %in% c("cell_ID")]

# add to seurat object
brain <- AddMetaData(brain,
                       metadata = metadata)

# save data
save(brain, file = paste0(run_name,"_qfn.seurat.Rdata"))


############################################################################################
###### Part X: In Development
############################################################################################
# ## Implementing nuclei counts using VistoSeq - https://lmweber.org/OSTA-book/image-segmentation-visium.html
# ## Implementing QC metric cutoffs defined by cell (nuclei) counts - https://lmweber.org/OSTA-book/quality-control.html
# library(scater)
# library(SummarizedExperiment)
# 
# # read in counts
# fnm <- file.path(data_path, "filtered_feature_bc_matrix")
# sce <- DropletUtils::read10xCounts(fnm,col.names = TRUE)
# 
# # read in image data
# img <- readImgData(
#   path = file.path(dir, "spatial"),
#   sample_id = run_name)
# 
# # read in spatial coordinates
# fnm <- file.path(data_path, "spatial", "tissue_positions.csv")
# xyz <- read.csv(fnm, header = TRUE)
# 
# # subset spatial coordinates by filtered spots
# row.names(xyz) <- xyz$barcode
# xyz <- xyz[sce@colData$Barcode,]
# 
# # construct observation & feature metadata
# rd <- S4Vectors::DataFrame(
#   symbol = rowData(sce)$Symbol)
# 
# # construct 'SpatialExperiment'
# spe <- SpatialExperiment(
#   assays = list(counts = assay(sce)),
#   rowData = rd, 
#   colData = DataFrame(xyz), 
#   spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
#   imgData = img,
#   sample_id = run_name)
# 
# # identify mitochondrial genes
# is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$symbol)
# table(is_mito)
# 
# # calculate per-spot QC metrics and store in colData
# spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# head(colData(spe))
# 
# 
# # add metadata to Giotto object
# addCellMetadata(
#   gobject
# )