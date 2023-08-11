# libs
library(Giotto)
library(GiottoData)
library(readxl)
library(corrplot)
library(edgeR)
library(ssizeRNA)
library(ggplot2)
library(magick)
library(MASS)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(readxl)
library(stringr)
source("/Users/zacc/github_repo/Giotto/R/general_help.R")
source("/Users/zacc/github_repo/Giotto/R/utilities.R")

############################################################################################
#### Inputs
############################################################################################
# path to results folder
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/cell_type_decon'
setwd(results_folder)
my_python_path = NULL # alternatively, "/local/python/path/python" if desired.
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = my_python_path)

############################################################################################
###### Part 1: Create Giotto Visium Object for each Visium experiment and join
############################################################################################
# A1 sample #
# provide path to visium folder
data_path = '/Users/zacc/USyd/spatial_transcriptomics/data/SANPIN_VisiumFFPE_Cytassist_results/221011/V52Y16-079-A1/VISIUM/V52Y16-079-A1/outs'
# name of images file within '/spatial' directory of visium folder
ctr_image = 'tissue_lowres_image.png'
## directly from visium folder
visium_brain = createGiottoVisiumObject(visium_dir = data_path,
                                        expr_data = 'raw',
                                        png_name = ctr_image,
                                        gene_column_index = 2,
                                        instructions = instrs)

# add polygon coordinates of tissue and features
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/giotto/A1/poly_coord_geomxA1SN.R")
my_polygon_coordinates <- poly_coord_geomxA1SN
coord1 <- my_polygon_coordinates[my_polygon_coordinates$name == "SN",]
coord_tissue_overlap <-  my_polygon_coordinates[my_polygon_coordinates$name == "tissue_overlap",] # redundant but needed for join
# We must transform the data.table or data.frame with coordinates into a Giotto polygon object
my_giotto_polygons_SN <- createGiottoPolygonsFromDfr(coord1, name = 'SN')
my_giotto_polygons_tissueoverlap <- createGiottoPolygonsFromDfr(coord_tissue_overlap, name = 'tissue_overlap') # redundant but needed for join
## then, add the polygons to the Giotto object
visium_brain <- addGiottoPolygons(gobject = visium_brain,gpolygons = list(my_giotto_polygons_SN,my_giotto_polygons_tissueoverlap))
# Select polygon cells
visium_brain <- addPolygonCells(visium_brain, polygon_name = 'SN')
visium_brain <- addPolygonCells(visium_brain, polygon_name = 'tissue_overlap')


# B1 sample #
# provide path to visium folder
data_path = '/Users/zacc/USyd/spatial_transcriptomics/data/SANPIN_VisiumFFPE_Cytassist_results/221011/V52Y16-079-B1/VISIUM/V52Y16-079-B1/outs'
# name of images file within '/spatial' directory of visium folder
pd_image = 'tissue_lowres_image.png'
## directly from visium folder
visium_brain_b1 = createGiottoVisiumObject(visium_dir = data_path,
                                           expr_data = 'raw',
                                           png_name = pd_image,
                                           gene_column_index = 2,
                                           instructions = instrs)

# add polygon coordinates of tissue and features
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/giotto/B1/poly_coord_geomxB1SN.R")
## We must transform the data.table or data.frame with coordinates into a Giotto polygon object
coord_SN <- my_polygon_coordinates
coord_tissue_overlap <- my_polygon_coordinates2
## We must transform the data.table or data.frame with coordinates into a Giotto polygon object
my_giotto_polygons_SN <- createGiottoPolygonsFromDfr(coord_SN , name = 'SN')
my_giotto_polygons_tissueoverlap <- createGiottoPolygonsFromDfr(coord_tissue_overlap, name = 'tissue_overlap')
## Then, add the polygons to the Giotto object
visium_brain_b1 <- addGiottoPolygons(gobject = visium_brain_b1,gpolygons = list(my_giotto_polygons_SN, my_giotto_polygons_tissueoverlap))
## Select polygon cells
visium_brain_b1 <- addPolygonCells(visium_brain_b1, polygon_name = 'SN')
visium_brain_b1 <- addPolygonCells(visium_brain_b1, polygon_name = 'tissue_overlap')


# join giotto objects
# joining with x_shift has the advantage that you can join both 2D and 3D data
# x_padding determines how much distance is between each dataset
# if x_shift = NULL, then the total shift will be guessed from the giotto image
testcombo = joinGiottoObjects(gobject_list = list(visium_brain, visium_brain_b1),
                              gobject_names = c('CTR', 'PD'),
                              join_method = 'shift', x_padding = 100)

# join info is stored in this slot
# simple list for now
testcombo@join_info

# check joined Giotto object
fDataDT(testcombo)
pDataDT(testcombo)
showGiottoImageNames(testcombo)
showGiottoSpatLocs(testcombo)
showGiottoExpression(testcombo)

# # this plots all the images by list_ID
# spatPlot2D(gobject = testcombo, cell_color = 'in_tissue',
#            show_image = T, image_name = c("CTR-image","PD-image"),
#            group_by = 'list_ID', point_alpha = 0.5,
#            point_size = 0.5, cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
#            save_param = list(save_name = "1a_plot"))
# 

############################################################################################
###### Part 2: Process Giotto Objects
############################################################################################
# subset on in-tissue spots
metadata = pDataDT(testcombo)
in_tissue_barcodes = metadata[in_tissue == 1 & tissue_overlap == 'no_polygon']$cell_ID
testcombo = subsetGiotto(testcombo, cell_ids = in_tissue_barcodes)

## filter
testcombo <- filterGiotto(gobject = testcombo,
                          expression_threshold = 1,
                          feat_det_in_min_cells = 50,
                          min_det_feats_per_cell = 500,
                          expression_values = c('raw'),
                          verbose = T)

## normalize
testcombo <- normalizeGiotto(gobject = testcombo, scalefactor = 6000)

## add gene & cell statistics
testcombo <- addStatistics(gobject = testcombo, expression_values = 'raw')
fDataDT(testcombo)
pDataDT(testcombo)

# ## visualize
# spatPlot2D(gobject = testcombo, group_by = 'list_ID', cell_color = 'nr_feats', color_as_factor = F, 
#            point_size = 0.75, save_param = list(save_name = "2a_plot"))

# split SN and VTA tissue sections
metadata = pDataDT(testcombo)
in_SN = metadata[in_tissue == 1 & SN %in% c("polygon 1","SN")]$cell_ID
in_VTA = metadata[in_tissue == 1 & SN %in% c("no_polygon")]$cell_ID
SN_vis = subsetGiotto(testcombo, cell_ids = in_SN)
LC_vis = subsetGiotto(testcombo, cell_ids = in_VTA)

## visualize
# spatPlot2D(gobject = SN_vis, group_by = 'list_ID', cell_color = 'nr_feats', color_as_factor = F, 
#            point_size = 1, save_param = list(save_name = "2b_plot"))
# spatPlot2D(gobject = LC_vis, group_by = 'list_ID', cell_color = 'nr_feats', color_as_factor = F, 
#            point_size = 1, save_param = list(save_name = "2c_plot"))
# 

##### Part C. Dimension Reduction and clustering
# Substantia Nigra
## highly variable features / genes (HVF)
SN_vis <- calculateHVF(gobject = SN_vis, save_plot = TRUE)
## run PCA on expression values (default)
gene_metadata = fDataDT(SN_vis)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID
## run PCA on expression values (default)
SN_vis <- runPCA(gobject = SN_vis, feats_to_use = featgenes)
# plots
screePlot(SN_vis, ncp = 30)
dimPlot2D(gobject = SN_vis,dim_reduction_to_use = "pca")
## run UMAP and tSNE on PCA space (default)
SN_vis <- runUMAP(SN_vis, dimensions_to_use = 1:10)
plotUMAP(gobject = SN_vis)
SN_vis <- runtSNE(SN_vis, dimensions_to_use = 1:10)
plotTSNE(gobject = SN_vis)
## sNN network (default)
SN_vis <- createNearestNetwork(gobject = SN_vis, dimensions_to_use = 1:10, k = 15)
## Leiden clustering
SN_vis <- doLeidenCluster(gobject = SN_vis, resolution = 0.4, n_iterations = 1000)

# Locus Ceoreluos
## highly variable features / genes (HVF)
LC_vis <- calculateHVF(gobject = LC_vis, save_plot = TRUE)
## run PCA on expression values (default)
gene_metadata = fDataDT(LC_vis)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID
## run PCA on expression values (default)
LC_vis <- runPCA(gobject = LC_vis, feats_to_use = featgenes)
# plots
screePlot(LC_vis, ncp = 30)
dimPlot2D(gobject = LC_vis,dim_reduction_to_use = "pca")
## run UMAP and tSNE on PCA space (default)
LC_vis <- runUMAP(LC_vis, dimensions_to_use = 1:10)
plotUMAP(gobject = LC_vis)
LC_vis <- runtSNE(LC_vis, dimensions_to_use = 1:10)
plotTSNE(gobject = LC_vis)
## sNN network (default)
LC_vis <- createNearestNetwork(gobject = LC_vis, dimensions_to_use = 1:10, k = 15)
## Leiden clustering
LC_vis <- doLeidenCluster(gobject = LC_vis, resolution = 0.4, n_iterations = 1000)


# plots
plotUMAP(gobject = SN_vis, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)
plotUMAP(gobject = SN_vis, cell_color = 'list_ID', point_size = 2.5)
spatDimPlot(gobject = SN_vis, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)


plotUMAP(gobject = LC_vis, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)
plotUMAP(gobject = LC_vis, cell_color = 'list_ID', point_size = 2.5)
spatDimPlot(gobject = LC_vis, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)



############################################################################################
###### Part 3: Cell type enrichment
############################################################################################
### Cell enrichment of DA/NE neurons, Glia subtypes and T-cells
# load Neuronal and Glia marker genes
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/giotto/compare_visium/gene_marker_list.cellenrich.R")

# test each cell enrichment marker selection method
gene_marker_list <- gene_marker_list_topfdrwithincell

# create PAGE matrix
PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = names(gene_marker_list),sign_list = gene_marker_list)

# Substantia Nigra
# runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
SN_vis = runPAGEEnrich(gobject = SN_vis, sign_matrix = PAGE_matrix_1)
# visulaize 
cell_types_PAGE = colnames(PAGE_matrix_1)
d <- 1:length(cell_types_PAGE)
split_list <- split(d, ceiling(seq_along(d)/4))
for (i in 1:length(split_list)){
  print(unlist(split_list[i]))
  spatCellPlot2D(gobject = SN_vis ,
                 spat_enr_names = 'PAGE',
                 cell_annotation_values = cell_types_PAGE[unlist(split_list[i])],
                 cow_n_col = 2,coord_fix_ratio = 1, point_size = 0.5, show_legend = T,
                 show_image = F, save_plot = T)
}

# Locus Ceoreluos
# runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
LC_vis = runPAGEEnrich(gobject = LC_vis, sign_matrix = PAGE_matrix_1)
# visulaize 
cell_types_PAGE = colnames(PAGE_matrix_1)
d <- 1:length(cell_types_PAGE)
split_list <- split(d, ceiling(seq_along(d)/4))
for (i in 1:length(split_list)){
  print(unlist(split_list[i]))
  spatCellPlot2D(gobject = LC_vis ,
                 spat_enr_names = 'PAGE',
                 cell_annotation_values = cell_types_PAGE[unlist(split_list[i])],
                 cow_n_col = 2,coord_fix_ratio = 1, point_size = 0.5, show_legend = T,
                 show_image = F, save_plot = T)
}



#######################################
### Marker genes of neuron subtypes ###
#######################################
GenewiseCounts <- as.matrix(get_expression_values(SN_vis, output = "matrix",values = "scaled"))
genes_interest <- c("GAL","NPY","ADRA1A","ADRA2A","ALDH1A1","CALB1","TH")

#row.names(GenewiseCounts )[grep("ADOR",row.names(GenewiseCounts ))]

# visulaize 
d <- 1:length(genes_interest)
split_list <- split(d, ceiling(seq_along(d)/4))
for (i in 1:length(split_list)){
  print(unlist(split_list[i]))
  spatFeatPlot2D(SN_vis, expression_values = 'normalized',
                 feats = genes_interest[unlist(split_list[i])], point_size = 0.5)
}

for (i in 1:length(split_list)){
  print(unlist(split_list[i]))
  spatFeatPlot2D(LC_vis, expression_values = 'normalized',
                 feats = genes_interest[unlist(split_list[i])], point_size = 0.5)
}
