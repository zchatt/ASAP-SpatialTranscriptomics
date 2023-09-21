# scripts for running RCTD and CSIDE method to identify ctDEG

library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(plyr)
library(data.table)

############################################################################################
#### Inputs
############################################################################################

analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis"
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/geomx_sep2023_norm.gx.Rdata"
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/RCTD_results'
tissue_analysis = "Midbrain" # tissue to analyse

# single-cell
cell_ranger_data = "/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger" # single-cell cell ranger output
filenames <- list.files("/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger", pattern="*UMAP.tsv", full.names=TRUE) # read in metadata / umap of cell IDs
cell_id <- do.call(rbind, lapply(filenames,fread))
#cell_id <- read.table(file = '/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger/da_UMAP.tsv', sep = '\t', header = TRUE) # read in metadata / umap of cell IDs
cell_id <- cell_id[!cell_id$NAME == "TYPE",]

############################################################################################
###### Part 1: RCTD & CSIDE analysis
############################################################################################
# setwd
setwd(results_folder)

### Data Preprocessing and running RCTD
if(!file.exists(file.path(results_folder,'myRCTDde_visium.rds'))) {
  
  ## spatial data
  # load geomx normalised data
  load(rdata)
  # select tissue
  norm.quantile
  run_index =  brain.merge@meta.data$tissue == tissue_analysis
  # load in counts matrix
  counts <- brain.merge@assays$SCT@counts[,run_index]  # load in counts matrix
  # load in coordinates
  offset = 500
  kth_image <- brain.merge@images
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
  coords = combined_coordinates[run_index,c("imagecol","imagerow")]
  #coords <- brain.merge@meta.data[run_index,]
  #coords <- coords[,c("array_col","array_row")]
  nUMI <- colSums(counts) # In this case, total counts per spot
  puck <- SpatialRNA(coords, counts, nUMI)
  barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
  plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 
  
  # #single-cell data
  sc_obj <- Read10X(cell_ranger_data) # single-cell cell ranger output
  # subset for cel-types
  sc_obj <- sc_obj[,cell_id$NAME]
  cell_types <- setNames(cell_id$Cell_Type, cell_id$NAME)
  #cell_types <- as.factor(cell_id$Cell_Type) # convert to factor data type
  cell_types <- as.factor(cell_types)
  nUMI <- colSums(sc_obj)
  
  ## create reference and run RCTD
  reference <- Reference(sc_obj, cell_types, nUMI)
  myRCTD <- create.RCTD(puck, reference, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD,file.path(results_folder,'myRCTDde_visium.rds'))
} else {
  myRCTD <- readRDS(file.path(results_folder,'myRCTDde_visium.rds'))
}

### Exploring the full mode results
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
print(head(norm_weights)) 

p1 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'SOX6_AGTR1'], ylimit = c(0,0.5), 
                           title ='plot of SOX6_AGTR1 weights', size = 0.5) 

ggsave("cside_1_sox6agtr1.png", p1, device = "png")

p2 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'CALB1_GEM'], ylimit = c(0,0.5), 
                           title ='plot of CALB1_GEM weights', size = 0.5) 

ggsave("cside_2_calb1gem.png", p2, device = "png")

## Running CSIDE - single
myRCTD@config$max_cores <- 2
myRCTD@config$doublet_mode = 'full'

explanatory.variable <- as.numeric(as.factor(mapvalues(brain.merge@meta.data$run_name,run_names,list_ID))) - 1
names(explanatory.variable) <-  barcodes

myRCTD <- run.CSIDE.single(myRCTD,
                           explanatory.variable,
                           cell_type_threshold = 10,
                           fdr = 0.05, 
                           doublet_mode = FALSE,
                           weight_threshold = 0.1) 

saveRDS(myRCTD,file.path(results_folder,'myRCTDde_visium.rds'))

make_all_de_plots(myRCTD,results_folder)

## CSIDE results
all_gene_list <- myRCTD@de_results$all_gene_list
sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val < 0.05,])
sig_gene_list <- lapply(all_gene_list, function(x) x[order(x$p_val),])


############################################################################################
###### Part X: In development
############################################################################################

### plotting cell-type distributions dorsal to ventral
norm_weights

x_coord <- myRCTD@spatialRNA@coords[explanatory.variable == 0,]
ct_data <- norm_weights[row.names(x_coord),]
# create bin
x_coord$x

library(ggplot2)
library(ggridges)
theme_set(theme_minimal()) 

ggplot(
  lincoln_weather, 
  aes(x = `Mean Temperature [F]`, y = `Month`, fill = stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Temp. [F]", option = "C") +
  labs(title = 'Temperatures in Lincoln NE')             




# ## Running CSIDE - regions/ multiple contrasts
# region_left <- barcodes[which(myRCTD@spatialRNA@coords$x < quantile(myRCTD@spatialRNA@coords$x, 1/3))]
# region_right <- barcodes[which(myRCTD@spatialRNA@coords$x > quantile(myRCTD@spatialRNA@coords$x, 2/3))]
# region_middle <- setdiff(barcodes, union(region_left, region_right))
# region_list <- list(region_left, region_right, region_middle)
# 
# class_num <- rep(0, length(barcodes)); names(class_num) <- barcodes
# class_num[region_middle] <- 1; class_num[region_right] <- 2 
# 
# p3  <- plot_class(myRCTD@spatialRNA, barcodes, factor(class_num),  title ='plot of regions')
# ggsave("cside_3_regionsplit.png", p3, device = "png")

# 
# ## build design matrix
# regions_list <- list()
# contrast <- mapvalues(brain.merge@meta.data$run_name,run_names,list_ID)
# for (i in 1:length(unique(contrast))) {
#   regions_list[[i]] <- barcodes[contrast == unique(contrast)[i]]
# }
# regions_list[[i+1]] <- NA
# 
# X <- build.designmatrix.regions(myRCTD,regions_list)
# barcodes <- rownames(X)
# myRCTD <- run.CSIDE(myRCTD, X, barcodes,cell_types = cell_types,
#                     cell_type_threshold = 10,
#                     fdr = 0.05, 
#                     doublet_mode = FALSE, 
#                     weight_threshold = 0.1,
#                     test_mode = 'categorical') 
# 
# saveRDS(myRCTD,file.path(results_folder,'myRCTDde_visium_SN_regions.rds'))
# 
# ## CSIDE results
# all_gene_list <- myRCTD@de_results$all_gene_list
# sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val_best < 0.05,])
# sig_gene_list <- lapply(all_gene_list, function(x) x[order(x$p_val_best),])
