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

# NOTE; taking excitatory, inhibitory and non-neurons from Kamath, NE and 5-HT from Webber, DA neuron populations from Siletti re-code
toMatch <- c("Ex_","Inh_","Astro_","Endo_","Ependyma_","Macro_","MG_","Olig_","OPC_","5HT","NE",
             names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Siletti"])))
group1 <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
sc_obj <- merge.combined.sct[,group1]
# remove cell-types with majority of cells within the Periaqueductal gray matter - Abaurre reference
sc_obj <- sc_obj[,!sc_obj@meta.data$cell_type_merge %in% c("CALB1_NPW_ONECUT1")]

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

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# combine with snRNAseq and plot
dplot1 <- do.call("rbind", res)
dplot <- as.data.frame(siletti_DA_truth)
colnames(dplot) <- c("Cell","Proportion")
dplot$group <- "Siletti re Reference"

dplot <- rbind(dplot,dplot1)
dplot$Proportion <- as.numeric(dplot$Proportion)

g1 <- ggplot(dplot, aes(x = group, y = Proportion, fill= Cell)) + 
  geom_bar(stat = "identity") + xlab("") + theme_bw()

# save data
ggsave("barplots_siletti_TH.CTR.png", g1)

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))
#save(meta_rctd,file = "meta_rctd_Siletti_DAn12.Rdata")



### ii) deconvolute using Kamath
set.seed(123)
# construct ST to deconvolute for SN tissue
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')

# construct references
# NOTE; taking all cells from Kamath
toMatch <- c(names(table(sc_obj@meta.data$cell_type_merge[sc_obj@meta.data$dataset_merge == "Kamath"])))
group1 <- sc_obj@meta.data$cell_type_merge %in% toMatch
sc_obj <- sc_obj[ ,group1]

# remove cell-types with < 25 cells AND/OR cells of interes
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
#cells_keep <- c("CALB1_CALCR","CALB1_CRYM_CCDC68","CALB1_GEM", "CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR",
#                "SOX6_AGTR1","SOX6_DDT","SOX6_GFRA2", "SOX6_PART1")
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# get count data and restrict to only genes within GeoMx
counts_sc <- sc_obj[["RNA"]]$counts
counts_sc <- counts_sc[row.names(counts_sc) %in% row.names(puck@counts),]

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# iii) compare to ground truth
meta_select <- meta[group2,]
roi_interest <- unique(meta_select$ROI)
res <- list()
for (i in 1:length(roi_interest)){
  # select CTR TH+
  norm_weights <- normalize_weights(myRCTD@results$weights)
  norm_weights_sub <- norm_weights[meta_select$Diagnosis == "CTR" & meta_select$segment == "Full ROI" & meta_select$ROI %in% roi_interest[i],]
  
  # stacked barplots of each
  m1 <- colMeans(norm_weights_sub)
  dplot1 <- as.data.frame(cbind(Cell = names(m1),Proportion = m1))
  dplot1$group <- paste0("CTR, Full ROI, ",roi_interest[i])
  res[[i]] <- dplot1
  
}

# combine with snRNAseq and plot
dplot1 <- do.call("rbind", res)
dplot <- as.data.frame(kamath_all_truth)
colnames(dplot) <- c("Cell","Proportion")
dplot$group <- "Kamath Reference"

dplot <- rbind(dplot,dplot1)
dplot$Proportion <- as.numeric(dplot$Proportion)

g1 <- ggplot(dplot, aes(x = group, y = Proportion, fill= Cell)) + 
  geom_bar(stat = "identity",colour="black") + xlab("") + theme_bw()

# save data
ggsave("barplots_kamath_Full.CTR.png", g1)

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))
#save(meta_rctd,file = "meta_rctd_Kamath_n63.Rdata")
#save(meta_rctd,file = "meta_rctd_Kamath_n10.Rdata")




#### iii) deconvolute using Kamath non-DA and Siletti DA
# construct ST to deconvolute for SN tissue
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')

# construct references
# NOTE; taking all cells from Kamath
toMatch <- c(names(table(sc_obj@meta.data$cell_type_merge[sc_obj@meta.data$dataset_merge %in% c("Kamath","Siletti")])))
group1 <- sc_obj@meta.data$cell_type_merge %in% toMatch
sc_obj <- sc_obj[ ,group1]

# remove cell-types with < 25 cells AND/OR cells of interest
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# remove cell-types with majority of cells within the Periaqueductal gray matter AND DA cells from Kamath et al.
sc_obj <- sc_obj[,!sc_obj@meta.data$cell_type_merge %in% c("CALB1_NPW_ONECUT1",
                                                           "CALB1_CALCR","CALB1_CRYM_CCDC68","CALB1_GEM",
                                                           "CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR",
                                                           "SOX6_AGTR1","SOX6_DDT","SOX6_GFRA2", "SOX6_PART1")]

# get count data and restrict to only genes within GeoMx
count1 <- sc_obj[["RNA"]]$counts.1
count1 <- count1[row.names(count1) %in% row.names(puck@counts),]

count2 <- sc_obj[["RNA"]]$counts.2
count2 <- count2[row.names(count2) %in% row.names(puck@counts),]

genes_select <- intersect(row.names(count1),row.names(count2))
counts_sc <- cbind(count1[genes_select,],count2[genes_select,])# merge 
table(colnames(sc_obj) == colnames(counts_sc))

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))
#save(meta_rctd,file = "meta_rctd_Kamath_nonDA.Siletti_DA_n65.Rdata")


#### i) deconvolute using Kamath non-DA and Siletti DA
# construct ST to deconvolute for SN tissue
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')

# construct references
# NOTE; taking all cells from Kamath
toMatch <- c(names(table(sc_obj@meta.data$cell_type_merge[sc_obj@meta.data$dataset_merge %in% c("Kamath","Siletti")])))
group1 <- sc_obj@meta.data$cell_type_merge %in% toMatch
sc_obj <- sc_obj[ ,group1]

# remove cell-types with < 25 cells AND/OR cells of interest
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# remove cell-types with majority of cells within the Periaqueductal gray matter AND DA cells from Kamath et al.
sc_obj <- sc_obj[,!sc_obj@meta.data$cell_type_merge %in% c("CALB1_NPW_ONECUT1",
                                                           "CALB1_CALCR","CALB1_CRYM_CCDC68","CALB1_GEM",
                                                           "CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR",
                                                           "SOX6_AGTR1","SOX6_DDT","SOX6_GFRA2", "SOX6_PART1")]

# get count data and restrict to only genes within GeoMx
count1 <- sc_obj[["RNA"]]$counts.1
count1 <- count1[row.names(count1) %in% row.names(puck@counts),]

count2 <- sc_obj[["RNA"]]$counts.2
count2 <- count2[row.names(count2) %in% row.names(puck@counts),]

genes_select <- intersect(row.names(count1),row.names(count2))
counts_sc <- cbind(count1[genes_select,],count2[genes_select,])# merge 
table(colnames(sc_obj) == colnames(counts_sc))

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))
#save(meta_rctd,file = "meta_rctd_Kamath_nonDA.Siletti_DA_n65.Rdata")



#### iv) deconvolute using Webber LC reference
# construct ST to deconvolute for SN tissue
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')

# construct references
toMatch <- c(names(table(sc_obj@meta.data$cell_type_merge[sc_obj@meta.data$dataset_merge == "Webber"])))
group1 <- sc_obj@meta.data$cell_type_merge %in% toMatch
sc_obj <- sc_obj[ ,group1]

# remove cell-types with < 25 cells
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# get count data and restrict to only genes within GeoMx
counts_sc <- round(10^sc_obj[["RNA"]]$counts,0)
counts_sc <- counts_sc[row.names(counts_sc) %in% row.names(puck@counts),]

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# # iii) compare to ground truth
# meta_select <- meta[group2,]
# roi_interest <- unique(meta_select$ROI)
# res <- list()
# for (i in 1:length(roi_interest)){
#   # select CTR TH+
#   norm_weights <- normalize_weights(myRCTD@results$weights)
#   norm_weights_sub <- norm_weights[meta_select$Diagnosis == "CTR" & meta_select$segment == "TH" & meta_select$ROI %in% roi_interest[i],]
# 
#   # stacked barplots of each
#   m1 <- colMeans(norm_weights_sub)
#   dplot1 <- as.data.frame(cbind(Cell = names(m1),Proportion = m1))
#   dplot1$group <- paste0("CTR, TH, ",roi_interest[i])
#   res[[i]] <- dplot1
# 
# }
# 
# # combine with snRNAseq and plot
# dplot1 <- do.call("rbind", res)
# dplot <- as.data.frame(webber_LC_truth)
# colnames(dplot) <- c("Cell","Proportion")
# dplot$group <- "Webber Reference"
# 
# dplot <- rbind(dplot,dplot1)
# dplot$Proportion <- as.numeric(dplot$Proportion)
# 
# g1 <- ggplot(dplot, aes(x = group, y = Proportion, fill= Cell)) +
#   geom_bar(stat = "identity",colour="black") + xlab("") + theme_bw()
# 
# # save data
# ggsave("barplots_webber_TH.LC.CTR.png", g1)

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))
#save(meta_rctd,file = "meta_rctd_Webber_n10.Rdata")








# NOTE: Below are scripts previously used when just deconvoluting using the Kamath et al. reference.
# There will be handy scripts for plotting and performing cell-type specific DEG analysis.











############################################################################################
###### Part 1: RCTD to obtain cell proportions
############################################################################################
# setwd
setwd(results_folder)

### Data Preprocessing and running RCTD
if(!file.exists(file.path(results_folder,'myRCTDde_gx.rds'))) {
  
  ## spatial data
  # load geomx normalised data
  load(rdata)
  # select samples
  #keep_index <- gxdat_s$Diagnosis == "CTR"
  keep_index <- gxdat_s$Diagnosis != "NTC"
  # load in counts matrix
  #counts <- gxdat_s@assays$GeoMx@counts[,keep_index]  # load in counts matrix
  counts <- gxdat_s@assays$RNA@counts[,keep_index]  # load in counts matrix
  # create metadata matrix 
  meta <- gxdat_s@meta.data[keep_index,]
  # confirm names
  table(colnames(counts) == row.names(meta))
  # load in coordinates
  coords = as.data.frame(gxdat_s@reductions$umap@cell.embeddings[keep_index,]) # we are using UMAP coordinates as dummy variables until we obtain DSP ROI locations
  colnames(coords) <- c("imagerow","imagecol")
  # process counts
  nUMI <- colSums(counts) # In this case, total counts per spot
  puck <- SpatialRNA(coords, counts, nUMI)
  barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
  plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')
  
  ## create reference and run RCTD
  reference <- Reference(counts_sc, cell_types, nUMIsc)
  myRCTD <- create.RCTD(puck, reference, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD,file.path(results_folder,'myRCTDde_gx.rds'))
} else {
  myRCTD <- readRDS(file.path(results_folder,'myRCTDde_gx.rds'))
}


############################################################################################
###### Part 2: RCTD plots and cell-type proportion differences
############################################################################################
### data
barcodes <- colnames(myRCTD@spatialRNA@counts)
norm_weights <- normalize_weights(myRCTD@results$weights)
myRCTD@config$doublet_mode <- "full"
meta <- as.data.frame(gxdat_s@meta.data)
table(row.names(meta) == row.names(norm_weights))
dplot <- cbind(meta,norm_weights)

# colour palettes
Diagnosis_col = c("grey","#00AFBB", "#E7B800","red")
dv200_bin_col = viridis(4)
segment_col = c("grey","black")
Brainregion_col = viridis(4)
area_bin_col = viridis(5)
ROI_col = viridis(7)

### UMAP plots
p1 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'NE'], 
                           title ='plot of SOX6_AGTR1_NOX4 weights', size = 0.5) 

ggsave("cside_1_sox6dtt.png", p1, device = "png")

p2 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'CALB1_RBP4'], ylimit = c(0,0.5), 
                           title ='plot of CALB1_RBP4 weights', size = 0.5) 

ggsave("cside_2_calb1rbp4.png", p2, device = "png")

### Violin plots
# 1. DaN key genes expression
tmp <- t(gxdat_s@assays$RNA@data[c("TH","SLC6A2","KCNJ6"),])
dplot <- cbind(dplot,tmp)
tmp <- dplot[,c("segment","TH","SLC6A2","KCNJ6")]

data_table <- tmp %>% pivot_longer(cols=c("TH","SLC6A2","KCNJ6"),
                                   names_to='Gene',
                                   values_to='Expression')
# make violin plot
v1 <- ggviolin( data_table, x = "Gene", y = "Expression", 
                fill = "segment", color = "segment") 


#  2. DaN cell proportions in TH +/-
dplot$DaN <- rowSums(dplot[,c("SOX6_AGTR1","SOX6_DDT", "SOX6_GFRA2","SOX6_PART1", "CALB1_CALCR","CALB1_CRYM_CCDC68",        
                              "CALB1_GEM","CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR")])
x_variable = "segment"
y_variable = "DaN"
x_lab = x
y_lab = paste0(y, " proportion")
colour_palette = segment_col
data_table <- dplot[dplot$Diagnosis == "CTR",c(x,y)]
# make violin plot
bxp <- ggviolin(
  data_table, x = x_variable, y = y_variable, 
  fill = x_variable, palette = colour_palette) + theme(legend.position = "none")
# perform t-test
data_table$x <- data_table[, x_variable]
data_table$y <- data_table[, y_variable]
stat.test <- data_table %>% t_test(y ~ x) %>% add_significance("p")
# add p-values to plot and save the result
v2 <- bxp + stat_pvalue_manual(stat.test,
                               y.position = max(data_table$y) * 1.4, step.increase = 0.1,
                               label = "p.signif") + ylab(y_lab) + xlab(x_lab)

arrange <- ggarrange(plotlist=list(v1,v2), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_DaN.confirm.png", arrange)


# 3. quantify neuronal proportions in TH + x region x Dx
# Violin plot function
violin_plot_function <- function(data_table, x_variable, y_variable, x_lab, y_lab, colour_palette) {
  # make violin plot
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette,yscale = "log10") + theme(legend.position = "none")
  
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  stat.test <- data_table %>% t_test(y ~ x) %>% add_significance("p.adj")
  print(stat.test)
  
  # add p-values to plot and save the result
  bxp <- bxp + stat_pvalue_manual(stat.test,
                                  y.position = max(data_table$y) * 1.4, step.increase = 0.1,
                                  label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  
  # return object
  return(bxp)
}

cell_types <- colnames(norm_weights)
contrast_matrix <- expand.grid(meta$segment,meta$ROI)
contrast_matrix <- contrast_matrix[!duplicated(contrast_matrix),]
dplot$Diagnosis <- as.factor(dplot$Diagnosis)
levels(dplot$Diagnosis) <- c("CTR","ILBD", "ePD", "lPD")

for (i in 1:nrow(contrast_matrix)){
  for (z in 1:length(cell_types)){
    print(paste0(contrast_matrix[i,2],sep=":",contrast_matrix[i,1]))
    print(z)
    x_variable = "Diagnosis"
    y_variable = cell_types[z]
    x_lab = paste0(contrast_matrix[i,2],sep=":",contrast_matrix[i,1])
    y_lab = paste0(y_variable, " proportion")
    colour_palette = segment_col
    data_table <- dplot[dplot$ROI == contrast_matrix[i,2] & dplot$segment == contrast_matrix[i,1],c(x_variable,y_variable)]
    
    colnames(data_table) <- gsub("-","_",colnames(data_table))
    y_variable <- gsub("-","_",y_variable)
    
    v1 <- violin_plot_function(data_table,x_variable ,y_variable,x_lab,y_lab, Diagnosis_col)
    arrange <- ggarrange(plotlist=list(v1), nrow=2, ncol=2, widths = c(2,2))
    ggsave(paste0("violin_cell.type_contrast",x_lab,z,".png"), arrange)
  }
}


### Heatmap plots
# aggregate data
x = norm_weights[meta$segment == "TH",]
y = meta[meta$segment == "TH",]
ct_agg_median_th <- as.data.frame(x) %>%
  group_by(y$ROI,y$Diagnosis) %>% 
  summarise_all("median")

x = norm_weights[meta$segment != "TH",]
y = meta[meta$segment != "TH",]
ct_agg_median_full <- as.data.frame(x) %>%
  group_by(y$ROI,y$Diagnosis) %>% 
  summarise_all("median")

# plot data frames
dplot <- ct_agg_median_th
dplot1 <- t(dplot[,3:ncol(dplot)])

# column annotations
anno_df = data.frame(
  Region = dplot$`y$ROI` ,
  Dx = dplot$`y$Diagnosis` 
)

ha = HeatmapAnnotation(df = anno_df,
                       col = list(Dx = c("CTR" = "grey","ILBD" = "#00AFBB", "ePD" = "#E7B800", "lPD"= "red" ),
                                  Region = c( "LC" = "#440154FF", "RN" = "#443A83FF",  "SND" = "#31688EFF", 
                                              "SNL" = "#21908CFF", "SNM" = "#35B779FF", "SNV" = "#8FD744FF", "VTA" = "#FDE725FF"))
)
# row colors
tmp <- row.names(dplot1)
tmp[!grepl("Ex|Inh|SOX6|CALB1",tmp)] <- "red3"  
tmp[grep("Ex|Inh|SOX6|CALB1",tmp)] <- "grey20"
row_annotation = c(col = tmp)

# draw heatmap
pdf("heatmap_cside_cellprop_TH.pdf")
hm <- Heatmap(dplot1, cluster_columns = FALSE,
              top_annotation = ha,
              column_split=sort(as.factor((dplot$`y$ROI`))),
              row_names_gp = gpar(col = row_annotation, fontsize = 7),
              heatmap_legend_param = list(title = "Cell Prop"))
draw(hm,
     column_title = "TH ROI's",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()



############################################################################################
###### Part 2b: Evaluate the cell-type proportion differences by Dx and Region
############################################################################################
# use cell-type proportions as expression data
exp_dat <- t(norm_weights)

# extract meta data
meta_dat <- as.data.frame(gxdat_s@meta.data)
table(row.names(meta_dat) == colnames(exp_dat))

# format numerical variables
meta_dat$DV200 <- as.numeric(meta_dat$DV200)
meta_dat$Age <- as.numeric(meta_dat$Age)
meta_dat$n_per.µm2.iNM <- as.numeric(meta_dat$n_per.µm2.iNM)
meta_dat$n_per.µm2.eNM <- as.numeric(meta_dat$n_per.µm2.eNM)
meta_dat$Diagnosis_stage <- as.character(meta_dat$Diagnosis)
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "CTR"] <- 0
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "ILBD"] <- 1
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "ePD"] <- 2
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "lPD"] <- 3
meta_dat$Diagnosis_stage <- as.numeric(meta_dat$Diagnosis_stage)

# format factors
factor_format <- c("Brainbank_ID","Sex","Diagnosis","Brainregion","ROI")
for (i in 1:length(factor_format)){
  meta_dat[,factor_format[i]] <- as.factor(meta_dat[,factor_format[i]])
}
#meta_dat$Brainregion <- factor(meta_dat$Brainregion, levels=c('A6','A9','RN','A10'))

# read in contrast matrix
cont_dat <- read_excel(contrast_path, sheet = "LIMMA_Voom")
cont_dat <- cont_dat[cont_dat$notes == "run",]

## Run Voom
res <- list() # list to collect results
for (val in 83:nrow(cont_dat)){
  # print model number
  print(cont_dat$model.number[val])
  
  # create targets df
  if (cont_dat$roi[val] != "NA"){
    print("segment & roi")
    targ <- meta_dat[meta_dat$segment %in% cont_dat$segment[val] & meta_dat$ROI %in% cont_dat$roi[val],]
  } else if (cont_dat$brainregion[val] != "NA"){
    print("segment & brainregion")
    targ <- meta_dat[meta_dat$segment %in% cont_dat$segment[val] & meta_dat$Brainregion %in% cont_dat$brainregion[val],]
  } else {
    print("segment")
    targ <- meta_dat[meta_dat$segment %in% cont_dat$segment[val],]
  }
  
  # create contrast arg list
  cont_list <- list()
  cont_vars <- factor_format
  for (z in 1:length(cont_vars)){
    contrasts(targ[,cont_vars[z]], contrasts = FALSE)
    cont_list[[z]] <- contrasts(targ[,cont_vars[z]], contrasts = FALSE)
  }
  names(cont_list) <- factor_format
  
  # create design matrix
  options(na.action='na.omit')
  design <- model.matrix(reformulate(cont_dat$model.matrix[val]),
                         data=targ, 
                         drop = FALSE,
                         contrasts.arg = cont_list)
  
  # make names
  colnames(design) <- make.names(colnames(design))
  design <- design[,!colnames(design) %in% c("DiagnosisILBD.n_per.µm2.iNM","DiagnosisILBD.n.iNM")]
  
  # create expression df & targ with design row.names
  targ <- targ[row.names(design),]
  y <- exp_dat[,row.names(targ)]
  
  # fit design
  v <- voom(y,design)
  vfit <- lmFit(v)
  
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=cont_dat$contrast[val],levels=design)
  #cont.matrix <- makeContrasts(A="DiagnosisCTR - DiagnosisILBD",levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant DEGs and assign to list
  tmp <- topTable(vfit2,number = Inf,sort.by="P")
  res[[val]] <- tmp[tmp$P.Value < 0.05,]
  
}

# ## SAVE results list
# cside_cellprop_res <- res
# save(cside_cellprop_res,cont_dat,file="cside_cellprop.rds")

# ## save results (p< 0.05) to excel spreadsheets
# # i. name each item by description
# names_1 <- paste0(cont_dat$description, "(",cont_dat$segment,")")
# names(voom_res) <- names_1

# plots to check results by group
# colour palettes
Diagnosis_col = c("grey","#00AFBB", "#E7B800","red")
Brainregion_col = viridis(3)
ROI_col = viridis(6)

# Violin plot function
violin_plot_function <- function(data_table, x_variable, y_variable, x_lab, y_lab, colour_palette) {
  # make violin plot
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none")
  
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  stat.test <- data_table %>% t_test(y ~ x) %>% add_significance("p.adj")
  print(stat.test)
  
  # add p-values to plot and save the result
  bxp <- bxp + stat_pvalue_manual(stat.test,
                                  y.position = max(data_table$y) * 1.4, step.increase = 0.1,
                                  label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  
  #return object
  return(bxp)
}

# example plots
dat <- t(norm_weights)

cell_name <- "Ex_LAMP5_BAIAP3"
cell_type <- dat[cell_name,]
dplot <- cbind(meta_dat,cell_type)
dplot <- dplot[dplot$ROI == "SNV" & dplot$segment == "TH",c("Diagnosis","cell_type")]
#dplot$Brainregion <- factor(dplot$Brainregion)

v1 <- violin_plot_function(dplot,"Diagnosis","cell_type","Diagnosis",paste0(cell_name, " proportion"), Diagnosis_col)




model <- lm(cell_type ~ Diagnosis, data = dplot)




############################################################################################
###### Part 3: Run C-SIDE 2 regions
############################################################################################
### data
norm_weights <- normalize_weights(myRCTD@results$weights)
myRCTD@config$doublet_mode <- "full"
meta <- as.data.frame(gxdat_s@meta.data)
table(row.names(meta) == row.names(norm_weights))
dplot <- cbind(meta,norm_weights)

# define thresholds
cell_type_threshold = 1
myRCTD@config$max_cores <- 8
weight_threshold = 0.1

# read in contrast matrix
cont_dat <- read_excel(contrast_path)
contrast_matrix <- cont_dat[cont_dat$model == "C-SIDE 2 region",]

# ## Running CSIDE - contrast per segment and roi
cside_res <- list()
for (i in 4:nrow(contrast_matrix)){
  # define subsets
  segment_run <- as.character(contrast_matrix[i,"segment"])
  roi_run <- as.character(contrast_matrix[i,"roi"])
  print(paste0(segment_run,".", roi_run,".",contrast_matrix$contrast[i]))
  
  # get meta data and barcodes
  meta_run <- meta[meta$segment == segment_run & meta$ROI == roi_run,]
  barcodes_run <- row.names(meta_run)
  
  # get list of cells that pass thresholds
  agg_cell_types <- aggregate_cell_types(myRCTD, barcodes_run,doublet_mode = F)
  cell_cell_type_threshold <- names(agg_cell_types)[agg_cell_types > cell_type_threshold]
  cells_weight_threshold <- names(colSums(norm_weights))[colSums(norm_weights) > weight_threshold]
  cell_run <- intersect(cells_weight_threshold,cell_cell_type_threshold)
  
  # define explanatory variables
  explanatory.variable <- as.numeric(as.factor(meta_run$Diagnosis == contrast_matrix$contrast[i])) - 1
  names(explanatory.variable) <-  barcodes_run 
  
  myRCTD_test <- run.CSIDE.single(myRCTD,
                                  explanatory.variable,
                                  cell_types = cell_run,
                                  cell_type_threshold = cell_type_threshold,
                                  weight_threshold = weight_threshold,
                                  fdr = 0.25, 
                                  doublet_mode = FALSE) 
  
  ## CSIDE results
  all_gene_list <- myRCTD_test@de_results$all_gene_list
  sig_gene_list <- lapply(all_gene_list, function(x) x[order(x$p_val),])
  cside_res[[i]] <- sig_gene_list
  names(cside_res[[i]]) <- paste0(segment_run,",", roi_run,",",contrast_matrix$contrast[i],",",names(cside_res[[i]]))
}

# ssave results list and contrasts
# save(cside_res, file = "cside_res_snvsnd.ctrpd.rds")


############################################################################################
###### Part 4: Run C-SIDE multiple regions
############################################################################################
load(rdata)
norm_weights <- normalize_weights(myRCTD@results$weights)
myRCTD@config$doublet_mode <- "full"
meta <- as.data.frame(gxdat_s@meta.data)
table(row.names(meta) == row.names(norm_weights))

# read in contrast matrix
cont_dat <- read_excel(contrast_path, sheet = "C-SIDE")
cont_dat <- cont_dat[cont_dat$notes == "run",]

# params
cell_type_threshold = 1
myRCTD@config$max_cores <- 12
weight_threshold = 0.1

# ## Running CSIDE - regions/ multiple contrasts
#cside_res <- list()
for (val in 12:nrow(cont_dat)){
  print(cont_dat$model.number[val])
  # define contrast items 
  contrast_items <- unlist(strsplit(cont_dat$contrast[val],","))
  
  # create targets df
  if (cont_dat$roi[val] != "NA"){
    print("segment & roi")
    targ <- meta_dat[meta_dat$segment %in% cont_dat$segment[val] & meta_dat$ROI %in% cont_dat$roi[val],]
    targ <- targ[targ[,colnames(targ) == cont_dat$contrast_column[val]] %in% contrast_items,]
    
  } else if (cont_dat$brainregion[val] != "NA"){
    print("segment & brainregion")
    targ <- meta_dat[meta_dat$segment %in% cont_dat$segment[val] & meta_dat$Brainregion %in% cont_dat$brainregion[val],]
    targ <- targ[targ[,colnames(targ) == cont_dat$contrast_column[val]] %in% contrast_items,]
    
  } else {
    print("segment")
    targ <- meta_dat[meta_dat$segment %in% cont_dat$segment[val],]
    targ <- targ[targ[,colnames(targ) == cont_dat$contrast_column[val]] %in% contrast_items,]
    
  }
  # get barcodes
  barcodes_run <- row.names(targ)
  
  # build item list and assignment matrix
  #item_list <- split(barcodes_run, f = droplevels(targ[,colnames(targ) == cont_dat$contrast_column[val]]))
  item_list <- split(barcodes_run, f = targ[,colnames(targ) == cont_dat$contrast_column[val]])
  X <- build.designmatrix.regions(myRCTD,item_list)
  barcodes2 <- rownames(X)
  
  # get list of cells that pass thresholds
  agg_cell_types <- aggregate_cell_types(myRCTD, barcodes_run,doublet_mode = F)
  cell_cell_type_threshold <- names(agg_cell_types)[agg_cell_types > cell_type_threshold]
  cells_weight_threshold <- names(colSums(norm_weights))[colSums(norm_weights) > weight_threshold]
  cell_run <- intersect(cells_weight_threshold,cell_cell_type_threshold)
  
  # run
  myRCTD_test <- run.CSIDE(myRCTD, X, barcodes_run, cell_run , 
                           test_mode = 'categorical', gene_threshold = 5e-5, 
                           doublet_mode = F, cell_type_threshold = cell_type_threshold, fdr = 0.25, weight_threshold = 0, 
                           log_fc_thresh = 0)
  ## CSIDE results
  all_gene_list <- myRCTD_test@de_results$all_gene_list
  #sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val_best < 0.05,])
  sig_gene_list <- lapply(all_gene_list, function(x) x[order(x$p_val_best),])
  cside_res[[val]] <- sig_gene_list 
}

# save results list
#save(cside_res,cont_dat,file="cside_deg.rds")


############################################################################################
###### Part 5: Plots of C-SIDE results
############################################################################################
load("cside_deg.rds")
# check gene in C-SIDE results

for (z in 1:15){
  tmp <- cside_res[[z]]
  print(z)
  for (i in 1:length(tmp)){
    print(names(tmp)[i])
    if (!is.na(tmp[[i]]["RGN","p_val_best"])){
      if (tmp[[i]]["RGN","p_val_best"] < 0.05){
        print(tmp[[i]]["PGK1",])
      }}
  }}


## Plots
# barplot for n sig genes x cell-type
x = unlist(lapply(sig_gene_list, nrow))

dplot <- as.data.frame(cbind(cell_type=names(x),sig_genes=x))
dplot$sig_genes <- as.numeric(dplot$sig_genes) 

g2 <- ggplot(dplot, aes(x=reorder(cell_type,-x),y=x)) +
  geom_bar(stat="identity", width=0.7)+  
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cell Type") + ylab("Sig CtDEG's (n)") 

# save plots
ggsave("ctDEG_results_volcano_full.ctr.a10a9rn.png", 
       ggarrange(g2,nrow=2,ncol=1)
       , device = "png")


# volcano plots
dat_tmp <- as.data.frame(all_gene_list$SOX6_DDT)

EnhancedVolcano(toptable = dat_tmp , x = "log_fc_best",y = "p_val_best",
                lab = rownames(dat_tmp),pCutoff = 0.05, FCcutoff = 0.5,
                title = "volcano_ctDEG_geomx_CTR",subtitle = "")

# define genesets
all_gene_sets = msigdbr(species = "human")
# using "H: hallmark gene sets"
# using "C2: curated gene sets"
# using "C5: ontology gene sets"
# using "C8: cell type signature gene sets"
msigdbr_df <- all_gene_sets %>%
  dplyr::filter(gs_cat == "C5")

# create list for fgseas
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# select res from list, perform FGSEA and plot results
for (i in 1:length(all_gene_list)){
  i = which(names(sig_gene_list) == "SOX6_DDT" )
  print(i)
  res <- as.data.frame(sig_gene_list[[i]])
  name1 <- names(sig_gene_list)[i]
  gene_stats <- res$log_fc
  names(gene_stats) <- row.names(res)
  
  # FGSEA
  fgseaRes <- fgsea(pathways = msigdbr_list,
                    stats    = gene_stats,
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes[order(pval), ]
  head(fgseaRes)
  
  # plot
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  pdf(paste0("FGSEA_C5celltype_cside_fullroi.ctr.a10a9rn",name1,".pdf"))
  print(plotGseaTable(msigdbr_list[topPathways], gene_stats, fgseaRes,
                      gseaParam=0.5, pathwayLabelStyle=list(size=6, color="black")))
  dev.off()
}

# quantitative expression plot
barcodes <- myRCTD_ctr.fulla9a10rn@internal_vars_de$barcodes
gene <- 'ATP6V0B'
Y <- myRCTD@spatialRNA@counts[gene, barcodes]
Yn <- Y / myRCTD_ctr.fulla9a10rn@spatialRNA@nUMI[barcodes]
dc_prop <- myRCTD_ctr.fulla9a10rn@internal_vars_de$my_beta[,'SOX6_DDT']
thresh <- median(dc_prop)
reg <- myRCTD_ctr.fulla9a10rn@internal_vars_de$X2[barcodes,2]
pred <- predict_CSIDE_all(myRCTD_ctr.fulla9a10rn, gene)


gene_fits <- myRCTD_ctr.fulla9a10rn@de_results$gene_fits
cell_types <- "SOX6_DDT"
my_beta <- myRCTD_ctr.fulla9a10rn@internal_vars_de$my_beta
cell_type_ind <- which(myRCTD_ctr.fulla9a10rn@internal_vars_de$cell_types == "SOX6_DDT")
X2_mat <- myRCTD_ctr.fulla9a10rn@internal_vars_de$X2[rownames(my_beta),]

pred_ct <- predict_CSIDE(cell_type_ind, gene_fits, gene, X2_mat)



ct_thresh <- 0
myRCTDde <- myRCTD_ctr.fulla9a10rn
cur_cell_types <- "SOX6_DDT"
p_df <- get_quant_df(myRCTDde, myRCTDde@de_results$gene_fits, myRCTDde@internal_vars_de$cell_types, 
                     cur_cell_types, gene, multi_region = T, prop_thresh = ct_thresh)

table(row.names(p_df) == row.names(meta[barcodes,]))

d1 <- cbind(meta[barcodes,],gene_exp = p_df$pred)
d1$Brainregion <- factor(d1$Brainregion, levels = c("A9","A10","RN"))
p1a <- ggplot(d1, aes(x=Brainregion, y=gene_exp, fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1),facet.by = "Brainregion", short.panel.labs = FALSE) + 
  scale_fill_brewer(palette='viridis') + theme_classic() + ylab(paste0(gene," expression")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) +
  stat_compare_means( label = "p.signif")


d2 <- cbind(meta[barcodes,],gene_exp = Y )
d2$Brainregion <- factor(d2$Brainregion, levels = c("A9","A10","RN"))
p1b <- ggplot(d2, aes(x=Brainregion, y=gene_exp, fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1),facet.by = "Brainregion", short.panel.labs = FALSE) + 
  scale_fill_brewer(palette='viridis') + theme_classic() + ylab(paste0(gene," expression")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) +
  stat_compare_means( label = "p.signif")

# save plots
ggsave("boxplot_ATP6V0B_full.ctr.a10a9rn.png", 
       ggarrange(p1b ,nrow=2,ncol=1), device = "png")


