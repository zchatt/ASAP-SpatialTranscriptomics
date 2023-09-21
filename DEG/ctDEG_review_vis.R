library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)


# 1. test RCTD and CSIDE method using snRNAseq data from Kamath et al -> csDEG
# 2. test CSIDE method using PAGE cell-enrichments from data from Kamath et al -> csDEG
# 3. test TOAST method using PAGE cell-enrichments from data from Kamath et al -> csDEG
## It would be ideal to use PAGE cell-enrichments for ctDEG detection due to the limited scRNAseq datasets

############################################################################################
#### Inputs
############################################################################################

analysis_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/SN_210823"

## read in Kamath et al cell ranger output
kamath_obj <- Read10X("/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger")

## read in metadata / umap of cell IDs
#kamath_meta <- read.table(file = '/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger/METADATA_PD.tsv', sep = '\t', header = TRUE)
umap_de <- read.table(file = '/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger/da_UMAP.tsv', sep = '\t', header = TRUE)
umap_de <- umap_de[-1,]

# load seurat object eg. made with "batch.cluster_vis.R"
load(paste0(basename(analysis_dir),"_merged.seurat.Rdata"))

############################################################################################
###### Part 1: RCTD & CSIDE analysis
############################################################################################

# setwd
setwd(analysis_dir)

# subset for DA neurons
kamath_obj_sub <- kamath_obj[,umap_de$NAME]

### Data Preprocessing and running RCTD
if(!file.exists(file.path(savedir,'myRCTD_visium_SN.rds'))) {
  # spatial data
  # select just control SN sample
  run_index =  brain.merge@meta.data$run_name == "V52Y16-079-A1" & brain.merge@meta.data$tissue == "SN" 
  # load in counts matrix
  counts <- brain.merge@assays$SCT@counts[,run_index]  # load in counts matrix
  # load in coordinates
  coords <- brain.merge@meta.data[run_index,]
  coords <- coords[,c("array_col","array_row")]
  nUMI <- colSums(counts) # In this case, total counts per spot
  puck <- SpatialRNA(coords, counts, nUMI)
  barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
  plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                       title ='plot of nUMI') 
  # single-cell data
  cell_types <- setNames(umap_de$Cell_Type, umap_de$NAME)
  cell_types <- as.factor(umap_de$Cell_Type) # convert to factor data type
  cell_types <- as.factor(cell_types)
  nUMI <- colSums(kamath_obj_sub)
 
  # create reference and run RCTD
  reference <- Reference(kamath_obj_sub, cell_types, nUMI)
  myRCTD <- create.RCTD(puck, reference, max_cores = 2)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD,file.path(savedir,'myRCTD_visium_SN.rds'))
}

### Exploring the full mode results
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
cell_types <- c('SOX6_AGTR1','CALB1_GEM')
print(head(norm_weights[,cell_types])) 

ggsave("qc_4_scatter.normeffect.png",
       ggarrange(p1,p2, ncol=2,nrow=2) + bgcolor("white"),
       device = "png")

p1 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'SOX6_AGTR1'], ylimit = c(0,0.5), 
                     title ='plot of SOX6_AGTR1 weights', size = 1) 

ggsave("cside_1_sox6agtr1.png", p1, device = "png")

p2 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'CALB1_GEM'], ylimit = c(0,0.5), 
                     title ='plot of CALB1_GEM weights', size = 1) 

ggsave("cside_2_calb1gem.png", p2, device = "png")

## Choose regions for ctDEG analysis
### Create SpatialRNA object
region_left <- barcodes[which(myRCTD@spatialRNA@coords$x < quantile(myRCTD@spatialRNA@coords$x, 1/3))]
region_right <- barcodes[which(myRCTD@spatialRNA@coords$x > quantile(myRCTD@spatialRNA@coords$x, 2/3))]
region_middle <- setdiff(barcodes, union(region_left, region_right))
region_list <- list(region_left, region_right, region_middle)

class_num <- rep(0, length(barcodes)); names(class_num) <- barcodes
class_num[region_middle] <- 1; class_num[region_right] <- 2 

p3  <- plot_class(myRCTD@spatialRNA, barcodes, factor(class_num),  title ='plot of regions')
ggsave("cside_3_regionsplit.png", p3, device = "png")

## Running CSIDE
# ctDEG
myRCTD@config$max_cores <- 2
myRCTD@config$doublet_mode = 'full'

# build design matrix
X <- build.designmatrix.regions(myRCTD,region_list)
barcodes <- rownames(X)
myRCTD <- run.CSIDE(myRCTD, X, barcodes,cell_types = cell_types,
                    cell_type_threshold = 10,
                    fdr = 0.01, 
                    doublet_mode = FALSE, 
                    weight_threshold = 0.1,
                    test_mode = 'categorical') 

saveRDS(myRCTD,file.path(savedir,'myRCTDde_visium_SN_regions.rds'))

## CSIDE results
all_gene_list <- myRCTD@de_results$all_gene_list
sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val_best < 0.05,])
sig_gene_list <- lapply(all_gene_list, function(x) x[order(x$p_val_best),])

############################################################################################
###### Part 2: CSIDE analysis using PAGE reference
############################################################################################
# replace RCTD cell-type weights with PAGE results
myRCTD@results$weights <- brain.merge@meta.data$

### Exploring the full mode results
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
cell_types <- c('SOX6_AGTR1','CALB1_GEM')
print(head(norm_weights[,cell_types])) 

ggsave("qc_4_scatter.normeffect.png",
       ggarrange(p1,p2, ncol=2,nrow=2) + bgcolor("white"),
       device = "png")

p1 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'SOX6_AGTR1'], ylimit = c(0,0.5), 
                           title ='plot of SOX6_AGTR1 weights', size = 1) 

ggsave("cside_1_sox6agtr1.png", p1, device = "png")

p2 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'CALB1_GEM'], ylimit = c(0,0.5), 
                           title ='plot of CALB1_GEM weights', size = 1) 

ggsave("cside_2_calb1gem.png", p2, device = "png")

## Choose regions for ctDEG analysis
### Create SpatialRNA object
region_left <- barcodes[which(myRCTD@spatialRNA@coords$x < quantile(myRCTD@spatialRNA@coords$x, 1/3))]
region_right <- barcodes[which(myRCTD@spatialRNA@coords$x > quantile(myRCTD@spatialRNA@coords$x, 2/3))]
region_middle <- setdiff(barcodes, union(region_left, region_right))
region_list <- list(region_left, region_right, region_middle)

class_num <- rep(0, length(barcodes)); names(class_num) <- barcodes
class_num[region_middle] <- 1; class_num[region_right] <- 2 

p3  <- plot_class(myRCTD@spatialRNA, barcodes, factor(class_num),  title ='plot of regions')
ggsave("cside_3_regionsplit.png", p3, device = "png")

## Running CSIDE
# ctDEG
myRCTD@config$max_cores <- 2
myRCTD@config$doublet_mode = 'full'

# build design matrix
X <- build.designmatrix.regions(myRCTD,region_list)
barcodes <- rownames(X)
myRCTD <- run.CSIDE.regions(myRCTD, X, barcodes,
                            cell_type_threshold = 10,
                            fdr = 0.01, 
                            doublet_mode = FALSE, 
                            weight_threshold = 0.8) 

saveRDS(myRCTD,file.path(savedir,'myRCTDde_visium_SN_regions.rds'))

############################################################################################
###### Part 3: TOAST using PAGE reference
############################################################################################


# run TOAST
library(TOAST)


Design_out <- makeDesign(design, Prop)
fitted_model <- fitModel(Design_out, Y_raw)
fitted_model$all_coefs # list all phenotype names
fitted_model$all_cell_types # list all cell type names
# coef should be one of above listed phenotypes
# cell_type should be one of above listed cell types
res_table <- csTest(fitted_model, coef = "age", 
                    cell_type = "Neuron", contrast_matrix = NULL)
head(res_table)
