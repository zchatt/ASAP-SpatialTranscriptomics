# Cell-type enrichment

# libraries
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/giotto_env.R")
library(preprocessCore)
library(Giotto)
library(Seurat)
library(ggpubr)
library(harmony)

############################################################################################
#### Inputs
############################################################################################

analysis_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/SN_210823"

# load cell-enrichment markers 
load("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/celltype_markers/tables/gene_marker_list.cellenrich.R")

############################################################################################
#### Part 1 : Convert Seurat to Giotto and calculate cell-type enrichment
############################################################################################
setwd(analysis_dir)

# load seurat object eg. made with "batch.cluster_vis.R"
load(paste0(basename(analysis_dir),"_merged.seurat.Rdata"))

## create Giotto object. Note; seuratToGiotto does not work with multiple slices
#  normalized scaled data
cd <- brain.merge@assays$SCT@scale.data
# Create a Giotto object
giotto_obj = createGiottoObject(expression = cd)
giotto_obj = setExpression(giotto_obj,cd,name = "custom")

# create PAGE matrix
PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = names(gene_marker_list[1:25]),sign_list = gene_marker_list[1:25])

# perform PAGE
giotto_obj = runPAGEEnrich(gobject = giotto_obj, 
                           sign_matrix = PAGE_matrix_1, 
                           expression_values = "custom",
                           min_overlap_genes = 1)
# extract PAGE results 
page_dat <- as.data.frame(get_spatial_enrichment(gobject = giotto_obj,
                                                 enrichm_name = "PAGE",
                                                 output = "data.table"))

row.names(page_dat) <- page_dat$cell_ID
page_dat <- page_dat[,!colnames(page_dat) %in% c("cell_ID")]

# add to seurat object
brain.merge <- AddMetaData(brain.merge,
                       metadata = page_dat)

# save
save(paste0(basename(analysis_dir),"_merged.seurat.Rdata"))
