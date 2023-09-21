# Cell-type enrichment - geomx

# libraries
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/giotto_env.R")
library(preprocessCore)
library(Giotto)
library(ggpubr)
library(harmony)
library(Seurat)

############################################################################################
#### Inputs
############################################################################################

analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis"
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/geomx_sep2023_seurat.Rdata"

# load cell-enrichment markers 
load("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/celltype_markers/tables/gene_marker_list.cellenrich.R")

############################################################################################
#### Part 1 : Convert Seurat to Giotto and calculate cell-type enrichment
############################################################################################
setwd(analysis_dir)

# load normalised geomx data
load(rdata)

## create Giotto object. Note; seuratToGiotto does not work with multiple slices
# Create a Giotto object
giotto_obj = createGiottoObject(raw_exprs = gxdat_s@assays$GeoMx@data,
                                custom_expr = scale(gxdat_s@assays$GeoMx@data),
                                instructions = my_instructions)

# create PAGE matrix
tmp <- lapply(gene_marker_list,function(x) x[1:100])
PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = names(tmp),sign_list = tmp)
PAGE_matrix_1 = PAGE_matrix_1[row.names(PAGE_matrix_1) %in% row.names(giotto_obj@custom_expr),]


# perform PAGE
giotto_obj = runPAGEEnrich(gobject = giotto_obj, 
                           sign_matrix = PAGE_matrix_1, 
                           min_overlap_genes = 1)

# extract PAGE results 
page_dat <- as.data.frame(giotto_obj@spatial_enrichment$PAGE)
row.names(page_dat) <- page_dat$cell_ID
page_dat <- page_dat[,!colnames(page_dat) %in% c("cell_ID")]

# add to seurat object
gxdat_s <- AddMetaData(gxdat_s,metadata = page_dat)

# save
save(gxdat_s,file=rdata)

############################################################################################
#### Part 2 : Plotting
############################################################################################

# extract PAGE results
page_dat
