# splitting R data from into sub-projects and removing non-shareable information 


#----------- Inputs -----------
rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/"


#----------- Split to project ----------- 
# setwd
setwd(analysis_dir)
# Load normalized and batch corrected rdata
load(rdata)
# extract expression values
exp_dat <- as.matrix(gxdat_s@assays$RNA$counts)


# Dataset 1: GeoMx human TH+ masked SNC, VTA
gxdat_s <- gxdat_s


# Dataset 2: GeoMx human unmasked ROIs




# Dataset 3: GeoMx human TH masked LC


