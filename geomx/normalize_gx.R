# Normalization of GeoMx using quantile normalization
# Please refer to "ASAP-SpatialTranscriptomics/geomx/normalization_review_gx.R" for methods in selecting quantile normalization

# Libraries
library(preprocessCore)

#########
## INPUT ##
#########

run_name <- "hu_brain_001"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/analysis"

#########
## RUN ##
#########
setwd(analysis_dir)
# Load qc data generated from geomx_lowlevel.R, should be in "/analysis" of the name "run_name_qc.gx.Rdata"
load(paste0(run_name, "_qc.gx.Rdata"))
# Quantile Normalization
norm.quantile <- normalize.quantiles(raw)
dimnames(norm.quantile) <- dimnames(exprs(target_demoData))
# Write to file
save(norm.quantile,file=paste0(run_name,"_norm.gx.Rdata"))
