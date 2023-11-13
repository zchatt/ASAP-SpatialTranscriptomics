# compare xenograft human cells with human and iPSCs


# 

# comparisons
# 1) hESC-derived dopaminergic transplants from Tiklová et al, 2020
  #- transplanted ventral midbrain (VM)-patterned human embryonic stem cells (hESCs)
  #and VM fetal tissue into adult rat brain. Performed snRNAseq on in vivo and in vitro. 
  # 660 analyzed cells before grafting (404 cells of hESC origin, 256 cells of fetal origin)
  # 746 cells grafted to striatum (683 cells of hESC origin, grafted rats n=2; 63 cells of fetal origin, grafted rats n=2)
  # described in GSE118412_series_matrix.txt.gz
# 2a) human SNpc DaNs from Kamath et al. 2022

# 2b) human SNpc DaNs from Siletti et al. 2023

# read in supercluster splatter seurat data
clust_splat <- readRDS("/Users/zacc/USyd/spatial_transcriptomics/analysis/asap_kirik/siletti/Supercluster_Splatter.rds")
clust_splat[["CellID"]] <- row.names(clust_splat@meta.data) 
# read in DA cell annotations
da_annots <- read.delim(file="/Users/zacc/USyd/spatial_transcriptomics/analysis/asap_kirik/siletti/DA_cellmeta.csv", 
                        sep=",", header=T, row.names = 1)

# read in all cell annotations
clust_splat_dan <- subset(clust_splat, subset = CellID %in% da_annots$CellID)


# 3) PPMI - 
# FOUNDIN-PD - The Foundational Data Initiative has generated single-cell RNAseq data on day-65 iPSC-derived dopaminergic neurons from the PPMI collection
# 96 cell lines

# need to subset on gbhpc
# rsync -au /Users/zacc/USyd/ASAP/snrna_seq/FOUNDIN-PD_150.8_SCRN_Seurat.tar.gz zac@10.65.82.197:/home/zac/FOUNDIN

readRDS("/Users/zacc/USyd/ASAP/snrna_seq/FOUNDIN-PD_150.8_SCRN_Seurat/iPSCsDopaALL_integratedAfterBroadCellType.RDS")





#/Users/zacc/USyd/spatial_transcriptomics/analysis/asap_kirik/PPMI_156_MJFF-1651.scRNAseq.data.h5ad # unsure what this is?
# /Users/zacc/PPMI_RNAseq_IR3_Analysis.tar.gz # whole blood RNAseq



# Azimuth


