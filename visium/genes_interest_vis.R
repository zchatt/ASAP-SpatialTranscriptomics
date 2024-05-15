### Analysis of genes or gene-sets of interest from Visium 10X datasets

library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(plyr)
library(data.table)
library(ggpubr)

############################################################################################
#### Inputs
############################################################################################

load("/Volumes/PRJ-2022ASAPs/Bioinformics/datasets/seurat_objects/SN_210823_merged.seurat.Rdata") # visium seurat object
working_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/SN_210823"
setwd(working_dir)

############################################################################################
#### Plot Gene of Interest
############################################################################################
# subset SN ventral tier
brain.merge_sn <- subset(brain.merge, subset = tissue == "SN" )

# genes of interest
gene_interest <- c("CACNA1C","CACNA1H","KCND3","CACNA2D1","KCNIP1","KCNA10","KCNC1","HCN1")
gene_interest <- c("TYR","MITF")

# plot each slice and combine
for (i in gene_interest){
  p1 <- SpatialFeaturePlot(brain.merge_sn, 
                           images = "slice1.1", features = i, pt.size.factor = 2,
                           alpha = c(1, 1), ncol = 1, image.alpha = 0, stroke = 0) + theme(aspect.ratio = 2) + ggtitle("PD")
  
  p2 <- SpatialFeaturePlot(brain.merge_sn,
                           images = "slice1", features = i, pt.size.factor = 2,
                           alpha = c(1, 1), ncol = 1, image.alpha = 0,  stroke = 0) + theme(aspect.ratio = 2)  + ggtitle("CTR")
  
  arrange <- ggarrange(plotlist=list(p2,p1), nrow=1, ncol=2)
  ggsave(paste0(i,"_seurat_spat.pdf"),arrange)
}
