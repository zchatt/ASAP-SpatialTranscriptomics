# Joining multiple datasets, batch correction and clustering

# libraries
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/giotto_env.R")
library(preprocessCore)
1library(Giotto)
library(Seurat)
library(ggpubr)
library(harmony)

############################################################################################
#### Inputs
############################################################################################
run_names = c("V52Y16-079-A1","V52Y16-079-B1")
results_folders = c("/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/V52Y16-079-A1",
                    "/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/V52Y16-079-B1")
analysis_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/SN_210823"

############################################################################################
#### Part 1 : Join and Batch correct
############################################################################################
setwd(analysis_dir)

#make syntactically valid names
synt_name = make.names(run_names)

# load seurat objects and assign to synt_name(s)
for (k in 1:length(run_names)){
  print(k)
  load(paste0(results_folders[k],"/",run_names[k],"_qfn.seurat.Rdata"))
  assign(synt_name[k],brain)
}

# merge objects
brain.merge <- merge(get(unlist(synt_name[1])),get(unlist(synt_name[2])))

if (length(synt_name) > 2) {
  for (k in 3:length(synt_name)){
    print(k)
    brain.merge <- merge(brain.merge,get(unlist(synt_name[k])))
  }
}

# set default assay
DefaultAssay(brain.merge) <- "SCT"

# clustering
var_features <- list()
for (k in 1:length(synt_name)){
  var_features[[k]] <- VariableFeatures(get(unlist(synt_name[k])))
}

VariableFeatures(brain.merge) <- c(unlist(var_features))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)

# add names to metadata
brain.merge@meta.data$run_name <- rep(run_names,unlist(lapply(synt_name,function(x) ncol(get(x)))))

## plot PC1/2
p1 <- DimPlot(object = brain.merge, reduction = "pca", pt.size = .1, group.by = "run_name")
p2 <- VlnPlot(object = brain.merge, features = "PC_1", pt.size = .1, group.by = "run_name")

ggsave("bclus_1_pca12.png",
       ggarrange(p1, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

ggsave("bclus_2_pca.exp.png",
       ggarrange(p2, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

## Batch correct using Harmony
brain.merge <- brain.merge %>% 
  RunHarmony("run_name", plot_convergence = TRUE)

## plot PC1/2
p1 <- DimPlot(object = brain.merge, reduction = "harmony", pt.size = .1, group.by = "run_name")
p2 <- VlnPlot(object = brain.merge, features = "harmony_1", pt.size = .1, group.by = "run_name")

ggsave("bclus_3_pca12.harmony.png",
       ggarrange(p1, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

ggsave("bclus_4_pca.exp.harmony.png",
       ggarrange(p2, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

############################################################################################
#### Part 2 : Dimension reduction
############################################################################################

brain.merge <- brain.merge %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

## plot
p1 <- DimPlot(brain.merge, reduction = "umap", group.by = "run_name", pt.size = .1, split.by = 'run_name')
p2 <- DimPlot(brain.merge, reduction = "umap", label = TRUE)
p3 <-SpatialDimPlot(brain.merge)

ggsave("bclus_5_runname.umap.png",
       ggarrange(p1, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

ggsave("bclus_6_clust.umap.png",
       ggarrange(p2, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

ggsave("bclus_7_clust.spatdim.png",
       ggarrange(p3, ncol=1,nrow=1) + bgcolor("white"),
       device = "png")

# save data
save(brain.merge, file = paste0(basename(analysis_dir),"_merged.seurat.Rdata"))

