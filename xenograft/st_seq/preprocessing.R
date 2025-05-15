# scripts for preprocessing Xenograft data generated from Visium

# libs
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(readxl)
library(ggpubr)
library(harmony)


#----------- Input -----------

# working dir
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft")

# read in meetadata - contains locations to  spaceranger output
meta_dat <- read_excel("/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft/xenograft_metadata.xlsx")
colnames(meta_dat) <- make.names(colnames(meta_dat))
meta_dat <- as.data.frame(meta_dat)
row.names(meta_dat) <- meta_dat$Sample_Name

#----------- create Seurat objects and merge -----------

# # directories of spaceranger "outs"
# S11W1_LB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_LB/spaceranger_count_cytassist/S11W1_LB/outs"
# S11W1_LT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_LT/spaceranger_count_cytassist/S11W1_LT/outs"
# S11W1_RB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_RB/outs"
# S11W1_RT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_RT/outs"
# 
# S12W2_LB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/SW12W2_LB/outs"
# S12W2_LT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S12W2_LT/outs"
# S12W2_RB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S12W2_RB/outs"
# S12W2_RT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S12W2_RT/outs"
# 
# # create lists
# list_sr <- c(S11W1_LB,S11W1_LT,S11W1_RB,S11W1_RT,
#              S12W2_LB,S12W2_LT,S12W2_RB,S12W2_RT)
# names(list_sr) <- c("S11W1_LB","S11W1_LT","S11W1_RB","S11W1_RT",
#                     "S12W2_LB","S12W2_LT","S12W2_RB","S12W2_RT")


# create lists of spaceranger out directories to parse
list_sr <- meta_dat$spaceranger.out
names(list_sr) <- meta_dat$Sample_Name

# load 10X from spaceranger and SCTransform
seurat_object_list <- list()
for (i in 1:length(list_sr)) {
  seurat_object_list[[i]] <- Load10X_Spatial(list_sr[i],
                                             slice = names(list_sr)[i])
  seurat_object_list[[i]]@meta.data$Sample_Name <- names(list_sr)[i]
}

# merge 
seurat_obj_merge <- merge(seurat_object_list[[1]],
                  y = seurat_object_list[2:length(seurat_object_list)], 
                  add.cell.ids = names(seurat_object_list), project = "xenograft")

# add metadata
md <- as.data.frame(seurat_obj_merge@meta.data)
md$row_id <- row.names(md)
tmp <- merge(md, meta_dat, by = "Sample_Name", all.x = TRUE, suffixes = c("", ".drop"))
tmp <- tmp[, !grepl("\\.drop$", colnames(tmp))]
row.names(tmp) <- tmp$row_id
tmp <- tmp[row.names(seurat_obj_merge@meta.data), ]

seurat_obj_merge <- AddMetaData(seurat_obj_merge,tmp)

# k-means clustering k=2 on each Sample_Name independently
set.seed(123)
md <- as.data.frame(seurat_obj_merge@meta.data)
res <- list()

for (i in 1:length(unique(md$Sample_Name))) {
  tmp <- md[md$Sample_Name == unique(md$Sample_Name)[i],]
  kmeans_result <- kmeans(tmp[,"nCount_Spatial"], centers=2)
  means <- tapply(tmp[,"nCount_Spatial"], kmeans_result$cluster, mean)
  
  # check if we need to swap cluster labels
  if (means[1] < means[2]) {
    # swap cluster assignments
    kmeans_result$cluster <- ifelse(kmeans_result$cluster == 1, 2, 1)
  }
  res[[i]] <- kmeans_result$cluster
}

seurat_obj_merge$kmeans_2 <- as.character(unlist(res))


#----------- Plot to compare to Lab images to ensure correct annotation of slices -----------

# plots
plot1 <- VlnPlot(seurat_obj_merge, features = "nCount_Spatial", pt.size = 0.1, group.by = "Sample_Name") + NoLegend()

ggsave("VlnPlot_nCount.png", plot1)


# Spatial plot of nCount
# Create the log-transformed metadata column if not already present
seurat_obj_merge$log10_nCount_Spatial <- log10(seurat_obj_merge$nCount_Spatial)

# Generate individual SpatialFeaturePlot objects (one per slice)
plots_list <- SpatialFeaturePlot(seurat_obj_merge, 
                                 features = "log10_nCount_Spatial", 
                                 combine = FALSE)

# If the returned list is not named, assign names from the seurat object images
if (is.null(names(plots_list))) {
  names(plots_list) <- names(seurat_obj_merge@images)
}

# Create a mapping from Sample_Name to Case_ID using the metadata
mapping <- unique(seurat_obj_merge@meta.data[, c("Sample_Name", "Case_ID")])
mapping_vec <- setNames(paste0(mapping$Sample_Name," - ",mapping$Case_ID), mapping$Sample_Name)

# Update each plot's title using the mapping and change title size
plots_list <- lapply(seq_along(plots_list), function(i) {
  nm <- names(plots_list)[i]
  plots_list[[i]] + 
    ggtitle(mapping_vec[[nm]]) +
    theme(plot.title = element_text(size = 8))
})


# Combine the updated plots into a 4x4 grid with one shared legend
combined_plot <- wrap_plots(plots_list, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

combined_plot

ggsave("SpatialFeaturePlot_nCount.png", combined_plot)





#----------- Run Seurat pre-processing pipeline -----------

# store mitochondrial percentage in object meta data
seurat_obj_merge <- PercentageFeatureSet(seurat_obj_merge, pattern = "^MT-", col.name = "percent.mt")

# Remove cells with zero nCount_RNA (UMIs)
seurat_obj_merge <- subset(seurat_obj_merge, subset = nCount_Spatial > 0)

# run seurat standard pipeline
seurat_obj_merge <- SCTransform(seurat_obj_merge, assay = "Spatial", verbose = FALSE)
DefaultAssay(seurat_obj_merge) <- "SCT"
seurat_obj_merge <- RunPCA(seurat_obj_merge, verbose = FALSE)
seurat_obj_merge <- FindNeighbors(seurat_obj_merge, dims = 1:30)
seurat_obj_merge <- FindClusters(seurat_obj_merge, verbose = FALSE)
seurat_obj_merge <- RunUMAP(seurat_obj_merge, dims = 1:30)

# save seurat object
SaveSeuratRds(seurat_obj_merge, "/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft/xenograft_seurat_250325.rds")


#----------- Plots -----------

plot_names <- unique(paste(seurat_obj_merge$Graft, seurat_obj_merge$Case_ID, seurat_obj_merge$Sample_Name, sep=" / "))

# nCounts
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_obj_merge, features = "nCount_Spatial", images = names(list_sr)[i],pt.size.factor = 2.5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("nCount_Spatial.png", arrange)


# mitochondrial
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_obj_merge, features = "percent.mt", images = names(list_sr)[i],pt.size.factor = 2.5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("percent.mt_Spatial.png", arrange)


# clusters
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialDimPlot(seurat_obj_merge, images = names(list_sr)[i],pt.size.factor = 2.5, image.alpha = 0) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(plot_names[i]) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("seuratclust_Spatial.png", arrange)

# k_means
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialDimPlot(seurat_obj_merge, group.by = "kmeans_2", images = names(list_sr)[i],pt.size.factor = 2.5,image.alpha = 0.5,
                                       cols = c("1" = 'red', "2" = 'grey')) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5)  +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("k2_Spatial.png", arrange)


#-------------- Harmony ------------
# run harmony
seurat_obj_merge <- RunHarmony(seurat_obj_merge,"run")

# dimensional reduction
seurat_obj_merge <- seurat_obj_merge %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 

# umap
seurat_obj_merge <- seurat_obj_merge %>% RunUMAP(reduction = "harmony",  dims = 1:30)

## run UMAP on total dataset harmony corrected
p1 <- DimPlot(seurat_obj_merge, reduction = "umap", group.by = "run", pt.size = .1) + ggtitle("Harmony UMAP - run")
p2 <- DimPlot(seurat_obj_merge, reduction = "umap", group.by = "grafted_area", label = TRUE) + ggtitle("Harmony UMAP - grafted_area")
p3 <- DimPlot(seurat_obj_merge, reduction = "umap", group.by = "Graft", label = TRUE) + ggtitle("Harmony UMAP - Graft")
#p4 <- DimPlot(seurat_obj_merge, reduction = "umap", label = TRUE,  pt.size = .1)
p4 <- FeaturePlot(seurat_obj_merge, features = "nCount_Spatial") + ggtitle("Harmony UMAP - nCount_Spatial")

arrange <- ggarrange(plotlist=list(p1,p2,p3,p4), nrow=2, ncol=2, widths = c(2,2))
ggsave("UMAP_harmony_total.png", arrange)

# save seurat object
SaveSeuratRds(seurat_obj_merge, "/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft/xenograft_seurat_250325.rds")



#----------- Highest expressed genes in k=2 of run 1 -----------

# read in data
seurat_obj_merge <- LoadSeuratRds("/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft/xenograft_seurat_250325.rds")

# subset for non-graft sites in run 1
seurat_sub <- subset(seurat_obj_merge, subset = kmeans_2 ==  2 & run == 1)

plot_names <- unique(paste(seurat_sub$Graft, seurat_sub$Case_ID, seurat_sub$Sample_Name, sep=" / "))
list_sr <- unique(seurat_sub$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_sub, features = "nCount_Spatial", images = list_sr[i],pt.size.factor = 2.5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("nCount_Spatial_run1_k2.png", arrange)

# get top genes by expression in k=2 (outside graft) from run 1
sum_exp <- rowSums(seurat_sub)
top_genes_k2 <- sum_exp[sum_exp > quantile(sum_exp, probs = seq(0.9, 0.99, by=0.01))[10]] # genes in 99% by expression


# plot top gene in k=2 by expression
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_sub, features = "GNAO1", images = list_sr[i],pt.size.factor = 2.5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("GNAO1_Spatial_run1_k2.png", arrange)


# plot top gene in k=2 by expression - total data
plot_names <- unique(paste(seurat_obj_merge$Graft, seurat_obj_merge$Case_ID, seurat_obj_merge$Sample_Name, sep=" / "))
list_sr <- unique(seurat_obj_merge$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_obj_merge, features = "GNAO1", images = list_sr[i],pt.size.factor = 2.5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("GNAO1_Spatial.png", arrange)

# calculate summed gene expression for top genes by expression in k=2 (outside graft) from run 1
seurat_obj_merge$top_gene_exp_run1.k2 <- colSums(scale(seurat_obj_merge@assays$SCT$counts)[names(top_genes_k2),])

hist(seurat_obj_merge$top_gene_exp_run1.k2)


# plot summed genes in k=2 by expression - total data
plot_names <- unique(paste(seurat_obj_merge$Graft, seurat_obj_merge$Case_ID, seurat_obj_merge$Sample_Name, sep=" / "))
list_sr <- unique(seurat_obj_merge$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_obj_merge, features = "top_gene_exp_run1.k2", images = list_sr[i],pt.size.factor = 2.5,image.alpha = 0.5) + 
    theme(legend.position = "none",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("TopRun1.k2_Spatial.png", arrange)


# k-means clustering k=2 on top run1 k2 genes expression on each Sample_Name independently
set.seed(123)
md <- as.data.frame(seurat_obj_merge@meta.data)
res <- list()

for (i in 1:length(unique(md$Sample_Name))) {
  tmp <- md[md$Sample_Name == unique(md$Sample_Name)[i],]
  kmeans_result <- kmeans(tmp[,"top_gene_exp_run1.k2"], centers=2)
  means <- tapply(tmp[,"top_gene_exp_run1.k2"], kmeans_result$cluster, mean)
  
  # check if we need to swap cluster labels
  if (means[1] < means[2]) {
    # swap cluster assignments
    kmeans_result$cluster <- ifelse(kmeans_result$cluster == 1, 2, 1)
  }
  res[[i]] <- kmeans_result$cluster
}

seurat_obj_merge$kmeans_2_top_run1.k2 <- as.character(unlist(res))


# Define graft as k=2 in top gene expression in run 1 k2, & k=1 in ncount
plot_names <- unique(paste(seurat_obj_merge$Graft, seurat_obj_merge$Case_ID, seurat_obj_merge$Sample_Name, sep=" / "))
list_sr <- unique(seurat_obj_merge$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialDimPlot(seurat_obj_merge, group.by = "kmeans_2_top_run1.k2", images = list_sr[i],pt.size.factor = 2.5,image.alpha = 0.5,
                                   cols = c("1" = 'red', "2" = 'grey')) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5)  +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("TopRun1.k2_kmeans_Spatial.png", arrange)


seurat_obj_merge$grafted_area <- seurat_obj_merge$kmeans_2_top_run1.k2 == "2" &  seurat_obj_merge$kmeans_2 == "1"


plot_names <- unique(paste(seurat_obj_merge$Graft, seurat_obj_merge$Case_ID, seurat_obj_merge$Sample_Name, sep=" / "))
list_sr <- unique(seurat_obj_merge$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialDimPlot(seurat_obj_merge, group.by = "grafted_area", images = list_sr[i],pt.size.factor = 2.5,image.alpha = 0.5,
                                   cols = c("TRUE" = 'red', "FALSE" = 'grey')) + 
    theme(legend.position = "none",legend.title = element_blank(),legend.key.size = unit(0.01, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5)  +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("Graft_Spatial.png", arrange)


#----------- Select human grafts and re-cluster -----------

p1 <- DimPlot(seurat_obj_merge, group.by = "run", label = TRUE) + ggtitle("UMAP - run")
p2 <- DimPlot(seurat_obj_merge, group.by = "grafted_area", label = TRUE) + ggtitle("UMAP - grafted_area")
p3 <- DimPlot(seurat_obj_merge, group.by = "Graft", label = TRUE) + ggtitle("UMAP - Graft")
p4 <- DimPlot(seurat_obj_merge, group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP - seurat_clusters")

arrange <- ggarrange(plotlist=list(p1,p2,p3), nrow=2, ncol=2, widths = c(2,2))
ggsave("UMAP_total.png", arrange)


# subset Cortex and re-cluster
seurat_cortex <- subset(seurat_obj_merge, subset = Graft == "Cortex")

# run harmony
seurat_cortex <- RunHarmony(seurat_cortex,"run")

# dimensional reduction
seurat_cortex <- seurat_cortex %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 

# umap
seurat_cortex <- seurat_cortex %>% RunUMAP(reduction = "harmony",  dims = 1:30)


# plot UMAPs 
p1 <- DimPlot(seurat_cortex, group.by = "run", label = TRUE) + ggtitle("Harmony UMAP - run")
p2 <- DimPlot(seurat_cortex, group.by = "grafted_area", label = TRUE) + ggtitle("UMAP - grafted_area")
p3 <- DimPlot(seurat_cortex, group.by = "Genotype", label = TRUE) + ggtitle("Harmony UMAP - Genotype")
p4 <- DimPlot(seurat_cortex, group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP - seurat_clusters")

arrange <- ggarrange(plotlist=list(p1,p2,p3,p4), nrow=2, ncol=2, widths = c(2,2))
ggsave("UMAP_cortex.png", arrange)

# table
tab_graft_cluster <- table(seurat_cortex$grafted_area,
      seurat_cortex$seurat_clusters)
frac_cluster_graft <- tab_graft_cluster[2,]/colSums(tab_graft_cluster)




library(ggplot2)
library(dplyr)
library(tidyr)


tab <- table(seurat_cortex$Genotype, seurat_cortex$seurat_clusters)
prop_tab <- prop.table(tab, margin = 1) * 100  # now each row sums to 100%
df_plot <- as.data.frame(prop_tab)
colnames(df_plot) <- c("Genotype", "Cluster", "Percentage")

# lot
ggplot(df_plot, aes(x = Genotype, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity") +
  ylab("Percentage of Cells") +
  xlab("Genotype") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  ggtitle("Stacked Barplot: % Cluster Composition per Genotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








# reassign grafts using cluster information
seurat_cortex$grafted_area[
  seurat_cortex$seurat_clusters %in% 
    names(frac_cluster_graft)[frac_cluster_graft > 0.5]] <- TRUE


p1 <- DimPlot(seurat_cortex, group.by = "run", label = TRUE) + ggtitle("Harmony UMAP - run")
p2 <- DimPlot(seurat_cortex, group.by = "grafted_area", label = TRUE) + ggtitle("UMAP - grafted_area")
p3 <- DimPlot(seurat_cortex, group.by = "Genotype", label = TRUE) + ggtitle("Harmony UMAP - Genotype")
p4 <- DimPlot(seurat_cortex, group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP - seurat_clusters")

arrange <- ggarrange(plotlist=list(p1,p2,p3,p4), nrow=2, ncol=2, widths = c(2,2))
ggsave("UMAP_cortex_clustassign.png", arrange)


#
plot_names <- unique(paste(seurat_cortex$Graft, seurat_cortex$Case_ID, seurat_cortex$Sample_Name, sep=" / "))
list_sr <- unique(seurat_cortex$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialDimPlot(seurat_cortex, group.by = "grafted_area", images = list_sr[i],pt.size.factor = 2.5,image.alpha = 0.5,
                                   cols = c("TRUE" = 'red', "FALSE" = 'grey')) + 
    theme(legend.position = "none",legend.title = element_blank(),legend.key.size = unit(0.01, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5)  +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=3, ncol=4, widths = c(2,2))
ggsave("Graft_Spatial_Cortex_clustassigned.png", arrange)







# subset grafts
seurat_graft <- subset(seurat_obj_merge, subset = grafted_area == "TRUE")

# run harmony
seurat_graft <- RunHarmony(seurat_graft,"run")

# dimensional reduction
seurat_graft <- seurat_graft %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 

# umap
seurat_graft <- seurat_graft %>% RunUMAP(reduction = "harmony",  dims = 1:30)


# plot UMAPs 
p1 <- DimPlot(seurat_graft, group.by = "run", label = TRUE) + ggtitle("Harmony UMAP - run")
p2 <- DimPlot(seurat_graft, group.by = "Graft", label = TRUE) + ggtitle("Harmony UMAP - Graft")
p3 <- DimPlot(seurat_graft, group.by = "Genotype", label = TRUE) + ggtitle("Harmony UMAP - Genotype")
p4 <- FeaturePlot(seurat_graft, features = "nCount_Spatial") + ggtitle("Harmony UMAP - nCount_Spatial")

arrange <- ggarrange(plotlist=list(p1,p2,p3,p4), nrow=2, ncol=2, widths = c(2,2))
ggsave("UMAP_graft.png", arrange)


# save seurat object
SaveSeuratRds(seurat_obj_merge, "/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft/xenograft_seurat_250325.rds")











# #----------- Technical reproducibility -----------
# 
# expression_data <- seurat_obj_merge@assays$SCT$data[,seurat_obj_merge$kmeans_2 == "1"]
# metadata <- seurat_obj_merge@meta.data[seurat_obj_merge$kmeans_2 == "1",]
# 
# # Convert data to data.table format
# dt_expression <- as.data.table(expression_data)
# setnames(dt_expression, colnames(expression_data))
# dt_metadata <- as.data.table(metadata)
# 
# # Add row identifiers for joining
# dt_expression$Gene_ID <- rownames(expression_data)
# 
# # Melt the expression data to long format
# dt_long_expression <- melt(dt_expression, id.vars = "Gene_ID", variable.name = "Sample_Name", value.name = "Expression")
# 
# # Join metadata to get Sample_Name information
# dt_long_expression <- merge(dt_long_expression, dt_metadata, by = "Sample_Name")
# 
# # Calculate mean expression per gene per Sample_Name
# result <- dt_long_expression[, .(Mean_Expression = mean(Expression)), by = .(Gene_ID, Sample_Name)]
# 
# # Optionally, cast back to wide format to have genes as rows and Sample_Names as columns
# dt_wide <- dcast(result, Gene_ID ~ Sample_Name, value.var = "Mean_Expression")
# 
# # Print the result
# print(dt_wide)
# 
# 
# 
# 
# 
# dt <- rbind(exp_dat)
# 
# dt2 <- exp_dat[lapply(.SD,mean), by=meta_dat$Sample_Name,]
# 
# # Melt the data frame to long format
# library(reshape2)
# long_exp_dat <- melt(exp_dat, id.vars="Gene")
# 
# # Join with metadata to include Sample_Name information
# long_exp_dat <- merge(long_exp_dat, metadata, by.x="variable", by.y="Sample")
# 
# # Calculate mean expression per gene per Sample_Name
# mean_expression_per_gene_per_Sample_Name <- long_exp_dat %>%
#   group_by(Sample_Name, Gene) %>%
#   summarize(mean_expression = mean(value), .groups = 'drop')
# 
# # Convert result back to a wide format if needed
# wide_expression_per_Sample_Name <- reshape2::dcast(mean_expression_per_gene_per_Sample_Name, Gene ~ Sample_Name, value.var="mean_expression")
# 
# 
# 
# mean_expression_per_Sample_Name <- colMeans(exp_dat) %>% 
#   data.frame(mean_expression=.) %>% 
#   cbind(meta_dat) %>%
#   group_by(Sample_Name) %>%
#   summarize(mean_expression=mean(mean_expression))
# 
# 
# Sample_Name
# 
# Case_ID
# 
# Graft