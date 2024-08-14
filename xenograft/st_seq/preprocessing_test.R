# libs
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)

#----------- metadata -----------
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft")
meta_dat <- read_excel("/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft/meta_data_test.xlsx")

#----------- spaceranger out directories -----------

# directories of spaceranger "outs"
S11W1_LB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_LB/spaceranger_count_cytassist/S11W1_LB/outs"
S11W1_LT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_LT/spaceranger_count_cytassist/S11W1_LT/outs"
S11W1_RB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_RB/outs"
S11W1_RT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S11W1_RT/outs"

S12W2_LB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/SW12W2_LB/outs"
S12W2_LT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S12W2_LT/outs"
S12W2_RB <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S12W2_RB/outs"
S12W2_RT <- "/Volumes/research-data/PRJ-ASAPbioin/10x_visium_data/mouse_xenografts/split_data/S12W2_RT/outs"

# create lists
list_sr <- c(S11W1_LB,S11W1_LT,S11W1_RB,S11W1_RT,
             S12W2_LB,S12W2_LT,S12W2_RB,S12W2_RT)
names(list_sr) <- c("S11W1_LB","S11W1_LT","S11W1_RB","S11W1_RT",
                    "S12W2_LB","S12W2_LT","S12W2_RB","S12W2_RT")

#----------- create Seurat objects and merge -----------

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
                  add.cell.ids = names(seurat_object_list), project = "xenograft_test")

# add metadata
md <- as.data.frame(seurat_obj_merge@meta.data)
md$row_id <- row.names(md)
tmp <- merge(md, meta_dat, by="Sample_Name", all.x = TRUE)
row.names(tmp) <- tmp$row_id
tmp <- tmp[row.names(seurat_obj_merge@meta.data),]
seurat_obj_merge@meta.data <- tmp

# k-means clustering k=2 on each Sample_Name independentlty
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

# # plots
# plot1 <- VlnPlot(seurat_obj_merge, features = "nCount_Spatial", pt.size = 0.1, group.by = "Sample_Name") + NoLegend()
# plot2 <- SpatialFeaturePlot(seurat_obj_merge, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)


#----------- Run Seurat pre-processing pipeline -----------

seurat_obj_merge <- SCTransform(seurat_obj_merge, assay = "Spatial", verbose = FALSE)
DefaultAssay(seurat_obj_merge) <- "SCT"
seurat_obj_merge <- RunPCA(seurat_obj_merge, verbose = FALSE)
seurat_obj_merge <- FindNeighbors(seurat_obj_merge, dims = 1:30)
seurat_obj_merge <- FindClusters(seurat_obj_merge, verbose = FALSE)
seurat_obj_merge <- RunUMAP(seurat_obj_merge, dims = 1:30)


#----------- Plots -----------
plot_names <- unique(paste(seurat_obj_merge$Graft, seurat_obj_merge$Case_ID, seurat_obj_merge$Sample_Name, sep=" / "))

# nCounts
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_obj_merge, features = "nCount_Spatial", images = names(list_sr)[i],pt.size.factor = 2.5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) 
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=2, widths = c(2,2))
ggsave("nCount_Spatial.png", arrange)

# clusters
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialDimPlot(seurat_obj_merge, images = names(list_sr)[i],pt.size.factor = 2.5, image.alpha = 0) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(plot_names[i]) 
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=2, widths = c(2,2))
ggsave("seuratclust_Spatial.png", arrange)

# k_means
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialDimPlot(seurat_obj_merge, group.by = "kmeans_2", images = names(list_sr)[i],pt.size.factor = 2.5,image.alpha = 0.5,
                                       cols = c("1" = 'red', "2" = 'grey')) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) 
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=2, widths = c(2,2))
ggsave("k2_Spatial.png", arrange)




#----------- Technical reproducibility -----------

expression_data <- seurat_obj_merge@assays$SCT$data[,seurat_obj_merge$kmeans_2 == "1"]
metadata <- seurat_obj_merge@meta.data[seurat_obj_merge$kmeans_2 == "1",]

# Convert data to data.table format
dt_expression <- as.data.table(expression_data)
setnames(dt_expression, colnames(expression_data))
dt_metadata <- as.data.table(metadata)

# Add row identifiers for joining
dt_expression$Gene_ID <- rownames(expression_data)

# Melt the expression data to long format
dt_long_expression <- melt(dt_expression, id.vars = "Gene_ID", variable.name = "Sample_Name", value.name = "Expression")

# Join metadata to get Sample_Name information
dt_long_expression <- merge(dt_long_expression, dt_metadata, by = "Sample_Name")

# Calculate mean expression per gene per Sample_Name
result <- dt_long_expression[, .(Mean_Expression = mean(Expression)), by = .(Gene_ID, Sample_Name)]

# Optionally, cast back to wide format to have genes as rows and Sample_Names as columns
dt_wide <- dcast(result, Gene_ID ~ Sample_Name, value.var = "Mean_Expression")

# Print the result
print(dt_wide)





dt <- rbind(exp_dat)

dt2 <- exp_dat[lapply(.SD,mean), by=meta_dat$Sample_Name,]

# Melt the data frame to long format
library(reshape2)
long_exp_dat <- melt(exp_dat, id.vars="Gene")

# Join with metadata to include Sample_Name information
long_exp_dat <- merge(long_exp_dat, metadata, by.x="variable", by.y="Sample")

# Calculate mean expression per gene per Sample_Name
mean_expression_per_gene_per_Sample_Name <- long_exp_dat %>%
  group_by(Sample_Name, Gene) %>%
  summarize(mean_expression = mean(value), .groups = 'drop')

# Convert result back to a wide format if needed
wide_expression_per_Sample_Name <- reshape2::dcast(mean_expression_per_gene_per_Sample_Name, Gene ~ Sample_Name, value.var="mean_expression")



mean_expression_per_Sample_Name <- colMeans(exp_dat) %>% 
  data.frame(mean_expression=.) %>% 
  cbind(meta_dat) %>%
  group_by(Sample_Name) %>%
  summarize(mean_expression=mean(mean_expression))


Sample_Name

Case_ID

Graft