# scripts for running RCTD and CSIDE method to identify ctDEG in Visium data

library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(plyr)
library(data.table)
library(msigdbr)
#library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(EnhancedVolcano)
library(viridis)
library(rstatix)
library(ComplexHeatmap)
library(readxl)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(semla)
library(tibble)
library(patchwork)
library(dplyr)
library(reshape2)
library(tidyr)


# source local
source("/Users/zacc/USyd/spatial_transcriptomics/reports/Abaurre_DA_manuscript/scripts/convenience.R")
source("/Users/zacc/github_repo/spacexr/R/CSIDE_plots.R")

register(SerialParam())


#------------ Inputs ------------

results_folder = "/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/vis_280125"
setwd(results_folder)

## ST data
# visium original (n=2)
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/SN_210823/SN_210823_merged.seurat.Rdata")
st_obj <- subset(x = brain.merge, subset = run_name == "V52Y16-079-A1" & tissue == "SN" )

# # extract counts, run Neuroestimator (gbhpc) and rebind
# counts <- st_obj@assays$Spatial@counts
# write.table(counts, file="/Users/zacc/USyd/spatial_transcriptomics/reports/Abaurre_DA_manuscript/V52Y16_079_A1_countmmatrix_121224.txt", sep="\t")


##  single-cell data
## merged dataset ##
load("/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/merged_kam.sil.web_seurat.Rdata")

#------------ pre-processing single-cell reference ------------

# i) taking excitatory, inhibitory and non-neurons from Kamath and DA neuron populations from Siletti re-code
toMatch <- c("Ex_","Inh_","Astro_","Endo_","Ependyma_","Macro_","MG_","Olig_","OPC_")
group1 <- unique(merge.combined.sct@meta.data$cell_type_merge)[grep(paste(toMatch,collapse="|"), unique(merge.combined.sct@meta.data$cell_type_merge))]

da_neurons <- names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Siletti"]))
group_total <- which(merge.combined.sct@meta.data$cell_type_merge %in% c(group1,da_neurons))
sc_obj <- merge.combined.sct[,group_total]

#  - save count matrix alone for Neuroestimator analysis -  run using Ubuntu with neuroestimator_run.R scripts 
counts <- sc_obj@assays$SCT$counts
write.table(counts, file= "abaurre_v1_countmmatrix.txt", sep="\t")

# ii) remove cell-types with <25 cells
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]

# iii) remove CALB1_NPW_ONECUT1 as from Periaqueductal gray matter
sc_obj <- sc_obj[,!sc_obj@meta.data$cell_type_merge %in% "CALB1_NPW_ONECUT1"]

# iv) get count data
DefaultAssay(sc_obj) <- "SCT"
counts_sc <- sc_obj[["SCT"]]$counts

# v) set cell types and quant nUMI and create reference
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
reference <- Reference(counts_sc, cell_types, nUMIsc)



#------------ pre-processing spatial dataset ------------

# i) load in counts matrix
counts <- st_obj@assays$Spatial@counts

# ii) load in coordinates
offset = 500
kth_image <- st_obj@images
num_elements <- length(kth_image)
coordinates_list <- vector("list", length = num_elements)
for (i in seq_along(kth_image)) {
  if ("coordinates" %in% slotNames(kth_image[[i]])) {
    coordinates_list[[i]] <- kth_image[[i]]@coordinates
    if (i > 1){
      coordinates_list[[i]]["imagecol"] <-  coordinates_list[[i]]["imagecol"] + max(coordinates_list[[i-1]]["imagecol"]) + offset
    }
  }
}
combined_coordinates <- do.call(rbind, coordinates_list)
coords = combined_coordinates[,c("imagecol","imagerow")]

# iii) inputs for RCTD
nUMI <- colSums(counts) # In this case, total counts per spot
puck <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 



#------------- Run RCTD to obtain cell proportions ------------- 

set.seed(100)

# i) create RCTD object
myRCTD <- create.RCTD(puck, reference, max_cores = 8)

# ii) re-run custom create to extract the DEG info, rank and get top 10 for each cell-type, write to file
myRCTD_logFC <- create.RCTD_logFC(puck, reference, max_cores = 8)

myRCTD_logFC_top10 <- lapply(myRCTD_logFC, function(x){
  x <- as.data.frame(x)
  tmp <- x[order(-as.numeric(x$logFC)),]
  tmp <- tmp[1:10,]
  return(tmp)
})

rctd_10 <- do.call(rbind,myRCTD_logFC_top10)
write.table(rctd_10,file="rctd_deg_top10.txt",quote = FALSE,sep = "\t", row.names = FALSE)

# iii) run RCTD
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# iv) add data to ST metadata
st_obj <- AddMetaData(st_obj,normalize_weights(myRCTD@results$weights))

# v) save objects
save(st_obj, myRCTD, file = "ST_objects_Abaurre_v1.RData")



#------ spatial heatmaps of DA enrichments (weights) -----------
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/vis_280125/ST_objects_Abaurre_v1.RData")

# set plot parameters
ratio = 1.9
meta_to_plot <- da_neurons[! da_neurons %in% c("CALB1_NPW_ONECUT1") ] # remove CALB1_NPW_ONECUT1 as from Periaqueductal gray matter

# loop through to plot enrichment of each cell-type
for (i in 1:length(meta_to_plot)){
  plot2 <- SpatialFeaturePlot(st_obj, features = meta_to_plot[i], images = "slice1",pt.size.factor = 3000,image.alpha = 0, 
  ) + 
    theme(legend.position = "right", aspect.ratio = ratio,legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(meta_to_plot[i])+ geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +
    theme(
      title = element_text(colour = "white",size=10),
      plot.background=element_rect(fill = "black"),
      panel.background = element_rect(fill = 'black'),
      
      # Change legend 
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "black", color = NA),
      legend.text=element_text(color="white",size=12,angle = 270)
    ) 
  
  ggsave(plot2,file=paste0("spatplot_rctd_weights_",meta_to_plot[i],"_black.png")) 
}



#------ spatial plot of DA neuron (winner-takes_all) -----------

# create scaled data for each DA cell-type and assign highest value as label of spot
da_weights <- as.data.frame(st_obj@meta.data)
da_weights_scale <- scale(da_weights[,colnames(da_weights) %in% da_neurons])

max_colnames_vector <- apply(da_weights_scale, 1, function(row) {
  colnames(da_weights_scale)[which.max(row)]
})

st_obj@meta.data$max_da_neuron <- max_colnames_vector

# plot spatial using cryptograph colors
plot2 <- SpatialDimPlot(
  st_obj, 
  "max_da_neuron", 
  images = "slice1", 
  pt.size.factor = 3000, 
  image.alpha = 0, 
  cols = DA_col_cytograph
) +
  theme(
    legend.position = "right",
    aspect.ratio = ratio,
    legend.title = element_blank(),
    legend.key.size = unit(0.8, 'cm'),
    # Title and background
    title = element_text(colour = "white", size = 10),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = 'black'),
    # Legend styling
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white", size = 6),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  ggtitle("DA cell-types") +
  geom_point(shape = 21, stroke = NA, fill = NA, size = 0.5)

# Save the plot
ggsave(plot2, file = "spatplot_rctd_wta_da_cytograph.png", bg = "black")


#------ spatial plot of DA neuron (winner-takes_all, TH mask) -----------

st_obj@meta.data$th_mask <- 1
st_obj@meta.data$th_mask[st_obj@assays$SCT@counts["TH",] < 1] <- -1

# Recreate the plot using ggplot2, setting alpha based on the new column
plot_data <- ggplot_build(plot2)
plot_df <- plot_data$data[[1]]
plot_df$alpha <- ifelse(st_obj@meta.data$th_mask[row.names(st_obj@meta.data) %in% row.names(st_obj@images$slice1@coordinates)] == -1, 0, 1) # Modify as needed

plot_df$alpha[st_obj@meta.data$th_mask[st_obj@assays$SCT@counts["TH",] > 0]] <- 1

# Create masked spatial plot
plot3 <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(aes(color = fill, alpha = alpha), size = 1) + 
  scale_color_identity() +  
  scale_alpha_identity() +  
  theme_minimal() +
  coord_fixed() +
  theme(
    legend.position = "right",
    aspect.ratio = ratio,
    legend.title = element_blank(),
    legend.key.size = unit(0.8, 'cm'),
    title = element_text(colour = "white", size = 10),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = 'black'),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white", size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  ggtitle("DA cell-types (TH+ masked)") +
  geom_point(shape = 21, stroke = NA, fill = NA, size = 0.5)

# Save the final plot
ggsave(plot3, file = "spatplot_rctd_wta_da.thmask_cytograph.png", bg = "black")



#------------ Add neuroestimator and manual brain region annotations ------------

# load neuroestimator results
load("/Users/zacc/USyd/spatial_transcriptomics/reports/Abaurre_DA_manuscript/neuroestimator_res_V52Y16_079_A1_121224.Rdata")
res <- as.data.frame(res)
row.names(res) <- gsub("\\.","-",row.names(res))
common_names <- intersect(row.names(res), colnames(st_obj))

# add neuroestimator to meta
st_obj <- st_obj[,common_names]
res <- res[common_names,]

st_obj$predicted_activity <- res

# add manual annotated regions with border labels to metadata
st_obj_new <- LoadSeuratRds("/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/processed_midbrain_sections/semla_mb_merged_hr.rds")

st_obj_new$region_annot[st_obj_new$region_annot == "region"] <- "WM"

regions_annot <- unique(st_obj_new$region_annot)
for (i in 1:length(regions_annot)){
  print(i)
  st_obj_new <- RegionNeighbors(st_obj_new, 
                                column_name = "region_annot", 
                                column_labels = regions_annot[i])
}

# get neighbouring spots to each region for plotting
meta1 <- st_obj_new@meta.data[st_obj_new@meta.data$sample_id == "case34.1",]
row.names(meta1) <- paste0(row.names(meta1),"_1")

common_names <- intersect(row.names(meta1), colnames(st_obj))

st_obj <- st_obj[,common_names]
meta1 <- meta1[common_names,]

st_obj$region_annot <- meta1$region_annot
st_obj$nb_to_WM <- meta1$nb_to_WM
st_obj$nb_to_SNV <- meta1$nb_to_SNV
st_obj$nb_to_SND <- meta1$nb_to_SND
st_obj$nb_to_RN <- meta1$nb_to_RN
st_obj$nb_to_VTA <- meta1$nb_to_VTA
st_obj$nb_to_SNM <- meta1$nb_to_SNM
st_obj$nb_to_SNL <- meta1$nb_to_SNL


# save
SaveSeuratRds(st_obj,file = "ST_visium_Abaurre_v1.rds")


#------------ Plot Neuroestimator ST ------------

st_obj <- LoadSeuratRds("ST_visium_Abaurre_v1.rds")

# set plot parameters
ratio = 1.9
meta_to_plot <- c("predicted_activity")
  
# loop through to plot enrichment of each cell-type
for (i in 1:length(meta_to_plot)){
  plot2 <- SpatialFeaturePlot(st_obj, features = meta_to_plot[i], images = "slice1",pt.size.factor = 3000,image.alpha = 0, 
  ) + 
    theme(legend.position = "right", aspect.ratio = ratio,legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(meta_to_plot[i])+ geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +
    theme(
      title = element_text(colour = "white",size=10),
      plot.background=element_rect(fill = "black"),
      panel.background = element_rect(fill = 'black'),
      
      # Change legend 
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "black", color = NA),
      legend.text=element_text(color="white",size=12,angle = 270)
    ) 
  
  ggsave(plot2,file=paste0("spatplot_",meta_to_plot[i],"_black.png")) 
}

#------------ Plot Gnee Expression interest in ST ------------

# set plot parameters
ratio = 1.9
meta_to_plot <- c("ALDH1A1","ANXA1","SNCA","VCAN","NDNF")

# loop through to plot enrichment of each cell-type
for (i in 1:length(meta_to_plot)){
  plot2 <- SpatialPlot(st_obj, features = meta_to_plot[i], images = "slice1",pt.size.factor = 3000,image.alpha = 0, 
  ) + 
    theme(legend.position = "right", aspect.ratio = ratio,legend.title = element_blank(),legend.key.size = unit(0.8, 'cm')) + 
    ggtitle(meta_to_plot[i])+ geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +
    theme(
      title = element_text(colour = "white",size=10),
      plot.background=element_rect(fill = "black"),
      panel.background = element_rect(fill = 'black'),
      
      # Change legend 
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "black", color = NA),
      legend.text=element_text(color="white",size=12,angle = 270)
    ) 
  
  ggsave(plot2,file=paste0("spatplot_",meta_to_plot[i],"_black.png")) 
}




#------------ Recipprical KS test to ID location ------------

# function for reciprical KS  test
ReciprocalKSTest <- function(obj, features, group.by, pval_cutoff = 0.05) {
  # Extract data
  data <- FetchData(obj, vars = c(group.by, features))
  colnames(data)[1] <- "group"  # Standardize column name for group.by
  
  # Convert features to long format
  long_data <- data %>%
    pivot_longer(cols = features, names_to = "feature", values_to = "value")
  
  # Prepare table to store results
  results <- list()
  
  for (feature in unique(long_data$feature)) {
    feature_data <- long_data %>% filter(feature == !!feature)  # Filter data for this feature
    
    # Calculate mean values for each group
    group_means <- feature_data %>%
      group_by(group) %>%
      summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_value))  # Order groups by mean values
    
    # Perform reciprocal KS tests
    ks_results <- data.frame()
    for (i in 1:(nrow(group_means) - 1)) {
      # Define the two sets to compare
      set1 <- feature_data %>% filter(group %in% group_means$group[1:i]) %>% pull(value)
      set2 <- feature_data %>% filter(group == group_means$group[i + 1]) %>% pull(value)
      
      # Perform the KS test
      if (length(set1) > 1 & length(set2) > 1) {
        p_value <- ks.test(set1, set2)$p.value
      } else {
        p_value <- NA  # Handle insufficient data
      }
      
      # Store the result
      ks_results <- rbind(
        ks_results,
        data.frame(
          Feature = feature,
          Group_Comparison = paste0(
            paste(group_means$group[1:i], collapse = "+"), " vs ", group_means$group[i + 1]
          ),
          P_Value = p_value
        )
      )
    }
    
    # Identify the first significant p-value
    first_significant <- ks_results %>%
      filter(P_Value < pval_cutoff) %>%
      slice(1)  # First row with p-value < cutoff
    
    # Add group means to the results
    group_means <- group_means %>%
      mutate(Feature = feature) %>%
      rename(Group = group, Mean_Value = mean_value)
    
    # Combine the results for this feature
    feature_results <- merge(ks_results, group_means, by = "Feature", all = TRUE)
    feature_results$First_Significant <- ifelse(
      feature_results$Group_Comparison == first_significant$Group_Comparison,
      TRUE, FALSE
    )
    results[[feature]] <- feature_results
  }
  
  # Combine all features into a single table
  final_results <- do.call(rbind, results)
  
  # Return the final table
  return(final_results)
}

# function for dotplot
DotPlotSingleDotOrdered <- function(obj, features, group.by, pval_results, pval_cutoff = 0.05) {
  # Extract data
  data <- FetchData(obj, vars = c(group.by, features))
  colnames(data)[1] <- "group"  # Standardize column name for group.by
  
  # Convert features to long format
  long_data <- data %>%
    pivot_longer(cols = features, names_to = "feature", values_to = "value")
  
  # Calculate mean values for each group and feature
  mean_data <- long_data %>%
    group_by(group, feature) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  
  # Process significant_regions by splitting Group_Comparison and keeping relevant matches
  significant_regions <- pval_results %>%
    mutate(Highlight_Regions = sapply(strsplit(Group_Comparison, " vs "), function(x) x[1])) %>%
    separate_rows(Highlight_Regions, sep = "\\+") %>%  # Split regions into separate rows
    dplyr::select(Feature, Highlight_Regions) %>%
    rename(feature = Feature, group = Highlight_Regions)
  
  # Add a column to mark highlights based on exact matches of group and feature 
  highlights <- mean_data
  highlights$highlight <- paste0(highlights$group,highlights$feature) %in% paste0(significant_regions$group,significant_regions$feature)
  
  # Reorder columns (groups) by their mean values
  highlights <- highlights %>%
    mutate(group = factor(group, levels = unique(group[order(mean_value, decreasing = TRUE)])))
  
  # Reorder rows (features) by their mean values
  highlights <- highlights %>%
    mutate(feature = factor(feature, levels = unique(feature[order(mean_value, decreasing = TRUE)])))
  
  # Create the DotPlot
  p <- ggplot(highlights, aes(x = group, y = feature, size = mean_value, fill = highlight)) +
    geom_point(shape = 21, color = "black", stroke = 0.5) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey"), name = "Enriched Region") +
    scale_size_continuous(range = c(0.5, 15), name = "Mean Proportion") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_blank(), # Remove y-axis title
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(),
      legend.position = "right"
    )
  
  return(p)
}

# Perform Reciprocal KS Tests
pval_results <- ReciprocalKSTest(obj = st_obj, features = da_neurons, group.by = "region_annot", pval_cutoff = 0.05)

# Select regions from KS test p < 0.05
pval_results <- pval_results[pval_results$P_Value < 0.05,]
pval_results <- pval_results[!duplicated(pval_results$Feature),]

# Save results as a table
write.csv(pval_results, "reciprocal_ks_results.csv", row.names = FALSE)

# Create DotPlot with Highlighted Significant Regions
pdf("dotplot_highlighted_significant.pdf", width = 12, height = 6)
DotPlotSingleDotOrdered(obj = st_obj, features = da_neurons, group.by = "region_annot", pval_results = pval_results, pval_cutoff = 0.05)
dev.off()


































#------------ Dot plot of brain region enrichment for each cell-type ------------

# Remove CALB1_NPW_ONECUT1
da_neurons <- da_neurons[! da_neurons %in% c("CALB1_NPW_ONECUT1") ] # remove CALB1_NPW_ONECUT1 as from Periaqueductal gray matter

# Function to create a DotPlot without overlapping dots
DotPlotSingleDot <- function(obj, features, group.by, split.by = NULL) {
  # Extract data
  data <- FetchData(obj, vars = c(group.by, features))
  colnames(data)[1] <- "group"  # Standardize column name for group.by
  
  # Convert features to long format
  long_data <- data %>%
    pivot_longer(cols = features, names_to = "feature", values_to = "value") %>%
    group_by(group, feature) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")  # Aggregate values
  
  # Highlight the group with the highest mean value for each feature
  highlights <- long_data %>%
    group_by(feature) %>%
    mutate(highlight = mean_value == max(mean_value)) %>%
    ungroup()
  
  # Create DotPlot
  p <- ggplot(highlights, aes(x = group, y = feature, size = mean_value, fill = highlight)) +
    geom_point(shape = 21, color = "black", stroke = 0.5) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) + # Highlight in red
    scale_size_continuous(range = c(0.5, 15)) + # Adjust dot size range
    theme_minimal() +
    xlab(group.by) + ylab("Feature") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  return(p)
}

# Example usage
pdf("dotplot_abaurre_manregion.pdf", width = 12, height = 6)
DotPlotSingleDot(obj = st_obj, features = da_neurons, group.by = "region_annot")
dev.off()









DotPlotSingleDotOrdered <- function(obj, features, group.by, pval_cutoff = 0.05) {
  # Extract data
  data <- FetchData(obj, vars = c(group.by, features))
  colnames(data)[1] <- "group"  # Standardize column name for group.by
  
  # Convert features to long format
  long_data <- data %>%
    pivot_longer(cols = features, names_to = "feature", values_to = "value")
  
  # Perform KS test for each feature
  pval_data <- long_data %>%
    group_by(feature) %>%
    summarize(
      # Compare each group to all others using KS test
      p_value = tryCatch({
        result <- ks.test(value[group == unique(group)[1]], value[group != unique(group)[1]])
        print(paste("Feature:", unique(feature), "- p-value:", result$p.value)) # Print p-value
        result$p.value
      }, error = function(e) NA), # Handle cases with insufficient data
      .groups = "drop"
    )
  
  # Calculate mean values and merge with p-values
  mean_data <- long_data %>%
    group_by(group, feature) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  
  highlights <- mean_data %>%
    left_join(pval_data, by = "feature") %>%
    group_by(feature) %>%
    mutate(
      highlight = ifelse(p_value < pval_cutoff & mean_value == max(mean_value), TRUE, FALSE)
    ) %>%
    ungroup()
  
  # Reorder columns (groups) by their mean values
  highlights <- highlights %>%
    mutate(group = factor(group, levels = unique(group[order(mean_value, decreasing = TRUE)])))

  # Reorder rows (features) by their mean values
  highlights <- highlights %>%
    mutate(feature = factor(feature, levels = unique(feature[order(mean_value, decreasing = TRUE)])))
  
  # Create DotPlot
  p <- ggplot(highlights, aes(x = group, y = feature, size = mean_value, fill = highlight)) +
    geom_point(shape = 21, color = "black", stroke = 0.5) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey"), name = "Enriched Region") +
    scale_size_continuous(range = c(0.5, 15), name = "Mean Proportion") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_blank(), # Remove y-axis title
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(),
      legend.position = "right"
    )
  
  return(p)
}

# Example usage
pdf("dotplot_ordered_highlighted.pdf", width = 12, height = 6)
DotPlotSingleDotOrdered(obj = st_obj, features = da_neurons, group.by = "region_annot", pval_cutoff = 0.05)
dev.off()



