# scripts for running RCTD and CSIDE method to identify ctDEG within GeoMx datasets

library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(plyr)
library(data.table)
library(msigdbr)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(EnhancedVolcano)
library(viridis)
library(rstatix)
library(ComplexHeatmap)
library(readxl)
library(dplyr)
library(nnet)
library(tidyr)
library(emmeans)
source("/Users/zacc/github_repo/spacexr/R/CSIDE_plots.R")
register(SerialParam())


#----------- Inputs -----------

analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geeomx_280125"
setwd(analysis_dir)

rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
setwd(analysis_dir)

##  single-cell data
load("/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/merged_kam.sil.web_seurat.Rdata")

#------------ pre-processing single-cell reference ------------

# i) taking excitatory, inhibitory and non-neurons from Kamath and DA neuron populations from Siletti re-code
toMatch <- c("Ex_","Inh_","Astro_","Endo_","Ependyma_","Macro_","MG_","Olig_","OPC_")
group1 <- unique(merge.combined.sct@meta.data$cell_type_merge)[grep(paste(toMatch,collapse="|"), unique(merge.combined.sct@meta.data$cell_type_merge))]

da_neurons <- names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Siletti"]))
group_total <- which(merge.combined.sct@meta.data$cell_type_merge %in% c(group1,da_neurons))
sc_obj <- merge.combined.sct[,group_total]

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

# vi) create reference with just DA neurons
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% da_neurons]
counts_sc <- sc_obj[["SCT"]]$counts

cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
reference_da <- Reference(counts_sc, cell_types, nUMIsc)


#-----------Format geomx spatial data -----------

# load geomx normalised data
load(rdata)

gxd# select samples
#keep_index <- gxdat_s$Diagnosis == "CTR"
keep_index <- gxdat_s$Diagnosis != "NTC"

# load in counts matrix
counts <- gxdat_s@assays$RNA@counts[,keep_index]  # load in counts matrix

# create metadata matrix 
meta <- gxdat_s@meta.data[keep_index,]

# confirm names
table(colnames(counts) == row.names(meta))

# load in coordinates
coords = as.data.frame(gxdat_s@reductions$umap@cell.embeddings[keep_index,]) # we are using UMAP coordinates as dummy variables until we obtain DSP ROI locations
colnames(coords) <- c("imagerow","imagecol")

# process counts
nUMI <- colSums(counts) # In this case, total counts per spot

# construct ST to deconvolute
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')



#------------ RCTD to obtain cell proportions ------------

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add data to ST metadata
gxdat_s@meta.data <- cbind(gxdat_s@meta.data,normalize_weights(myRCTD@results$weights))

# save total object
gxdat_s <- UpdateSeuratObject(gxdat_s)
save(gxdat_s, myRCTD, file = "ST_objects_Abaurre_gx_v1.RData")

# format object for sharing
gxdat_s_ctr <- subset(gxdat_s, subset = Diagnosis == "CTR") # subset controls
gxdat_s_ctr <- subset(gxdat_s_ctr, subset = ROI != "LC") # subset LC
gx_id <- read_excel("/Users/zacc/USyd/ASAP/CRN_upload/GeoMx cohort.xlsx", sheet = "zac") # recode id's
gxdat_s_ctr$MJFF_number <- mapvalues(gxdat_s_ctr$Brainbank_ID,gx_id$`ID-Never share this anyone else`,gx_id$`MJFF number`) # recode

remove_exact <- c("Brainbank_ID", "Brainbank_ID_sub", "Brainbank") # Define the exact columns to remove
start_col <- which(colnames(gxdat_s_ctr@meta.data) == "ROI Area [µm²].iNM") # indice start for NM
end_col <- which(colnames(gxdat_s_ctr@meta.data) == "NM_mean.optical.density.eNM") # indice end for NM

if (length(start_col) == 1 & length(end_col) == 1) {
  remove_range <- colnames(gxdat_s_ctr@meta.data)[start_col:end_col]
} else {
  remove_range <- c()  # Ensure this is empty if columns not found
}
remove_cols <- c(remove_exact, remove_range) # Combine all columns to remove
gxdat_s_ctr@meta.data <- gxdat_s_ctr@meta.data[, !colnames(gxdat_s_ctr@meta.data) %in% remove_cols] # Remove the specified columns from meta.data

# save total object
save(gxdat_s_ctr, myRCTD, file = "ST_objects_Abaurre_gx_ctr_v1.RData")



## Alternative  - run RCTD with only DA
##  run RCTD
myRCTD <- create.RCTD(puck, reference_da, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add data to ST metadata
meta_da <- cbind(gxdat_s@meta.data,normalize_weights(myRCTD@results$weights))




#------------ Multinomial Logistic Regression + Tukey’s Test to identify enriched location ------------
load("ST_objects_Abaurre_gx_ctr_v1.RData")

# define DA neurons
da_neurons <- da_neurons[!da_neurons %in% "CALB1_NPW_ONECUT1"]

# function to reciprocally perform Multinomial Logistic Regression + Tukey’s Test on highest (mean) ROI #1 v else, ROI #1 + ROI #2 v else etc.
ReciprocalTukeyTest <- function(obj, features, group.by, pval_cutoff = 0.05) {
  # Extract data
  data <- FetchData(obj, vars = c(group.by, features, "Age", "Sex", "PMD.hs"))
  colnames(data)[1] <- "group"  # Standardize column name for group.by
  
  # Convert features to long format
  long_data <- data %>%
    pivot_longer(cols = features, names_to = "feature", values_to = "value")
  
  # Ensure categorical variables are factors, numeric variables are numeric
  long_data$group <- as.factor(long_data$group)
  long_data$Sex <- as.factor(long_data$Sex)
  long_data$Age <- as.numeric(long_data$Age)
  long_data$PMD.hs <- as.numeric(long_data$PMD.hs)
  
  # Prepare table to store results
  results <- list()
  
  for (feature in unique(long_data$feature)) {
    print(feature)
    feature_data <- long_data %>% filter(feature == !!feature)
    
    # Calculate mean values for each group (ROI)
    group_means <- feature_data %>%
      group_by(group) %>%
      dplyr::summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_value))
    
    # Initialize storage for pairwise cumulative comparisons
    tukey_rec_results <- data.frame()
    
    for (i in 1:(nrow(group_means) - 1)) {
      print(i)
      # Define cumulative groups for pairwise comparison
      included_groups <- group_means$group[1:i]  # First i groups combined
      next_group <- group_means$group[i + 1]  # The next group to compare
      
      # Create a new binary variable: Included (1) vs. Next Group (0)
      feature_data <- feature_data %>%
        mutate(
          Group_Binary = ifelse(group %in% included_groups, "Included", "Next_Group")
        )
      
      # Convert to factor
      feature_data$Group_Binary <- as.factor(feature_data$Group_Binary)
      
      # Fit multinomial logistic regression for `Included vs. Next_Group`, controlling for Age, Sex, and PMD.hs
      formula <- as.formula("Group_Binary ~ value + Age + Sex + PMD.hs")
      model <- multinom(formula, data = feature_data, trace = FALSE)
      
      # Perform Tukey’s post-hoc test comparing `Included vs. Next_Group`
      tukey_results <- emmeans(model, pairwise ~ Group_Binary, adjust = "tukey")$contrasts %>%
        as.data.frame()
      
      # Extract p-value for `Included vs. Next_Group`
      p_value <- tukey_results$p.value[1]  # Extract first p-value
      
      # Store the result
      tukey_rec_results <- rbind(
        tukey_rec_results,
        data.frame(
          Feature = feature,
          Group_Comparison = paste0(paste(included_groups, collapse = "+"), " vs. ", next_group),
          P_Value = p_value
        )
      )
    }
    
    # Identify the first significant comparison (p < cutoff)
    first_significant <- tukey_rec_results %>%
      dplyr::filter(P_Value < pval_cutoff) %>%
      dplyr::slice(1)
    
    # Add group means to the results
    group_means <- group_means %>%
      mutate(Feature = feature) %>%
      dplyr::rename(Group = group, Mean_Value = mean_value)
    
    # Combine the results for this feature
    feature_results <- merge(tukey_rec_results, group_means, by = "Feature", all = TRUE)
    feature_results$First_Significant <- ifelse(
      feature_results$Group_Comparison == first_significant$Group_Comparison,
      TRUE, FALSE
    )
    
    results[[feature]] <- feature_results
  }
  
  final_results <- do.call(rbind, results)
  
  return(final_results)
}


## TH
# format column names and object
colnames(gxdat_s_ctr@meta.data) <- make.names(colnames(gxdat_s_ctr@meta.data))
obj <- subset(gxdat_s_ctr, subset = segment == "TH")

# Run reciprocal Tukey’s test while controlling for Age and Sex
tukey_results <- ReciprocalTukeyTest(obj = obj, features = da_neurons, group.by = "ROI", pval_cutoff = 0.05)

# Save results as a table
tukey_results <- tukey_results[tukey_results$P_Value < 0.05,]
tukey_results <- tukey_results[!duplicated(tukey_results$Feature),]
write.csv(tukey_results, "reciprocal_logreg_tuk_results_th.csv", row.names = FALSE)

# 
# 
# ## Full ROI
# # format column names and object
# colnames(gxdat_s_ctr@meta.data) <- make.names(colnames(gxdat_s_ctr@meta.data))
# obj <- subset(gxdat_s_ctr, subset = segment != "TH")
# 
# # Run reciprocal Tukey’s test while controlling for Age and Sex
# tukey_results <- ReciprocalTukeyTest(obj = obj, features = da_neurons, group.by = "ROI", pval_cutoff = 0.05)
# 
# # Save results as a table
# tukey_results <- tukey_results[tukey_results$P_Value < 0.05,]
# tukey_results <- tukey_results[!duplicated(tukey_results$Feature),]
# write.csv(tukey_results, "reciprocal_logreg_tuk_results_full.csv", row.names = FALSE)


#------------ Group-mean to ID location ------------

# function for dotplot - note ordering is same as Visium
DotPlotSingleDotOrdered <- function(obj, features, group.by) {
  # Custom row order for features
  custom_row_order <- c(
    "CALB1_CRYM",
    "SOX6_AGTR1_NOX4",
    "SOX6_SMOC1_LPL",
    "CALB1_SEMA3D_RSPO3",
    "CALB1_VIP_NPPC",
    "CALB1_CRYM_CALCR",
    "CALB1_NEUROD6_PPP1R17",
    "GAD2_EBF2_NPSR1",
    "SOX6_GFRA2_TBC1D8B",
    "CALB1_CBLN4_PAX5",
    "CALB1_PRLHR_TRHR",
    "GAD2_CALCRL_KCNK13"
  )
  
  # Custom column order for groups
  custom_col_order <- c("SNM", "SNV", "SNL", "SND", "VTA", "RN")
  
  # Extract data
  data <- FetchData(obj, vars = c(group.by, features))
  colnames(data)[1] <- "group"  # Standardize column name for group.by
  
  # Convert features to long format
  long_data <- data %>%
    pivot_longer(cols = features, names_to = "feature", values_to = "value")
  
  # Calculate mean values for each group and feature
  mean_data <- long_data %>%
    group_by(group, feature) %>%
    dplyr::summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  
  # Identify the group with the highest mean for each feature
  highest_means <- mean_data %>%
    group_by(feature) %>%
    dplyr::filter(mean_value == max(mean_value)) %>%  # Keep only the highest mean per feature
    ungroup()
  
  # Add a column to mark highlights based on the highest mean per feature
  highlights <- mean_data %>%
    mutate(highlight = paste0(group, feature) %in% paste0(highest_means$group, highest_means$feature))
  
  # Reorder rows (features) based on the custom row order
  highlights <- highlights %>%
    mutate(feature = factor(feature, levels = custom_row_order, ordered = TRUE))
  
  # Reorder columns (groups) based on the custom column order
  highlights <- highlights %>%
    mutate(group = factor(group, levels = custom_col_order, ordered = TRUE))
  
  # Create the DotPlot
  p <- ggplot(highlights, aes(x = group, y = feature, size = mean_value, fill = highlight)) +
    geom_point(shape = 21, color = "black", stroke = 0.5) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey"), name = "Highest Mean") +
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

## TH-masked
obj <- subset(gxdat_s_ctr, subset = segment == "TH")

# Create DotPlot with Highlighted Significant Regions
pdf("dotplot_highlighted_highest_gx_TH.pdf", width = 12, height = 6)
DotPlotSingleDotOrdered(obj = obj, features = da_neurons, group.by = "ROI")
dev.off()


## Full ROI
obj <- subset(gxdat_s_ctr, subset = segment != "TH")

# Create DotPlot with Highlighted Significant Regions
pdf("dotplot_highlighted_highest_gx_full.pdf", width = 12, height = 6)
DotPlotSingleDotOrdered(obj = obj, features = da_neurons, group.by = "ROI")
dev.off()




#------------ Boxplots of t------------

## plot
dplot <- proportions_total
g2 <- ggplot(dplot, aes(x = Dx, y = proportion, fill = cell_type_merge )) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ cell_type_merge, nrow = 1, scales = "free_x") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
    strip.text.x = element_text(size = 5, angle = 90, hjust = 0, vjust = 0.5)  # Align text start to the right
  ) +
  labs(x = "Diagnosis", y = "Proportion of Cell Type", fill = "Cell", title = "") +
  scale_fill_manual(values = DA_col)






